#include "cov/nbo.hpp"
#include "cov/wavefunction_io.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Production-parser observation only. The external checker owns all verdicts.
int main(int argc, char** argv) {
    try {
        if(argc!=7 && argc!=10) throw std::runtime_error("usage: cov_nbo_probe canonical.fchk output.nbo archive.47 aonbo.37 nbomo.49 output-directory [naomo aonao naonbo]");
        const auto out=std::filesystem::u8path(argv[6]);
        if(std::filesystem::exists(out)) throw std::runtime_error("refusing to overwrite probe evidence");
        std::filesystem::create_directories(out);
        cov::WavefunctionParseOptions parse_options;
        parse_options.auto_enrich_gaussian_log=false;
        auto canonical=cov::parse_wavefunction(std::filesystem::u8path(argv[1]),parse_options);
        const auto original_orbitals=canonical.orbitals;
        cov::NboReadOptions options;
        if(std::string(argv[3])!="-")options.archive47=std::filesystem::u8path(argv[3]);
        if(std::string(argv[4])!="-")options.aonbo=std::filesystem::u8path(argv[4]);
        if(std::string(argv[5])!="-")options.nbomo=std::filesystem::u8path(argv[5]);
        if(argc==10) {
            if(std::string(argv[7])!="-")options.naomo=std::filesystem::u8path(argv[7]);
            if(std::string(argv[8])!="-")options.aonao=std::filesystem::u8path(argv[8]);
            if(std::string(argv[9])!="-")options.naonbo=std::filesystem::u8path(argv[9]);
        }
        auto dataset=cov::read_nbo(std::filesystem::u8path(argv[2]),options);
        const auto association=cov::associate_nbo(dataset,canonical);
        bool canonical_preserved=original_orbitals.size()==canonical.orbitals.size();
        for(std::size_t i=0;canonical_preserved && i<original_orbitals.size();++i) {
            const auto& before=original_orbitals[i];const auto& after=canonical.orbitals[i];
            canonical_preserved=before.coefficients==after.coefficients &&
                before.gaussian_source_coefficients==after.gaussian_source_coefficients &&
                (before.energy_hartree==after.energy_hartree ||
                 (std::isnan(before.energy_hartree)&&std::isnan(after.energy_hartree))) &&
                before.occupation==after.occupation && before.spin==after.spin &&
                before.source_orbital_index==after.source_orbital_index && before.symmetry==after.symmetry;
        }
        {std::ofstream f(out/"canonical-identity.json");f<<std::setprecision(17);
         f<<"{\"preserved\":"<<(canonical_preserved?"true":"false")<<",\"orbitals\":[";
         for(std::size_t i=0;i<canonical.orbitals.size();++i){if(i)f<<',';const auto& o=canonical.orbitals[i];
             f<<"{\"index\":"<<i<<",\"source_index\":"<<o.source_orbital_index
              <<",\"spin\":"<<static_cast<int>(o.spin)<<",\"occupation\":"<<o.occupation<<",\"energy_hartree\":";
             if(std::isfinite(o.energy_hartree))f<<o.energy_hartree;else f<<"null";f<<'}';}
         f<<"]}";if(!f)throw std::runtime_error("canonical identity write failed");}
        if(!canonical_preserved)throw std::runtime_error("association modified canonical orbitals");
        {std::ofstream f(out/"dataset.json");f<<cov::serialize_nbo_json(dataset);if(!f)throw std::runtime_error("dataset write failed");}
        if(!association.compatible) {
            std::cerr<<association.status<<": "<<association.detail<<'\n';
            return 2;
        }
        const auto local=cov::make_nbo_wavefunction(dataset,canonical);
        std::ofstream f(out/"render-contract.json");f<<std::setprecision(17);
        f<<"{\"schema\":1,\"basis_count\":"<<local.basis_count<<",\"orbitals\":[";
        for(std::size_t i=0;i<local.orbitals.size();++i) {
            if(i)f<<',';
            const auto& o=local.orbitals[i];
            f<<"{\"index\":"<<i<<",\"spin\":"<<static_cast<int>(o.spin)<<",\"source_index\":"<<o.source_orbital_index<<",\"occupation\":"<<o.occupation<<",\"canonical_energy\":";
            if(std::isfinite(o.energy_hartree))f<<o.energy_hartree;else f<<"null";
            f<<",\"gaussian_coefficients\":[";
            for(std::size_t k=0;k<o.gaussian_source_coefficients.size();++k){if(k)f<<',';f<<o.gaussian_source_coefficients[k];}
            f<<"],\"internal_coefficients\":[";
            for(std::size_t k=0;k<o.coefficients.size();++k){if(k)f<<',';f<<o.coefficients[k];}
            f<<"]}";
        }
        f<<"]}";
        if(!f)throw std::runtime_error("render contract write failed");
        std::cout<<"observed "<<local.orbitals.size()<<" NBO orbitals; external checks required\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
