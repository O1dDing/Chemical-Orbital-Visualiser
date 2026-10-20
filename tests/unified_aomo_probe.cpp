#include "cov/nbo_salc.hpp"
#include "cov/wavefunction_io.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
std::string quote(const std::string& value) {
    std::string out="\"";
    for (const auto ch:value) {
        switch(ch) {
        case '\\':out+="\\\\";break;
        case '"':out+="\\\"";break;
        case '\n':out+="\\n";break;
        case '\r':out+="\\r";break;
        case '\t':out+="\\t";break;
        default:out+=ch;
        }
    }
    return out+'"';
}
void values(std::ostream& out,const std::vector<double>& data) {
    out<<'[';
    for(std::size_t i=0;i<data.size();++i) {
        if(i)out<<',';
        if(!std::isfinite(data[i]))throw std::runtime_error("Nonfinite scientific evidence");
        out<<data[i];
    }
    out<<']';
}
}

// This observation adapter exercises the same scientific model and selection
// path used by the UI. Acceptance is owned by an independent input-file oracle
// and visible native tests; this executable does not certify its own numbers.
int main(int argc,char** argv) {try {
    if(argc!=4)throw std::runtime_error("Usage: cov_unified_aomo_probe canonical.fchk analysis_directory output.json");
    const auto canonical=cov::parse_wavefunction(argv[1]);
    const auto original=canonical;
    const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1)throw std::runtime_error("Exactly one explicit analysis candidate is required");
    const auto integrated=cov::read_nbo_integration(canonical,found.candidates.front());
    const auto model=cov::build_nbo_salc_model(canonical,integrated);
    const auto output=std::filesystem::path(argv[3]);
    if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output);
    if(!out)throw std::runtime_error("Cannot write unified orbital evidence");
    out<<std::setprecision(17)<<"{\"schema\":1,\"salc\":"<<cov::serialize_nbo_salc_json(model)
       <<",\"canonical\":{\"basis_count\":"<<canonical.basis_count<<",\"orbitals\":[";
    for(std::size_t i=0;i<canonical.orbitals.size();++i) {
        const auto& mo=canonical.orbitals[i];
        if(i)out<<',';
        out<<"{\"index\":"<<i<<",\"global_index\":"<<i
           <<",\"source_index\":"<<mo.source_orbital_index
           <<",\"within_spin_index\":"<<mo.source_orbital_index
           <<",\"spin\":"<<quote(mo.spin==cov::Spin::Beta?"beta":"alpha")
           <<",\"energy_hartree\":"<<mo.energy_hartree<<",\"occupation\":"<<mo.occupation<<'}';
    }
    out<<"]},\"side_fields\":[";
    std::size_t selection_failures=0;
    for(std::size_t i=0;i<model.orbitals.size();++i) {
        const auto selection=cov::nbo_salc_selection(model,i);
        const auto view=cov::make_nbo_selection_view(integrated,canonical,selection);
        if(i)out<<',';
        out<<"{\"index\":"<<i<<",\"id\":"<<quote(model.orbitals[i].id)
           <<",\"available\":"<<(view.available?"true":"false")
           <<",\"detail\":"<<quote(view.detail)<<",\"metric_norm2\":";
        values(out,view.metric_norm2);
        out<<",\"coefficients\":";
        if(view.available&&view.wavefunction.orbitals.size()==1)values(out,view.wavefunction.orbitals.front().coefficients);
        else {out<<"[]";++selection_failures;}
        out<<'}';
    }
    bool preserved=canonical.ao_overlap==original.ao_overlap&&canonical.orbitals.size()==original.orbitals.size();
    for(std::size_t i=0;i<canonical.orbitals.size();++i) {
        const auto& now=canonical.orbitals[i];const auto& before=original.orbitals[i];
        preserved=preserved&&now.coefficients==before.coefficients
            &&now.gaussian_source_coefficients==before.gaussian_source_coefficients
            &&now.energy_hartree==before.energy_hartree&&now.occupation==before.occupation
            &&now.spin==before.spin&&now.source_orbital_index==before.source_orbital_index;
    }
    out<<"],\"selection_failures\":"<<selection_failures<<",\"canonical_preserved\":"<<(preserved?"true":"false")<<'}';
    out.close();if(!out)throw std::runtime_error("Evidence write failed");
    std::cout<<"unified model "<<model.status<<", side orbitals="<<model.orbitals.size()
             <<", canonical preserved="<<preserved<<", selection failures="<<selection_failures<<'\n';
    return preserved&&selection_failures==0?0:2;
}catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
