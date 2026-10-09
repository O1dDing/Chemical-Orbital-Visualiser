#include "cov/nbo.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
void write(const std::filesystem::path& p,const std::string& s){std::ofstream f(p);f<<s;if(!f)throw std::runtime_error("fixture write failed");}
template<class F> void rejects(F f,const char* why){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,why);}
}
int main(int argc,char** argv){try{
    const auto dir=std::filesystem::temp_directory_path()/"cov_nbo_smoke";std::filesystem::create_directories(dir);
    const std::string text=" *********************************** NBO 7.0 ***********************************\n Cite this program [NBO 7.0.10 (8-Feb-2021)] as:\n NATURAL POPULATIONS: Natural atomic orbital occupancies\n 1 H 1 s Val( 1s) 1.00000 -0.50000\n 2 H 2 s Val( 1s) 1.00000 -0.50000\n Summary of Natural Population Analysis:\n H 1 0.00000 0.00000 1.00000 0.00000 1.00000\n H 2 0.00000 0.00000 1.00000 0.00000 1.00000\n (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (2.00000) BD ( 1) H 1- H 2\n ( 50.00%) 0.7071* H 1 s(100.00%)\n ( 50.00%) 0.7071* H 2 s(100.00%)\n 2. (0.00000) BD*( 1) H 1- H 2\n SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS\n Threshold for printing: 0.50 kcal/mol\n 1. BD ( 1) H 1- H 2  2. BD*( 1) H 1- H 2 1.25 0.80 0.010\n NATURAL BOND ORBITALS (Summary):\n 1. BD ( 1) H 1- H 2 2.00000 -0.50000\n 2. BD*( 1) H 1- H 2 0.00000 0.50000\n Total Lewis 2.00000 (100.00%)\n NBO dipole matrix\n 1. BD ( 1) H 1- H 2 0.00 0.00 -0.00\n Wiberg bond index matrix in the NAO basis:\n Atom 1 2\n 1. H 0.0000 1.0000\n 2. H 1.0000 0.0000\n Wiberg bond index, Totals by atom:\n";
    auto p=dir/"fixture.nbo";write(p,text);auto d=cov::read_nbo(p);
    require(d.naos.size()==2&&d.populations.size()==2&&d.orbitals.size()==2,"NPA/NAO/NBO counts");require(d.orbitals[0].components.size()==2,"local components");require(d.orbitals[0].diagonal_fock_hartree&&*d.orbitals[0].diagonal_fock_hartree==-0.5,"diagonal Fock semantics");require(d.e2.size()==1&&d.e2[0].printing_threshold==0.5&&d.e2[0].donor==1&&d.e2[0].acceptor==2,"E2 threshold/identity");require(d.wiberg.size()==4,"Wiberg full matrix");require(d.orbitals[0].source.line_begin>0,"source line provenance");
    write(p,text+text);rejects([&]{cov::read_nbo(p);},"ambiguous analysis accepted");cov::NboReadOptions select;select.analysis_segment=1;require(cov::read_nbo(p,select).orbitals.size()==2,"explicit selection");
    write(p,"NATURAL POPULATIONS: Natural atomic orbital occupancies\n 1 H 1 s Val( 1s) 1.0 -0.5\n");auto missing=cov::read_nbo(p);require(missing.e2.empty()&&!missing.e2_sections[0].printing_threshold&&!missing.e2_sections[0].missing_reason.empty(),"missing E2 became zero");
    write(p,"CMO: NBO Analysis of Canonical Molecular Orbitals\n Leading (> 5%) NBO Contributions to Molecular Orbitals\n MO 1 (occ): orbital energy = -1.0 a.u.\n 0.99*[ 1]: CR (1) H 1\n");auto cmo=cov::read_nbo(p);require(cmo.matrices.empty()&&cmo.cmo_summaries.size()==1&&!cmo.warnings.empty(),"CMO summary promoted to full decomposition");
    write(p,"******* Alpha spin orbitals *******\n (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (1.00000) LP ( 1) H 1 s(100.00%)\n******* Beta  spin orbitals *******\n (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (0.00000) LP ( 1) H 1 s(100.00%)\n");auto spins=cov::read_nbo(p);require(spins.orbitals.size()==2&&spins.orbitals[0].spin==cov::NboSpin::Alpha&&spins.orbitals[1].spin==cov::NboSpin::Beta,"spin identity reset");
    write(p," (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (1.98000) 3C ( 1) B 1- H 2- B 3\n ( 30.00%) 0.5477* B 1 s(50.00%)p 1.0(50.00%)\n ( 40.00%) 0.6325* H 2 s(100.00%)\n ( 30.00%) 0.5477* B 3 s(50.00%)p 1.0(50.00%)\n");auto three=cov::read_nbo(p);require(three.orbitals.size()==1&&three.orbitals[0].atoms.size()==3&&three.orbitals[0].components.size()==3,"general three-center parsing");
    cov::Wavefunction w;w.source=cov::WavefunctionSource::Fchk;w.basis_count=1;w.atoms.push_back({"He",2,0,0,0,2});w.shells.push_back({0,0,1,0,0,0});w.primitives.push_back({1.0,1.0});w.gaussian_ao_transform.push_back({0,1,1});w.ao_overlap={1};w.alpha_electrons=w.beta_electrons=1;w.electron_counts_provenance=cov::DataProvenance::Producer;w.orbital_occupation_model=cov::OrbitalOccupationModel::CanonicalShared;cov::MolecularOrbital mo;mo.occupation=2;mo.source_orbital_index=0;mo.coefficients={1};mo.gaussian_source_coefficients={1};w.orbitals.push_back(mo);
    cov::NboDataset one;cov::NboArchive a;a.atoms=w.atoms;a.basis_count=1;a.density_is_bond_order=true;a.centers={1};a.labels={1};a.ncomp={1};a.nprim={1};a.nptr={1};a.exponents={1};a.cs={1};for(const auto& k:{"OVERLAP","LCAOMO","DENSITY"}){cov::NboMatrix m;m.kind=k;m.rows=m.columns=1;m.values={std::string(k)=="DENSITY"?2.0:1.0};a.matrices.push_back(m);}one.archive=a;cov::NboOrbital no;no.id=1;no.occupation=2;one.orbitals={no};cov::NboMatrix b;b.kind="AONBO";b.rows=b.columns=1;b.values={1};one.matrices={b};require(cov::associate_nbo(one,w).compatible,"strict synthetic association");auto nw=cov::make_nbo_wavefunction(one,w);require(nw.orbitals.size()==1&&std::isnan(nw.orbitals[0].energy_hartree)&&w.orbitals[0].energy_hartree==0,"independent orbital energy semantics");auto bad=w;bad.atoms[0].x=0.1;require(!cov::associate_nbo(one,bad).compatible,"geometry mismatch accepted");bad=w;bad.primitives[0].exponent=2;require(!cov::associate_nbo(one,bad).compatible,"basis mismatch accepted");bad=w;bad.atoms[0].nuclear_charge=1;require(!cov::associate_nbo(one,bad).compatible,"ECP charge mismatch accepted");bad=w;bad.orbitals[0].occupation=1;require(!cov::associate_nbo(one,bad).compatible,"density mismatch accepted");
    // Some post-analysis FCHK files retain only alpha C but do retain both
    // producer densities. This supports beta NBOs without inventing beta C.
    cov::NboDataset open=one;open.archive->open_shell=true;open.archive->matrices.clear();open.matrices.clear();open.orbitals.clear();for(auto spin:{cov::NboSpin::Alpha,cov::NboSpin::Beta}){for(const auto& k:{"LCAOMO","DENSITY"}){cov::NboMatrix m;m.kind=k;m.spin=spin;m.rows=m.columns=1;m.values={1};open.archive->matrices.push_back(m);}auto nb=b;nb.spin=spin;open.matrices.push_back(nb);auto orbital=no;orbital.spin=spin;orbital.occupation=1;open.orbitals.push_back(orbital);}auto s=a.matrices[0];open.archive->matrices.push_back(s);auto density_w=w;density_w.total_density_packed={2};density_w.spin_density_packed={0};density_w.total_density_provenance=density_w.spin_density_provenance=cov::DataProvenance::Producer;auto evidence=cov::associate_nbo(open,density_w);require(evidence.compatible&&evidence.canonical_evidence.size()==2&&evidence.canonical_evidence[0].direct_fchk_coefficients&&!evidence.canonical_evidence[1].direct_fchk_coefficients&&evidence.canonical_evidence[1].density_verified,"beta density-supported provenance");require(cov::make_nbo_wavefunction(open,density_w).orbitals.size()==2,"density-supported beta rendering");density_w.spin_density_provenance=cov::DataProvenance::Derived;require(!cov::associate_nbo(open,density_w).compatible,"derived spin density admitted as independent evidence");
    write(p," (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (1.98000) 3C ( 1) B 1- B 2- H 3\n 2. (0.01000) 3Cn( 1) B 1- B 2- H 3\n 3. (0.01000) 3C*( 1) B 1- B 2- H 3\n");auto kinds=cov::read_nbo(p);require(kinds.orbitals.size()==3&&kinds.orbitals[1].kind=="3Cn","three-center nonbonding suffix identity");
    write(p,"NATURAL POPULATIONS: Natural atomic orbital occupancies\n NAO Atom No lang Type(AO) Occupancy Spin\n 1 C 1 s Val( 2s) 1.05 0.05\n Summary of Natural Population Analysis:\n C 1 -0.50000 2.00000 4.49000 0.01000 6.50000 1.00000\n ******* Alpha spin orbitals *******\n NATURAL POPULATIONS: Natural atomic orbital occupancies\n NAO Atom No lang Type(AO) Occupancy Energy\n 1 C 1 s Val( 2s) 0.55 -0.50\n");auto spin_columns=cov::read_nbo(p);require(spin_columns.populations.size()==1&&spin_columns.populations[0].spin==cov::NboSpin::Total&&spin_columns.populations[0].spin_density==1.0,"total NPA extra spin-density column");require(spin_columns.naos.size()==2&&spin_columns.naos[0].spin_density==0.05&&!spin_columns.naos[0].energy_hartree&&spin_columns.naos[1].energy_hartree==-0.5&&!spin_columns.naos[1].spin_density,"NAO spin density must not become energy");
    const std::string parents=" (Occupancy) Bond orbital / Coefficients / Hybrids\n 1. (1.00000) LP ( 1) H 1 s(100.00%)\n 2. (1.00000) LP ( 2) H 1 p(100.00%)\n";
    const std::string nlmo_head=" NLMO / Occupancy / Percent from Parent NBO\n";
    write(p,parents+nlmo_head+" 1. (1.00000) 99.0000% LP ( 2) H 1\n");
    auto parent=cov::read_nbo(p);require(parent.nlmos.size()==1&&parent.nlmos[0].parent_nbo==2&&
        parent.nlmos[0].parent_percent==99&&parent.nlmos[0].source.line_begin>0,
        "cached parent key changed actual parent identity or provenance");
    write(p,parents+nlmo_head+" 1. (1.00000) 99.0000% LP ( 3) H 1\n 2. (1.00000) 98.0000% unknown label\n");
    auto unknown_parent=cov::read_nbo(p);require(unknown_parent.nlmos.size()==2&&
        !unknown_parent.nlmos[0].parent_nbo&&!unknown_parent.nlmos[1].parent_nbo,"missing parent was guessed");
    write(p,parents+" 3. (1.00000) LP ( 1) H 1 p(100.00%)\n"+nlmo_head+
        " 1. (1.00000) 99.0000% LP ( 1) H 1\n");
    rejects([&]{cov::read_nbo(p);},"same-spin duplicate parent key lost ambiguity rejection");
    write(p,parents+nlmo_head+" 1. (1.00000) 99.0000% unknown label\n 1. (1.00000) 99.0000% unknown label\n");
    rejects([&]{cov::read_nbo(p);},"duplicate NLMO identity hidden by missing parent");
    const auto parent_block=parents+nlmo_head+" 1. (1.00000) 99.0000% LP ( 1) H 1\n";
    write(p," Alpha spin orbitals\n"+parent_block+" Beta spin orbitals\n"+parent_block);
    auto parent_spins=cov::read_nbo(p);require(parent_spins.nlmos.size()==2&&
        parent_spins.nlmos[0].parent_nbo==1&&parent_spins.nlmos[1].parent_nbo==1&&
        parent_spins.nlmos[0].spin!=parent_spins.nlmos[1].spin,"parent spin isolation lost");
    write(p,parents+" NATURAL BOND ORBITALS (Summary):\n Molecular unit 1\n"
        " 1. LP ( 1) H 1 1.00000 -0.50000\n Total Lewis 1.00000 (100.00%)\n"
        " Molecular unit 2\n 2. LP ( 2) H 1 1.00000 -0.25000\n"
        " Total Lewis 1.00000 (100.00%)\n NATURAL LOCALIZED MOLECULAR ORBITAL (NLMO) ANALYSIS:\n");
    auto units=cov::read_nbo(p);require(units.orbitals.size()==2&&
        units.orbitals[0].diagonal_fock_hartree==-.5&&units.orbitals[1].diagonal_fock_hartree==-.25,
        "molecular unit subtotal truncated the remaining NBO energy summary");
    const auto matrix_path=dir/"numeric-matrix.37";
    write(matrix_path,"NBOs in the AO basis:\n 999 integer metadata\n 1 2 3\n"
        "\t+.5D+1 1. -0\f\n 2e0\v-0.25d+1 +3E-2\r\n 0.0 0.0\n");
    const auto numeric_matrix=cov::read_nbo_matrix_rectangular(matrix_path,"AONBO",2,3,false);
    require(numeric_matrix.size()==1&&numeric_matrix[0].values==
        std::vector<double>{5.0,-0.0,-2.5,1.0,2.0,0.03}&&
        std::signbit(numeric_matrix[0].values[1])&&
        numeric_matrix[0].source.line_begin==4&&numeric_matrix[0].source.line_end==5,
        "matrix numeric scanner altered values, orientation, signed zero or source lines");
    for(const std::string bad_token:{"1e", ".", "+-1", "1,2", "NaN", "Inf", "0x1p0", "1.0-2.0", "1D+", "1e999bad"}) {
        write(matrix_path,"NBOs in the AO basis:\n 1.0 0.0\n "+bad_token+" 1.0\n");
        rejects([&]{cov::read_nbo_matrix(matrix_path,"AONBO",2,false);},
            "malformed or adjacent matrix number bypassed payload rejection");
    }
    write(matrix_path,"NBOs in the AO basis:\n 1e999bad\n .5 0.0 0.0 1.0\n");
    require(cov::read_nbo_matrix(matrix_path,"AONBO",2,false)[0].values[0]==.5,
        "malformed pre-payload header must be validated before numeric conversion");
    for(const auto* payload:{"1.0 0.0\n 0 1\n","1.0 0.0 0.0 1.0 2.0\n",
                             "1.0 0.0 0.0 1e999\n"}) {
        write(matrix_path,std::string("NBOs in the AO basis:\n")+payload);
        rejects([&]{cov::read_nbo_matrix(matrix_path,"AONBO",2,false);},
            "integer interruption, boundary overflow or numeric overflow was accepted");
    }
    if(argc>1){auto real=cov::read_nbo(argv[1]);require(!real.orbitals.empty()&&!real.naos.empty(),"real sample parse empty");std::cout<<"real_sample_orbitals="<<real.orbitals.size()<<" e2="<<real.e2.size()<<"\n";}
    std::cout<<"NBO smoke: passed provenance, ambiguity, spins, three-center, thresholds, CMO, association rejection and independent rendering contracts\n";return 0;
}catch(const std::exception& e){std::cerr<<"NBO smoke failed: "<<e.what()<<'\n';return 1;}}
