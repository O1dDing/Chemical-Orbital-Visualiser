#include "cov/atomic_orbital_symmetry.hpp"
#include "cov/local_angular_generators.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
int failures=0;
void check(bool passed,const std::string& message) {
    if(!passed) {std::cerr<<message<<'\n';++failures;}
}
cov::Wavefunction atom(int l,bool pure,int shells=1) {
    cov::Wavefunction w;w.atoms.push_back({"X",10,0,0,0});
    for(int i=0;i<shells;++i) {
        cov::Shell shell;shell.basis_offset=w.basis_count;
        shell.angular_momentum=static_cast<std::uint8_t>(l);shell.pure=pure?1:0;
        w.basis_count+=cov::shell_basis_count(shell);w.shells.push_back(shell);
    }
    const auto n=w.basis_count;w.ao_overlap.assign(n*n,0);
    for(std::size_t i=0;i<n;++i)w.ao_overlap[i*n+i]=1;
    return w;
}
double inner(const cov::Wavefunction& w,const std::vector<double>& a,const std::vector<double>& b) {
    double result=0;const std::size_t n=w.basis_count;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)result+=a[i]*w.ao_overlap[i*n+j]*b[j];
    return result;
}
void add(cov::Wavefunction& w,std::vector<double> c,double energy=-1,cov::Spin spin=cov::Spin::Alpha) {
    cov::MolecularOrbital mo;mo.coefficients=std::move(c);mo.energy_hartree=energy;
    mo.spin=spin;mo.occupation=2;w.orbitals.push_back(std::move(mo));
}
void rotated_full_basis(cov::Wavefunction& w,double energy=-1) {
    const std::size_t n=w.basis_count;
    std::vector<std::vector<double>> columns(n,std::vector<double>(n,0));
    for(std::size_t i=0;i<n;++i)columns[i][i]=1;
    for(std::size_t a=0;a<n;++a)for(std::size_t b=a+1;b<n;++b) {
        const double angle=.137*(a+1)*(b+2),c=std::cos(angle),s=std::sin(angle);
        for(std::size_t row=0;row<n;++row) {
            const double x=columns[a][row],y=columns[b][row];
            columns[a][row]=c*x+s*y;columns[b][row]=-s*x+c*y;
        }
    }
    for(auto& c:columns)add(w,c,energy);
}
void cartesian_metric(cov::Wavefunction& w) {
    // Equal radial envelope: normalized Cartesian d diagonals have <xx|yy>=1/3.
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)w.ao_overlap[i*6+j]=i==j?1:1.0/3;
}
std::vector<std::vector<double>> angular_columns(const cov::Wavefunction& w,int l) {
    const auto generators=cov::local_angular_generators(2,false,l,{1,0,0,0,1,0,0,0,1});
    std::vector<std::vector<double>> columns;
    for(std::size_t j=0;j<generators.columns;++j) {
        std::vector<double> c(w.basis_count);
        for(std::size_t i=0;i<c.size();++i)c[i]=generators.coefficients[i*generators.columns+j];
        // Gram-Schmidt in actual S; source columns are never normalized by the analyser.
        for(const auto& previous:columns) {
            const double overlap=inner(w,previous,c);
            for(std::size_t i=0;i<c.size();++i)c[i]-=overlap*previous[i];
        }
        const double norm=std::sqrt(inner(w,c,c));for(auto& v:c)v/=norm;
        columns.push_back(std::move(c));
    }
    return columns;
}
}

int main() {
    for(int l=1;l<=4;++l) {
        auto w=atom(l,l>=2);rotated_full_basis(w);
        const auto source=w.orbitals;const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(result.available,"rotated angular basis unavailable l="+std::to_string(l));
        check(result.symmetry_group=="SO(3)","single atom assigned a finite group");
        for(std::size_t i=0;i<result.orbitals.size();++i) {
            const auto& a=result.orbitals[i];
            check(a.angular_momentum_verified&&a.angular_momentum==l,"rotated pure l not recognized");
            check(a.partner_block_verified&&a.containing_span_closed,"full rotated 2l+1 partner span not recognized");
            check(a.radial_copy_ordinal==1&&a.principal_n==0,"radial copy fabricated principal n");
            check(a.inversion_parity==(l%2?-1:1),"parity inconsistent with l");
            check(source[i].coefficients==w.orbitals[i].coefficients&&source[i].energy_hartree==w.orbitals[i].energy_hartree&&
                  source[i].occupation==w.orbitals[i].occupation,"analysis altered source MO");
        }
    }
    {
        auto w=atom(1,false);w.ao_overlap={2,0,0,0,2,0,0,0,2};
        for(int m=0;m<3;++m) {std::vector<double> c(3);c[m]=1/std::sqrt(2.0);add(w,c);}
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(result.available&&result.orbitals[0].partner_block_verified,"nonidentity S normalization rejected");
        check(std::abs(result.orbitals[0].source_squared_s_norm-1)<1e-12,"S norm replaced by coefficient squared norm");
        w.orbitals[0].coefficients[0]=1;
        const auto bad=cov::analyse_atomic_orbital_symmetry(w);
        check(!bad.orbitals[0].angular_momentum_verified&&bad.orbitals[0].status=="source_metric_not_normalized","unnormalized source hard-assigned");
    }
    {
        auto w=atom(0,false);cov::Shell shell;shell.basis_offset=1;shell.angular_momentum=1;
        w.shells.push_back(shell);w.basis_count=4;w.ao_overlap.assign(16,0);
        for(int i=0;i<4;++i)w.ao_overlap[4*i+i]=1;
        add(w,{std::sqrt(.5),std::sqrt(.5),0,0});
        const auto result=cov::analyse_atomic_orbital_symmetry(w);const auto& a=result.orbitals[0];
        check(result.available&&a.decomposition_verified&&!a.angular_momentum_verified&&a.status=="mixed_angular_momentum","different-l mixture hard-assigned");
        check(a.components.size()==2&&std::abs(a.components[0].weight-.5)<1e-12&&std::abs(a.components[1].weight-.5)<1e-12,"different-l S weights incorrect");
    }
    {
        auto w=atom(2,false);cartesian_metric(w);
        auto d=angular_columns(w,2),s=angular_columns(w,0);
        for(auto c:d)add(w,c);add(w,s.front(),-.5);
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(result.available,"Cartesian d angular decomposition unavailable");
        if(result.orbitals.size()==6) {
            for(int i=0;i<5;++i)check(result.orbitals[i].angular_momentum==2&&result.orbitals[i].partner_block_verified,"Cartesian pure d combination rejected");
            check(result.orbitals[5].angular_momentum==0&&result.orbitals[5].partner_block_verified,"Cartesian trace direction mislabeled d");
        }
        w.orbitals.clear();add(w,{1,0,0,0,0,0});
        const auto mixed=cov::analyse_atomic_orbital_symmetry(w);
        check(mixed.available&&mixed.orbitals[0].decomposition_verified&&!mixed.orbitals[0].angular_momentum_verified,"Cartesian xx lower-l pollution hard-assigned d");
        if(mixed.orbitals[0].components.size()==2) {
            check(std::abs(mixed.orbitals[0].components[0].weight-5.0/9)<1e-12,"Cartesian xx s weight should be 5/9");
            check(std::abs(mixed.orbitals[0].components[1].weight-4.0/9)<1e-12,"Cartesian xx d weight should be 4/9");
        }
    }
    {
        auto w=atom(1,false);add(w,{1,0,0});
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(result.orbitals[0].angular_momentum_verified&&!result.orbitals[0].partner_block_verified,"one p column promoted to complete shell");
    }
    {
        auto w=atom(1,false,2);add(w,{1,0,0,0,0,0});add(w,{0,0,0,0,1,0});add(w,{0,0,1,0,0,0});
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(const auto& a:result.orbitals)check(a.angular_momentum_verified&&!a.containing_span_closed&&a.generator_leakage_squared>.5,"radially inconsistent 3D p span accepted");
    }
    {
        auto w=atom(1,false,2);rotated_full_basis(w);
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(const auto& a:result.orbitals)check(a.angular_momentum_verified&&a.containing_span_closed&&a.representation_multiplicity==2&&
            !a.partner_block_verified&&a.radial_copy_ordinal==0,"mixed repeated radial copies assigned individual partners");
    }
    {
        auto w=atom(1,false,2);const double overlap=.35;
        for(int m=0;m<3;++m) {w.ao_overlap[m*6+3+m]=overlap;w.ao_overlap[(3+m)*6+m]=overlap;}
        for(int m=0;m<3;++m) {std::vector<double> c(6);c[m]=1;add(w,c,-1);}
        for(int m=0;m<3;++m) {
            std::vector<double> c(6);c[m]=-overlap/std::sqrt(1-overlap*overlap);
            c[3+m]=1/std::sqrt(1-overlap*overlap);add(w,c,.2);
        }
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(result.available&&result.angular_reference_ranks[1]==6,"overlapping radial references lost angular rank");
        for(std::size_t i=0;i<6;++i)check(result.orbitals[i].partner_block_verified&&
            result.orbitals[i].radial_copy_ordinal==(i<3?1u:2u),"overlapping radial copies not identified in S");
    }
    {
        auto w=atom(1,false);add(w,{1,0,0},-1,cov::Spin::Alpha);add(w,{0,1,0},-1,cov::Spin::Beta);add(w,{0,0,1},-1,cov::Spin::Beta);
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(const auto& a:result.orbitals)check(!a.partner_block_verified,"different spins counted as spatial partners");
    }
    {
        auto w=atom(1,false);rotated_full_basis(w);w.atoms[0].atomic_number=18;w.atoms[0].nuclear_charge=8;
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(const auto& a:result.orbitals)check(a.angular_momentum_verified&&a.principal_n==0,"ECP inferred missing principal/core shells");
        w.atoms.push_back({"H",1,1,0,0});
        const auto molecule=cov::analyse_atomic_orbital_symmetry(w);
        check(!molecule.applicable&&!molecule.available&&molecule.status=="not_single_atom","multiatom molecule accepted as atom");
    }
    {
        auto w=atom(1,false);rotated_full_basis(w);w.ao_overlap[0]=2;
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        check(!result.available&&result.status=="atomic_metric_not_rotation_invariant","anisotropic AO metric accepted as spherical");
        w=atom(1,false);w.shells[0].basis_offset=1;
        check(!cov::analyse_atomic_orbital_symmetry(w).available,"invalid AO shell range accepted");
    }
    {
        // Even an already closed higher p copy cannot be numbered if an earlier
        // mixed source could contain another p occurrence.
        auto w=atom(1,false,2);
        cov::Shell shell;shell.basis_offset=6;shell.angular_momentum=0;
        w.shells.push_back(shell);w.basis_count=7;w.ao_overlap.assign(49,0);
        for(int i=0;i<7;++i)w.ao_overlap[7*i+i]=1;
        add(w,{std::sqrt(.5),0,0,0,0,0,std::sqrt(.5)},-2);
        add(w,{0,0,0,1,0,0,0},-1);add(w,{0,0,0,0,1,0,0},-1);add(w,{0,0,0,0,0,1,0},-1);
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(std::size_t i=1;i<4;++i)check(result.orbitals[i].partner_block_verified&&result.orbitals[i].radial_copy_ordinal==0,
            "incomplete mixed lower occurrence ignored in radial ordinal");
    }
    {
        auto w=atom(0,false);add(w,{1},-2);add(w,{1},-1);
        const auto result=cov::analyse_atomic_orbital_symmetry(w);
        for(const auto& a:result.orbitals)check(a.partner_block_verified&&a.radial_copy_ordinal==0,
            "duplicated source s columns counted as separate radial occurrences");
    }
    if(failures) {std::cerr<<failures<<" atomic symmetry checks failed\n";return 1;}
    std::cout<<"atomic_orbital_symmetry_smoke passed\n";
}
