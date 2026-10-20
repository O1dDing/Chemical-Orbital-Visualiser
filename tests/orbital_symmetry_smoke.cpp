#include "cov/orbital_symmetry.hpp"

#include <cmath>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

cov::Atom atom(const char* symbol, int z, double x, double y, double zc) {
    return {symbol,z,x,y,zc};
}

void identity_overlap(cov::Wavefunction& wf) {
    wf.ao_overlap.assign(static_cast<std::size_t>(wf.basis_count)*wf.basis_count,0.0);
    for (std::size_t i=0;i<wf.basis_count;++i) wf.ao_overlap[i*wf.basis_count+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;
}

void add_s_shell(cov::Wavefunction& wf, std::size_t atom_index) {
    const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back({1.0f,1.0f});
    cov::Shell shell;
    shell.atom_index=static_cast<std::uint32_t>(atom_index);
    shell.primitive_offset=primitive_offset;
    shell.primitive_count=1;
    shell.basis_offset=wf.basis_count;
    shell.angular_momentum=0;
    shell.pure=0;
    wf.shells.push_back(shell);
    ++wf.basis_count;
}

void add_p_shell(cov::Wavefunction& wf, std::size_t atom_index) {
    const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back({1.0f,1.0f});
    cov::Shell shell;
    shell.atom_index=static_cast<std::uint32_t>(atom_index);
    shell.primitive_offset=primitive_offset;
    shell.primitive_count=1;
    shell.basis_offset=wf.basis_count;
    shell.angular_momentum=1;
    shell.pure=0;
    wf.shells.push_back(shell);
    wf.basis_count+=3;
}

void add_mo(cov::Wavefunction& wf, double energy, std::vector<double> c) {
    cov::MolecularOrbital mo;
    mo.energy_hartree=energy;
    mo.spin=cov::Spin::Alpha;
    mo.occupation=0.0f;
    for (double value:c) mo.coefficients.push_back(value);
    wf.orbitals.push_back(std::move(mo));
}

bool all_label(const cov::Wavefunction& wf,
               std::size_t begin,
               std::size_t end,
               const std::string& expected) {
    for (std::size_t i=begin;i<end;++i) {
        if (wf.orbitals[i].symmetry!=expected ||
            wf.orbitals[i].symmetry_provenance!=cov::DataProvenance::Derived) {
            std::cerr << "expected " << expected << " at MO " << i
                      << ", got '" << wf.orbitals[i].symmetry << "'\n";
            return false;
        }
    }
    return true;
}

cov::Wavefunction d3h_ring() {
    cov::Wavefunction wf;
    const double r=1.6;
    for (int i=0;i<3;++i) {
        const double a=2.0*pi*static_cast<double>(i)/3.0;
        wf.atoms.push_back(atom("H",1,r*std::cos(a),r*std::sin(a),0.0));
        add_s_shell(wf,static_cast<std::size_t>(i));
    }
    const double s3=std::sqrt(3.0), s6=std::sqrt(6.0), s2=std::sqrt(2.0);
    add_mo(wf,-1.0,{1.0/s3,1.0/s3,1.0/s3});
    add_mo(wf,-0.4,{2.0/s6,-1.0/s6,-1.0/s6});
    add_mo(wf,-0.4,{0.0,1.0/s2,-1.0/s2});
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction d5h_ring() {
    cov::Wavefunction wf;
    const double r=2.0;
    for (int i=0;i<5;++i) {
        const double a=2.0*pi*static_cast<double>(i)/5.0;
        wf.atoms.push_back(atom("C",6,r*std::cos(a),r*std::sin(a),0.0));
        add_s_shell(wf,static_cast<std::size_t>(i));
    }
    std::vector<double> a1(5,1.0/std::sqrt(5.0));
    std::vector<double> c1(5),s1(5),c2(5),s2(5);
    const double norm=std::sqrt(2.0/5.0);
    for (int i=0;i<5;++i) {
        const double a=2.0*pi*static_cast<double>(i)/5.0;
        c1[i]=norm*std::cos(a); s1[i]=norm*std::sin(a);
        c2[i]=norm*std::cos(2.0*a); s2[i]=norm*std::sin(2.0*a);
    }
    add_mo(wf,-1.0,a1);
    add_mo(wf,-0.5,c1); add_mo(wf,-0.5,s1);
    add_mo(wf, 0.1,c2); add_mo(wf, 0.1,s2);
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction even_ring(int n, bool perpendicular_p, bool change_representation) {
    cov::Wavefunction wf;
    for (int i=0;i<n;++i) {
        const double angle=2*pi*i/n;
        wf.atoms.push_back(atom("C",6,2*std::cos(angle),2*std::sin(angle),0));
        if (perpendicular_p) add_p_shell(wf,i);
        else add_s_shell(wf,i);
    }
    // Analytic discrete Fourier irreps in an orthonormal site basis. This
    // representation fixture does not claim to be a Gaussian wavefunction.
    const auto add_fourier=[&](int k,bool sine,double energy) {
        std::vector<double> values(wf.basis_count);
        const double scale=std::sqrt((k==0 || k==n/2)?1.0/n:2.0/n);
        for (int i=0;i<n;++i) {
            const double angle=2*pi*k*i/n;
            values[perpendicular_p?3*i+2:i]=scale*(sine?std::sin(angle):std::cos(angle));
        }
        add_mo(wf,energy,values);
    };
    add_fourier(0,false,-1.0);
    for (int k=1;k<n/2;++k) {
        add_fourier(k,false,-1+.3*k);add_fourier(k,true,-1+.3*k);
    }
    add_fourier(n/2,false,1.0);
    identity_overlap(wf);
    if (change_representation) {
        const std::array<double,3> axis{1/std::sqrt(14.0),2/std::sqrt(14.0),3/std::sqrt(14.0)};
        const double c=std::cos(.47),s=std::sin(.47);
        std::array<std::array<double,3>,3> rotation{};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            rotation[i][j]=(i==j?c:0)+(1-c)*axis[i]*axis[j];
        rotation[0][1]-=s*axis[2];rotation[1][0]+=s*axis[2];
        rotation[0][2]+=s*axis[1];rotation[2][0]-=s*axis[1];
        rotation[1][2]-=s*axis[0];rotation[2][1]+=s*axis[0];
        for (auto& a:wf.atoms) {
            const std::array<double,3> old{a.x,a.y,a.z};
            std::array<double,3> moved{2,-3,4};
            for (int i=0;i<3;++i) for (int j=0;j<3;++j) moved[i]+=rotation[i][j]*old[j];
            a.x=moved[0];a.y=moved[1];a.z=moved[2];
        }
        if (perpendicular_p) for (auto& mo:wf.orbitals) for (int a=0;a<n;++a) {
            const double z=mo.coefficients[3*a+2];
            for (int i=0;i<3;++i) mo.coefficients[3*a+i]=rotation[i][2]*z;
        }
        // The E subspace is unchanged by a different canonical basis or by
        // the arbitrary overall phase of any individual real MO.
        auto first=wf.orbitals[1].coefficients,second=wf.orbitals[2].coefficients;
        for (std::size_t i=0;i<first.size();++i) {
            wf.orbitals[1].coefficients[i]=std::cos(.713)*first[i]+std::sin(.713)*second[i];
            wf.orbitals[2].coefficients[i]=-std::sin(.713)*first[i]+std::cos(.713)*second[i];
        }
        for (std::size_t i=0;i<wf.orbitals.size();i+=2)
            for (auto& value:wf.orbitals[i].coefficients) value=-value;
        std::reverse(wf.atoms.begin(),wf.atoms.end());
        for (auto& shell:wf.shells) shell.atom_index=n-1-shell.atom_index;
    }
    return wf;
}

cov::Wavefunction td_ligand_sigma() {
    cov::Wavefunction wf;
    wf.atoms={
        atom("Cl",17, 1, 1, 1), atom("Cl",17, 1,-1,-1),
        atom("Cl",17,-1, 1,-1), atom("Cl",17,-1,-1, 1),
    };
    for (std::size_t i=0;i<4;++i) add_s_shell(wf,i);
    add_mo(wf,-1.0,{0.5,0.5,0.5,0.5});
    add_mo(wf,-0.2,{1/std::sqrt(2.0),-1/std::sqrt(2.0),0,0});
    add_mo(wf,-0.2,{1/std::sqrt(6.0),1/std::sqrt(6.0),-2/std::sqrt(6.0),0});
    add_mo(wf,-0.2,{1/std::sqrt(12.0),1/std::sqrt(12.0),1/std::sqrt(12.0),-3/std::sqrt(12.0)});
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction oh_central_p() {
    cov::Wavefunction wf;
    wf.atoms={
        atom("Ti",22,0,0,0),
        atom("F",9, 2,0,0),atom("F",9,-2,0,0),
        atom("F",9,0, 2,0),atom("F",9,0,-2,0),
        atom("F",9,0,0, 2),atom("F",9,0,0,-2),
    };
    add_p_shell(wf,0);
    add_mo(wf,-0.3,{1,0,0});
    add_mo(wf,-0.3,{0,1,0});
    add_mo(wf,-0.3,{0,0,1});
    identity_overlap(wf);
    return wf;
}

} // namespace

int main() {
    // A square site basis plus an off-plane, different element has C4v.
    // The alternating site combination is B1 relative to atom-supported
    // mirrors, including after a general rigid rotation and atom reversal.
    for(bool rotated:{false,true}){
        auto wf=even_ring(4,false,rotated);std::array<double,3> centre{};
        for(const auto& a:wf.atoms){centre[0]+=a.x/4;centre[1]+=a.y/4;centre[2]+=a.z/4;}
        const std::array<double,3> a{wf.atoms[0].x-centre[0],wf.atoms[0].y-centre[1],wf.atoms[0].z-centre[2]},b{wf.atoms[1].x-centre[0],wf.atoms[1].y-centre[1],wf.atoms[1].z-centre[2]};
        const std::array<double,3> normal{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
        wf.atoms.push_back(atom("N",7,centre[0]+normal[0]/2,centre[1]+normal[1]/2,centre[2]+normal[2]/2));
        const auto result=cov::derive_orbital_symmetry(wf);
        if(result.point_group!="C4v"||result.orbitals_labelled!=4){std::cerr<<"C4v finite classifier incomplete\n";return 1;}
        if(!all_label(wf,0,1,"A1")||!all_label(wf,1,3,"E")||!all_label(wf,3,4,"B1"))return 1;
    }
    {
        auto wf=d3h_ring();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="D3h" || r.orbitals_labelled!=3u ||
            !all_label(wf,0,1,"A1'") || !all_label(wf,1,3,"E'")) {
            std::cerr << "D3h derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    for (const int n:{4,6}) for (const bool perpendicular_p:{false,true})
        for (const bool change_representation:{false,true}) {
        auto wf=even_ring(n,perpendicular_p,change_representation);
        const auto result=cov::derive_orbital_symmetry(wf);
        const std::vector<std::string> expected=n==4
            ?(perpendicular_p?std::vector<std::string>{"A2u","Eg","Eg","B2u"}
                             :std::vector<std::string>{"A1g","Eu","Eu","B1g"})
            :(perpendicular_p?std::vector<std::string>{"A2u","E1g","E1g","E2u","E2u","B2g"}
                             :std::vector<std::string>{"A1g","E1u","E1u","E2g","E2g","B1u"});
        if (result.point_group!="D"+std::to_string(n)+"h" || result.orbitals_labelled!=expected.size()) {
            std::cerr << "even Dnh representation was not fully classified\n";return EXIT_FAILURE;
        }
        for (std::size_t i=0;i<expected.size();++i)
            if (!all_label(wf,i,i+1,expected[i])) return EXIT_FAILURE;
        for (const auto& assignment:wf.derived_orbital_symmetry_assignments) {
            if (assignment.point_group!=result.point_group || !assignment.axes_available ||
                assignment.orbital_indices.empty() || assignment.subspace_retention<.99999 ||
                assignment.maximum_character_error>1e-5) {
                std::cerr << "Dnh label lost its subspace/frame evidence\n";return EXIT_FAILURE;
            }
        }
    }
    {
        auto wf=even_ring(6,true,false);
        wf.orbitals.erase(wf.orbitals.begin()+2); // Incomplete E1g subspace.
        (void)cov::derive_orbital_symmetry(wf);
        if (!wf.orbitals[1].symmetry.empty()) {
            std::cerr << "non-closed E subspace received a strict 1D label\n";return EXIT_FAILURE;
        }
    }
    {
        auto wf=d5h_ring();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="D5h" || r.orbitals_labelled!=5u ||
            !all_label(wf,0,1,"A1'") || !all_label(wf,1,3,"E1'") ||
            !all_label(wf,3,5,"E2'")) {
            std::cerr << "D5h derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=td_ligand_sigma();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="Td" || r.orbitals_labelled!=4u ||
            !all_label(wf,0,1,"A1") || !all_label(wf,1,4,"T2")) {
            std::cerr << "Td derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=oh_central_p();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="Oh" || r.orbitals_labelled!=3u ||
            !all_label(wf,0,3,"T1u")) {
            std::cerr << "Oh AO-angular transform/derived-irrep regression"
                      << " point_group=" << r.point_group
                      << " groups_examined=" << r.groups_examined
                      << " groups_labelled=" << r.groups_labelled
                      << " orbitals_labelled=" << r.orbitals_labelled
                      << " worst_retention=" << r.worst_subspace_retention << '\n';
            for (std::size_t i=0;i<wf.orbitals.size();++i) {
                std::cerr << "  MO " << i << " label='" << wf.orbitals[i].symmetry
                          << "' provenance=" << static_cast<int>(wf.orbitals[i].symmetry_provenance)
                          << '\n';
            }
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=d3h_ring();
        wf.orbitals[0].symmetry="producer";
        wf.orbitals[0].symmetry_provenance=cov::DataProvenance::Producer;
        (void)cov::derive_orbital_symmetry(wf);
        if (wf.orbitals[0].symmetry!="producer" ||
            wf.orbitals[0].symmetry_provenance!=cov::DataProvenance::Producer) {
            std::cerr << "producer MO symmetry was overwritten\n";
            return EXIT_FAILURE;
        }
    }
    std::cout << "orbital symmetry smoke test passed\n";
    return EXIT_SUCCESS;
}
