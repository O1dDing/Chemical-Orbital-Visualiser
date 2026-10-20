#include "cov/symmetry.hpp"
#include "cov/point_group_irreps.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

cov::Atom atom(const char* symbol, int z, double x, double y, double zc) {
    return {symbol, z, x, y, zc};
}

bool expect_group(const cov::Wavefunction& wf,
                  const std::string& expected,
                  const char* label) {
    const auto result = cov::analyse_molecular_symmetry(wf);
    if (result.point_group != expected) {
        std::cerr << label << ": expected " << expected
                  << ", got " << result.point_group
                  << " (ops=" << result.operations.size()
                  << ", tol=" << result.tolerance_bohr << ")\n";
        return false;
    }
    if (result.operations.empty()) {
        std::cerr << label << ": no validated symmetry operations retained\n";
        return false;
    }
    if(!result.linear){
        const auto table=cov::finite_point_group_irreps(result);
        if(!table.valid||result.operations.size()!=cov::finite_point_group_order(expected)){
            std::cerr<<label<<": finite geometry operations are incomplete or inconsistent: "<<table.reason<<'\n';return false;
        }
        for(const auto& op:result.operations)if(op.max_mapping_error_bohr>result.tolerance_bohr){
            std::cerr<<label<<": group completion exceeded the original geometry tolerance\n";return false;
        }
    }
    return true;
}

cov::Wavefunction h2() {
    cov::Wavefunction wf;
    wf.atoms = {
        atom("H", 1, 0.0, 0.0, -0.7),
        atom("H", 1, 0.0, 0.0,  0.7),
    };
    return wf;
}

cov::Wavefunction trigonal_planar() {
    cov::Wavefunction wf;
    const double r = 1.4;
    for (int i = 0; i < 3; ++i) {
        const double a = 2.0 * pi * static_cast<double>(i) / 3.0;
        wf.atoms.push_back(atom("H", 1, r * std::cos(a), r * std::sin(a), 0.0));
    }
    return wf;
}

cov::Wavefunction pentagonal_planar() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Fe", 26, 0.0, 0.0, 0.0));
    const double r = 2.0;
    for (int i = 0; i < 5; ++i) {
        const double a = 2.0 * pi * static_cast<double>(i) / 5.0;
        wf.atoms.push_back(atom("C", 6, r * std::cos(a), r * std::sin(a), 0.0));
    }
    return wf;
}

cov::Wavefunction tetrahedral() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Zn", 30, 0.0, 0.0, 0.0));
    wf.atoms.push_back(atom("Cl", 17,  1.0,  1.0,  1.0));
    wf.atoms.push_back(atom("Cl", 17,  1.0, -1.0, -1.0));
    wf.atoms.push_back(atom("Cl", 17, -1.0,  1.0, -1.0));
    wf.atoms.push_back(atom("Cl", 17, -1.0, -1.0,  1.0));
    return wf;
}

cov::Wavefunction octahedral() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Ti", 22, 0.0, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9,  1.8, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, -1.8, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0,  1.8, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0, -1.8, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0, 0.0,  1.8));
    wf.atoms.push_back(atom("F", 9, 0.0, 0.0, -1.8));
    return wf;
}

cov::Wavefunction dihedral_rounded(bool staggered) {
    cov::Wavefunction wf;
    wf.atoms={atom("C",6,0,0,-.7654),atom("C",6,0,0,.7653)};
    const auto rounded=[](double x){return std::round(x*10000)/10000;};
    for(int side:{-1,1})for(int j=0;j<3;++j){
        const double a=2*pi*j/3+(side==1&&staggered?pi/3:0);
        wf.atoms.push_back(atom("H",1,rounded(1.0185*std::cos(a)),rounded(1.0185*std::sin(a)),side*1.1636));
    }
    return wf;
}

cov::Wavefunction perpendicular_terminal_planes() {
    cov::Wavefunction wf;
    wf.atoms={atom("C",6,0,0,0),atom("C",6,0,0,1.3024),atom("C",6,0,0,-1.3024)};
    for(int sign:{-1,1}){
        wf.atoms.push_back(atom("H",1,sign*.6557,sign*.6557,1.8669));
        wf.atoms.push_back(atom("H",1,sign*.6557,-sign*.6557,-1.8669));
    }
    return wf;
}

cov::Wavefunction tetrahedral_water_planes() {
    cov::Wavefunction wf;wf.atoms.push_back(atom("Mg",12,0,0,0));
    for(int coordinate=0;coordinate<3;++coordinate)for(int sign:{-1,1}){
        std::array<double,3> p{};p[coordinate]=sign*2.09;wf.atoms.push_back(atom("O",8,p[0],p[1],p[2]));
        for(int hydrogen_sign:{-1,1}){
            p[coordinate]=sign*2.6988;p[(coordinate+1)%3]=hydrogen_sign*.7934;
            wf.atoms.push_back(atom("H",1,p[0],p[1],p[2]));
        }
    }return wf;
}

cov::Wavefunction rounded_polygon(int n) {
    cov::Wavefunction wf;const auto rounded=[](double x){return std::round(x*10000)/10000;};
    for(int j=0;j<n;++j){const double a=2*pi*j/n;wf.atoms.push_back(atom("C",6,rounded(1.6115*std::cos(a)),rounded(1.6115*std::sin(a)),0));}
    return wf;
}

cov::Wavefunction rounded_icosahedron() {
    cov::Wavefunction wf;const double phi=(1+std::sqrt(5.0))/2;
    for(double a:{-1.0,1.0})for(double b:{-phi,phi}){
        wf.atoms.push_back(atom("B",5,0,a,std::round(b*10000)/10000));
        wf.atoms.push_back(atom("B",5,a,std::round(b*10000)/10000,0));
        wf.atoms.push_back(atom("B",5,std::round(b*10000)/10000,0,a));
    }return wf;
}

cov::Wavefunction rotated_reordered(cov::Wavefunction w) {
    const double a=.371,b=.537,c=std::cos(a),s=std::sin(a),d=std::cos(b),t=std::sin(b);
    for(auto& p:w.atoms){const double x=d*p.x+t*p.z,y=p.y,z=-t*p.x+d*p.z;
        p.x=c*x-s*y+3.1;p.y=s*x+c*y-1.2;p.z=z+.9;}
    std::reverse(w.atoms.begin(),w.atoms.end());return w;
}

} // namespace

int main() {
    if (!expect_group(h2(), "Dinfh", "H2")) return EXIT_FAILURE;
    if (!expect_group(trigonal_planar(), "D3h", "D3h planar")) return EXIT_FAILURE;
    if (!expect_group(pentagonal_planar(), "D5h", "D5h planar")) return EXIT_FAILURE;
    if (!expect_group(tetrahedral(), "Td", "Td tetrahedron")) return EXIT_FAILURE;
    if (!expect_group(octahedral(), "Oh", "Oh octahedron")) return EXIT_FAILURE;
    const std::vector<std::pair<cov::Wavefunction,std::string>> regressions{
        {perpendicular_terminal_planes(),"D2d"},{dihedral_rounded(true),"D3d"},
        {dihedral_rounded(false),"D3h"},{tetrahedral_water_planes(),"Th"},
        {rounded_polygon(7),"D7h"},{rounded_icosahedron(),"Ih"}};
    for(const auto& [w,group]:regressions){
        if(!expect_group(w,group,"consistent finite geometry"))return EXIT_FAILURE;
        if(!expect_group(rotated_reordered(w),group,"3D rotated/reordered finite geometry"))return EXIT_FAILURE;
    }

    cov::Wavefunction producer = tetrahedral();
    producer.point_group_detected = "PRODUCER-GROUP";
    producer.point_group_used = "PRODUCER-GROUP";
    producer.point_group_provenance = cov::DataProvenance::Producer;
    cov::derive_point_group_from_geometry(producer);
    if (producer.point_group_detected != "PRODUCER-GROUP" ||
        producer.point_group_provenance != cov::DataProvenance::Producer) {
        std::cerr << "producer symmetry provenance was overwritten\n";
        return EXIT_FAILURE;
    }

    cov::Wavefunction derived = tetrahedral();
    cov::derive_point_group_from_geometry(derived);
    if (derived.point_group_detected != "Td" ||
        derived.point_group_provenance != cov::DataProvenance::Derived ||
        derived.point_group_detected_provenance != cov::DataProvenance::Derived ||
        !derived.point_group_used.empty() ||
        derived.point_group_used_provenance != cov::DataProvenance::Unavailable) {
        std::cerr << "geometry-derived symmetry was not stored with derived provenance\n";
        return EXIT_FAILURE;
    }
    auto used_only=tetrahedral();
    used_only.point_group_used="C1";
    used_only.point_group_used_provenance=cov::DataProvenance::Producer;
    used_only.point_group_provenance=cov::DataProvenance::Producer;
    cov::derive_point_group_from_geometry(used_only);
    if (used_only.point_group_detected!="Td" || used_only.point_group_used!="C1" ||
        used_only.point_group_detected_provenance!=cov::DataProvenance::Derived ||
        used_only.point_group_used_provenance!=cov::DataProvenance::Producer) {
        std::cerr << "geometry-derived full group was blocked by producer-used subgroup\n";
        return EXIT_FAILURE;
    }

    std::cout << "symmetry smoke test passed\n";
    return EXIT_SUCCESS;
}
