#include "cov/molecular_point_group_frame.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>

namespace {
cov::Wavefunction read(const std::string& file){
    std::ifstream input(std::string(COV_SOURCE_DIR)+"/tests/fixtures/point_groups/"+file+".xyz");
    std::size_t n=0;input>>n;std::string line;std::getline(input,line);std::getline(input,line);
    const std::map<std::string,unsigned> elements{{"H",1},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},
        {"Ne",10},{"Mg",12},{"Si",14},{"P",15},{"S",16},{"Ar",18},{"Cr",24},{"Fe",26},{"Co",27},{"U",92}};
    cov::Wavefunction w;for(std::size_t i=0;i<n;++i){cov::Atom atom;input>>atom.symbol>>atom.x>>atom.y>>atom.z;
        if(!input||!elements.contains(atom.symbol))throw std::runtime_error("invalid public fixture "+file);
        atom.atomic_number=elements.at(atom.symbol);atom.x*=1.889726124626;atom.y*=1.889726124626;atom.z*=1.889726124626;w.atoms.push_back(atom);}
    if(!n)throw std::runtime_error("missing public fixture "+file);return w;
}
bool check(const cov::Wavefunction& w,const std::string& name,const std::string& expected){
    const auto geometry=cov::analyse_molecular_symmetry(w);const auto full=cov::complete_molecular_point_group(w,geometry);
    const auto table=cov::molecular_point_group_irreps(w,full);
    std::cout<<name<<": group="<<geometry.point_group<<" operations="<<full.operations.size()<<" expected="<<expected<<" table="<<table.valid<<" "<<table.reason<<'\n';
    return geometry.point_group==expected&&full.operations.size()==cov::finite_point_group_order(expected)&&table.valid;
}
}
int main(){try{
    bool okay=true;
    const std::pair<const char*,const char*> cases[]={{"allene-D2d","D2d"},{"ethane-staggered-D3d","D3d"},{"ethane-eclipsed-D3h","D3h"},{"neopentane-T","T"},{"neopentane-Td","Td"},{"mg-aqua-Th","Th"},{"octamethyl-POSS-O","O"},{"tropylium-D7h","D7h"},{"dodecaborate-Ih","Ih"},{"cyclopentadienyl-D5h","D5h"},{"corannulene-C5v","C5v"}};
    for(const auto& [file,group]:cases){auto w=read(file);okay=check(w,file,group)&&okay;
        const double c=std::cos(.371),s=std::sin(.371);for(auto& a:w.atoms){const double x=a.x,y=a.y;a.x=c*x-s*y+3.1;a.y=s*x+c*y-1.2;a.z+=.9;}
        std::reverse(w.atoms.begin(),w.atoms.end());okay=check(w,std::string(file)+" rotated/reordered",group)&&okay;}
    std::ifstream representatives(std::string(COV_SOURCE_DIR)+"/tests/fixtures/point_groups/representatives/cases.tsv");
    std::string file,group;
    std::size_t representative_count=0;
    while(representatives>>file>>group) {
        auto w=read("representatives/"+file);okay=check(w,file,group)&&okay;
        const double c=std::cos(.713),s=std::sin(.713);
        for(auto& a:w.atoms){const double y=a.y,z=a.z;a.x+=2.4;a.y=c*y-s*z-.9;a.z=s*y+c*z+1.3;}
        std::reverse(w.atoms.begin(),w.atoms.end());okay=check(w,file+" rotated/reordered",group)&&okay;
        ++representative_count;
    }
    if(representative_count<30)throw std::runtime_error("incomplete extended geometry fixture inventory");
    // These stored atoms exposed loss of valid operations when the initial
    // frame was perturbed. Accept every operation only at the original tolerance.
    for(const auto& [fixture,expected,atom]:std::vector<std::tuple<std::string,std::string,std::size_t>>{
        {"ferrocene-staggered-D5d","D5d",8},{"corannulene-C5v","C5v",24},
        {"fullerene-C140-I","I",33}}) {
        const auto original=read("representatives/"+fixture);
        const auto tolerance=cov::analyse_molecular_symmetry(original).tolerance_bohr;
        for(double scale:{.25,.5,.75})for(int axis=0;axis<3;++axis) {
            auto moved=original;auto& a=moved.atoms.at(atom);
            if(axis==0)a.x+=scale*tolerance;if(axis==1)a.y+=scale*tolerance;if(axis==2)a.z+=scale*tolerance;
            const auto analysis=cov::analyse_molecular_symmetry(moved);
            okay=analysis.point_group==expected&&analysis.operations.size()==cov::finite_point_group_order(expected)&&
                analysis.max_mapping_error_bohr<=analysis.tolerance_bohr&&okay;
        }
    }
    for(const auto& [file,group]:std::vector<std::pair<std::string,std::string>>{{"neopentane-T","Td"},{"mg-aqua-Th","Oh"},{"octamethyl-POSS-O","Oh"}}){auto w=read(file);
        w.atoms.erase(std::remove_if(w.atoms.begin(),w.atoms.end(),[](const auto& a){return a.atomic_number==1;}),w.atoms.end());
        okay=check(w,file+" explicit non-H analysis",group)&&okay;}
    auto perturbed=read("ethane-staggered-D3d");perturbed.atoms[0].x+=1e-6;okay=check(perturbed,"within recorded geometry tolerance","D3d")&&okay;
    perturbed.atoms[0].x+=.03;if(cov::analyse_molecular_symmetry(perturbed).point_group=="D3d"){std::cerr<<"above-tolerance perturbation silently symmetrized\n";okay=false;}
    return okay?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
