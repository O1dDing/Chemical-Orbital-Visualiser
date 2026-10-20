#include "cov/point_group_irreps.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <numeric>
#include <random>

namespace {
using Mat=std::array<double,9>;
using Vec=std::array<double,3>;
constexpr Mat e{1,0,0,0,1,0,0,0,1};
constexpr double pi=std::numbers::pi;
bool okay=true;
void require(bool condition,const std::string& detail){if(!condition){okay=false;std::cerr<<detail<<'\n';}}
Mat product(const Mat& a,const Mat& b){Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
Mat transpose(const Mat& a){return {a[0],a[3],a[6],a[1],a[4],a[7],a[2],a[5],a[8]};}
double error(const Mat& a,const Mat& b){double d=0;for(int i=0;i<9;++i)d=std::max(d,std::abs(a[i]-b[i]));return d;}
double determinant(const Mat& a){return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);}
Mat neg(Mat a){for(auto& x:a)x=-x;return a;}
Mat rz(double theta){return {std::cos(theta),-std::sin(theta),0,std::sin(theta),std::cos(theta),0,0,0,1};}
constexpr Mat h{1,0,0,0,1,0,0,0,-1},v{1,0,0,0,-1,0,0,0,1},c2x{1,0,0,0,-1,0,0,0,-1};
cov::MolecularSymmetry closure(const std::string& name,const std::vector<Mat>& generators){
    cov::MolecularSymmetry s;s.point_group=name;s.operations.emplace_back();
    for(std::size_t i=0;i<s.operations.size();++i)for(const auto& generator:generators){const auto a=product(s.operations[i].matrix,generator);
        if(std::none_of(s.operations.begin(),s.operations.end(),[&](const auto& op){return error(op.matrix,a)<1e-8;})){cov::SymmetryOperation op;op.matrix=a;s.operations.push_back(op);}
        if(s.operations.size()>2000){require(false,"test generator closure exceeded expected finite bound");return s;}
    }return s;
}
const cov::PointGroupIrrepRow* row(const cov::PointGroupIrrepTable& t,const std::string& label){for(const auto& r:t.rows)if(r.label==label)return &r;return nullptr;}

// Every group is checked against an independent physical vector character,
// a deliberately mixed/repeated character and the full regular character.
cov::PointGroupIrrepTable check(const cov::MolecularSymmetry& s,const cov::PointGroupIrrepOptions& options={}){
    auto t=cov::finite_point_group_irreps(s,options);require(t.valid,s.point_group+": "+t.reason);
    if(!t.valid)return t;
    require(s.operations.size()==cov::finite_point_group_order(s.point_group),s.point_group+": parameterized group order");
    std::vector<double> vector_char,axial_char,quadratic_char,regular(s.operations.size()),mixed(s.operations.size());std::size_t dimension=0;
    regular[t.identity_operation]=double(s.operations.size());
    for(const auto& op:s.operations){const double trace=op.matrix[0]+op.matrix[4]+op.matrix[8];const auto square=product(op.matrix,op.matrix);
        vector_char.push_back(trace);axial_char.push_back(determinant(op.matrix)*trace);
        quadratic_char.push_back(.5*(trace*trace+square[0]+square[4]+square[8])-1);}
    const auto p=cov::decompose_point_group_characters(t,vector_char);require(p.valid&&p.dimension==3,s.point_group+": vector decomposition: "+p.reason);
    const auto axial=cov::decompose_point_group_characters(t,axial_char);require(axial.valid&&axial.dimension==3,s.point_group+": axial vector decomposition: "+axial.reason);
    const auto quadratic=cov::decompose_point_group_characters(t,quadratic_char);require(quadratic.valid&&quadratic.dimension==5,s.point_group+": traceless quadratic decomposition: "+quadratic.reason);
    const auto r=cov::decompose_point_group_characters(t,regular);require(r.valid&&r.dimension==s.operations.size(),s.point_group+": regular decomposition: "+r.reason);
    for(std::size_t j=0;j<t.rows.size();++j){const auto count=1+j%3;dimension+=count*t.rows[j].dimension;for(std::size_t k=0;k<mixed.size();++k)mixed[k]+=double(count)*t.rows[j].characters[k];}
    const auto m=cov::decompose_point_group_characters(t,mixed);require(m.valid&&m.dimension==dimension,s.point_group+": mixed repeated decomposition: "+m.reason);
    if(m.valid)for(std::size_t j=0;j<t.rows.size();++j)require(m.multiplicities[j]==1+j%3,s.point_group+": mixed copy count lost");
    mixed[(t.identity_operation+1)%mixed.size()]+=.13;
    require(!cov::decompose_point_group_characters(t,mixed).valid,s.point_group+": inconsistent characters accepted");
    return t;
}

void transformed(const cov::MolecularSymmetry& s,const cov::PointGroupIrrepOptions& axes){
    // Generic rigid rotation; explicit axes rotate with the actual operations.
    const double a=.713,b=.437;const Mat ry{std::cos(b),0,std::sin(b),0,1,0,-std::sin(b),0,std::cos(b)};
    const Mat q=product(rz(a),ry);cov::MolecularSymmetry rotated=s;
    for(auto& op:rotated.operations)op.matrix=product(product(q,op.matrix),transpose(q));
    auto options=axes;const auto apply=[&](const Vec& x){Vec y{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)y[i]+=q[3*i+j]*x[j];return y;};
    if(options.principal_axis)options.principal_axis=apply(*options.principal_axis);
    if(options.reference_axis)options.reference_axis=apply(*options.reference_axis);
    const auto original=check(s,axes),t=check(rotated,options);
    if(original.valid&&t.valid){require(original.rows.size()==t.rows.size(),s.point_group+": rigid-rotation row count");for(std::size_t j=0;j<t.rows.size();++j){require(original.rows[j].label==t.rows[j].label,s.point_group+": rigid-rotation label");for(std::size_t k=0;k<s.operations.size();++k)require(std::abs(original.rows[j].characters[k]-t.rows[j].characters[k])<1e-7,s.point_group+": rigid-rotation character");}}
    std::vector<std::size_t> permutation(s.operations.size());std::iota(permutation.begin(),permutation.end(),0);std::mt19937 random(1701);std::shuffle(permutation.begin(),permutation.end(),random);
    cov::MolecularSymmetry reordered=rotated;for(std::size_t i=0;i<permutation.size();++i)reordered.operations[i]=rotated.operations[permutation[i]];
    const auto u=check(reordered,options);if(t.valid&&u.valid)for(std::size_t j=0;j<t.rows.size();++j)for(std::size_t k=0;k<permutation.size();++k)require(std::abs(u.rows[j].characters[k]-t.rows[j].characters[permutation[k]])<1e-7,s.point_group+": operation permutation");
}

cov::MolecularSymmetry cubic(const std::string& name){
    cov::MolecularSymmetry s;s.point_group=name;std::array<int,3> permutation{0,1,2};
    do{for(int sx:{-1,1})for(int sy:{-1,1})for(int sz:{-1,1}){Mat a{};a[permutation[0]]=sx;a[3+permutation[1]]=sy;a[6+permutation[2]]=sz;
        const bool tetrahedral=sx*sy*sz==1,proper=determinant(a)>0;
        if((name=="Oh")||(name=="O"&&proper)||(name=="Td"&&tetrahedral)||(name=="T"&&tetrahedral&&proper)||(name=="Th"&&((tetrahedral&&proper)||(!tetrahedral&&!proper)))){cov::SymmetryOperation op;op.matrix=a;s.operations.push_back(op);}
    }}while(std::next_permutation(permutation.begin(),permutation.end()));return s;
}
double dot(const Vec& a,const Vec& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec cross(const Vec& a,const Vec& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
Vec unit(Vec a){const double n=std::sqrt(dot(a,a));for(auto& x:a)x/=n;return a;}
Mat pair_frame(Vec a,Vec b){a=unit(a);const double p=dot(a,b);for(int i=0;i<3;++i)b[i]-=p*a[i];b=unit(b);const auto c=cross(a,b);return {a[0],b[0],c[0],a[1],b[1],c[1],a[2],b[2],c[2]};}
cov::MolecularSymmetry icosahedral(bool inversion){
    const double phi=(1+std::sqrt(5.0))/2;std::vector<Vec> vertices;
    for(double a:{-1.0,1.0})for(double b:{-phi,phi}){vertices.push_back({0,a,b});vertices.push_back({a,b,0});vertices.push_back({b,0,a});}
    const auto frame=pair_frame(vertices[0],vertices[1]);const double target=dot(vertices[0],vertices[1]);cov::MolecularSymmetry s;s.point_group=inversion?"Ih":"I";
    for(const auto& a:vertices)for(const auto& b:vertices){if(std::abs(dot(a,b)-target)>1e-8)continue;const auto matrix=product(pair_frame(a,b),transpose(frame));
        bool maps=true;for(const auto& x:vertices){Vec y{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)y[i]+=matrix[3*i+j]*x[j];if(std::none_of(vertices.begin(),vertices.end(),[&](const auto& z){Vec d{};for(int k=0;k<3;++k)d[k]=y[k]-z[k];return dot(d,d)<1e-10;})){maps=false;break;}}
        if(maps){cov::SymmetryOperation op;op.matrix=matrix;s.operations.push_back(op);if(inversion){op.matrix=neg(matrix);s.operations.push_back(op);}}
    }return s;
}
void expect(const cov::MolecularSymmetry& s,const cov::PointGroupIrrepOptions& axes,
            const std::vector<std::pair<std::string,std::size_t>>& p,
            const std::vector<std::pair<std::string,std::size_t>>& axial,
            const std::vector<std::pair<std::string,std::size_t>>& d){
    const auto table=cov::finite_point_group_irreps(s,axes);require(table.valid,s.point_group+": known physical characters table");if(!table.valid)return;
    for(int type=0;type<3;++type){std::vector<double> characters;
        for(const auto& op:s.operations){const auto& r=op.matrix;const double trace=r[0]+r[4]+r[8];const auto r2=product(r,r);
            characters.push_back(type==0?trace:type==1?determinant(r)*trace:.5*(trace*trace+r2[0]+r2[4]+r2[8])-1);}
        const auto decomposition=cov::decompose_point_group_characters(table,characters);require(decomposition.valid,s.point_group+": physical character decomposition");if(!decomposition.valid)continue;
        const auto& expected=type==0?p:type==1?axial:d;
        for(std::size_t j=0;j<table.rows.size();++j){std::size_t count=0;for(const auto& item:expected)if(item.first==table.rows[j].label)count=item.second;
            require(decomposition.multiplicities[j]==count,s.point_group+": known "+(type==0?"p":type==1?"R":"d")+" multiplicity for "+table.rows[j].label);}
    }
}
} // namespace

int main(){
    cov::PointGroupIrrepOptions axes;axes.principal_axis=Vec{0,0,1};axes.reference_axis=Vec{1,0,0};
    for(const auto n:{2,3,4,5,6,7,8,11,12,17,65}){
        const auto r=rz(2*pi/n),s=product(rz(pi/n),h);const auto number=std::to_string(n);
        for(const auto& g:std::vector<std::pair<std::string,std::vector<Mat>>>{{"C"+number,{r}},{"C"+number+"v",{r,v}},{"C"+number+"h",{r,h}},{"D"+number,{r,c2x}},{"D"+number+"h",{r,c2x,h}},{"D"+number+"d",{s,c2x}},{"S"+number,{product(r,h)}}}){
            const auto group=closure(g.first,g.second);transformed(group,axes);check(group);
        }
    }
    for(const auto& name:{"C1","Cs","Ci"})transformed(closure(name,std::string(name)=="C1"?std::vector<Mat>{}:std::string(name)=="Cs"?std::vector<Mat>{h}:std::vector<Mat>{neg(e)}),{});
    for(const auto& name:{"T","Td","Th","O","Oh"})transformed(cubic(name),{});
    transformed(icosahedral(false),{});transformed(icosahedral(true),{});
    expect(cubic("T"),{},{{"T",1}},{{"T",1}},{{"E",1},{"T",1}});
    expect(cubic("Th"),{},{{"Tu",1}},{{"Tg",1}},{{"Eg",1},{"Tg",1}});
    expect(cubic("Td"),{},{{"T2",1}},{{"T1",1}},{{"E",1},{"T2",1}});
    expect(cubic("O"),{},{{"T1",1}},{{"T1",1}},{{"E",1},{"T2",1}});
    expect(cubic("Oh"),{},{{"T1u",1}},{{"T1g",1}},{{"Eg",1},{"T2g",1}});
    expect(icosahedral(false),{},{{"T1",1}},{{"T1",1}},{{"H",1}});
    expect(icosahedral(true),{},{{"T1u",1}},{{"T1g",1}},{{"Hg",1}});
    expect(closure("C2v",{rz(pi),v}),axes,{{"A1",1},{"B1",1},{"B2",1}},{{"A2",1},{"B1",1},{"B2",1}},{{"A1",2},{"A2",1},{"B1",1},{"B2",1}});
    expect(closure("D2h",{rz(pi),c2x,h}),axes,{{"B1u",1},{"B2u",1},{"B3u",1}},{{"B1g",1},{"B2g",1},{"B3g",1}},{{"Ag",2},{"B1g",1},{"B2g",1},{"B3g",1}});
    expect(closure("D3h",{rz(2*pi/3),c2x,h}),axes,{{"A2''",1},{"E'",1}},{{"A2'",1},{"E''",1}},{{"A1'",1},{"E'",1},{"E''",1}});
    expect(closure("D4h",{rz(pi/2),c2x,h}),axes,{{"A2u",1},{"Eu",1}},{{"A2g",1},{"Eg",1}},{{"A1g",1},{"B1g",1},{"B2g",1},{"Eg",1}});
    expect(closure("D3d",{product(rz(pi/3),h),c2x}),axes,{{"A2u",1},{"Eu",1}},{{"A2g",1},{"Eg",1}},{{"A1g",1},{"Eg",2}});
    expect(closure("D5d",{product(rz(pi/5),h),c2x}),axes,{{"A2u",1},{"E1u",1}},{{"A2g",1},{"E1g",1}},{{"A1g",1},{"E1g",1},{"E2g",1}});
    expect(closure("D2d",{product(rz(pi/2),h),c2x}),axes,{{"B2",1},{"E",1}},{{"A2",1},{"E",1}},{{"A1",1},{"B1",1},{"B2",1},{"E",1}});
    expect(closure("D4d",{product(rz(pi/4),h),c2x}),axes,{{"B2",1},{"E1",1}},{{"A2",1},{"E3",1}},{{"A1",1},{"E2",1},{"E3",1}});
    expect(closure("S4",{product(rz(pi/2),h)}),axes,{{"B",1},{"E",1}},{{"A",1},{"E",1}},{{"A",1},{"B",2},{"E",1}});
    expect(closure("S6",{product(rz(pi/3),h)}),axes,{{"Au",1},{"Eu",1}},{{"Ag",1},{"Eg",1}},{{"Ag",1},{"Eg",2}});
    const auto c2v=check(closure("C2v",{rz(pi),v}),axes);const auto b1=row(c2v,"B1"),b2=row(c2v,"B2");
    require(b1&&b2,"C2v complete B labels");if(b1&&b2){const auto g=closure("C2v",{rz(pi),v});for(std::size_t k=0;k<g.operations.size();++k){require(std::abs(b1->characters[k]-g.operations[k].matrix[0])<1e-7,"C2v B1 transforms as x");require(std::abs(b2->characters[k]-g.operations[k].matrix[4])<1e-7,"C2v B2 transforms as y");}}
    const auto c3=check(closure("C3",{rz(2*pi/3)}));require(row(c3,"E")&&row(c3,"E")->character_norm==2,"C3 real paired E norm");
    const auto d4d=check(closure("D4d",{product(rz(pi/4),h),c2x}),axes);require(row(d4d,"B1")&&row(d4d,"E3"),"D4d complete B/E rows");
    const auto ih=check(icosahedral(true));require(row(ih,"Hg")&&row(ih,"Hg")->dimension==5&&row(ih,"Gu")&&row(ih,"Gu")->dimension==4,"Ih high-dimensional rows");
    auto broken=closure("C4v",{rz(pi/2),v});broken.operations.pop_back();require(!cov::finite_point_group_irreps(broken).valid,"missing operation accepted");
    broken=closure("C4v",{rz(pi/2),v});broken.operations[1].matrix[0]+=.01;require(!cov::finite_point_group_irreps(broken).valid,"nonorthogonal operation accepted");
    broken=closure("C4v",{rz(pi/2),v});broken.operations[1]=broken.operations[0];require(!cov::finite_point_group_irreps(broken).valid,"duplicate operation accepted");
    broken=closure("C4v",{rz(pi/2),v});broken.point_group="D4";require(!cov::finite_point_group_irreps(broken).valid,"proper/improper groups confused");
    require(cov::finite_point_group_order("C127v")==254&&cov::finite_point_group_order("D127d")==508,"arbitrary n order API");
    for(const auto& name:{"Dinfh","Cinfv","Kh","C0","D1","C3q","S4v","x"})require(cov::finite_point_group_order(name)==0,std::string(name)+": unsupported finite order");
    if(okay)std::cout<<"complete finite point-group characters, axes, paired norms, decompositions and invariance passed\n";
    return okay?EXIT_SUCCESS:EXIT_FAILURE;
}
