#pragma once
#include "cov/symmetry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace cov::orbital_characters {
using Vector=std::array<double,3>;
using Matrix=std::array<double,9>;
inline double dot(const Vector& a,const Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Matrix multiply(const Matrix& a,const Matrix& b){Matrix c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
inline double matrix_error(const Matrix& a,const Matrix& b){double e=0;for(int i=0;i<9;++i)e=std::max(e,std::abs(a[i]-b[i]));return e;}
// Geometry may supply generators. Complete the same validated atom action,
// without inventing an input-axis binding or using operation kind metadata.
inline std::vector<SymmetryOperation> d2h_operations(const Wavefunction& w,const std::vector<SymmetryOperation>& generators){
    SymmetryOperation e;e.atom_permutation.resize(w.atoms.size());std::iota(e.atom_permutation.begin(),e.atom_permutation.end(),0);
    std::vector<SymmetryOperation> ops{e};
    for(std::size_t i=0;i<ops.size();++i)for(const auto& g:generators){
        if(g.atom_permutation.size()!=w.atoms.size())return {};
        SymmetryOperation p;p.matrix=multiply(ops[i].matrix,g.matrix);p.atom_permutation.resize(w.atoms.size());
        for(std::size_t a=0;a<w.atoms.size();++a){if(g.atom_permutation[a]>=w.atoms.size())return {};p.atom_permutation[a]=ops[i].atom_permutation[g.atom_permutation[a]];}
        if(std::any_of(ops.begin(),ops.end(),[&](const auto& q){return p.atom_permutation==q.atom_permutation&&matrix_error(p.matrix,q.matrix)<2e-5;}))continue;
        if(ops.size()>=8)return {};ops.push_back(std::move(p));
    }
    return ops.size()==8?ops:std::vector<SymmetryOperation>{};
}
struct D2hAxes {std::array<Vector,3> xyz{};bool resolved=false;std::string detail;};
inline D2hAxes d2h_axes(const Wavefunction& w,const std::vector<SymmetryOperation>& ops,const Vector& centre,double tolerance){
    D2hAxes out;std::vector<Vector> axes;
    for(const auto& op:ops){const auto& m=op.matrix;const double determinant=m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);
        if(std::abs(determinant-1)>2e-5||std::abs(m[0]+m[4]+m[8]+1)>2e-5)continue;
        Vector axis{};double best=0;for(int c=0;c<3;++c){Vector v{};for(int r=0;r<3;++r)v[r]=m[3*r+c]+(r==c?1.:0.);const double n=dot(v,v);if(n>best){best=n;axis=v;}}
        if(best<=1e-20)continue;for(auto& x:axis)x/=std::sqrt(best);
        if(std::none_of(axes.begin(),axes.end(),[&](const auto& v){return std::abs(dot(v,axis))>1-2e-5;}))axes.push_back(axis);
    }
    out.detail="D2h axes unresolved; only axis-independent Ag/Au names are available";
    if(axes.size()!=3)return out;for(int a=0;a<3;++a)for(int b=0;b<a;++b)if(std::abs(dot(axes[a],axes[b]))>2e-5)return out;
    std::array<double,3> moment{};for(int k=0;k<3;++k)for(const auto& a:w.atoms){const Vector r{a.x-centre[0],a.y-centre[1],a.z-centre[2]};const double projection=dot(r,axes[k]);moment[k]+=projection*projection;}
    std::array<int,3> order{0,1,2};std::sort(order.begin(),order.end(),[&](int a,int b){return moment[a]<moment[b];});
    const double margin=std::max(tolerance*tolerance*double(w.atoms.size()),1e-10*std::max(1.,moment[order[2]]));
    if(moment[order[1]]-moment[order[0]]<=margin||moment[order[2]]-moment[order[1]]<=margin)return out;
    // A recorded geometrical convention, covariant under rigid rotation and
    // atom permutation. Degenerate moments never get arbitrary B1/B2/B3 names.
    out.xyz={axes[order[2]],axes[order[0]],axes[order[1]]};out.resolved=true;
    std::ostringstream detail;detail<<"D2h: x/y/z are C2 axes with maximal/minimal/intermediate unweighted nuclear second moments; B1u=z, B2u=y, B3u=x; input-frame axes=";
    for(const auto& a:out.xyz)detail<<'('<<a[0]<<','<<a[1]<<','<<a[2]<<')';out.detail=detail.str();return out;
}
inline std::string match_d2h(const D2hAxes& axes,const std::vector<SymmetryOperation>& ops,const std::vector<double>& chars,std::size_t dimension,double tolerance){
    if(ops.size()!=8||chars.size()!=ops.size()||dimension!=1)return {};
    struct Row{const char* name;std::array<int,3> powers;};
    const Row rows[]={{"Ag",{0,0,0}},{"B1g",{1,1,0}},{"B2g",{1,0,1}},{"B3g",{0,1,1}},{"Au",{1,1,1}},{"B1u",{0,0,1}},{"B2u",{0,1,0}},{"B3u",{1,0,0}}};
    for(const auto& row:rows){if(!axes.resolved&&std::string(row.name)!="Ag"&&std::string(row.name)!="Au")continue;bool same=true;
        for(std::size_t i=0;i<ops.size();++i){double expected=1;
            if(!axes.resolved){const auto& m=ops[i].matrix;if(std::string(row.name)=="Au")expected=m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);}
            else for(int k=0;k<3;++k){Vector transformed{};for(int r=0;r<3;++r)for(int c=0;c<3;++c)transformed[r]+=ops[i].matrix[3*r+c]*axes.xyz[k][c];const double sign=dot(axes.xyz[k],transformed);if(std::abs(std::abs(sign)-1)>2e-5){same=false;break;}if(row.powers[k])expected*=sign;}
            if(!std::isfinite(chars[i])||std::abs(chars[i]-expected)>tolerance){same=false;break;}
        }if(same)return row.name;
    }return {};
}
} // namespace cov::orbital_characters
