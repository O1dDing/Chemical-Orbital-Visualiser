#pragma once

#include "cov/point_group_irreps.hpp"
#include "cov/d2h_orbital_characters.hpp"

namespace cov {
// One recorded coordinate convention for canonical and fragment analyses.
// A tied geometrical convention falls back to the documented input-axis
// convention of the table builder, not to an unsupported/unknown point group.
inline PointGroupIrrepTable molecular_point_group_irreps(
    const Wavefunction& w,const MolecularSymmetry& geometry) {
    PointGroupIrrepOptions options;
    using namespace orbital_characters;
    if(geometry.point_group=="D2h") {
        const auto axes=d2h_axes(w,geometry.operations,geometry.centre_bohr,geometry.tolerance_bohr);
        if(axes.resolved){options.principal_axis=axes.xyz[2];options.reference_axis=axes.xyz[0];}
    }
    if(geometry.point_group=="C2v") {
        const SymmetryOperation* rotation=nullptr;
        std::vector<const SymmetryOperation*> mirrors;
        const auto determinant=[](const Matrix& m){return m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);};
        for(const auto& op:geometry.operations){const auto& m=op.matrix;const double tr=m[0]+m[4]+m[8];
            if(std::abs(determinant(m)-1)<2e-5&&std::abs(tr+1)<2e-5)rotation=&op;
            if(std::abs(determinant(m)+1)<2e-5&&std::abs(tr-1)<2e-5)mirrors.push_back(&op);}
        if(rotation&&mirrors.size()==2){std::size_t support[2]{};
            for(int j=0;j<2;++j)for(std::size_t a=0;a<mirrors[j]->atom_permutation.size();++a)if(mirrors[j]->atom_permutation[a]==a)++support[j];
            if(support[0]!=support[1]){const auto& yz=*mirrors[support[0]>support[1]?0:1];Vector x{},z{};double nx=0,nz=0;
                for(int c=0;c<3;++c){Vector a{},b{};for(int r=0;r<3;++r){a[r]=(r==c?1.:0.)-yz.matrix[3*r+c];b[r]=(r==c?1.:0.)+rotation->matrix[3*r+c];}
                    if(dot(a,a)>nx){nx=dot(a,a);x=a;}if(dot(b,b)>nz){nz=dot(b,b);z=b;}}
                if(nx>1e-20&&nz>1e-20){for(auto& a:x)a/=std::sqrt(nx);for(auto& a:z)a/=std::sqrt(nz);options.principal_axis=z;options.reference_axis=x;}
            }
        }
    }
    // Cnv mirror classes are named by their actual atom support where the
    // geometry distinguishes them, so an input rotation cannot swap B1/B2.
    if(geometry.point_group.size()>2&&geometry.point_group.front()=='C'&&geometry.point_group.back()=='v'&&geometry.point_group!="C2v"){
        const auto base=finite_point_group_irreps(geometry);
        if(base.valid){const auto z=base.principal_axis;std::size_t best=0,least=w.atoms.size()+1;const SymmetryOperation* mirror=nullptr;
            for(const auto& op:geometry.operations){const auto& m=op.matrix;const double determinant=m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);
                if(std::abs(determinant+1)>2e-5||std::abs(m[0]+m[4]+m[8]-1)>2e-5)continue;
                std::size_t support=0;for(std::size_t i=0;i<op.atom_permutation.size();++i)if(op.atom_permutation[i]==i)++support;
                least=std::min(least,support);if(support>best){best=support;mirror=&op;}}
            if(mirror&&best>least){Vector x{};double norm=0;for(int c=0;c<3;++c){Vector v{};for(int r=0;r<3;++r)v[r]=mirror->matrix[3*r+c]+(r==c?1.:0.)-2*z[r]*z[c];
                    if(dot(v,v)>norm){norm=dot(v,v);x=v;}}
                if(norm>1e-20){for(auto& a:x)a/=std::sqrt(norm);options.principal_axis=z;options.reference_axis=x;}}
        }
    }
    // In even Dnh rings, C2 through atoms and C2 through gaps are different
    // classes. Prefer the atom-supported class when geometry distinguishes it;
    // this preserves B names under a rigid rotation of the input coordinates.
    if(geometry.point_group.size()>2&&geometry.point_group.front()=='D'&&geometry.point_group!="D2h"){
        const auto base=finite_point_group_irreps(geometry);
        if(base.valid){const auto z=base.principal_axis;std::size_t best=0;const SymmetryOperation* reference=nullptr;
            for(const auto& op:geometry.operations){const auto& m=op.matrix;const double determinant=m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);
                if(std::abs(determinant-1)>2e-5||std::abs(m[0]+m[4]+m[8]+1)>2e-5)continue;
                Vector transformed{};for(int r=0;r<3;++r)for(int c=0;c<3;++c)transformed[r]+=m[3*r+c]*z[c];
                if(dot(transformed,z)>-1+2e-5)continue;
                std::size_t support=0;for(std::size_t i=0;i<op.atom_permutation.size();++i)if(op.atom_permutation[i]==i)++support;
                if(support>best){best=support;reference=&op;}}
            if(reference){Vector x{};double norm=0;for(int c=0;c<3;++c){Vector v{};for(int r=0;r<3;++r)v[r]=reference->matrix[3*r+c]+(r==c?1.:0.);if(dot(v,v)>norm){norm=dot(v,v);x=v;}}
                if(norm>1e-20){for(auto& a:x)a/=std::sqrt(norm);options.principal_axis=z;options.reference_axis=x;}}
        }
    }
    return finite_point_group_irreps(geometry,options);
}

inline MolecularSymmetry complete_molecular_point_group(
    const Wavefunction& w,const MolecularSymmetry& geometry,std::size_t maximum_order=256) {
    auto result=geometry;
    const auto expected=finite_point_group_order(geometry.point_group);
    if(geometry.linear||!expected||expected>maximum_order){result.operations.clear();return result;}
    using namespace orbital_characters;
    SymmetryOperation e;e.atom_permutation.resize(w.atoms.size());std::iota(e.atom_permutation.begin(),e.atom_permutation.end(),0);
    result.operations={e};
    for(std::size_t i=0;i<result.operations.size();++i)for(const auto& g:geometry.operations){
        if(g.atom_permutation.size()!=w.atoms.size()){result.operations.clear();return result;}
        SymmetryOperation p;p.matrix=multiply(result.operations[i].matrix,g.matrix);p.atom_permutation.resize(w.atoms.size());
        for(std::size_t a=0;a<w.atoms.size();++a){if(g.atom_permutation[a]>=w.atoms.size()){result.operations.clear();return result;}p.atom_permutation[a]=result.operations[i].atom_permutation[g.atom_permutation[a]];}
        if(std::any_of(result.operations.begin(),result.operations.end(),[&](const auto& q){return p.atom_permutation==q.atom_permutation&&matrix_error(p.matrix,q.matrix)<2e-5;}))continue;
        if(result.operations.size()>=expected){result.operations.clear();return result;}
        for(std::size_t a=0;a<w.atoms.size();++a){const auto& x=w.atoms[a];const auto& y=w.atoms[p.atom_permutation[a]];
            Vector source{x.x-geometry.centre_bohr[0],x.y-geometry.centre_bohr[1],x.z-geometry.centre_bohr[2]},target{y.x-geometry.centre_bohr[0],y.y-geometry.centre_bohr[1],y.z-geometry.centre_bohr[2]},diff{};
            for(int r=0;r<3;++r){diff[r]=-target[r];for(int c=0;c<3;++c)diff[r]+=p.matrix[3*r+c]*source[c];}
            p.max_mapping_error_bohr=std::max(p.max_mapping_error_bohr,std::sqrt(dot(diff,diff)));
            if(x.atomic_number!=y.atomic_number||p.max_mapping_error_bohr>geometry.tolerance_bohr){result.operations.clear();return result;}
        }
        result.operations.push_back(std::move(p));
    }
    if(result.operations.size()!=expected)result.operations.clear();return result;
}
} // namespace cov
