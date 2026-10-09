#include "cov/point_group_irreps.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <tuple>

namespace cov {
namespace {
using Vec = std::array<double, 3>;
using Mat = std::array<double, 9>;
constexpr Mat identity{1,0,0,0,1,0,0,0,1};
constexpr double pi = std::numbers::pi;

double dot(const Vec& a, const Vec& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec cross(const Vec& a, const Vec& b) { return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
bool unit(Vec& a) { const double n=std::sqrt(dot(a,a)); if(!std::isfinite(n)||n<1e-12)return false; for(auto& x:a)x/=n;return true; }
Vec canonical_axis(Vec a) {
    unit(a);
    // An unoriented physical axis acquires a reproducible input-frame sign.
    int k=0;for(int j=1;j<3;++j)if(std::abs(a[j])>std::abs(a[k])+1e-10)k=j;
    if(a[k]<0)for(auto& x:a)x=-x;return a;
}
Mat multiply(const Mat& a,const Mat& b) { Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c; }
Mat transpose(const Mat& a) { return {a[0],a[3],a[6],a[1],a[4],a[7],a[2],a[5],a[8]}; }
double error(const Mat& a,const Mat& b) { double e=0;for(int i=0;i<9;++i)e=std::max(e,std::abs(a[i]-b[i]));return e; }
double det(const Mat& a) { return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]); }
double trace(const Mat& a) { return a[0]+a[4]+a[8]; }
Mat negative(Mat a) { for(auto& x:a)x=-x;return a; }
Mat rotation(const Vec& z,double theta) {
    const double c=std::cos(theta),s=std::sin(theta);const Mat k{0,-z[2],z[1],z[2],0,-z[0],-z[1],z[0],0};Mat a{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[3*i+j]=c*identity[3*i+j]+(1-c)*z[i]*z[j]+s*k[3*i+j];return a;
}
Mat mirror(const Vec& normal) { Mat a=identity;for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[3*i+j]-=2*normal[i]*normal[j];return a; }
Vec fixed_axis(const Mat& a) {
    Vec best{};double size=0;
    // Nullspace of R-I, including half-turns where the skew part vanishes.
    const Vec r[3]{{a[0]-1,a[1],a[2]},{a[3],a[4]-1,a[5]},{a[6],a[7],a[8]-1}};
    for(int i=0;i<3;++i)for(int j=i+1;j<3;++j){auto v=cross(r[i],r[j]);if(dot(v,v)>size){size=dot(v,v);best=v;}}
    return canonical_axis(best);
}
struct Group { std::string name;char family=0;std::size_t n=0;char suffix=0;std::size_t order=0; };
Group parse(std::string_view text) {
    std::string s;for(unsigned char c:text)if(!std::isspace(c)&&c!='_')s+=char(std::tolower(c));
    if(s=="c1")return {"C1",'a',1,0,1};
    if(s=="cs"||s=="c1h"||s=="c1v"||s=="s1")return {"Cs",'a',1,'s',2};
    if(s=="ci"||s=="s2")return {"Ci",'a',1,'i',2};
    for(const auto& p:std::vector<std::pair<std::string,std::size_t>>{{"t",12},{"td",24},{"th",24},{"o",24},{"oh",48},{"i",60},{"ih",120}})
        if(s==p.first){std::string name=s;name[0]=char(std::toupper(name[0]));return {name,'p',0,0,p.second};}
    if(s.size()<2||(s[0]!='c'&&s[0]!='d'&&s[0]!='s'))return {};
    std::size_t n=0,k=1;for(;k<s.size()&&std::isdigit(static_cast<unsigned char>(s[k]));++k){if(n>(std::numeric_limits<std::size_t>::max()-9)/10)return {};n=n*10+std::size_t(s[k]-'0');}
    if(n<2||k==1||s.size()>k+1)return {};
    const char suffix=k<s.size()?s[k]:0,family=s[0];std::size_t factor=1;
    if(family=='c'){if(suffix&&suffix!='h'&&suffix!='v')return {};factor=suffix?2:1;}
    if(family=='d'){if(suffix&&suffix!='h'&&suffix!='d')return {};factor=suffix?4:2;}
    if(family=='s'){if(suffix)return {};factor=n%2?2:1;}
    if(n>std::numeric_limits<std::size_t>::max()/factor)return {};
    std::string name(1,char(std::toupper(family)));name+=std::to_string(n);if(suffix)name+=suffix;
    return {name,family,n,suffix,n*factor};
}

// Index only the first matrix row to avoid cubic closure searches. Adjacent
// bins are examined because validated floating-point matrices can straddle a
// rounding boundary. Equality is always checked on all nine matrix entries.
class MatrixIndex {
    using Key=std::array<long long,3>;
    const std::vector<SymmetryOperation>& ops;
    double tol;
    std::map<Key,std::vector<std::size_t>> buckets;
    Key key(const Mat& a)const{return {std::llround(a[0]/tol),std::llround(a[1]/tol),std::llround(a[2]/tol)};}
public:
    MatrixIndex(const std::vector<SymmetryOperation>& operations,double tolerance):ops(operations),tol(tolerance){for(std::size_t i=0;i<ops.size();++i)buckets[key(ops[i].matrix)].push_back(i);}
    std::size_t find(const Mat& a)const{
        const auto b=key(a);
        for(int x=-1;x<=1;++x)for(int y=-1;y<=1;++y)for(int z=-1;z<=1;++z){const auto it=buckets.find({b[0]+x,b[1]+y,b[2]+z});if(it!=buckets.end())for(auto i:it->second)if(error(a,ops[i].matrix)<tol)return i;}
        return ops.size();
    }
};

bool validate_operations(PointGroupIrrepTable& table,const std::vector<SymmetryOperation>& ops,double tol) {
    MatrixIndex index(ops,tol);std::size_t e=ops.size();
    for(std::size_t i=0;i<ops.size();++i){const auto& a=ops[i].matrix;
        if(!std::all_of(a.begin(),a.end(),[](double x){return std::isfinite(x);})||error(multiply(a,transpose(a)),identity)>tol||std::abs(std::abs(det(a))-1)>tol){table.reason="operation matrix is not finite and orthogonal";return false;}
        if(index.find(a)!=i){table.reason="duplicate operation matrices";return false;}
        if(error(a,identity)<tol)e=i;
    }
    if(e==ops.size()){table.reason="identity operation is missing";return false;}
    table.identity_operation=e;
    for(const auto& a:ops)for(const auto& b:ops)if(index.find(multiply(a.matrix,b.matrix))==ops.size()){table.reason="operation matrices are not a closed finite group";return false;}
    return true;
}

struct BaseRow { std::string label;std::size_t dimension;double norm;int frequency;int sign; };
std::vector<BaseRow> cyclic_rows(std::size_t n) {
    std::vector<BaseRow> rows{{"A",1,1,0,1}};
    if(n%2==0)rows.push_back({"B",1,1,int(n/2),1});
    const std::size_t ecount=(n-1)/2;
    for(std::size_t j=1;j<=ecount;++j)rows.push_back({ecount==1?"E":"E"+std::to_string(j),2,2,int(j),1});
    return rows;
}
std::vector<BaseRow> dihedral_rows(std::size_t n,bool d2_names) {
    std::vector<BaseRow> rows{{"A1",1,1,0,1},{"A2",1,1,0,-1}};
    if(n%2==0){rows.push_back({"B1",1,1,int(n/2),1});rows.push_back({"B2",1,1,int(n/2),-1});}
    const std::size_t ecount=(n-1)/2;
    for(std::size_t j=1;j<=ecount;++j)rows.push_back({ecount==1?"E":"E"+std::to_string(j),2,1,int(j),0});
    if(n==2&&d2_names){rows[0].label="A";rows[1].label="B1";rows[2].label="B3";rows[3].label="B2";}
    return rows;
}

bool select_axes(PointGroupIrrepTable& table,const Group& g,const std::vector<SymmetryOperation>& ops,const PointGroupIrrepOptions& options) {
    Vec z{};const double tol=options.tolerance;
    if(options.principal_axis){z=*options.principal_axis;if(!unit(z)){table.reason="principal axis is zero or nonfinite";return false;}}
    else {
        // Find an actual minimal-angle generator. Choosing its unsigned axis
        // by input Cartesian preference handles the equivalent D2 axes.
        double score=-1;std::vector<Vec> candidates;
        const bool improper=(g.family=='s'||(g.family=='d'&&g.suffix=='d'&&g.n%2==0));
        for(const auto& op:ops){Mat a=op.matrix;
            if(improper){if(det(a)>0)continue;a=negative(a);}else if(det(a)<0)continue;
            if(error(a,identity)<tol)continue;
            const double c=(trace(a)-1)/2;
            // For improper generators -S, the axis eigenvalue is +1 but the
            // angle is pi-theta: select the target, not its higher powers.
            const double target=improper?-std::cos(2*pi/(g.family=='s'?double(g.n):2*double(g.n))):std::cos(2*pi/double(g.n));
            const double candidate_score=-std::abs(c-target);
            if(candidate_score>score+tol){score=candidate_score;candidates.clear();}
            if(std::abs(candidate_score-score)<tol)candidates.push_back(fixed_axis(a));
        }
        if(candidates.empty()||score<-tol*4){table.reason="point-group principal generator is missing";return false;}
        z=candidates.front();for(const auto& a:candidates)if(std::make_tuple(std::abs(a[2]),std::abs(a[1]),std::abs(a[0]),a)>std::make_tuple(std::abs(z[2]),std::abs(z[1]),std::abs(z[0]),z))z=a;
        if(!unit(z)){table.reason="principal rotation axis is unresolved";return false;}
    }
    Vec seed=options.reference_axis.value_or(Vec{1,0,0});
    const double zz=dot(seed,z);for(int i=0;i<3;++i)seed[i]-=zz*z[i];
    if(!unit(seed)){
        if(options.reference_axis){table.reason="reference axis is parallel to the principal axis";return false;}
        seed={0,1,0};const double p=dot(seed,z);for(int i=0;i<3;++i)seed[i]-=p*z[i];if(!unit(seed)){table.reason="reference axis is unresolved";return false;}
    }
    Vec x=seed;
    const bool dihedral=(g.family=='d'||(g.family=='c'&&g.suffix=='v'));
    if(dihedral){
        double best=-2;bool found=false;
        for(const auto& op:ops){Vec candidate{};
            if(g.family=='c'){
                if(det(op.matrix)>0||std::abs(trace(op.matrix)-1)>tol)continue;
                const auto normal=fixed_axis(negative(op.matrix));if(std::abs(dot(normal,z))>tol)continue;
                candidate=cross(normal,z);
            }else{
                if(det(op.matrix)<0||std::abs(trace(op.matrix)+1)>tol)continue;
                candidate=fixed_axis(op.matrix);if(std::abs(dot(candidate,z))>tol)continue;
            }
            if(!unit(candidate))continue;if(dot(candidate,seed)<0)for(auto& a:candidate)a=-a;
            const double score=dot(candidate,seed);
            if(score>best+1e-10||(std::abs(score-best)<=1e-10&&candidate>x)){best=score;x=candidate;found=true;}
        }
        if(!found){table.reason="perpendicular twofold axis or vertical mirror is missing";return false;}
        // Supplied axes name physical operations, not merely nearby planes.
        if(options.reference_axis&&std::abs(best-1)>tol){table.reason="reference axis does not name a group twofold axis or vertical plane";return false;}
    }
    table.principal_axis=z;table.reference_axis=x;
    std::ostringstream detail;detail<<g.name<<": z=("<<z[0]<<','<<z[1]<<','<<z[2]<<"), x=("<<x[0]<<','<<x[1]<<','<<x[2]<<"); ";
    if(g.family=='c'&&g.suffix=='v')detail<<"sigma_v(xz) is the reference mirror; B1 is even on it";
    else if(g.family=='d')detail<<"C2(x) is the reference perpendicular twofold; even-n B1 is even on it; D2 B1/B2/B3 transform as z/y/x";
    else detail<<"rotation powers are measured about z";
    if(!options.reference_axis&&dihedral)detail<<"; equivalent axes resolved by nearest input Cartesian x, then deterministic Cartesian tie";
    table.axis_detail=detail.str();return true;
}

bool axial(PointGroupIrrepTable& table,const Group& g,const std::vector<SymmetryOperation>& ops,const PointGroupIrrepOptions& options) {
    if(!select_axes(table,g,ops,options))return false;
    const Vec z=table.principal_axis,x=table.reference_axis;
    std::size_t n=g.n;bool dih=g.family=='d'||(g.family=='c'&&g.suffix=='v');bool extend=false,invert=false;
    bool sn_parity=false;
    Mat generator=rotation(z,2*pi/double(n)),second=identity,central=identity;
    if(g.family=='c'&&g.suffix=='v')second=mirror(cross(z,x));
    if(g.family=='d')second=rotation(x,pi);
    if(g.suffix=='h'){
        extend=true;invert=n%2==0;central=invert?negative(identity):mirror(z);
    }else if(g.family=='d'&&g.suffix=='d'){
        if(n%2){extend=true;invert=true;central=negative(identity);}
        else{n*=2;generator=multiply(rotation(z,2*pi/double(n)),mirror(z));}
    }else if(g.family=='s'){
        if(n%2){extend=true;central=mirror(z);generator=rotation(z,2*pi/double(n));}
        else if(n%4==2){n/=2;extend=true;invert=true;sn_parity=true;central=negative(identity);generator=rotation(z,2*pi/double(n));}
        else generator=multiply(generator,mirror(z));
    }
    MatrixIndex index(ops,options.tolerance);std::vector<std::tuple<std::size_t,int,int>> codes(ops.size());std::vector<bool> seen(ops.size());
    Mat power=identity;
    for(std::size_t k=0;k<n;++k){for(int f=0;f<(dih?2:1);++f)for(int h=0;h<(extend?2:1);++h){Mat a=power;if(f)a=multiply(a,second);if(h)a=multiply(a,central);const auto i=index.find(a);
        if(i==ops.size()||seen[i]){table.reason="advertised point group does not match its matrix generators and axis convention";return false;}seen[i]=true;codes[i]={k,f,h};}
        power=multiply(power,generator);
    }
    if(error(power,identity)>options.tolerance||std::count(seen.begin(),seen.end(),true)!=std::ptrdiff_t(ops.size())){table.reason="matrix generator order or group completeness disagrees with the point group";return false;}
    auto base=dih?dihedral_rows(n,g.family=='d'&&g.n==2&&g.suffix!='d'):cyclic_rows(n);
    for(int parity=0;parity<(extend?2:1);++parity)for(const auto& row:base){
        std::string label=row.label;
        if(extend){
            if(invert)label+=parity?"u":"g";else label+=parity?"''":"'";
        }
        PointGroupIrrepRow out{label,row.dimension,row.norm,{}};out.characters.resize(ops.size());
        for(std::size_t i=0;i<ops.size();++i){const auto [k,f,h]=codes[i];double value=row.dimension==2?2*std::cos(2*pi*double(row.frequency)*double(k)/double(n)):std::cos(2*pi*double(row.frequency)*double(k)/double(n));
            if(f)value=row.dimension==2?0:value*row.sign;if(h&&parity)value=-value;
            out.characters[i]=std::abs(value)<1e-12?0:value;
        }
        table.rows.push_back(std::move(out));
    }
    if(sn_parity)table.axis_detail+="; S(4k+2) uses proper C(2k+1) frequencies and inversion g/u";
    return true;
}

bool polyhedral(PointGroupIrrepTable& table,const Group& g,const std::vector<SymmetryOperation>& ops,double tol) {
    const bool direct=g.name=="Th"||g.name=="Oh"||g.name=="Ih";
    const std::string base=direct?g.name.substr(0,1):g.name;
    // Chemical class conventions, calibrated against A. Gelessus' university
    // tables: https://symmetry.constructor.university/ . I/Ih are also derived
    // in MIT 18.702: https://math.mit.edu/classes/18.702/summary-feb22.pdf .
    // O columns E,8C3,6C2(edge),6C4,3C2(face); Td E,8C3,3C2,6S4,6sigma_d.
    struct Row { const char* name;std::size_t dim;double norm;std::vector<double> values; };
    std::vector<Row> rows;std::vector<std::size_t> expected;
    const double phi=(1+std::sqrt(5.0))/2;
    if(base=="T"){rows={{"A",1,1,{1,1,1}},{"E",2,2,{2,-1,2}},{"T",3,1,{3,0,-1}}};expected={1,8,3};}
    else if(base=="Td"){rows={{"A1",1,1,{1,1,1,1,1}},{"A2",1,1,{1,1,1,-1,-1}},{"E",2,1,{2,-1,2,0,0}},{"T1",3,1,{3,0,-1,1,-1}},{"T2",3,1,{3,0,-1,-1,1}}};expected={1,8,3,6,6};}
    else if(base=="O"){rows={{"A1",1,1,{1,1,1,1,1}},{"A2",1,1,{1,1,-1,-1,1}},{"E",2,1,{2,-1,0,0,2}},{"T1",3,1,{3,0,-1,1,-1}},{"T2",3,1,{3,0,1,-1,-1}}};expected={1,8,6,6,3};}
    else if(base=="I"){rows={{"A",1,1,{1,1,1,1,1}},{"T1",3,1,{3,phi,1-phi,0,-1}},{"T2",3,1,{3,1-phi,phi,0,-1}},{"G",4,1,{4,-1,-1,1,0}},{"H",5,1,{5,0,0,-1,1}}};expected={1,12,12,20,15};}
    else{table.reason="unsupported polyhedral group";return false;}
    std::vector<Mat> face_halfturns;
    if(base=="O")for(const auto& op:ops){Mat a=op.matrix;if(det(a)<0)continue;if(std::abs(trace(a)-1)<tol)face_halfturns.push_back(multiply(a,a));}
    std::vector<std::size_t> classes(ops.size()),counts(expected.size()*(direct?2:1));
    std::vector<bool> odd(ops.size());
    for(std::size_t i=0;i<ops.size();++i){Mat a=ops[i].matrix;odd[i]=det(a)<0;
        if(odd[i]&&direct)a=negative(a);
        else if(odd[i]&&base!="Td"){table.reason="improper element in a proper polyhedral group";return false;}
        const double t=trace(a);std::size_t c=expected.size();
        if(error(a,identity)<tol)c=0;
        else if(base=="Td"&&odd[i]){if(std::abs(t+1)<tol)c=3;else if(std::abs(t-1)<tol)c=4;}
        else if(base=="T"||base=="Td"){if(std::abs(t)<tol)c=1;else if(std::abs(t+1)<tol)c=2;}
        else if(base=="O"){
            if(std::abs(t)<tol)c=1;else if(std::abs(t-1)<tol)c=3;else if(std::abs(t+1)<tol){c=2;for(const auto& b:face_halfturns)if(error(a,b)<tol){c=4;break;}}
        }else if(base=="I"){
            if(std::abs(t-phi)<tol)c=1;else if(std::abs(t-(1-phi))<tol)c=2;else if(std::abs(t)<tol)c=3;else if(std::abs(t+1)<tol)c=4;
        }
        if(c==expected.size()){table.reason="operation does not belong to an advertised polyhedral rotation class";return false;}
        classes[i]=c;++counts[c+(direct&&odd[i]?expected.size():0)];
    }
    for(std::size_t i=0;i<counts.size();++i)if(counts[i]!=expected[i%expected.size()]){table.reason="polyhedral operation class populations are incomplete or mislabeled";return false;}
    for(int parity=0;parity<(direct?2:1);++parity)for(const auto& row:rows){PointGroupIrrepRow out{std::string(row.name)+(direct?(parity?"u":"g"):""),row.dim,row.norm,{}};
        for(std::size_t i=0;i<ops.size();++i)out.characters.push_back(row.values[classes[i]]*(direct&&odd[i]&&parity?-1:1));table.rows.push_back(std::move(out));}
    table.axis_detail=g.name+": chemical classes from actual matrix traces; O face C2 distinguished by squares of C4; inversion parity from -R; T E is a real conjugate pair";
    return true;
}

bool validate_rows(PointGroupIrrepTable& table,const std::vector<SymmetryOperation>& ops,double tol) {
    double completeness=0;
    for(std::size_t a=0;a<table.rows.size();++a){const auto& row=table.rows[a];
        if(row.characters.size()!=ops.size()||std::abs(row.characters[table.identity_operation]-double(row.dimension))>tol){table.reason="irrep identity character does not equal dimension";return false;}
        completeness+=double(row.dimension*row.dimension)/row.character_norm;
        for(std::size_t b=0;b<=a;++b){const auto& other=table.rows[b];double inner=0;for(std::size_t k=0;k<ops.size();++k){const double value=row.characters[k];if(!std::isfinite(value)||std::abs(value)>double(row.dimension)+tol){table.reason="invalid irrep character";return false;}inner+=value*other.characters[k];}
            inner/=double(ops.size());if(std::abs(inner-(a==b?row.character_norm:0))>tol){table.reason="irrep character orthogonality or real-pair norm failed";return false;}}
    }
    if(std::abs(completeness-double(ops.size()))>tol){table.reason="irrep dimension completeness failed";return false;}
    // The sum d/norm * chi must reconstruct the real regular representation.
    for(std::size_t k=0;k<ops.size();++k){double value=0;for(const auto& row:table.rows)value+=double(row.dimension)/row.character_norm*row.characters[k];if(std::abs(value-(k==table.identity_operation?double(ops.size()):0))>tol*std::max(1.0,double(ops.size()))){table.reason="irrep regular-character completeness failed";return false;}}
    return true;
}
} // namespace

std::size_t finite_point_group_order(std::string_view group) { return parse(group).order; }

PointGroupIrrepTable finite_point_group_irreps(const MolecularSymmetry& symmetry,const PointGroupIrrepOptions& options) {
    PointGroupIrrepTable table;const auto g=parse(symmetry.point_group);table.point_group=g.name.empty()?symmetry.point_group:g.name;
    if(symmetry.linear||!g.order){table.reason="continuous, atomic, unsupported, or invalid point group; finite service excluded";return table;}
    if(!std::isfinite(options.tolerance)||options.tolerance<1e-10||options.tolerance>1e-2){table.reason="invalid operation tolerance";return table;}
    if(symmetry.operations.size()!=g.order){table.reason="operation count does not match finite point-group order";return table;}
    if(!validate_operations(table,symmetry.operations,options.tolerance))return table;
    bool good=false;
    if(g.family=='a'){
        table.rows.push_back({g.suffix=='s'?"A'":g.suffix=='i'?"Ag":"A",1,1,std::vector<double>(g.order,1)});
        if(g.order==2){const auto& a=symmetry.operations[1-table.identity_operation].matrix;const bool valid=g.suffix=='i'?error(a,negative(identity))<options.tolerance:det(a)<0&&std::abs(trace(a)-1)<options.tolerance;
            if(!valid){table.reason="nonaxial point-group operation does not match its name";return table;}
            std::vector<double> chars(g.order,-1);chars[table.identity_operation]=1;table.rows.push_back({g.suffix=='s'?"A''":"Au",1,1,std::move(chars)});}
        table.axis_detail=g.name+": actual identity/reflection/inversion matrices";good=true;
    }else if(g.family=='p')good=polyhedral(table,g,symmetry.operations,options.tolerance);
    else good=axial(table,g,symmetry.operations,options);
    if(good&&validate_rows(table,symmetry.operations,options.tolerance)){table.valid=true;table.reason="complete finite real chemical character table validated";}
    return table;
}

PointGroupCharacterDecomposition decompose_point_group_characters(const PointGroupIrrepTable& table,const std::vector<double>& characters,double tolerance) {
    PointGroupCharacterDecomposition result;
    if(!table.valid||table.rows.empty()){result.reason="finite character table is not valid";return result;}
    if(!std::isfinite(tolerance)||tolerance<=0||characters.size()!=table.rows.front().characters.size()||!std::all_of(characters.begin(),characters.end(),[](double x){return std::isfinite(x);})){result.reason="invalid character vector or tolerance";return result;}
    const double dimension=characters[table.identity_operation];
    if(dimension<0||dimension>double(std::numeric_limits<std::size_t>::max())||std::abs(dimension-std::round(dimension))>tolerance){result.reason="identity character is not a nonnegative integer dimension";return result;}
    for(const auto& row:table.rows){double inner=0;for(std::size_t i=0;i<characters.size();++i)inner+=row.characters[i]*characters[i];const double coefficient=inner/(double(characters.size())*row.character_norm);
        if(!std::isfinite(coefficient)||coefficient<-tolerance||coefficient>double(std::numeric_limits<std::size_t>::max())||std::abs(coefficient-std::round(coefficient))>tolerance){result.reason="characters do not contain integer complete real irrep copies";return result;}
        const auto copies=std::size_t(std::max(0.0,std::round(coefficient)));result.multiplicities.push_back(copies);
        if(copies>(std::numeric_limits<std::size_t>::max()-result.dimension)/row.dimension){result.reason="decomposition dimension overflow";return result;}result.dimension+=copies*row.dimension;
    }
    if(std::abs(double(result.dimension)-dimension)>tolerance){result.reason="decomposed irrep dimensions do not match the identity character";return result;}
    for(std::size_t k=0;k<characters.size();++k){double reconstruction=0;for(std::size_t j=0;j<table.rows.size();++j)reconstruction+=double(result.multiplicities[j])*table.rows[j].characters[k];result.maximum_reconstruction_error=std::max(result.maximum_reconstruction_error,std::abs(reconstruction-characters[k]));}
    if(result.maximum_reconstruction_error>tolerance){result.reason="irrep characters do not reconstruct every operation character";return result;}
    result.valid=true;result.reason="integer real-irrep multiplicities, dimensions and all operation characters reconstructed";return result;
}
} // namespace cov
