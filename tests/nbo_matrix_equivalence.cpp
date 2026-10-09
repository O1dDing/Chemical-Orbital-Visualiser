// Frozen pre-optimization reader from bbadfa3, retained only as a regression oracle.
#include "cov/nbo_integration.hpp"
#include "cov/wavefunction_io.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
namespace reference {
using namespace cov;
const std::string num=R"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)";
std::string trim(std::string s) { auto a=s.find_first_not_of(" \t\r\n"); return a==s.npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1); }
std::string upper(std::string s) { for(auto& c:s)c=static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return s; }
double number(std::string s) { for(auto& c:s)if(c=='D'||c=='d')c='E'; double x=std::stod(s); if(!std::isfinite(x))throw std::runtime_error("NBO nonfinite number"); return x; }
std::vector<double> numbers(const std::string& s) { static const std::regex r(num); std::vector<double> v; for(std::sregex_iterator i(s.begin(),s.end(),r),e;i!=e;++i)v.push_back(number(i->str())); return v; }
std::vector<std::string> lines(const std::filesystem::path& p) { std::ifstream f(p); if(!f)throw std::runtime_error("Cannot read NBO file: "+p.string()); std::vector<std::string> l; std::string s; while(std::getline(f,s))l.push_back(s); return l; }
NboSource source(const std::filesystem::path& p,std::string block,std::size_t line,std::string raw={},std::string version={},std::size_t seg=0) { return {p.string(),std::move(block),std::move(version),std::move(raw),line,line,seg}; }
std::vector<NboMatrix> read_matrix_file(const std::filesystem::path& p,const std::string& kind,std::size_t n,bool open,std::size_t cols=0) {
    if(!cols)cols=n;
    auto ls=lines(p);
    // W37/PLOT appends occupancies and orbital/atom descriptors to each spin
    // block. Only the first NBAS*NBAS real values belong to the matrix.
    // OPEN files delimit their independent payloads with ALPHA/BETA SPIN.
    const std::regex numeric("^\\s*"+num+"(?:\\s+"+num+")*\\s*$");
    const std::map<std::string,std::string> headings{{"AONBO","NBOs in the AO basis:"},{"NBOMO","MOs in the NBO basis:"},{"AONAO","NAOs in the AO basis:"},{"NAOMO","MOs in the NAO basis:"},{"NAONBO","NBOs in the NAO basis:"},{"AONHO","NHOs in the AO basis:"},{"AONLMO","NLMOs in the AO basis:"},{"AOPNAO","PNAOs in the AO basis:"},{"NAONHO","NHOs in the NAO basis:"},{"NHONBO","NBOs in the NHO basis:"},{"NBONLMO","NLMOs in the NBO basis:"},{"NLMOMO","MOs in the NLMO basis:"},{"NAONLMO","NLMOs in the NAO basis:"},{"AOMO","MOs in the AO basis:"}};
    const auto hi=headings.find(kind);if(hi==headings.end())throw std::runtime_error("Unsupported explicit NBO matrix kind");const auto& heading=hi->second;
    std::size_t header=ls.size();for(std::size_t i=0;i<ls.size();++i)if(ls[i].find(heading)!=ls[i].npos){if(header!=ls.size())throw std::runtime_error("Repeated matrix header");header=i;}
    if(header==ls.size())throw std::runtime_error("NBO matrix type header absent or mismatched: "+kind);
    // NBO uses a shared NAO basis for both spins; AONAO can therefore have
    // one explicitly unpolarized payload while NAOMO/NAONBO have two.
    if(open&&(kind=="AONAO"||kind=="AOPNAO")){bool spin_header=false;for(const auto& line:ls){const auto u=upper(trim(line));if(u=="ALPHA SPIN"||u=="BETA SPIN"||u=="BETA  SPIN")spin_header=true;}if(!spin_header)open=false;}
    std::vector<std::size_t> starts;if(open){for(std::size_t i=header+1;i<ls.size();++i){auto u=upper(trim(ls[i]));if(u=="ALPHA SPIN"||u=="BETA  SPIN"||u=="BETA SPIN")starts.push_back(i+1);}if(starts.size()!=2||upper(ls[starts[0]-1]).find("ALPHA")==std::string::npos||upper(ls[starts[1]-1]).find("BETA")==std::string::npos)throw std::runtime_error("OPEN matrix requires exactly one alpha and beta block");}else starts={header+1};
    const std::size_t count=n*cols;std::vector<NboMatrix> out;
    for(std::size_t s=0;s<starts.size();++s){std::vector<double> values;std::size_t first=0,last=0;const auto end=s+1<starts.size()?starts[s+1]-1:ls.size();for(std::size_t i=starts[s];i<end&&values.size()<count;++i){if(!std::regex_match(ls[i],numeric)){if(!values.empty())throw std::runtime_error("Interrupted NBO matrix payload");continue;}if(ls[i].find_first_of(".EeDd")==std::string::npos){if(values.empty())continue;throw std::runtime_error("Integer metadata encountered before complete NBO matrix");}if(!first)first=i+1;auto row=numbers(ls[i]);if(values.size()+row.size()>count)throw std::runtime_error("NBO matrix payload boundary not aligned");values.insert(values.end(),row.begin(),row.end());last=i+1;}if(values.size()!=count)throw std::runtime_error(kind+" incomplete matrix payload");NboMatrix m;m.kind=kind;m.spin=open?(s==0?NboSpin::Alpha:NboSpin::Beta):NboSpin::Total;m.rows=n;m.columns=cols;m.values.resize(count);m.source=source(p,kind,first);m.source.line_end=last;for(std::size_t j=0;j<cols;++j)for(std::size_t i=0;i<n;++i)m.values[i*cols+j]=values[j*n+i];out.push_back(std::move(m));}return out;
}
}
int main(int argc,char** argv) {try {
    if(argc!=3)throw std::runtime_error("Usage: cov_nbo_matrix_equivalence canonical.fchk package_directory");
    const auto canonical=cov::parse_wavefunction(argv[1]);
    const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1)throw std::runtime_error("Expected one candidate");
    const auto integrated=cov::read_nbo_integration(canonical,found.candidates[0]);
    std::set<std::string> seen;std::size_t matrices=0,numbers=0;
    for(const auto& input:integrated.dataset.matrices) {
        if(!seen.insert(input.kind+":"+input.source.path).second)continue;
        const bool open=integrated.dataset.archive->open_shell;
        const auto before=reference::read_matrix_file(input.source.path,input.kind,input.rows,open,input.columns);
        const auto after=cov::read_nbo_matrix_rectangular(input.source.path,input.kind,input.rows,input.columns,open);
        if(before.size()!=after.size())throw std::runtime_error("Spin block count changed");
        for(std::size_t b=0;b<before.size();++b) {
            const auto& x=before[b];const auto& y=after[b];
            if(x.rows!=y.rows||x.columns!=y.columns||x.spin!=y.spin||x.kind!=y.kind||
               x.source.path!=y.source.path||x.source.block!=y.source.block||
               x.source.line_begin!=y.source.line_begin||x.source.line_end!=y.source.line_end||
               x.values.size()!=y.values.size())throw std::runtime_error("Matrix identity or source boundary changed");
            for(std::size_t i=0;i<x.values.size();++i)
                if(std::bit_cast<std::uint64_t>(x.values[i])!=std::bit_cast<std::uint64_t>(y.values[i]))
                    throw std::runtime_error("Matrix value bits changed");
            ++matrices;numbers+=x.values.size();
        }
    }
    if(!matrices)throw std::runtime_error("No matrices checked");
    std::cout<<"{\"matrices\":"<<matrices<<",\"values\":"<<numbers<<",\"bitwise_equal\":true,\"source_boundaries_equal\":true}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
