#include "cov/ui.hpp"

#include <cstring>
#include <cstdlib>
#include <iostream>

int main() {
    using cov::ui::Language;
    using cov::ui::Text;

    const char* en_file=cov::ui::tr(Text::MoldenPath,Language::English);
    const char* en_idle=cov::ui::tr(Text::IdleHint,Language::English);
    const char* en_mo=cov::ui::tr(Text::MoldenMO,Language::English);
    const char* zh_file=cov::ui::tr(Text::MoldenPath,Language::ChineseSimplified);
    const char* en_tracking=cov::ui::tr(Text::FrameTracking,Language::English);
    const char* zh_fallback=cov::ui::tr(
        Text::ConservativeFallback,Language::ChineseSimplified);

    if (!en_file || std::strcmp(en_file,"Calculation file or folder")!=0 ||
        !en_idle || std::strstr(en_idle,"FCHK")==nullptr || std::strstr(en_idle,"NBO")==nullptr ||
        !en_mo || std::strcmp(en_mo,"Source MO")!=0 ||
        !zh_file || std::strcmp(zh_file,"计算文件或文件夹")!=0 ||
        !en_tracking || std::strcmp(en_tracking,"Frame continuity") != 0 ||
        !zh_fallback || std::strcmp(zh_fallback,"已使用简化匹配") != 0) {
        std::cerr << "Calculation-input UI wording regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "Calculation-input UI smoke test passed\n";
    return EXIT_SUCCESS;
}
