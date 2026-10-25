#include "cov/file_dialog.hpp"

#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <iterator>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
struct DirectoryGuard {
    std::filesystem::path path;
    ~DirectoryGuard() { std::error_code e; std::filesystem::remove_all(path,e); }
};
void executable(const std::filesystem::path& path) {
    std::ofstream file(path);
    file << "#!/bin/sh\n"
            "printf '%s\\n' \"$@\" > \"$COV_DIALOG_CAPTURE\"\n"
            "if [ \"$COV_DIALOG_EXIT\" != 0 ]; then exit \"$COV_DIALOG_EXIT\"; fi\n"
            "printf '%s\\n' \"$COV_DIALOG_SELECTION\"\n";
    file.close();
    require(chmod(path.c_str(),0700)==0,"executable picker stub");
}
std::string captured(const std::filesystem::path& path) {
    std::ifstream file(path);
    return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}
void check() {
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    DirectoryGuard directory{std::filesystem::temp_directory_path() /
        ("cov_dialog_"+std::to_string(getpid())+"_"+std::to_string(stamp))};
    std::filesystem::create_directories(directory.path);
    const auto capture=directory.path/"arguments.txt";
    const auto selection=directory.path/"mol ' $ 漢字.covnbopkg";
    setenv("PATH",directory.path.c_str(),1);
    setenv("DISPLAY",":99",1);
    setenv("COV_DIALOG_CAPTURE",capture.c_str(),1);
    setenv("COV_DIALOG_SELECTION",selection.c_str(),1);
    setenv("COV_DIALOG_EXIT","0",1);
#if defined(__APPLE__)
    executable(directory.path/"osascript");
#else
    executable(directory.path/"zenity");
#endif
    const auto selected=cov::open_wavefunction_file_dialog(cov::ui::Language::French,true);
    require(selected.selected() && selected.path==selection,"literal Unicode path from picker");
    const auto nbo=captured(capture);
    require(nbo.find("Ouvrir des données NBO")!=std::string::npos,"French NBO title");
#if !defined(__APPLE__)
    require(nbo.find("*.log *.out *.nbo *.47")!=std::string::npos,"NBO file filter");
    std::filesystem::remove(directory.path/"zenity");
    executable(directory.path/"kdialog");
#endif
    const auto calculation=cov::open_wavefunction_file_dialog(cov::ui::Language::Japanese,false);
    require(calculation.selected(),"fallback picker opens calculation");
    const auto normal=captured(capture);
    require(normal.find("計算ファイルを開く")!=std::string::npos,"Japanese calculation title");
#if !defined(__APPLE__)
    require(normal.find("*.fchk")!=std::string::npos && normal.find("*.covnbopkg")!=std::string::npos,
            "calculation filter includes package input");
    require(normal.find("--title\n")!=std::string::npos,"fallback title argument");
#endif
    setenv("COV_DIALOG_EXIT","1",1);
    require(cov::open_wavefunction_file_dialog(cov::ui::Language::French).cancelled,"cancel is not failure");
    setenv("COV_DIALOG_EXIT","23",1);
    const auto failure=cov::open_wavefunction_file_dialog(cov::ui::Language::ChineseSimplified);
    require(failure.error.find("文件选择窗口")!=std::string::npos && failure.error.find("23")!=std::string::npos,
            "localized failure keeps exit status");
}
}
int main() {
    try { check(); std::cout << "file_dialog_smoke passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
#else
int main() { return 0; }
#endif
