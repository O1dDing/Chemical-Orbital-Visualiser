#include "cov/formchk.hpp"

#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* name) {
    if (!condition) throw std::runtime_error(name);
}

struct DirectoryGuard {
    std::filesystem::path path;
    ~DirectoryGuard() {
        std::error_code error;
        std::filesystem::remove_all(path,error);
    }
};

std::string read_line(std::ifstream& file) {
    std::string line;
    std::getline(file,line);
    return line;
}

void check_exit_and_literals() {
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    DirectoryGuard directory{std::filesystem::temp_directory_path() /
        ("cov_formchk_spawn_"+std::to_string(getpid())+"_"+std::to_string(stamp))};
    std::filesystem::create_directories(directory.path);
    const auto checkpoint=directory.path / "wave ' $ ü 漢字.chk";
    const auto executable=directory.path / "fake ' $ ü 漢字 formchk";
    const auto capture=directory.path / "literal arguments.txt";
    { std::ofstream file(checkpoint); file << "binary placeholder"; }
    {
        std::ofstream file(executable);
        file << "#!/bin/sh\n"
                "printf '%s\\n' \"$1\" \"$2\" > \"$COV_FORMCHK_CAPTURE\"\n"
                "exit 23\n";
    }
    require(chmod(executable.c_str(),0700)==0,"make fake formchk executable");
    require(setenv("COV_FORMCHK",executable.c_str(),1)==0,"set formchk executable");
    require(setenv("COV_FORMCHK_CAPTURE",capture.c_str(),1)==0,"set capture path");
    std::string failure;
    try { (void)cov::parse_gaussian_chk_via_formchk(checkpoint); }
    catch (const std::runtime_error& error) { failure=error.what(); }
    require(failure.find("exit code 23")!=std::string::npos,"real exit status");
    std::ifstream captured(capture);
    require(static_cast<bool>(captured),"formchk received literal arguments");
    require(read_line(captured)==checkpoint.string(),"literal CHK path");
    const auto temporary=read_line(captured);
    require(temporary.find("cov_formchk_")!=std::string::npos &&
            temporary.find(".fchk")!=std::string::npos,"temporary output argument");
    require(!std::filesystem::exists(temporary),"temporary output cleaned");

    require(setenv("COV_FORMCHK","no-such-formchk-executable-for-cov-smoke",1)==0,
            "set missing executable");
    failure.clear();
    try { (void)cov::parse_gaussian_chk_via_formchk(checkpoint); }
    catch (const std::runtime_error& error) { failure=error.what(); }
    require(failure.find("Could not start Gaussian formchk")!=std::string::npos,
            "missing executable reported");
}

} // namespace

int main() {
    try {
        check_exit_and_literals();
        std::cout << "formchk_spawn_smoke: literal argv and exit status verified\n";
    } catch (const std::exception& error) {
        std::cerr << "formchk_spawn_smoke: " << error.what() << '\n';
        return 1;
    }
}
#else
int main() { return 0; }
#endif
