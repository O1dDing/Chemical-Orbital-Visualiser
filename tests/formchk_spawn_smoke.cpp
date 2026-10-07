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
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>

#include <chrono>
#include <cstdlib>
#include <cwchar>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* name) {
    if (!condition) throw std::runtime_error(name);
}

struct DirectoryGuard {
    std::filesystem::path path;
    ~DirectoryGuard() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

std::filesystem::path own_executable() {
    std::vector<wchar_t> buffer(MAX_PATH);
    for (;;) {
        const DWORD length = GetModuleFileNameW(nullptr, buffer.data(),
                                               static_cast<DWORD>(buffer.size()));
        require(length != 0, "locate fake formchk executable");
        if (length < buffer.size()) return std::filesystem::path(buffer.data());
        buffer.resize(buffer.size() * 2);
    }
}

void write_argument(std::ofstream& file, const wchar_t* argument) {
    const DWORD length = static_cast<DWORD>(std::wcslen(argument));
    file.write(reinterpret_cast<const char*>(&length), sizeof(length));
    file.write(reinterpret_cast<const char*>(argument), length * sizeof(wchar_t));
}

std::wstring read_argument(std::ifstream& file) {
    DWORD length = 0;
    file.read(reinterpret_cast<char*>(&length), sizeof(length));
    require(static_cast<bool>(file) && length < 32768, "read argument length");
    std::wstring argument(length, L'\0');
    file.read(reinterpret_cast<char*>(argument.data()), length * sizeof(wchar_t));
    require(static_cast<bool>(file), "read literal argument");
    return argument;
}

void check_exit_and_literals() {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    DirectoryGuard directory{std::filesystem::temp_directory_path() /
        (L"cov_formchk_spawn_" + std::to_wstring(GetCurrentProcessId()) + L"_" +
         std::to_wstring(stamp))};
    std::filesystem::create_directories(directory.path);
    const auto checkpoint = directory.path / L"wave %USERNAME% \u6F22\u5B57.chk";
    const auto executable = directory.path / L"fake %USERNAME% \u6F22\u5B57 formchk.exe";
    const auto capture = directory.path / L"literal arguments.bin";
    { std::ofstream file(checkpoint); file << "binary placeholder"; }
    std::filesystem::copy_file(own_executable(), executable);
    require(_wputenv_s(L"COV_FORMCHK", executable.c_str()) == 0,
            "set fake formchk executable");
    require(_wputenv_s(L"COV_FORMCHK_CAPTURE", capture.c_str()) == 0,
            "set capture path");

    std::string failure;
    try { (void)cov::parse_gaussian_chk_via_formchk(checkpoint); }
    catch (const std::runtime_error& error) { failure = error.what(); }
    require(failure.find("exit code 23") != std::string::npos, "real exit status");
    std::ifstream captured(capture, std::ios::binary);
    require(static_cast<bool>(captured), "formchk received literal arguments");
    require(read_argument(captured) == checkpoint.wstring(), "literal CHK path");
    const auto temporary = read_argument(captured);
    require(temporary.find(L"cov_formchk_") != std::wstring::npos &&
            temporary.find(L".fchk") != std::wstring::npos,
            "temporary output argument");
    require(!std::filesystem::exists(temporary), "temporary output cleaned");

    require(_wputenv_s(L"COV_FORMCHK", L"no-such-formchk-executable-for-cov-smoke.exe") == 0,
            "set missing executable");
    failure.clear();
    try { (void)cov::parse_gaussian_chk_via_formchk(checkpoint); }
    catch (const std::runtime_error& error) { failure = error.what(); }
    require(failure.find("Could not start Gaussian formchk") != std::string::npos,
            "missing executable reported");
}

} // namespace

int wmain(int argc, wchar_t* argv[]) {
    if (argc == 3) {
        const wchar_t* capture = _wgetenv(L"COV_FORMCHK_CAPTURE");
        if (!capture || !*capture) return 99;
        std::ofstream file(std::filesystem::path(capture), std::ios::binary);
        if (!file) return 99;
        write_argument(file, argv[1]);
        write_argument(file, argv[2]);
        file.close();
        if (!file) return 99;
        { std::ofstream output{std::filesystem::path(argv[2])}; output << "temporary output"; }
        return 23;
    }
    try {
        check_exit_and_literals();
        std::cout << "formchk_spawn_smoke: literal argv and exit status verified\n";
    } catch (const std::exception& error) {
        std::cerr << "formchk_spawn_smoke: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
#endif
