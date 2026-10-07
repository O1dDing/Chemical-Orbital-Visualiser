#include "cov/formchk.hpp"

#include "cov/fchk_overlap.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#else
#include <cerrno>
#include <cstring>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

extern "C" { extern char** environ; }
#endif

namespace cov {
namespace {

#ifdef _WIN32
std::wstring quote_windows_argument(const std::wstring& value) {
    std::wstring quoted = L"\"";
    std::size_t backslashes = 0;
    for (const wchar_t character : value) {
        if (character == L'\\') {
            ++backslashes;
            continue;
        }
        quoted.append(backslashes * (character == L'"' ? 2 : 1), L'\\');
        backslashes = 0;
        if (character == L'"') quoted += L'\\';
        quoted += character;
    }
    quoted.append(backslashes * 2, L'\\');
    quoted += L'"';
    return quoted;
}

DWORD run_formchk(const std::wstring& executable, const std::filesystem::path& input,
                  const std::filesystem::path& output) {
    // CreateProcessW receives literal arguments: cmd.exe would expand %NAME%
    // even inside quotes, changing valid checkpoint and executable paths.
    std::wstring command = quote_windows_argument(executable) + L" " +
                           quote_windows_argument(input.wstring()) + L" " +
                           quote_windows_argument(output.wstring());
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(nullptr, command.data(), nullptr, nullptr, FALSE,
                        CREATE_NO_WINDOW, nullptr, nullptr, &startup, &process)) {
        throw std::runtime_error("Could not start Gaussian formchk: Windows error " +
                                 std::to_string(GetLastError()) +
                                 ". Install Gaussian formchk or set COV_FORMCHK to its executable path.");
    }
    CloseHandle(process.hThread);
    const DWORD wait_result = WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exit_code = 0;
    const bool obtained = wait_result == WAIT_OBJECT_0 &&
                          GetExitCodeProcess(process.hProcess, &exit_code) != FALSE;
    CloseHandle(process.hProcess);
    if (!obtained)
        throw std::runtime_error("Could not obtain Gaussian formchk exit status");
    return exit_code;
}
#else
int run_formchk(const std::string& executable, const std::string& input,
                const std::string& output) {
    if (executable.find('\0') != std::string::npos ||
        input.find('\0') != std::string::npos ||
        output.find('\0') != std::string::npos)
        throw std::invalid_argument("formchk path contains a NUL byte");

    char* argv[] = {const_cast<char*>(executable.c_str()),
                    const_cast<char*>(input.c_str()),
                    const_cast<char*>(output.c_str()), nullptr};
    pid_t child = -1;
    const int launched = posix_spawnp(&child, executable.c_str(), nullptr,
                                      nullptr, argv, environ);
    if (launched != 0)
        throw std::runtime_error("Could not start Gaussian formchk: " +
                                 std::string(std::strerror(launched)) +
                                 ". Install Gaussian formchk or set COV_FORMCHK to its executable path.");

    int status = 0;
    pid_t waited;
    do { waited = waitpid(child, &status, 0); }
    while (waited == -1 && errno == EINTR);
    if (waited == -1)
        throw std::runtime_error("Could not wait for Gaussian formchk: " +
                                 std::string(std::strerror(errno)));
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    if (WIFSIGNALED(status))
        throw std::runtime_error("Gaussian formchk terminated by signal " +
                                 std::to_string(WTERMSIG(status)));
    throw std::runtime_error("Gaussian formchk did not report an exit status");
}
#endif

struct TemporaryFileGuard {
    std::filesystem::path path;
    ~TemporaryFileGuard() {
        std::error_code ec;
        if (!path.empty()) std::filesystem::remove(path, ec);
    }
};

} // namespace

Wavefunction parse_gaussian_chk_via_formchk(const std::filesystem::path& chk_path,
                                            const FchkParseOptions& options) {
    std::error_code ec;
    if (!std::filesystem::exists(chk_path, ec) || ec) {
        throw std::runtime_error("Gaussian CHK file does not exist: " + chk_path.string());
    }

#ifdef _WIN32
    // std::getenv and path(string) can lose characters outside the active
    // Windows code page; the executable path must stay wide through launch.
    const wchar_t* configured = _wgetenv(L"COV_FORMCHK");
    const std::wstring executable = configured && *configured
        ? std::wstring(configured) : L"formchk.exe";
#else
    const char* configured = std::getenv("COV_FORMCHK");
    const std::string executable = configured && *configured
        ? std::string(configured) : "formchk";
#endif

    const auto stamp = std::chrono::high_resolution_clock::now()
                           .time_since_epoch().count();
    TemporaryFileGuard output{
        std::filesystem::temp_directory_path() /
        ("cov_formchk_" + std::to_string(stamp)
#ifndef _WIN32
         + "_" + std::to_string(getpid())
#endif
         + ".fchk")
    };

#ifdef _WIN32
    const DWORD code = run_formchk(executable, chk_path, output.path);
#else
    const int code = run_formchk(executable, chk_path.string(), output.path.string());
#endif
    if (code != 0) {
        throw std::runtime_error(
            "Gaussian formchk failed with exit code " + std::to_string(code) +
            ". Install Gaussian formchk or set COV_FORMCHK to its executable path.");
    }

    if (!std::filesystem::exists(output.path, ec) || ec ||
        std::filesystem::file_size(output.path, ec) == 0 || ec) {
        throw std::runtime_error(
            "Gaussian formchk reported success but did not produce a usable FCHK file");
    }

    Wavefunction wf=parse_fchk(output.path, options);
    (void)enrich_fchk_overlap_from_file(wf, output.path);
    return wf;
}

} // namespace cov
