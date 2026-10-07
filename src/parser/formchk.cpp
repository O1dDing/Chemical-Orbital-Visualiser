#include "cov/formchk.hpp"

#include "cov/fchk_overlap.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>

#ifndef _WIN32
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
std::string quoted_shell_argument(const std::string& value, const char* label) {
    if (value.find('"') != std::string::npos ||
        value.find('\r') != std::string::npos ||
        value.find('\n') != std::string::npos) {
        throw std::runtime_error(std::string(label) +
                                 " contains characters unsafe for formchk invocation");
    }
    return '"' + value + '"';
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

    std::string executable;
    if (const char* configured = std::getenv("COV_FORMCHK"); configured && *configured) {
        executable = configured;
    } else {
#ifdef _WIN32
        executable = "formchk.exe";
#else
        executable = "formchk";
#endif
    }

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
    const std::string command =
        quoted_shell_argument(executable, "formchk executable") + " " +
        quoted_shell_argument(chk_path.string(), "CHK path") + " " +
        quoted_shell_argument(output.path.string(), "temporary FCHK path");

    const int code = std::system(command.c_str());
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
