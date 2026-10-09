#include "cov/file_dialog.hpp"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#include <commdlg.h>

#include <array>
#include <sstream>

namespace cov {

FileDialogResult open_wavefunction_file_dialog(ui::Language language, bool nbo_input) {
    FileDialogResult result;

    std::array<wchar_t, 32768> buffer{};
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFile = buffer.data();
    ofn.nMaxFile = static_cast<DWORD>(buffer.size());
    const auto local = [&](const wchar_t* en, const wchar_t* zh,
                           const wchar_t* ja, const wchar_t* fr) {
        switch(language) {
            case ui::Language::ChineseSimplified: return zh;
            case ui::Language::Japanese: return ja;
            case ui::Language::French: return fr;
            default: return en;
        }
    };
    std::wstring filters;
    const auto add_filter = [&](const wchar_t* label, const wchar_t* pattern) {
        filters += label; filters.push_back(L'\0');
        filters += pattern; filters.push_back(L'\0');
    };
    if(nbo_input) {
        add_filter(local(L"NBO data", L"NBO 数据", L"NBO データ", L"Données NBO"),
                   L"*.covnbopkg;*.log;*.out;*.nbo;*.47");
    } else {
        add_filter(local(L"Calculation files", L"计算文件", L"計算ファイル", L"Fichiers de calcul"),
                   L"*.fchk;*.fch;*.chk;*.molden;*.molden.input;*.molden.inp;*.covnbopkg");
    }
    add_filter(L"Gaussian FCHK / FCH", L"*.fchk;*.fch");
    add_filter(L"Gaussian CHK (formchk)", L"*.chk");
    add_filter(L"Molden", L"*.molden;*.molden.input;*.molden.inp");
    add_filter(local(L"All files", L"所有文件", L"すべてのファイル", L"Tous les fichiers"), L"*.*");
    filters.push_back(L'\0');
    ofn.lpstrFilter = filters.c_str();
    ofn.nFilterIndex = 1;
    ofn.lpstrTitle = nbo_input
        ? local(L"Open NBO data", L"打开 NBO 数据", L"NBO データを開く", L"Ouvrir des données NBO")
        : local(L"Open calculation", L"打开计算文件", L"計算ファイルを開く", L"Ouvrir un calcul");
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST |
                OFN_NOCHANGEDIR | OFN_DONTADDTORECENT;

    if (GetOpenFileNameW(&ofn) != FALSE) {
        result.path = std::filesystem::path(buffer.data());
        return result;
    }

    const DWORD error = CommDlgExtendedError();
    if (error == 0) {
        result.cancelled = true;
        return result;
    }

    std::ostringstream message;
    message << error;
    const char* detail = language == ui::Language::ChineseSimplified ? "无法打开文件窗口。错误码：" :
        language == ui::Language::Japanese ? "ファイル選択を開けません。エラー：" :
        language == ui::Language::French ? "La fenêtre de sélection ne s’ouvre pas. Code : " :
        "The file picker could not open. Error: ";
    result.error = detail + message.str();
    return result;
}

} // namespace cov

#else

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

extern "C" { extern char** environ; }

namespace cov {
namespace {

struct DialogRun {
    bool available = true;
    int exit_code = -1;
    std::string output;
    std::string error;
};

DialogRun run_dialog(const std::vector<std::string>& arguments) {
    DialogRun result;
    int descriptors[2];
    if (pipe(descriptors) != 0) {
        result.error = "Could not create file dialog pipe: " + std::string(std::strerror(errno));
        return result;
    }

    posix_spawn_file_actions_t actions;
    const int init_error = posix_spawn_file_actions_init(&actions);
    int action_error = init_error;
    if (action_error == 0)
        action_error = posix_spawn_file_actions_addclose(&actions,descriptors[0]);
    if (action_error == 0)
        action_error = posix_spawn_file_actions_adddup2(&actions,descriptors[1],STDOUT_FILENO);
    if (action_error == 0 && descriptors[1] != STDOUT_FILENO)
        action_error = posix_spawn_file_actions_addclose(&actions,descriptors[1]);
    if (action_error != 0) {
        close(descriptors[0]);
        close(descriptors[1]);
        if (init_error == 0) posix_spawn_file_actions_destroy(&actions);
        result.error = "Could not configure file dialog: " + std::string(std::strerror(action_error));
        return result;
    }

    std::vector<char*> argv;
    argv.reserve(arguments.size()+1);
    for (const auto& argument : arguments)
        argv.push_back(const_cast<char*>(argument.c_str()));
    argv.push_back(nullptr);
    pid_t child = -1;
    const int launched = posix_spawnp(&child,arguments[0].c_str(),&actions,
                                      nullptr,argv.data(),environ);
    posix_spawn_file_actions_destroy(&actions);
    close(descriptors[1]);
    if (launched != 0) {
        close(descriptors[0]);
        if (launched == ENOENT) result.available = false;
        else result.error = "Could not start file dialog: " + std::string(std::strerror(launched));
        return result;
    }

    std::array<char,4096> buffer{};
    bool too_long = false;
    for (;;) {
        const ssize_t received = read(descriptors[0],buffer.data(),buffer.size());
        if (received == 0) break;
        if (received < 0) {
            if (errno == EINTR) continue;
            result.error = "Could not read file dialog selection: " + std::string(std::strerror(errno));
            break;
        }
        if (result.output.size()+static_cast<std::size_t>(received) > 65536u)
            too_long = true;
        else if (!too_long)
            result.output.append(buffer.data(),static_cast<std::size_t>(received));
    }
    close(descriptors[0]);
    int status = 0;
    pid_t waited;
    do { waited = waitpid(child,&status,0); }
    while (waited == -1 && errno == EINTR);
    if (waited == -1)
        result.error = "Could not wait for file dialog: " + std::string(std::strerror(errno));
    else if (WIFEXITED(status)) result.exit_code = WEXITSTATUS(status);
    else result.error = "File dialog terminated unexpectedly";
    if (too_long) result.error = "File dialog returned an oversized path";
    return result;
}

} // namespace

FileDialogResult open_wavefunction_file_dialog(ui::Language, bool) {
    FileDialogResult result;
#if defined(__APPLE__)
    const std::vector<std::vector<std::string>> commands{{
        "osascript", "-e",
        "POSIX path of (choose file with prompt \"Open wavefunction\")"
    }};
#else
    if (!std::getenv("DISPLAY") && !std::getenv("WAYLAND_DISPLAY")) {
        result.supported = false;
        return result;
    }
    const std::vector<std::vector<std::string>> commands{
        {"zenity","--file-selection","--title=Open wavefunction"},
        {"kdialog","--getopenfilename",".",
         "Wavefunctions (*.fchk *.fch *.chk *.molden *.molden.input *.molden.inp);;All files (*)"}
    };
#endif
    for (const auto& command : commands) {
        auto launched = run_dialog(command);
        if (!launched.available) continue;
        if (!launched.error.empty()) {
            result.error = std::move(launched.error);
            return result;
        }
        if (launched.exit_code == 1) {
            result.cancelled = true;
            return result;
        }
        if (launched.exit_code != 0) {
            result.error = "File dialog exited with status " + std::to_string(launched.exit_code);
            return result;
        }
        if (!launched.output.empty() && launched.output.back() == '\n')
            launched.output.pop_back();
        if (!launched.output.empty() && launched.output.back() == '\r')
            launched.output.pop_back();
        if (launched.output.empty()) {
            result.error = "File dialog returned an empty path";
            return result;
        }
        result.path = std::filesystem::path(std::move(launched.output));
        return result;
    }
    result.supported = false;
    return result;
}

} // namespace cov

#endif
