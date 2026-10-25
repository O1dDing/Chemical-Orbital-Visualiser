#!/usr/bin/env python3
"""Build release archives from existing, inspected binaries; never publish them.

Run from a clean checkout. The build tree must belong to that checkout, and
--commit must be the full source commit selected by the workflow.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import struct
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import zipfile


ROOT = Path(__file__).resolve().parents[1]
LABELS = {"windows": "Windows-x64", "macos": "macOS-universal", "linux": "Linux-x86_64"}
REGRESSION_URL = (
    "https://raw.githubusercontent.com/cclib/cclib/"
    "21daa960123d28aa21eeaaacc2a1dea39e136829/data/FChk/basicGaussian16/dvb_sp.fchk"
)


def run(*args):
    return subprocess.run([str(a) for a in args], check=True, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout


def copy_file(source, target):
    if not source.is_file():
        raise RuntimeError(f"Required package input is missing: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def tar_permissions(info):
    info.uid = info.gid = 0
    info.uname = info.gname = "root"
    info.mode = 0o755 if info.isdir() or info.name.endswith("/cov") else 0o644
    return info


def common_files(target, build):
    for name in ("LICENSE", "THIRD_PARTY_NOTICES.md"):
        copy_file(ROOT / name, target / name)
    shutil.copytree(ROOT / "licenses", target / "licenses")
    copy_file(build / "_deps/glfw-src/LICENSE.md", target / "licenses/glfw/LICENSE.md")
    copy_file(build / "_deps/imgui-src/LICENSE.txt", target / "licenses/imgui/LICENSE.txt")
    shutil.copytree(ROOT / "examples", target / "examples")
    # The shipped one-job example imports this small process helper.
    copy_file(ROOT / "tests/validation_process.py", target / "tests/validation_process.py")
    # Ship user help, without development archives or acceptance reports.
    for pattern in ("UI*.md", "AOMO_NBO*.md", "NBO_ONE_JOB*.md"):
        for source in sorted((ROOT / "docs").glob(pattern)):
            destination = target / "docs" / source.name
            copy_file(source, destination)
            text = destination.read_text(encoding="utf-8")
            # Each shipped language has a platform quickstart. Its section
            # headings differ from the repository README, so discard anchors.
            text = re.sub(r"\((\.\./README(?:\.(?:zh-CN|ja|fr))?\.md)#[^)]*\)", r"(\1)", text)
            destination.write_text(text, encoding="utf-8")


def quickstart(platform, version, language="en"):
    if language != "en":
        return translated_quickstart(platform, version, language)
    start = {
        "windows": "Extract the ZIP, then open `cov.exe`.",
        "macos": "Extract the ZIP, then open `Chemical Orbital Visualiser.app`. "
                 "You can move the app to Applications. This app is signed locally "
                 "and is not notarised. If macOS blocks it, use Finder's Open action "
                 "or Privacy & Security to allow the downloaded app.",
        "linux": "Extract the archive, then run `./cov` from the extracted folder.",
    }[platform]
    requirements = {
        "windows": "Windows x64 and OpenGL 2.1 or newer. Orbital rendering uses NVIDIA "
                   "CUDA when available and otherwise uses the CPU. "
                   "GPU operation needs an NVIDIA driver compatible with CUDA "
                   "12.8; installing the CUDA toolkit is optional. UI fonts come from Windows.",
        "macos": "macOS 12 or newer, on Apple silicon or Intel, with OpenGL 2.1. "
                 "The app includes the native Metal module and CPU fallback. "
                 "UI fonts come from macOS.",
        "linux": "Linux x86_64 with glibc 2.35 or newer, an X11 display (or XWayland), "
                 "and OpenGL 2.1. This CPU build is made on Ubuntu 22.04. The C++ runtime "
                 "is linked into the executable; display and system libraries come from "
                 "your distribution. On Ubuntu/Debian, install the runtime libraries "
                 "and UI fonts with:\n\n"
                 "```sh\nsudo apt-get install libgl1 libx11-6 libxrandr2 libxinerama1 "
                 "libxcursor1 libxi6 fonts-dejavu-core fonts-noto-cjk\n```",
    }[platform]
    docs = "Chemical Orbital Visualiser.app/Contents/Resources/docs" if platform == "macos" else "docs"
    examples = "Chemical Orbital Visualiser.app/Contents/Resources/examples" if platform == "macos" else "examples"
    browse_help = " In Finder, choose Show Package Contents to browse the app's included help and examples." if platform == "macos" else ""
    return (
        f"# Chemical Orbital Visualiser (COV) {version}\n\n"
        "**English** · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · [Français](README.fr.md)\n\n"
        "Explore orbital energies, occupations, energy-level diagrams, orbital "
        "composition, and three-dimensional orbital shapes. With matching NBO files, "
        "inspect charge, bonding, and connections between atomic, localised, and molecular orbitals.\n\n"
        f"{start}\n\nOpen a Gaussian FCHK/FCH or compatible Molden file, or open the "
        "calculation folder to load its wavefunction and matching NBO files together. "
        f"The `{examples}` folder includes a small hydrogen molecule. Browse levels, select "
        "an orbital, and export diagrams or data from the interface.\n\n"
        f"{requirements}\n\n"
        "Gaussian CHK conversion requires your separately installed `formchk`. "
        "Gaussian and NBO programs are not included. NBO shapes and composition "
        "require the matching report, `.47` archive, and orbital matrices.\n\n"
        f"User help is in `{docs}`.{browse_help} Library notices and the COV licence are included.\n"
    )


def translated_quickstart(platform, version, language):
    docs = "Chemical Orbital Visualiser.app/Contents/Resources/docs" if platform == "macos" else "docs"
    examples = "Chemical Orbital Visualiser.app/Contents/Resources/examples" if platform == "macos" else "examples"
    navigation = {
        "zh-CN": "[English](README.md) · **简体中文** · [日本語](README.ja.md) · [Français](README.fr.md)",
        "ja": "[English](README.md) · [简体中文](README.zh-CN.md) · **日本語** · [Français](README.fr.md)",
        "fr": "[English](README.md) · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · **Français**",
    }[language]
    introduction = {
        "zh-CN": "查看轨道能量、电子占据、能级图、轨道组成和三维轨道形状。加载匹配的 NBO 文件后，还可查看电荷、成键，以及原子轨道、局域轨道和分子轨道之间的联系。",
        "ja": "軌道エネルギー、電子占有数、エネルギー準位図、軌道組成、三次元の軌道形状を確認できます。対応する NBO ファイルがあれば、電荷、結合、原子軌道・局在軌道・分子軌道の関係も確認できます。",
        "fr": "Explorez les énergies orbitalaires, les occupations électroniques, les diagrammes de niveaux d’énergie, la composition et les formes tridimensionnelles des orbitales. Avec les fichiers NBO correspondants, examinez les charges, les liaisons et les relations entre orbitales atomiques, localisées et moléculaires.",
    }[language]
    start = {
        "zh-CN": {
            "windows": "解压 ZIP，然后打开 `cov.exe`。",
            "macos": "解压 ZIP，然后打开 `Chemical Orbital Visualiser.app`。也可将应用移到 Applications。此应用采用本地签名，未经 Apple 公证。如 macOS 阻止打开，可在 Finder 中使用“打开”，或在“隐私与安全性”中允许此应用。",
            "linux": "解压归档，然后在解压后的目录中运行 `./cov`。",
        },
        "ja": {
            "windows": "ZIP を展開し、`cov.exe` を開きます。",
            "macos": "ZIP を展開し、`Chemical Orbital Visualiser.app` を開きます。アプリケーションフォルダへ移動することもできます。このアプリはローカル署名済みで、Apple の公証は受けていません。macOS にブロックされた場合は Finder の「開く」、または「プライバシーとセキュリティ」から実行を許可してください。",
            "linux": "アーカイブを展開し、展開先のフォルダから `./cov` を実行します。",
        },
        "fr": {
            "windows": "Extrayez le ZIP, puis ouvrez `cov.exe`.",
            "macos": "Extrayez le ZIP, puis ouvrez `Chemical Orbital Visualiser.app`. Vous pouvez déplacer l’application dans Applications. Elle est signée localement et n’est pas notariée par Apple. Si macOS la bloque, utilisez l’action Ouvrir du Finder ou autorisez-la dans Confidentialité et sécurité.",
            "linux": "Extrayez l’archive, puis lancez `./cov` depuis le dossier extrait.",
        },
    }[language][platform]
    usage = {
        "zh-CN": f"打开 Gaussian FCHK/FCH 或兼容的 Molden 文件，也可打开计算目录，一起加载波函数及其匹配的 NBO 文件。`{examples}` 中包含一个小型氢分子示例。浏览能级、选择轨道，并在界面中导出能级图或数据。",
        "ja": f"Gaussian FCHK/FCH または対応する Molden ファイルを開くか、計算フォルダを開いて波動関数と対応する NBO ファイルをまとめて読み込みます。`{examples}` フォルダには小さな水素分子の例があります。準位を参照して軌道を選択し、画面から図やデータをエクスポートできます。",
        "fr": f"Ouvrez un fichier Gaussian FCHK/FCH ou Molden compatible, ou ouvrez le dossier de calcul pour charger ensemble la fonction d’onde et les fichiers NBO correspondants. Le dossier `{examples}` contient un petit exemple de molécule d’hydrogène. Parcourez les niveaux, sélectionnez une orbitale, puis exportez des diagrammes ou des données depuis l’interface.",
    }[language]
    requirements = {
        "zh-CN": {
            "windows": "需要 Windows x64 和 OpenGL 2.1 或更新版本。轨道渲染在 NVIDIA CUDA 可用时使用 CUDA，否则使用 CPU。GPU 运行需要兼容 CUDA 12.8 的 NVIDIA 驱动；无需另行安装 CUDA 工具包。界面字体使用 Windows 系统字体。",
            "macos": "需要 macOS 12 或更新版本、Apple 芯片或 Intel 处理器，以及 OpenGL 2.1。应用包含原生 Metal 模块和 CPU 回退。界面字体使用 macOS 系统字体。",
            "linux": "需要 Linux x86_64、glibc 2.35 或更新版本、X11 显示环境（或 XWayland）和 OpenGL 2.1。此 CPU 版本在 Ubuntu 22.04 上构建。C++ 运行库已链接到程序中；显示及系统库由发行版提供。在 Ubuntu/Debian 上，可使用以下命令安装运行库和界面字体：",
        },
        "ja": {
            "windows": "Windows x64 と OpenGL 2.1 以降が必要です。軌道描画には利用可能な場合に NVIDIA CUDA を使用し、それ以外は CPU を使用します。GPU での実行には CUDA 12.8 対応の NVIDIA ドライバーが必要です。CUDA ツールキットの追加インストールは不要です。画面には Windows のシステムフォントを使用します。",
            "macos": "macOS 12 以降、Apple シリコンまたは Intel、OpenGL 2.1 が必要です。アプリにはネイティブ Metal モジュールと CPU へのフォールバックが含まれます。画面には macOS のシステムフォントを使用します。",
            "linux": "Linux x86_64、glibc 2.35 以降、X11 ディスプレイ（または XWayland）、OpenGL 2.1 が必要です。この CPU ビルドは Ubuntu 22.04 で作成されています。C++ ランタイムは実行ファイルにリンクされ、ディスプレイとシステムのライブラリはディストリビューションが提供します。Ubuntu/Debian では、次のコマンドで実行時ライブラリと画面用フォントをインストールできます。",
        },
        "fr": {
            "windows": "Windows x64 et OpenGL 2.1 ou ultérieur sont nécessaires. Le rendu orbital utilise NVIDIA CUDA lorsqu’il est disponible, et le CPU sinon. L’utilisation du GPU nécessite un pilote NVIDIA compatible avec CUDA 12.8 ; l’installation du toolkit CUDA est facultative. Les polices de l’interface proviennent de Windows.",
            "macos": "macOS 12 ou ultérieur, un processeur Apple silicon ou Intel et OpenGL 2.1 sont nécessaires. L’application inclut le module Metal natif et le repli sur le CPU. Les polices de l’interface proviennent de macOS.",
            "linux": "Linux x86_64 avec glibc 2.35 ou ultérieur, un affichage X11 (ou XWayland) et OpenGL 2.1 sont nécessaires. Cette version CPU est construite sur Ubuntu 22.04. Le runtime C++ est lié à l’exécutable ; les bibliothèques d’affichage et système proviennent de votre distribution. Sur Ubuntu/Debian, installez les bibliothèques et les polices avec :",
        },
    }[language][platform]
    if platform == "linux":
        requirements += "\n\n```sh\nsudo apt-get install libgl1 libx11-6 libxrandr2 libxinerama1 libxcursor1 libxi6 fonts-dejavu-core fonts-noto-cjk\n```"
    dependencies = {
        "zh-CN": "Gaussian CHK 转换需要您另行安装的 `formchk`。本包不含 Gaussian 和 NBO 程序。NBO 形状和组成需要匹配的报告、`.47` 归档和轨道矩阵。",
        "ja": "Gaussian CHK の変換には、別途インストールした `formchk` が必要です。Gaussian と NBO のプログラムは含まれていません。NBO の形状と組成には、対応するレポート、`.47` アーカイブ、軌道行列が必要です。",
        "fr": "La conversion des fichiers Gaussian CHK nécessite votre propre installation de `formchk`. Les programmes Gaussian et NBO ne sont pas inclus. Les formes et la composition NBO nécessitent le rapport, l’archive `.47` et les matrices orbitalaires correspondants.",
    }[language]
    help_text = {
        "zh-CN": f"用户帮助位于 `{docs}`。随包附有库声明和 COV 许可证。",
        "ja": f"ユーザーヘルプは `{docs}` にあります。ライブラリの告知と COV のライセンスも同梱しています。",
        "fr": f"L’aide se trouve dans `{docs}`. Les notices des bibliothèques et la licence de COV sont incluses.",
    }[language]
    if platform == "macos":
        help_text += {
            "zh-CN": " 在 Finder 中选择“显示包内容”，即可浏览随应用附带的帮助和示例。",
            "ja": " Finder で「パッケージの内容を表示」を選ぶと、同梱のヘルプと例を参照できます。",
            "fr": " Dans le Finder, choisissez Afficher le contenu du paquet pour parcourir l’aide et les exemples inclus.",
        }[language]
    return f"# Chemical Orbital Visualiser (COV) {version}\n\n{navigation}\n\n{introduction}\n\n{start}\n\n{usage}\n\n{requirements}\n\n{dependencies}\n\n{help_text}\n"


def write_quickstarts(target, platform, version):
    for language, suffix in (("en", ""), ("zh-CN", ".zh-CN"), ("ja", ".ja"), ("fr", ".fr")):
        (target / f"README{suffix}.md").write_text(quickstart(platform, version, language), encoding="utf-8")


def numeric_version(version):
    match = re.fullmatch(r"v(\d+\.\d+\.\d+)(?:-[A-Za-z0-9.]+)?", version)
    if not match:
        raise RuntimeError("Invalid release version")
    return match.group(1)


def inspect_windows(binary, audit):
    with binary.open("rb") as stream:
        if stream.read(2) != b"MZ":
            raise RuntimeError("Not a Windows PE executable")
        stream.seek(0x3C)
        pe_offset = struct.unpack("<I", stream.read(4))[0]
        stream.seek(pe_offset)
        if stream.read(6) != b"PE\0\0\x64\x86":
            raise RuntimeError("Windows executable is not x64")
    vswhere = Path(os.environ["ProgramFiles(x86)"]) / "Microsoft Visual Studio/Installer/vswhere.exe"
    installation = run(vswhere, "-latest", "-products", "*", "-requires",
                       "Microsoft.VisualStudio.Component.VC.Tools.x86.x64", "-property", "installationPath").strip()
    candidates = sorted(Path(installation).glob("VC/Tools/MSVC/*/bin/Hostx64/x64/dumpbin.exe"))
    if not candidates:
        raise RuntimeError("MSVC dependency inspector was not found")
    dependencies = run(candidates[-1], "/DEPENDENTS", binary)
    (audit / "windows-dependencies.txt").write_text(dependencies, encoding="utf-8")
    system_dlls = {
        "kernel32.dll", "user32.dll", "gdi32.dll", "shell32.dll", "comdlg32.dll",
        "advapi32.dll", "ole32.dll", "oleaut32.dll", "opengl32.dll", "ws2_32.dll",
        "bcrypt.dll", "crypt32.dll", "version.dll", "shlwapi.dll", "winmm.dll",
        "imm32.dll", "dwmapi.dll", "setupapi.dll", "ntdll.dll", "secur32.dll",
        "comctl32.dll", "ucrtbase.dll",
    }
    dlls = re.findall(r"^\s+([\w.-]+\.dll)\s*$", dependencies, re.MULTILINE | re.IGNORECASE)
    if not dlls:
        raise RuntimeError("Could not read Windows startup DLL dependencies")
    unknown = [name for name in dlls if name.lower() not in system_dlls
               and not name.lower().startswith(("api-ms-win-", "ext-ms-win-"))]
    if unknown:
        raise RuntimeError(f"Unpackaged Windows startup dependencies: {unknown}")
    cuda = Path(os.environ["CUDA_PATH"]) / "bin/cuobjdump.exe"
    cubins = run(cuda, "--list-elf", binary)
    (audit / "cuda-architectures.txt").write_text(cubins, encoding="utf-8")
    for architecture in (60, 61, 70, 75, 80, 86, 89, 90, 100, 120):
        if not re.search(rf"\bsm_{architecture}\b", cubins):
            raise RuntimeError(f"CUDA sm_{architecture} cubin is missing")
    ptx = run(cuda, "--dump-ptx", binary)
    if not re.search(r"\.target\s+sm_120\b", ptx):
        raise RuntimeError("CUDA compute_120 PTX is missing")
    (audit / "cuda-ptx-target.txt").write_text("compute_120 PTX present\n", encoding="utf-8")


def inspect_linux(binary, audit):
    header = binary.read_bytes()[:20]
    if header[:6] != b"\x7fELF\x02\x01" or struct.unpack("<H", header[18:20])[0] != 62:
        raise RuntimeError("Linux executable must be little-endian ELF x86_64")
    dependencies = run("ldd", binary)
    (audit / "linux-dependencies.txt").write_text(dependencies, encoding="utf-8")
    if "not found" in dependencies or re.search(r"lib(?:stdc\+\+|gcc_s)\.", dependencies):
        raise RuntimeError("Missing Linux library or C++ runtime was not linked statically")
    # Keep the documented baseline honest, even if runner compiler defaults change.
    symbols = run("readelf", "--version-info", binary)
    versions = [tuple(map(int, v.split("."))) for v in re.findall(r"\bGLIBC_(\d+(?:\.\d+)+)", symbols)]
    if not versions or max(versions) > (2, 35):
        raise RuntimeError("Linux glibc requirement exceeds the Ubuntu 22.04 baseline")
    (audit / "linux-symbol-versions.txt").write_text(symbols, encoding="utf-8")
    # No fetched/build-tree shared libraries may accidentally escape packaging.
    for line in dependencies.splitlines():
        match = re.search(r"=>\s+(\S+)", line)
        if match and not match.group(1).startswith(("/lib/", "/lib64/", "/usr/lib/", "/usr/lib64/")):
            raise RuntimeError(f"Linux executable depends on a non-system library: {line}")


def runtime_notices(platform, target):
    if platform == "linux":
        copy_file(Path("/usr/share/doc/libstdc++6/copyright"), target / "licenses/gcc/copyright")
        copy_file(Path("/usr/share/common-licenses/GPL-3"), target / "licenses/gcc/GPL-3")
    elif platform == "windows":
        toolkit = Path(os.environ["CUDA_PATH"])
        candidates = [toolkit / "EULA.txt", toolkit / "doc/EULA.txt", toolkit / "LICENSE.txt"]
        notice = next((source for source in candidates if source.is_file()), None)
        if notice:
            copy_file(notice, target / "licenses/nvidia/CUDA-LICENSE.txt")
        else:
            # Partial toolkit installs may omit documentation. Retain the
            # versioned NVIDIA notice for the statically linked CUDA runtime.
            url = "https://docs.nvidia.com/cuda/archive/12.8.1/eula/index.html"
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            if b"License Agreement" not in data or b"cudart_static.lib" not in data:
                raise RuntimeError("CUDA 12.8 license notice download is invalid")
            location = target / "licenses/nvidia/CUDA-12.8-EULA.html"
            location.parent.mkdir(parents=True, exist_ok=True)
            location.write_bytes(data)


def prepare_macos(app, audit, version):
    macos = app / "Contents/MacOS"
    viewer = macos / "cov"
    module = macos / "libcov_compute_metal.dylib"
    for binary in (viewer, module):
        architectures = run("lipo", "-archs", binary).split()
        if set(architectures) != {"arm64", "x86_64"}:
            raise RuntimeError(f"Not a universal arm64+x86_64 binary: {binary}: {architectures}")
        load_commands = run("otool", "-l", binary)
        minimum_versions = re.findall(r"\bminos (\d+(?:\.\d+)+)|\bversion (\d+(?:\.\d+)+)\s+sdk", load_commands)
        if not minimum_versions or any(tuple(map(int, (modern or legacy).split(".")[:2])) > (12, 0)
                                       for modern, legacy in minimum_versions):
            raise RuntimeError(f"macOS deployment target exceeds 12.0: {binary}")
        rpaths = re.findall(r"cmd LC_RPATH\s+cmdsize \d+\s+path (.*?) \(offset", load_commands)
        for rpath in dict.fromkeys(rpaths):
            if rpath != "@loader_path":
                run("install_name_tool", "-delete_rpath", rpath, binary)
        if "@loader_path" not in rpaths:
            run("install_name_tool", "-add_rpath", "@loader_path", binary)
    run("install_name_tool", "-id", "@rpath/libcov_compute_metal.dylib", module)
    for binary in (viewer, module):
        dependencies = run("otool", "-L", binary)
        (audit / f"macos-{binary.name}-dependencies.txt").write_text(dependencies, encoding="utf-8")
        for line in dependencies.splitlines():
            if " (compatibility version" not in line:
                continue
            name = line.strip().split(" (compatibility version", 1)[0]
            if name == "@rpath/libcov_compute_metal.dylib" and binary == module:
                continue  # The dylib's own install name.
            if not name.startswith(("/usr/lib/", "/System/Library/Frameworks/")):
                raise RuntimeError(f"Unpackaged macOS dependency: {name}")
        (audit / f"macos-{binary.name}-architectures.txt").write_text(
            run("lipo", "-archs", binary), encoding="utf-8")
    info = {
        "CFBundleName": "Chemical Orbital Visualiser", "CFBundleDisplayName": "Chemical Orbital Visualiser",
        "CFBundleIdentifier": "org.chemicalorbitalvisualiser.cov", "CFBundleExecutable": "cov",
        "CFBundlePackageType": "APPL", "CFBundleShortVersionString": numeric_version(version),
        "CFBundleVersion": numeric_version(version), "COVReleaseVersion": version,
        "LSMinimumSystemVersion": "12.0", "NSHighResolutionCapable": True,
    }
    with (app / "Contents/Info.plist").open("wb") as stream:
        plistlib.dump(info, stream)
    # Rpath edits invalidate the compiler's signatures. Sign inside out, then
    # check the complete app. This is ad-hoc signing, not notarisation.
    run("codesign", "--force", "--sign", "-", module)
    run("codesign", "--force", "--sign", "-", app)
    run("codesign", "--verify", "--deep", "--strict", "--verbose=2", app)


def package(args):
    build = args.build.resolve()
    source_root = ROOT.resolve(strict=True)
    actual_commit = run("git", "-C", source_root, "rev-parse", "HEAD").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit) or actual_commit != args.commit:
        raise RuntimeError("Requested source commit does not match the checkout")
    if run("git", "-C", source_root, "status", "--porcelain", "--untracked-files=normal").strip():
        raise RuntimeError("Release packaging requires a clean source checkout")
    base_version = numeric_version(args.version)
    cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
    if not re.search(rf"^CMAKE_PROJECT_VERSION:STATIC={re.escape(base_version)}$", cache, re.MULTILINE):
        raise RuntimeError("Release version does not match the configured project version")
    home = re.search(r"^CMAKE_HOME_DIRECTORY:INTERNAL=(.+)$", cache, re.MULTILINE)
    # Resolve both paths: macOS /var aliases /private/var, and Windows runner
    # temporary paths may use short names. Comparing only one resolved side
    # rejects a build that belongs to the very same directory.
    if not home or Path(home.group(1).strip()).resolve(strict=True) != source_root:
        raise RuntimeError("Build tree does not belong to this source checkout")
    expected = {"COV_BUILD_VIEWER": "ON", "COV_ENABLE_VALIDATION": "OFF",
                "COV_ENABLE_CUDA": "ON" if args.platform == "windows" else "OFF"}
    if args.platform == "macos":
        expected["COV_ENABLE_METAL"] = "ON"
    for key, value in expected.items():
        if not re.search(rf"^{key}:BOOL={value}$", cache, re.MULTILINE):
            raise RuntimeError(f"Unexpected build configuration: {key} must be {value}")
    destination = args.output.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    name = f"Chemical-Orbital-Visualiser-{args.version}-{LABELS[args.platform]}"
    archive = destination / (name + (".tar.gz" if args.platform == "linux" else ".zip"))
    if archive.exists():
        raise RuntimeError(f"Refusing to overwrite an existing package: {archive}")
    audit = build / "release-audit"
    audit.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cov-package-") as temporary:
        stage = Path(temporary) / name
        stage.mkdir()
        write_quickstarts(stage, args.platform, args.version)
        if args.platform == "macos":
            app = stage / "Chemical Orbital Visualiser.app"
            resources = app / "Contents/Resources"
            common_files(resources, build)
            write_quickstarts(resources, args.platform, args.version)
            copy_file(build / "cov", app / "Contents/MacOS/cov")
            copy_file(build / "libcov_compute_metal.dylib", app / "Contents/MacOS/libcov_compute_metal.dylib")
            prepare_macos(app, audit, args.version)
            run("ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", stage, archive)
        else:
            common_files(stage, build)
            runtime_notices(args.platform, stage)
            binary = stage / ("cov.exe" if args.platform == "windows" else "cov")
            copy_file(build / "Release/cov.exe" if args.platform == "windows" else build / "cov", binary)
            if args.platform == "windows":
                inspect_windows(binary, audit)
                with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as output:
                    for file in sorted(stage.rglob("*")):
                        if file.is_file():
                            output.write(file, file.relative_to(stage.parent))
            else:
                binary.chmod(0o755)
                inspect_linux(binary, audit)
                with tarfile.open(archive, "w:gz") as output:
                    output.add(stage, arcname=name, filter=tar_permissions)
    digest = hashlib.sha256()
    with archive.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    record = {"source_commit": actual_commit, "version": args.version, "platform": args.platform,
              "archive": archive.name, "sha256": digest.hexdigest()}
    (audit / "package.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Created {archive}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("fetch-regression")
    fetch.add_argument("--build", type=Path, required=True)
    release = commands.add_parser("package")
    release.add_argument("--platform", choices=LABELS, required=True)
    release.add_argument("--build", type=Path, required=True)
    release.add_argument("--version", required=True)
    release.add_argument("--commit", required=True)
    release.add_argument("--output", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    if args.command == "fetch-regression":
        args.build.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(REGRESSION_URL, timeout=60) as response:
            data = response.read()
        if b"Number of atoms" not in data or b"Number of basis functions" not in data:
            raise RuntimeError("Pinned regression download is not an FCHK input")
        (args.build / "real_gaussian16_dvb_sp.fchk").write_bytes(data)
    else:
        package(args)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        if isinstance(error, subprocess.CalledProcessError):
            print(error.stdout or "", file=sys.stderr)
        print(f"Release packaging failed: {error}", file=sys.stderr)
        sys.exit(1)
