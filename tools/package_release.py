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
            text = re.sub(r"\(\.\./README(?:\.zh-CN)?\.md(?:#[^)]*)?\)", "(../README.md)", text)
            destination.write_text(text, encoding="utf-8")


def quickstart(platform, version):
    start = {
        "windows": "Extract the ZIP, then open `cov.exe`.",
        "macos": "Extract the ZIP, then open `Chemical Orbital Visualiser.app`. "
                 "You can move the app to Applications. This preview is signed locally "
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
                 "The app includes the native Metal module and CPU fallback. UI fonts "
                 "come from macOS; Japanese text can use Noto Sans CJK installed as "
                 "`/Library/Fonts/NotoSansCJK-Regular.ttc`.",
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
    return (
        f"# Chemical Orbital Visualiser (COV) {version}\n\n"
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
        f"User help is in `{docs}`. Library notices and the COV licence are included.\n"
    )


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
        "CFBundlePackageType": "APPL", "CFBundleShortVersionString": "0.4.0",
        "CFBundleVersion": "0.4.0d4", "COVReleaseVersion": version,
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
    actual_commit = run("git", "-C", ROOT, "rev-parse", "HEAD").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit) or actual_commit != args.commit:
        raise RuntimeError("Requested source commit does not match the checkout")
    if run("git", "-C", ROOT, "status", "--porcelain", "--untracked-files=normal").strip():
        raise RuntimeError("Release packaging requires a clean source checkout")
    if not re.fullmatch(r"v\d+\.\d+\.\d+(?:-[A-Za-z0-9.]+)?", args.version):
        raise RuntimeError("Invalid release version")
    cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
    home = re.search(r"^CMAKE_HOME_DIRECTORY:INTERNAL=(.+)$", cache, re.MULTILINE)
    if not home or Path(home.group(1).strip()).resolve() != ROOT:
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
        readme = quickstart(args.platform, args.version)
        (stage / "README.md").write_text(readme, encoding="utf-8")
        if args.platform == "macos":
            app = stage / "Chemical Orbital Visualiser.app"
            resources = app / "Contents/Resources"
            common_files(resources, build)
            (resources / "README.md").write_text(readme, encoding="utf-8")
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
