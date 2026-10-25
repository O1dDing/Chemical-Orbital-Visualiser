"""Archive/provenance checks with synthetic inputs; platform tools run in CI."""

import argparse
import importlib.util
import json
from pathlib import Path
import plistlib
import struct
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import zipfile


spec = importlib.util.spec_from_file_location("package_release", Path(__file__).parents[1] / "package_release.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)
COMMIT = "a" * 40


class PackageChecks(unittest.TestCase):
    def fixture(self, directory, platform):
        root = Path(directory) / "source with spaces"
        build = root / "build-release"
        for filename in ("LICENSE", "THIRD_PARTY_NOTICES.md", "licenses/eigen/COPYING",
                         "examples/h2.molden", "tests/validation_process.py", "docs/UI.md",
                         "docs/AOMO_NBO.md", "docs/NBO_ONE_JOB.md", "docs/UI.zh-CN.md",
                         "docs/UI.ja.md", "docs/UI.fr.md",
                         "build-release/_deps/glfw-src/LICENSE.md", "build-release/_deps/imgui-src/LICENSE.txt"):
            path = root / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("example\n", encoding="utf-8")
        (root / "docs/acceptance-report.md").write_text("not for release")
        (root / "docs/UI.md").write_text("[quickstart](../README.zh-CN.md#下载)", encoding="utf-8")
        for language in ("zh-CN", "ja", "fr"):
            (root / f"docs/UI.{language}.md").write_text(
                f"[quickstart](../README.{language}.md#start)", encoding="utf-8")
        cuda = "ON" if platform == "windows" else "OFF"
        (build / "CMakeCache.txt").write_text(
            f"CMAKE_HOME_DIRECTORY:INTERNAL={root}\nCMAKE_PROJECT_VERSION:STATIC=0.4.0\nCOV_BUILD_VIEWER:BOOL=ON\n"
            f"COV_ENABLE_VALIDATION:BOOL=OFF\nCOV_ENABLE_CUDA:BOOL={cuda}\n"
            f"COV_ENABLE_METAL:BOOL={'ON' if platform == 'macos' else 'OFF'}\n", encoding="utf-8")
        binary = build / ("Release/cov.exe" if platform == "windows" else "cov")
        binary.parent.mkdir(exist_ok=True)
        binary.write_bytes(b"synthetic executable")
        if platform == "macos":
            (build / "libcov_compute_metal.dylib").write_bytes(b"synthetic Metal module")
        args = argparse.Namespace(build=build, platform=platform, commit=COMMIT,
                                  version="v0.4.0", output=root / "dist")
        return root, args

    def checkout(self, *args):
        if "rev-parse" in args:
            return COMMIT + "\n"
        if "status" in args:
            return ""
        raise AssertionError(args)

    def test_windows_archive_and_source_record(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "windows")
            with patch.object(release, "ROOT", root), patch.object(release, "run", self.checkout), \
                    patch.object(release, "inspect_windows"), patch.object(release, "runtime_notices"):
                release.package(args)
                archive = args.output / "Chemical-Orbital-Visualiser-v0.4.0-Windows-x64.zip"
                with zipfile.ZipFile(archive) as bundle:
                    members = bundle.namelist()
                    prefix = archive.stem + "/"
                    self.assertIn(prefix + "cov.exe", members)
                    self.assertEqual(bundle.read(prefix + "cov.exe"), b"synthetic executable")
                    self.assertIn(prefix + "licenses/glfw/LICENSE.md", members)
                    self.assertIn(prefix + "tests/validation_process.py", members)
                    self.assertIn(prefix + "examples/h2.molden", members)
                    self.assertFalse(any("acceptance-report" in name for name in members))
                    self.assertIn("../README.zh-CN.md", bundle.read(prefix + "docs/UI.md").decode())
                    for suffix in ("", ".zh-CN", ".ja", ".fr"):
                        text = bundle.read(prefix + f"README{suffix}.md").decode("utf-8")
                        self.assertIn("v0.4.0", text)
                        self.assertIn("cov.exe", text)
                        self.assertIn("CUDA 12.8", text)
                    for language in ("zh-CN", "ja", "fr"):
                        text = bundle.read(prefix + f"docs/UI.{language}.md").decode("utf-8")
                        self.assertEqual(text, f"[quickstart](../README.{language}.md)")
                record = json.loads((args.build / "release-audit/package.json").read_text())
                self.assertEqual(record["source_commit"], COMMIT)
                self.assertEqual(len(record["sha256"]), 64)
                with self.assertRaisesRegex(RuntimeError, "overwrite"):
                    release.package(args)

    def test_linux_tar_preserves_executable_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            with patch.object(release, "ROOT", root), patch.object(release, "run", self.checkout), \
                    patch.object(release, "inspect_linux"), patch.object(release, "runtime_notices"):
                release.package(args)
            archive = args.output / "Chemical-Orbital-Visualiser-v0.4.0-Linux-x86_64.tar.gz"
            with tarfile.open(archive) as bundle:
                binary = bundle.getmember(archive.name.removesuffix(".tar.gz") + "/cov")
                self.assertEqual(binary.mode & 0o777, 0o755)

    def test_macos_archive_includes_localised_help_and_resources(self):
        def archive_tools(*command):
            if command[0] != "ditto":
                return self.checkout(*command)
            stage, archive = command[-2:]
            with zipfile.ZipFile(archive, "w") as output:
                for file in stage.rglob("*"):
                    if file.is_file():
                        output.write(file, file.relative_to(stage.parent))
            return ""

        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "macos")
            with patch.object(release, "ROOT", root), patch.object(release, "run", archive_tools), \
                    patch.object(release, "prepare_macos"):
                release.package(args)
            archive = args.output / "Chemical-Orbital-Visualiser-v0.4.0-macOS-universal.zip"
            with zipfile.ZipFile(archive) as bundle:
                prefix = archive.stem + "/"
                app = prefix + "Chemical Orbital Visualiser.app/Contents/"
                self.assertEqual(bundle.read(app + "MacOS/cov"), b"synthetic executable")
                self.assertEqual(bundle.read(app + "MacOS/libcov_compute_metal.dylib"), b"synthetic Metal module")
                for suffix in ("", ".zh-CN", ".ja", ".fr"):
                    filename = f"README{suffix}.md"
                    self.assertEqual(bundle.read(prefix + filename), bundle.read(app + "Resources/" + filename))
                    self.assertIn(app + f"Resources/docs/UI{suffix}.md", bundle.namelist())

    def test_mismatched_project_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            args.version = "v0.3.0"
            with patch.object(release, "ROOT", root), patch.object(release, "run", self.checkout):
                with self.assertRaisesRegex(RuntimeError, "project version"):
                    release.package(args)
            self.assertFalse(args.output.exists())

    def test_localised_quickstarts_keep_platform_requirements(self):
        for language in ("en", "zh-CN", "ja", "fr"):
            for platform in ("windows", "macos", "linux"):
                with self.subTest(language=language, platform=platform):
                    text = release.quickstart(platform, "v0.4.0", language)
                    self.assertIn("OpenGL 2.1", text)
                    self.assertIn("formchk", text)
                    self.assertIn(".47", text)
                    self.assertNotIn("preview", text)
                    for suffix in ("", ".zh-CN", ".ja", ".fr"):
                        if suffix != ("" if language == "en" else f".{language}"):
                            self.assertIn(f"README{suffix}.md", text)
                    if platform == "linux":
                        self.assertIn("2.35", text)
                        self.assertIn("sudo apt-get install", text)
                    elif platform == "macos":
                        self.assertIn("Metal", text)
                        self.assertIn("Contents/Resources/docs", text)

    def test_commit_mismatch_stops_before_packaging(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            args.commit = "b" * 40
            with patch.object(release, "ROOT", root), patch.object(release, "run", self.checkout):
                with self.assertRaisesRegex(RuntimeError, "commit"):
                    release.package(args)
            self.assertFalse(args.output.exists())

    def test_source_directory_alias_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            alias = root / ".." / root.name
            self.assertNotEqual(alias, root.resolve())
            with patch.object(release, "ROOT", alias), patch.object(release, "run", self.checkout), \
                    patch.object(release, "inspect_linux"), patch.object(release, "runtime_notices"):
                release.package(args)
            self.assertTrue(next(args.output.glob("*.tar.gz")).is_file())

    def test_different_source_directory_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            other = root.parent / "different checkout"
            other.mkdir()
            cache = args.build / "CMakeCache.txt"
            cache.write_text(cache.read_text(encoding="utf-8").replace(str(root), str(other)), encoding="utf-8")
            with patch.object(release, "ROOT", root), patch.object(release, "run", self.checkout):
                with self.assertRaisesRegex(RuntimeError, "does not belong"):
                    release.package(args)
            self.assertFalse(args.output.exists())

    def test_linux_rejects_runtime_or_baseline_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "cov"
            header = bytearray(20)
            header[:6] = b"\x7fELF\x02\x01"
            header[18:20] = struct.pack("<H", 62)
            binary.write_bytes(header)
            audit = Path(directory)
            for dependency in ("libfoo.so => not found", "libstdc++.so.6 => /lib/libstdc++.so.6"):
                with patch.object(release, "run", return_value=dependency):
                    with self.assertRaises(RuntimeError):
                        release.inspect_linux(binary, audit)
            with patch.object(release, "run", side_effect=["libc.so.6 => /lib/libc.so.6", "GLIBC_2.36"]):
                with self.assertRaisesRegex(RuntimeError, "baseline"):
                    release.inspect_linux(binary, audit)

    def test_macos_rejects_missing_architecture(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(release, "run", return_value="arm64"):
                with self.assertRaisesRegex(RuntimeError, "universal"):
                    release.prepare_macos(Path(directory) / "COV.app", Path(directory), "v0.4.0")

    def test_macos_stable_bundle_metadata(self):
        def native_tools(*args):
            if args[0] == "lipo":
                return "arm64 x86_64"
            if args[0] == "otool" and args[1] == "-l":
                return "cmd LC_BUILD_VERSION\nminos 12.0\ncmd LC_RPATH\ncmdsize 40\npath @loader_path (offset 12)"
            if args[0] == "otool" and args[1] == "-L":
                return "\t/usr/lib/libSystem.B.dylib (compatibility version 1.0.0)"
            return ""

        with tempfile.TemporaryDirectory() as directory:
            app = Path(directory) / "COV.app"
            (app / "Contents/MacOS").mkdir(parents=True)
            with patch.object(release, "run", native_tools):
                release.prepare_macos(app, Path(directory), "v0.4.0")
            with (app / "Contents/Info.plist").open("rb") as stream:
                info = plistlib.load(stream)
            self.assertEqual(info["CFBundleShortVersionString"], "0.4.0")
            self.assertEqual(info["CFBundleVersion"], "0.4.0")
            self.assertEqual(info["COVReleaseVersion"], "v0.4.0")

    def test_dirty_checkout_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, args = self.fixture(directory, "linux")
            def dirty_checkout(*command):
                return " M src/main.cpp\n" if "status" in command else self.checkout(*command)
            with patch.object(release, "ROOT", root), patch.object(release, "run", dirty_checkout):
                with self.assertRaisesRegex(RuntimeError, "clean source checkout"):
                    release.package(args)
            self.assertFalse(args.output.exists())


if __name__ == "__main__":
    unittest.main()
