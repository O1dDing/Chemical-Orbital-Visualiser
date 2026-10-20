"""Archive/provenance checks with synthetic inputs; platform tools run in CI."""

import argparse
import importlib.util
import json
from pathlib import Path
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
                         "docs/AOMO_NBO.md", "docs/NBO_ONE_JOB.md",
                         "build-release/_deps/glfw-src/LICENSE.md", "build-release/_deps/imgui-src/LICENSE.txt"):
            path = root / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("example\n", encoding="utf-8")
        (root / "docs/acceptance-report.md").write_text("not for release")
        (root / "docs/UI.md").write_text("[quickstart](../README.zh-CN.md#下载)", encoding="utf-8")
        cuda = "ON" if platform == "windows" else "OFF"
        (build / "CMakeCache.txt").write_text(
            f"CMAKE_HOME_DIRECTORY:INTERNAL={root}\nCOV_BUILD_VIEWER:BOOL=ON\n"
            f"COV_ENABLE_VALIDATION:BOOL=OFF\nCOV_ENABLE_CUDA:BOOL={cuda}\n", encoding="utf-8")
        binary = build / ("Release/cov.exe" if platform == "windows" else "cov")
        binary.parent.mkdir(exist_ok=True)
        binary.write_bytes(b"synthetic executable")
        args = argparse.Namespace(build=build, platform=platform, commit=COMMIT,
                                  version="v0.4.0-pre.4", output=root / "dist")
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
                archive = args.output / "Chemical-Orbital-Visualiser-v0.4.0-pre.4-Windows-x64.zip"
                with zipfile.ZipFile(archive) as bundle:
                    members = bundle.namelist()
                    prefix = archive.stem + "/"
                    self.assertIn(prefix + "cov.exe", members)
                    self.assertIn(prefix + "licenses/glfw/LICENSE.md", members)
                    self.assertIn(prefix + "tests/validation_process.py", members)
                    self.assertIn(prefix + "examples/h2.molden", members)
                    self.assertFalse(any("acceptance-report" in name for name in members))
                    self.assertIn("../README.md", bundle.read(prefix + "docs/UI.md").decode())
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
            archive = args.output / "Chemical-Orbital-Visualiser-v0.4.0-pre.4-Linux-x86_64.tar.gz"
            with tarfile.open(archive) as bundle:
                binary = bundle.getmember(archive.name.removesuffix(".tar.gz") + "/cov")
                self.assertEqual(binary.mode & 0o777, 0o755)

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
                    release.prepare_macos(Path(directory) / "COV.app", Path(directory), "v0.4.0-pre.4")


if __name__ == "__main__":
    unittest.main()
