# Build Chemical Orbital Visualiser

[简体中文](BUILD.zh-CN.md)

These instructions apply to this source tree. Published Windows downloads have their own platform requirements; changing a source build does not change an existing download.

## Requirements

- CMake 3.28 or newer
- C++20 compiler, Ninja, and Git
- Internet access for CMake's pinned Eigen, GLFW, and Dear ImGui source dependencies
- For the viewer build: OpenGL development files and, on Linux, the GLFW X11 or Wayland development dependencies; running the viewer needs an OpenGL 2.1 compatibility context

On Windows, use a Visual Studio 2022 x64 developer command prompt or developer PowerShell with the Desktop development with C++ workload. On Debian or Ubuntu, use a C++20 toolchain and the matching OpenGL/GLFW development packages. macOS builds need Xcode command-line tools. FreeBSD is a source portability target for the CPU build.

## Presets

Run these commands from the source root:

```text
cmake --preset cpu-core
cmake --build --preset cpu-core --parallel
ctest --preset cpu-core
```

`cpu-core` builds the core libraries and tests without the desktop viewer or a GPU SDK. To build the desktop viewer with CPU grid computation:

```text
cmake --preset portable
cmake --build --preset portable --parallel
```

Run `build/portable/cov` on Linux or macOS, or `build\portable\cov.exe` on Windows. For example:

```bash
./build/portable/cov --compute-backend=cpu ./examples/h2.molden
```

In Windows PowerShell:

```powershell
.\build\portable\cov.exe --compute-backend=cpu .\examples\h2.molden
```

The `cuda`, `metal`, and `webgpu` presets enable their named backends. `metal` appears only on macOS. Each viewer preset uses a separate build directory. The CUDA preset needs CUDA Toolkit 12.8; running its backend needs a compatible NVIDIA driver. Its default device-code targets are `60, 61, 70, 75, 80, 86, 89, 90, 100, 120` plus `compute_120` PTX. Use CUDA 12.8 for this list; CUDA 13 does not build the older targets. To build only one supported architecture, for example:

```bash
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=86-real
cmake --build --preset cuda --parallel
```

On macOS, `cmake --preset metal` followed by `cmake --build --preset metal` builds the OpenGL 2.1 viewer and its Metal compute module. To enable optional WebGPU, extract wgpu-native **v29.0.1.1**, then run `cmake --preset webgpu -DCOV_WGPU_NATIVE_ROOT=/absolute/path/to/wgpu-native` and `cmake --build --preset webgpu`. The Windows WebGPU runtime DLL must also be beside `cov.exe` at launch.

## AMD and Intel compute modules

HIP and SYCL use separate compiler toolchains, so build each as a standalone shared library. The main viewer remains a normal C++ build. On Linux, from the source root:

```bash
# AMD HIP/ROCm; choose a gfx target supported by the installed SDK.
cmake -S src/compute/hip -B build/hip -DCOV_ROOT="$PWD" \
  -DCMAKE_HIP_COMPILER=/opt/rocm/llvm/bin/clang++ -DCMAKE_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

# Intel oneAPI DPC++; the module selects Intel Level Zero GPUs.
source /opt/intel/oneapi/setvars.sh --force
cmake -S src/compute/sycl -B build/sycl -DCOV_ROOT="$PWD" -DCMAKE_CXX_COMPILER=icpx
cmake --build build/sycl --parallel
```

On Windows, use a Visual Studio 2022 x64 developer prompt after installing the matching SDK:

```bat
cmake -S src/compute/hip -B build/hip -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER="C:/Program Files/AMD/ROCm/7.2/bin/hipcc.exe" -DCOV_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
cmake -S src/compute/sycl -B build/sycl -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER=icx-cl
cmake --build build/sycl --parallel
```

Put the resulting `cov_compute_hip.dll` or `cov_compute_sycl.dll` (Windows), or `libcov_compute_hip.so` or `libcov_compute_sycl.so` (Linux), beside `cov.exe`/`cov`. Alternatively, set `COV_COMPUTE_MODULE_DIR` to the **absolute** directory containing the modules. Keep the corresponding HIP or oneAPI runtime available to the system loader.

`--compute-backend=auto` prefers the native backend matching the display GPU, then other available native backends, and falls back to CPU computation. Explicit choices are `cpu`, `cuda`, `hip`, `sycl`, `metal`, and `webgpu`; `--compute-device=N` selects a zero-based GPU index for a GPU backend. WebGPU is opt-in. OpenCL is reserved for future work and has no build option in this tree.

## Inputs and parser

COV reads Gaussian FCHK/FCH and Molden wavefunctions. Local CHK conversion requires a separately installed `formchk`; set `COV_FORMCHK` to its path if automatic discovery fails. A Gaussian `.log` or `.out` file can add information to the corresponding wavefunction, but does not replace the FCHK/FCH input.

The current input limit is 100 atoms. Cartesian and real-spherical `s/p/d/f/g` basis functions are supported. For a Molden file, the expanded basis count must match the coefficients of each molecular orbital. An out-of-range coefficient index is treated as an input or shell-convention error; the parser does not guess a repair.

## Rendering architecture

The CPU parses the input and packs the selected orbital for grid computation. The CPU or selected compute backend evaluates a three-dimensional grid; the viewer uploads it to an OpenGL 3D texture. CUDA can also write directly to the texture. The renderer draws positive and negative isosurfaces from that texture. Changing the isovalue changes the display without recalculating the orbital grid.

For each grid point, the evaluator computes

\[
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r).
\]

It does not allocate a full grid-points-by-basis-functions matrix. Available grid resolutions are 64³, 128³, 256³ and 512³.

For Gaussian/NBO file preparation, see the [one-job template](NBO_ONE_JOB.md).
