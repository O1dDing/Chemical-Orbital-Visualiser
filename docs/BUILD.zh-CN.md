# 从源码编译 Chemical Orbital Visualiser

[English](BUILD.md)

本文适用于当前源码树。已发布的 Windows 下载包有各自的平台要求；修改源码构建配置不会改变已有下载包。

## 环境要求

- CMake 3.28 或更新版本
- 支持 C++20 的编译器、Ninja 和 Git
- 联网获取 CMake 固定版本的 Eigen、GLFW、Dear ImGui 源码依赖
- 构建桌面程序需要 OpenGL 开发文件；Linux 还需要 GLFW 对应的 X11 或 Wayland 开发依赖。运行时需要 OpenGL 2.1 兼容上下文

Windows 使用装有“使用 C++ 的桌面开发”工作负载的 Visual Studio 2022 x64 开发者命令提示符或 PowerShell。Debian/Ubuntu 需要 C++20 工具链及相应的 OpenGL/GLFW 开发包；macOS 使用 Xcode 16.2 或更新版本，包含 Metal 编译器。FreeBSD 源码构建使用 CPU 后端，并需要 `mesa-libs`、`libglvnd` 开发包。

## 构建预设

以下命令均在源码根目录运行：

```text
cmake --preset cpu-core
cmake --build --preset cpu-core --parallel
ctest --preset cpu-core
```

`cpu-core` 构建核心库和测试，不构建桌面程序，也不需要 GPU SDK。构建使用 CPU 网格计算的桌面程序：

```text
cmake --preset portable
cmake --build --preset portable --parallel
```

Linux/macOS 的程序路径为 `build/portable/cov`，Windows 为 `build\portable\cov.exe`。例如：

```bash
./build/portable/cov --compute-backend=cpu ./examples/h2.molden
```

Windows PowerShell 使用：

```powershell
.\build\portable\cov.exe --compute-backend=cpu .\examples\h2.molden
```

`cuda`、`metal`、`webgpu` 预设分别启用对应后端；`metal` 只在 macOS 显示。每个桌面预设使用独立构建目录。CUDA 预设需要 CUDA Toolkit 12.8；运行 CUDA 后端还需要兼容的 NVIDIA 驱动。默认设备代码目标为 `60、61、70、75、80、86、89、90、100、120`，另含 `compute_120` PTX。这组目标请使用 CUDA 12.8；CUDA 13 不再编译其中较旧的目标。如只需一个受支持架构，可指定：

```bash
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=86-real
cmake --build --preset cuda --parallel
```

macOS 运行 `cmake --preset metal` 和 `cmake --build --preset metal`，构建 OpenGL 2.1 桌面程序及 Metal 计算模块。可选 WebGPU 需要解压 wgpu-native **v29.0.1.1**，再运行 `cmake --preset webgpu -DCOV_WGPU_NATIVE_ROOT=/absolute/path/to/wgpu-native` 和 `cmake --build --preset webgpu`。Windows 启动时还需将 WebGPU 运行库 DLL 放在 `cov.exe` 旁。

## AMD 和 Intel 计算模块

HIP 与 SYCL 使用各自的编译器，应分别构建为共享库；主程序仍由普通 C++ 工具链编译。Linux 下在源码根目录运行：

```bash
# AMD HIP/ROCm；选择已安装 SDK 支持的 gfx 目标。
cmake -S src/compute/hip -B build/hip -DCOV_ROOT="$PWD" \
  -DCMAKE_HIP_COMPILER=/opt/rocm/llvm/bin/clang++ -DCMAKE_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

# Intel oneAPI DPC++；模块仅选择 Intel Level Zero GPU。
source /opt/intel/oneapi/setvars.sh --force
cmake -S src/compute/sycl -B build/sycl -DCOV_ROOT="$PWD" -DCMAKE_CXX_COMPILER=icpx
cmake --build build/sycl --parallel
```

Windows 下安装相应 SDK 后，在 Visual Studio 2022 x64 开发者命令提示符运行：

```bat
cmake -S src/compute/hip -B build/hip -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER="C:/Program Files/AMD/ROCm/7.2/bin/hipcc.exe" -DCOV_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
cmake -S src/compute/sycl -B build/sycl -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER=icx-cl
cmake --build build/sycl --parallel
```

将生成的 `cov_compute_hip.dll` 或 `cov_compute_sycl.dll`（Windows），或 `libcov_compute_hip.so` 或 `libcov_compute_sycl.so`（Linux），放在 `cov.exe`/`cov` 旁。也可将 `COV_COMPUTE_MODULE_DIR` 设为模块所在目录的**绝对路径**。Windows 上将对应运行库 DLL 放在模块旁，或将 `COV_COMPUTE_RUNTIME_DIR` 设为 SDK 运行库所在目录的绝对路径。Linux 上的 HIP 或 oneAPI 运行库仍须可由系统加载。

`--compute-backend=auto` 优先选择与显示 GPU 匹配的原生后端，再尝试其他可用原生后端，最后回退到 CPU。可显式指定 `cpu`、`cuda`、`hip`、`sycl`、`metal`、`webgpu`；`--compute-device=N` 指定 GPU 后端的零起始设备编号。WebGPU 需显式选择。OpenCL 保留以后开发，本源码树没有对应构建选项。

## 输入文件

COV 读取 Gaussian FCHK/FCH 和 Molden 波函数。打开 CHK 需要单独安装 `formchk`；找不到转换器时，将 `COV_FORMCHK` 设为可执行文件的完整路径。Gaussian `.log` 或 `.out` 可补充波函数信息，不能替代 FCHK/FCH。

每个输入最多 100 个原子。支持 Cartesian 和实球谐 `s/p/d/f/g` 基函数。Molden 展开后的基函数数目须与轨道系数一致。

Gaussian/NBO 文件准备见[单次作业模板](NBO_ONE_JOB.zh-CN.md)。
