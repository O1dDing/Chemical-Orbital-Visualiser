# Chemical Orbital Visualiser (COV)

[English](README.md) · **简体中文** · [日本語](README.ja.md) · [Français](README.fr.md)

查看轨道能量、电子占据情况，以及原子轨道与分子轨道之间的联系。

COV 可以用交互式能级图查看轨道能量和电子占据，并导出图片与数据。软件支持 Gaussian FCHK/FCH 和 Molden 文件。

v0.4 预览版加入了 NBO 分析。配合计算生成的 NBO 文件，可以查看原子轨道和局域轨道怎样组成分子轨道，沿着图中的连线查看各项贡献，以及电荷和成键信息。

选中轨道后，还可以查看它的三维形状。预览版支持 Windows、macOS 和 Linux，各平台均可使用 CPU 计算。界面支持中文、英文、日文和法文。

## 功能

- **查看能级图。** 查看轨道能量和电子占据情况、切换能量单位，并将能级图导出为 PNG 或 SVG。
- **查看对称性与组成 — v0.4 预览版。** 查看当前视图的对称性及轨道组名称、受限开壳层对应轨道的共享能量，以及适用配合物的金属/配体组成。简要视图可同时筛除内层和配体背景及匹配的 AO/SALC 贡献。
- **查看图中细节 — v0.4 预览版。** 简并轨道组和等价键采用一致的成键显示，并提供浮动轨道详情和标签开关。
- **追踪轨道联系 — v0.4 预览版。** 查看原子轨道和局域轨道怎样组成分子轨道，点击连线查看各项贡献。
- **查看电荷与成键 — v0.4 预览版。** 如果 NBO 文件中包含相关数据，可以查看 NPA 电荷、自旋布居、Wiberg 键级指数和给体–受体相互作用。
- **导出图像和数据。** 将能级图保存为 PNG 或 SVG，并将相应数据导出为 CSV 或 JSON。
- **浏览轨道。** 跳转到 HOMO 或 LUMO、搜索轨道列表，并显示已占据、未占据、内层或价层轨道。
- **查看三维轨道。** 选择轨道、调整其等值面，并旋转或缩放分子视图。

## 下载

| 版本 | 包含的功能 | 下载 |
|---|---|---|
| [稳定版 v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | 轨道能量、占据情况、能级图和三维视图 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [预览版 v0.4.0-pre.4](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.4) | NBO 分析、轨道组成与联系、对称性及成键视图 | [Windows ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Windows-x64.zip) · [macOS ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-macOS-universal.zip) · [Linux tar.gz](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Linux-x86_64.tar.gz) |

v0.3.0 下载包保留旧产品名称，需要 NVIDIA RTX 50 系列显卡。预览版的三个平台包来自同一源码版本。

## 开始使用

1. 下载并解压对应平台的程序包。Windows 运行 `cov.exe`；macOS 打开 `Chemical Orbital Visualiser.app`；Linux 运行 `./cov`。
2. 打开 Gaussian FCHK/FCH 或兼容的 Molden 文件，也可以将文件拖入窗口。
3. 在轨道列表和能级图中查看能量与占据情况。选择一个能级以查看对应轨道，或导出能级图。
4. 使用 v0.4 预览版时，可以打开计算文件夹，一并加载其中的波函数和 NBO 文件。如果文件夹中有多次计算，请选择要打开的一次。

已安装 `formchk` 时，COV 可以将 Gaussian CHK 转换为 FCHK。FCHK/FCH 或 Molden 文件可用于查看 MO 的能量、占据和形状；NBO 的轨道形状与组成分析还需要配套的报告、`.47` 档案和轨道矩阵。

## 输入文件与运行要求

- **波函数：** Gaussian `.fchk` / `.fch`，或兼容的 `.molden` / `.mol` / `.input` 文件。已安装 `formchk` 时也可打开 Gaussian `.chk`。
- **NBO 文件 — v0.4 预览版：** [准备计算文件](docs/NBO_ONE_JOB.zh-CN.md)说明怎样生成文件；[使用 NBO 结果](docs/AOMO_NBO.zh-CN.md)列出各视图需要哪些文件。
- **分子大小：** 每个输入最多包含 100 个原子。
- **预览版平台：** Windows x64；Apple Silicon 或 Intel 的 macOS 12 及更新版本；Linux x86_64，程序包以 Ubuntu 22.04 为基线。
- **预览版图形与计算：** OpenGL 2.1 或更新版本。Windows 包含 CUDA 12.8，并可回退至 CPU。使用 CUDA 需要兼容的 NVIDIA 驱动。macOS 提供 Metal 和 CPU 计算；Linux 使用 CPU 计算。打开预览版不要求 NVIDIA 显卡。其他可选 GPU 模块可从源码编译。
- **稳定版 v0.3.0：** Windows 包需要 NVIDIA RTX 50 系列显卡、兼容的驱动，以及 OpenGL 2.1 或更新版本。

## 文档

- [使用 COV](docs/UI.zh-CN.md)
- [使用 NBO 结果](docs/AOMO_NBO.zh-CN.md) — v0.4 预览版
- [准备计算文件](docs/NBO_ONE_JOB.zh-CN.md) — 含源码树中的单次作业模板
- [源码编译](docs/BUILD.zh-CN.md)
- [稳定版发布说明](docs/releases/v0.3.0.md) · [预览版发布说明](docs/releases/v0.4.0-pre.4.zh-CN.md) · [旧版 v0.3 预览](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0-pre-archive)

[反馈问题或建议功能](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose)。

## 许可证

[Apache License 2.0](LICENSE)。随附程序库的许可证列于[第三方声明](THIRD_PARTY_NOTICES.md)。
