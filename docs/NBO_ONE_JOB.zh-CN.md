# 准备 COV 计算文件

[English](NBO_ONE_JOB.md) · **简体中文** · [日本語](NBO_ONE_JOB.ja.md) · [Français](NBO_ONE_JOB.fr.md)

FCHK/FCH 或 Molden 文件可用于查看分子轨道的能量、占据和形状。NBO 的轨道形状和组成分析还需要配套的 NBO 报告、`.47` 档案和轨道矩阵。各视图需要哪些文件，见[使用 NBO 结果](AOMO_NBO.zh-CN.md)。

## 从计算到 COV

先得到波函数，再进行 NBO 分析。Gaussian 可以在同一个作业中调用 NBO。`formchk` 将保留正则 MO 的检查点转换为 FCHK；GenNBO 则可以读取保存的 `.47` 档案，生成新的 NBO 报告及所需矩阵，不必重做 SCF。

FCHK 和 NBO 数据应对应相同的几何、基组和电子态。GenNBO 重新分析时可能使用不同选项，因此报告和矩阵应使用该次生成的一套。需要时可以先优化结构、计算频率；只是查看已有波函数时，并不要求增加这些步骤。

[单次作业模板](../examples/nbo-one-job/cov_one_job.py)依次运行 Gaussian/NBO、`formchk` 和 GenNBO，再将结果整理给 COV。COV 应用本身负责打开这些结果。

## 使用 Windows 模板

模板从 COV 源码树运行，需要 Python 3.11 或更新版本、Gaussian 16W A.03 和 NBO7 i8。Gaussian 与 NBO 需已安装，`gaunbo6.bat` 指向相应 NBO 安装目录。A.03 使用 `Pop=NBO6Read` 调用 NBO7。其他 Gaussian 修订版使用各自接口，本模板检查 A.03。

模板使用源码树中的[进程管理文件](../tests/validation_process.py)。

复制 [water.json](../examples/nbo-one-job/water.json)，填写分子、几何、方法、基组、电荷和多重度。示例为固定几何水分子，PBE1PBE / 6-31G(d)，电荷 0，多重度 1。自定义基组或 ECP 的完整输入段放在 `basis_ecp_tail` 中。

在源码根目录运行，使用本机安装路径和一个尚不存在的输出目录：

```powershell
python examples/nbo-one-job/cov_one_job.py `
  --recipe examples/nbo-one-job/water.json `
  --output work/water-job `
  --gaussian-bin 'C:/Gaussian/Gaussian 16 W' `
  --nbo-bin 'C:/NBO/bin'
```

这条命令只预览输入、检查路径，不创建文件或开始计算。在同一命令后加 `--run` 即提交作业。模板使用给定几何进行单点计算；优化和频率计算需另行准备输入。

最多使用三个 Gaussian 线程，整棵进程树限制为三个核心、20 GiB。Gaussian、`formchk`、GenNBO 的时限分别为 900、300、600 秒。`--memory-limit-gib` 和 `--gaussian-timeout` 可调整相应限制；`--help` 列出选项。

## 打开结果

在 COV 中打开 `work/water-job/cov-package`，或将其中的 `drop.covnbopkg` 拖入 COV。该目录包含 `canonical.fchk`、一份 `analysis.nbo` 报告、`.47` 档案和请求的轨道矩阵。

`cov-package` 中的报告、档案和矩阵来自同一次 GenNBO 分析。`canonical.fchk` 和 `analysis.nbo` 是模板采用的名称，COV 不要求固定使用这两个文件名。`.covnbopkg` 只是打开清单，不提供额外的轨道数据。

模板保留原始检查点和 Gaussian 档案，运行记录放在输入目录外。作业失败时也保留输出；再次尝试使用新目录。

## 矩阵输出

下表是模板的输出设置。NBO 允许更改这些编号，文件中的矩阵标题说明实际内容。模板输出较齐全的一组文件，单个 COV 视图可能只用到其中一部分。

| 输出选项 | 文件 | 内容 |
| --- | --- | --- |
| `AONBO=W37` | `FILE.37` | AO 基底中的 NBO 系数 |
| `NBOMO=W49` | `FILE.49` | NBO 基底中的 MO 系数 |
| `SAO=W50` | `FILE.50` | AO 重叠矩阵 |
| `NAOMO=W51` | `FILE.51` | NAO 基底中的 MO 系数 |
| `AONAO=W52` | `FILE.52` | AO 基底中的 NAO 系数 |
| `NAONBO=W53` | `FILE.53` | NAO 基底中的 NBO 系数 |
| `AONHO=W54` | `FILE.54` | AO 基底中的 NHO 系数 |
| `AONLMO=W55` | `FILE.55` | AO 基底中的 NLMO 系数 |
| `AOPNAO=W56` | `FILE.56` | AO 基底中的正交化前 NAO 系数 |
| `NAONHO=W57` | `FILE.57` | NAO 基底中的 NHO 系数 |
| `NAONLMO=W58` | `FILE.58` | NAO 基底中的 NLMO 系数 |
| `NHONBO=W59` | `FILE.59` | NHO 基底中的 NBO 系数 |
| `NBONLMO=W60` | `FILE.60` | NBO 基底中的 NLMO 系数 |
| `NLMOMO=W61` | `FILE.61` | NLMO 基底中的 MO 系数 |
| `AOMO=W62` | `FILE.62` | AO 基底中的正则 MO 系数 |

`ARCHIVE` 写出 `.47` 档案。`PRINT=3` 包含 NLMO 和 Wiberg 输出，模板还请求 `E2PERT=0.0`。完整 NBO 选项见 [cov_one_job.py](../examples/nbo-one-job/cov_one_job.py)。

`.47` 不一定包含 Fock 数据。这个模板也为能量视图准备输入，因此档案缺少重叠、密度、正则 MO 或 Fock 段时会停止。COV 中不依赖缺失数据的视图仍可使用较少的文件。开壳层视图需要对应的自旋数据；模板不会替分子选择限制或非限制方法。

## 保存与补导出

原始输入、计算日志和检查点便于以后补导出，但打开已有导出文件时，并不要求同时提供它们。正则 MO 的 FCHK 应从保留正则轨道的检查点生成；SaveNBOs 或 SaveNLMOs 可能会替换检查点中的轨道。

| 缺少什么？ | 怎样补上？ |
| --- | --- |
| 没有 FCHK，但原始 Gaussian CHK 还在 | 用 `formchk` 导出，或由 COV 调用已安装的转换器。 |
| 缺少报告段落或轨道矩阵，所需数据仍在 `.47` 中 | 增加输出选项后运行 GenNBO，使用这次生成的报告和矩阵。通常不需要新的 SCF。 |
| 档案没有 Fock 等数据 | 原软件支持时，从原始计算状态补导出。GenNBO 不能补出输入中没有的数据。 |
| 没有 `.covnbopkg` 打开清单 | 直接打开目录或选择文件。 |
| 文件来自不同计算 | 找回配套的波函数和 NBO 输出，改文件名不会改变内容。 |
| 只剩坐标 | 需要重新进行电子结构计算，获得波函数。 |

SALC 由 COV 根据导入的轨道和分子对称性构造，不需要另做一次 SCF。形状与能量所需的数据不同，详见[使用 NBO 结果](AOMO_NBO.zh-CN.md)。

矩阵输出和 `.47` 重新分析见 [NBO 7 手册](https://nbo.chem.wisc.edu/nboman.pdf) B.2.4–B.2.6 与 B.7。
