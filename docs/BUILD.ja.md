# Chemical Orbital Visualiser をソースからビルドする

[English](BUILD.md) · [简体中文](BUILD.zh-CN.md) · **日本語** · [Français](BUILD.fr.md)

この説明は現在のソースツリーに対応します。配布パッケージと動作要件は[ダウンロードの説明](../README.ja.md#ダウンロード)にあります。

## 必要な環境

- CMake 3.28 以降
- C++20 対応コンパイラー、Ninja、Git
- CMake が指定バージョンの Eigen、GLFW、Dear ImGui のソースを取得するためのインターネット接続
- デスクトップ表示には OpenGL 開発ファイル。Linux では GLFW の X11 または Wayland 開発依存関係も必要です。実行時には OpenGL 2.1 互換コンテキストが必要です。

Windows では、「C++ によるデスクトップ開発」を導入した Visual Studio 2022 の x64 開発者コマンドプロンプトまたは開発者 PowerShell を使います。Debian/Ubuntu では C++20 ツールチェーンと対応する OpenGL/GLFW 開発パッケージを使います。macOS では Metal コンパイラーを含む Xcode 16.2 以降を使います。FreeBSD のソースビルドは CPU バックエンドを使い、`mesa-libs` と `libglvnd` の開発パッケージが必要です。

## プリセット

ソースのルートディレクトリで実行します。

```text
cmake --preset cpu-core
cmake --build --preset cpu-core --parallel
ctest --preset cpu-core
```

`cpu-core` はコアライブラリーとテストをビルドし、デスクトップ表示や GPU SDK は必要としません。CPU で軌道グリッドを計算するデスクトップ版をビルドする場合は、次を使います。

```text
cmake --preset portable
cmake --build --preset portable --parallel
```

Linux/macOS では `build/portable/cov`、Windows では `build\portable\cov.exe` を実行します。例：

```bash
./build/portable/cov --compute-backend=cpu ./examples/h2.molden
```

Windows PowerShell：

```powershell
.\build\portable\cov.exe --compute-backend=cpu .\examples\h2.molden
```

`cuda`、`metal`、`webgpu` のプリセットはそれぞれ対応するバックエンドを有効にします。`metal` は macOS にのみ表示されます。各デスクトップ版プリセットは個別のビルドディレクトリを使います。CUDA プリセットには CUDA Toolkit 12.8 が必要で、実行には互換性のある NVIDIA ドライバーも必要です。既定のデバイスコード対象は `60, 61, 70, 75, 80, 86, 89, 90, 100, 120` と `compute_120` PTX です。この一覧には CUDA 12.8 を使います。CUDA 13 では一部の古い対象をビルドできません。一つの対応アーキテクチャに限定する場合の例：

```bash
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=86-real
cmake --build --preset cuda --parallel
```

macOS では `cmake --preset metal` と `cmake --build --preset metal` で OpenGL 2.1 のデスクトップ版と Metal 計算モジュールをビルドします。任意の WebGPU を有効にするには、wgpu-native **v29.0.1.1** を展開し、`cmake --preset webgpu -DCOV_WGPU_NATIVE_ROOT=/absolute/path/to/wgpu-native`、`cmake --build --preset webgpu` を実行します。Windows では起動時に WebGPU のランタイム DLL も `cov.exe` の隣に置きます。

## AMD と Intel の計算モジュール

HIP と SYCL は別々のコンパイラーツールチェーンを使うため、それぞれ独立した共有ライブラリーとしてビルドします。デスクトップ本体は通常の C++ ビルドです。Linux ではソースのルートで次を実行します。

```bash
# AMD HIP/ROCm：インストール済み SDK が対応する gfx 対象を選択します。
cmake -S src/compute/hip -B build/hip -DCOV_ROOT="$PWD" \
  -DCMAKE_HIP_COMPILER=/opt/rocm/llvm/bin/clang++ -DCMAKE_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

# Intel oneAPI DPC++：Intel Level Zero GPU を選択するモジュールです。
source /opt/intel/oneapi/setvars.sh --force
cmake -S src/compute/sycl -B build/sycl -DCOV_ROOT="$PWD" -DCMAKE_CXX_COMPILER=icpx
cmake --build build/sycl --parallel
```

Windows では、対応する SDK のインストール後に Visual Studio 2022 の x64 開発者プロンプトで実行します。

```bat
cmake -S src/compute/hip -B build/hip -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER="C:/Program Files/AMD/ROCm/7.2/bin/hipcc.exe" -DCOV_HIP_ARCHITECTURES=gfx1100
cmake --build build/hip --parallel

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
cmake -S src/compute/sycl -B build/sycl -G Ninja -DCOV_ROOT="%CD%" -DCMAKE_CXX_COMPILER=icx-cl
cmake --build build/sycl --parallel
```

生成された `cov_compute_hip.dll` または `cov_compute_sycl.dll`（Windows）、`libcov_compute_hip.so` または `libcov_compute_sycl.so`（Linux）を `cov.exe`/`cov` の隣に置きます。環境変数 `COV_COMPUTE_MODULE_DIR` にモジュールのあるディレクトリの**絶対パス**を設定する方法もあります。Windows では対応するランタイム DLL をモジュールの隣に置くか、`COV_COMPUTE_RUNTIME_DIR` に SDK ランタイムディレクトリの絶対パスを設定します。Linux では HIP/oneAPI のランタイムをシステムローダーから利用できるようにします。

`--compute-backend=auto` は表示 GPU に対応するネイティブバックエンドを優先し、ほかの利用可能なネイティブバックエンドを試した後、CPU に切り替えます。`cpu`、`cuda`、`hip`、`sycl`、`metal`、`webgpu` を明示的に指定できます。`--compute-device=N` は GPU バックエンドのデバイスを 0 始まりの番号で選択します。WebGPU は明示的な選択が必要です。OpenCL は今後の開発用で、このソースツリーにビルドオプションはありません。

## 入力ファイルと読み込み

COV は Gaussian の FCHK/FCH と Molden の波動関数を読み込みます。CHK のローカル変換には別途インストールした `formchk` が必要です。自動検出で見つからなければ、`COV_FORMCHK` に実行ファイルのパスを設定します。Gaussian の `.log` または `.out` は対応する波動関数の情報を補えますが、FCHK/FCH の代わりにはなりません。

入力は最大 100 原子です。デカルトおよび実球面調和関数の `s/p/d/f/g` 基底関数に対応します。Molden では、展開後の基底関数数と各分子軌道の係数が一致する必要があります。範囲外の係数インデックスは入力または殻の規約のエラーとして扱い、推測による修復は行いません。

## 描画の構成

CPU が入力を読み込み、選択した軌道をグリッド計算用にまとめます。CPU または選択した計算バックエンドが 3 次元グリッドを計算し、デスクトップ表示側が OpenGL の 3D テクスチャーに転送します。CUDA はテクスチャーに直接書き込むこともできます。描画処理は、そのテクスチャーから正と負の等値面を描きます。等値面の値を変更しても軌道グリッドの再計算は不要です。

各グリッド点で次を計算します。

$$
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r).
$$

グリッド点数と基底関数数の積の大きさを持つ行列は確保しません。グリッド解像度は 64³、128³、256³、512³ です。

Gaussian/NBO ファイルの準備は[単一ジョブのテンプレート](NBO_ONE_JOB.ja.md)を参照してください。
