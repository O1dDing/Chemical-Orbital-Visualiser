# Chemical Orbital Visualiser (COV)

[English](README.md) · [简体中文](README.zh-CN.md) · **日本語** · [Français](README.fr.md)

軌道のエネルギーと電子占有数、原子軌道と分子軌道のつながりを調べられます。

COV は軌道のエネルギーと占有数を操作可能なエネルギー準位図に表示し、図とデータを書き出せます。Gaussian の FCHK/FCH ファイルと Molden ファイルを読み込みます。

v0.4 プレビュー版では NBO 解析にも対応しました。同じ計算の NBO 出力を使って、原子軌道や局在化軌道が分子軌道にどのように寄与するかを確認し、準位図上で軌道間のつながりをたどり、電荷や結合に関する情報を調べられます。

軌道を選ぶと、その形状を 3D で表示できます。プレビュー版は Windows、macOS、Linux に対応し、各環境で CPU 計算を利用できます。画面表示は英語、簡体字中国語、日本語、フランス語に対応しています。

## 主な機能

- **エネルギー準位図を読む。** 軌道のエネルギーと電子占有数を確認し、エネルギーの単位を切り替え、図を PNG または SVG で書き出せます。
- **対称性と組成を見る — v0.4 プレビュー版。** 現在の表示に対応する対称性と軌道群の名前、制限開殻計算の対応軌道に共通するエネルギー、対象となる錯体の金属・配位子組成を確認できます。簡略表示では、内殻・配位子の背景と対応する AO/SALC 寄与をまとめて絞り込めます。
- **図の詳細を調べる — v0.4 プレビュー版。** 縮退軌道群と等価な結合の結合表示をそろえ、軌道詳細のフローティング表示やラベルの切り替えを利用できます。
- **軌道間のつながりをたどる — v0.4 プレビュー版。** 原子軌道や局在化軌道の分子軌道への寄与を表示します。つながりを選ぶと、その寄与を詳しく確認できます。
- **電荷と結合を調べる — v0.4 プレビュー版。** NBO ファイルにデータが含まれていれば、NPA 電荷、スピン分布、Wiberg 結合指数、供与体–受容体相互作用を確認できます。
- **図とデータを書き出す。** エネルギー準位図を PNG または SVG で保存し、対応するデータを CSV または JSON で書き出せます。
- **軌道を探す。** HOMO や LUMO に移動し、軌道一覧を検索できます。占有軌道、空軌道、内殻軌道、価電子軌道だけを表示することもできます。
- **軌道を 3D で見る。** 軌道を選び、等値面を調整し、分子の表示を回転・拡大縮小できます。

## ダウンロード

| バージョン | 主な内容 | ダウンロード |
|---|---|---|
| [安定版 v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | 軌道のエネルギーと占有数、エネルギー準位図、3D 表示 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [プレビュー版 v0.4.0-pre.4](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.4) | NBO 解析、軌道の組成とつながり、対称性と結合の表示 | [Windows ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Windows-x64.zip) · [macOS ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-macOS-universal.zip) · [Linux tar.gz](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Linux-x86_64.tar.gz) |

v0.3.0 の配布ファイルには旧製品名が残っており、NVIDIA RTX 50 シリーズが必要です。プレビュー版の 3 環境用パッケージは同じソース版から作成しています。

## はじめに

1. お使いの環境用のパッケージをダウンロードして展開します。Windows は `cov.exe`、macOS は `Chemical Orbital Visualiser.app`、Linux は `./cov` を起動します。
2. Gaussian の FCHK/FCH ファイルか互換性のある Molden ファイルを開くか、ウィンドウにドラッグします。
3. 軌道一覧とエネルギー準位図でエネルギーと占有数を確認します。準位を選んで対応する軌道を表示するか、図を書き出します。
4. v0.4 プレビュー版では、計算フォルダーを開くと波動関数と NBO ファイルをまとめて読み込めます。フォルダーに複数の計算がある場合は、開く計算を選んでください。

`formchk` がインストールされていれば、COV は Gaussian CHK を FCHK に変換できます。FCHK/FCH と Molden ファイルには、MO のエネルギー、占有数、形状のデータが含まれます。NBO 軌道の形状や組成の解析には、対応するレポート、`.47` アーカイブ、軌道行列も必要です。

## 入力ファイルと動作要件

- **波動関数：** Gaussian の `.fchk` / `.fch`、または互換性のある `.molden` / `.mol` / `.input` ファイル。Gaussian の `.chk` は、インストール済みの `formchk` を使って開けます。
- **NBO ファイル — v0.4 プレビュー版：** [計算ファイルの準備](docs/NBO_ONE_JOB.ja.md)では作成方法を、[NBO の結果を使う](docs/AOMO_NBO.ja.md)では各表示に必要なファイルを説明しています。
- **分子の大きさ：** 入力ファイルあたり最大 100 原子。
- **プレビュー版の環境：** Windows x64、Apple Silicon または Intel の macOS 12 以降、Linux x86_64（パッケージの基準は Ubuntu 22.04）。
- **プレビュー版の表示と計算：** OpenGL 2.1 以降。Windows 版は CUDA 12.8 を含み、CPU にも切り替えられます。CUDA の利用には互換性のある NVIDIA ドライバーが必要です。macOS 版は Metal と CPU、Linux 版は CPU で計算します。プレビュー版の起動に NVIDIA GPU は必要ありません。追加の GPU モジュールはソースからビルドできます。
- **安定版 v0.3.0：** Windows パッケージには NVIDIA RTX 50 シリーズ、互換性のあるドライバー、OpenGL 2.1 以降が必要です。

## ドキュメント

- [COV の使い方（英語）](docs/UI.md)
- [NBO の結果を使う](docs/AOMO_NBO.ja.md) — v0.4 プレビュー版
- [計算ファイルの準備](docs/NBO_ONE_JOB.ja.md) — ソースツリーの実行テンプレートを含みます
- [ソースからのビルド（英語）](docs/BUILD.md)
- [安定版のリリースノート](docs/releases/v0.3.0.md) · [プレビュー版のリリースノート](docs/releases/v0.4.0-pre.4.ja.md) · [過去の v0.3 プレビュー版](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0-pre-archive)

[問題の報告・機能の提案](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose)。

## ライセンス

[Apache License 2.0](LICENSE)。同梱ライブラリのライセンスは[サードパーティーに関するお知らせ](THIRD_PARTY_NOTICES.md)に記載しています。
