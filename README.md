# Chemical Orbital Visualiser (COV)

**English** · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · [Français](README.fr.md)

Explore orbital energies, occupations and the connections between atomic and molecular orbitals.

COV displays orbital energies and occupations in interactive energy-level diagrams, with figures and data available for export. It reads Gaussian FCHK/FCH and Molden files.

The v0.4 preview adds NBO analysis. With the calculation's NBO output, you can see how atomic and localised orbitals contribute to molecular orbitals, follow their connections in the diagram, and inspect charge and bonding information.

Select an orbital to view its shape in 3D. The preview runs on Windows, macOS and Linux, with CPU computation available on each platform. The interface is available in English, Simplified Chinese, Japanese and French.

## What you can do

- **Read energy-level diagrams.** See orbital energies and electron occupations, change energy units, and export a diagram as PNG or SVG.
- **Read symmetry and composition — v0.4 preview.** See symmetry and group names for the current view, shared energies for matching restricted-open-shell orbitals, and metal/ligand composition for applicable complexes. A brief view can filter core and ligand background together with matching AO/SALC contributions.
- **Inspect diagram details — v0.4 preview.** Use consistent bonding displays for degenerate groups and equivalent bonds, floating orbital details and label toggles.
- **Follow orbital connections — v0.4 preview.** See how atomic and localised orbitals contribute to molecular orbitals. Select a connection to inspect its contribution.
- **Inspect charge and bonding — v0.4 preview.** Read NPA charges, spin populations, Wiberg bond indices and donor–acceptor interactions when the NBO files contain them.
- **Export figures and data.** Save energy diagrams as PNG or SVG, and the corresponding data as CSV or JSON.
- **Browse orbitals.** Jump to HOMO or LUMO, search the orbital list, and show occupied, virtual, core or valence orbitals.
- **View orbitals in 3D.** Select an orbital, adjust its isosurface, and rotate or zoom the molecular view.

## Download

| Version | Includes | Downloads |
|---|---|---|
| [Stable v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | Orbital energies, occupations, energy-level diagrams and 3D views | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [Preview v0.4.0-pre.4](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.4) | NBO analysis, orbital composition and connections, symmetry and bonding views | [Windows ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Windows-x64.zip) · [macOS ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-macOS-universal.zip) · [Linux archive](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.4/Chemical-Orbital-Visualiser-v0.4.0-pre.4-Linux-x86_64.tar.gz) |

The v0.3.0 download retains the earlier product name and requires an NVIDIA RTX 50-series GPU. The three preview packages use the same source version.

## Get started

1. Download and extract the package for your platform. On Windows, run `cov.exe`; on macOS, open `Chemical Orbital Visualiser.app`; on Linux, run `./cov`.
2. Open a Gaussian FCHK/FCH or compatible Molden file, or drag it into the window.
3. Browse the energies and occupations in the orbital list and energy diagram. Select a level to view its orbital, or export the diagram.
4. In the v0.4 preview, open the calculation folder to load its wavefunction and NBO files together. If the folder contains several calculations, choose the one to open.

With `formchk` installed, COV can convert Gaussian CHK files to FCHK. FCHK/FCH and Molden files provide MO energies, occupations and shapes. NBO orbital shapes and composition analysis also need the matching report, `.47` archive and orbital matrices.

## Input and requirements

- **Wavefunctions:** Gaussian `.fchk` / `.fch`, or compatible `.molden` / `.mol` / `.input` files. Gaussian `.chk` can be opened with an installed `formchk`.
- **NBO files — v0.4 preview:** [Prepare calculation files](docs/NBO_ONE_JOB.md) explains how to produce them; [Using NBO results](docs/AOMO_NBO.md) lists the files for each view.
- **Molecule size:** up to 100 atoms per input.
- **Preview platforms:** Windows x64; macOS 12 or newer on Apple Silicon or Intel; Linux x86_64 with Ubuntu 22.04 as the package baseline.
- **Preview graphics and computation:** OpenGL 2.1 or newer. Windows includes CUDA 12.8 with CPU fallback. CUDA use needs a compatible NVIDIA driver. macOS includes Metal and CPU computation; Linux uses CPU computation. An NVIDIA GPU is not required to open the preview. Optional GPU modules can be built from source.
- **Stable v0.3.0:** the Windows package requires an NVIDIA RTX 50-series GPU, a compatible driver and OpenGL 2.1 or newer.

## Documentation

- [Using COV](docs/UI.md)
- [Using NBO results](docs/AOMO_NBO.md) — v0.4 preview
- [Prepare calculation files](docs/NBO_ONE_JOB.md) — includes the source-tree job template
- [Building from source](docs/BUILD.md)
- [Stable release notes](docs/releases/v0.3.0.md) · [Preview release notes](docs/releases/v0.4.0-pre.4.md) · [Older v0.3 previews](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0-pre-archive)

[Report a problem or suggest a feature](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose).

## Licence

[Apache License 2.0](LICENSE). Bundled library licences are listed in [Third-party notices](THIRD_PARTY_NOTICES.md).
