# Prepare calculation files for COV

**English** · [简体中文](NBO_ONE_JOB.zh-CN.md) · [日本語](NBO_ONE_JOB.ja.md) · [Français](NBO_ONE_JOB.fr.md)

FCHK/FCH or Molden files provide molecular orbital energies, occupations and shapes. NBO orbital shapes and composition analysis need the corresponding NBO report, `.47` archive and orbital matrices. The files required for each view are listed in [Using NBO results](AOMO_NBO.md).

## From a calculation to COV

The wavefunction comes first. Gaussian can then run NBO as part of the same job. `formchk` exports the canonical checkpoint as FCHK; GenNBO can read the saved `.47` archive to produce another NBO report and the requested matrices without repeating the SCF calculation.

The FCHK and NBO data describe the same geometry, basis and electronic state. A GenNBO reanalysis may use different analysis options, so its report and matrices belong together. Optimization and frequency calculations can precede this step when needed; they are not required just to display an existing wavefunction.

The [one-job template](../examples/nbo-one-job/cov_one_job.py) submits Gaussian/NBO, `formchk` and GenNBO in sequence, then collects their output for COV. The COV application itself opens the results.

## Use the Windows template

The template runs from the COV source tree with Python 3.11 or newer, Gaussian 16W A.03 and NBO7 i8. Gaussian and NBO must already be installed, with `gaunbo6.bat` pointing to that NBO installation. A.03 uses `Pop=NBO6Read` to call NBO7. Other Gaussian revisions use their own interface setup; this template checks for A.03.

The template uses [tests/validation_process.py](../tests/validation_process.py) from the source tree.

Copy [water.json](../examples/nbo-one-job/water.json) and set the molecule, geometry, method, basis, charge and multiplicity. The example is a fixed-geometry PBE1PBE / 6-31G(d) water single point, charge 0, multiplicity 1. `basis_ecp_tail` accepts the explicit basis/ECP block when needed.

Run from the source-tree root, using your installation paths and a new output directory:

```powershell
python examples/nbo-one-job/cov_one_job.py `
  --recipe examples/nbo-one-job/water.json `
  --output work/water-job `
  --gaussian-bin 'C:/Gaussian/Gaussian 16 W' `
  --nbo-bin 'C:/NBO/bin'
```

This previews the input and checks the paths without creating files or starting a calculation. Add `--run` to the same command to submit it. The output directory must not already exist. The template keeps the supplied geometry and uses one single-point step; optimization and frequency calculations need their own inputs.

The template permits up to three Gaussian threads. The process tree is limited to three cores and 20 GiB; timeouts are 900 seconds for Gaussian, 300 for `formchk` and 600 for GenNBO. `--memory-limit-gib` and `--gaussian-timeout` adjust the documented limits; `--help` lists the options.

## Open the result

Open `work/water-job/cov-package` in COV, or drag its `drop.covnbopkg` into COV. The folder contains `canonical.fchk`, one `analysis.nbo` report, the `.47` archive and the requested orbital matrices.

The report, archive and matrices in `cov-package` come from the same GenNBO analysis. `canonical.fchk` and `analysis.nbo` are template filenames; COV does not require those exact names. The `.covnbopkg` file is an opening list, not an extra source of orbital data.

The template preserves the original checkpoint and Gaussian archive. Run records stay outside the input folder. Failed jobs retain their output; a further attempt uses a new directory.

## Matrix output

These are the template's output settings. The numbers can be changed in NBO; the matrix headings identify the contents. The template requests a broad set of files, although a particular COV view may need only some of them.

| Output option | File | Contents |
| --- | --- | --- |
| `AONBO=W37` | `FILE.37` | NBO coefficients in the AO basis |
| `NBOMO=W49` | `FILE.49` | MO coefficients in the NBO basis |
| `SAO=W50` | `FILE.50` | AO overlap matrix |
| `NAOMO=W51` | `FILE.51` | MO coefficients in the NAO basis |
| `AONAO=W52` | `FILE.52` | NAO coefficients in the AO basis |
| `NAONBO=W53` | `FILE.53` | NBO coefficients in the NAO basis |
| `AONHO=W54` | `FILE.54` | NHO coefficients in the AO basis |
| `AONLMO=W55` | `FILE.55` | NLMO coefficients in the AO basis |
| `AOPNAO=W56` | `FILE.56` | Pre-orthogonal NAO coefficients in the AO basis |
| `NAONHO=W57` | `FILE.57` | NHO coefficients in the NAO basis |
| `NAONLMO=W58` | `FILE.58` | NLMO coefficients in the NAO basis |
| `NHONBO=W59` | `FILE.59` | NBO coefficients in the NHO basis |
| `NBONLMO=W60` | `FILE.60` | NLMO coefficients in the NBO basis |
| `NLMOMO=W61` | `FILE.61` | MO coefficients in the NLMO basis |
| `AOMO=W62` | `FILE.62` | Canonical MO coefficients in the AO basis |

`ARCHIVE` writes the `.47` file. `PRINT=3` includes NLMO and Wiberg output; the template also requests `E2PERT=0.0`. Its full NBO keylist is in [cov_one_job.py](../examples/nbo-one-job/cov_one_job.py).

A `.47` file does not necessarily contain Fock data. This template stops if the archive lacks overlap, density, canonical MO or Fock sections, because it prepares inputs for the energy views as well. COV can still use a smaller file set for views that do not need the missing data. Open-shell views need the corresponding spin data; the template does not choose restricted or unrestricted methods for the molecule.

## Keep and recover files

The original input, calculation log and checkpoint are useful for later exports. They are not all required to open an existing set of exported files. The canonical FCHK comes from a checkpoint retaining the canonical MOs; SaveNBOs or SaveNLMOs can replace those orbitals in a checkpoint.

| What is missing? | How it can be recovered |
| --- | --- |
| FCHK, with the original Gaussian CHK still available | Export with `formchk`, or let COV call the installed converter. |
| A report section or orbital matrix, with the needed data still in `.47` | Rerun GenNBO with the additional output options; use that run's report and matrices together. No new SCF is normally needed. |
| Fock or other data absent from the archive | Export from the original calculation state if the software supports it. GenNBO cannot supply data absent from its input. |
| `.covnbopkg` opening list | Open the folder or select the files directly. |
| Files belonging to different calculations | Find the matching wavefunction and NBO outputs. Renaming files does not change their contents. |
| Only coordinates remain | A new electronic-structure calculation is needed to obtain the wavefunction. |

SALCs are constructed in COV from the imported orbitals and molecular symmetry; they do not require another SCF calculation. Their shapes and energies have different data requirements, described in [Using NBO results](AOMO_NBO.md).

Matrix output and `.47` reanalysis are described in the [NBO 7 manual](https://nbo.chem.wisc.edu/nboman.pdf), sections B.2.4–B.2.6 and B.7.
