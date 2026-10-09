# Explore AO–MO and NBO relationships

**English** · [简体中文](AOMO_NBO.zh-CN.md) · [日本語](AOMO_NBO.ja.md) · [Français](AOMO_NBO.fr.md)

Chemical Orbital Visualiser (COV) brings orbital energies, occupations and composition into one view. In the v0.4.0 preview, you can open an existing Gaussian FCHK and its corresponding NBO results, then move from an MO, atom or bond to related orbitals and contributions. Inspect a selected orbital in 3D when you want to see its shape.

## Which files do I need?

COV reads results that already exist; it does not run Gaussian or NBO. Drop the files together, drop their directory, or enter a path. If the directory contains several calculations, choose the one you want. The advanced path controls help when the files are stored in different places.

| What you want to see | Files to provide |
| --- | --- |
| MO energies, occupations, energy diagram, Gaussian AO–MO coefficient links and 3D MO | Gaussian FCHK/FCH; Molden supports basic MO views. |
| Printed NPA charges, Wiberg indices and E(2) values | FCHK and an NBO report containing those sections. |
| NAO shapes, NPA charge and Wiberg colouring | Base files and AONAO, including the corresponding density. |
| MO composition in NAOs | Base files, AONAO and NAOMO; NAONBO provides another transform. |
| NBO shapes and canonical MO–NBO links | Base files and AONBO; MO links use the corresponding canonical coefficients in the FCHK. |
| NHO shapes and components | Base files, AONAO, AONHO and NAONHO; NHO–NBO links also use NHONBO and AONBO. |
| NLMO shapes, main NBO and tail | Base files, AONBO, AONLMO, NBONLMO and matching report entries. MO links use the corresponding canonical coefficients in the FCHK. |
| PNAO shapes | Base files, AONAO and AOPNAO. |
| SALC shapes | Base files and AONAO, with a geometry and selected orbital space that support the symmetry combination. |
| NAO and SALC energies calculated from the operator | The corresponding orbitals and Fock data, including the spin components needed for that calculation. |
| E(2) donor–acceptor orbitals | Base files, AONBO and report entries identifying the donor and acceptor. |

Base files mean the matching FCHK, NBO report and `.47` archive. The FCHK and NBO data must describe the same wavefunction. The report and orbital matrices must come from the same NBO analysis; a GenNBO reanalysis has its own report and matrices. The archive supplies the structure, basis, overlap, density and canonical MO coefficients. Spin colouring also needs the corresponding spin-resolved density and report data.

NBOMO and NLMOMO provide additional transforms when supplied. COV can also calculate MO links from the local orbital coefficients, overlap and canonical MO coefficients, so those two files are not always needed. If the FCHK lacks canonical coefficients for one spin, links to those MOs are unavailable.

Fock data are not required simply to read printed E(2) values or show the corresponding donor and acceptor orbitals. When supplied, they also allow COV to compare the printed coupling and energy difference with the matrix. NAO/NBO energies printed in the report can remain available without an archive Fock matrix. SALC shapes and their calculated energies have separate requirements.

Matrix file numbers come from the NBO output settings; COV identifies their contents from the matrix headings. The names `canonical.fchk` and `analysis.nbo` are used by our template, not required filenames. A `.covnbopkg` file lists the files to open; opening the directory or selecting the files also works.

For the calculation sequence, output options and missing files, see [Prepare calculation files](NBO_ONE_JOB.md). The current NBO association uses Gaussian FCHK/FCH data; Molden remains an option for ordinary MO views. If NBO data are missing, those ordinary MO views remain available.

## Follow an orbital through the diagram

Choose an MO in the browser or energy diagram. The diagram focuses on valence and nearby unoccupied levels; **All** in the browser shows the full imported orbital list. Click an MO to highlight its available links to Gaussian atomic orbitals (AO), natural atomic orbitals (NAO) and atom groups. Click an AO, NAO or link to inspect a contribution or find related MOs. You can select several terms and show their signed partial sum or overlay them in 3D. Folding groups tidies the diagram without removing their members.

NAOs form an orthogonal representation, so squared NAO coefficients can describe weights in that representation. Squared coefficients of the original, generally nonorthogonal Gaussian AOs are not atomic populations. When NAO or SALC energies are available, they can use the same numerical axis as the MOs. Their values are operator expectation values in the molecular environment; the central values are canonical MO energies. The illustrative side layout arranges orbitals without using their energies. A group of orbitals is a SALC only when the available data supports that symmetry meaning.

A local basis may contain fewer orbitals than the Gaussian AO basis. The diagram then shows the available contributions and the uncovered part, without rescaling the contributions to 100%.

**Overview** shows a compact view. **Research analysis** adds detail, and **Full basis** includes core and Rydberg orbitals. These presets change what is shown.

## Start from the molecule

Click an atom or bond in the 3D view to open its related values and orbital links. With the report and matching matrices available, atoms can be coloured by NPA charge or spin population; bonds can show Wiberg indices, bonding and antibonding orbitals, and coordination or multicentre relationships. A Wiberg index is a continuous value, not an integer bond order.

NHO views can show directional orbital lobes and their angular-momentum components. Select an E(2) interaction to see its donor and acceptor orbitals together; the E(2) value is a perturbation estimate, not a bond or reaction energy. NLMO views can separate the full orbital, its main NBO component and the remaining tail when the needed data is available.

A small contribution may be hidden by the current isosurface threshold. Use **Fit component** to adjust the display threshold or **Reveal bonds** to make the molecule easier to see. Keep the same threshold when comparing the sizes of different orbitals.

## Export

**Export images** saves PNG and SVG figures. **Analysis data (advanced) → Export analysis data** separately saves JSON/CSV and the associated analysis files. Interactive 3D exploration remains in COV.
