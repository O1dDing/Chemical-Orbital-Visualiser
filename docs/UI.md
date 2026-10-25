# Using Chemical Orbital Visualiser

**English** · [简体中文](UI.zh-CN.md) · [日本語](UI.ja.md) · [Français](UI.fr.md)

This guide covers COV v0.4.0. The [download instructions](../README.md#download) cover the published packages.

Chemical Orbital Visualiser (COV) brings orbital energies, occupations and relationships into one view. Open an existing Gaussian FCHK/FCH or Molden file, choose an orbital from the list, compare levels in the energy diagram, and export a figure or data table. The molecule and selected orbital can also be inspected in 3D.

## Open a calculation

Use **Open File**, drag a file into the window, or enter its path. COV reads calculation results; it does not run Gaussian. A Gaussian CHK file can be converted locally if `formchk` is installed. If COV cannot find the converter, create an FCHK file with `formchk` first and open that file instead.

If COV cannot find `formchk`, set `COV_FORMCHK` to its full executable path and restart COV.

Opening a calculation folder loads its wavefunction and NBO files together. If the folder contains several calculations, choose the one to open. The current input limit is 100 atoms per molecule.

## Browse orbitals and energies

Search the orbital list, jump to HOMO or LUMO, or use the core, valence and virtual filters. **All** returns to the complete imported list; filtering does not remove orbitals from the calculation. Details for the selected orbital show its energy, occupation, spin and any available symmetry or chemical labels.

Choose the energy unit that suits your work: Ha, eV, J/mol, kJ/mol, cal/mol or kcal/mol. The original energy is kept in Hartree; changing the displayed unit does not change the calculation.

Orbital names such as `2a₁` appear wherever the symmetry and occurrence order are available. Otherwise, COV shows the original MO number and spin. Both names and source numbers are searchable. Levels with close energies may be grouped at their mean energy; each member keeps its own energy and identity in the details. The grouping tolerance starts at `1e-5 Ha` and can be adjusted.

## Compare levels in the diagram

The central diagram shows occupied valence levels and nearby unoccupied levels, with electron occupations. Deep core and distant virtual levels can be hidden to keep the figure readable; they remain in the imported orbital list. Selecting an orbital for 3D inspection does not add it to or reorder the compact diagram.

The energy axis can be linear or nonlinear. The nonlinear view spreads crowded levels apart and is labelled **Nonlinear energy axis**, with energy values and units on the ticks. Grouped levels use the members’ mean energy.

**Export images** saves PNG and SVG figures. **Analysis data (advanced) → Export analysis data** separately saves JSON and CSV data. Figures use the displayed orbital names and the current interface language. The data includes orbital numbers, energies, occupations, spin, grouping and which levels appeared in the figure.

## Inspect the molecule and orbital

Rotate the scene and adjust atom, bond, molecule and orbital visibility. The default **Enhanced ball-and-stick** view emphasizes the framework; **Stick + delocalisation** reduces atom clutter. Hydrogen visibility can be changed separately. The positive and negative phases of an orbital surface use different colours.

Some bond and contact styles depend on information available in the calculation. A drawn line or colour is a display aid, not a measured bond order by itself.

The details window can be moved, resized and reopened. The interface is available in English, 简体中文, 日本語 and Français.

For AO–MO and NBO views, see the [AO–MO / NBO guide](AOMO_NBO.md).
