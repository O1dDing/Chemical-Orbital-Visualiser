# MO composition and NBO data

For file preparation and everyday use, see [Using NBO results](AOMO_NBO.md) and the [one-job template](NBO_ONE_JOB.md). This page describes the matrix relationships and exports used by the source implementation.

## Wavefunction and analysis identity

Use the canonical FCHK for the calculation step represented by the NBO archive. GenNBO produces a new analysis of that archived wavefunction; use its report and matrices together. Canonical and local orbitals retain separate identities, energies, spins and selections.

Association uses ordered atoms, geometry, basis functions, overlap, density and the available canonical coefficient columns. File dimensions alone do not establish a match. Missing beta canonical columns remain unavailable; an archive transform does not create a beta FCHK orbital. Source files and analysis segments are retained in exported records.

## Matrix relationships

For compatible AO representations, let A be AONAO, S the AO overlap, C the canonical coefficients and T NAOMO. The importer checks AᵀSA, AT = C and T = AᵀSC, plus available NAONBO/AONBO/NBOMO relationships. It reads matrix roles from their headings.

NAO group weights sum squared coefficients over the selected atoms or shells. Original nonorthogonal AO coefficient squares are coefficients, rather than atomic populations. Electron contributions use the corresponding occupation and density data. ECP-replaced core counts and explicit-electron populations remain separate.

A local basis can have fewer columns than the AO basis. Its available contributions and uncovered part remain separate; omitted directions are not filled or renormalized. Core, valence, Rydberg and angular labels come from the actual report. Unknown or missing values remain unavailable.

NAO and fragment side energies are operator expectation values in the molecular environment. Canonical MO energies, local Fock diagonals and E(2) interaction estimates are distinct quantities. SALC construction also uses the geometry, selected fragment and supported symmetry operations.

## Selection and export

The composition view uses the central diagram's actual MO members. Selecting an orbital for 3D inspection does not add it to that diagram. Atom and shell groups retain their membership; adding a ligand group is an explicit action.

NBO exports include report tables, `.nbo.json`, `.view.json`, `.view.svg` and `.view.png`. MO composition exports include `.focus.json`, `.focus.csv`, `.focus.groups.csv`, `.focus.svg` and `.focus.png`. The full numerical data retain the source identities, matrices and unavailable fields; view records retain the selected members and groups.

## Source checks

The local plan protocol and scene capture are described in [native checks](native-validation.md). NBO controls use `nbo.*` semantic targets; composition controls use `nbo.focus.*`.

The source tree includes `tests/nbo_integration_suite.py` and `tests/nbo_focus_suite.py` to read existing producer data and compare imports, exports and selected GPU textures. Their corresponding negative-input scripts check rejected or incomplete inputs. Independent reference libraries are supplied separately. These checks use fresh result directories and do not start chemistry jobs.

Method references are listed in [NBO and symmetry references](SCIENTIFIC-SOURCES.md).
