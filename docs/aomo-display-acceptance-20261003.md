# AOMO display and scientific consistency acceptance — 2026-10-03

This change repairs shared scientific/display paths rather than recognizing
molecule names or MO numbers. Source canonical coefficients, energies,
occupations and immutable orbital references are preserved. Overview filtering,
partner display and current-view numbering are derived state.

## Scientific contracts

- Nuclear geometry, density-verified electronic symmetry and local coordination
  symmetry are separate scopes. Actual retained operation matrices define the
  electronic naming scope. In the saved inputs, OLD-245/250 support Cs electronic
  scopes under C3v geometry; OLD-258 supports C4v under D4d geometry. An independent
  literal-density/analytic-overlap calculation agrees within 2e-14.
- Repeated multidimensional SALC spaces are resolved into complete copies with
  explicit source maps. The same transformation is applied to full Fock and
  density matrices and signed contributions. Copy membership does not by itself
  prove quantitative energy degeneracy. Approximate dominant names keep their
  approximation marker and do not acquire a false certified partner identity.
- Ordered occupied NBO donor to unoccupied acceptor evidence uses the actual
  same-spin operator and complete internal ligand pi/pi-star spaces. Printed E(2)
  count is not a vote. Two occupied endpoints remain occupied-space mixing.
- Weak subsidiary relations use frozen-operator energy sensitivity and principal
  angles, with complete-space, frontier and near-resonance protections. The
  `cov.pi-display.2026-10-03.v1` budgets (0.025 eV and sin-squared angle 0.02)
  are display calibration budgets, not universal chemical bond-energy cutoffs or
  self-consistent deletion energies. Unavailable operators cannot qualify.
- Whole-group bonding descriptions use a fixed structural scope and the
  spectrum of its metric-normalized Hermitian overlap operator. Rotating partners
  cannot flip their shared result. A mixed-sign spectrum is mixed; an unresolved
  numerical direction is not silently called nonbonding. This descriptor is not
  a bond energy or a statement about every local bond.
- Default structural bonds are separate from one localized Lewis assignment and
  from interaction/contact edges. Continuous Wiberg values are not individually
  rounded to integer bond multiplicities. Explicit Lewis display remains optional.

## Actual numerical and integration checks

| Scope | Result |
| --- | --- |
| Automated checks | 57/57 passed, including the saved FCHK regression |
| Naming corpus | 106/106 packages; 122 raw/spatial models; 3,880 certified multidimensional copies |
| Actual-input covariance | Five representative models; maximum space leakage 7.44e-9, group-weight change 1.78e-11, energy-trace change 1.26e-12 Ha |
| Repeated-copy mapping | Eight real packages, including E/T/G/H spaces, passed |
| Pi direction | Carbonyl and cyanide controls; nine actual E(2) print deletion/duplication runs preserve direction and source-group mapping |
| Weak-display calibration | 106 packages, 409 scopes, 24 eligible, six subsidiary scopes across five molecules qualify; main valence-d pi-star channels remain |
| Original conflicting groups | 53 source groups / 121 members: 33 positive and 20 negative whole-group spectra; no member contradiction |
| Previously missing bonding rows | 681 original columns: 326 positive, 345 negative, two numerical unresolved, eight single-atom not-applicable |
| Three-mode round trips | 106 packages; 16,665 full central and 18,507 full side objects agree with independently read source counts |
| Dense all-links drawing | 1,771,210 nonzero edges; clipping before tessellation reduced the measured dash workload from 970,832,802 to 4,300,355 |
| Actual same-window load/attach | Four packages, 32 stages; complete model restoration, correct canonical fallback and no stale NBO identity |
| Native four-language copy/export | 154 actual commands, zero failures; CSV/JSON/current diagram agree and ordinary picture output is only PNG/SVG |

OLD-245 received missing actual alpha/beta operators from its original frozen RO
density. No new SCF solution replaced its canonical input; all 21 current side
objects now have quantitative energy. OLD-273 Cr 4p and 4s display energies are
1.6133941228591366 and 0.5152546269578446 Ha, respectively, with both source spin
values and summed populations retained.

NaCl and MgF2 saved controls exercise higher-p inclusion from retained-MO
contribution rather than a metal-d prerequisite. Boric acid, an aluminate,
formate/formic acid and a silicic control distinguish equivalent delocalized
bonds from genuinely inequivalent bonds. Additional calculations are test inputs,
not hardcoded branches.

Overview folding removes corresponding unused complete AO/SALC groups together
with core/background MO groups, but preserves side spaces used by retained MOs.
Research/full modes restore original members. Hiding numbers and ignoring H
labels are presentation options only. An energy display group is labelled as an
energy group rather than being presented as a certified irreducible partner set.

Ordinary Windows checks exercised the movable details window above the main
panel, number controls, and startup from Chinese/space/chemical-symbol paths.
The application embeds a UTF-8 process manifest on MSVC Windows builds and reads
startup arguments through the Unicode API. Real Gaussian formchk conversions
returned byte-identical FCHK data without a visible command window.

## Reproducibility and limits

Focused numerical probes are built with the test targets; the independent
actual-display collector is described in [display-forensics.md](display-forensics.md).
The native collector records final rendered text and source/selection identity;
it does not claim a full pixel or isosurface proof for every orbital. Full-library
mode coverage is separate from selecting every object in every mode.

Raw library outputs are intentionally not bundled into picture exports or this
repository. PG-032's normalized full analysis JSON was 122,820,654 bytes, and a
single selected MO scope was 4,719,731 bytes; a full scientific analysis is still
larger than a current-view object export. All 6,480 selected relations round-trip
exactly, and smaller complete relation sets retain every field and ordering.

Passing these checks does not imply every source orbital has a strict symmetry
label. The 122 naming models retain 262 unclassified canonical and 1,227
unclassified side entries when the actual evidence is insufficient. Approximate
entries remain distinguishable. The two unresolved overlap columns are near
the numerical floor, and isolated atoms have no interatomic bonding scope.
No fake names, occupations, energies or donor directions fill those states.
