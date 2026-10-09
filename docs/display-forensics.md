# Inspect the actual COV display

This optional development tool records the native UI after `ImGui::Render()`.
It uses the same diagram, details window, input handling and OpenGL rendering as
the application. It does not recreate details from a separate chemistry report.
It is absent from the ordinary build and from normal picture exports.

Build with `COV_ENABLE_VALIDATION=ON`. Then, with Python 3:

```text
python tests/cov_display_forensics.py --exe <cov_validation.exe> --library <library-directory> --output <new-directory> --cases PG-001 PG-006 --all-nodes
```

Omit `--all-nodes` for an inventory. Use `--selected` with stable hit IDs for a
targeted replay. `--node-limit` limits collection explicitly; omitted objects
remain listed as not attempted. `--language` selects the real language control.
The default mesh resolution is 48 to make repeated selection affordable; use
`--resolution` to change it. These runs do not establish isosurface accuracy.

The runner creates a hidden native GLFW window, never moves the desktop cursor,
and does not focus another application. Each case is a separate process. A sweep
uses one process for all its objects. Nodes outside the scroll viewport are
reached through real wheel input. Overlapping bars use the application's chooser.
Navigation reveals the union of the selected bar and its label before clicking
the real bar. Wheel requests use whole input ticks so a half-pixel layout edge
cannot stall indefinitely after the UI rounds its scroll position.
The details window is scrolled from top to bottom, with overlapping captured
pages. No button result, selected orbital, scroll offset or display string is
forced by the collector.

Each case produces a primary `.display.json.gz` bundle. Repeated JSON containers
are stored once under content hashes and referenced by individual frames. The
Python `load_bundle()` helper restores the decoded files exactly and checks the
shared-object hashes. Uncompressed native files are development staging; `--zip`
optionally archives that staging as well and can be substantially larger. Normal
users exporting a picture receive none of these files.

## What a capture contains

- The actual loaded canonical path and attached NBO source. A failed load that
  leaves the previous molecule on screen cannot pass the identity check.
- The actual diagram snapshot: mode, filters, names, members, energies,
  occupations, node identities, bar and label rectangles, and visible edge data.
- The selected and rendered orbital identities, including set, source, spin and
  coefficient source, plus the current details route.
- Final rendered text runs decoded from font-atlas glyph quads. This includes
  directly drawn graph labels as well as ordinary widgets, Chinese, Greek and
  subscript glyphs. Separately formatted metadata is not substituted for them.
- Draw order, window rectangles, command/viewport clipping and conservative
  later-window overlap. An open details window is not proof it is visible.

`--capture` adds actual framebuffer BMPs; use it for a small number of visual
checks. It is off by default. No coefficient/Fock/density matrices, full input
analysis, volume texture dumps or per-frame history are requested. The native
`--validation-forensic` switch also bypasses the older large diagnostic writers.
Ordinary PNG/SVG export remains independent of this tool.

## Native plan additions

Use these only with `--validation-forensic`:

```text
inspect "current-view"
inspect-details "selected-details"
```

`inspect` saves one settled `.view.json`. `inspect-details` requires an open
details window and writes `selected-details-p000.view.json`, etc. It reports a
failure when scrolling is blocked or the page limit is reached. Native title-bar
targets use the outer window clip, while content targets use the content clip.
An `executed` pointer command alone is never counted as a successful inspection.

## Completeness and limits

Use `--expand-pi` to open an applicable interaction-details header through the
actual native input path before capturing its pages. `--preset research` or
`--preset full` selects the actual AO/MO view combo; the default preserves the
initial view. These options do not replace selection with an analysis export.
An absent header is recorded as inapplicable, while a header that remains closed
after an attempted expansion fails the requested inspection. Mode-network
headers retain their own open/closed targets for explicit expansion replays.

The collector also checks the frozen complete-NAO composition denominator,
agreement between selected details and graph objects, final displayed counts,
common-mode/source-edge identity and the separation of a multi-group network
from a two-endpoint relation. These checks supplement chemical controls; they
do not make a successful native session a scientific validation by itself.

The report distinguishes included, attempted, captured, failed and unvisited
objects. Its scope is the source nodes included in that initial view, not every
orbital in the source file or every view mode. Folded group headers are listed
separately. Capture checks include loaded-file identity, selection identity,
rendered-scene identity and actual clean details text on every page. It flags
duplicate displayed endpoint pairs and differing rendered bonding labels within
a displayed group, without declaring the chemistry correct or incorrect.

The decoder handles standard indexed ImGui font quads. Spaces and run boundaries
are inferred from geometry. Missing source characters cannot be reconstructed:
the actual fallback glyph is recorded. Truncation, ambiguous or unmatched glyphs,
unsupported callbacks and clipping are explicit. A later translucent window is
reported as possible coverage. Same-list shapes, the 3D scene and final pixel
legibility require a framebuffer check; the tool does not claim to solve OCR or
prove all pixels unobscured.

The automated checks cover real ImGui font rendering, CJK and scientific glyphs,
clipping, fallback glyphs, navigation, plan generation, wrong-file rejection,
selection/scene mismatch and covered details. Existing library files and ordinary
application exports are not modified by collection.

## Current-view data export

Picture export produces PNG and SVG only. Optional data export uses
`cov_aomo_objects_v5` CSV and `cov_aomo_unified_view_v5` JSON for the current
display objects. The CSV records display identity, original references, energy,
occupation and membership; the JSON supplies source mappings and shared tables.
The `object_table_file` value names the actual sibling CSV, and
`quantitative_links_ref` resolves within the JSON. A source/derived spatial
identity is distinct from its translated or view-numbered label.

Source projections are stored once and referenced by source identity instead of
being repeated for every MO/relation combination. This is a lossless schema
change, not an extra output threshold. A selected MO scope retains its original
global index. Expanded relation consumers use `routed_mo_relations`; absent
weights and occupied-space-mixing states remain explicit. RO display rows retain
both spin originals, the validated mean energy and summed occupation. Missing
energy is never serialized as zero.
