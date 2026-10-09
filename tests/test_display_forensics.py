"""Fixture tests for native plan contracts and honest display coverage."""
from __future__ import annotations

import json
import gzip
from pathlib import Path
import tempfile
import unittest

import cov_display_forensics as f


def node(index: int, *, group: bool = False, lane: int = 1, clipped: bool = False) -> dict:
    return {"id": f"canonical_mo:{index}", "hit_id": f"aomo.node.canonical_mo:{index}",
            "choice_hit_id": f"aomo.node.choice.canonical_mo:{index}", "canonical_index": index,
            "group_header": group, "lane": lane, "scroll_clipped": clipped,
            "clickable_in_viewport": not clipped, "display_group_id": "group-A",
            "logical_rect": [10, index * 40, 72, index * 40 + 22]}


def text_run(text: str, *, clipping: str = "unclipped", occlusion: str = "no_later_window_overlap",
             x: int = 10, y: int = 10, width: int = 150) -> dict:
    return {"text": text, "owner": "轨道详情###cov.orbital.details", "clipping": clipping,
            "occlusion": occlusion, "bounds": [x, y, x + width, y + 20], "baseline_origin_y": y}


def fixture_view(index: int, canonical: Path, *, details: bool = False, extra_trace: list | None = None) -> dict:
    active = {"kind": "canonical", "canonical_index": index}
    trace = [{"kind": "forensic.input", "data": {"canonical_path": canonical.as_posix(), "error": "", "nbo_source": "FILE.47"}},
             {"kind": "forensic.aomo", "data": {"focused_canonical_index": index, "active_view": active,
                                                  "selected_side_node_id": "", "nodes": [node(index)]}}]
    if details:
        trace.append({"kind": "forensic.diagram.details", "data": {"open": True, "content_visible": True,
                                                                    "active_view": active}})
    trace.extend(extra_trace or [])
    return {"state": {"applied_mo": index, "scene_matches_applied": True}, "draw_trace": trace,
            "rendered": {"status": "captured", "counters": {"ambiguous_glyphs": 0,
                         "unmatched_font_texture_quads": 0, "callbacks_not_decoded": 0, "invalid_index_groups": 0},
                         "text_runs": [text_run("Actual final text")], "glyphs": []}}


class PlanTests(unittest.TestCase):
    def test_paths_and_unicode_have_native_quoted_syntax(self):
        self.assertEqual(f.quoted(r"E:\样本\OLD-012"), '"E:/样本/OLD-012"')
        self.assertEqual(f.quoted('choice "A"'), '"choice \\"A\\""')
        with self.assertRaises(ValueError):
            f.quoted("bad\npath")

    def test_inventory_requests_real_language_input_and_low_scene_resolution(self):
        plan = f.inventory_plan(Path("E:/library/packages/OLD-012"))
        self.assertIn("scene 0.60 0.65 0.35 2.2 0.03 48", plan.commands)
        self.assertIn('drop "E:/library/packages/OLD-012"', plan.commands)
        start = plan.commands.index('click "language"')
        self.assertEqual(plan.commands[start:start + 4], ['click "language"', 'key "Home"',
                                                         'key "Down"', 'key "Enter"'])
        self.assertTrue(plan.text().startswith("COV_VALIDATION 1\n"))
        self.assertFalse(any("export" in command for command in plan.commands))

    def test_offscreen_nodes_are_attempted_and_group_headers_are_reported_separately(self):
        nodes = [node(1, clipped=True), node(2, group=True), node(3)]
        chosen = f.choose_nodes(nodes)
        self.assertEqual([n["canonical_index"] for n in chosen], [1, 3])
        plan = f.sweep_plan(Path("E:/packages/X"), chosen, language="current", inventory_nodes=nodes)
        self.assertIn('seek "aomo.node.canonical_mo:1.reveal"', plan.commands)
        for attempt in plan.attempts:
            commands = plan.commands[attempt["command_start"]:attempt["command_end"]]
            self.assertLess(next(i for i, c in enumerate(commands) if c.startswith("inspect ")),
                            commands.index('click "diagram.details"'))
            self.assertEqual(commands[-1], 'click "diagram.details.close.bottom"')
        coverage = f.summarize_coverage(nodes, [], False)
        self.assertEqual(coverage["included_source_nodes"], 2)
        self.assertEqual(coverage["included_group_headers"], 1)
        self.assertFalse(coverage["complete"])

    def test_selection_does_not_invent_absent_hit_ids(self):
        with self.assertRaisesRegex(ValueError, "absent"):
            f.choose_nodes([node(1)], ["aomo.node.not-in-inventory"])

    def test_overlap_chooser_uses_actual_inventory_geometry(self):
        a, b = node(1), node(2)
        b["logical_rect"] = a["logical_rect"][:]
        plan = f.sweep_plan(Path("E:/packages/X"), [a], language="current", inventory_nodes=[a, b])
        self.assertTrue(plan.attempts[0]["overlap_choice"])
        self.assertIn('click "aomo.node.choice.canonical_mo:1"', plan.commands)
        self.assertFalse(f.needs_choice(a, [a]))

    def test_artifact_names_avoid_source_number_assumptions_and_collisions(self):
        a, b = f.artifact_id("aomo.node:α"), f.artifact_id("aomo.node:β")
        self.assertNotEqual(a, b)
        self.assertRegex(a, r"^[A-Za-z0-9_-]+$")


class EvidenceTests(unittest.TestCase):
    def test_input_path_mismatch_is_rejected_despite_matching_mo_index(self):
        intended = Path("E:/library/A/canonical.fchk")
        view = fixture_view(1, Path("E:/library/B/canonical.fchk"))
        self.assertFalse(f.input_identity(view, intended)["matched"])
        self.assertTrue(f.selection_match(node(1), view))
        self.assertFalse(f.input_identity({"draw_trace": []}, intended)["matched"])

    def evaluate(self, directory: Path, *, failed_action: bool = False, missing_pages: bool = False,
                 wrong_input: bool = False, wrong_selection: bool = False) -> dict:
        canonical = Path("E:/library/X/canonical.fchk")
        plan = f.sweep_plan(Path("E:/library/X"), [node(1)], language="current")
        attempt = plan.attempts[0]
        selected_index = 2 if wrong_selection else 1
        loaded = Path("E:/library/OTHER/canonical.fchk") if wrong_input else canonical
        name = attempt["artifact_id"]
        f.write_json(directory / (name + "-selected.view.json"), fixture_view(selected_index, loaded))
        if not missing_pages:
            f.write_json(directory / (name + "-details-p000.view.json"), fixture_view(selected_index, loaded, details=True))
        actions = [{"command": i, "status": "executed", "detail": ""}
                   for i in range(attempt["command_start"], attempt["command_end"])]
        if failed_action:
            actions[0].update(status="failed", detail="target clipped or unreachable by wheel input")
        return f.evaluate_attempt(attempt, directory, actions, canonical)

    def test_complete_actual_capture_and_failed_seek_have_distinct_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good = root / "good"
            good.mkdir()
            success = self.evaluate(good)
            self.assertEqual(success["status"], "captured")
            self.assertTrue(f.summarize_coverage([node(1)], [success], True)["complete"])
            bad = root / "bad"
            bad.mkdir()
            failure = self.evaluate(bad, failed_action=True)
            self.assertEqual(failure["status"], "failed")
            coverage = f.summarize_coverage([node(1)], [failure], True)
            self.assertEqual(coverage["captured_fraction"], 0)
            self.assertFalse(coverage["complete"])

    def test_missing_pages_stale_input_and_wrong_selection_never_count_as_captured(self):
        for flag in ("missing_pages", "wrong_input", "wrong_selection"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as tmp:
                result = self.evaluate(Path(tmp), **{flag: True})
                self.assertEqual(result["status"], "failed")
                self.assertTrue(result["failures"])

    def test_node_limit_does_not_report_full_coverage(self):
        nodes = [node(1), node(2), node(3)]
        chosen = f.choose_nodes(nodes, limit=1)
        attempts = [{"hit_id": chosen[0]["hit_id"], "status": "captured"}]
        coverage = f.summarize_coverage(nodes, attempts, True)
        self.assertEqual(coverage["captured_fraction"], 1 / 3)
        self.assertEqual(len(coverage["not_attempted"]), 2)
        self.assertFalse(coverage["complete"])

    def test_covered_details_and_decoder_gaps_do_not_claim_complete_text(self):
        canonical = Path("E:/library/X/canonical.fchk")
        view = fixture_view(1, canonical, details=True)
        view["rendered"]["text_runs"] = [text_run("real but covered", occlusion="possible_later_window_cover")]
        self.assertFalse(f.rendered_completeness(view, True)["complete"])
        view["rendered"]["text_runs"] = [text_run("clean")]
        for field in ("unmatched_font_texture_quads", "callbacks_not_decoded", "ambiguous_glyphs"):
            with self.subTest(field=field):
                view["rendered"]["counters"][field] = 1
                self.assertFalse(f.rendered_completeness(view, True)["complete"])
                view["rendered"]["counters"][field] = 0
        view["rendered"]["status"] = "truncated"
        self.assertFalse(f.rendered_completeness(view, True)["complete"])

    def test_scene_identity_mismatch_is_a_failed_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.evaluate(root)
            selected = root / result["selected_capture"]
            value = json.loads(selected.read_text(encoding="utf-8"))
            value["state"]["scene_matches_applied"] = False
            f.write_json(selected, value)
            plan = f.sweep_plan(Path("E:/library/X"), [node(1)], language="current")
            attempt = plan.attempts[0]
            actions = [{"command": i, "status": "executed"} for i in range(attempt["command_start"], attempt["command_end"])]
            revised = f.evaluate_attempt(attempt, root, actions, Path("E:/library/X/canonical.fchk"))
            self.assertEqual(revised["status"], "failed")
            self.assertTrue(any("scene" in error for error in revised["failures"]))

    def test_clean_transcript_preserves_covered_or_clipped_text_separately(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            view = fixture_view(1, Path("E:/canonical.fchk"))
            view["rendered"]["text_runs"] = [text_run("清晰文本"),
                text_run("剪裁文本", clipping="command_partial"),
                text_run("遮挡文本", occlusion="possible_later_window_cover")]
            f.write_json(root / "inventory.view.json", view)
            f.write_rendered_report(root, root / "rendered.txt")
            clean = (root / "rendered.txt").read_text(encoding="utf-8")
            qualified = (root / "rendered-qualified.txt").read_text(encoding="utf-8")
            self.assertIn("清晰文本", clean)
            self.assertNotIn("剪裁文本", clean)
            self.assertIn("剪裁文本", qualified)
            self.assertIn("遮挡文本", qualified)

    def test_duplicate_gap_endpoint_observation_is_per_actual_page(self):
        gap = {"kind": "details.energy-gap", "data": {"gap_kind": "pi_partner", "lower_orbitals": [1, 2],
                                                       "upper_orbitals": [3], "label": "actual label"}}
        self.assertEqual(f.duplicate_gap_pairs({"draw_trace": [gap]}), [])
        duplicate = f.duplicate_gap_pairs({"draw_trace": [gap, gap]})
        self.assertEqual(duplicate[0]["record_count"], 2)
        self.assertEqual(duplicate[0]["lower_orbitals"], [1, 2])

    def test_group_bonding_difference_uses_clean_final_runs_not_semantic_text(self):
        view = {"rendered": {"text_runs": [text_run("各原子对的综合成键作用: 反键")]},
                "draw_trace": [{"kind": "draw.text", "data": {"value": "fabricated bonding"}}]}
        observed = f.observed_bonding_labels(view)
        self.assertEqual(observed[0]["value"], "反键")
        attempts = [{"hit_id": node(1)["hit_id"], "status": "captured", "observed_bonding_labels": observed},
                    {"hit_id": node(2)["hit_id"], "status": "captured", "observed_bonding_labels": [{"value": "成键"}]}]
        flags = f.group_label_differences([node(1), node(2)], attempts)
        self.assertEqual(flags[0]["group_id"], "group-A")
        self.assertIn("not assessed", flags[0]["scientific_verdict"])


class CompactBundleTests(unittest.TestCase):
    def test_roundtrip_preserves_decoded_views_bounds_errors_and_reserved_keys(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native = root / "native"
            native.mkdir()
            shared = {"nodes": [{"label": "真实显示αβ" * 80, "bounds": [-12.5, 10, 800.25, 44],
                                  "errors": [None, "遮挡或未完成"], "$ref": "source key", "r": 0,
                                  "literal_array": [2, 0]}]}
            original = {"a.view.json": {"state": {"frame": 1}, "actual": shared},
                        "b.view.json": {"state": {"frame": 2}, "actual": shared}}
            for name, value in original.items():
                f.write_json(native / name, value)
            f.write_json(native / "session.json", {"failed_commands": 1, "capture_status": "failed"})
            (native / "actions.jsonl").write_text('{"status":"failed","detail":"target covered"}\n', encoding="utf-8")
            (native / "plan.txt").write_text('COV_VALIDATION 1\ninspect "a"\n', encoding="utf-8")
            bundle = root / "case.display.json.gz"
            stats = f.compact_bundle(native, bundle)
            decoded = f.load_bundle(bundle)
            for name, value in original.items():
                self.assertEqual(decoded[name], value)
                self.assertTrue((native / name).is_file())
            self.assertEqual(decoded["session.json"]["failed_commands"], 1)
            self.assertEqual(decoded["actions.jsonl"][0]["detail"], "target covered")
            self.assertEqual(decoded["plan.txt"], 'COV_VALIDATION 1\ninspect "a"\n')
            self.assertGreater(stats["shared_objects"], 0)
            self.assertEqual(stats["view_files"], 2)
            self.assertEqual(stats["compressed_bytes"], bundle.stat().st_size)
            decoded["a.view.json"]["actual"]["nodes"][0]["label"] = "edited locally"
            self.assertEqual(decoded["b.view.json"], original["b.view.json"])

    def test_corrupt_shared_object_hash_is_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native = root / "native"
            native.mkdir()
            for name in ("a.view.json", "b.view.json"):
                f.write_json(native / name, {"actual": {"label": "long actual text " * 80}})
            bundle = root / "case.display.json.gz"
            f.compact_bundle(native, bundle)
            with gzip.open(bundle, "rt", encoding="utf-8") as stream:
                value = json.load(stream)
            value["objects"][0]["sha256"] = "0" * 64
            corrupt = root / "corrupt.display.json.gz"
            with gzip.open(corrupt, "wt", encoding="utf-8") as stream:
                json.dump(value, stream)
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                f.load_bundle(corrupt)


class ScientificDisplayContractTests(unittest.TestCase):
    def test_research_uses_actual_combo_in_inventory_and_sweep(self):
        for plan in (f.inventory_plan(Path("E:/packages/X"), preset="research"),
                     f.sweep_plan(Path("E:/packages/X"), [node(1)], preset="research")):
            self.assertLess(plan.commands.index('click "aomo.preset"'),
                            plan.commands.index('click "aomo.preset.research"'))
            self.assertLess(plan.commands.index('click "aomo.preset.research"'),
                            plan.commands.index('inspect "inventory"'))
    def test_expansion_precedes_actual_capture(self):
        plan = f.sweep_plan(Path("E:/packages/X"), [node(1)], expand_pi=True)
        expand = plan.commands.index('expand "details.pi.toggle"')
        capture = next(i for i, c in enumerate(plan.commands) if c.startswith('inspect-details '))
        self.assertLess(expand, capture)
        self.assertTrue(plan.attempts[0]['pi_expansion_requested'])

    def test_conditional_percentages_cannot_pass_complete_norm_check(self):
        ledger = dict.fromkeys(f.COMPOSITION_BUCKETS, 0.0)
        ledger.update(available=True, complete=True, weight_sum=1.0,
                      centre_current_p=.002, ligand_valence=.035)
        view = {"draw_trace": [{"kind": "details.composition", "data": ledger}]}
        self.assertEqual(f.display_semantics(view)["status"], "failed")
        ledger['ligand_other'] = .963
        self.assertEqual(f.display_semantics(view)["status"], "passed")

    def test_direction_needs_common_mode_and_same_source_edge(self):
        mode = {"verified": True, "direction_verified": True,
                "two_endpoint_relation": True,
                "direction": "centre_to_ligand", "shared_mode_ids": ["mode-a"],
                "shared_fragment_contraction_norm_hartree": .04, "matched_edge_ids": []}
        view = {"draw_trace": [{"kind": "details.energy-gap", "data": {"mode": mode}}]}
        self.assertEqual(f.display_semantics(view)["status"], "failed")
        mode['matched_edge_ids'] = ['edge-1-to-2']
        self.assertEqual(f.display_semantics(view)["status"], "passed")
        view['targets'] = [{'id': 'details.pi.toggle.closed'}]
        self.assertEqual(f.display_semantics(view, True)["status"], "failed")

    def test_network_cannot_be_repeated_or_turned_into_pair(self):
        network = dict(verified=True, mode_id="mode-a", channel_id="source-a",
                       nodes=[{"members": [2]}, {"members": [5]}, {"members": [8]}])
        record = {"kind": "details.pi-network", "data": {"network": network}}
        self.assertEqual(f.display_semantics({"draw_trace": [record]})["status"], "passed")
        self.assertEqual(f.display_semantics({"draw_trace": [record, record]})["status"], "failed")
        mode = dict(verified=True, shared_mode_ids=["mode-a"],
                    shared_fragment_contraction_norm_hartree=.02, two_endpoint_relation=False)
        self.assertEqual(f.display_semantics({"draw_trace": [
            {"kind": "details.energy-gap", "data": {"mode": mode}}]})["status"], "failed")

    def test_joined_backbond_requires_occupied_metal_and_applicable_group(self):
        assessment = dict(available=True, applicable=True, members=[51, 52, 53],
                          occupied_metal_backbond_supported=False)
        view = {"draw_trace": [{"kind": "details.pi.joined-summary",
                               "data": {"label": "π 回馈", "assessment": assessment}}]}
        self.assertEqual(f.display_semantics(view)["status"], "failed")
        assessment["occupied_metal_backbond_supported"] = True
        self.assertEqual(f.display_semantics(view)["status"], "passed")
        assessment["applicable"] = False
        self.assertEqual(f.display_semantics(view)["status"], "failed")

    def test_disconnected_atom_bucket_closes_without_relabelling_as_ligand(self):
        ledger = dict.fromkeys(f.COMPOSITION_BUCKETS, 0.0)
        ledger.update(available=True, complete=True, weight_sum=1.0,
                      centre_current_s=.3, ligand_valence=.4, other_atoms=.3)
        view = {"draw_trace": [{"kind": "details.composition", "data": ledger}]}
        self.assertEqual(f.display_semantics(view)["status"], "passed")
        ledger["ligand_valence"] += .3
        self.assertEqual(f.display_semantics(view)["status"], "failed")

    def test_explicit_source_roundoff_is_not_double_counted_by_checker(self):
        ledger = dict.fromkeys(f.COMPOSITION_BUCKETS, 0.0)
        ledger.update(available=True, complete=True, weight_sum=1-2.34e-8,
                      normalization_error=2.34e-8, ligand_valence=1-2.34e-8,
                      unresolved=2.34e-8)
        view = {"draw_trace": [{"kind": "details.composition", "data": ledger}]}
        self.assertEqual(f.display_semantics(view)["status"], "passed")
        ledger['unresolved'] = .01
        self.assertEqual(f.display_semantics(view)["status"], "failed")


if __name__ == "__main__":
    unittest.main()
