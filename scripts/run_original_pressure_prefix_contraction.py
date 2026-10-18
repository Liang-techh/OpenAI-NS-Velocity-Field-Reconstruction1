"""Build and inspect the saved-data pressure-prefix contraction packet.

The runner authenticates the immutable parent/native packets before decoding
them.  ``--write-report`` and ``--inspect-report`` are data-only operations;
the two numerical kernels are imported only by ``--compute-new``.  Every new
component is checkpointed before it is inspected so an interrupted run cannot
silently replay or overwrite saved arithmetic.
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
RUNNER = "scripts/run_original_pressure_prefix_contraction.py"
PARENT = "lei_ren_part1_paper_compliant_current_original_pressure_window_transport.json"
PARENT_BYTES = 832779
PARENT_SHA = "e425bbdff9f36ffb8c69052b377ec15440752cdaac8ec0919cb251ecd1b442e3"
NATIVE = "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
NATIVE_BYTES = 42177499
NATIVE_SHA = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
KERNEL = "src/openai_ns_reconstruction/original_pressure_prefix_log_transport.py"
KERNEL_SHA = "a4cdbef23e3435abd51f317e6fea8e0b1be6df308796a6970b0cf259d41d237b"
BOUNDARY = "src/openai_ns_reconstruction/original_first_switch_boundary_enclosure.py"
FIXTURE = "original_pressure_prefix_log_transport_fixture.json"
FIXTURE_BYTES = 36075
FIXTURE_SHA = "361cfa3b063fb1ee69fa85eca09156b04468b166fa939282e0ecf6f9a301c96e"
FIXTURE_RECORDED_SHA = "2eb16efa3ce8fa2a371f5906223bfbba75e252418a3d5bb3e8043d0dcceea20e"
BOUNDARY_EVIDENCE = {
    "experiments/root_st073/original_first_switch_boundary_initial_unaccepted_fixture.json":
        (22358, "9b91add827778ab6608e99214e67f2951e407e0a7fbc84f39fd7a97a275e4625"),
    "experiments/root_st073/original_first_switch_boundary_symmetric_error_regression.json":
        (52923, "b2cbee3e19ae08d96d6cfeb6f06221fa7deaf3a0ab5bd0eb0ccefadefeae22f4"),
    "experiments/root_st073/original_first_switch_boundary_initial_unaccepted_source.txt":
        (19291, "850b942665363d8a239120c86fb56d8a92521ec10fdfed725315a7be56147dd2"),
    "experiments/root_st073/original_first_switch_boundary_root_review.json":
        (3560, "8612a4b4e7db1e922edc53a6e365c4d68870536f94828e989882dd29f21745c2"),
}
REPORT_NAME = "lei_ren_part1_paper_compliant_current_original_pressure_prefix_contraction.json"
REPORT = BASE / REPORT_NAME
CHECKPOINT_DIR = ROOT / "tmp/original_pressure_prefix_contraction"
CHECKPOINTS = {
    "source": CHECKPOINT_DIR / "new_source.json",
    "boundary": CHECKPOINT_DIR / "new_first_B.json",
    "transport": CHECKPOINT_DIR / "new_transport.json",
    "join": CHECKPOINT_DIR / "new_join.json",
}
CHARTS = ("first_switch", "second_switch", "post_power")
ORDERS = (0, 1)
ZERO = (0, 0, 0, 0)
REMOVED = ((141114, 184687), (141115, 184688))
OPAQUE = ((183935, 184677), (183936, 184678))
FALSE_SCOPE = (
    "full_pressure_numeric", "complete_pressure_numeric", "global_closure",
    "global_five_moment_closure", "full_goal_complete",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return digest(path.read_bytes())


def _safe_relative(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute() and ".." not in path.parts,
            "unsafe relative path: " + name)
    return path


def _input_path(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute(), "absolute inherited input path: " + name)
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py",
                "unexpected inherited traversal: " + name)
        return (BASE / path).resolve()
    if path.parts and path.parts[0] == "experiments":
        return (ROOT / path).resolve()
    if path.parts and path.parts[0] in ("src", "scripts", "tmp"):
        return (ROOT / path).resolve()
    return (BASE / path.name).resolve()


def _implementation_path(name: str) -> Path:
    return (ROOT / _safe_relative(name)).resolve()


def _endpoint(value: object, label: str = "endpoint") -> tuple[int, int, int, int]:
    require(isinstance(value, (list, tuple)) and len(value) == 4,
            label + " must be a four-component exact MPF tuple")
    require(all(type(part) is int for part in value), label + " is not integer")
    result = tuple(value)
    sign, mantissa, exponent, bits = result
    require(sign in (0, 1) and mantissa >= 0 and bits == mantissa.bit_length(),
            label + " is noncanonical")
    require((mantissa == 0 and result == ZERO) or (mantissa > 0 and mantissa & 1),
            label + " has a non-normalized mantissa")
    return result


def _cmp_abs(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    lm, rm = left[1], right[1]
    if lm == 0 or rm == 0:
        return (lm > 0) - (rm > 0)
    lt, rt = lm.bit_length() + left[2], rm.bit_length() + right[2]
    if lt != rt:
        return (lt > rt) - (lt < rt)
    if left[2] >= right[2]:
        lv, rv = lm << (left[2] - right[2]), rm
    else:
        lv, rv = lm, rm << (right[2] - left[2])
    return (lv > rv) - (lv < rv)


def _cmp_exact(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    if left[1] == 0:
        return 0 if right[1] == 0 else (1 if right[0] else -1)
    if right[1] == 0:
        return -1 if left[0] else 1
    if left[0] != right[0]:
        return -1 if left[0] else 1
    result = _cmp_abs(left, right)
    return -result if left[0] else result


def _record(value: object, label: str = "record") -> dict[str, object]:
    require(isinstance(value, Mapping), label + " must be an object")
    lower = _endpoint(value.get("lower_exact_mpf_tuple"), label + ".lower")
    upper = _endpoint(value.get("upper_exact_mpf_tuple"), label + ".upper")
    require(_cmp_exact(lower, upper) <= 0, label + " is reversed")
    require(value.get("directed") is True, label + " must be directed")
    return {"lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper), "directed": True}


def _record_equal(left: object, right: object) -> bool:
    try:
        return _record(left) == _record(right)
    except (TypeError, ValueError):
        return False


def _records_count(value: object) -> int:
    if isinstance(value, Mapping):
        count = 0
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            _record(value)
            count = 1
        return count + sum(_records_count(item) for item in value.values())
    if isinstance(value, list):
        return sum(_records_count(item) for item in value)
    return 0


def _row_record(row: object, label: str, order: int | None = None) -> dict[str, object]:
    require(isinstance(row, Mapping), label + " row is not an object")
    if order is not None:
        require(row.get("ordinary_Z_order") == order, label + " order changed")
    return _record(row.get("value", row), label + ".value")


def _pair(rows: object, label: str) -> list[dict[str, object]]:
    require(isinstance(rows, list) and len(rows) >= 2, label + " must contain Z0/Z1 rows")
    return [_row_record(rows[order], label, order) for order in ORDERS]


def _authenticated_parent_files() -> tuple[dict[str, dict], dict[str, bytes]]:
    expected = {PARENT: (PARENT_BYTES, PARENT_SHA), NATIVE: (NATIVE_BYTES, NATIVE_SHA)}
    raws = {}
    for name, (size, expected_sha) in expected.items():
        raw = (BASE / name).read_bytes()
        require(len(raw) == size and digest(raw) == expected_sha,
                "authenticated report changed before decode: " + name)
        raws[name] = raw
    values = {name: json.loads(raw) for name, raw in raws.items()}
    require(all(isinstance(value, dict) for value in values.values()),
            "authenticated report is not an object")
    for report in values.values():
        for field, resolver in (("input_hashes", _input_path),
                                ("implementation_hashes", _implementation_path)):
            pins = report.get(field)
            require(isinstance(pins, dict), "missing " + field)
            for name, expected_sha in pins.items():
                require(isinstance(name, str) and isinstance(expected_sha, str),
                        "malformed " + field + " pin")
                path = resolver(name)
                require(path.is_file() and sha(path) == expected_sha,
                        field + " changed: " + name)
    return values, raws


def _preflight_file(relative: str, expected_sha: str | None = None) -> dict[str, object]:
    path = ROOT / _safe_relative(relative)
    raw = path.read_bytes()
    require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
            "LF/EOF preflight failed: " + relative)
    tree = ast.parse(raw.decode("utf-8"), filename=relative)
    actual = digest(raw)
    if expected_sha is not None:
        require(actual == expected_sha, "reviewed source changed: " + relative)
    return {"path": relative, "bytes": len(raw), "sha256": actual,
            "ast_parse": True, "source_executed": False,
            "ast_nodes": len(list(ast.walk(tree)))}


def _fixture_audit() -> dict[str, object]:
    path = BASE / FIXTURE
    raw = path.read_bytes()
    require(len(raw) == FIXTURE_BYTES and digest(raw) == FIXTURE_SHA,
            "saved prefix transport fixture changed")
    value = json.loads(raw)
    require(isinstance(value, dict) and value.get("passed") is True,
            "saved prefix transport fixture did not pass")
    require(value.get("fixture_sha256") == FIXTURE_RECORDED_SHA
            and value.get("module_sha256") == KERNEL_SHA,
            "prefix fixture/kernel pin mismatch")
    checks = value.get("checks")
    require(isinstance(checks, dict) and checks and all(item is True for item in checks.values()),
            "prefix fixture checks are incomplete")
    return {"path": "experiments/root_st073/" + FIXTURE, "bytes": FIXTURE_BYTES,
            "sha256": FIXTURE_SHA, "passed": True, "module_sha256": KERNEL_SHA,
            "source_executed": False, "checks": copy.deepcopy(checks)}


def _boundary_evidence_audit() -> dict[str, object]:
    records = {}
    decoded = {}
    for name, (size, expected_sha) in BOUNDARY_EVIDENCE.items():
        path = ROOT / name
        raw = path.read_bytes()
        require(len(raw) == size and digest(raw) == expected_sha,
                "boundary review evidence changed: " + name)
        records[name] = {"bytes": size, "sha256": expected_sha,
                         "source_executed": False}
        if name.endswith(".json"):
            decoded[name] = json.loads(raw)
    initial = decoded[next(name for name in decoded if name.endswith("initial_unaccepted_fixture.json"))]
    regression = decoded[next(name for name in decoded if name.endswith("symmetric_error_regression.json"))]
    review = decoded[next(name for name in decoded if name.endswith("root_review.json"))]
    source_path = ROOT / "experiments/root_st073/original_first_switch_boundary_initial_unaccepted_source.txt"
    source_raw = source_path.read_bytes()
    require(b"\r" not in source_raw and source_raw.endswith(b"\n")
            and len(list(ast.walk(ast.parse(source_raw.decode("utf-8"))))) > 0,
            "initial boundary source AST premise changed")
    require(initial.get("passed") is True and initial.get("module_sha256") !=
            digest((ROOT / BOUNDARY).read_bytes()),
            "initial boundary fixture acceptance receipt changed")
    require(regression.get("passed") is True
            and regression.get("module_sha256") == digest((ROOT / BOUNDARY).read_bytes())
            and regression.get("preflight_before_new_math", {}).get("passed") is True
            and regression.get("regression_contract", {}).get("new_input_cases_only") is True
            and regression.get("regression_contract", {}).get("both_negative_and_positive_errors_checked_by_exact_tuples") is True
            and regression.get("regression_contract", {}).get("actual_saved_data_produced") is False,
            "boundary regression receipt changed")
    require(review.get("passed") is True
            and review.get("final_module_sha256") == digest((ROOT / BOUNDARY).read_bytes())
            and review.get("initial_fixture_accepted") is False
            and review.get("new_regression_passed") is True
            and review.get("old_fixture_replayed") is False
            and review.get("actual_saved_data_produced") is False
            and review.get("changed_function_names") == ["bound_boundary"]
            and review.get("original_exact_proof_reused", {}).get("passed") is True
            and isinstance(review.get("unchanged_function_AST_sha256"), dict)
            and review.get("unchanged_function_AST_sha256"),
            "boundary root review receipt changed")
    return {"files": records, "review_schema": review.get("schema"),
            "final_module_sha256": review.get("final_module_sha256"),
            "initial_fixture_accepted": False, "new_regression_passed": True,
            "old_fixture_replayed": False, "actual_saved_data_produced": False,
            "changed_function_names": ["bound_boundary"],
            "original_exact_proof_reused": True}


def preflight() -> dict[str, object]:
    runner = _preflight_file(RUNNER)
    kernel = _preflight_file(KERNEL, KERNEL_SHA)
    boundary = _preflight_file(BOUNDARY)
    fixture = _fixture_audit()
    evidence = _boundary_evidence_audit()
    values, raws = _authenticated_parent_files()
    return {"runner": runner, "kernel": kernel, "boundary": boundary,
            "fixture": fixture, "boundary_evidence": evidence,
            "parent": values[PARENT], "native": values[NATIVE],
            "parent_raw": raws[PARENT], "native_raw": raws[NATIVE]}


def _validate_source_bindings(parent: Mapping[str, object], native: Mapping[str, object]) -> None:
    seam = native.get("actual_macro_switch_RQy_seam")
    require(isinstance(seam, Mapping) and seam.get("passed") is True,
            "RQy seam packet did not pass")
    guards = seam.get("source_guards", {})
    require(isinstance(guards, Mapping) and guards.get("passed") is True
            and guards.get("loop_common_N") is True
            and guards.get("native_J_distinct") == {"frozen_macro": 66987, "first_switch": 9}
            and guards.get("common_N") == 18
            and guards.get("nonreset_M") == {"frozen_macro": True, "first_switch": True},
            "RQy seam controls changed")
    proof = seam.get("proof_rows")
    require(isinstance(proof, list) and len(proof) == 2, "RQy proof row count changed")
    for order, row in enumerate(proof):
        require(row.get("ordinary_Z_order") == order
                and row.get("frozen_macro_upper_m1") == 141114 + order
                and row.get("first_switch_lower_s0") == 184687 + order,
                "RQy proof handles changed")
        identity = row.get("source_identity")
        require(identity == {"common_physical_radius_cutpoint": True,
                             "common_physical_Q_y_source": True,
                             "same_N": True, "nonreset_M": True,
                             "distinct_native_J": True, "passed": True},
                "RQy proof identity changed")
    symbolic = seam.get("symbolic_controls")
    require(isinstance(symbolic, Mapping) and symbolic.get("passed") is True
            and symbolic.get("ordinary_Z_orders") == [0, 1]
            and all(row.get("passed") is True and row.get("difference") == "0"
                    for row in symbolic.get("rows", [])),
            "RQy symbolic controls changed")

    source = native.get("actual_restoration_pressure_source", {})
    graph = source.get("graph_patch", {})
    require(graph.get("base_node_count") == 185704
            and graph.get("base_prefix_unchanged") is True
            and graph.get("base_canonical_nodes_sha256") ==
            "2d6efeb13efe91dd95d4f24fa03eb4626c7685d04e5f29f2bdc6532e101994e3"
            and graph.get("full_node_count") == 186086
            and graph.get("full_canonical_nodes_sha256") ==
            "dfbba64e1f106fa42c233af28a096f8dfad1e35580163c54cad9eb523d21271e",
            "restoration source graph patch changed")
    require(graph.get("segment_starts") == {
        "q4": 0, "accepted_merged": 183700, "accepted_seam": 184687,
        "accepted_first": 184703, "accepted_post": 184898,
        "accepted_long_reference": 185308, "accepted_correlation": 185690,
        "new_suffix": 185704}, "source segment starts changed")
    corr = native.get("actual_restoration_correlated_pressure", {})
    cgraph = corr.get("graph_patch", {})
    require(cgraph.get("base_node_count") == 186086
            and cgraph.get("base_canonical_nodes_sha256") == graph.get("full_canonical_nodes_sha256")
            and cgraph.get("full_node_count") == 186100
            and cgraph.get("full_canonical_nodes_sha256") ==
            "83b8594b5b775520cd38c9da1e6633e800ebf170dfb48b6b12eb34f903514fc0",
            "correlated source graph patch changed")
    first_source = native.get("actual_macro_first_switch_pressure_source", {}).get("first_switch_source", {})
    b_source = first_source.get("B_source", {})
    require(b_source.get("passed") is True
            and b_source.get("rows") == [184811, 184812]
            and b_source.get("formula") ==
            "B=[R^2*Q/L+R^2*Q^2/2-2*Z^2*R^2*U^2/L^2]_{s=0}^{s=1}"
            and b_source.get("physical_radius_power") == 2
            and b_source.get("physical_Jacobian_applied") is False
            and b_source.get("endpoint_measure_has_no_J") is True
            and b_source.get("ordinary_Z_orders") == [0, 1]
            and first_source.get("reduced_volume_D") == [184880, 184881]
            and first_source.get("reduced_S1") == [184882, 184883],
            "first-switch B source binding changed")
    nested_partition = parent.get("actual_pressure_window_transport", {}).get("source_partition")
    require(isinstance(nested_partition, Mapping)
            and nested_partition == parent.get("source_partition")
            and nested_partition.get("rows")
            and len(nested_partition["rows"]) == 2,
            "parent source partition is missing")


def _source_component(parent: Mapping[str, object], native: Mapping[str, object]) -> dict[str, object]:
    _validate_source_bindings(parent, native)
    preserved = parent.get("preserved", {})
    require(preserved.get("cumulative_terms") == 212
            and preserved.get("signed_core_lower_terms") == 6
            and preserved.get("reserved_n1_own_swirl_registration_count") == 1
            and preserved.get("all_inherited_errors") is True,
            "parent preservation ledger changed")
    old = copy.deepcopy(parent["actual_pressure_window_transport"]["source_partition"])
    rows = []
    removed = []
    for order, row in enumerate(old["rows"]):
        require(row.get("ordinary_Z_order") == order, "parent source row order changed")
        atoms = row.get("remaining_signed_atoms")
        require(isinstance(atoms, Mapping), "parent remaining atoms missing")
        pair = REMOVED[order]
        require(atoms.get(str(pair[0])) == -1 and atoms.get(str(pair[1])) == 1,
                "expected seam atoms are not present")
        reduced = {str(k): v for k, v in atoms.items()
                   if int(k) not in pair}
        require(len(reduced) == 5, "seam removal did not leave five atoms")
        require(set(reduced) == {str(handle) for handle in
                                 (183935 + order, 184677 + order, 184884 + order,
                                  185280 + order, 185292 + order)}
                and all(value == 1 for value in reduced.values()),
                "remaining signed atom set changed")
        rows.append({"ordinary_Z_order": order,
                     "adopted_reduced_root": row["adopted_reduced_root"],
                     "removed_seam_atoms": {str(pair[0]): -1, str(pair[1]): 1},
                     "remaining_signed_atoms": reduced,
                     "remaining_atoms_opaque": True})
        removed.append({"ordinary_Z_order": order, "frozen_macro_upper_m1": pair[0],
                        "first_switch_lower_s0": pair[1], "signs": [-1, 1],
                        "same_source_seam": True})
    seam = native["actual_macro_switch_RQy_seam"]
    joined = native["actual_micro_macro_correlated_m15"]["actual_joined_pressure_function"]
    correlated = native["actual_micro_macro_correlated_m15"]
    require(correlated.get("old_212_terms_and_errors_unchanged") is True
            and correlated.get("old_six_numeric_terms_and_all_error_caps_unchanged") is True
            and correlated.get("signed_function_values_evaluated") is False
            and correlated.get("signed_six_integrals_evaluated") is False
            and correlated.get("accepted_ancestors_replayed") is False
            and correlated.get("numerical_phase_or_quadrature_evaluated") is False,
            "old signed terms/error ledger changed")
    require(joined.get("absolute_M15_formula") == "M15(m=1)=accepted_H+joined_shifted_tail"
            and joined.get("shifted_H_anchor_digest") ==
            "ec934fbfce5aa21708541d5d65a87a920c3239e599f333370aef1dc650985e17"
            and joined.get("joined_shifted_tail") == [184683, 184684]
            and joined.get("pressure_increment") == [184685, 184686]
            and joined.get("pressure_increment_numerical_enclosure") is False
            and joined.get("original_first_lower_M15_retained") is True
            and joined.get("shifted_H_anchor_retained") is True,
            "accepted H anchor receipt changed")
    source_graph = native["actual_restoration_pressure_source"]["graph_patch"]
    corr_graph = native["actual_restoration_correlated_pressure"]["graph_patch"]
    return {
        "schema": "original_pressure_prefix_contraction_source_v1",
        "passed": True,
        "parent_source_partition": old,
        "removed_same_seam_pairs": removed,
        "remaining_rows": rows,
        "remaining_atoms_per_Z": 5,
        "opaque_atoms": [{"ordinary_Z_order": order,
                          "handles": list(OPAQUE[order]),
                          "reason": "unresolved first-micro lower RQy and outer-micro/macro half-S1"}
                         for order in ORDERS],
        "first_switch_B_source": copy.deepcopy(
            native["actual_macro_first_switch_pressure_source"]["first_switch_source"]["B_source"]),
        "seam_receipt": {
            "proof_rows": copy.deepcopy(seam["proof_rows"]),
            "symbolic_controls": copy.deepcopy(seam["symbolic_controls"]),
            "common_Qy_radius": True, "same_N": True,
            "nonreset_M": True, "distinct_native_J": True,
        },
        "accepted_H_anchor_receipt": {
            "absolute_M15_formula": joined["absolute_M15_formula"],
            "shifted_H_anchor_digest": joined["shifted_H_anchor_digest"],
            "joined_shifted_tail": copy.deepcopy(joined["joined_shifted_tail"]),
            "pressure_increment": copy.deepcopy(joined["pressure_increment"]),
            "pressure_increment_source_defined": True,
            "pressure_increment_numerical_enclosure": False,
            "do_not_add_raw_increment_to_H": True,
            "first_lower_RQy_remains_opaque": True,
            "identity": "absolute_M15_at_Rm/Rm^2 = H/Rm^2 + outer_halfS1/Rm^2 + new_partial_rows",
            "pressure_identity": "pressure_increment_adopted/Rm^2 = firstlowerRQy/Rm^2 + outerhalfS1/Rm^2 + new_partial_rows",
            "source_only": True,
        },
        "graph_patch_receipt": {
            "restoration_base_node": 185704,
            "restoration_base_sha256": source_graph["base_canonical_nodes_sha256"],
            "restoration_full_node": 186086,
            "restoration_full_sha256": source_graph["full_canonical_nodes_sha256"],
            "correlated_full_node": 186100,
            "correlated_full_sha256": corr_graph["full_canonical_nodes_sha256"],
            "segment_starts": copy.deepcopy(source_graph["segment_starts"]),
        },
        "source_contract": {
            "old_certificate_preserved_verbatim": True,
            "old_partition_preserved_verbatim": True,
            "only_two_same_seam_pairs_removed": True,
            "standalone_endpoints_and_error_caps_preserved": True,
            "opaque_atoms_retained_explicitly": True,
            "source_execution": False,
            "accepted_ancestors_replayed": False,
            "cumulative_terms_preserved": 212,
            "signed_core_lower_terms_preserved": 6,
            "reserved_n1_own_swirl_registration_count_preserved": 1,
            "new_reserved_term": False,
        },
    }


def _setup(audit: Mapping[str, object] | None = None) -> dict[str, object]:
    if audit is None:
        values, raws = _authenticated_parent_files()
        parent, native = values[PARENT], values[NATIVE]
    else:
        parent, native = audit["parent"], audit["native"]
        raws = {PARENT: audit["parent_raw"], NATIVE: audit["native_raw"]}
    require(parent.get("source_only") is False
            and parent.get("accepted_ancestors_replayed") is False,
            "parent scope changed")
    _validate_source_bindings(parent, native)
    partition = parent["source_partition"]
    first = native["actual_first_switch_background_I1_numeric"]
    first_fin = native["actual_first_switch_finite_N_I1_numeric"]
    switch_bg = native["actual_switch_post_background_I1_numeric"]
    switch_fin = native["actual_switch_post_finite_N_I1_numeric"]
    rest_bg = native["actual_restoration_background_I1_numeric"]
    for item, label in ((first, "first background"), (first_fin, "first finite-N"),
                        (switch_bg, "switch background"), (switch_fin, "switch finite-N"),
                        (rest_bg, "restoration background")):
        require(item.get("passed") is True, label + " packet did not pass")
    bg_contract = first.get("source_contract", {})
    require(bg_contract.get("normalized_by_Wswitch") is True
            and bg_contract.get("incoming_memory_retained") is True
            and bg_contract.get("signed_values_are_enclosures") is True
            and bg_contract.get("finite_N_added") is False,
            "first background source contract changed")
    geometry_receipt = first.get("geometry_measure", {})
    require(geometry_receipt.get("radius_Z_independent") is True
            and geometry_receipt.get("length_Z_independent") is True
            and first.get("actual_phase_Z_independent") is True
            and first.get("Q_b_proof", {}).get("ordinary_product_binomial") is True,
            "first boundary ordinary-Z geometry premise changed")
    fin_contract = first_fin.get("source_contract", {})
    require(fin_contract.get("physical_delta_U_local_only") is True
            and fin_contract.get("no_incoming_U_added") is True
            and fin_contract.get("incoming_physical_delta_M_from_macro") is True
            and fin_contract.get("physical_Pstar_applied_once") is True
            and fin_contract.get("phase_Z_independent") is True
            and fin_contract.get("nonreset_mean_positive_kernel") is True,
            "first finite-N source contract changed")
    d_full = first_fin["D_full"]
    require(d_full.get("passed") is True
            and d_full.get("computed_from_full_I1_once") is True
            and d_full.get("separate_D_error_not_added") is True
            and d_full.get("normalized_by_Wswitch") is True, "first D source changed")
    s1_rows = {}
    for chart in ("second_switch", "post_switch"):
        item = switch_fin["charts"][chart]
        full = item["S1_full"]
        require(item.get("passed") is True and full.get("passed") is True
                and full.get("computed_from_B_and_D_once") is True
                and full.get("separate_S1_error_not_added") is True,
                "S1 source changed: " + chart)
        s1_rows[chart] = _pair(full["rows"], chart + ".S1_full")
    geometry = first["geometry_measure"]
    post_geometry = switch_bg["charts"]["second_switch"]["geometry"]
    require(_record_equal(geometry["loghbS"], post_geometry["loghbS"]),
            "second-switch loghbS does not match first-switch source")
    switch_contract = switch_bg.get("source_contract", {})
    require(switch_contract.get("normalized_second") == "Wswitch"
            and switch_contract.get("normalized_post") == "Wmacro"
            and switch_contract.get("incoming_memory_retained") is True
            and switch_contract.get("finite_N_phase_not_selected") is True
            and switch_contract.get("source_executed") is False,
            "switch background scale contract changed")
    rest_scale = rest_bg.get("physical_scale", {})
    require(rest_scale.get("Wpost") == "Rm^2"
            and rest_scale.get("radius_materialized") is False,
            "post-restore Rm scale contract changed")
    first_w = _record(geometry["logWswitch"], "first logWswitch")
    hb = _record(geometry["loghbS"], "loghbS")
    macro_w = native["actual_macro_background_I1_numeric"]["geometric_measure"]["Wmacro_record"]
    require(_record_equal(macro_w, {"lower_exact_mpf_tuple": [0, 625, 4, 10],
                                    "upper_exact_mpf_tuple": [0, 625, 4, 10],
                                    "directed": True}), "Wmacro is not exact 10000")
    log_rm = _record(rest_bg["charts"]["post_restore"]["log_R_upper"], "log_R_upper")
    u_rows = _pair(first["native_fields"]["U"], "first U")
    q_rows = _pair(first["Q_b"], "first Q_b")
    proof = first["Q_b_proof"]
    z = _record(proof["Z"][0], "Q_b_proof.Z")
    delta = _record(proof["delta"], "Q_b_proof.delta")
    errors = first_fin["physical_error_log_bounds"]
    du = [_record(row["log_absolute_upper_bound"], "delta_U.log_absolute_upper_bound")
          for row in errors["delta_U"][:2]]
    dq = [_record(row["log_absolute_upper_bound"], "delta_Q.log_absolute_upper_bound")
          for row in errors["delta_Q"][:2]]
    require(len(u_rows) == len(q_rows) == len(du) == len(dq) == 2,
            "first-switch boundary row shape changed")
    r_upper = _record(geometry["R_upper"], "first R_upper")
    b_source = copy.deepcopy(native["actual_macro_first_switch_pressure_source"]["first_switch_source"]["B_source"])
    old_partial = parent["actual_pressure_window_transport"]["partial_window_rows"]
    parent_partial = [_record(row, "parent partial") for row in old_partial]
    require(len(parent_partial) == 2, "parent partial rows changed")
    return {
        "values": {PARENT: parent, NATIVE: native}, "raws": raws,
        "parent": parent, "native": native, "source_component": _source_component(parent, native),
        "parent_partial": parent_partial, "D_first": _pair(d_full["rows"], "first D_full"),
        "S1_second": s1_rows["second_switch"], "S1_post": s1_rows["post_switch"],
        "U_first": u_rows, "Q_first": q_rows, "z": z, "delta": delta,
        "delta_U_logs": du, "delta_Q_logs": dq, "R_upper": r_upper,
        "logW_first": first_w, "loghbS": hb, "logRm": log_rm,
        "B_source_binding": b_source,
    }


def _to_interval(c: object, record: Mapping[str, object], make_mpf: object, label: str):
    canonical = _record(record, label)
    lo = make_mpf(tuple(canonical["lower_exact_mpf_tuple"]))
    hi = make_mpf(tuple(canonical["upper_exact_mpf_tuple"]))
    return c.mpf([lo, hi])


def _serialize_interval(value: object, label: str = "interval") -> dict[str, object]:
    rows = getattr(value, "_mpi_", None)
    require(rows is not None, label + " is not an interval")
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    require(len(rows) == 2, label + " has malformed endpoints")
    lower, upper = _endpoint(tuple(rows[0]), label + ".lower"), _endpoint(tuple(rows[1]), label + ".upper")
    require(_cmp_exact(lower, upper) <= 0, label + " is reversed")
    return {"lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper), "directed": True}


def _write_json(path: Path, value: Mapping[str, object]) -> None:
    require(not path.exists(), "saved output exists; overwrite refused: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def _checkpoint(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    value = json.loads(path.read_bytes())
    require(isinstance(value, dict), "checkpoint is not an object: " + str(path))
    return value


def _checkpoint_meta(payload: Mapping[str, object], seconds: float, audit: Mapping[str, object]) -> dict[str, object]:
    return {
        "checkpoint_schema": "original_pressure_prefix_contraction_checkpoint_v1",
        "payload": payload, "interval_precision_bits": 1600,
        "execution_seconds": seconds, "computation_runner_sha256": sha(ROOT / RUNNER),
        "prefix_kernel_sha256": KERNEL_SHA,
        "boundary_kernel_sha256": audit["boundary"]["sha256"],
        "parent_report": PARENT, "parent_report_sha256": PARENT_SHA,
        "native_report": NATIVE, "native_report_sha256": NATIVE_SHA,
        "accepted_ancestors_replayed": False, "source_execution": False,
    }


def _compute_new(setup: Mapping[str, object], audit: Mapping[str, object]) -> dict[str, object]:
    require(not REPORT.exists(), "new report exists; production/replay refused")
    source_cp = _checkpoint(CHECKPOINTS["source"])
    if source_cp is not None:
        _validate_checkpoint_meta(source_cp, audit)
        require(source_cp.get("payload") == setup["source_component"],
                "existing source checkpoint changed")
    if source_cp is None:
        source_cp = _checkpoint_meta(setup["source_component"], 0.0, audit)
        _write_json(CHECKPOINTS["source"], source_cp)
    boundary_cp = _checkpoint(CHECKPOINTS["boundary"])
    transport_cp = _checkpoint(CHECKPOINTS["transport"])
    join_cp = _checkpoint(CHECKPOINTS["join"])
    for name, value in (("boundary", boundary_cp), ("transport", transport_cp),
                        ("join", join_cp)):
        if value is not None:
            _validate_checkpoint_meta(value, audit)
            payload = value.get("payload")
            require(isinstance(payload, Mapping), "existing " + name + " payload is missing")
            expected_schema = {
                "boundary": "original_first_switch_boundary_enclosure_v1",
                "transport": "original_pressure_prefix_log_transport_v1",
                "join": "original_pressure_prefix_contraction_v1",
            }[name]
            require(payload.get("schema") == expected_schema and payload.get("passed") is True,
                    "existing " + name + " checkpoint changed")
            if name == "boundary":
                for key, expected in (("input_U_rows", setup["U_first"]),
                                      ("input_Q_rows", setup["Q_first"]),
                                      ("input_delta_U_logs", setup["delta_U_logs"]),
                                      ("input_delta_Q_logs", setup["delta_Q_logs"]),
                                      ("input_z", setup["z"]),
                                      ("input_delta", setup["delta"]),
                                      ("input_R_upper", setup["R_upper"]),
                                      ("input_logRm", setup["logRm"]),
                                      ("source_binding", setup["B_source_binding"])):
                    require(payload.get(key) == expected,
                            "existing boundary source binding changed: " + key)
            elif name == "transport":
                require(payload.get("input_quantities") == {
                    "first_switch": "D_first", "second_switch": "S1_second", "post_power": "S1_post"}
                        and payload.get("input_rows") == {
                            "first_switch": setup["D_first"],
                            "second_switch": setup["S1_second"],
                            "post_power": setup["S1_post"]}
                        and payload.get("logRm") == setup["logRm"],
                        "existing transport source binding changed")
    # Numeric imports are intentionally confined to this compute path and
    # happen only after every existing checkpoint has passed its data audit.
    sys.path.insert(0, str(ROOT / "src"))
    from mpmath import make_mpf
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_first_switch_boundary_enclosure as boundary
    from openai_ns_reconstruction import original_pressure_prefix_log_transport as transport
    context = MPIntervalContext()
    context.prec = 1600
    if boundary_cp is None:
        started = time.monotonic()
        U = tuple(_to_interval(context, item, make_mpf, "U") for item in setup["U_first"])
        Q = tuple(_to_interval(context, item, make_mpf, "Q") for item in setup["Q_first"])
        dU = tuple(_to_interval(context, item, make_mpf, "delta_U") for item in setup["delta_U_logs"])
        dQ = tuple(_to_interval(context, item, make_mpf, "delta_Q") for item in setup["delta_Q_logs"])
        result = boundary.bound_boundary(
            context, U, Q, _to_interval(context, setup["z"], make_mpf, "z"),
            _to_interval(context, setup["delta"], make_mpf, "delta"), dU, dQ,
            _to_interval(context, setup["R_upper"], make_mpf, "R_upper"),
            _to_interval(context, setup["logRm"], make_mpf, "logRm"), radius_bits=1024)
        raw = result["rawIVresult"]
        payload = {
            "schema": "original_first_switch_boundary_enclosure_v1", "passed": True,
            "input_quantity": "B_first", "normalization": "own_Rm^2",
            "input_U_rows": copy.deepcopy(setup["U_first"]),
            "input_Q_rows": copy.deepcopy(setup["Q_first"]),
            "input_delta_U_logs": copy.deepcopy(setup["delta_U_logs"]),
            "input_delta_Q_logs": copy.deepcopy(setup["delta_Q_logs"]),
            "input_z": copy.deepcopy(setup["z"]), "input_delta": copy.deepcopy(setup["delta"]),
            "input_R_upper": copy.deepcopy(setup["R_upper"]),
            "input_logRm": copy.deepcopy(setup["logRm"]),
            "source_binding": copy.deepcopy(setup["B_source_binding"]),
            "raw_rows": [_serialize_interval(row, "B_first[Z" + str(i) + "]")
                         for i, row in enumerate(raw)],
            "metadata": copy.deepcopy(result["metadata"]),
            "source_contract": {"first_switch_D_is_not_S1": True,
                                "B_first_already_Rm2_normalized": True,
                                "B_first_not_divided_by_Wswitch_again": True,
                                "positive_symmetric_bound_is_not_value": True,
                                "new_values_are_directed_function_enclosures": True},
        }
        boundary_cp = _checkpoint_meta(payload, time.monotonic() - started, audit)
        _write_json(CHECKPOINTS["boundary"], boundary_cp)
    if transport_cp is None:
        started = time.monotonic()
        rows = {
            "first_switch": tuple(_to_interval(context, item, make_mpf, "D_first")
                                   for item in setup["D_first"]),
            "second_switch": tuple(_to_interval(context, item, make_mpf, "S1_second")
                                    for item in setup["S1_second"]),
            "post_power": tuple(_to_interval(context, item, make_mpf, "S1_post")
                                 for item in setup["S1_post"]),
        }
        logW = {
            "first_switch": _to_interval(context, setup["logW_first"], make_mpf, "logW_first"),
            "second_switch": context.ln(context.mpf(10000))
            + _to_interval(context, setup["loghbS"], make_mpf, "loghbS"),
            "post_power": context.ln(context.mpf(10000)),
        }
        result = transport.transport(context, rows, logW,
                                     _to_interval(context, setup["logRm"], make_mpf, "logRm"),
                                     radius_bits=1024, geometry_Z_independent=True)
        payload = {
            "schema": "original_pressure_prefix_log_transport_v1", "passed": True,
            "input_quantities": {"first_switch": "D_first", "second_switch": "S1_second",
                                  "post_power": "S1_post"},
            "input_rows": {chart: [_serialize_interval(row, chart) for row in rows[chart]]
                           for chart in CHARTS},
            "logW_by_chart": {chart: _serialize_interval(value, chart + ".logW")
                              for chart, value in logW.items()},
            "logW_provenance": {
                "first_switch": "saved actual_first_switch_background_I1_numeric.geometry_measure.logWswitch",
                "second_switch": "ln(10000)+saved actual_switch_post_background_I1_numeric.charts.second_switch.geometry.loghbS",
                "post_power": "ln(10000) from saved Wmacro=10000; no R_upper^2",
                "logRm": "saved actual_restoration_background_I1_numeric.charts.post_restore.log_R_upper",
            },
            "logRm": copy.deepcopy(setup["logRm"]),
            "raw_rows": [_serialize_interval(row, "transport[Z" + str(i) + "]")
                         for i, row in enumerate(result["rawIVresult"])],
            "metadata": copy.deepcopy(result["metadata"]),
            "source_contract": {"D_first_used_once": True, "S1_second_used_once": True,
                                "S1_post_used_once": True, "geometry_Z_independent": True,
                                "source_half_factor_applied": False,
                                "unmaterializable_weights_exponentiated": False},
        }
        transport_cp = _checkpoint_meta(payload, time.monotonic() - started, audit)
        _write_json(CHECKPOINTS["transport"], transport_cp)
    if join_cp is None:
        started = time.monotonic()
        b_rows = boundary_cp["payload"]["raw_rows"]
        t_rows = transport_cp["payload"]["raw_rows"]
        parent_rows = setup["parent_partial"]
        b_iv = tuple(_to_interval(context, row, make_mpf, "B row") for row in b_rows)
        t_iv = tuple(_to_interval(context, row, make_mpf, "transport row") for row in t_rows)
        p_iv = tuple(_to_interval(context, row, make_mpf, "parent row") for row in parent_rows)
        half = context.mpf(1) / 2
        joined = tuple(p_iv[order] + half * (b_iv[order] + t_iv[order]) for order in ORDERS)
        payload = {
            "schema": "original_pressure_prefix_contraction_v1", "passed": True,
            "normalization": "own_Rm^2", "ordinary_Z_orders": [0, 1],
            "parent_partial_window_rows": copy.deepcopy(parent_rows),
            "first_B_rows": copy.deepcopy(b_rows),
            "transport_rows": copy.deepcopy(t_rows),
            "half_factor": {"numerator": 1, "denominator": 2,
                            "applied_to_new_B_plus_transport": True},
            "new_partial_rows": [_serialize_interval(row, "new_partial[Z" + str(i) + "]")
                                 for i, row in enumerate(joined)],
            "opaque_atoms": copy.deepcopy(setup["source_component"]["opaque_atoms"]),
            "functional_identity": {
                "formula": "adopted_root/Rm^2 = twoOpaque/Rm^2 + new_partial_rows",
                "twoOpaque_per_Z": [list(pair) for pair in OPAQUE],
                "geometry": "frozen", "full_pressure_numeric": False,
                "global_closure": False, "independent_m15_claim": False,
            },
            "source_contract": {"parent_partial_rows_preserved": True,
                                "new_B_and_transport_halved_once": True,
                                "parent_half_S1_not_reapplied": True,
                                "opaque_terms_not_numeric": True,
                                "new_values_are_directed_function_enclosures": True},
        }
        join_cp = _checkpoint_meta(payload, time.monotonic() - started, audit)
        _write_json(CHECKPOINTS["join"], join_cp)
    return {"source": source_cp, "boundary": boundary_cp,
            "transport": transport_cp, "join": join_cp}


def _validate_checkpoint_meta(value: Mapping[str, object], audit: Mapping[str, object]) -> None:
    require(value.get("checkpoint_schema") ==
            "original_pressure_prefix_contraction_checkpoint_v1",
            "checkpoint schema changed")
    require(value.get("prefix_kernel_sha256") == KERNEL_SHA
            and value.get("boundary_kernel_sha256") == audit["boundary"]["sha256"],
            "checkpoint kernel pin changed")
    require(value.get("parent_report") == PARENT and value.get("parent_report_sha256") == PARENT_SHA
            and value.get("native_report") == NATIVE and value.get("native_report_sha256") == NATIVE_SHA,
            "checkpoint source pin changed")
    require(value.get("computation_runner_sha256") == sha(ROOT / RUNNER)
            and value.get("interval_precision_bits") == 1600
            and isinstance(value.get("execution_seconds"), (int, float))
            and value.get("execution_seconds") >= 0,
            "checkpoint execution metadata changed")
    require(value.get("accepted_ancestors_replayed") is False
            and value.get("source_execution") is False, "checkpoint audit changed")


def _validate_components(checkpoints: Mapping[str, Mapping[str, object]],
                         setup: Mapping[str, object], audit: Mapping[str, object]) -> int:
    for name, value in checkpoints.items():
        _validate_checkpoint_meta(value, audit)
        require(isinstance(value.get("payload"), Mapping), "checkpoint payload missing: " + name)
    source = checkpoints["source"]["payload"]
    require(source == setup["source_component"], "source component changed")
    boundary = checkpoints["boundary"]["payload"]
    require(boundary.get("schema") == "original_first_switch_boundary_enclosure_v1"
            and boundary.get("passed") is True and boundary.get("input_quantity") == "B_first",
            "boundary component changed")
    for key, expected in (("input_U_rows", setup["U_first"]),
                          ("input_Q_rows", setup["Q_first"]),
                          ("input_delta_U_logs", setup["delta_U_logs"]),
                          ("input_delta_Q_logs", setup["delta_Q_logs"]),
                          ("input_z", setup["z"]), ("input_delta", setup["delta"]),
                          ("input_R_upper", setup["R_upper"]), ("input_logRm", setup["logRm"])):
        require(boundary.get(key) == expected, "boundary input changed: " + key)
    require(boundary.get("source_binding") == setup["B_source_binding"],
            "boundary source binding changed")
    for row in boundary.get("raw_rows", []): _record(row, "boundary raw row")
    require(len(boundary.get("raw_rows", [])) == 2
            and boundary.get("metadata", {}).get("schema") ==
            "original_first_switch_boundary_enclosure_v1"
            and boundary.get("source_contract", {}).get("B_first_already_Rm2_normalized") is True
            and boundary.get("source_contract", {}).get("B_first_not_divided_by_Wswitch_again") is True,
            "boundary output metadata changed")
    transport = checkpoints["transport"]["payload"]
    require(transport.get("schema") == "original_pressure_prefix_log_transport_v1"
            and transport.get("passed") is True
            and transport.get("input_quantities") == {
                "first_switch": "D_first", "second_switch": "S1_second", "post_power": "S1_post"},
            "transport component changed")
    expected_rows = {"first_switch": setup["D_first"],
                     "second_switch": setup["S1_second"], "post_power": setup["S1_post"]}
    for chart in CHARTS:
        actual = transport.get("input_rows", {}).get(chart)
        require(isinstance(actual, list) and len(actual) == 2
                and actual == expected_rows[chart], "transport input rows changed: " + chart)
    for row in transport.get("raw_rows", []): _record(row, "transport raw row")
    require(len(transport.get("raw_rows", [])) == 2,
            "transport raw rows changed")
    metadata = transport.get("metadata", {})
    require(metadata.get("schema") == "original_pressure_prefix_log_transport_v1"
            and metadata.get("geometry_Z_independent") is True
            and metadata.get("radius_bits") == 1024
            and metadata.get("source_contract", {}).get("unmaterializable_weight_exponentiated") is False,
            "transport metadata changed")
    require(isinstance(transport.get("logW_provenance"), Mapping),
            "transport logW provenance missing")
    join = checkpoints["join"]["payload"]
    require(join.get("schema") == "original_pressure_prefix_contraction_v1"
            and join.get("passed") is True
            and join.get("parent_partial_window_rows") == setup["parent_partial"]
            and join.get("first_B_rows") == boundary.get("raw_rows")
            and join.get("transport_rows") == transport.get("raw_rows")
            and len(join.get("new_partial_rows", [])) == 2,
            "join component changed")
    for row in join["new_partial_rows"]: _record(row, "new partial row")
    return sum(_records_count(checkpoints[name]["payload"]) for name in checkpoints)


def _merge_pins(*maps: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for mapping in maps:
        for name, value in mapping.items():
            require(name not in result or result[name] == value,
                    "conflicting source pin: " + name)
            result[name] = value
    return result


def _expected_inputs(setup: Mapping[str, object]) -> dict[str, str]:
    parent, native = setup["parent"], setup["native"]
    result = _merge_pins(parent["input_hashes"], native["input_hashes"])
    for name, value in ((PARENT, PARENT_SHA), (NATIVE, NATIVE_SHA),
                        (FIXTURE, FIXTURE_SHA)):
        require(name not in result or result[name] == value, "input pin conflict: " + name)
        result[name] = value
    for name, (_, value) in BOUNDARY_EVIDENCE.items():
        require(name not in result or result[name] == value, "input pin conflict: " + name)
        result[name] = value
    return result


def _expected_implementations(setup: Mapping[str, object], audit: Mapping[str, object]) -> dict[str, str]:
    result = _merge_pins(setup["parent"]["implementation_hashes"],
                         setup["native"]["implementation_hashes"])
    for name, value in ((KERNEL, KERNEL_SHA), (BOUNDARY, audit["boundary"]["sha256"]),
                        (RUNNER, sha(ROOT / RUNNER))):
        require(name not in result or result[name] == value,
                "implementation pin conflict: " + name)
        result[name] = value
    return result


def _package(setup: Mapping[str, object], audit: Mapping[str, object]) -> dict[str, object]:
    require(not REPORT.exists(), "report overwrite refused")
    checkpoints = {name: _checkpoint(path) for name, path in CHECKPOINTS.items()}
    require(all(value is not None for value in checkpoints.values()), "new checkpoint missing")
    count = _validate_components(checkpoints, setup, audit)
    payloads = {name: checkpoints[name]["payload"] for name in CHECKPOINTS}
    parent = setup["parent"]
    report = copy.deepcopy(parent)
    checkpoint_hashes = {name: {"bytes": len(CHECKPOINTS[name].read_bytes()),
                                "sha256": sha(CHECKPOINTS[name])}
                         for name in CHECKPOINTS}
    elapsed = sum(float(checkpoints[name]["execution_seconds"]) for name in CHECKPOINTS)
    report.update({
        "schema": "original_pressure_prefix_contraction_v1",
        "stage": "current_original_pressure_prefix_contraction",
        "parent_report": PARENT, "parent_report_sha256": PARENT_SHA,
        "immutable_parent_packet_sha256": PARENT_SHA,
        "input_hashes": _expected_inputs(setup),
        "implementation_hashes": _expected_implementations(setup, audit),
        "actual_pressure_prefix_contraction": payloads["join"],
        "actual_first_switch_boundary_enclosure": payloads["boundary"],
        "actual_pressure_prefix_log_transport": payloads["transport"],
        "prefix_source_component": payloads["source"],
        "new_partial_rows": copy.deepcopy(payloads["join"]["new_partial_rows"]),
        "prefix_functional_identity": copy.deepcopy(payloads["join"]["functional_identity"]),
        "authenticated_source_reports": {
            **copy.deepcopy(parent.get("authenticated_source_reports", {})),
            PARENT: {"bytes": PARENT_BYTES, "sha256": PARENT_SHA},
            NATIVE: {"bytes": NATIVE_BYTES, "sha256": NATIVE_SHA},
            FIXTURE: {"bytes": FIXTURE_BYTES, "sha256": FIXTURE_SHA},
            **copy.deepcopy(audit["boundary_evidence"]["files"]),
        },
        "execution_audit": {
            "accepted_ancestors_replayed": False, "source_execution": False,
            "new_saved_data_production_attempts": 1,
            "checkpoint_hashes": checkpoint_hashes,
            "computation_runner_sha256": sha(ROOT / RUNNER),
            "prefix_kernel_sha256": KERNEL_SHA,
            "boundary_kernel_sha256": audit["boundary"]["sha256"],
            "canonical_exact_tuple_records": count,
        },
        "source_only": False, "accepted_ancestors_replayed": False,
        "execution_seconds": elapsed, "interval_precision_bits": 1600,
    })
    fixtures = copy.deepcopy(parent.get("fresh_fixtures", {}))
    fixtures[FIXTURE] = {
        "bytes": FIXTURE_BYTES, "sha256": FIXTURE_SHA,
        "module_sha256": KERNEL_SHA, "source_executed": False,
    }
    report["fresh_fixtures"] = fixtures
    scope = dict(parent.get("scope", {}))
    scope.update({"pressure_prefix_contraction_partial_group_computed": True,
                  "partial_pressure_root_group_numeric": True,
                  "new_values_are_directed_function_enclosures": True,
                  "point_centers_selected": False})
    for key in FALSE_SCOPE: scope[key] = False
    report["scope"] = scope
    report["limits"] = list(parent.get("limits", [])) + [
        "Prefix contraction covers D_first, S1_second, and S1_post only; two opaque atoms remain per ordinary Z row.",
        "The first-switch B boundary and all transported rows are directed function enclosures, not selected centers.",
        "The parent source certificate, standalone endpoints, and finite-N error caps are preserved verbatim.",
        "No full pressure, global closure, or independent m15 claim is made.",
    ]
    return report


def package() -> None:
    audit = preflight()
    setup = _setup(audit)
    report = _package(setup, audit)
    _write_json(REPORT, report)
    raw = REPORT.read_bytes()
    print(json.dumps({"report": str(REPORT), "bytes": len(raw), "sha256": digest(raw)}))


def inspect() -> None:
    audit = preflight()
    setup = _setup(audit)
    require(REPORT.exists(), "prefix contraction report is missing")
    checkpoints = {name: _checkpoint(path) for name, path in CHECKPOINTS.items()}
    require(all(value is not None for value in checkpoints.values()), "new checkpoint missing")
    count = _validate_components(checkpoints, setup, audit)
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    require(isinstance(report, dict), "report is not an object")
    require(report.get("input_hashes") == _expected_inputs(setup)
            and report.get("implementation_hashes") == _expected_implementations(setup, audit),
            "report source pins changed")
    require(report.get("actual_pressure_window_transport") == setup["parent"].get("actual_pressure_window_transport")
            and report.get("source_partition") == setup["parent"].get("source_partition")
            and report.get("preserved") == setup["parent"].get("preserved"),
            "immutable parent packet changed")
    payloads = {name: checkpoints[name]["payload"] for name in CHECKPOINTS}
    require(report.get("actual_pressure_prefix_contraction") == payloads["join"]
            and report.get("actual_first_switch_boundary_enclosure") == payloads["boundary"]
            and report.get("actual_pressure_prefix_log_transport") == payloads["transport"]
            and report.get("prefix_source_component") == payloads["source"]
            and report.get("new_partial_rows") == payloads["join"]["new_partial_rows"],
            "report components do not match checkpoints")
    execution = report.get("execution_audit")
    expected_hashes = {name: {"bytes": len(CHECKPOINTS[name].read_bytes()),
                              "sha256": sha(CHECKPOINTS[name])}
                       for name in CHECKPOINTS}
    require(isinstance(execution, Mapping)
            and execution.get("checkpoint_hashes") == expected_hashes
            and execution.get("accepted_ancestors_replayed") is False
            and execution.get("source_execution") is False,
            "report execution audit changed")
    require(report.get("source_only") is False
            and report.get("accepted_ancestors_replayed") is False
            and report.get("interval_precision_bits") == 1600,
            "report scope metadata changed")
    scope = report.get("scope")
    require(isinstance(scope, Mapping)
            and scope.get("pressure_prefix_contraction_partial_group_computed") is True
            and all(scope.get(key) is False for key in FALSE_SCOPE),
            "report scope overclaimed")
    print(json.dumps({"passed": True, "report": REPORT_NAME,
                      "sha256": digest(raw), "bytes": len(raw),
                      "canonical_exact_tuple_records": count,
                      "input_pins": len(report["input_hashes"]),
                      "implementation_pins": len(report["implementation_hashes"]),
                      "inspection_math_providers_imported": False,
                      "accepted_ancestors_replayed": False, "source_execution": False}))


def produce() -> None:
    audit = preflight()
    setup = _setup(audit)
    checkpoints = _compute_new(setup, audit)
    print(json.dumps({"checkpoints": {name: {"sha256": sha(path)}
                                       for name, path in CHECKPOINTS.items()},
                      "kernel_sha256": KERNEL_SHA,
                      "boundary_kernel_sha256": audit["boundary"]["sha256"],
                      "source_execution": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--compute-new", action="store_true")
    actions.add_argument("--write-report", action="store_true")
    actions.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    if args.compute_new:
        produce()
    elif args.write_report:
        package()
    else:
        inspect()
