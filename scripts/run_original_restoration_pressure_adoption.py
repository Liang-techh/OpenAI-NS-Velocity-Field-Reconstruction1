"""Package the source-only restoration pressure adoption receipt.

The producer consumes the authenticated saved endpoint report and the reviewed
source adopter.  Inspection and packaging read only JSON, hashes, and ASTs;
they never import the adopter or any numerical provider.
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
RUNNER = "scripts/run_original_restoration_pressure_adoption.py"
PARENT = "lei_ren_part1_paper_compliant_current_original_post_restore_logu_C1.json"
PARENT_SHA = "1d778c10c272acafee88bd42178ad4ae8198d231de1f8b12d185a09da75f8072"
PARENT_BYTES = 657774
ENDPOINT_REL = "experiments/root_st073/" + PARENT
NAME = "lei_ren_part1_paper_compliant_current_original_restoration_pressure_adoption.json"
REPORT = BASE / NAME
CHECKPOINT = ROOT / "tmp/original_restoration_pressure_adoption/new_source.json"
ADOPTER = "src/openai_ns_reconstruction/original_restoration_pressure_adoption.py"
REVIEWED_ADOPTER_SHA = "275e76bb4f3964078ee83e6c01913daad8ce53f8f11696d7ac6e27e380d3ad4e"
GRAPH_COUNTS = {"prefix": 185704, "pressure_source": 186086, "correlated_pressure_source": 186100}
GRAPH_SHAS = {
    "prefix": "2d6efeb13efe91dd95d4f24fa03eb4626c7685d04e5f29f2bdc6532e101994e3",
    "pressure_source": "dfbba64e1f106fa42c233af28a096f8dfad1e35580163c54cad9eb523d21271e",
    "correlated_pressure_source": "83b8594b5b775520cd38c9da1e6633e800ebf170dfb48b6b12eb34f903514fc0",
}

FALSE_SCOPE = (
    "signed_pressure_join_completed",
    "independent_reserved_full_pulse_m15_matched",
    "global_five_moment_closure",
    "exact_heat_exterior",
    "admissible_stress_and_flat_remainder",
    "true_n_dependent_recursion",
    "realizable_oscillatory_correction",
    "independent_Cartesian_NS_validation",
    "legacy_packets_recertified",
    "actual_phase_or_inverse_values_installed",
    "full_goal_complete",
    "complete_pressure_numeric",
    "full_pressure_numeric",
    "global_closure",
    "terminal_functional_match",
)


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return digest(path.read_bytes())


def input_path(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute(), "Absolute input path")
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py", "Input traversal")
        return (BASE / path).resolve()
    return ROOT / path if path.parts[0] == "experiments" else BASE / path.name


def parent() -> dict:
    raw = (BASE / PARENT).read_bytes()
    require(len(raw) == PARENT_BYTES and digest(raw) == PARENT_SHA, "Parent changed before decode")
    value = json.loads(raw)
    require(isinstance(value, dict), "Parent report is not an object")
    for name, expected in value["input_hashes"].items():
        require(sha(input_path(name)) == expected, "Input changed: " + name)
    for name, expected in value["implementation_hashes"].items():
        require(sha(ROOT / name) == expected, "Implementation changed: " + name)
    return value


def preflight() -> None:
    for name in (ADOPTER, RUNNER):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)
    require(REVIEWED_ADOPTER_SHA != "REVIEW_PENDING", "Reviewed adopter SHA is pending root review")
    require(sha(ROOT / ADOPTER) == REVIEWED_ADOPTER_SHA, "Reviewed adopter source changed")


def write_json(path: Path, value: Mapping[str, object]) -> None:
    require(not path.exists(), "Saved source adoption exists; replay/overwrite refused: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def _mpf(value: object) -> tuple[int, int, int, int]:
    require(isinstance(value, (list, tuple)) and len(value) == 4, "Malformed exact MPF tuple")
    require(all(type(item) is int for item in value), "Noninteger exact MPF tuple")
    sign, mantissa, exponent, bitcount = value
    require(sign in (0, 1) and mantissa >= 0 and bitcount == mantissa.bit_length(), "Noncanonical exact MPF tuple")
    require((mantissa == 0 and tuple(value) == (0, 0, 0, 0)) or (mantissa > 0 and mantissa & 1), "Noncanonical MPF")
    return tuple(value)


def _cmp_abs(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    _, lm, le, _ = left
    _, rm, re, _ = right
    if lm == 0 or rm == 0:
        return (lm > 0) - (rm > 0)
    ltop, rtop = lm.bit_length() + le, rm.bit_length() + re
    if ltop != rtop:
        return (ltop > rtop) - (ltop < rtop)
    if le >= re:
        lv, rv = lm << (le - re), rm
    else:
        lv, rv = lm, rm << (re - le)
    return (lv > rv) - (lv < rv)


def _cmp_mpf(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    if left[1] == 0:
        if right[1] == 0:
            return 0
        return 1 if right[0] else -1
    if right[1] == 0:
        return -1 if left[0] else 1
    if left[0] != right[0]:
        return -1 if left[0] else 1
    result = _cmp_abs(left, right)
    return -result if left[0] else result


def exact_records(tree: object) -> int:
    if isinstance(tree, Mapping):
        count = 0
        if "lower_exact_mpf_tuple" in tree or "upper_exact_mpf_tuple" in tree:
            require("lower_exact_mpf_tuple" in tree and "upper_exact_mpf_tuple" in tree, "Incomplete exact enclosure")
            lower, upper = _mpf(tree["lower_exact_mpf_tuple"]), _mpf(tree["upper_exact_mpf_tuple"])
            require(_cmp_mpf(lower, upper) <= 0, "Reversed directed exact enclosure")
            count = 1
        return count + sum(exact_records(value) for value in tree.values() if isinstance(value, (Mapping, list)))
    if isinstance(tree, list):
        return sum(exact_records(value) for value in tree if isinstance(value, (Mapping, list)))
    return 0


def validate_payload(payload: object) -> dict:
    require(isinstance(payload, dict), "Adopter payload is not an object")
    require(payload.get("schema") == "original_restoration_pressure_adoption_v1", "Adopter schema changed")
    require(payload.get("passed") is True and payload.get("source_adoption") is True, "Source adoption did not pass")
    require(payload.get("namespace_identity_supported") is True, "Adopter namespace identity is not supported")
    require(payload.get("source_execution") is False, "Adopter executed a source")
    require(payload.get("accepted_ancestors_replayed") is False, "Adopter replayed an ancestor")
    for key in ("full_joined_pressure_numeric", "global_five_moment_closure", "complete_pressure_numeric", "full_goal_complete"):
        require(payload.get(key) is False, "Adopter overclaimed: " + key)
    raw = payload["raw_uncontracted_pressure"]
    require(raw.get("roots") == [186084, 186085] and raw.get("joined_contracted") is None, "Blocked raw roots changed")
    require(raw.get("status") == "BLOCKED_UNPROVEN_RESTORATION_POSTRESTORE_SOURCE_CONGRUENCE", "Blocked raw status changed")
    adopted = payload["adopted_pressure_join"]
    require(adopted.get("roots") == [186096, 186099], "Adopted contracted roots changed")
    require(adopted.get("ordinary_Z_orders") == [0, 1], "Adopted Z orders changed")
    require(adopted.get("external_Rm_retained") == [183696, 183697], "External Rm handles changed")
    namespace = payload["namespace"]
    require(namespace.get("raw_uncontracted_roots") == [186084, 186085], "Namespace raw roots changed")
    require(namespace.get("adopted_contracted_roots") == [186096, 186099], "Namespace adopted roots changed")
    require(namespace.get("external_Rm_handles") == [183696, 183697], "Namespace Rm handles changed")
    for graph_name, count in GRAPH_COUNTS.items():
        graph = namespace.get(graph_name)
        require(graph == {"node_count": count, "sha256": GRAPH_SHAS[graph_name]}, "Source graph identity changed: " + graph_name)
    require(set(payload.get("seam_identity_digests", {})) == {"reference_to_restoration", "restoration_to_postrestore"}, "Source seam digest targets changed")
    boundary_graph = namespace.get("endpoint_boundary", {}).get("source_graph")
    require(boundary_graph == {"node_count": GRAPH_COUNTS["pressure_source"], "sha256": GRAPH_SHAS["pressure_source"]}, "Endpoint source graph identity changed")
    endpoint_binding = payload.get("source_bindings", {}).get("endpoint_report", {})
    require(endpoint_binding == {"path": ENDPOINT_REL, "bytes": PARENT_BYTES, "sha256": PARENT_SHA}, "Endpoint source binding changed")
    normalized = payload["external_Rm_normalized_C1"]
    require(normalized.get("status") == "BOUND", "Normalized external Rm status changed")
    require(normalized.get("source_handles") == [183696, 183697] and normalized.get("coefficient") == -1, "Normalized Rm source changed")
    require(normalized.get("physical_scale") == "Wpost=Rm^2" and normalized.get("normalization") == "(Rm*Qy)/Wpost=Qy/Rm", "Rm normalization changed")
    enclosure = normalized.get("enclosure", {})
    require(enclosure.get("directed") is True and enclosure.get("ordinary_Z_orders") == [0, 1], "C1 enclosure metadata changed")
    require(enclosure.get("enclosure_not_actual_zero_or_value") is True, "C1 enclosure was selected as an actual value")
    require(payload.get("api_path") == ADOPTER and payload.get("api_sha256") == sha(ROOT / ADOPTER), "Adopter API pin changed")
    require(exact_records(payload) > 0, "Adopter payload has no exact tuple metadata")
    return payload


def produce() -> None:
    require(not CHECKPOINT.exists() and not REPORT.exists(), "Saved source adoption exists; replay refused")
    preflight()
    parent()
    endpoint_raw = (BASE / PARENT).read_bytes()
    sys.path.insert(0, str(ROOT / "src"))
    from openai_ns_reconstruction import original_restoration_pressure_adoption as provider
    started = time.monotonic()
    payload = provider.adopt(ROOT, endpoint_raw)
    write_json(CHECKPOINT, {
        "payload": payload,
        "execution_seconds": time.monotonic() - started,
        "adopter_module_sha256": sha(ROOT / ADOPTER),
        "computation_runner_sha256": sha(ROOT / RUNNER),
        "parent_report_sha256": PARENT_SHA,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
    })
    validate_payload(payload)
    print(json.dumps({"checkpoint": str(CHECKPOINT), "sha256": sha(CHECKPOINT)}))


def package() -> None:
    require(not REPORT.exists(), "Report overwrite refused")
    preflight()
    old = parent()
    calc = json.loads(CHECKPOINT.read_bytes())
    require(calc.get("parent_report_sha256") == PARENT_SHA and calc.get("source_executed") is False, "Checkpoint/source mismatch")
    require(calc.get("accepted_ancestors_replayed") is False, "Checkpoint replayed an ancestor")
    payload = validate_payload(calc.get("payload"))
    require(payload["external_Rm_normalized_C1"]["enclosure"] == old["actual_post_restore_logu_C1"]["honest_normalized_C1_enclosure"], "Inherited C1 enclosure changed")
    require(calc.get("adopter_module_sha256") == REVIEWED_ADOPTER_SHA, "Checkpoint adopter pin mismatch")
    require(calc.get("computation_runner_sha256") == sha(ROOT / RUNNER), "Checkpoint runner pin mismatch")
    inputs = dict(old["input_hashes"])
    inputs[PARENT] = PARENT_SHA
    implementations = dict(old["implementation_hashes"])
    implementations[ADOPTER] = REVIEWED_ADOPTER_SHA
    implementations[RUNNER] = sha(ROOT / RUNNER)
    report = copy.deepcopy(old)
    report.update({
        "schema": "original_restoration_pressure_adoption_v1",
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA,
        "input_hashes": inputs,
        "implementation_hashes": implementations,
        "actual_restoration_pressure_adoption": payload,
        "source_only": True,
        "accepted_ancestors_replayed": False,
        "restoration_correlated_pressure_source_adopted": payload["source_adoption"],
        "execution_audit": {
            "accepted_ancestors_replayed": False,
            "new_saved_data_production_attempts": 1,
            "checkpoint_sha256": sha(CHECKPOINT),
            "computation_runner_sha256": calc["computation_runner_sha256"],
        },
    })
    scope = dict(old.get("scope", {}))
    scope["restoration_correlated_pressure_source_adopted"] = payload["source_adoption"]
    for key in FALSE_SCOPE:
        scope[key] = False
    report["scope"] = scope
    report["limits"] = list(old.get("limits", [])) + [
        "Source-linked restoration pressure adoption only; no accepted producer or numerical pressure replay.",
        "Contracted roots 186096 and 186099 are adopted while blocked raw roots 186084 and 186085 remain attached.",
        "The external -Rm*Qy boundary is retained and normalized by Wpost=Rm^2; the C1 enclosure is not an actual selected zero.",
    ]
    write_json(REPORT, report)
    print(json.dumps({"report": str(REPORT), "bytes": REPORT.stat().st_size, "sha256": sha(REPORT)}))


def inspect() -> None:
    preflight()
    old = parent()
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    expected_inputs = dict(old["input_hashes"])
    expected_inputs[PARENT] = PARENT_SHA
    expected_implementations = dict(old["implementation_hashes"])
    expected_implementations[ADOPTER] = REVIEWED_ADOPTER_SHA
    expected_implementations[RUNNER] = sha(ROOT / RUNNER)
    require(report.get("input_hashes") == expected_inputs and report.get("implementation_hashes") == expected_implementations, "Source pins changed")
    require(report.get("actual_post_restore_logu_C1") == old.get("actual_post_restore_logu_C1"), "Inherited endpoint packet changed")
    require(report.get("preserved") == old.get("preserved"), "Inherited preservation ledger changed")
    payload = validate_payload(report.get("actual_restoration_pressure_adoption"))
    require(payload["external_Rm_normalized_C1"]["enclosure"] == old["actual_post_restore_logu_C1"]["honest_normalized_C1_enclosure"], "Inherited C1 enclosure changed")
    require(report.get("source_only") is True and report.get("accepted_ancestors_replayed") is False, "Source-only receipt changed")
    require(report.get("restoration_correlated_pressure_source_adopted") is payload["source_adoption"], "Adoption outcome flag changed")
    scope = report.get("scope", {})
    require(scope.get("restoration_correlated_pressure_source_adopted") is payload["source_adoption"], "Scope adoption outcome changed")
    require(all(scope.get(key) is False for key in FALSE_SCOPE), "Downstream scope overclaimed")
    require(report.get("execution_audit", {}).get("accepted_ancestors_replayed") is False, "Execution audit replayed an ancestor")
    require(report.get("execution_audit", {}).get("computation_runner_sha256") == sha(ROOT / RUNNER), "Runner audit pin changed")
    print(json.dumps({
        "passed": True,
        "report": NAME,
        "sha256": digest(raw),
        "bytes": len(raw),
        "canonical_exact_tuple_records": exact_records(payload),
        "input_pins": len(expected_inputs),
        "implementation_pins": len(expected_implementations),
        "restoration_correlated_pressure_source_adopted": payload["source_adoption"],
        "inspection_math_providers_imported": False,
        "accepted_ancestors_replayed": False,
    }))


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
