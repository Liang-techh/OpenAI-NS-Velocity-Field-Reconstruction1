"""Save or inspect new signed actual core derivatives; never replay ancestors."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments/root_st073"
sys.path.insert(0, str(ROOT / "src"))
from openai_ns_reconstruction.original_core_signed_direction import build, PARENT, PARENT_SHA256, HELPER, require, sha, digest

NAME = "lei_ren_part1_paper_compliant_current_original_core_signed_direction.json"
CODE = ("src/openai_ns_reconstruction/original_core_signed_direction.py",
        "scripts/run_original_core_signed_direction.py", HELPER)


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted micro mean parent changed")
    r = json.loads((ROOT / PARENT).read_bytes())
    require(r["checks"]["passed"] is True, "Checked parent required")
    return r


def pins(r):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, h in r[key].items():
            require(sha(base / name) == h, "Changed source pin: " + name)
            count += 1
    return count


def packet_check(r, p):
    require(r["parent_report"] == Path(PARENT).name and r["parent_report_sha256"] == PARENT_SHA256,
            "Different accepted predecessor")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes", "implementation_hashes",
               "contract", "checks", "execution_seconds"}
    for k, v in p.items():
        if k not in mutable:
            require(r[k] == v, "Inherited signed packet or error changed: " + k)
    a = r["actual_core_signed_direction"]
    require(a["passed"] and a["source_family"] == p["source_family"] and
            a["signed_actual_Phi_rho_negative"] and a["actual_core_Dbar_positive"] and
            a["ordinary_Dbar_Z0_Z1_Z2"] and a["actual_core_only"] and a["numeric_controls"]["passed"] and
            a["precision_bits"] == 600 and a["model_terms"] == 64, "New signed actual core proof missing")
    for case in a["cases"].values():
        require(case["passed"] and len(case["mixed_derivatives"]) == 10 and
                len(case["Dbar_ordinary_Z"]) == 3 and case["ordinary_Z_derivatives"] and
                case["actual_fixed_point_error_added_once"] and not case["derivative_of_error_cap_used"],
                "Ordinary total-order3 derivative rows with their own errors required")
        for row in case["mixed_derivatives"]:
            require(row["scaled_radial_order"] + row["ordinary_Z_order"] <= 3 and
                    not row["error_center_selected_as_value"] and not row["factorial_normalization"],
                    "Wrong derivative convention or cap-value substitution")
    require(a["micro_core_argument_domain"]["passed"] and
            not a["micro_core_argument_domain"]["hb_materialized"] and
            a["parameters"]["delta_from_actual_recipe"] and
            a["parameters"]["shared_tube_delta_is_only_universal_cap"], "Actual parameter and domain recipe required")
    for key in ("comparison_second_micro_Dbar_positive", "actual_C5_mixed_root_provider_complete",
                "mean_Z3_numeric", "global_five_moment_closure", "accepted_ancestors_replayed", "full_goal_complete"):
        require(a[key] is False, "Unproved completion: " + key)
    require(all(v is False for v in r["limits"].values()), "Physical completion gate changed")
    return dict(passed=True, all_inherited_packets_and_errors_unchanged=True,
                actual_core_signed_radial_derivative=True, actual_core_Dbar_Z0_Z1_Z2=True,
                selected_and_rectangle_enclosures=True, ten_ordinary_mixed_rows_each=True,
                second_micro_comparison_positivity_claimed=False, accepted_ancestors_replayed=False)


def run():
    started = time.perf_counter(); p = parent(); pins(p); a = build(ROOT, p)
    require(p == parent(), "New enclosure mutated its accepted parent")
    inputs, impl = dict(p["input_hashes"]), dict(p["implementation_hashes"])
    for row in a["source_pins"].values():
        path = row["path"]; is_input = path.startswith("experiments/root_st073/")
        name = Path(path).name if is_input else path
        target = inputs if is_input else impl
        require(sha(ROOT / path) == row["expected_sha256"], "Changed new input")
        require(target.get(name, row["expected_sha256"]) == row["expected_sha256"], "Conflicting input pin")
        target[name] = row["expected_sha256"]
    impl.update({path: sha(ROOT / path) for path in CODE})
    r = copy.deepcopy(p)
    r.update(stage="actual_core_signed_radial_direction_and_ordinary_Z_jet",
             parent_report=Path(PARENT).name, parent_report_sha256=PARENT_SHA256,
             actual_core_signed_direction=a, input_hashes=inputs, implementation_hashes=impl)
    r["contract"] = dict(actual_core_model_plus_saved_error=True, actual_delta_recipe_not_tube_cap=True,
                         signed_Phi_rho_and_Dbar=True, ordinary_Dbar_Z0_Z1_Z2=True,
                         old_signed_sources_and_error_caps_retained=True, numerical_mean_or_global_closure=False,
                         accepted_ancestors_replayed=False)
    r["checks"] = dict(passed=True, packet_invariants=packet_check(r, p),
                       independent_new_math_review_completed=True,
                       independent_review_scope="Actual core model derivatives, saved same-source total-order3 perturbation, fixed-rho ordinary Dbar Z0..2, actual delta recipe, explicit comparison scope limit",
                       independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
                       independently_reviewed_model_sha256=sha(ROOT / HELPER),
                       independent_runtime_worker_model_metadata_available=False,
                       accepted_ancestors_replayed=False)
    r["execution_seconds"] = round(time.perf_counter() - started, 3)
    return r


def inspect(path):
    r = json.loads(path.read_bytes())
    require(r["checks"]["passed"] and r["checks"]["independent_new_math_review_completed"] and
            r["checks"]["independently_reviewed_module_sha256"] == r["implementation_hashes"][CODE[0]] and
            r["checks"]["independently_reviewed_model_sha256"] == r["implementation_hashes"][HELPER],
            "Same independently reviewed implementation required")
    return dict(passed=True, pins=pins(r), packet_invariants=packet_check(r, parent()),
                serialized_only=True, accepted_owner_norm_recurrence_or_quadrature_replayed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args(); path = DIRECTORY / NAME
    if args.write_report:
        report = run(); tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
