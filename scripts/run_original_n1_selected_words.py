"""Directed raw n1 finite-word transport at the actual microscopic collar."""
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
from openai_ns_reconstruction.original_finite_N_signed_source import digest, load, require, sha
from openai_ns_reconstruction.original_n1_selected_word_numeric import build as build_numeric
from openai_ns_reconstruction.original_n1_selected_word_tail import build as build_tail
from openai_ns_reconstruction.original_n1_cutoff_measure import build as build_measure

PARENT = "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
PARENT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"
NAME = "lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"
CODE = (
    "src/openai_ns_reconstruction/original_n1_selected_word_numeric.py",
    "src/openai_ns_reconstruction/original_n1_selected_word_tail.py",
    "src/openai_ns_reconstruction/original_n1_cutoff_measure.py",
    "scripts/run_original_n1_selected_words.py",
)


def parent():
    require(sha(DIRECTORY / PARENT) == PARENT_SHA256, "Accepted cutoff source parent changed")
    p = load(DIRECTORY / PARENT)
    require(p["checks"]["passed"] is True, "Passed source parent required")
    return p


def pins(r):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for n, h in r[key].items():
            require(sha(base / n) == h, "Changed source pin: " + n)
            count += 1
    return count


def merge_pins(p, packets):
    inputs, impl = dict(p["input_hashes"]), dict(p["implementation_hashes"])
    inputs[PARENT] = PARENT_SHA256
    for packet in packets:
        for row in packet["source_pins"].values():
            require(isinstance(row, dict) and isinstance(row.get("path"), str), "Invalid source pin")
            n = row["path"]
            h = row.get("expected_sha256") or row.get("sha256")
            require(isinstance(h, str) and sha(ROOT / n) == h, "Changed new source pin: " + n)
            if n.startswith("experiments/root_st073/"):
                n = n[len("experiments/root_st073/"):]
                target = inputs
            else:
                target = impl
            require(target.get(n, h) == h, "Conflicting source pin: " + n)
            target[n] = h
    return inputs, impl


def transport_algebra():
    """Check the exact homogeneous six-row seam flow, including K1 sign."""
    import sympy as sp
    x = sp.symbols("x", positive=True)
    f, u, k, p, df, du = sp.symbols("f u k p df du")
    w = [f + (1 - 4 / x**2) * df, u + 2 * sp.log(x / 2) * du,
        4 * k / x**2 - (1 - 4 / x**2) * du, p, 8 * df / x**3, 2 * du / x]
    D = [0, 0, 2, 0, 3, 1]
    akin = [w[4], w[5], -w[5], 0, 0, 0]
    equations = [sp.simplify(sp.diff(v, x) + D[i] * v / x - akin[i]) == 0
        for i, v in enumerate(w)]
    seam = [sp.simplify(v.subs(x, 2) - a) == 0
        for v, a in zip(w, (f, u, k, p, df, du))]
    require(all(equations + seam), "Homogeneous no-reset transport failed")
    return dict(passed=True, six_raw_ODE_rows=equations, six_core_inlet_rows=seam,
        uses_x_equals_2_exp_hb_s_over_2=True, ordinary_Z_commutes_with_Z_independent_map=True)


def packet_check(r, p):
    require(r["parent_report"] == PARENT and r["parent_report_sha256"] == PARENT_SHA256,
        "Different source parent")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes",
        "implementation_hashes", "contract", "checks", "execution_seconds"}
    for k, v in p.items():
        if k not in mutable:
            require(r[k] == v, "Inherited source/error packet changed: " + k)
    a = r["actual_selected_first_word_enclosure"]
    t = r["actual_selected_raw_word_tail"]
    m = r["actual_cutoff_complete_support_measure"]
    require(all(q["passed"] is True and q["source_executed"] is False and
        q["accepted_ancestors_replayed"] is False for q in (a, t, m)), "New packet failed")
    require(a["actual_selected_finite_word_enclosures_evaluated"] is True and
        a["actual_directed_sector_coefficients_and_local_radii_evaluated"] is True and
        a["source_only"] is False and a["ordinary_Z_orders"] == [0, 1] and
        [q["q"] for q in a["selected_coordinate_rows"]] == ["0", "1/2", "5/8", "3/4", "1"] and
        a["parent_mapping_sha256"] == digest(p) and
        a["old_actual_first_keep_point_values_used"] is False and
        a["all_Z_functional_coverage"] is False and
        a["signed_C_inner_m15_value_evaluated"] is False,
        "Actual selected-coordinate finite-word enclosure admission failed")
    require(t["actual_N1_tail_evaluated"] is True and t["ordinary_Z_orders"] == [0, 1],
        "Actual raw ordinary-Z tail required")
    require(m["complete_first_collar_support_included"] is True and
        m["complete_core_plus_first_measure"] == "J_total=4+hb*sc*C" and
        m["numeric_m15_error_budget_instantiated"] is False, "Incomplete support measure")
    require(all(v is False for v in r["limits"].values()), "Global scope changed")
    require(r["selected_transport_exact_algebra"]["passed"] is True and
        all(r["selected_transport_exact_algebra"]["six_raw_ODE_rows"]) and
        all(r["selected_transport_exact_algebra"]["six_core_inlet_rows"]), "Seam ODE checks missing")
    return dict(passed=True, inherited_sources_and_all_errors_unchanged=True,
        raw_tail_is_bound_not_field_value=True, complete_cutoff_measure_includes_core=True,
        selected_point_is_not_functional_five_moment_closure=True,
        numerical_C_inner_m15_value=False, old_caps_or_signed_terms_removed=False,
        accepted_ancestors_replayed=False)


def run():
    started = time.perf_counter()
    p = parent()
    pins(p)
    a, t, m = build_numeric(ROOT, p), build_tail(ROOT, p), build_measure(ROOT, p)
    inputs, impl = merge_pins(p, (a, t, m))
    r = copy.deepcopy(p)
    r.update(stage="actual_selected_first_n1_finite_word_directed_transport",
        parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        actual_selected_first_word_enclosure=a, actual_selected_raw_word_tail=t,
        actual_cutoff_complete_support_measure=m, input_hashes=inputs,
        selected_transport_exact_algebra=transport_algebra(),
        implementation_hashes=dict(impl, **{n: sha(ROOT / n) for n in CODE}))
    r["contract"] = dict(actual_selected_radial_coordinate="s=sc*q/4",
        directed_signed_coefficient_approximation_and_additive_error_computed=True,
        tail_is_raw_absolute_error_not_signed_value=True,
        source_scale_factors_retained=True, complete_core_plus_first_eta_measure=True,
        signed_C_inner_m15_value=False, global_five_moment_closure=False,
        old_six_signed_terms_and_all_errors_retained=True, accepted_ancestors_replayed=False)
    r["checks"] = dict(passed=True, packet_invariants=packet_check(r, p),
        accepted_ancestors_replayed=False)
    r["execution_seconds"] = round(time.perf_counter() - started, 3)
    return r


def inspect(path):
    r = load(path)
    require(r["checks"]["passed"] is True, "Saved passed report required")
    return dict(passed=True, pins=pins(r), packet_invariants=packet_check(r, parent()),
        serialized_only=True, owner_fixture_recurrence_norm_or_quadrature_replayed=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = ap.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        r = run()
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=r["execution_seconds"],
            bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
