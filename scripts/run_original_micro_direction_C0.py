"""Save or inspect actual two-micro C0 radii; never replay ancestors."""
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
from openai_ns_reconstruction.original_micro_direction_C0 import build, PARENT, PARENT_SHA256, HELPER, require, sha, digest

NAME = "lei_ren_part1_paper_compliant_current_original_micro_direction_C0.json"
CODE = ("src/openai_ns_reconstruction/original_micro_direction_C0.py",
        "scripts/run_original_micro_direction_C0.py", HELPER)


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
    require(r['parent_report'] == Path(PARENT).name and r['parent_report_sha256'] == PARENT_SHA256, 'Same signed core predecessor required')
    mutable = {'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes','contract','checks','execution_seconds'}
    for k,v in p.items():
        if k not in mutable: require(r[k] == v, 'Inherited signed packet or error changed: '+k)
    a=r['actual_micro_direction_C0']
    require(a['passed'] and a['source_family'] == p['source_family'] and a['ordinary_Z_orders'] == [0] and
            a['both_micro_actual_a_and_E_positive'] and a['actual_finiteN_and_joint_mean_C0_enclosed'], 'Actual two-micro C0 gates missing')
    require(a['comparison_direction']['two_micro_comparison_positive'] and
            a['comparison_direction']['source_bound_from_own_comparison_histories'] and
            not a['comparison_direction']['core_and_comparison_fields_identified_on_second_micro'] and
            a['comparison_algebra']['F0_dressing_and_epsilon_Lambda_identity'], 'Comparison is its own source, with exact scale cancellation')
    e=a['actual_C0_errors']
    require(e['passed'] and e['ordinary_Z_order'] == 0 and e['common_coarse_error_radius']=='2^-1024' and
            e['eta_positive_and_below_one'] and e['F0_positive_and_at_most_one'] and
            e['joint_mean_incoming_history_retained'] and e['zero_is_error_box_center_not_function_value'] and
            len(e['logarithmic_error_radii']) == 5, 'Source-bound relative/physical/mean error radii required')
    require(a['universal_loop_controls']['passed'] and
            a['universal_loop_controls']['phase_normalization']['numerical_quadrature'] is False and
            a['directed_relative_errors']['N_materialized'] is False, 'Phase normalization and nonmaterialized frequency required')
    gates=a['actual_root_domain_gates']
    require(gates['passed'] and gates['s']==[0,2] and gates['hb_strictly_positive_and_below_one'] and
            gates['active_implies_Delta_less_eta'] and gates['flat_implies_A_and_B_exact_zero'], 'Actual root branch/domain gates required')
    require(a['loop_cut_variable_binding']['helper_cut']=='sigma_flat(1-Delta/eta)' and
            a['loop_cut_variable_binding']['source_argument_itself_not_assumed_below_one'], 'Helper cut denotes the normalized sigma value')
    for k in ('full_C5_numeric_containment','numerical_Cy_Z3','numerical_QZZ','forward_Rp_discrepancy_enclosed',
              'global_five_moment_closure','accepted_ancestors_replayed','full_goal_complete'):
        require(a[k] is False, 'Unproved higher derivative or physical closure: '+k)
    require(all(v is False for v in r['limits'].values()), 'Physical completion gate changed')
    return dict(passed=True, both_comparison_directions_positive=True, actual_a_and_E_positive=True,
                actual_relative_physical_joint_mean_C0_radii=True, original_joint_memory_retained=True,
                all_inherited_packets_and_errors_unchanged=True, no_second_core_replacement=True,
                accepted_ancestors_replayed=False, higher_Z_derivative_numeric_enclosure=False)


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
    r.update(stage="actual_two_micro_comparison_direction_and_finite_N_C0_radii",
             parent_report=Path(PARENT).name, parent_report_sha256=PARENT_SHA256,
             actual_micro_direction_C0=a, input_hashes=inputs, implementation_hashes=impl)
    r['contract'] = dict(actual_two_micro_comparison_direction_positive=True, actual_root_a_and_E_positive=True,
                         relative_and_physical_finite_N_C0_radii=True, joint_mean_C0_memory_retained=True,
                         old_signed_sources_and_error_caps_retained=True, numerical_Z3_or_global_closure=False,
                         accepted_ancestors_replayed=False)
    r["checks"] = dict(passed=True, packet_invariants=packet_check(r, p),
                       independent_new_math_review_completed=True,
                       independent_review_scope="Actual own-comparison direction bounds, exact F0/Lambda cancellation, source-bound loop C0 theorem, physical finiteN and nonreset mean radii; higher Z rows open",
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
