"""Current-source switch mixed4 check using unchanged admitted row algorithms.

Only source-graph ownership and provider/receipt names change in the local
verifier replay. Original formula fixtures are consumed by current hashes;
their finite scalar parameters are never treated as production field data.
"""
import ast
import copy
import json
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_microswitch_mixed_C4_check as template
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import (
    CompliantActualSwitchMixedC4, ActualSwitchSourceDispatcher,
    R100_source_join, mixed_source_bindings, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def replay_rows(admitted,holder):
    path=Path(template.__file__)
    tree=ast.parse(path.read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="run")
    expected=ast.dump(ast.parse(
        "json.loads((HERE/'lei_ren_part1_paper_compliant_long_reshape_mixed_C4.json').read_bytes())['shared_exact_axial_source']",
        mode="eval").body)
    changed=0
    for node in ast.walk(fn):
        if (isinstance(node,ast.Assign) and any(ast.unparse(v)=="shared" for v in node.targets)
                and ast.dump(node.value)==expected):
            node.value=ast.parse("raw['shared_exact_axial_source']",mode="eval").body
            changed+=1
    if changed!=1:
        raise ValueError("Original row verifier source-graph ownership assignment changed")
    if (len(fn.body)<3 or not isinstance(fn.body[-1],ast.Return)
            or not isinstance(fn.body[-3],ast.Expr)
            or "write_text" not in ast.unparse(fn.body[-3])
            or not isinstance(fn.body[-2],ast.Expr)
            or not ast.unparse(fn.body[-2]).startswith("print(")):
        raise ValueError("Original verifier publication tail changed")
    # Publication belongs to this checker after all added source gates pass.
    fn.body=fn.body[:-3]+[fn.body[-1]]
    fn.name="_verify_current_rows"
    env=dict(template.__dict__)
    env["__file__"]=__file__
    env["NAME"]=PREFIX+"actual_switch_mixed_C4.json"
    def factory():
        holder["provider"]=CompliantActualSwitchMixedC4()
        return holder["provider"]
    env["CompliantMicroswitchMixedC4"]=factory
    fixture_keys=("comparison_fixture","controls_fixture","physical_fixture","factored_scale_fixture","functional_identities")
    for name in fixture_keys:
        if not admitted[name]["passed"]:
            raise ValueError("Admitted original generic formula evidence required: "+name)
        env[name]=lambda name=name:copy.deepcopy(admitted[name])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(path),"exec"),env)
    return env["_verify_current_rows"]()


def run():
    with mp.workdps(400):
        original_name=PREFIX+"microswitch_mixed_C4_check.json"
        admitted=json.loads((HERE/original_name).read_bytes())
        _verify_hashes(admitted)
        if not admitted["all_passed"]:
            raise ValueError("Original unchanged mixed4 operator/fixture evidence required")
        name=PREFIX+"actual_switch_mixed_C4.json"
        raw=json.loads((HERE/name).read_bytes())
        _verify_hashes(raw)
        holder={}
        result=replay_rows(admitted,holder)
        provider=holder["provider"]
        if raw["current_actual_mixed4_source_bindings"]!=mixed_source_bindings(provider):
            raise ValueError("Current original-math/source provider binding changed")
        if raw["actual_R100_functional_mixed4_join"]!=R100_source_join():
            raise ValueError("Original R100 functional mixed4 join changed")
        if raw["shared_exact_axial_source"]!=provider.shared_axial_source:
            raise ValueError("Current defining axial graph differs from actual provider")
        graph=raw["shared_exact_axial_source"]
        if (not graph["finite_width_source_integral_enclosures_available"]
                or graph["formal_integrals_numerically_reconstructed"]
                or not graph["finite_width_source_integral_bindings"]["known_comparison_direction_preserved"]):
            raise ValueError("Source-enclosure graph promoted or known direction replaced")
        # The two source references are shared exactly with centered E.
        args=graph["V110"]["args"]
        if graph["V100"]["args"]!=args[:-1] or graph["E_V110_minus_4Z"]["args"]!=args[1:]:
            raise ValueError("Centered E no longer shares the actual V110 defining integrals")
        trace_count=0
        current=provider.switch.bridge.actual([-1,1],provider.switch.bridge.r/100)
        inlet=raw["actual_R100_inlet"]["actual_parent_axial5_packet"]
        expected=encode(pack(current))
        fields=(
            "F_actual_over_F0_axial5_coefficients","F_actual_true_axial5_divided_by_F0",
            "Uz_actual_axial5_coefficients","actual_moment_shape_axial5_coefficients",
            "actual_Q_axial4_coefficients","pressure_axis_axial5_coefficients",
            "pressure_increment_true_axial5_divided_by_R_F0_squared")
        for key in fields:
            if canonical_source(inlet[key])!=canonical_source(expected[key]):
                raise ValueError("Current actual R100 trace data changed: "+key)
            trace_count+=1
        for flag in ("actual_point_moment_history_recovered","downstream_reshape_actual_feedback_installed",
                     "full_implicit_leading_inputs_recomputed","global_completed_tensor_admissibility",
                     "temporal_recursion"):
            if raw[flag] is not False:
                raise ValueError("Unbuilt scope promoted: "+flag)
        # Exercise the existing dispatcher evaluation with this checked provider.
        # Receipt publication follows every check, so no provisional PASS exists.
        dispatcher=ActualSwitchSourceDispatcher()
        dispatcher.providers["actual_switch_mixed_C4"]=provider
        dispatcher.family=provider.family
        dispatcher.source=provider.source
        dispatched=0
        for chart,coordinate,key in (
                ("switch_first",[0,1],"whole_first"),
                ("switch_second",[1,2],"whole_second"),
                ("switch_power",[0,1],"whole_postswitch_power")):
            packet=dispatcher.evaluate(chart,[-1,1],coordinate)
            if (not packet["current_finite_width_actual_history_used"]
                    or packet["acceptance_receipt"]!=PREFIX+"actual_switch_mixed_C4_check.json"):
                raise ValueError("Dispatcher returned a legacy switch receipt")
            for group,rows in packet["physical_mixed_grids"].items():
                if canonical_source(encode(pack(rows)))!=canonical_source(raw[key][group]):
                    raise ValueError("Dispatcher did not return the current mixed4 grids: "+chart)
                for grid in rows.values():
                    for value in grid.values():
                        lo,hi=endpoints(value)
                        if not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current dispatched derivative")
                        dispatched+=1
        hashes=dict(result["input_hashes"])
        hashes[name]=sha(name)
        hashes[original_name]=sha(original_name)
        hashes[Path(template.__file__).name]=sha(Path(template.__file__).name)
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result.update(
            datum_enclosure_sha256=provider.core.datum.datum_sha,
            original_generic_fixture_receipt_consumed=original_name,
            original_row_verifier_replayed=True,
            verifier_replay_changes=["current provider and receipt owner","current exact source-graph owner",
                                     "hash-bound generic fixture reuse","publication only after added gates"],
            current_input_hashes_checked=len(hashes),
            current_actual_mixed4_input_AST_bindings_checked=len(
                raw["current_actual_mixed4_source_bindings"]["original_mixed4_actual_input_AST_bindings"]),
            current_actual_field_AST_bindings_checked=len(
                raw["actual_R100_functional_mixed4_join"]["actual_R100_trace_AST_bindings"]["current_actual_field_AST_bindings"]),
            actual_R100_field_moment_Q_pressure_trace_groups_checked=trace_count,
            actual_R100_functional_mixed4_join=raw["actual_R100_functional_mixed4_join"],
            actual_R100_functional_mixed4_join_certified=True,
            actual_finite_width_bridge_histories_installed_in_switch_mixed4=True,
            actual_R100_R110_feedback_mixed4_available=True,
            current_actual_source_dispatcher_charts_checked=3,
            current_actual_source_dispatcher_derivative_rows_checked=dispatched,
            actual_point_moment_history_recovered=False,
            downstream_reshape_actual_feedback_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            global_completed_tensor_admissibility=False,
            temporal_recursion=False,all_passed=True,input_hashes=hashes)
        Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
        print("PASS actual finite-width switch mixed4, current dispatcher and functional R100 join",flush=True)
        return result


if __name__=="__main__":
    run()
