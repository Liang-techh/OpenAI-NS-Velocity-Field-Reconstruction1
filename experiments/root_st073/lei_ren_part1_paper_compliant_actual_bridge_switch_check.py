"""Focused actual R100 source transfer, own-history and original-power checks."""
import json
from pathlib import Path
import mpmath as mp

from lei_ren_part1_paper_compliant_actual_bridge_switch import (
    CompliantActualBridgeSwitch, source_bindings, sha, PREFIX, HERE)
from lei_ren_part1_paper_compliant_actual_bridge_integrals import named_moments
from lei_ren_part1_paper_compliant_first_switch_leading_check import _direction_913
from lei_ren_part1_paper_compliant_inner_switch_profiles import power_transport
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    IntervalTaylor, MTH, MZ, MTHZ, MZT, MP)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode


def contains(a,b):
    return (len(a.coefficients)==len(b.coefficients)
            and all(endpoints(x)[0]<=endpoints(y)[0] and endpoints(x)[1]>=endpoints(y)[1]
                    for x,y in zip(a.coefficients,b.coefficients)))


def independent_postpower(c):
    """Direct positive Volterra integrals, independent of closed-form weights."""
    initial=dict(H="0.9",mean="-0.1",K="0.2",A="0.3",B="0.8",C="0.4")
    jet=lambda v:IntervalTaylor.constant(c,c.mpf(v),5)
    interval=power_transport(c,c.mpf(103)/110,{n:jet(v) for n,v in initial.items()},
                             jet("1.3"),jet("0.7"))
    rows={MTH:("H",2,"2phi"),MZ:("mean",1,"V"),MTHZ:("K",2,"2phiV"),
          "A":("A",1,"V2"),"B":("B",2,"phi2"),MP:("C",1,"phi2")}
    result=[]
    with mp.workdps(120):
        Y=mp.log(mp.mpf(110)/103)
        phi=lambda t:mp.mpf("1.3")*mp.exp(-mp.mpf(".4")*t)
        V=mp.mpf(".7")
        source={"2phi":lambda t:2*phi(t),"V":lambda t:V,
                "2phiV":lambda t:2*phi(t)*V,"V2":lambda t:V*V,
                "phi2":lambda t:phi(t)**2}
        for name,(key,rate,expression) in rows.items():
            actual=mp.exp(-rate*Y)*mp.mpf(initial[key])+mp.quad(
                lambda t:mp.exp(-rate*(Y-t))*source[expression](t),[0,Y])
            row=interval[MZT]["axial" if name=="A" else "swirl"] if name in ("A","B") else interval[name]
            lo,hi=endpoints(row[0])
            miss=max(mp.mpf(0),lo-actual,actual-hi)
            if miss>mp.mpf("1e-70"):
                raise ArithmeticError("Independent postpower integral: "+name)
            result.append(dict(moment=name,direct_integral=actual,enclosure=row[0],
                               positive_enclosure_miss=miss))
    return result


def run():
    with mp.workdps(400):
        name=PREFIX+"actual_bridge_switch.json"
        record=json.loads((HERE/name).read_bytes())
        _verify_hashes(record)
        provider=CompliantActualBridgeSwitch()
        c=provider.ctx
        assert record["source_bindings"]==source_bindings(provider)
        assert (record["actual_five_defect_family_sha256"],record["implicit_source_sha256"],
                record["datum_enclosure_sha256"])==(
                provider.bridge.family,provider.bridge.source,provider.bridge.core.datum.datum_sha)
        jet=lambda values:IntervalTaylor(c,[read_interval(c,v) for v in values])
        finite=joins=modes=0
        for label,packet in record["actual_R100_R110_packets"].items():
            Z=read_interval(c,packet["Z"])
            upstream=provider.upstream.packet(Z,1,"macro")
            inlet=packet["actual_R100_inlet"]
            assert upstream["coordinate_R100_endpoint"]
            assert contains(jet(inlet["F_actual_over_F0_axial5_coefficients"]),IntervalTaylor(c,upstream["actual_phi_axial5"]))
            assert contains(jet(inlet["Uz_actual_axial5_coefficients"]),IntervalTaylor(c,upstream["actual_raw_V_axial5"]))
            mapped=named_moments({n:IntervalTaylor(c,r) for n,r in upstream["actual_own_six_moments_axial5"].items()})
            stored=inlet["actual_moment_shape_axial5_coefficients"]
            for key,row in mapped.items():
                parts=row.items() if isinstance(row,dict) else [(None,row)]
                for part,value in parts:
                    own=jet(stored[key][part] if part is not None else stored[key])
                    assert contains(own,value),label+" actual own moment "+key
                    joins+=6
            assert contains(jet(inlet["pressure_axis_axial5_coefficients"]),IntervalTaylor(c,upstream["pressure_axis_axial5"]))
            assert contains(jet(inlet["pressure_increment_true_axial5_divided_by_R_F0_squared"]),IntervalTaylor(
                c,upstream["actual_pressure_increment_axial5_divided_by_R_F0_squared"]))
            joins+=24
            signed=packet["signed_actual_switch_source_integrals"]
            assert contains(jet(signed["actual_incoming_phi_axial5"]),IntervalTaylor(c,upstream["actual_phi_axial5"]))
            assert contains(jet(signed["actual_incoming_V_axial5"]),IntervalTaylor(c,upstream["actual_raw_V_axial5"]))
            comparison=provider.bridge.comparison(Z,provider.bridge.r/100)
            assert comparison["known_comparison_direction_preserved"]
            inp=provider.bridge.inputs(Z)
            direct=_direction_913(c,Z,provider.bridge.delta,comparison["phi"],comparison["v"],
                                  comparison["moments"],inp["p0"],inp["F0_ratios"],inp["F0_squared_ratios"])
            for key,row in direct.items():
                declared=[jet(r) for r in signed["comparison_direction_modes"][key]]
                assert all(r.order==5 for r in declared)
                assert contains(sum(declared,declared[0]*0),row),label+" independent prescribed direction "+key
                modes+=6
            for endpoint in ("actual_R100_inlet","actual_first_switch_exit","actual_R2_inlet","actual_R110_inlet"):
                state=packet[endpoint]
                lists=[state[k] for k in (
                    "F_actual_over_F0_axial5_coefficients","Uz_actual_axial5_coefficients",
                    "pressure_axis_axial5_coefficients",
                    "pressure_increment_true_axial5_divided_by_R_F0_squared")]
                moments=state["actual_moment_shape_axial5_coefficients"]
                for value in moments.values():
                    lists.extend(value.values() if isinstance(value,dict) else [value])
                assert len(lists)==10
                for values in lists:
                    assert len(values)==6
                    for value in values:
                        lo,hi=endpoints(read_interval(c,value))
                        assert mp.isfinite(lo) and mp.isfinite(hi)
                        finite+=1
                assert len(state["actual_Q_axial4_coefficients"])==5
                assert state["original_P0_retained"] and state["all_actual_moments_inherited"]
                assert not state["comparison_moments_substituted"]
            assert packet["actual_R110_inlet"]["exact_power_moment_transport_from_actual_R2"]
            for flag in ("actual_point_moment_history_recovered","existing_mixed4_provider_replaced",
                         "full_implicit_leading_inputs_recomputed","temporal_recursion"):
                assert packet[flag] is False
            assert packet["same_source_R100_functional_join_installed"]
        assert endpoints(read_interval(c,record["actual_R100_R110_packets"]["whole_Z"]["Z"]))==(mp.mpf(-1),mp.mpf(1))
        for proof in record["bridge_product_cap_proofs"]:
            assert proof["passed"]
            assert endpoints(read_interval(c,proof["source_log_upper"]))[1]<=endpoints(c.ln(read_interval(c,proof["cap"])))[0]
        for proof in record["signed_switch_product_cap_proofs"]:
            assert proof["exact_source_not_replaced"]
            assert endpoints(read_interval(c,proof["log_magnitude_upper"]))[1]<=endpoints(read_interval(c,proof["log_cap"]))[0]
        fixtures=independent_postpower(c)
        hashes=dict(record["input_hashes"])
        hashes[name]=sha(name)
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        dep=PREFIX+"first_switch_leading_check.py"
        hashes[dep]=sha(dep)
        result=dict(
            all_passed=True,
            actual_five_defect_family_sha256=provider.bridge.family,
            implicit_source_sha256=provider.bridge.source,
            datum_enclosure_sha256=provider.bridge.core.datum.datum_sha,
            current_input_hashes_checked=len(hashes),
            original_callable_identities_checked=len(record["source_bindings"]["original_callable_identities"]),
            comparison_continuation_AST_bindings_checked=len(record["source_bindings"]["comparison_continuation_AST_bindings"]),
            actual_R100_field_moment_pressure_coefficients_transferred=joins,
            independent_comparison_direction_coefficients_checked=modes,
            finite_actual_R100_R2_R110_axial_coefficients_checked=finite,
            independent_six_postpower_integrals=fixtures,
            whole_real_Z_R100_R110_composition_checked=True,
            R100_R110_actual_feedback_composed=True,
            actual_point_moment_history_recovered=False,
            existing_mixed4_provider_replaced=False,
            full_implicit_leading_inputs_recomputed=False,
            global_completed_tensor_admissibility=False,
            temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix(".json").write_text(json.dumps(_encode(result),indent=2)+"\n",encoding="utf8")
        print("PASS actual finite-width R100-to-R110 own-history composition",flush=True)
        return result


if __name__=="__main__":
    run()
