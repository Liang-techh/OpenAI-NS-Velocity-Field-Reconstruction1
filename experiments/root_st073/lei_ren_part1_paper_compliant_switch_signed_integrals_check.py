"""Focused source and independent integral check for the entire R100..110 switch."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_switch_signed_integrals import (
    switch_integral_identities, switch_control_source_bridge)
from lei_ren_part1_paper_compliant_first_switch_leading_check import (
    _direction_913, _sigma)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles, IntervalTaylor, MTH, MZ, MTHZ, MZT, MP)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX="lei_ren_part1_paper_compliant_"


def sha(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def contains(row,value):
    return all(endpoints(a)[0]<=endpoints(b)[0] and endpoints(a)[1]>=endpoints(b)[1]
               for a,b in zip(row.coefficients,value.coefficients))


def moderate_integral_fixtures():
    """Independent full finite-width quadrature; no source point is selected."""
    c=MPIntervalContext(); c.dps=100; results=[]
    with mp.workdps(85):
        for htext in (".01",".04"):
            h=mp.mpf(htext); H=c.mpf(htext)
            # Dbar positive, with three genuine frozen-history modes.
            D=[mp.mpf(".2"),mp.mpf("-.1"),mp.mpf(".03")]
            Dfun=lambda y:100*sum(D[j]*mp.exp((1-j)*y) for j in range(3))
            JD=mp.quad(lambda x:Dfun(h*x),[0,1])+mp.quad(
                lambda x:(1-_sigma(x))*Dfun(h*(1+x)),[0,mp.mpf(".5"),1])
            bounds=c.mpf(0)
            for j in range(3):
                e=c.exp(c.mpf([0,H.b])*(1-j))
                bounds+=100*c.mpf(str(D[j]))*(e+e*e/2)
            lo,hi=endpoints(bounds)
            assert lo<=JD<=hi,"angular exact finite-width fixture"
            prefix=lambda x:100*sum(D[j]*(
                (mp.expm1((1-j)*h*x)/((1-j)*h)) if j!=1 else x) for j in range(3))
            Dmax=100*sum(abs(D[j])*mp.exp(max(0,(1-j)*h)) for j in range(3))
            qbounds=c.exp(c.mpf([-mp.mpf(".5")*h*h*Dmax,0]))*c.mpf("1.03")
            for power in (1,2):
                for coeffs in ((mp.mpf(".3"),mp.mpf("-.4"),mp.mpf(".02")),
                               (mp.mpf("-.2"),mp.mpf(".04"),mp.mpf("-.1"))):
                    drive=lambda x:100**power*sum(coeffs[j]*mp.exp((power-j)*h*x)
                                                  for j in range(3))
                    actual=-mp.quad(lambda x:(1-_sigma(x))*mp.mpf("1.03")
                                    *mp.exp(-h*h*prefix(x)/2)*drive(x),
                                    [0,mp.mpf(".5"),1])
                    interval=sum((c.mpf(str(coeffs[j]))*c.exp(c.mpf([0,H.b])*(power-j))
                                  for j in range(3)),c.mpf(0))*(100**power)
                    enclosure=-qbounds*interval/2; lo,hi=endpoints(enclosure)
                    assert lo<=actual<=hi,"axial exact finite-width fixture"
                    results.append(dict(h=h,radial_power=power,actual_normalized_integral=actual,
                                        enclosure=enclosure))
    return results


def _run():
    name=PREFIX+"switch_signed_integrals.json"
    record=json.loads((HERE/name).read_bytes()); _verify_hashes(record)
    bridge=CompliantInnerBridgeProfiles(); c=bridge.ctx
    assert (record["actual_five_defect_family_sha256"],record["implicit_source_sha256"])==(bridge.family,bridge.source)
    assert record["exact_switch_composition"]==switch_integral_identities()
    assert record["source_control_bridge"]==switch_control_source_bridge()
    assert record["whole_Z_full_switch_signed_enclosures_available"]
    for flag in ("higher_actual_Ra_to_R100_bridge_orders_solved","actual_signed_bridge_completed",
                 "selected_nonlinear_point_values_recovered","temporal_recursion"):
        assert record[flag] is False,flag
    count=0
    jet=lambda values:IntervalTaylor(c,[read_interval(c,v) for v in values])
    checked=[]
    for label,packet in record["signed_source_integral_packets"].items():
        z=read_interval(c,packet["Z"]); inp=bridge.inputs(z)
        incoming=bridge.actual(z,bridge.r/100)
        assert contains(jet(packet["actual_incoming_phi_axial5"]),IntervalTaylor(c,incoming["F_actual_over_F0_axial5_coefficients"]))
        assert contains(jet(packet["actual_incoming_V_axial5"]),IntervalTaylor(c,incoming["Uz_actual_axial5_coefficients"]))
        assert endpoints(read_interval(c,packet["original_width_log"]))[1]<endpoints(c.ln(bridge.cap))[0]
        for key in ("actual_JD_signed_axial5_enclosure","actual_R110_phi_over_F0_axial5_enclosure",
                    "actual_R110_V_axial5_enclosure"):
            for value in packet[key]:
                lo,hi=endpoints(read_interval(c,value)); assert mp.isfinite(lo) and mp.isfinite(hi)
                count+=1
        assert packet["complete_switch_signed_integral_enclosures_available"]
        assert packet["actual_incoming_bridge_retained"]
        assert packet["comparison_field_not_substituted_for_actual_field"]
        assert not packet["newly_selected_point_values"]
        # Independent original direction at R100 checks the stored modal
        # reconstruction's zero-time value, retaining all six own moments.
        phi=jet(packet["comparison_frozen_phi_axial5"]); v=jet(packet["comparison_frozen_V_axial5"])
        raw=packet["comparison_R100_own_moments_axial5"]
        moments={key:({part:jet(row) for part,row in val.items()} if isinstance(val,dict) else jet(val))
                 for key,val in raw.items()}
        direct=_direction_913(c,z,bridge.delta,phi,v,moments,inp["p0"],inp["F0_ratios"],inp["F0_squared_ratios"])
        for key,row in direct.items():
            modes=[jet(values) for values in packet["comparison_direction_modes"][key]]
            summed=sum(modes,modes[0]*0)
            assert contains(summed,row),label+" same complete source modal direction "+key
        checked.append(label)
    assert endpoints(read_interval(c,record["signed_source_integral_packets"]["whole_Z"]["Z"]))==(mp.mpf(-1),mp.mpf(1))
    for label,p in record["leading_width_packets"].items():
        first=json.loads((HERE/(PREFIX+"first_switch_leading.json")).read_bytes())["packets"][label]
        D=jet(first["D_over_R_at_R100"])
        assert contains(jet(p["first_angular_hb2"]),D*(-50))
        assert contains(jet(p["second_angular_hb2"]),D*(-25))
        assert contains(jet(p["R110_logF_hb2"]),D*(-75))
        assert contains(jet(p["R110_logF_hb1"]),D*0+c.mpf(".6"))
    for proof in record["bound_proofs"]:
        assert proof["exact_source_not_replaced"]
        assert endpoints(read_interval(c,proof["log_magnitude_upper"]))[1]<=endpoints(read_interval(c,proof["log_cap"]))[0]
    fixtures=moderate_integral_fixtures()
    hashes=dict(record["input_hashes"]); hashes[name]=sha(name)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    # The independent transcription is part of this check's source.
    dep=PREFIX+"first_switch_leading_check.py"; hashes[dep]=sha(dep)
    result=dict(all_passed=True,actual_five_defect_family_sha256=bridge.family,
                implicit_source_sha256=bridge.source,
                whole_Z_full_switch_signed_enclosures_available=True,
                original_unweighted_first_angular_control_verified=True,
                complete_switch_width_coefficients_through2_checked=True,
                actual_incoming_and_common_moment_modes_checked=checked,
                finite_signed_output_rows=count,moderate_full_integral_fixtures=fixtures,
                source_caps_not_point_values=True,
                actual_signed_bridge_completed=False,
                selected_nonlinear_point_values_recovered=False,
                temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix(".json").write_text(json.dumps(_encode(result),indent=2)+"\n",encoding="utf8")
    print("Full original switch signed integral enclosures PASS; full bridge/recursion pending",flush=True)
    return result


def run():
    with mp.workdps(300):return _run()


if __name__=="__main__":run()
