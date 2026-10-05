"""Independent receipt check for the first R=100 micro-chart leading term."""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles, IntervalTaylor, derivative, dress, square,
    MTH, MZ, MTHZ, MZT, MP,
)
from lei_ren_part1_paper_compliant_macro_signed_integrals import (
    COMPARISON_NAME, COMPARISON_CHECK_NAME, _packet_for_coordinate,
    _packet_phi, _packet_v, _packet_moments_correct, _verify_hashes,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_compliant_"
NAME = PREFIX + "first_switch_leading.json"
SOURCE = PREFIX + "first_switch_leading.py"
COMPARISON_SOURCE = PREFIX + "comparison_point_integrals.json"


def _hash(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _jet(c, row):
    from lei_ren_part1_paper_compliant_inner_bridge_profiles import IntervalTaylor
    return IntervalTaylor(c, [read_interval(c, value) for value in row])


def _contains(stored, expected):
    return all(endpoints(a)[0] <= endpoints(b)[0]
               and endpoints(a)[1] >= endpoints(b)[1]
               for a, b in zip(stored.coefficients, expected.coefficients))


def _sigma(x):
    if x <= 0: return mp.mpf(0)
    if x >= 1: return mp.mpf(1)
    if x > mp.mpf("0.5"): return 1 - _sigma(1-x)
    odds = 1/(1-x)**2 - 1/x**2
    e = mp.exp(odds)
    return e/(1+e)


def _direction_913(c, zvalue, delta, phi, value, moments, pressure, ratios, ratios2):
    """Independent transcription of the four original (9.13) direction rows."""
    z=IntervalTaylor.variable(c,zvalue,phi.order); d=1-square(z); L=1-square(z)*delta
    m=moments[MZ]; a=moments[MZT]["axial"]
    f=dress(phi,ratios); h=dress(moments[MTH],ratios); k=dress(moments[MTHZ],ratios)
    b=dress(moments[MZT]["swirl"],ratios2); p=dress(moments[MP],ratios2)
    W=1-(z*m)*(1-delta)-d*derivative(m)
    angular=h*(1-delta/2)-(z*derivative(h))*((1-delta)/2)-d*derivative(k)+(z*k)*(2*delta-1)
    return dict(
        D_over_R=(-W+angular/(f*2))/L,
        drive_hydro=(-W*value+(m-z*derivative(m))*((1-delta)/2)+(z*a)*(2*delta)-d*derivative(a))/(2*L),
        drive_pressure=((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L),
        drive_swirl=(-(z*b)*(2*delta)+d*derivative(b)+(z*p)*(2*(1+delta))-d*derivative(p))/(2*L))


def run():
    bridge = CompliantInnerBridgeProfiles(); c = bridge.ctx
    record = json.loads((HERE/NAME).read_bytes())
    comparison = json.loads((HERE/COMPARISON_SOURCE).read_bytes())
    comparison_check = json.loads((HERE/COMPARISON_CHECK_NAME).read_bytes())
    _verify_hashes(record); _verify_hashes(comparison); _verify_hashes(comparison_check)
    failures=[]
    def require(ok, msg):
        if not ok: failures.append(msg)
    require(record["actual_five_defect_family_sha256"]==bridge.family,"family hash")
    require(record["implicit_source_sha256"]==bridge.source,"implicit source hash")
    require(record["datum_enclosure_sha256"]==bridge.core.datum.datum_sha,"pressure datum hash")
    require(record["leading_hb2_first_switch_only"] is True,"scope flag")
    require(record["full_first_switch_resolved"] is False,"finite-width scope")
    require(record["actual_signed_bridge_completed"] is False,"bridge scope")
    require(record["exact_first_switch"]=="first micro chart: R=100*exp(hb*s), 0<=s<=1","original switch radius")

    with mp.workdps(200):
        pulse = mp.quad(lambda s: 1-_sigma(s), [0, mp.mpf("0.5"), 1])
    weight = read_interval(c, record["packets"][".5"]["exact_pulse_weight"])
    wlo, whi = endpoints(weight)
    require(wlo <= pulse <= whi, "independent pulse quadrature")
    require(abs(pulse-mp.mpf("0.5")) < mp.mpf("1e-150"), "exact symmetric pulse integral")

    from lei_ren_part1_paper_compliant_first_switch_leading import switch_control_source_bridge
    require(record["source_control_bridge"]==switch_control_source_bridge(),"original unweighted angular control")
    tree=ast.parse((HERE/(PREFIX+"microswitch_mixed_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="switch_controls")
    first=next(n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test)=="branch == 'first'")
    require(ast.unparse(first.body[0])=="a = tinyD","independent first angular branch")
    checked=[]
    for label in (".5","0","exact_shared_root"):
        packets=comparison["comparison_point_packets"][label]["macro"]
        endpoint=_packet_for_coordinate(c,packets,mp.mpf(1))
        z=read_interval(c,endpoint["Z"]);phi=_packet_phi(c,endpoint);v=_packet_v(c,endpoint)
        moments=_packet_moments_correct(c,endpoint);inputs=bridge.inputs(z)
        current=_direction_913(c,z,bridge.delta,phi,v,moments,inputs["p0"],inputs["F0_ratios"],inputs["F0_squared_ratios"])
        saved=record["packets"][label]
        expected_f=current["D_over_R"]*(-c.mpf(50))
        require(endpoints(read_interval(c,saved["exact_angular_first_chart_weight"]))==(mp.mpf(1),mp.mpf(1)),label+" unweighted angular source")
        saved_f=_jet(c,saved["J_logF_hb2"])
        require(_contains(saved_f,expected_f),label+" log-F coefficient")
        for component, source, radius_power, scale_log in (
            ("hydro","drive_hydro",100,c.mpf(0)),
            ("pressure","drive_pressure",100,2*bridge.core.logP),
            ("swirl","drive_swirl",10000,c.mpf([endpoints(-2*bridge.core.logC-2*bridge.core.Lambda*bridge.core.Gbar)[0],
                                                   endpoints(-2*bridge.core.logC)[1]]))):
            item=saved["J_V_"+component]
            coeff=_jet(c,item["normalized_axial_coefficients"])
            expected=current[source]*(-c.mpf("0.5")*radius_power)
            require(_contains(coeff,expected),label+" "+component+" coefficient")
            lo,hi=endpoints(read_interval(c,item["positive_scale_log"]))
            elo,ehi=endpoints(scale_log)
            require(lo<=elo and hi>=ehi,label+" "+component+" full source scale")
            require(item.get("positive_scale_log_is_enclosure") is True,label+" source log enclosure")
        checked.append(label)

    result=dict(all_passed=not failures,failures=failures,
        leading_hb2_first_switch_checked=True,original_angular_first_weight1_checked=True,independent_sigma_integral=pulse,
        points_checked=checked,full_F0_squared_source_log_range_checked=True,full_first_switch_resolved=False,
        actual_signed_bridge_completed=False,
        input_hashes={NAME:_hash(NAME),SOURCE:_hash(SOURCE),
                      Path(__file__).name:_hash(Path(__file__).name),
                      COMPARISON_SOURCE:_hash(COMPARISON_SOURCE),
                      COMPARISON_NAME:_hash(COMPARISON_NAME),
                      COMPARISON_CHECK_NAME:_hash(COMPARISON_CHECK_NAME),
                      PREFIX+"macro_signed_integrals.py":_hash(PREFIX+"macro_signed_integrals.py"),
                      PREFIX+"flat_pulse_derivatives.py":_hash(PREFIX+"flat_pulse_derivatives.py"),
                      PREFIX+"microswitch_mixed_C4.py":_hash(PREFIX+"microswitch_mixed_C4.py")})
    (HERE/(Path(__file__).stem+".json")).write_text(json.dumps(_encode(result),indent=2)+"\n",encoding="utf8")
    print("First-switch leading checker:","PASS" if result["all_passed"] else "FAIL",flush=True)
    if failures: raise RuntimeError("; ".join(failures))
    return result


if __name__=="__main__": run()
