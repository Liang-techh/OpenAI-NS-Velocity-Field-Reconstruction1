"""Leading h_b^2 signed correction across the first R=100 micro chart.

The switch is R=100 exp(h_b s), 0<=s<=1.  At leading h_b^2, its smooth
axial cutoff integrates exactly to 1/2, while the angular weight is1; the actual/comparison quotient and current
direction are their R=100 zero-width endpoint values.  All pressure and swirl
amplitudes remain logarithmic/factored.  Higher orders and the complete
finite-width switch are deliberately unresolved.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles, direction,
)
from lei_ren_part1_paper_compliant_macro_signed_integrals import (
    COMPARISON_NAME, COMPARISON_CHECK_NAME, _packet_for_coordinate,
    _packet_phi, _packet_v, _packet_moments_correct, _as_output, _verify_hashes,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_compliant_"


def _hash(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _factored(value, log_scale):
    return dict(positive_scale_log=log_scale,
                normalized_axial_coefficients=list(value.coefficients),
                exact_representation="exp(positive_scale_log)*normalized_axial_coefficients")


def _tree(value):
    if isinstance(value, dict):
        return {key: _tree(row) for key, row in value.items()}
    return list(value.coefficients)



def switch_control_source_bridge():
    """Bind the original unequal angular/axial first-switch weights."""
    path=HERE/(PREFIX+"microswitch_mixed_C4.py")
    tree=ast.parse(path.read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="switch_controls")
    branch=next(n for n in ast.walk(fn) if isinstance(n,ast.If)
                and ast.unparse(n.test)=="branch == 'first'")
    if len(branch.body)!=1 or not isinstance(branch.body[0],ast.Assign):
        raise ValueError("Original angular first-switch branch changed")
    if ast.unparse(branch.body[0].targets[0])!="a" or ast.unparse(branch.body[0].value)!="tinyD":
        raise ValueError("First angular shear must keep unweighted tinyD")
    if len(branch.orelse)!=1 or not isinstance(branch.orelse[0],ast.Assign):
        raise ValueError("Original second angular switch branch changed")
    second="[sum((tinyD[j]*(math.comb(k,j)*complement[k-j]) for j in range(k+1)),tinyD[0]*0)+sig[k]*c.mpf('.8') for k in range(4)]"
    if ast.dump(branch.orelse[0].value)!=ast.dump(ast.parse(second,mode="eval").body):
        raise ValueError("Original second angular shear composition changed")
    expected={
        "sig":"[cutoff[k]*math.factorial(k) for k in range(5)]",
        "complement":"[1-sig[0]]+[-row for row in sig[1:]]",
        "tinyD":"[scale_width(row,k+1) for k,row in enumerate(Dbar_y)]",
        "logF":"[scale_width(row,1)*(-c.mpf('.5')) for row in a]",
    }
    proofs={"first_angular_shear_is_unweighted_tinyD":True}
    for target,expression in expected.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
               and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(v)==wanted for v in found)!=1:
            raise ValueError("Original angular switch assignment changed: "+target)
        proofs[target+"_source_binding"]=True
    post=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="postpower")
    for target,expression in {"length":"c.ln(c.mpf(110)/100)-2*self.h","zeta":"length*fraction","theta":"c.exp(-zeta)","phi":"phi2*theta**c.mpf('.4')"}.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        found=[n.value for n in ast.walk(post) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(v)==wanted for v in found)!=1:
            raise ValueError("Original complete postpower changed: "+target)
        proofs["postpower_"+target+"_source_binding"]=True
    proofs["second_angular_control_source_binding"]=True
    return dict(identities=proofs,input_hashes={path.name:_hash(path.name)},
                first_angular_weight="1",first_axial_weight="1-sigma(s)",
                paper_source="Lei-Ren Part I v2: original 100-to-110 switch prescription")


class FirstSwitchLeading:
    def __init__(self):
        self.bridge = CompliantInnerBridgeProfiles()
        self.ctx = c = self.bridge.ctx
        self.comparison = json.loads((HERE / COMPARISON_NAME).read_bytes())
        self.comparison_check = json.loads((HERE / COMPARISON_CHECK_NAME).read_bytes())
        _verify_hashes(self.comparison)
        _verify_hashes(self.comparison_check)
        if not self.comparison_check.get("all_passed"):
            raise ValueError("Accepted comparison-history checker required")
        if (self.comparison["actual_five_defect_family_sha256"] != self.bridge.family
                or self.comparison["implicit_source_sha256"] != self.bridge.source
                or self.comparison["datum_enclosure_sha256"] != self.bridge.core.datum.datum_sha):
            raise ValueError("First switch must inherit the exact current source family")
        self.pulse_weight = c.mpf("0.5")
        self.pressure_log = 2 * self.bridge.core.logP
        self.swirl_log = -2 * self.bridge.core.logC - 2 * self.bridge.core.Lambda * self.bridge.core.Gbar
        self.hashes = dict(self.comparison.get("input_hashes", {}))
        self.hashes.update({
            COMPARISON_NAME: _hash(COMPARISON_NAME),
            COMPARISON_CHECK_NAME: _hash(COMPARISON_CHECK_NAME),
            PREFIX + "flat_pulse_derivatives.py": _hash(PREFIX + "flat_pulse_derivatives.py"),
            PREFIX + "microswitch_mixed_C4.py": _hash(PREFIX + "microswitch_mixed_C4.py"),
            Path(__file__).name: _hash(Path(__file__).name),
        })

    def packet(self, label):
        c = self.ctx
        macro = self.comparison["comparison_point_packets"][label]["macro"]
        at_r100 = _packet_for_coordinate(c, macro, mp.mpf(1))
        z = read_interval(c, at_r100["Z"])
        phi = _packet_phi(c, at_r100)
        value = _packet_v(c, at_r100)
        moments = _packet_moments_correct(c, at_r100)
        inputs = self.bridge.inputs(z)
        current = direction(c, z, self.bridge.delta, phi, value, moments,
                            inputs["p0"], inputs["F0_ratios"], inputs["F0_squared_ratios"])

        # First chart keeps a=hb*Dbar: the angular equation has NO cutoff.
        # d_s log F=-hb^2*Dbar/2; only the axial shear is switched off.
        jf = current["D_over_R"] * (-c.mpf(100) / 2)
        # d_s V = -h_b^2 (1-sigma(s))*G + O(h_b^3), with R=100 at this order.
        hydro = current["drive_hydro"] * (-c.mpf(100) * self.pulse_weight)
        pressure = current["drive_pressure"] * (-c.mpf(100) * self.pulse_weight)
        swirl = current["drive_swirl"] * (-c.mpf(10000) * self.pulse_weight)
        return dict(
            Z=z,
            source_radius="R=100*exp(hb*s)",
            phase_bounds=[0, 1],
            exact_pulse_weight=self.pulse_weight,
            exact_angular_first_chart_weight=c.mpf(1),
            endpoint_moments_R100=_tree(moments),
            D_over_R_at_R100=_as_output(current["D_over_R"]),
            J_logF_hb2=_as_output(jf),
            J_V_hydro=_factored(hydro, c.mpf(0)),
            J_V_pressure=_factored(pressure, self.pressure_log),
            J_V_swirl=_factored(swirl, self.swirl_log),
            J_V=dict(exact_representation="sum(exp(scale_log_i)*normalized_i(Z))",
                     terms=dict(hydro=_factored(hydro, c.mpf(0)),
                                pressure=_factored(pressure, self.pressure_log),
                                swirl=_factored(swirl, self.swirl_log)),
                     amplitude_scales_materialized=False),
            formal_endpoint_changes=dict(
                log_F_100exp_hb_over_F100="hb^2*J_logF_hb2 + O(hb^3)",
                V_100exp_hb_minus_V100="hb^2*J_V + O(hb^3)",
                actual_over_comparison_quotient_at_leading_order="1",
                pressure_scale_log=self.pressure_log,
                swirl_scale_log=self.swirl_log),
            source_equation="V_s=-hb^2*(1-sigma(s))*(phi_actual/barphi)*G(R,Z)",
            angular_source_equation="logF_s=-hb^2*Dbar(R,Z)/2; no sigma factor",
            source_fields_at_leading_order="comparison moments at exact R=100 endpoint; actual/barphi=1",
            higher_hb_orders_resolved=False,
            full_first_switch_resolved=False,
        )

    def report(self):
        return dict(
            actual_five_defect_family_sha256=self.bridge.family,
            implicit_source_sha256=self.bridge.source,
            datum_enclosure_sha256=self.bridge.core.datum.datum_sha,
            comparison_source_namespace=self.comparison["original_comparison_source"]["source_namespace"],
            exact_first_switch="first micro chart: R=100*exp(hb*s), 0<=s<=1",
            pulse_weight_identity="integral_0^1(1-sigma(s))ds=1/2",
            source_control_bridge=switch_control_source_bridge(),
            packets={label: self.packet(label) for label in (".5", "0", "exact_shared_root")},
            leading_hb2_first_switch_only=True,
            full_first_switch_resolved=False,
            actual_signed_bridge_completed=False,
            input_hashes=self.hashes,
        )


def run():
    with mp.workdps(300):
        result = FirstSwitchLeading().report()
    result["input_hashes"][Path(__file__).name] = _hash(Path(__file__).name)
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(_encode(result), indent=2) + "\n", encoding="utf8")
    print("First-switch leading h_b^2 signed packet generated", flush=True)
    return result


if __name__ == "__main__":
    run()
