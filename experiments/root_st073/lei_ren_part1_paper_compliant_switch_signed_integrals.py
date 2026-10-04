"""Signed full-source R100-to-R110 switch integrals and width coefficients.

The original angular first chart has unit weight; only its axial shear
uses 1-sigma. The frozen comparison's complete six moments determine exact
exponential direction modes. All short-switch integrals are bounded as
signed axial jets, with hb retained as the original positive source.
These are enclosures of source functions, not selected point values.
"""
from functools import wraps
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_first_switch_leading import (
    FirstSwitchLeading, switch_control_source_bridge)
from lei_ren_part1_paper_compliant_macro_signed_integrals import (
    _mode_rows, _as_output, _verify_hashes)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    IntervalTaylor, symmetric, logarithm)
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import scaled_positive_source
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX="lei_ren_part1_paper_compliant_"


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def switch_integral_identities():
    hb,Y,J,D=s.symbols("hb Y JD Dbar100",real=True)
    first=-hb**2*D/2
    second=-hb/s.Integer(5)-hb**2*D/4
    post=-s.Rational(2,5)*(Y-2*hb)
    total=-s.Rational(2,5)*Y+s.Rational(3,5)*hb-s.Rational(3,4)*hb**2*D
    if s.expand(first+second+post-total)!=0:
        raise ArithmeticError("Original angular switch composition changed")
    exact=-hb/s.Integer(5)-hb**2*J/2-s.Rational(2,5)*(Y-2*hb)
    if s.expand(exact-(-s.Rational(2,5)*Y+s.Rational(3,5)*hb-hb**2*J/2))!=0:
        raise ArithmeticError("Exact JD-to-R110 identity changed")
    return dict(first_angular_weight1=True,second_angular_weight_half=True,
                first_axial_weight_half=True,second_and_post_axial_shear_zero=True,
                exact_R2_logF="-hb/5-hb^2*JD/2",
                exact_R110_logF="-2*log(110/100)/5+3*hb/5-hb^2*JD/2",
                JD="int_0^1 Dbar(100exp(hb*s))ds+int_0^1 (1-sigma(s))*Dbar(100exp(hb*(1+s)))ds",
                JD_at_zero_width="3*Dbar100/2",
                logF_hb2="-3*Dbar100/4",
                normalized_F110_ratio_hb2="9/50-3*Dbar100/4",
                width_independent_of_Z=True,
                actual_source_width_not_a_cap=True)


def exponential_range(c,rate,width_cap):
    """Enclosure of exp(rate*hb*s), 0<=s<=1 and 0<hb<=cap."""
    return c.exp(c.mpf([0,endpoints(width_cap)[1]])*rate)


def angular_kernel_ranges(c,width_cap):
    """Exact positive-weight bounds for the two short angular integrals."""
    weights=[]
    for rate in (1,0,-1):
        first=exponential_range(c,rate,width_cap)
        # The actual second-chart weight integrates to1/2, but its radius
        # includes the preceding hb chart. No cap is substituted for hb.
        preceding=exponential_range(c,rate,width_cap)
        second=preceding*first/2
        weights.append(first+second)
    return weights


def weighted_range(rows,weights):
    out=rows[0]*0
    for row,weight in zip(rows,weights):
        out+=row*weight
    return out


def source_term(row,log_scale):
    return dict(positive_source_log=log_scale,
                signed_axial_coefficients=_as_output(row),
                exact_source_representation="exp(source_log)*signed_coefficient(Z)")


def source_precision(fn):
    @wraps(fn)
    def call(*args,**kwargs):
        with mp.workdps(300):return fn(*args,**kwargs)
    return call


class CompliantSwitchSignedIntegrals:
    @source_precision
    def __init__(self):
        self.leading=FirstSwitchLeading()
        self.bridge=self.leading.bridge; self.ctx=c=self.bridge.ctx
        self.family=self.bridge.family; self.source=self.bridge.source
        self.hashes=dict(self.leading.hashes); self.hashes.update(self.bridge.hashes)
        name=PREFIX+"first_switch_leading_check.json"
        check=json.loads((HERE/name).read_bytes())
        _verify_hashes(check)
        if not check["all_passed"] or not check["original_angular_first_weight1_checked"]:
            raise ValueError("Correct original angular first-switch weight required")
        self.hashes.update(check["input_hashes"]); self.hashes[name]=digest(name)
        self.hashes[PREFIX+"macro_signed_integrals.py"]=digest(PREFIX+"macro_signed_integrals.py")
        from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
        terminal=read_interval(c,self.bridge.records["global_exit_certificate"]["terminal_comparison_D_floor"])
        if endpoints(terminal)[0]<=0:
            raise ValueError("Accepted positive comparison Dbar on100..110 required")
        if endpoints(2*self.bridge.cap)[1]>=endpoints(c.ln(c.mpf(110)/100))[0]:
            raise ValueError("Original second switch must precede110")
        self.terminal_Dbar_floor=terminal
        self.proofs=[]
        self.source_control_bridge=switch_control_source_bridge()
        self.hashes.update(self.source_control_bridge["input_hashes"])
        name=PREFIX+"microswitch_mixed_C4_check.json"
        receipt=json.loads((HERE/name).read_bytes()); _verify_hashes(receipt)
        if (not receipt["all_passed"]
            or receipt["actual_five_defect_family_sha256"]!=self.family
            or receipt["implicit_source_sha256"]!=self.source
            or not receipt["complete_postswitch_power_installed"]):
            raise ValueError("Matching full original switch/power source required")
        self.hashes.update(receipt["input_hashes"]); self.hashes[name]=digest(name)

    @source_precision
    def packet(self,Z):
        c=self.ctx; z=c.mpf(Z); hcap=self.bridge.cap
        inputs=self.bridge.inputs(z)
        comparison=self.bridge.comparison(z,self.bridge.r/100)
        phi=comparison["phi"].truncate(5); v=comparison["v"].truncate(5)
        moments={key:({part:row.truncate(5) for part,row in value.items()}
                      if isinstance(value,dict) else value.truncate(5))
                 for key,value in comparison["moments"].items()}
        modes=_mode_rows(self.bridge,z,inputs,phi,v,moments)
        modes={key:[row.truncate(5) for row in values] for key,values in modes.items()}
        incoming=self.bridge.actual(z,self.bridge.r/100)
        actualphi=IntervalTaylor(c,incoming["F_actual_over_F0_axial5_coefficients"])
        actualv=IntervalTaylor(c,incoming["Uz_actual_axial5_coefficients"])
        quotient=actualphi/phi
        angular=weighted_range(modes["D_over_R"],angular_kernel_ranges(c,hcap))*100
        # Every prefix of chart1 has length<=1. Bound its axial jets before
        # the exponential; the admitted same-source theorem gives Dbar>0.
        Drange=weighted_range(modes["D_over_R"],
                             [exponential_range(c,rate,hcap) for rate in (1,0,-1)])*100
        bound=IntervalTaylor(c,[abs(row) for row in Drange.coefficients])
        scaled=scaled_positive_source(c,2*self.bridge.logh,bound,self.proofs)/2
        ell=IntervalTaylor(c,[c.mpf([-endpoints(scaled[0])[1],0])]
                          +[symmetric(c,row) for row in scaled.coefficients[1:]])
        ratio=quotient*ell.exp()
        axial={}
        for part,power,scale in (
            ("hydro",1,c.mpf(0)),
            ("pressure",1,2*self.bridge.core.logP),
            ("swirl",2,-2*self.bridge.core.logC-2*self.bridge.core.Lambda*self.bridge.core.Gbar)):
            # Coefficientwise integration against the positive original
            # weight1-sigma, of total mass1/2, retains signed bounds.
            drive=weighted_range(modes["drive_"+part],
                                 [exponential_range(c,power-j,hcap) for j in range(3)])*(100**power)
            coeff=-(ratio*drive)/2
            axial[part]=source_term(coeff,2*self.bridge.logh+scale)
        constant=lambda value:IntervalTaylor.constant(c,value,5)
        linear=scaled_positive_source(c,self.bridge.logh,constant(c.mpf("0.6")),self.proofs)
        quadratic=scaled_positive_source(c,2*self.bridge.logh,-angular/2,self.proofs)
        log_correction=linear+quadratic
        phi110=actualphi*log_correction.exp()*(c.mpf(100)/110)**c.mpf(".4")
        v110=actualv
        for part,term in axial.items():
            row=IntervalTaylor(c,term["signed_axial_coefficients"])
            v110+=scaled_positive_source(c,term["positive_source_log"],row,self.proofs)
        return dict(Z=z,source_domain="R100 through both microscopic charts and postpower toR110",
                    original_width_source="hb=cstar*K^-100=epsilon_b; independent of Z",
                    original_width_log=self.bridge.logh,
                    comparison_frozen_phi_axial5=_as_output(phi),
                    comparison_frozen_V_axial5=_as_output(v),
                    comparison_R100_own_moments_axial5={key:({part:_as_output(row) for part,row in value.items()} if isinstance(value,dict) else _as_output(value)) for key,value in moments.items()},
                    actual_incoming_phi_axial5=_as_output(actualphi),
                    actual_incoming_V_axial5=_as_output(actualv),
                    comparison_direction_modes={key:[_as_output(row) for row in rows] for key,rows in modes.items()},
                    actual_JD_signed_axial5_enclosure=_as_output(angular),
                    actual_first_chart_logF_prefix_axial5_enclosure=_as_output(ell),
                    actual_incoming_over_frozen_comparison_quotient_axial5=_as_output(quotient),
                    actual_axial_switch_increment_terms=axial,
                    actual_R110_phi_over_F0_axial5_enclosure=_as_output(phi110),
                    actual_R110_V_axial5_enclosure=_as_output(v110),
                    exact_logF110_source_terms=dict(
                        zeroth="-2*log(110/100)/5",
                        linear=source_term(constant(c.mpf(".6")),self.bridge.logh),
                        quadratic=source_term(-angular/2,2*self.bridge.logh)),
                    complete_switch_signed_integral_enclosures_available=True,
                    source_modes_use_complete_comparison_moments=True,
                    actual_incoming_bridge_retained=True,
                    comparison_field_not_substituted_for_actual_field=True,
                    numerical_cap_is_only_error_bound=True,
                    newly_selected_point_values=False,
                    actual_signed_bridge_completed=False,temporal_recursion=False)

    @source_precision
    def leading_packet(self,label):
        first=self.leading.packet(label); c=self.ctx
        D=IntervalTaylor(c,first["D_over_R_at_R100"])
        zero=D*0; firstj=D*(-50); secondj=D*(-25); total=firstj+secondj
        return dict(Z=first["Z"],first_angular_hb2=_as_output(firstj),
                    second_angular_hb1=_as_output(zero-c.mpf(".2")),
                    second_angular_hb2=_as_output(secondj),
                    postpower_angular_hb1=_as_output(zero+c.mpf(".8")),
                    R110_logF_constant=-c.mpf(".4")*c.ln(c.mpf(110)/100),
                    R110_logF_hb1=_as_output(zero+c.mpf(".6")),
                    R110_logF_hb2=_as_output(total),
                    R110_normalized_F_ratio_hb2=_as_output(total+c.mpf(".18")),
                    R110_logUtheta_constant=c.mpf(".1")*c.ln(c.mpf(110)/100),
                    R110_axial_hb2_terms=first["J_V"]["terms"],
                    all_width_orders_through2=True,
                    finite_width_bridge_higher_order_feedback_solved=False)

    def report(self):
        packets={label:self.packet(z) for label,z in (("whole_Z",[-1,1]),("0","0"),(".5",".5"))}
        root=self.leading.comparison["comparison_point_packets"]["exact_shared_root"]["macro"][0]["Z"]
        from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
        packets["exact_shared_root"]=self.packet(read_interval(self.ctx,root))
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    datum_enclosure_sha256=self.bridge.core.datum.datum_sha,
                    exact_switch_composition=switch_integral_identities(),
                    consumed_positive_comparison_Dbar_100_to_110_floor=self.terminal_Dbar_floor,
                    source_control_bridge=self.source_control_bridge,
                    signed_source_integral_packets=packets,
                    leading_width_packets={label:self.leading_packet(label)
                                           for label in ("0",".5","exact_shared_root")},
                    bound_proofs=self.proofs,source_width_cap=self.bridge.cap,
                    whole_Z_full_switch_signed_enclosures_available=True,
                    higher_actual_Ra_to_R100_bridge_orders_solved=False,
                    actual_signed_bridge_completed=False,
                    selected_nonlinear_point_values_recovered=False,
                    temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(300):
        result=CompliantSwitchSignedIntegrals().report()
    result["input_hashes"][Path(__file__).name]=digest(Path(__file__).name)
    Path(__file__).with_suffix(".json").write_text(json.dumps(_encode(result),indent=2)+"\n",encoding="utf8")
    print("Original complete switch signed integrals and R110 width coefficients generated",flush=True)
    return result


if __name__=="__main__":run()
