"""Install actual finite-width histories into the original switch mixed4 math.

The original evaluate/postpower algorithms are inherited unchanged. A current
source-functional bridge/switch provider supplies every actual parent history.
Physical derivatives preserve hb/amplitude factors until their final sums.
"""
import ast
import copy
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_actual_bridge_switch import (
    CompliantActualBridgeSwitch, source_precision, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_bridge_integrals import bridge_control_bindings
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import (
    CompliantMicroswitchMixedC4, IntervalTaylor, scaled_positive_source)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_source_dispatcher import CompliantSourceDispatcher
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack


def assignment_source_bindings(module,method,expected):
    tree=ast.parse((HERE/(PREFIX+module+".py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    result={}
    for target,expression in expected.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
              and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(row)==wanted for row in rows)!=1:
            raise ValueError("Current source assignment changed: "+module+"."+method+"."+target)
        result[target]=True
    return result


def current_actual_field_bindings():
    return assignment_source_bindings("actual_bridge_integrals","packet",{
        "(angular, axial)":"self.contributions(p,chart,q)",
        "(ell, delta_phi, delta_V)":"self.field_changes(p,angular,axial)",
        "exp_ell":"ell.exp()",
        "exp_minus_one":"IntervalTaylor(c,[c.expm1(ell[0])]+list(exp_ell.coefficients[1:]))",
        "point_delta_phi":"phi0*exp_minus_one",
        "phi":"phi0+point_delta_phi",
        "V":"V0+delta_V",
    })


def actual_R100_trace_bindings():
    inputs=assignment_source_bindings("inner_switch_profiles","inputs",{
        "incoming":"self.bridge.actual(Z,self.bridge.r/100)",
        "phi":"jet(incoming['F_actual_over_F0_axial5_coefficients'])",
        "v":"jet(incoming['Uz_actual_axial5_coefficients'])",
        "stored":"incoming['actual_moment_shape_axial5_coefficients']",
        "source":"self.bridge.inputs(Z)",
    })
    phase=assignment_source_bindings("inner_switch_profiles","phase",{
        "start":"slo==shi==0",
        "(phi, v, moments)":"(inp['phi'],inp['v'],inp['moments'])",
        "theta":"c.mpf(1)",
        "ell":"IntervalTaylor.constant(c,0,5)",
    })
    tree=ast.parse((HERE/(PREFIX+"inner_switch_profiles.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="packet")
    expected={
        "pressure_axis_axial5_coefficients":"list(source['p0'].truncate(5).coefficients)",
        "pressure_increment_true_axial5_divided_by_R_F0_squared":"list(dress(moments[MP],source['F0_squared_ratios']).coefficients)",
        "actual_moment_shape_axial5_coefficients":"coefficient_lists(moments)",
        "F_actual_true_axial5_divided_by_F0":"list(dress(phi,source['F0_ratios']).coefficients)",
    }
    for keyword,expression in expected.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        values=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call)
                for kw in node.keywords if kw.arg==keyword]
        if len(values)!=1 or ast.dump(values[0])!=wanted:
            raise ValueError("Same actual R100 pressure/moment output changed: "+keyword)
    return dict(original_inputs_AST_bindings=inputs,exact_phase0_AST_bindings=phase,
                same_moment_pressure_amplitude_output_AST_bindings=list(expected),
                current_actual_field_AST_bindings=current_actual_field_bindings())


def actual_axial_graph(history,original,hashes):
    """Same defining integrals, with current finite-width source ownership.

    Graph references are exact source recipes, not the values of any box.
    The admitted original bridge and switch controls bind those integrands.
    Both V and centered E reuse identical references at every axial order.
    """
    graph=copy.deepcopy(original)
    if (not graph["caps_are_not_source_integrals"]
            or graph["formal_signed_integrals"]["I_bridge"]["sign"]!=-1
            or graph["formal_signed_integrals"]["I_first_switch"]["sign"]!=-1):
        raise ValueError("Original signed axial source graph required")
    bindings=bridge_control_bindings()
    if not bindings["actual_V_source_integral"]["verified"]:
        raise ValueError("Original bridge V integral must be source-bound")
    if graph["formal_signed_integrals"]["I_bridge"]["exact_original_V_source"]!=bindings[
            "actual_V_source_integral"]["expression"]:
        raise ValueError("Original bridge integrand differs from current finite-width source")
    first=history.switch.phase([-1,1],1)
    if graph["formal_signed_integrals"]["I_first_switch"]["exact_original_V_source"]!=first["exact_Uz_source"]:
        raise ValueError("Original switch integrand differs from actual current source")
    new_bindings={name:hashes[name] for name in (
        PREFIX+"actual_bridge_integrals.py",PREFIX+"actual_bridge_integrals_check.json",
        PREFIX+"actual_bridge_switch.py",PREFIX+"actual_bridge_switch_check.json")}
    namespace=hashlib.sha256(json.dumps(dict(
        family=history.bridge.family,source=history.bridge.source,
        datum=history.bridge.core.datum.datum_sha,
        integrands=graph["formal_signed_integrals"],current_source_bindings=new_bindings),
        sort_keys=True).encode("utf8")).hexdigest()
    old=graph["shared_source_namespace"]
    def rebind(value):
        if isinstance(value,dict):
            return {key:(namespace if key=="shared_source" and row==old else rebind(row))
                    for key,row in value.items()}
        if isinstance(value,list):
            return [rebind(row) for row in value]
        return value
    graph=rebind(graph)
    graph["shared_source_namespace"]=namespace
    graph["source_bindings"].update(new_bindings)
    graph["V110_enclosure_binding"]="current finite-width actual bridge R100 -> original first switch -> constant raw V through R110"
    graph["E_enclosure_binding"]="the same current raw V110 source minus4Z; pressure primitive C is separate"
    graph["finite_width_source_integral_bindings"]=dict(
        bridge=dict(provider="CompliantActualBridgeIntegrals",method="packet",
                    source_call="packet(Z,1,'macro')",signed_terms="signed_axial_integral_terms",
                    own_moments="actual_own_six_moments_axial5",
                    exact_source_width="hb=epsilon_b=cstar*K^-100"),
        switch=dict(provider="CompliantActualBridgeSwitch",method="packet",
                    signed_terms="signed_actual_switch_source_integrals.actual_axial_switch_increment_terms",
                    actual_histories="actual_R100_inlet through actual_R110_inlet"),
        known_comparison_direction_preserved=True,
        source_controls_AST_bindings=bindings,
        current_actual_field_AST_bindings=current_actual_field_bindings(),
        consumed_current_actual_provider_bindings=history.upstream.bindings,
        original_exact_chi="1-(1-hb)*sigma(y/hb); hb=cstar*K^-100")
    graph["finite_width_source_integral_enclosures_available"]=True
    graph["formal_integrals_numerically_reconstructed"]=False
    return graph


def R100_source_join():
    """Original late bridge and flat first-switch controls in common y units."""
    import lei_ren_part1_paper_compliant_bridge_mixed_C4 as bridge_module
    import lei_ren_part1_paper_compliant_microswitch_mixed_C4 as switch_module
    if bridge_module.phase_physical is not switch_module.phase_physical:
        raise ValueError("Both traces must use the same physical primitive operator")
    path=HERE/(PREFIX+"bridge_mixed_C4.py")
    tree=ast.parse(path.read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="evaluate")
    assignments={
        "barphi":"[algebra.lift(comparison['phi'])]+[algebra.lift(0)]*3",
        "barV":"[algebra.lift(comparison['v'])]+[algebra.lift(0)]*3",
        "barlog":"[algebra.lift(0)]*3",
    }
    for target,expression in assignments.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
              and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(row)==wanted for row in rows)!=1:
            raise ValueError("Original late frozen comparison changed: "+target)
    h=s.symbols("hb",positive=True)
    D=s.symbols("D0:4");g=s.symbols("g0:4");q=s.symbols("quotient")
    bridge_logF=[-h*row/2 for row in D]
    switch_logF=[-h**(k+2)*row/2 for k,row in enumerate(D)]
    B=bridge_module.exponential_derivatives(bridge_logF)
    S=switch_module.exponential_derivatives(switch_logF)
    count=0
    def zero(value):
        nonlocal count
        if s.expand(value)!=0:
            raise ArithmeticError("R100 original control/coordinate join changed")
        count+=1
    for k in range(4):
        zero(switch_logF[k]-h**(k+1)*bridge_logF[k])
        By=bridge_logF[k]+(s.Rational(1,2) if k==0 else 0)
        Ss=switch_logF[k]+(h/2 if k==0 else 0)
        zero(Ss-h**(k+1)*By)
        Vy=-h*q*sum(s.binomial(k,j)*B[j]*g[k-j] for j in range(k+1))
        Vs=-q*sum(s.binomial(k,j)*S[j]*h**(k-j+2)*g[k-j] for j in range(k+1))
        zero(Vs-h**(k+1)*Vy)
    for k in range(5):
        zero(S[k]-h**k*B[k])
    return dict(original_frozen_comparison_AST_assignments=3,
                symbolic_original_control_coordinate_identities=count,
                exact_common_physical_primitive_operator=True,
                source_geometry="R=100; bridge chi=hb; first sigma(0) and all positive endpoint jets vanish",
                positive_derivative_units="D_s^k=hb^k*D_y^k, k=1..4",
                actual_trace_data="same upstream.packet(Z,1,'macro') phi/raw V/H/M/K/A/B/C/P0",
                actual_R100_trace_AST_bindings=actual_R100_trace_bindings(),
                comparison_direction="same finite-width frozen comparison endpoint and its own histories",
                no_interval_overlap_used_as_identity=True)


def mixed_source_bindings(provider):
    callables={name:getattr(CompliantActualSwitchMixedC4,name) is getattr(CompliantMicroswitchMixedC4,name)
               for name in ("evaluate","postpower","width")}
    if not all(callables.values()):
        raise ValueError("Original switch mixed4 algorithms changed")
    tree=ast.parse((HERE/(PREFIX+"microswitch_mixed_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="evaluate")
    expected={
        "parent":"self.switch.phase(Z,phase)",
        "inp":"self.switch.bridge.inputs(Z)",
        "comparison":"self.switch.bridge.comparison(c.mpf(Z),self.switch.bridge.r/R)",
        "phi":"jet(parent['F_actual_over_F0_axial5_coefficients'])",
        "v":"jet(parent['Uz_actual_axial5_coefficients'])",
        "raw":"parent['actual_moment_shape_axial5_coefficients']",
        "quotient":"phi/comparison['phi'].truncate(5)",
        "directions":"comparison_radial_directions(c,Z,self.core.delta,comparison['phi'],comparison['v'],comparison['moments'],inp['p0'],inp['F0_ratios'],inp['F0_squared_ratios'])",
        "Dbar":"[rate_rows(directions['D_over_R'],c.mpf(1),k)*R for k in range(4)]",
        "controls":"switch_controls(c,branch,sigma_jets(c,phase-start),Dbar,drive_y,algebra.lift(quotient),algebra.width)",
    }
    proofs={}
    for target,expression in expected.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
              and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(row)==wanted for row in rows)!=1:
            raise ValueError("Original actual mixed4 input changed: "+target)
        proofs[target]=True
    return dict(original_mixed4_callable_identities=callables,
                original_mixed4_actual_input_AST_bindings=proofs,
                actual_switch_object_is_current_history=provider.switch is provider.history.switch,
                canonical_pressure_and_ratios_from_current_actual_bridge=True,
                current_known_comparison_axial6_to_direction_axial5=True,
                source_width_not_a_parameter_value=True)


class CompliantActualSwitchMixedC4(CompliantMicroswitchMixedC4):
    @source_precision
    def __init__(self):
        # Consume the existing math gates, then install current actual history.
        super().__init__()
        original_name=PREFIX+"microswitch_mixed_C4_check.json"
        original_check=json.loads((HERE/original_name).read_bytes())
        _verify_hashes(original_check)
        if (not original_check["all_passed"]
                or not original_check["phase1_R2_and_R110_functional_mixed4_joins_certified"]
                or (original_check["actual_five_defect_family_sha256"],original_check["implicit_source_sha256"])!=(
                    self.family,self.source)):
            raise ValueError("Admitted unchanged original mixed4 math and fixtures required")
        self.hashes.update(original_check["input_hashes"])
        self.hashes[original_name]=sha(original_name)
        self.history=CompliantActualBridgeSwitch()
        name=PREFIX+"actual_bridge_switch_check.json"
        check=json.loads((HERE/name).read_bytes())
        _verify_hashes(check)
        if (not check["all_passed"] or not check["R100_R110_actual_feedback_composed"]
                or (check["actual_five_defect_family_sha256"],check["implicit_source_sha256"],
                    check["datum_enclosure_sha256"])!=(
                    self.history.bridge.family,self.history.bridge.source,self.history.bridge.core.datum.datum_sha)):
            raise ValueError("Current same-source actual finite-width history required")
        self.switch=self.history.switch
        self.core=self.switch.core
        self.ctx=c=self.switch.ctx
        self.family=self.switch.family
        self.source=self.switch.source
        self.hashes.update(self.history.hashes)
        self.hashes.update(check["input_hashes"])
        self.hashes[name]=sha(name)
        self.proofs=[]
        self.logh=self.switch.logh
        self.invP2=c.exp(-2*self.core.logP)
        self.logF0=c.mpf([endpoints(-self.core.logC-self.core.Lambda*self.core.Gbar)[0],
                          endpoints(-self.core.logC)[1]])
        self.h=scaled_positive_source(c,self.logh,IntervalTaylor.constant(c,1,5),self.proofs)[0]
        self.distance=c.ln(100/self.switch.bridge.r)
        self.continuation_distance=c.ln(c.mpf("4.1")/4)
        if (endpoints(self.distance)[0]<=endpoints(2*self.switch.cap)[1]
                or endpoints(self.continuation_distance)[0]<=endpoints(2*self.switch.cap)[1]):
            raise ValueError("Known comparison smoothing must end before R100")
        self.shared_axial_source=actual_axial_graph(self.history,self.shared_axial_source,self.hashes)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.source_bindings=mixed_source_bindings(self)
        self.R100_join=R100_source_join()

    @source_precision
    def report(self):
        result=super().report()
        result.update(
            datum_enclosure_sha256=self.core.datum.datum_sha,
            current_actual_mixed4_source_bindings=self.source_bindings,
            actual_R100_functional_mixed4_join=self.R100_join,
            actual_R100_functional_mixed4_join_certified=True,
            actual_finite_width_bridge_histories_installed_in_switch_mixed4=True,
            actual_R100_R110_feedback_mixed4_available=True,
            switch_mixed4_source_namespace=self.shared_axial_source["shared_source_namespace"],
            current_actual_histories_and_canonical_P0_preserved=True,
            existing_mixed4_algorithms_reused_unchanged=True,
            actual_point_moment_history_recovered=False,
            downstream_reshape_actual_feedback_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            global_completed_tensor_admissibility=False)
        return result


class ActualSwitchSourceDispatcher(CompliantSourceDispatcher):
    """Existing chart evaluation with current receipt ownership for switches."""
    CHARTS=("switch_first","switch_second","switch_power")

    def provider(self,chart):
        if chart not in self.CHARTS:
            raise ValueError("This current-history dispatcher owns the three switch charts only")
        if "actual_switch_mixed_C4" not in self.providers:
            name=PREFIX+"actual_switch_mixed_C4_check.json"
            check=json.loads((HERE/name).read_bytes())
            _verify_hashes(check)
            if not check["all_passed"] or not check["actual_R100_R110_feedback_mixed4_available"]:
                raise ValueError("Current actual switch mixed4 acceptance required")
            self.family=check["actual_five_defect_family_sha256"]
            self.source=check["implicit_source_sha256"]
            self.hashes.update(check["input_hashes"])
            self.hashes[name]=sha(name)
            self.providers["actual_switch_mixed_C4"]=CompliantActualSwitchMixedC4()
        return self.providers["actual_switch_mixed_C4"]

    def evaluate(self,chart,Z,coordinate):
        result=super().evaluate(chart,Z,coordinate)
        result["acceptance_receipt"]=PREFIX+"actual_switch_mixed_C4_check.json"
        result["current_finite_width_actual_history_used"]=True
        return result


def run():
    with mp.workdps(400):
        result=CompliantActualSwitchMixedC4().report()
    Path(__file__).with_suffix(".json").write_text(json.dumps(encode(pack(result)),indent=2)+"\n",encoding="utf8")
    print("Current actual finite-width R100 histories installed in original switch mixed4",flush=True)
    return result


if __name__=="__main__":
    run()
