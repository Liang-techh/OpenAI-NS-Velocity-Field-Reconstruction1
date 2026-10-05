"""Same-source finite-width actual bridge histories through original R100..110.

This local adapter reuses the admitted switch methods unchanged. Comparison
histories still prescribe Dbar/Ebar; actual own histories recover Q/pressure.
No existing mixed4/implicit provider is replaced by this axial companion.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_actual_bridge_integrals import (
    CompliantActualBridgeIntegrals, named_moments, value_enclosure)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles, direction, coefficient_lists, dress,
    MTH, MZ, MTHZ, MZT, MP)
from lei_ren_part1_paper_compliant_inner_switch_profiles import CompliantInnerSwitchProfiles
from lei_ren_part1_paper_compliant_switch_signed_integrals import (
    CompliantSwitchSignedIntegrals, source_precision, switch_control_source_bridge)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX="lei_ren_part1_paper_compliant_"


def sha(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def recontext(c,value):
    """Outward numeric transport; no defining function or parameter changes."""
    if isinstance(value,IntervalTaylor):
        return IntervalTaylor(c,value.coefficients)
    if isinstance(value,dict):
        return {key:recontext(c,row) for key,row in value.items()}
    if isinstance(value,list):
        return [recontext(c,row) for row in value]
    return value


class CoreContext:
    """Same canonical core; only Taylor context objects are reconciled."""
    def __init__(self,core,c):
        self.original=core
        self.ctx=c

    def __getattr__(self,name):
        return getattr(self.original,name)

    def axis_inputs(self,Z):
        return recontext(self.ctx,self.original.axis_inputs(Z))


class ActualR100BridgeAdapter(CompliantInnerBridgeProfiles):
    def __init__(self,upstream):
        self.upstream=upstream
        self.ctx=c=upstream.ctx
        native=upstream.bridge
        self.core=CoreContext(upstream.core,c)
        self.family=self.core.family
        self.source=self.core.source
        self.delta=self.core.delta
        self.r=c.mpf(native.r)
        self.logh=c.mpf(upstream.logh)
        self.cap=upstream.cap
        self.records=native.records
        self.hashes=dict(upstream.hashes)
        self.cap_proofs=[]
        self.cache={}
        self.comparison_cache={}
        if (native.family,native.source,native.core.datum.datum_sha)!=(
                self.family,self.source,self.core.datum.datum_sha):
            raise ValueError("Same original actual/comparison/core family required")

    def inputs(self,Z):
        return self.upstream.prepare(Z)["inputs"]

    def actual(self,Z,theta):
        c=self.ctx
        if c.mpf(theta)._mpi_!=(self.r/100)._mpi_:
            raise ValueError("This functional adapter supports the named R100 interface only")
        incoming=self.upstream.packet(Z,1,"macro")
        if not incoming["coordinate_R100_endpoint"]:
            raise ValueError("Exact original macro endpoint required")
        phi=IntervalTaylor(c,incoming["actual_phi_axial5"])
        moments=named_moments({name:IntervalTaylor(c,row)
                              for name,row in incoming["actual_own_six_moments_axial5"].items()})
        return dict(
            F_actual_over_F0_axial5_coefficients=list(phi.coefficients),
            F_actual_true_axial5_divided_by_F0=list(dress(phi,self.inputs(Z)["F0_ratios"]).coefficients),
            Uz_actual_axial5_coefficients=incoming["actual_raw_V_axial5"],
            actual_moment_shape_axial5_coefficients=coefficient_lists(moments),
            actual_Q_axial4_coefficients=incoming["actual_radial_Q_axial4"],
            pressure_axis_axial5_coefficients=incoming["pressure_axis_axial5"],
            pressure_increment_true_axial5_divided_by_R_F0_squared=incoming[
                "actual_pressure_increment_axial5_divided_by_R_F0_squared"],
            R100_actual_integral_and_own_feedback_source=incoming,
            exact_R100_functional_source="upstream.packet(Z,1,'macro')",
            comparison_moments_substituted=False)

    def comparison(self,Z,theta):
        """Continue the known frozen comparison's actual finite-width history."""
        c=self.ctx
        key=c.mpf(Z)._mpi_
        if key not in self.comparison_cache:
            macro=self.upstream.comparison.macro(Z,1)
            if not macro["exact_endpoint_at_R100"]:
                raise ValueError("Known comparison must start at the same exact R100")
            self.comparison_cache[key]=(
                value_enclosure(macro["phi"]),value_enclosure(macro["V"]),
                {name:value_enclosure(row) for name,row in macro["moments"].items()})
        phi,V,incoming=self.comparison_cache[key]
        theta=c.mpf(theta)
        if theta._mpi_==(self.r/100)._mpi_:
            ratio=c.mpf(1)
        else:
            lower=endpoints(self.r/110)[0]
            upper=endpoints(self.r/100)[1]
            if endpoints(theta)[0]<lower or endpoints(theta)[1]>upper:
                raise ValueError("Known comparison continuation only on R100..R110")
            ratio=intersection(c,theta*(100/self.r),
                               c.mpf([endpoints(c.mpf(100)/110)[0],1]))
        targets=dict(H=phi,M=V,K=phi*V,A=V*V,B=phi*phi/2,C=phi*phi)
        own={name:targets[name]+(row-targets[name])*ratio**(
                2 if name in ("H","K","B") else 1) for name,row in incoming.items()}
        moments=named_moments(own)
        inp=self.inputs(Z)
        rows=direction(c,c.mpf(Z),self.delta,phi,V,moments,
                       inp["p0"],inp["F0_ratios"],inp["F0_squared_ratios"])
        if phi.order!=6 or V.order!=6 or any(row.order!=5 for row in rows.values()):
            raise ValueError("Known source requires axial6 before directional differentiation")
        return dict(phi=phi,v=V,moments=moments,direction=rows,
                    comparison_R100_history_retained=True,
                    exact_comparison_theta_source="100/R",
                    known_comparison_direction_preserved=True)


class ActualInnerSwitch(CompliantInnerSwitchProfiles):
    def __init__(self,bridge):
        self.bridge=bridge
        self.core=bridge.core
        self.ctx=c=bridge.ctx
        self.family=bridge.family
        self.source=bridge.source
        self.logh=bridge.logh
        self.cap=bridge.cap
        self.hashes=dict(bridge.hashes)
        ledger=bridge.records["K1_ledger"]
        self.a_upper=read_interval(c,ledger["decreasing_K_smallness_bounds"]["initial_switch_a"])
        if endpoints(self.a_upper)[1]>=endpoints(c.mpf(4)/5)[0]:
            raise ValueError("Original microscopic shear bound required")
        if endpoints(self.logh)[1]>=endpoints(c.ln(c.ln(c.mpf(110)/100)/2))[0]:
            raise ValueError("Original R2 must precede R110")
        self.cache={}


class ActualSignedSwitch(CompliantSwitchSignedIntegrals):
    @source_precision
    def __init__(self,bridge):
        # Preserve all original signed-switch gates; use the new named inlet.
        super().__init__()
        self.bridge=bridge
        self.ctx=bridge.ctx
        self.family=bridge.family
        self.source=bridge.source
        self.hashes.update(bridge.hashes)


def source_bindings(provider):
    """Callable identities plus exact endpoint geometry, rather than overlap."""
    h,Y,q=s.symbols("hb Y q",real=True)
    if s.expand((2*h+q*(Y-2*h)).subs(q,1)-Y)!=0:
        raise ArithmeticError("Original R100 geometry changed")
    checks={name:getattr(ActualInnerSwitch,name) is getattr(CompliantInnerSwitchProfiles,name)
            for name in ("inputs","packet","phase","post","inlet")}
    checks["signed_switch_packet_callable"]=(
        ActualSignedSwitch.packet is CompliantSwitchSignedIntegrals.packet)
    checks["original_width_product_callable"]=(
        ActualR100BridgeAdapter.width_product is CompliantInnerBridgeProfiles.width_product)
    if not all(checks.values()):
        raise ValueError("Original admitted switch callables changed")
    tree=ast.parse(Path(__file__).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="actual")
    wanted=ast.dump(ast.parse('self.upstream.packet(Z,1,"macro")',mode="eval").body)
    calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call)]
    if sum(ast.dump(n)==wanted for n in calls)!=1:
        raise ValueError("Named original actual R100 source call changed")
    comp=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="comparison")
    continuation={}
    expected={
        "targets":"dict(H=phi,M=V,K=phi*V,A=V*V,B=phi*phi/2,C=phi*phi)",
        "own":"{name:targets[name]+(row-targets[name])*ratio**(2 if name in ('H','K','B') else 1) for name,row in incoming.items()}",
        "ratio":"intersection(c,theta*(100/self.r),c.mpf([endpoints(c.mpf(100)/110)[0],1]))",
        "lower":"endpoints(self.r/110)[0]",
        "upper":"endpoints(self.r/100)[1]",
    }
    for target,expression in expected.items():
        wanted=ast.dump(ast.parse(expression,mode="eval").body)
        matches=[n.value for n in ast.walk(comp) if isinstance(n,ast.Assign)
                 and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(n)==wanted for n in matches)!=1:
            raise ValueError("Known comparison continuation changed: "+target)
        continuation[target]=True
    wanted=ast.dump(ast.parse("self.upstream.comparison.macro(Z,1)",mode="eval").body)
    if sum(ast.dump(n)==wanted for n in ast.walk(comp) if isinstance(n,ast.Call))!=1:
        raise ValueError("Known comparison R100 endpoint call changed")
    guard=ast.dump(ast.parse('not macro["exact_endpoint_at_R100"]',mode="eval").body)
    if sum(ast.dump(n.test)==guard for n in ast.walk(comp) if isinstance(n,ast.If))!=1:
        raise ValueError("Known comparison endpoint geometry guard changed")
    continuation["same_exact_R100_macro_call"]=True
    continuation["exact_endpoint_at_R100_consumed"]=True
    ledger=provider.bridge.records["K1_ledger"]
    if not ledger["h_b_equals_epsilon_b_by_definition"] or not ledger["positive_width_not_materialized"]:
        raise ValueError("Exact original positive width definition required")
    return dict(
        original_callable_identities=checks,
        exact_macro_q1_radius_identity="2hb+(log(100/Ra)-2hb)=log(100/Ra); R=100",
        actual_R100_source_call_AST_bound=True,
        comparison_continuation_AST_bindings=continuation,
        six_moment_map=dict(H=MTH,M=MZ,K=MTHZ,A=MZT+"['axial']",B=MZT+"['swirl']",C=MP),
        P0_and_amplitude_ratios_source="same upstream.prepare(Z)['inputs'] from original bridge.inputs(Z)",
        same_width_source="hb=epsilon_b=cstar*K^-100",
        original_switch_controls=switch_control_source_bridge(),
        known_comparison_direction_preserved=True,
        actual_histories_never_used_to_redefine_prescribed_shear=True)


class CompliantActualBridgeSwitch:
    @source_precision
    def __init__(self):
        self.upstream=CompliantActualBridgeIntegrals()
        name=PREFIX+"actual_bridge_integrals_check.json"
        check=json.loads((HERE/name).read_bytes())
        _verify_hashes(check)
        if (not check["all_passed"]
                or not check["full_finite_width_signed_bridge_integral_enclosures_available"]
                or not check["actual_own_six_moment_feedback_enclosures_available"]
                or (check["actual_five_defect_family_sha256"],check["implicit_source_sha256"],
                    check["datum_enclosure_sha256"])!=(
                    self.upstream.core.family,self.upstream.core.source,self.upstream.core.datum.datum_sha)):
            raise ValueError("Accepted same-family finite-width actual bridge and own feedback required")
        self.bridge=ActualR100BridgeAdapter(self.upstream)
        self.ctx=self.bridge.ctx
        self.switch=ActualInnerSwitch(self.bridge)
        self.signed=ActualSignedSwitch(self.bridge)
        self.hashes=dict(self.upstream.hashes)
        self.hashes.update(check["input_hashes"])
        self.hashes[name]=sha(name)
        self.hashes.update(self.signed.hashes)
        for module in (
                "inner_switch_profiles","inner_bridge_profiles","switch_signed_integrals",
                "first_switch_leading","comparison_point_integrals"):
            self.hashes[PREFIX+module+".py"]=sha(PREFIX+module+".py")
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.bindings=source_bindings(self)

    @source_precision
    def packet(self,Z):
        c=self.ctx
        start=self.switch.phase(Z,0)
        first=self.switch.phase(Z,1)
        second=self.switch.phase(Z,2)
        inlet=self.switch.inlet(Z)
        signed=self.signed.packet(Z)
        if not all(len(inlet[key])==6 for key in (
                "F_actual_over_F0_axial5_coefficients","Uz_actual_axial5_coefficients",
                "pressure_axis_axial5_coefficients",
                "pressure_increment_true_axial5_divided_by_R_F0_squared")):
            raise ValueError("Complete R110 axial5 output required")
        return dict(Z=c.mpf(Z),actual_R100_inlet=start,actual_first_switch_exit=first,
                    actual_R2_inlet=second,actual_R110_inlet=inlet,
                    signed_actual_switch_source_integrals=signed,
                    actual_finite_width_bridge_feedback_consumed=True,
                    same_source_R100_functional_join_installed=True,
                    original_own_six_moments_preserved_through_R110=True,
                    source_function_enclosures_only=True,
                    actual_point_moment_history_recovered=False,
                    existing_mixed4_provider_replaced=False,
                    full_implicit_leading_inputs_recomputed=False,
                    temporal_recursion=False)

    def report(self):
        packets={name:self.packet(Z) for name,Z in (("whole_Z",[-1,1]),("0",0),(".5",".5"))}
        return dict(
            actual_five_defect_family_sha256=self.bridge.family,
            implicit_source_sha256=self.bridge.source,
            datum_enclosure_sha256=self.bridge.core.datum.datum_sha,
            source_bindings=self.bindings,source_axial_domain=[-1,1],
            actual_R100_R110_packets=packets,
            bridge_product_cap_proofs=self.bridge.cap_proofs,
            signed_switch_product_cap_proofs=self.signed.proofs,
            R100_R110_actual_feedback_composed=True,
            arbitrary_Z_actual_R110_velocity_moment_pressure_axial5_available=True,
            known_comparison_direction_preserved=True,
            actual_point_moment_history_recovered=False,
            existing_mixed4_provider_replaced=False,
            full_implicit_leading_inputs_recomputed=False,
            global_completed_tensor_admissibility=False,
            exact_production_point_parameters_selected=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(400):
        result=CompliantActualBridgeSwitch().report()
    Path(__file__).with_suffix(".json").write_text(json.dumps(_encode(result),indent=2)+"\n",encoding="utf8")
    print("Actual finite-width R100 histories composed through original switches to R110",flush=True)
    return result


if __name__=="__main__":
    run()
