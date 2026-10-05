"""Current finite-width R110 histories in original long-reshape mixed4.

Only local provider wiring and source ownership change. Original profile
kernels and physical primitive derivative algorithms remain unchanged.
Current restoration/implicit/global admission is deliberately still open.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_long_reshape_profiles as profile_module
import lei_ren_part1_paper_compliant_long_reshape_mixed_C4 as mixed_module
import lei_ren_part1_paper_compliant_inner_bridge_profiles as bridge_module
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import (
    CompliantActualSwitchMixedC4, assignment_source_bindings, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack

PROFILE=profile_module.CompliantLongReshapeProfiles
MIXED=mixed_module.CompliantLongReshapeMixedC4


def accepted(name, family, source, flag):
    receipt=json.loads((HERE/name).read_bytes())
    _verify_hashes(receipt)
    if (not receipt["all_passed"] or not receipt[flag]
            or (receipt["actual_five_defect_family_sha256"],receipt["implicit_source_sha256"])!=(family,source)):
        raise ValueError("Current same-family acceptance required: "+name)
    return receipt


def bridge_source_definitions(provider):
    """Expose original exact definitions without editing accepted R100 packets."""
    tree=ast.parse(Path(bridge_module.__file__).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="actual")
    definitions=[kw.value for n in ast.walk(fn) if isinstance(n,ast.Call)
                 for kw in n.keywords if kw.arg=="source_integral_definitions"]
    if len(definitions)!=1:
        raise ValueError("Unique original bridge definition schema required")
    entries={}
    for kw in definitions[0].keywords:
        if kw.arg in ("F","V","chi"):
            if not isinstance(kw.value,ast.Constant) or not isinstance(kw.value.value,str):
                raise ValueError("Original defining expression is no longer a literal")
            entries[kw.arg]=kw.value.value
    current=provider.history.upstream.bindings
    for key in ("F","V"):
        binding=current["actual_"+key+"_source_integral"]
        if not binding["verified"] or entries.get(key)!=binding["expression"]:
            raise ValueError("Current finite-width provider changed original "+key+" source")
    if entries.get("chi")!=provider.shared_axial_source["finite_width_source_integral_bindings"]["original_exact_chi"]:
        raise ValueError("Current original microscopic control differs")
    return entries


def profile_source_bindings(provider):
    rows=assignment_source_bindings("long_reshape_profiles","inputs",{
        "inlet":"self.switch.inlet(Z)",
        "original":"jet(inlet['actual_R110_log_shape_axial5_coefficients'])",
        "phi":"jet(inlet['F_actual_over_F0_axial5_coefficients'])",
        "v":"jet(inlet['Uz_actual_axial5_coefficients'])",
        "raw":"inlet['actual_moment_shape_axial5_coefficients']",
        "normalized":"dict(theta=moments[MTH]/(phi*2),theta_z=moments[MTHZ]/(phi*2),pressure=moments[MP]/square(phi),swirl=moments[MZT]['swirl']/square(phi),mean=moments[MZ],axial=moments[MZT]['axial'])",
    })
    normalization=assignment_source_bindings("inner_switch_profiles","inlet",{
        "result":"self.post(Z,110)",
        "phi":"IntervalTaylor(c,result['F_actual_over_F0_axial5_coefficients'])",
        "B":"G*(-self.core.Lambda)+logarithm(phi)+logarithm(1+square(z))+c.ln(c.mpf(220))/2",
    })
    transfer=assignment_source_bindings("long_reshape_profiles","evaluate",{
        "inp":"self.inputs(Z)",
        "mean":"inp['moments']['mean']*theta+inp['v']*(1-theta)",
        "axial":"inp['moments']['axial']*theta+square(inp['v'])*(1-theta)",
    })
    checks={name:getattr(CompliantActualLongReshapeProfiles,name) is getattr(PROFILE,name)
            for name in ("inputs","evaluate","report")}
    checks["current_original_inlet_callable"]=(type(provider.reshape.switch).inlet is profile_module.CompliantInnerSwitchProfiles.inlet)
    checks["physical_mixed_evaluate_callable"]=(CompliantActualLongReshapeMixedC4.original_evaluate is MIXED.evaluate)
    checks["original_local_power_callable"]=(CompliantActualLongReshapeMixedC4.power_before_R110 is MIXED.power_before_R110)
    if not all(checks.values()) or provider.reshape.switch is not provider.history.switch:
        raise ValueError("Original long reshape algorithms/current inlet identity required")
    return dict(unchanged_original_callables=checks,
                actual_R110_input_AST_bindings=rows,actual_R110_normalization_AST_bindings=normalization,
                constant_V_transport_AST_bindings=transfer,
                original_B_C2_and_shear_theorems_consume_same_defining_source=True,
                current_FV_control_bindings=provider.history.upstream.bindings)


class CompliantActualLongReshapeProfiles(PROFILE):
    @source_precision
    def __init__(self, switch_mixed):
        # Retain the original source-general A/T/B/cutoff/reference gates.
        super().__init__()
        prior=(self.family,self.source,self.core.datum.datum_sha)
        self.switch_mixed=switch_mixed
        self.switch=switch_mixed.switch
        self.core=self.switch.core
        self.ctx=c=self.switch.ctx
        self.family=self.switch.family
        self.source=self.switch.source
        if prior!=(self.family,self.source,self.core.datum.datum_sha):
            raise ValueError("Original core/source/datum must be unchanged")
        name=PREFIX+"actual_switch_mixed_C4_check.json"
        self.current_switch_check=accepted(name,self.family,self.source,
                                          "actual_R100_R110_feedback_mixed4_available")
        if self.current_switch_check["datum_enclosure_sha256"]!=self.core.datum.datum_sha:
            raise ValueError("Current switch pressure datum differs")
        self.hashes.update(switch_mixed.hashes)
        self.hashes.update(self.current_switch_check["input_hashes"])
        self.hashes[name]=sha(name)
        self.records["current_actual_switch_mixed4_check"]=self.current_switch_check
        self.cache={};self.proofs=[]
        self.A=read_interval(c,self.core.records["physical_norm_family"]["A_upper"])
        self.T=400*self.A
        join=self.records["reference_join_bounds"]
        if endpoints(read_interval(c,join["B_C2_upper"]))!=endpoints(2*self.A):
            raise ValueError("Same original B C2 bound required")
        radius=read_interval(c,join["B_C2_upper"])*8/self.T
        self.q_source=c.mpf([endpoints(c.mpf(".1")-radius)[0],
                            endpoints(c.mpf(".1")+radius)[1]])
        if (endpoints(self.q_source)[0]<endpoints(c.mpf(".05"))[1]
                or endpoints(self.q_source)[1]>endpoints(c.mpf(".15"))[0]):
            raise ValueError("Original positive kernel rate theorem changed")
        self.logref=10*(self.core.logC+self.core.logP)
        if endpoints(self.logref-self.T)[0]<=8:
            raise ValueError("Original reference ordering changed")
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)


class CompliantActualLongReshapeMixedC4(MIXED):
    original_evaluate=MIXED.evaluate

    @source_precision
    def __init__(self):
        self.switch_mixed=CompliantActualSwitchMixedC4()
        self.history=self.switch_mixed.history
        self.reshape=CompliantActualLongReshapeProfiles(self.switch_mixed)
        self.ctx=c=self.reshape.ctx
        self.family=self.reshape.family;self.source=self.reshape.source
        self.hashes=dict(self.reshape.hashes);self.proofs=[]
        # This receipt certifies unchanged algorithms/finite fixtures only.
        name=PREFIX+"long_reshape_mixed_C4_check.json"
        self.original_check=accepted(name,self.family,self.source,
                    "actual_long_reshape_all_mixed_derivatives_total_order_le4_available")
        self.hashes.update(self.original_check["input_hashes"]);self.hashes[name]=sha(name)
        self.invP2=c.exp(-2*self.reshape.core.logP)
        self.shared_axial_source=self.formal_axial_source()
        self.exact_bridge_source_definitions=bridge_source_definitions(self)
        self.source_bindings=profile_source_bindings(self)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def formal_axial_source(self):
        # Preserve the accepted current signed source namespace, not boxes.
        graph=copy.deepcopy(self.switch_mixed.shared_axial_source)
        if (not graph["caps_are_not_source_integrals"]
                or not graph["finite_width_source_integral_enclosures_available"]
                or graph["formal_integrals_numerically_reconstructed"]):
            raise ValueError("Current exact signed source graph required")
        return graph

    @source_precision
    def evaluate(self,Z,phase):
        packet=self.original_evaluate(Z,phase)
        packet.update(current_actual_R110_feedback_in_long_reshape=True,
                      upstream_actual_switch_mixed4_available=True,
                      Rsh_reference_mixed4_join_certified=False,
                      current_reference_restoration_installed=False,
                      current_actual_moment_patch_installed=False)
        return packet

    @source_precision
    def report(self):
        result=super().report()
        result.update(
            datum_enclosure_sha256=self.reshape.core.datum.datum_sha,
            current_actual_long_reshape_source_bindings=self.source_bindings,
            exact_bridge_source_definitions=self.exact_bridge_source_definitions,
            current_actual_R110_feedback_in_long_reshape=True,
            current_actual_long_reshape_mixed4_available=True,
            upstream_actual_switch_mixed4_available=True,
            current_actual_source_namespace=self.shared_axial_source["shared_source_namespace"],
            original_selected_A_T_logref_retained=True,
            source_enclosures_not_point_histories=True,
            Rsh_reference_mixed4_join_certified=False,
            current_reference_restoration_installed=False,
            current_actual_moment_patch_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False)
        return result


@source_precision
def run():
    result=CompliantActualLongReshapeMixedC4().report()
    Path(__file__).with_suffix(".json").write_text(json.dumps(encode(pack(result)),indent=2)+"\n",encoding="utf8")
    print("Current finite-width R110 histories installed in original long-reshape mixed4",flush=True)
    return result


if __name__=="__main__":
    run()
