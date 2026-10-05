"""Current finite-width/restored histories in the original implicit five-bump patch.

Original transport, nonlinear map, common Jacobian, bumps and pressure
formulas are inherited unchanged. Legacy receipts admit only unchanged
algorithms and generic fixtures; this provider has its own acceptance.
"""
import json
from pathlib import Path

import lei_ren_part1_paper_compliant_actual_moment_patch as original
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import (
    CompliantActualReferenceRestoreMixedC4, accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import SharedFiveMomentRepair, pack
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

BASE=original.CompliantActualMomentPatch


def current_patch_source_bindings():
    callables={name:getattr(CompliantActualFeedbackMomentPatch,name) is getattr(BASE,name)
        for name in ("actual_data","coefficients","evaluate")}
    if not all(callables.values()):
        raise ValueError("Original actual patch transport/implicit/field algorithms changed")
    transport=assignment_source_bindings("actual_moment_patch","actual_data",{
        "inp":"self.reference.inputs(Z)",
        "initial":"inp['source_centered_Rsh']",
        "E":"inp['E']",
        "invAm2":"square(1+square(z))*c.exp(c.mpf('1.2')-2*self.core.logP)",
        "exact":"[a+b for a,b in zip(dominant,tails)]",
        "transported":"self.reference.defects(Z)['actual_normalized_five_defect_axial5_coefficients']",
    })
    inverse=assignment_source_bindings("actual_moment_patch","coefficients",{
        "data":"self.actual_data(Z)",
        "inverse":"implicit_axial_jets(self.ctx,self.W,data['actual_defects'],data['invAm2'],self.scales)",
    })
    pressure=assignment_source_bindings("actual_moment_patch","evaluate",{
        "p0":"self.reference.inputs(Z)['original_axis_pressure']",
        "P":"p0+square(am)*pressure_moment",
        "(inverse, data)":"self.coefficients(Z)",
    })
    return dict(unchanged_original_callables=callables,
        current_transport_source_AST_bindings=transport,
        current_implicit_inputs_AST_bindings=inverse,current_pressure_AST_bindings=pressure,
        same_implicit_Jacobian_function_for_orders_2_through_5=True,
        caps_and_legacy_point_coefficients_not_selected=True)


class CompliantActualFeedbackMomentPatch(BASE):
    @source_precision
    def __init__(self,require_checked=True):
        self.reference_mixed=CompliantActualReferenceRestoreMixedC4()
        self.reference=self.reference_mixed.reference;self.core=self.reference.core
        self.ctx=c=self.reference.ctx;self.family=self.reference.family;self.source=self.reference.source
        self.repair=SharedFiveMomentRepair();self.hashes=dict(self.reference_mixed.hashes);self.cache={}
        if (self.repair.admit["actual_five_defect_family_sha256"]!=self.family
                or self.repair.datum.source_sha!=self.source
                or self.repair.admit["implicit_source_sha256"]!=self.source
                or self.repair.datum.datum_sha!=self.core.datum.datum_sha
                or self.repair.admit["datum_enclosure_sha256"]!=self.core.datum.datum_sha):
            raise ValueError("Current patch requires same family/source and both pressure datum hashes")
        for part,flag in (
            ("actual_reference_restore_mixed_C4_check","current_actual_reference_restore_mixed4_available"),
            ("actual_Rsh_source_join_check","current_Rsh_source_functional_join_certified")):
            name=PREFIX+part+".json";record=accepted(name,self.family,self.source,flag)
            if record["datum_enclosure_sha256"]!=self.core.datum.datum_sha:
                raise ValueError("Current upstream patch pressure datum differs")
            if part=="actual_reference_restore_mixed_C4_check" and not (
                record["current_correlated_E_installed"] and record["current_original_restore_end_exact_4Z"]
                and record["current_E_source_C2_theorem_verified_before_intersection"]):
                raise ValueError("Current original E/restoration source class not admitted")
            self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)
        name=PREFIX+"five_moment_repair_check.json";generic=json.loads((HERE/name).read_bytes())
        _verify_hashes(generic)
        if (generic["actual_five_defect_family_sha256"]!=self.family or not all(generic[k] for k in (
            "actual_implicit_functional_five_moment_closure_independently_checked",
            "corrected_partial_moments_and_same_pressure_independently_checked",
            "residual_zero_containment_not_used_as_closure_proof"))):
            raise ValueError("Unchanged original functional five-bump map fixture missing")
        self.hashes.update(generic["input_hashes"]);self.hashes[name]=sha(name)
        name=PREFIX+"actual_moment_patch_check.json"
        self.generic_moment_check=accepted(name,self.family,self.source,
            "actual_patch_coefficient_axial5_enclosures_available")
        if (not self.generic_moment_check["source_transport_and_implicit_structure"]["passed"]
                or not self.generic_moment_check["independent_coefficient_fixture"]["passed"]
                or not self.generic_moment_check["residual_zero_containment_not_used_as_closure_proof"]):
            raise ValueError("Unchanged transport/Jacobian/finite physical coefficient fixture required")
        self.hashes.update(self.generic_moment_check["input_hashes"]);self.hashes[name]=sha(name)
        self.hashes.update(self.repair.hashes)
        convert=lambda v:c.mpf(endpoints(v))
        self.W={name:([[convert(v) for v in row] for row in values] if name=="L"
            else [convert(v) for v in values]) for name,values in self.repair.W.items()}
        self.scales=[convert(v) for v in self.repair.scales]
        self.oldtail=read_interval(c,self.repair.admit["actual_tail_C1_positive_cap"])
        self.oldrho=read_interval(c,self.repair.admit["actual_v1_minus_4Z_minus_j_C1_upper"])
        self.e=read_interval(c,self.repair.admit["complete_actual_functional_e_upper"])
        admit=self.repair.admit
        if (admit["actual_tail_cap_definition"]!="10^-800/Pstar^2; positive and checked against all exact log bounds"
                or "inherited summed C2 boundrho" not in admit["proof"]["C1"]
                or "independent Taylor boxes enlarge this correlated class" not in admit["proof"]["boxes"]
                or "true source functions, not every arbitrary function" not in admit["proof"]["boxes"]):
            raise ValueError("Original normalized correlated-C1 source class semantics changed")
        self.source_class_admission=dict(
            tail_cap_definition=admit["actual_tail_cap_definition"],
            C2_to_correlated_C1_source_proof=admit["proof"]["C1"],
            correlated_functions_not_arbitrary_independent_boxes=admit["proof"]["boxes"],
            current_C2_bound_controls_original_summed_C1_norm=True)
        self.bindings=current_patch_source_bindings()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        if require_checked:
            name=PREFIX+"actual_feedback_moment_patch_check.json"
            record=accepted(name,self.family,self.source,"current_actual_moment_patch_installed")
            if (record["datum_enclosure_sha256"]!=self.core.datum.datum_sha
                    or not record["current_actual_patch_implicit_axial5_recomputed"]
                    or not record["current_five_functional_terminal_identities_connected"]):
                raise ValueError("Current independent patch acceptance required")
            self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)

    @source_precision
    def report(self):
        result=super().report()
        result.update(datum_enclosure_sha256=self.core.datum.datum_sha,
            current_patch_source_bindings=self.bindings,
            original_correlated_source_class_admission=self.source_class_admission,
            current_exact_centered_E_source=self.reference_mixed.shared_axial_source["E_V110_minus_4Z"],
            current_reference_restoration_histories_used=True,
            current_Rsh_source_functional_join_certified=True,
            current_actual_moment_patch_installed=True,
            current_actual_patch_implicit_axial5_recomputed=True,
            current_five_functional_terminal_identities_connected=True,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,
            generic_legacy_fixtures_reused_not_current_source_admission=True)
        return result


@source_precision
def run():
    result=CompliantActualFeedbackMomentPatch(require_checked=False).report()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(encode(pack(result)),indent=2)+"\n").encode("utf8"))
    print("Current five defects and common implicit axial5 Jacobian connected to original patch",flush=True)
    return result


if __name__=="__main__":
    run()
