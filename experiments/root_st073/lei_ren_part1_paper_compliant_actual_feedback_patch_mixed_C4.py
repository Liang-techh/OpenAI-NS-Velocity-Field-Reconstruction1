"""Current admitted implicit patch in unchanged physical mixed4 algorithms."""
import json
from pathlib import Path

import lei_ren_part1_paper_compliant_actual_patch_mixed_C4 as original
from lei_ren_part1_paper_compliant_actual_feedback_moment_patch import (
    CompliantActualFeedbackMomentPatch, accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

BASE=original.CompliantActualPatchMixedC4


def current_mixed_source_bindings():
    callables={name:getattr(CompliantActualFeedbackPatchMixedC4,name) is getattr(BASE,name)
        for name in ("gamma","evaluate")}
    if not all(callables.values()):
        raise ValueError("Original physical patch mixed4 methods changed")
    incoming=assignment_source_bindings("actual_patch_mixed_C4","evaluate",{
        "parent":"self.patch.evaluate(x,Z)",
        "(inverse, data)":"self.patch.coefficients(Z)",
        "initial":"dict(mass=jet('Mz_over_R_axial5'),theta=jet('Mtheta_over_sqrt2_Rm_1p5_Am_axial5'),mixed=jet('Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5'),energy=jet('Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5'),pressure=jet('Mp_over_Am2_axial5'))",
        "p0":"jet('original_P0_axial5')",
        "mixed":"patch_mixed(c,x,Z,self.patch.core.delta,am,data['invAm2'],self.invP2,H,V,g,initial,p0)",
    })
    normalization=assignment_source_bindings("actual_patch_mixed_C4","gamma",{
        "argument":"(c.mpf(x)-c.mpf(center))/self.radius",
        "derivatives":"[rows[k]*(math.factorial(k)/(self.radius**(k+1)*self.N)) for k in range(5)]",
    })
    physical=assignment_source_bindings("actual_patch_mixed_C4","patch_mixed",{
        "P":"[p0+square(am)*pressure[0]]+[square(am)*row for row in pressure[1:]]",
        "Q":"[(2*z*V[k]-(z*mass[k])*(1-delta)-d*derivative(mass[k]))/L for k in range(5)]",
        "Ur":"[sum((Q[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),Q[0]*0) for k in range(5)]",
        "grid":"lambda rows:{'x'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(rows) for n in range(5-k)}",
        "ygrid":"lambda rows:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(log_radial_derivatives(c,x,rows)) for n in range(5-k)}",
    })
    return dict(unchanged_original_mixed_callables=callables,
        current_inherited_histories_and_P0_AST_bindings=incoming,
        original_beta_derivative_normalization_AST_bindings=normalization,
        physical_pressure_radial_prefactor_and_grid_AST_bindings=physical,
        current_Rm_join_uses_same_five_source_transport_identities=True,
        current_Rh_join_uses_same_current_unique_implicit_function=True,
        interval_overlap_not_used_as_functional_join_proof=True)


class CompliantActualFeedbackPatchMixedC4(BASE):
    @source_precision
    def __init__(self):
        self.patch=CompliantActualFeedbackMomentPatch()
        self.ctx=c=self.patch.ctx;self.family=self.patch.family;self.source=self.patch.source
        self.hashes=dict(self.patch.hashes)
        name=PREFIX+"flat_pulse_derivatives_check.json";flat=json.loads((HERE/name).read_bytes())
        _verify_hashes(flat)
        if (not flat["all_passed"] or not flat["original_radial_shape_derivatives_C4_available"]
                or flat["actual_five_defect_family_sha256"]!=self.family):
            raise ValueError("Original compact beta C4 source gate required")
        self.hashes.update(flat["input_hashes"]);self.hashes[name]=sha(name)
        name=PREFIX+"actual_patch_mixed_C4_check.json"
        self.generic_check=accepted(name,self.family,self.source,
            "actual_patch_all_mixed_derivatives_total_order_le4_available")
        if (not self.generic_check["exact_coordinate_and_join_checks"]["passed"]
                or not self.generic_check["independent_physical_fixture"]["passed"]):
            raise ValueError("Unchanged physical/coordinate/support fixture evidence required")
        self.hashes.update(self.generic_check["input_hashes"]);self.hashes[name]=sha(name)
        self.radius=c.mpf(1)/40;self.N=c.mpf(endpoints(self.patch.repair.normalization))
        self.invP2=c.exp(-2*self.patch.core.logP)
        self.bindings=current_mixed_source_bindings()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    @source_precision
    def report(self):
        result=super().report()
        result.update(datum_enclosure_sha256=self.patch.core.datum.datum_sha,
            current_patch_mixed_source_bindings=self.bindings,
            current_actual_moment_patch_installed=True,
            current_actual_patch_mixed4_available=True,
            current_actual_patch_implicit_axial5_recomputed=True,
            current_Rm_and_Rh_functional_mixed4_joins_certified=True,
            current_five_functional_terminal_identities_connected=True,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,
            generic_legacy_fixture_reused_not_current_source_admission=True)
        return result


@source_precision
def run():
    result=CompliantActualFeedbackPatchMixedC4().report()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(encode(pack(result)),indent=2)+"\n").encode("utf8"))
    print("Current original five-bump patch physical x/Z and logR/Z mixed4 generated",flush=True)
    return result


if __name__=="__main__":
    run()
