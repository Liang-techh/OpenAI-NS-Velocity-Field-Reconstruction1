"""One current matched source chain: switches, reshape, restore and patch.

This explicit-chart API returns directed source enclosures with their
original derivative coordinates and formal physical units. It does not
select production points or claim the remaining core-to-heat assembly.
The accepted legacy dispatcher and route table are not modified.
"""
import json
from pathlib import Path

import sympy as s

from lei_ren_part1_paper_compliant_source_dispatcher import CompliantSourceDispatcher, ROUTES
from lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4 import CompliantActualFeedbackPatchMixedC4
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted, sha, HERE, PREFIX
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import return_call_binding
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

SPECS={
    "switch_first":("actual_switch_mixed_C4","CompliantActualSwitchMixedC4","actual_switch_mixed_C4_check"),
    "switch_second":("actual_switch_mixed_C4","CompliantActualSwitchMixedC4","actual_switch_mixed_C4_check"),
    "switch_power":("actual_switch_mixed_C4","CompliantActualSwitchMixedC4","actual_switch_mixed_C4_check"),
    "reshape":("actual_long_reshape_mixed_C4","CompliantActualLongReshapeMixedC4","actual_long_reshape_mixed_C4_check"),
    "inner_reference":("actual_reference_restore_mixed_C4","CompliantActualReferenceRestoreMixedC4","actual_reference_restore_mixed_C4_check"),
    "axial_restore":("actual_reference_restore_mixed_C4","CompliantActualReferenceRestoreMixedC4","actual_reference_restore_mixed_C4_check"),
    "restore_buffer":("actual_reference_restore_mixed_C4","CompliantActualReferenceRestoreMixedC4","actual_reference_restore_mixed_C4_check"),
    "actual_patch":("actual_feedback_patch_mixed_C4","CompliantActualFeedbackPatchMixedC4","actual_feedback_patch_mixed_C4_check"),
}
CURRENT_ROUTES={chart:(stem,cls,*ROUTES[chart][2:5],receipt,ROUTES[chart][6])
    for chart,(stem,cls,receipt) in SPECS.items()}
GATES={
    "actual_switch_mixed_C4_check":("actual_R100_R110_feedback_mixed4_available","actual_R2_moments_velocity_pressure_retained"),
    "actual_long_reshape_mixed_C4_check":("current_actual_long_reshape_mixed4_available","R110_postswitch_power_reshape_local_mixed4_join_certified"),
    "actual_reference_restore_mixed_C4_check":("current_actual_reference_restore_mixed4_available","current_correlated_E_installed","current_original_restore_end_exact_4Z"),
    "actual_feedback_patch_mixed_C4_check":("current_actual_patch_mixed4_available","current_Rm_and_Rh_functional_mixed4_joins_certified","current_Rm_source_transport_and_P0_retained"),
    "actual_Rsh_source_join_check":("current_Rsh_source_functional_join_certified",),
}
GRID_KEYS=("physical_velocity_pressure_phase_Z_mixed4","physical_five_primitive_phase_Z_mixed4",
    "physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4",
    "physical_velocity_pressure_x_Z_mixed4","physical_five_primitive_x_Z_mixed4")
SCOPES=("full_implicit_leading_inputs_recomputed","actual_point_moment_history_recovered",
    "global_completed_tensor_admissibility","full_inner_interfaces_certified",
    "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed","temporal_recursion")


def Rm_functional_identity():
    """Current reference and zero-bump patch share physical functions near Rm."""
    y,z=s.symbols("y z",real=True);R0,Pstar=s.symbols("Rm Pstar",positive=True)
    delta=s.symbols("delta",real=True)
    em,eh,ek,aa,eb,ep,P0=[s.Function(n)(z) for n in ("em","eh","ek","aa","eb","ep","P0")]
    am=s.exp(-s.Rational(3,5))/(1+z*z);x=s.exp(y)
    d=(em,ek,eh,aa/(Pstar*am)**2-eb/2,ep/2)
    V=4*z;H=s.Rational(5,8)+eh*x**(-s.Rational(8,5))
    K=4*z*H+ek*x**(-s.Rational(8,5));m=4*z+em/x
    A=16*z*z+8*z*em/x+aa/x
    b=s.Rational(5,6)+eb*x**(-s.Rational(6,5));p=5+ep*x**(-s.Rational(1,5))
    u=am*x**s.Rational(1,10)
    Q=(2*z*V-(1-delta)*z*m-(1-z*z)*s.diff(m,z))/(1-delta*z*z)
    left=dict(Utheta_over_Pstar=u,Uz=V,Ur_over_sqrt_Rm_over_2=s.sqrt(x)*Q,
        P_over_Pstar2=P0+u*u*p/2,
        Mtheta_over_sqrt2_Rm_1p5_Pstar=am*x**s.Rational(8,5)*H,
        Mtheta_z_over_sqrt2_Rm_1p5_Pstar=am*x**s.Rational(8,5)*K,
        Mz_over_Rm=x*m,Mztheta_over_Rm_Pstar2=x*(A/Pstar**2-u*u*b/2),
        Mp_over_Pstar2=u*u*p/2)
    mass=4*z+d[0]/x
    theta=d[2]+s.Rational(5,8)*x**s.Rational(8,5)
    mixed=4*z*theta+d[1]
    energy=d[3]-s.Rational(5,12)*x**s.Rational(6,5)
    pressure=d[4]+s.Rational(5,2)*x**s.Rational(1,5)
    G=x*mass
    Qp=(2*z*V-(1-delta)*z*mass-(1-z*z)*s.diff(mass,z))/(1-delta*z*z)
    right=dict(Utheta_over_Pstar=am*x**s.Rational(1,10),Uz=V,
        Ur_over_sqrt_Rm_over_2=s.sqrt(x)*Qp,P_over_Pstar2=P0+am*am*pressure,
        Mtheta_over_sqrt2_Rm_1p5_Pstar=am*theta,
        Mtheta_z_over_sqrt2_Rm_1p5_Pstar=am*mixed,Mz_over_Rm=G,
        Mztheta_over_Rm_Pstar2=(8*z*G-16*z*z*x)/Pstar**2+am*am*energy,
        Mp_over_Pstar2=am*am*pressure)
    rows={}
    for name,expression in left.items():
        difference=s.simplify(expression-right[name])
        if difference!=0:raise ArithmeticError("Current external Rm source identity failed: "+name)
        rows[name]=["y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)]
    bindings=dict(
        original_reference_postrestore_terminal_call=return_call_binding(
            "reference_restore_mixed_C4","postrestore",
            "self.packet(Z,self.reference.terminal(Z,offset),[c.mpf(0)]*5,'restore_exit_to_Rm')"),
        original_actual_patch_transported_defects=assignment_source_bindings(
            "actual_moment_patch","actual_data",{
                "transported":"self.reference.defects(Z)['actual_normalized_five_defect_axial5_coefficients']"}),
        original_Rm_defect_coordinate_transform=assignment_source_bindings(
            "reference_restore_profiles","defects",{
                "packet":"self.terminal(Z,offset)",
                "x":"offset+6",
                "invAm2":"square(1+square(z))*c.exp(c.mpf('1.2')-2*self.core.logP)",
                "rows":"[centered['mean_error']*c.exp(x),centered['mixed_error']*c.exp(c.mpf('1.6')*x),centered['angular_error']*c.exp(c.mpf('1.6')*x),(centered['axial_square']*invAm2)*c.exp(x)-centered['swirl_error']*(c.exp(c.mpf('1.2')*x)/2),centered['pressure_error']*(c.exp(c.mpf('.2')*x)/2)]",
            }),
        original_actual_patch_inlet_histories=assignment_source_bindings(
            "actual_moment_patch","evaluate",{
                "mass":"z*4+defects[0]/x",
                "theta":"defects[2]+x**c.mpf('1.6')*c.mpf('.625')",
                "mixed":"(z*theta)*4+defects[1]",
                "energy":"defects[3]-x**c.mpf('1.2')*c.mpf(5)/12",
                "pressure_moment":"defects[4]+x**c.mpf('.2')*c.mpf('2.5')",
            }))
    return dict(exact_physical_function_identities=list(left),
        mixed4_rows_implied_by_exact_open_neighborhood_identity=rows,
        total_implied_physical_mixed4_rows=135,current_source_AST_bindings=bindings,
        reference_side_coordinate="y=log(R/Rm)=offset+6",
        patch_side_coordinate="x=R/Rm=exp(y)",
        original_first_support_lower_edge="49/40 > 1",
        exact_five_defect_source_transform=["em","ek","eh","aa/(Pstar*am)^2-eb/2","ep/2"],
        same_P0_and_source_functions_arbitrary_smooth_Z=True,
        physical_function_identity_precedes_any_enclosure_intersection=True,
        reference_grid_Utheta_and_angular_primitives_convert_by_fixed_exact_am_Z0=True,
        fixed_normalization_multipliers_not_differentiated_in_Z=True,
        interval_overlap_is_not_functional_proof=True,passed=True)


class CurrentMatchedSourceDispatcher(CompliantSourceDispatcher):
    CHARTS=tuple(CURRENT_ROUTES)

    def __init__(self,require_checked=True):
        super().__init__()
        self.anchor=None;self.datum_sha=None;self.require_checked=require_checked
        self.join_proof=None

    def _receipt(self,stem):
        if stem not in self.receipts:
            name=PREFIX+stem+".json"
            record=json.loads((HERE/name).read_bytes());_verify_hashes(record)
            if not record.get("all_passed") or any(not record[k] for k in GATES[stem]):
                raise ValueError("Current matched source owner not accepted: "+stem)
            if (record["actual_five_defect_family_sha256"],record["implicit_source_sha256"],record["datum_enclosure_sha256"])!=(
                self.family,self.source,self.datum_sha):
                raise ValueError("Current matched provider family/source/datum differs")
            for path,digest in record["input_hashes"].items():
                if path in self.hashes and self.hashes[path]!=digest:
                    raise ValueError("Inconsistent current source dependency: "+path)
            self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)
            self.receipts[stem]=record
        return self.receipts[stem]

    @source_precision
    def _load(self):
        if self.anchor is not None:return
        self.anchor=CompliantActualFeedbackPatchMixedC4()
        reference=self.anchor.patch.reference_mixed
        long=reference.long_mixed;switch=long.switch_mixed
        self.family=self.anchor.family;self.source=self.anchor.source
        self.datum_sha=self.anchor.patch.core.datum.datum_sha
        self.providers=dict(actual_switch_mixed_C4=switch,actual_long_reshape_mixed_C4=long,
            actual_reference_restore_mixed_C4=reference,actual_feedback_patch_mixed_C4=self.anchor)
        self.hashes.update(self.anchor.hashes)
        for stem in GATES:self._receipt(stem)
        generic=reference.original_check["structural_checks"]
        transport=self.anchor.patch.generic_moment_check["source_transport_and_implicit_structure"]
        if (not generic["passed"] or generic["exact_actual_Rm_primitive_inlet_identities"]!=5
                or not transport["passed"] or transport["exact_actual_dominant_tail_transport_identities"]!=5):
            raise ValueError("Current common Rm source transport needs unchanged universal proof")
        self.join_proof=Rm_functional_identity()
        if not all(self.provider_graph_identity().values()):
            raise ValueError("One common current matched provider graph required")
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        if self.require_checked:
            name=PREFIX+"current_matched_source_dispatcher_check.json"
            record=accepted(name,self.family,self.source,"all_current_matched_charts_exercised")
            if record["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current dispatcher acceptance datum differs")
            self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)

    def provider_graph_identity(self):
        patch=self.providers["actual_feedback_patch_mixed_C4"]
        reference=self.providers["actual_reference_restore_mixed_C4"]
        long=self.providers["actual_long_reshape_mixed_C4"]
        switch=self.providers["actual_switch_mixed_C4"]
        return dict(switch_is_nested_current_long_switch=switch is long.switch_mixed,
            long_is_nested_current_reference_long=long is reference.long_mixed,
            reference_is_nested_current_patch_reference=reference is patch.patch.reference_mixed,
            patch_reference_is_same_restoration_profile=patch.patch.reference is reference.reference,
            same_current_core_object=patch.patch.core is long.reshape.core is switch.core,
            same_current_history_object=reference.reference.history is long.history is switch.history,
            same_exact_signed_axial_source_graph=reference.shared_axial_source==long.shared_axial_source==switch.shared_axial_source,
            same_canonical_pressure_datum=patch.patch.core.datum.datum_sha==self.datum_sha)

    def provider(self,chart):
        if chart not in CURRENT_ROUTES:
            raise ValueError("Current matched source dispatcher owns exactly eight charts: "+str(chart))
        self._load()
        stem,_,_,_,_,receipt,_=CURRENT_ROUTES[chart]
        self._receipt(receipt)
        return self.providers[stem]

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        provider=self.provider(chart)
        stem,cls,method,description,domain,receipt,extra=CURRENT_ROUTES[chart]
        fn=getattr(provider,method)
        if chart=="actual_patch":packet=fn(coordinate,Z)
        elif extra is not None:packet=fn(Z,coordinate,extra)
        else:packet=fn(Z,coordinate)
        grids={key:packet[key] for key in GRID_KEYS if key in packet}
        if not grids:raise ValueError("Current chart returned no physical mixed4 grid")
        derivative=("original phase s,Z; D_y^k=hb^-k D_s^k in factored source ledgers" if chart.startswith("switch_") and chart!="switch_power"
            else "ordinary x=R/Rm,Z and y=logR,Z" if chart=="actual_patch" else "ordinary y=logR,Z; selector is not derivative coordinate")
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_coverage_coordinate=description,
            source_coordinate_domain=domain,derivative_coordinate=derivative,
            physical_mixed_grids=grids,source_packet=packet,
            original_scale_metadata=packet.get("exact_formal_prefactors",dict(rule="Original grid component units and same positive source radius/amplitude; no numerical common rescaling")),
            source_provider=PREFIX+stem+"."+cls,acceptance_receipt=PREFIX+receipt+".json",
            current_matched_source_history_used=True,same_current_provider_graph=True,
            output_kind="directed source-field enclosures with original formal scales; no production point selection",
            **dict.fromkeys(SCOPES,False))

    def manifest(self):
        self._load()
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,
            ordered_current_chart_registry={name:dict(provider=PREFIX+spec[0]+"."+spec[1],
                method=spec[2],coverage_coordinate=spec[3],domain=spec[4],acceptance_receipt=PREFIX+spec[5]+".json")
                for name,spec in CURRENT_ROUTES.items()},
            current_route_specific_acceptance_gates=GATES,provider_graph_identity=self.provider_graph_identity(),
            current_external_Rm_functional_join=self.join_proof,
            current_R110_functional_join_receipt=PREFIX+"actual_long_reshape_mixed_C4_check.json",
            current_Rsh_functional_join_receipt=PREFIX+"actual_Rsh_source_join_check.json",
            current_Rz_and_restore_exit_join_receipt=PREFIX+"actual_reference_restore_mixed_C4_check.json",
            current_R110_Rsh_Rz_restore_exit_Rm_functional_mixed4_joins_composed=True,
            current_Rh_outgoing_unique_implicit_reference_closure_available=True,
            current_Rh_external_neighbor_join_certified=False,
            all_current_matched_charts_callable=True,legacy_route_table_unchanged=True,
            all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
            all_absolute_radii_automatically_dispatched=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))


@source_precision
def run():
    field=CurrentMatchedSourceDispatcher(require_checked=False)
    result=field.manifest();packets={}
    for chart in CURRENT_ROUTES:
        c=field.provider(chart).ctx
        domain=([1,field.anchor.ctx.exp(1)] if chart=="actual_patch" else [-7,-6] if chart=="restore_buffer"
            else [1,2] if chart=="switch_second" else [0,1])
        if chart=="actual_patch":domain=c.mpf([1,endpoints(c.exp(1))[1]])
        packet=field.evaluate(chart,[-1,1],domain)
        packets[chart]=packet
        print("Current matched source chart: "+chart,flush=True)
    result.update(whole_current_chart_evaluations=packets,all_current_matched_charts_exercised=True,
        input_hashes=dict(field.hashes))
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(encode(pack(result)),indent=2)+"\n").encode("utf8"))
    return result


if __name__=="__main__":
    run()
