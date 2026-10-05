"""Current native pulse charts after the checked current Rp source join.

All six charts use one unchanged native pulse object. The original
functional ODE/coordinate certificate is retained in its narrower scope;
a uniform native mixed-four interface certificate is not asserted here.
"""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_actual_Rp_source_join import (
    CurrentRpPulseSourceDispatcher, CurrentPrePulseSourceDispatcher,
    SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_source_dispatcher import ROUTES
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField
from lei_ren_part1_paper_compliant_pulse_high_jets import CompliantPulseHighJets
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_pulse_interface_certificate import functional_identities
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

PULSE_CHARTS=("pulse_entrance","pulse_main","pulse_exit","pulse_gap",
    "pulse_gap_end","pulse_end")
EXPECTED={
    "pulse_entrance":("entrance","[0,.02/mu]"),
    "pulse_main":("main","[.02,10]"),
    "pulse_exit":("main","[10,11]"),
    "pulse_gap":("gap","[11,12]"),
    "pulse_gap_end":("gap_from_end","[-1/mu,-4]"),
    "pulse_end":("end","[-4,0]")}
CHARTS=CurrentPrePulseSourceDispatcher.CHARTS+PULSE_CHARTS
RECEIPT=PREFIX+"current_native_pulse_source_dispatcher_check.json"
UNIFORM="uniform_pulse_C4_chart_interface_certificate_available"


class CurrentNativePulseSourceDispatcher(CurrentRpPulseSourceDispatcher):
    CHARTS=CHARTS

    def __init__(self,require_checked=True):
        # The already accepted external Rp join is required for every owner.
        super().__init__(require_checked=True)
        self.require_chain_checked=require_checked
        self.chain_loaded=False;self.chain_acceptance_loaded=False

    @source_precision
    def _load_chain(self):
        if self.chain_loaded:return
        super()._load_Rp()
        field=self.native_pulse
        routes={}
        for chart,(method,domain) in EXPECTED.items():
            route=ROUTES[chart]
            if (route[:3]!=("pulse_mixed_C4","CompliantPulseMixedC4",method)
                    or route[4]!=domain or route[6] is not None):
                raise ValueError("Native pulse route/domain changed: "+chart)
            routes[chart]=dict(provider=PREFIX+"pulse_mixed_C4.CompliantPulseMixedC4",
                method=method,coverage_coordinate=route[3],domain=domain,
                acceptance_receipt=RECEIPT)
        methods=dict(
            entrance=field.entrance.__func__ is CompliantAxialPulseField.entrance,
            main=field.main.__func__ is CompliantPulseHighJets.main,
            gap=field.gap.__func__ is CompliantAxialPulseField.gap,
            gap_from_end=field.gap_from_end.__func__ is CompliantAxialPulseField.gap_from_end,
            end=field.end.__func__ is CompliantPulseHighJets.end,
            gap_packet=field._gap.__func__ is CompliantPulseHighJets._gap,
            high_packet=field._high_packet.__func__ is CompliantPulseMixedC4._high_packet,
            external_Rp_receipt_loaded=self.Rp_acceptance_loaded,
            same_native_class=type(field) is CompliantPulseMixedC4)
        if not all(methods.values()):raise ValueError("Unchanged same-object native pulse methods required")
        name=PREFIX+"pulse_interface_certificate.json"
        record=accepted(name,self.family,self.source,
            "exact_functional_main_gap_and_gap_end_identities_certified")
        proof=functional_identities()
        if (proof!=record["source_bound_functional_pulse_identities"]
                or not all(proof.values())
                or not record["functional_axial_derivatives_through5_identified"]
                or not record["exact_uncapped_selected_sources_used"]
                or record["full_pulse_C4_installed"]
                or record["temporal_recursion"]):
            raise ValueError("Original functional interface scope/source differs")
        for path,digest in record["input_hashes"].items():
            if path in self.hashes and self.hashes[path]!=digest:
                raise ValueError("Current pulse functional source conflict: "+path)
            self.hashes[path]=digest
        self.hashes[name]=sha(name);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.chain_routes=routes;self.chain_methods=methods;self.chain_functional_proof=proof
        if self.require_chain_checked:
            record=accepted(RECEIPT,self.family,self.source,
                "current_native_pulse_source_ownership_certified")
            if record["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current pulse chain acceptance datum differs")
            self.hashes.update(record["input_hashes"]);self.hashes[RECEIPT]=sha(RECEIPT)
            self.chain_acceptance_loaded=True
        self.chain_loaded=True

    def provider(self,chart):
        if chart in PULSE_CHARTS:
            self._load_chain();return self.native_pulse
        return super().provider(chart)

    def _packet(self,chart,packet,domain,coordinate):
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            source_provider=PREFIX+"pulse_mixed_C4.CompliantPulseMixedC4",
            acceptance_receipt=RECEIPT,source_coverage_coordinate=coordinate,
            source_coordinate_domain=domain,
            derivative_coordinate="ordinary logR,Z; full radial velocity factors differentiated",
            physical_mixed_grids={"physical_mixed_derivatives_total_order_le4":
                packet["physical_mixed_derivatives_total_order_le4"]},
            source_packet=packet,original_scale_metadata=packet["formal_log_Utheta_over_Pstar"],
            same_single_current_native_pulse_object_used=True,
            current_Rp_external_pulse_join_certified=self.Rp_acceptance_loaded,
            current_native_pulse_source_ownership_certified=self.chain_acceptance_loaded,
            **{UNIFORM:False},
            output_kind="current native pulse source enclosures; no production point selection",
            **dict.fromkeys(SCOPES,False))

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in PULSE_CHARTS:return super().evaluate(chart,Z,coordinate)
        field=self.provider(chart);c=field.ctx;v=c.mpf(coordinate)
        lo,hi=endpoints(v);method,domain=EXPECTED[chart]
        if chart=="pulse_entrance":
            upper=endpoints(c.mpf(".02")/field.mu*field.mu)[1]
            valid=lo>=0 and endpoints(field.mu*v)[1]<=upper
        elif chart=="pulse_gap_end":
            valid=hi<=-4 and endpoints(v+1/field.mu)[0]>=0
        else:
            bounds={"pulse_main":(".02","10"),"pulse_exit":("10","11"),
                "pulse_gap":("11","12"),"pulse_end":("-4","0")}[chart]
            # Endpoint enclosures are allowed; chart coverage is never widened
            # to a different route by relying on the broader native method.
            valid=lo>=endpoints(c.mpf(bounds[0]))[0] and hi<=endpoints(c.mpf(bounds[1]))[1]
        if not valid:raise ValueError("Current native pulse route requires "+chart+" "+domain)
        packet=getattr(field,method)(Z,v)
        route=self.chain_routes[chart]
        return self._packet(chart,packet,domain,route["coverage_coordinate"])

    @source_precision
    def gap_overlap(self,Z):
        """Legal main-coordinate overlap completes the reciprocal endpoint box."""
        field=self.provider("pulse_gap")
        packet=field.gap(Z,["12","12.0001"])
        return self._packet("pulse_gap_overlap",packet,"[12,12.0001]",
            "xi=mu*log(R/Rp); supplemental overlap at the gap coordinate change")

    @source_precision
    def coverage(self):
        self._load_chain();field=self.native_pulse;c=field.ctx
        start=-endpoints(1/field.mu)[0]
        image=13+field.mu*c.mpf(start)
        if endpoints(image)[1]>endpoints(c.mpf("12.0001"))[0]:
            raise ArithmeticError("Original gap coordinate domain not fully covered")
        if endpoints(c.mpf(start)+1/field.mu)[0]<0:
            raise ArithmeticError("Gap-end numerical box not legal for native guard")
        domains=dict(pulse_entrance=[0,endpoints(c.mpf(".02")/field.mu)[1]],
            pulse_main=[".02",10],pulse_exit=[10,11],pulse_gap=[11,12],
            pulse_gap_end=[start,-4],pulse_end=[-4,0])
        return dict(chart_boxes=domains,legal_gap_end_start=start,
            gap_end_start_main_coordinate_enclosure=image,
            supplemental_gap_main_box=["12","12.0001"],
            reciprocal_lower_endpoint_used_only_to_select_legal_coverage=True,
            actual_exact_source_endpoint_not_replaced=True,
            gap_coordinate_source_identity="xi=13+mu*s",
            exact_gap_boundary_covered_by_main_overlap=True,
            original_six_chart_domain_union_fully_enclosed=True)

    def manifest(self):
        result=super().manifest();self._load_chain()
        result["ordered_current_chart_registry"].update(self.chain_routes)
        result.update(current_downstream_chart_owner_count=20,
            current_native_pulse_chart_names=list(PULSE_CHARTS),
            current_native_pulse_method_identity=self.chain_methods,
            retained_original_functional_pulse_identities=self.chain_functional_proof,
            current_native_pulse_coverage=self.coverage(),
            current_native_pulse_source_ownership_proved=True,
            current_native_pulse_source_ownership_certified=self.chain_acceptance_loaded,
            all_current_pulse_charts_installed=True,
            all_current_pulse_charts_callable=True,
            native_pulse_internal_ODE_coordinate_identities_retained=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            quantitative_flat_velocity_interface_bound_ledger_available=False,
            all_profile_source_charts_callable=False,
            full_current_core_to_heat_physical_assembly=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))
        return result


@source_precision
def build():
    field=CurrentNativePulseSourceDispatcher(require_checked=False)
    result=field.manifest();packets={}
    for chart,domain in field.coverage()["chart_boxes"].items():
        packets[chart]=field.evaluate(chart,[-1,1],domain)
        print("Current native pulse source chart: "+chart,flush=True)
    result["whole_current_native_pulse_charts"]=packets
    result["whole_current_gap_coordinate_overlap"]=field.gap_overlap([-1,1])
    result["input_hashes"]=dict(field.hashes)
    return encode(pack(result))


def run():
    result=build()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current native pulse source chain generated: 20 owners, six charts, legal full-domain overlap",flush=True)
    return result


if __name__=="__main__":
    run()
