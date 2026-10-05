"""Current downstream source fields in physical Cartesian/time coordinates.

Inject current fourteen-chart ownership into unchanged coordinate operators.
The result is signed source bounds, not chosen nonlinear point fields.
No legacy core or unadmitted pulse/heat owner is silently loaded.
"""
import json
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_global_physical_assembly as original
from lei_ren_part1_paper_compliant_global_physical_assembly_check import source_identities
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import (
    CurrentPrePulseSourceDispatcher, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

BASE=original.CompliantGlobalPhysicalAssembly
CHARTS=CurrentPrePulseSourceDispatcher.CHARTS
OPEN=("full_point_physical_field_evaluation","full_background_NS_validation",
    "physical_energy_integral_certified","independently_bounded_flat_remainder",
    "all_33_current_source_charts_physical_spatial4_time1_mapped")


class CurrentDownstreamPhysicalAssembly(BASE):
    @source_precision
    def __init__(self,require_checked=True):
        self.dispatch=CurrentPrePulseSourceDispatcher()
        manifest=self.dispatch.manifest()
        self.source_owners=manifest["ordered_current_chart_registry"]
        self.hashes=dict(self.dispatch.hashes);self.family=self.dispatch.family;self.source=self.dispatch.source
        self.datum_sha=self.dispatch.datum_sha
        self.core=self.dispatch.anchor.patch.core
        self.ctx=c=self.core.ctx;self.delta=self.core.delta
        self.pre=self.dispatch.rh_reference;self.params=self.pre.params
        if (not manifest["current_Rh_to_Rp_source_chain_certified"]
                or manifest["current_downstream_chart_owner_count"]!=14
                or self.pre.datum.datum_sha!=self.datum_sha
                or self.core.datum.datum_sha!=self.datum_sha
                or self.params is not self.pre.datum.parameters):
            raise ValueError("Current nested source/core/pressure/parameter graph required")
        self.logP=c.mpf(endpoints(self.params.logPstar));self.logC=self.core.logC
        self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw
        self.parameter_source_bindings=dict(
            current_mapper=assignment_source_bindings("current_downstream_physical_assembly","__init__",{
                "self.delta":"self.core.delta","self.params":"self.pre.params",
                "self.logP":"c.mpf(endpoints(self.params.logPstar))","self.logC":"self.core.logC",
                "self.logRref":"c.ln(110)+10*(self.logC+self.logP)",
                "self.logRp":"self.logRref+self.logP+1+self.params.Tw"}),
            core=assignment_source_bindings("core_physical_field","__init__",{
                "self.delta":"self.outer.delta","self.logC":"read(norm,'selected_logCstar')"}),
            outer_initial=assignment_source_bindings("outer_initial","__init__",{
                "self.params":"self.datum.parameters","self.delta":"self.repair.delta",
                "self.logC":"read_interval(c,self.repair.records['compliant_physical_norm_family']['selected_logCstar'])"}),
            pre=assignment_source_bindings("pre_pulse_mixed_C4","__init__",{
                "self.datum":"self.initial.datum","self.delta":"self.initial.delta"}),
            repair_delta=assignment_source_bindings("five_moment_repair","__init__",{
                "self.delta":"c.mpf(endpoints(self.datum.parameters.delta))"}))
        if not all(self.current_provider_graph().values()):
            raise ValueError("Current physical Cstar/parameter source graph differs")
        self.operator_bindings={name:getattr(type(self),name) is getattr(BASE,name)
            for name in ("radius","normalized_sources")}
        if not all(self.operator_bindings.values()):
            raise ValueError("Original coordinate source/unit algorithms changed")
        name=PREFIX+"global_physical_assembly_check.json"
        record=accepted(name,self.family,self.source,"exact_source_divergence_identity_checked")
        for key in ("independent_implicit_physical_fixture","independent_micro_scale_fixture"):
            if not record[key]["passed"]:raise ValueError("Independent unchanged coordinate fixture required")
        if (record["independent_implicit_physical_fixture"]["independent_implicit_root_cartesian_derivatives"]!=140
                or record["independent_implicit_physical_fixture"]["independent_fixed_x_time_derivatives"]!=4
                or record["independent_micro_scale_fixture"]["independent_signed_micro_and_macro_source_scale_rows"]!=120):
            raise ValueError("Full fixed-unit, moving-basis and microscopic fixture scope required")
        for path,digest in record["input_hashes"].items():
            if path in self.hashes and self.hashes[path]!=digest:
                raise ValueError("Current physical operator/source dependency conflict: "+path)
        self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)
        self.original_operator_fixture=record
        self.source_identity_proof=source_identities()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.physical_acceptance_loaded=False
        if require_checked:
            name=PREFIX+"current_downstream_physical_assembly_check.json"
            record=accepted(name,self.family,self.source,
                "current_downstream_cartesian_spatial4_time1_certified")
            if record["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current physical acceptance datum differs")
            self.hashes.update(record["input_hashes"]);self.hashes[name]=sha(name)
            self.physical_acceptance_loaded=True

    def current_provider_graph(self):
        return dict(nested_current_core_object_used=self.core is self.dispatch.anchor.patch.core,
            same_current_outer_reference_used=self.pre is self.dispatch.rh_reference,
            exact_pre_parameter_object_used=self.params is self.pre.params is self.pre.datum.parameters,
            same_Cstar_defining_record=(
                self.core.records["physical_norm_family"]==
                self.pre.initial.repair.records["compliant_physical_norm_family"]),
            original_positive_mu_and_Tw_preserved=(
                endpoints(self.params.mu)[0]>0 and endpoints(self.params.Tw)[0]>0),
            same_current_core_and_pre_datum_definition=self.core.datum.definition==self.pre.datum.definition,
            same_current_family_source_and_datum=(
                self.family,self.source,self.datum_sha)==(
                self.dispatch.family,self.dispatch.source,self.dispatch.datum_sha),
            no_legacy_core_owner_initialized="core_physical_field" not in self.dispatch.providers)

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau="-1",theta="0",axis=False):
        if chart not in CHARTS or axis:
            raise ValueError("Current physical adapter owns fourteen downstream charts, no core/axis or pulse/heat yet")
        packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[chart]["provider"],
            current_source_acceptance_receipt=self.source_owners[chart]["acceptance_receipt"],
            current_downstream_cartesian_spatial4_time1_proved=True,
            current_downstream_cartesian_spatial4_time1_certified=self.physical_acceptance_loaded,
            current_source_dispatcher_used=True,current_core_parameter_source_used=True,
            current_Rp_external_pulse_join_certified=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))
        return packet

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,
            current_provider_graph_identity=self.current_provider_graph(),
            current_physical_parameter_source_bindings=self.parameter_source_bindings,
            unchanged_coordinate_source_methods=self.operator_bindings,
            unchanged_original_full_evaluate_operator_called=True,
            exact_source_and_divergence_identity=self.source_identity_proof,
            reused_independent_coordinate_fixture=self.original_operator_fixture["independent_implicit_physical_fixture"],
            reused_independent_micro_scale_fixture=self.original_operator_fixture["independent_micro_scale_fixture"],
            current_downstream_cartesian_spatial4_time1_proved=True,
            current_downstream_cartesian_spatial4_time1_certified=self.physical_acceptance_loaded,
            current_Rp_external_pulse_join_certified=False,
            current_core_axis_physical_owner_installed=False,
            all_profile_source_charts_callable=False,
            output_kind="signed physical Cartesian/time source bounds for fourteen current downstream charts",
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        c=self.ctx
        domains={chart:[0,1] for chart in CHARTS}
        domains.update(switch_second=[1,2],restore_buffer=[-7,-6],
            actual_patch=[1,endpoints(c.exp(1))[1]],Rh_reference=[-5,0],O2_buffer=[0,11])
        result=self.manifest();packets={}
        for chart in CHARTS:
            packets[chart]=self.evaluate(chart,[-1,1],domains[chart],theta=None)
            print("Current downstream physical source mapped: "+chart,flush=True)
        result.update(whole_current_downstream_physical_maps=packets,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentDownstreamPhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current fourteen-chart physical spatial4/fixed-x time1 source assembly generated",flush=True)
    return result


if __name__=="__main__":
    run()
