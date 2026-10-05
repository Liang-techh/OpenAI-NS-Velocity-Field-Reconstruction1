"""Current twenty-owner Cartesian/time source assembly including native pulses.

The fourteen-chart physical acceptance is retained without regeneration.
Six current native pulse charts and their reciprocal-boundary overlap use
the unchanged original physical operators and the same current source graph.
These are signed/logarithmic source bounds, not production point fields.
"""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_global_physical_assembly as original
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    CompliantGlobalPhysicalAssembly as BASE, UZ, UT, UR, P)
from lei_ren_part1_paper_compliant_current_downstream_physical_assembly import (
    CurrentDownstreamPhysicalAssembly, OPEN)
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher, PULSE_CHARTS, CHARTS, UNIFORM,
    SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum, CompliantOuterParameters
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

RECEIPT=PREFIX+"current_pulse_physical_assembly_check.json"
PREVIOUS=PREFIX+"current_downstream_physical_assembly.json"


def original_pulse_operator_bindings():
    """Bind pulse branches separately from unrelated post-pulse assignments."""
    tree=ast.parse((HERE/(PREFIX+"global_physical_assembly.py")).read_text(encoding="utf8"))
    result={}
    for method,wanted in {
            "radius":{"offset":("v","v/self.params.mu","13/self.params.mu+v")},
            "normalized_sources":{
                "offset":("c.mpf(value)","c.mpf(value)/self.params.mu","13/self.params.mu+c.mpf(value)"),
                "logs[4]":("self.logP-(c.mpf('.5')+self.params.mu)*offset",),
                "amplitudes[UT]":("{4:mp.mpf(1)}",),
                "amplitudes[UZ]":("{4:mp.mpf(1)}",)}}.items():
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        branches=[n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test)=="chart in PULSE"]
        if len(branches)!=1:raise ValueError("Unique original pulse physical branch required")
        branch=branches[0];rows={}
        for target,expressions in wanted.items():
            actual=[n.value for n in ast.walk(branch) if isinstance(n,ast.Assign)
                and any(ast.unparse(v)==target for v in n.targets)]
            for expression in expressions:
                expected=ast.dump(ast.parse(expression,mode="eval").body)
                if sum(ast.dump(v)==expected for v in actual)!=1:
                    raise ValueError("Original pulse operator source changed: "+method+"."+target)
            rows[target]=list(expressions)
        if method=="radius":
            expected=ast.dump(ast.parse("self.logRp+offset",mode="eval").body)
            if sum(isinstance(n.value,ast.Tuple) and ast.dump(n.value.elts[0])==expected
                    for n in ast.walk(branch) if isinstance(n,ast.Return))!=1:
                raise ValueError("Original pulse radius return changed")
        else:
            expected=ast.dump(ast.parse("amplitudes[UR].update({4:mp.mpf(1)})",mode="eval").body)
            if sum(ast.dump(n)==expected for n in ast.walk(branch) if isinstance(n,ast.Call))!=1:
                raise ValueError("Original full radial pulse unit changed")
        result[method]=rows
    result["original_component_units"]=assignment_source_bindings("global_physical_assembly","normalized_sources",{
        "logs":"[c.mpf(0),self.logP,c.mpf(0),c.mpf(0),c.mpf(0),logR,c.ln(2)]",
        "amplitudes":"{UZ:{},UT:{},UR:{5:mp.mpf('.5'),6:mp.mpf('-.5')},P:{1:mp.mpf(2)}}"})
    return dict(original_pulse_branch_assignments=result,
        full_theta_axial_radial_units_and_pressure_Pstar_squared_retained=True,passed=True)


def pressure_publication_bindings():
    tree=ast.parse((HERE/(PREFIX+"pulse_radial_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="pressure_moment")
    expected=ast.dump(ast.parse("p+p0",mode="eval").body)
    entries=[v.value for n in ast.walk(fn) if isinstance(n,ast.Return) and isinstance(n.value,ast.Call)
        for v in n.value.keywords if v.arg=="P_over_Pstar_squared"]
    if len(entries)!=1 or ast.dump(entries[0])!=expected:
        raise ValueError("Native absolute pressure publication changed")
    tree=ast.parse((HERE/(PREFIX+"pulse_mixed_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="transport_mixed")
    for target,expression in (("primitive","[point['pressure']['P_over_Pstar_squared']]"),("physical","pressure")):
        entries=[value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(v)==target for v in n.targets) and isinstance(n.value,ast.Dict)
            for key,value in zip(n.value.keys,n.value.values)
            if isinstance(key,ast.Constant) and key.value=="P_over_Pstar_squared"]
        if len(entries)!=1 or ast.dump(entries[0])!=ast.dump(ast.parse(expression,mode="eval").body):
            raise ValueError("Native pressure mixed grid publication changed")
    assignment_source_bindings("pulse_mixed_C4","transport_mixed",{"pressure":"primitive['P_over_Pstar_squared']"})
    return dict(native_absolute_p_plus_P0_publication=True,
        same_pressure_primitive_to_physical_grid=True,passed=True)


def pulse_factor_rebase_proof():
    """Combine the common +/-offset before enclosing gigantic source logs."""
    lp,rp,t,mu,l2=s.symbols("logP logRp offset mu log2",real=True)
    old={UT:lp-(s.Rational(1,2)+mu)*t,UZ:lp-(s.Rational(1,2)+mu)*t,
        UR:lp-(s.Rational(1,2)+mu)*t+(rp+t-l2)/2,P:2*lp}
    new={UT:lp-mu*t+rp/2-(rp+t)/2,UZ:lp-mu*t+rp/2-(rp+t)/2,
        UR:lp-mu*t+rp/2-l2/2,P:2*lp}
    count=0
    for n in range(5):
        for label in old:
            if s.expand((old[label]-n*(rp+t)/2)-(new[label]-n*(rp+t)/2))!=0:
                raise ArithmeticError("Pulse fixed physical unit rebase changed its source")
            count+=1
    assignments=assignment_source_bindings("current_pulse_physical_assembly","normalized_sources",{
        "(grids, logs, amplitudes)":"BASE.normalized_sources(self,chart,Z,value,packet,provider,logR)",
        "logs[2]":"self.logRp/2","logs[4]":"self.logP-scaled",
        "amplitudes[UT]":"{2:mp.mpf(1),4:mp.mpf(1),5:mp.mpf('-.5')}",
        "amplitudes[UZ]":"{2:mp.mpf(1),4:mp.mpf(1),5:mp.mpf('-.5')}",
        "amplitudes[UR]":"{2:mp.mpf(1),4:mp.mpf(1),6:mp.mpf('-.5')}"})
    coordinates=[assignment_source_bindings("current_pulse_physical_assembly","normalized_sources",{
        "scaled":expression}) for expression in ("self.params.mu*c.mpf(value)","c.mpf(value)","13+self.params.mu*c.mpf(value)")]
    return dict(exact_component_factor_identities_including_radial_orders0_through4=count,
        original_full_source_derivative_grids_retained=True,
        giant_offset_cancellation_removed_before_interval_bounds=True,
        no_new_derivative_of_fixed_basepoint_units=True,
        actual_rebase_assignments=assignments,actual_scaled_coordinate_bindings=coordinates,passed=True)


def native_parameter_source_bridge(owner):
    """Typed constructor/AST source identity precedes copy diagnostics."""
    pulse=owner.pulse;future=pulse.high.base.future;initial=future.angular.initial
    datums=(pulse.pulse.initial.datum,initial.datum,owner.pre.datum,owner.core.datum)
    if any(type(v) is not CompliantPressureDatum or type(v.parameters) is not CompliantOuterParameters
            or v.parameters.Md!="40" or v.parameters.precision!=160 for v in datums):
        raise ValueError("One pinned compliant parameter constructor source required")
    for key in ("definition","source_sha","datum_sha","input_hashes"):
        if any(getattr(v,key)!=getattr(datums[0],key) for v in datums):
            raise ValueError("Native/current parameter defining source differs: "+key)
    if any(type(v).normalized_jets is not CompliantPressureDatum.normalized_jets for v in datums):
        raise ValueError("Native/current analytic pressure callable differs")
    constructor=class_assignment("pressure_source","CompliantPressureDatum","__init__",
        "self.parameters","CompliantOuterParameters(Md,precision)")
    name="lei_ren_part1_paper_logarithmic_outer_parameters.py"
    tree=ast.parse((HERE/name).read_text(encoding="utf8"))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=="LogarithmicOuterParameters")
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=="__init__")
    wanted={"self.logPstar":"c.exp(self.md)+11","self.log_mu":"c.ln(c.mpf('.001'))-4*self.logPstar",
        "choice":"-4*self.logPstar-30","limit":"-200*c.ln(10)",
        "self.log_delta":"choice","self.mu":"c.exp(self.log_mu)","self.delta":"c.exp(self.log_delta)"}
    for target,expression in wanted.items():
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(v)==target for v in n.targets)]
        expected=ast.dump(ast.parse(expression,mode="eval").body)
        if sum(ast.dump(v)==expected for v in values)!=1:
            raise ValueError("Actual mu/delta parameter formula changed: "+target)
    test=ast.dump(ast.parse("endpoints(choice)[1]<endpoints(limit)[0]",mode="eval").body)
    branch=[n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.dump(n.test)==test]
    if len(branch)!=1 or not any(isinstance(n,ast.Assign) and ast.unparse(n.value)=="choice"
            and any(ast.unparse(v)=="self.log_delta" for v in n.targets) for n in branch[0].body):
        raise ValueError("Selected exact delta min branch changed")
    for datum in datums:
        p=datum.parameters;c=p.ctx
        if endpoints(-4*p.logPstar-30)[1]>=endpoints(-200*c.ln(10))[0]:
            raise ValueError("Md40 selected delta branch not justified")
    paths={}
    for stem,values in {
            "future_swirl_energy":{"self.params":"self.repair.params","self.mu":"self.repair.mu","self.delta":"self.repair.delta"},
            "outer_angular_repair":{"self.params":"self.heat.params","self.mu":"self.params.mu","self.delta":"self.heat.delta"},
            "exact_heat_component":{"self.params":"self.angular.params","self.delta":"self.angular.delta"},
            "outer_angular_candidate":{"self.params":"self.buffer.params","self.mu":"self.params.mu","self.delta":"self.initial.delta"},
            "outer_buffer":{"self.params":"self.initial.params"},
            "outer_initial":{"self.params":"self.datum.parameters","self.delta":"self.repair.delta"},
            "five_moment_repair":{"self.delta":"c.mpf(endpoints(self.datum.parameters.delta))"},
            "pre_pulse_mixed_C4":{"self.params":"self.initial.params","self.delta":"self.initial.delta"}}.items():
        paths[stem]=assignment_source_bindings(stem,"__init__",values)
    graph=dict(native_future_parameter_objects= future.params is future.repair.params is future.heat.params is future.angular.params is initial.params is initial.datum.parameters,
        native_future_mu_object=future.mu is future.params.mu is pulse.high.base.mu,
        native_future_delta_object=future.delta is future.repair.delta is future.heat.delta is future.angular.delta is initial.delta,
        native_initial_parameter_object=pulse.pulse.initial.params is pulse.pulse.initial.datum.parameters,
        native_delta_endpoint_copy=endpoints(future.delta)==endpoints(initial.datum.parameters.delta),
        native_pulse_mu_endpoint_copy=endpoints(pulse.mu)==endpoints(future.mu),
        native_pulse_delta_endpoint_copy=endpoints(pulse.delta)==endpoints(future.delta),
        current_pre_parameter_object=owner.params is owner.pre.params is owner.pre.datum.parameters,
        direct_current_mapper_delta_copy=endpoints(owner.delta)==endpoints(owner.params.delta),
        retained_core_delta_copy_compatibility=endpoints(owner.core.delta)==endpoints(owner.delta))
    if not all(graph.values()):raise ValueError("Actual native/current parameter origin or copy differs: "+str(graph))
    return dict(typed_common_parameter_constructor=True,constructor_source_assignment=constructor,
        same_exact_parameter_definition_and_source_hashes=True,actual_mu_delta_defining_assignments=wanted,
        exact_Md40_delta_choice_branch_proved=True,original_source_object_paths=paths,
        native_current_parameter_object_and_copy_graph=graph,
        endpoint_copy_checks_are_diagnostics_after_defining_source_identity=True,
        cross_instance_parameter_objects_not_required_identical=True,passed=True,
        input_hashes={name:sha(name),PREFIX+"pressure_source.py":sha(PREFIX+"pressure_source.py")})


class _GapOverlapView:
    """Private coverage view; forwards the unchanged same-object native owner."""
    def provider(self,chart):
        if chart!="pulse_gap":raise ValueError("Overlap view only owns gap coverage")
        return self.dispatch.provider(chart)

    def evaluate(self,chart,Z,coordinate):
        field=self.provider(chart);c=field.ctx
        if endpoints(c.mpf(coordinate))!=endpoints(c.mpf(["12","12.0001"])):
            raise ValueError("Only the accepted supplemental gap coordinate box is allowed")
        return self.dispatch.gap_overlap(Z)


class CurrentPulsePhysicalAssembly(BASE):
    @source_precision
    def __init__(self,require_checked=True):
        self.dispatch=CurrentNativePulseSourceDispatcher()
        manifest=self.dispatch.manifest();self.source_owners=manifest["ordered_current_chart_registry"]
        self.hashes=dict(self.dispatch.hashes);self.family=self.dispatch.family;self.source=self.dispatch.source
        self.datum_sha=self.dispatch.datum_sha
        self.core=self.dispatch.anchor.patch.core;self.ctx=c=self.core.ctx
        self.pre=self.dispatch.rh_reference;self.params=self.pre.params;self.pulse=self.dispatch.native_pulse
        self.delta=c.mpf(endpoints(self.params.delta))
        self.logP=c.mpf(endpoints(self.params.logPstar));self.logC=self.core.logC
        self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw
        self.logRv=self.logRp+13/self.params.mu
        if (manifest["current_downstream_chart_owner_count"]!=20
                or not manifest["current_native_pulse_source_ownership_certified"]
                or not manifest["current_Rp_external_pulse_join_certified"]
                or self.params is not self.pre.datum.parameters
                or self.core.datum.datum_sha!=self.datum_sha
                or manifest[UNIFORM]):
            raise ValueError("Accepted current twenty-owner graph and limited interface scope required")
        self.parameter_source_bindings=dict(
            current_mapper=assignment_source_bindings("current_pulse_physical_assembly","__init__",{
                "self.delta":"c.mpf(endpoints(self.params.delta))","self.params":"self.pre.params",
                "self.pulse":"self.dispatch.native_pulse",
                "self.logP":"c.mpf(endpoints(self.params.logPstar))","self.logC":"self.core.logC",
                "self.logRref":"c.ln(110)+10*(self.logC+self.logP)",
                "self.logRp":"self.logRref+self.logP+1+self.params.Tw",
                "self.logRv":"self.logRp+13/self.params.mu"}),
            original_pulse_radius_and_units=original_pulse_operator_bindings(),
            native_mu_delta=assignment_source_bindings("pulse_radial_C4","__init__",{
                "self.mu":"box(self.high.base.mu)","self.delta":"box(self.high.base.future.delta)"}),
            native_pulse_mu=assignment_source_bindings("axial_amplitude_selection","__init__",{
                "self.mu":"self.future.mu"}),
            native_pressure=assignment_source_bindings("pulse_radial_C4","pressure_moment",{
                "rows":"self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),5)['normalized_pressure_coefficients']",
                "p0":"IntervalTaylor(c,[c.mpf(endpoints(v)) for v in rows])"}),
            native_pressure_publication=pressure_publication_bindings())
        self.rebase_proof=pulse_factor_rebase_proof()
        self.native_parameter_bridge=native_parameter_source_bridge(self)
        self.hashes.update(self.native_parameter_bridge["input_hashes"])
        self.operator_bindings=dict(radius=getattr(type(self),"radius") is BASE.radius,
            cartesian_source_row=BASE.evaluate.__globals__["cartesian_source_row"] is original.cartesian_source_row,
            fixed_x_time_source_row=BASE.evaluate.__globals__["time_source_row"] is original.time_source_row,
            original_normalized_source_rows_retained=self.rebase_proof["original_full_source_derivative_grids_retained"])
        if not all(self.operator_bindings.values()):raise ValueError("Original physical source operators changed")
        name=PREFIX+"current_downstream_physical_assembly_check.json"
        receipt=accepted(name,self.family,self.source,"current_downstream_cartesian_spatial4_time1_certified")
        previous=json.loads((HERE/PREVIOUS).read_bytes())
        if (receipt["datum_enclosure_sha256"]!=self.datum_sha
                or previous["current_source_owner_registry"]!={k:self.source_owners[k] for k in CHARTS if k not in PULSE_CHARTS}
                or tuple(previous["current_physical_chart_owners"])!=CHARTS[:14]
                or not receipt["independent_current_patch_fixed_unit_fixture"]["passed"]
                or not receipt["unchanged_independent_cartesian_coordinate_fixture"]["passed"]
                or not receipt["unchanged_independent_microscopic_scale_fixture"]["passed"]):
            raise ValueError("Accepted unchanged fourteen-chart source/map evidence differs")
        for path,digest in receipt["input_hashes"].items():
            if path in self.hashes and self.hashes[path]!=digest:
                raise ValueError("Current physical evidence source conflict: "+path)
            self.hashes[path]=digest
        self.hashes[name]=sha(name);self.hashes[PREVIOUS]=sha(PREVIOUS)
        self.previous_physical_receipt=receipt;self.previous_physical_manifest=previous
        self.source_identity_proof=previous["exact_source_and_divergence_identity"]
        if not all(self.current_provider_graph().values()):raise ValueError("Current physical source graph changed")
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,"current_pulse_cartesian_spatial4_time1_certified")
            if receipt["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current pulse physical acceptance datum differs")
            self.hashes.update(receipt["input_hashes"]);self.hashes[RECEIPT]=sha(RECEIPT)
            self.physical_acceptance_loaded=True

    def current_provider_graph(self):
        return dict(**CurrentDownstreamPhysicalAssembly.current_provider_graph(self),
            same_single_current_native_pulse_object=self.pulse is self.dispatch.native_pulse,
            all_native_pulse_provider_objects_identical=all(self.dispatch.provider(k) is self.pulse for k in PULSE_CHARTS),
            checked_current_native_source_ownership=self.dispatch.chain_acceptance_loaded,
            checked_current_Rp_external_join=self.dispatch.Rp_acceptance_loaded,
            same_native_and_current_parameter_source_proof=self.native_parameter_bridge["passed"],
            same_native_and_current_pressure_Cstar_source_graph=all(self.dispatch.Rp_graph.values()),
            original_pulse_ODE_and_coordinate_identities_retained=all(self.dispatch.chain_functional_proof.values()))

    def normalized_sources(self,chart,Z,value,packet,provider,logR):
        grids,logs,amplitudes=BASE.normalized_sources(self,chart,Z,value,packet,provider,logR)
        if chart not in PULSE_CHARTS:return grids,logs,amplitudes
        c=self.ctx
        if chart=="pulse_entrance":scaled=self.params.mu*c.mpf(value)
        elif chart in ("pulse_main","pulse_exit","pulse_gap"):scaled=c.mpf(value)
        else:scaled=13+self.params.mu*c.mpf(value)
        logs=list(logs);logs[2]=self.logRp/2;logs[4]=self.logP-scaled
        amplitudes[UT]={2:mp.mpf(1),4:mp.mpf(1),5:mp.mpf('-.5')}
        amplitudes[UZ]={2:mp.mpf(1),4:mp.mpf(1),5:mp.mpf('-.5')}
        amplitudes[UR]={2:mp.mpf(1),4:mp.mpf(1),6:mp.mpf('-.5')}
        return grids,tuple(logs),amplitudes

    def _annotate(self,packet,owner):
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[owner]["provider"],
            current_source_acceptance_receipt=self.source_owners[owner]["acceptance_receipt"],
            current_pulse_physical_acceptance_receipt=RECEIPT,
            current_pulse_cartesian_spatial4_time1_proved=True,
            current_pulse_cartesian_spatial4_time1_certified=self.physical_acceptance_loaded,
            current_twenty_downstream_physical_source_ownership_certified=self.physical_acceptance_loaded,
            current_source_dispatcher_used=True,current_core_parameter_source_used=True,
            current_Rp_external_pulse_join_certified=self.dispatch.Rp_acceptance_loaded,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))
        if owner in PULSE_CHARTS:
            packet.update(pulse_correlated_radius_velocity_factors_combined_before_bounds=True,
                pulse_source_log_base_meanings=["unused","logPstar","logRp/2","unused",
                    "logPstar-mu*offset","logR","log2"],
                pulse_original_fixed_units_exactly_rebased=True)
        return packet

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau="-1",theta="0",axis=False):
        if chart not in CHARTS or axis:
            raise ValueError("Current physical adapter owns twenty downstream charts; core/axis/bridge/heat admission remains open")
        packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)
        return self._annotate(packet,chart)

    @source_precision
    def evaluate_gap_overlap(self,Z,log_tau="-1",theta="0"):
        # A private view changes only the coverage call. Physical operators
        # still see the original pulse_gap coordinate and the same provider.
        view=object.__new__(type(self));view.__dict__=dict(self.__dict__)
        coverage=_GapOverlapView();coverage.dispatch=self.dispatch;view.dispatch=coverage
        packet=BASE.evaluate(view,"pulse_gap",Z,["12","12.0001"],
            log_tau=log_tau,theta=theta,axis=False)
        packet.update(chart="pulse_gap_overlap",original_chart_for_physical_operators="pulse_gap",
            supplemental_coverage_only=True,same_current_gap_owner_used=True)
        return self._annotate(packet,"pulse_gap")

    def manifest(self):
        old=self.previous_physical_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,
            current_provider_graph_identity=self.current_provider_graph(),
            current_physical_parameter_source_bindings=self.parameter_source_bindings,
            current_native_parameter_defining_source_bridge=self.native_parameter_bridge,
            original_radius_and_cartesian_time_operators_retained=self.operator_bindings,
            unchanged_original_full_evaluate_operator_called=True,
            pulse_fixed_unit_rebase_proof=self.rebase_proof,
            exact_source_and_divergence_identity=self.source_identity_proof,
            current_native_pulse_source_coverage=self.dispatch.coverage(),
            retained_fourteen_chart_physical_evidence=dict(report=PREVIOUS,
                receipt=PREFIX+"current_downstream_physical_assembly_check.json",
                report_sha256=sha(PREVIOUS),owners=list(CHARTS[:14]),
                total_source_contributions=old["total_current_spatial_and_time_source_contributions_checked"],
                same_registry_sources_and_datum=True),
            reused_independent_coordinate_fixture=old["unchanged_independent_cartesian_coordinate_fixture"],
            reused_independent_micro_scale_fixture=old["unchanged_independent_microscopic_scale_fixture"],
            reused_independent_patch_fixed_unit_fixture=old["independent_current_patch_fixed_unit_fixture"],
            current_pulse_cartesian_spatial4_time1_proved=True,
            current_pulse_cartesian_spatial4_time1_certified=self.physical_acceptance_loaded,
            current_twenty_downstream_physical_source_ownership_certified=self.physical_acceptance_loaded,
            current_Rp_external_pulse_join_certified=self.dispatch.Rp_acceptance_loaded,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            current_core_axis_physical_owner_installed=False,all_profile_source_charts_callable=False,
            output_kind="signed physical Cartesian/time source bounds for twenty current downstream owners",
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for chart,domain in self.dispatch.coverage()["chart_boxes"].items():
            packets[chart]=self.evaluate(chart,[-1,1],domain,theta=None)
            print("Current native pulse physically mapped: "+chart,flush=True)
        result.update(whole_current_pulse_physical_maps=packets,
            whole_current_gap_overlap_physical_map=self.evaluate_gap_overlap([-1,1],theta=None),
            input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentPulsePhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current pulse physical source assembly generated: twenty owners, Cartesian spatial4/fixed-x time1",flush=True)
    return result


if __name__=="__main__":
    run()
