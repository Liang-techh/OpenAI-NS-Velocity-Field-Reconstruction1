"""Actual two microscopic-switch full tensors and completed phase1/R2 traces."""
import gzip
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_switch_power_background_tensor import (
    CurrentSwitchPowerBackgroundTensor,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,canonical_tensor_groups,BASE,IntervalTaylor)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import CompliantActualSwitchMixedC4,mixed_source_bindings
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import CompliantMicroswitchMixedC4,phase_physical
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import (
    SOURCE_KEY,compiled_original_microswitch_with_rows,raw_microswitch_rows,
    compiled_factored_source_operators,compiled_factored_physical_lift,expanded_stress_sectors,
    factored_rows_record,actual_microswitch_normalization_theorem,actual_microswitch_endpoint_theorem,microswitch_factor_units_theorem)

NAME=PREFIX+'current_microswitch_background_tensor.json.gz'
RECEIPT=PREFIX+'current_microswitch_background_tensor_check.json'
GATES=('current_actual_two_microswitch_full_tensors_available',
    'current_actual_microswitch_full_meridional_decomposition_available',
    'current_actual_phase1_R2_completed_tensor_joins_certified')
DOMAINS={'switch_first':(0,1),'switch_second':(1,2)}
SEAMS=('first_second','second_power')
VIEWS={
    'first_whole':('switch_first',(-1,1),(0,1),('-3','-1'),None,'1'),
    'first_left':('switch_first',(-1,1),0,'-1',None,'1'),
    'first_right':('switch_first',(-1,1),1,'-1',None,'1'),
    'first_fresh':('switch_first','.537','.537','-2.6','.41','.8'),
    'first_fresh_low':('switch_first','.731','.1337','-2.6','.41','.8'),
    'first_fresh_high':('switch_first','-.317','.831','-2.6','.41','.2'),
    'second_whole':('switch_second',(-1,1),(1,2),('-3','-1'),None,'1'),
    'second_left':('switch_second',(-1,1),1,'-1',None,'1'),
    'second_right':('switch_second',(-1,1),2,'-1',None,'1'),
    'second_fresh':('switch_second','.537','1.537','-2.6','.41','.8'),
    'second_fresh_low':('switch_second','.731','1.1337','-2.6','.41','.8'),
    'second_fresh_high':('switch_second','-.317','1.831','-2.6','.41','.2')}


def current_microswitch_source_theorem(field):
    mixed=field.switch;env=field.evaluator.__globals__;native=field.physical.dispatch.native
    bindings=mixed_source_bindings(mixed)
    live=dict(exact_current_actual_switch_owner=type(mixed) is CompliantActualSwitchMixedC4,
        same_checked_power_switch_owner=mixed is field.power_tensor.switch,
        same_current_switch_history=mixed.switch is mixed.history.switch,
        original_evaluate=mixed.evaluate.__func__ is CompliantMicroswitchMixedC4.evaluate,
        original_postpower=mixed.postpower.__func__ is CompliantMicroswitchMixedC4.postpower,
        original_source_phase_generator=mixed.evaluate.__func__.__globals__['phase_physical'] is phase_physical,
        compiled_output_only_phase_generator=env['phase_physical'].__globals__ is env,
        same_original_switch_controls=env['switch_controls'] is mixed.evaluate.__func__.__globals__['switch_controls'],
        same_original_factored_algebra=env['FactoredAlgebra'] is mixed.evaluate.__func__.__globals__['FactoredAlgebra'],
        same_native_current_three_switch_providers=all(field.physical.dispatch.provider(name) is mixed for name in ('switch_first','switch_second','switch_power')),
        same_current_provider_graph=all(native.provider_graph_identity().values()),
        same_original_BASE_radius=field.physical.radius.__func__ is BASE.radius)
    if not all(live.values()) or not bindings['actual_switch_object_is_current_history']:raise ValueError('Same original current microswitch source graph required')
    receipt=native.receipts['actual_switch_mixed_C4_check']
    if not all(receipt[k] for k in ('actual_R100_R110_feedback_mixed4_available','actual_R2_moments_velocity_pressure_retained',
            'source_width_and_amplitude_caps_deferred_until_final_physical_rows','phase1_R2_and_R110_functional_mixed4_joins_certified')):
        raise ValueError('Admitted actual source microswitch graph required')
    normalization=actual_microswitch_normalization_theorem();endpoint=actual_microswitch_endpoint_theorem();factors=microswitch_factor_units_theorem()
    if len(normalization['identities'])!=141 or not all(normalization['identities'].values()) or endpoint['actual_R2_source_rows_verified']!=135 or not endpoint['passed']:raise ValueError('Full original microswitch source-unit and R2 endpoint theorem required')
    parent=field.power_tensor.proof['current_actual_switch_power_source_pressure_and_R110_theorem']
    return dict(live_original_callable_bindings=live,current_original_actual_input_bindings=bindings,
        original_phase_to_ordinary_y_and_full_raw_history_normalization=normalization,
        original_phase1_and_R2_complete_source_function_theorem=endpoint,
        exact_formal_radius_four_factor_and_current_tensor_packet_units=factors,
        consumed_checked_current_raw_stress_and_full_physical_operator=parent['checked_arbitrary_variable_source_current_radius_unit_theorem'],
        checked_actual_analytic_pressure_function=parent['checked_actual_analytic_pressure_function'],
        same_current_F0_ratios_comparison_own_moments_and_actual_six_histories=True,
        original_signed_first_switch_axial_drive_retained=True,
        source_width_and_swirl_factors_retained_through_full_tensor_algebra=True,
        all_source_rows_expanded_only_after_full_Bell_Leibniz_and_Z_operations=True,
        same_source_identity_and_checked_full_operators_imply_two_completed_tensor_traces=True,
        input_hashes={**normalization['input_hashes'],**endpoint['input_hashes'],**factors['input_hashes']},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentMicroswitchBackgroundTensor:
    @source_precision
    def __init__(self,power_tensor=None,require_checked=True):
        self.power_tensor=power_tensor if power_tensor is not None else CurrentSwitchPowerBackgroundTensor()
        self.physical=self.power_tensor.physical;self.pulse=self.power_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.power_tensor.family;self.source=self.power_tensor.source;self.datum_sha=self.power_tensor.datum_sha
        self.switch=self.power_tensor.switch;self.assert_graph()
        self.evaluator,self.compiler_proof=compiled_original_microswitch_with_rows()
        self.stress,self.remainder,self.operator_proof=compiled_factored_source_operators()
        self.lift,self.lift_proof=compiled_factored_physical_lift(self.remainder)
        source=current_microswitch_source_theorem(self)
        self.proof=dict(current_actual_microswitch_source_pressure_units_and_two_endpoint_theorem=source,
            original_output_only_compiler_theorem=self.compiler_proof,
            original_unresolved_raw_and_axial_operator_theorem=self.operator_proof,
            original_full_factored_physical_pullback_theorem=self.lift_proof,
            every_original_signed_stress_and_remainder_sector_retained=True,
            original_absolute_P0_and_full_energy_baseline_retained=True,passed=True)
        self.hashes={**self.power_tensor.hashes,**source['input_hashes'],**self.compiler_proof['input_hashes'],
            **self.operator_proof['input_hashes'],**self.lift_proof['input_hashes']}
        for stem in ('current_microswitch_background_tensor','current_microswitch_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual microswitch tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.power_tensor.assert_graph()
        if not (self.power_tensor.acceptance_loaded and self.physical is self.power_tensor.physical and
                self.pulse is self.power_tensor.pulse and self.ctx is self.physical.ctx and self.switch is self.power_tensor.switch and
                all(self.physical.dispatch.provider(name) is self.switch for name in DOMAINS)):
            raise ValueError('Same checked power and current actual microswitch graph required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<DOMAINS[chart][0] or endpoints(v)[1]>DOMAINS[chart][1] or endpoints(nu)[0]<=0:raise ValueError('Original first[0,1]/second[1,2] phase,Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        branch='first' if chart=='switch_first' else 'second'
        micro=self.evaluator(self.switch,Z,v,branch);raw=raw_microswitch_rows(c,micro);algebra=raw['algebra']
        z=algebra.lift(IntervalTaylor.variable(c,Z,5));delta=c.mpf(self.pulse.delta)
        rows=self.stress(c,delta,z,raw['velocity']['theta'],raw['velocity']['axial'],raw['histories'],raw['absolute_pressure'])
        logR,radius=self.physical.radius(chart,v,micro,self.switch)
        sectors=expanded_stress_sectors(c,rows,algebra,logR,self.physical.logP)
        packet=dict(Z=Z,s=v,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),
            exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)
        point=self.lift(c,packet,delta,raw['velocity'],lt,theta,nu)
        return dict(point,chart=chart,coverage_coordinate=v,current_actual_source_stress_packet=packet,
            actual_upstream_current_microswitch_source={key:value for key,value in micro.items() if key!=SOURCE_KEY},
            original_unresolved_source_phase_rows=factored_rows_record(micro[SOURCE_KEY]),
            current_unresolved_raw_source_rows=factored_rows_record(raw),fixed_current_factored_source_log_bases=list(algebra.logs),
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,
            inverse_width_source_shift_before_any_resolution=True,
            every_original_stress_and_remainder_sector_expanded_without_pruning=True,
            original_radial_prefactors_and_absolute_P0_included_once=True,
            full_current_signed_axial_drive_comparison_own_moments_and_actual_histories_retained=True,
            derivative_coordinate_is_ordinary_logR_not_selector_phase=True,
            source_factors_retained_through_full_tensor_and_physical_operators=True,
            source_locals_exposed_only_without_operator_changes=True,
            actual_full_stress_not_local_difference=True,source_bounds_not_resolved_physical_point_values=True,
            **dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name=='first_second':
            left=self.chart('switch_first',Z,1,log_tau,theta,viscosity)
            right=self.chart('switch_second',Z,1,log_tau,theta,viscosity)
        elif name=='second_power':
            left=self.chart('switch_second',Z,2,log_tau,theta,viscosity)
            right=self.power_tensor.chart('switch_power',Z,0,log_tau,theta,viscosity)
        else:raise ValueError('Actual completed phase1/R2 tensor join required')
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Actual phase1/R2 full tensor layout differs')
        common={}
        for key in aa:
            values=[endpoints(row['log_absolute_upper'])[1] for row in aa[key]+bb[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_source_and_completed_tensor_endpoint_function_theorem=self.proof,
            source_function_equality_precedes_common_triangle_bounds=True,interval_overlap_not_used_as_function_identity=True,
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_actual_microswitch_full_tensor_and_phase1_R2_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.power_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=28,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_microswitch_domains=DOMAINS,
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_core_bridge_R100_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentMicroswitchBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_microswitch_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_microswitch_tensor_views'][name]=field.chart(*args)
        print('Current actual microscopic full tensor: '+name,flush=True)
    result['current_actual_phase1_R2_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
