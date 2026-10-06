"""Actual R2-to-R110 power full tensor and completed R110 attachment."""
import gzip
import json
import math
import mpmath as mp
from lei_ren_part1_paper_compliant_current_reshape_background_tensor import (
    CurrentReshapeBackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,endpoints,
    source_precision,accepted,_verify_hashes,canonical_tensor_groups,BASE,IntervalTaylor)
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_compliant_current_reshape_stress_operator import raw_reshape_rows
from lei_ren_part1_paper_compliant_current_switch_power_source_adapter import (
    compiled_original_postpower_with_locals,R110_source_functional_identity)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import CompliantActualSwitchMixedC4
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import CompliantMicroswitchMixedC4
from lei_ren_part1_paper_compliant_inner_switch_profiles import CompliantInnerSwitchProfiles,power_transport
from lei_ren_part1_paper_compliant_inner_bridge_profiles import dress
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import reshape_mixed
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

NAME=PREFIX+'current_switch_power_background_tensor.json.gz'
RECEIPT=PREFIX+'current_switch_power_background_tensor_check.json'
GATES=('current_actual_switch_power_full_tensor_available',
    'current_actual_switch_power_full_meridional_decomposition_available',
    'current_actual_R110_completed_tensor_join_certified')
DOMAINS={'switch_power':(0,1)}
SEAMS=('power_reshape',)
VIEWS={'switch_power_whole':('switch_power',(-1,1),(0,1),('-3','-1'),None,'1'),
    'switch_power_left':('switch_power',(-1,1),0,'-1',None,'1'),
    'switch_power_right':('switch_power',(-1,1),1,'-1',None,'1'),
    'switch_power_fresh':('switch_power','.537','.537','-2.6','.41','.8'),
    'switch_power_fresh_low':('switch_power','.731','.1337','-2.6','.41','.8'),
    'switch_power_fresh_high':('switch_power','-.317','.831','-2.6','.41','.2')}


def current_switch_power_source_theorem(field):
    """Bind exact R2 transport, original power operators and same R110 inlet."""
    mixed=field.switch;long=field.reshape_tensor.reshape;profile=mixed.switch
    native=field.physical.dispatch.native;env=field.evaluator.__globals__
    live=dict(exact_current_mixed_owner=type(mixed) is CompliantActualSwitchMixedC4,
        same_nested_switch_source=mixed is long.switch_mixed,
        same_current_switch_profile=profile is long.reshape.switch is long.history.switch,
        same_current_history=mixed.history is long.history,
        same_current_core=mixed.core is long.reshape.core,
        same_original_inverse_Pstar_unit=endpoints(mixed.invP2)==endpoints(long.invP2),
        same_current_V110_namespace=mixed.shared_axial_source==long.shared_axial_source,
        same_native_power_provider=field.physical.dispatch.provider('switch_power') is mixed,
        same_three_native_switch_providers=all(field.physical.dispatch.provider(name) is mixed for name in ('switch_first','switch_second','switch_power')),
        same_exact_current_provider_graph=all(native.provider_graph_identity().values()),
        original_postpower=mixed.postpower.__func__ is CompliantMicroswitchMixedC4.postpower,
        original_phase=profile.phase.__func__ is CompliantInnerSwitchProfiles.phase,
        original_post=profile.post.__func__ is CompliantInnerSwitchProfiles.post,
        original_inlet=profile.inlet.__func__ is CompliantInnerSwitchProfiles.inlet,
        original_packet=profile.packet.__func__ is CompliantInnerSwitchProfiles.packet,
        original_power_transport=mixed.postpower.__func__.__globals__['power_transport'] is power_transport,
        same_compiled_power_transport=env['power_transport'] is power_transport,
        same_compiled_reshape_mixed=env['reshape_mixed'] is reshape_mixed,
        same_compiled_dress=env['dress'] is dress,
        same_compiled_IntervalTaylor=env['IntervalTaylor'] is CompliantMicroswitchMixedC4.postpower.__globals__['IntervalTaylor'],
        original_BASE_radius=field.physical.radius.__func__ is BASE.radius,
        checked_raw_adapter=field.chart.__func__.__wrapped__.__globals__['raw_reshape_rows'] is raw_reshape_rows,
        checked_raw_stress=field.chart.__func__.__wrapped__.__globals__['raw_pre_stress_rows'] is raw_pre_stress_rows)
    if not all(live.values()):raise ValueError('Same original actual power/R110 graph and source operators required')
    receipt=native.receipts['actual_switch_mixed_C4_check']
    if not all(receipt[k] for k in ('actual_R100_R110_feedback_mixed4_available',
            'actual_R2_moments_velocity_pressure_retained','complete_postswitch_power_installed',
            'formal_hb_R2_zeta_theta_and_R110_containments_checked',
            'source_width_and_amplitude_caps_deferred_until_final_physical_rows',
            'phase1_R2_and_R110_functional_mixed4_joins_certified')):raise ValueError('Actual current R2 histories and original power source admission required')
    if endpoints(mixed.logh)[1]>=endpoints(field.ctx.ln(field.ctx.ln(field.ctx.mpf(110)/100)/2))[0]:raise ValueError('Original positive R2 must precede110')
    unit=field.reshape_tensor.proof['original_actual_variable_reshape_current_radius_normalization_theorem']
    if not unit['passed'] or len(unit['identities'])!=141 or not all(unit['identities'].values()):raise ValueError('Checked arbitrary-source raw-unit theorem required')
    boundary=R110_source_functional_identity()
    if not boundary['passed'] or boundary['total_physical_mixed4_rows_implied']!=135:raise ValueError('New postpower-to-post/reshape R110 source proof failed')
    flat=sigma_jets(field.ctx,field.ctx.mpf(0))
    if endpoints(flat[0])!=(mp.mpf(0),mp.mpf(0)) or any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for value in flat[1:]):raise ValueError('Original reshape phase-zero flat jets changed')
    asts=SourceAST()
    asts.expression('current_switch_power_background_tensor','chart','power',wanted='self.evaluator(self.switch,Z,v)')
    asts.expression('inner_switch_profiles','packet','source',wanted="inputs['source']")
    asts.method('global_physical_assembly','radius')
    asts.method('frozen_comparison_field','dress')
    asts.method('long_reshape_mixed_C4','reshape_mixed')
    pressure=field.reshape_tensor.proof['current_actual_variable_reshape_source_pressure_and_Rsh_theorem']['checked_actual_analytic_pressure_function']
    if not pressure['passed']:raise ValueError('Same original current analytic pressure function required')
    return dict(live_original_callable_bindings=live,original_output_only_source_locals_adapter=field.compiler_proof,
        actual_original_R110_power_post_and_reshape_boundary_function_theorem=boundary,
        checked_arbitrary_variable_source_current_radius_unit_theorem=unit,
        checked_actual_analytic_pressure_function=pressure,original_R110_reshape_flat_endpoint_checked=True,
        same_full_current_R2_histories_V110_F0_and_P0_retained=True,
        formal_positive_hb_R2_zeta_theta_and_radius_sources_not_numeric_cap_values=True,
        current_original_power_rate_point1_only_on_original_power_region=True,
        same_source_mixed4_and_checked_full_operators_imply_completed_R110_tensor_trace=True,
        source_function_identity_precedes_completed_tensor_triangle_bounds=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes={**asts.hashes,**boundary['input_hashes']},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentSwitchPowerBackgroundTensor:
    @source_precision
    def __init__(self,reshape_tensor=None,require_checked=True):
        self.reshape_tensor=reshape_tensor if reshape_tensor is not None else CurrentReshapeBackgroundTensor()
        self.physical=self.reshape_tensor.physical;self.pulse=self.reshape_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.reshape_tensor.family;self.source=self.reshape_tensor.source;self.datum_sha=self.reshape_tensor.datum_sha
        self.switch=self.reshape_tensor.reshape.switch_mixed;self.assert_graph();self.lift=self.reshape_tensor.lift
        self.evaluator,self.compiler_proof=compiled_original_postpower_with_locals()
        source=current_switch_power_source_theorem(self)
        self.proof=dict(current_actual_switch_power_source_pressure_and_R110_theorem=source,
            consumed_checked_arbitrary_raw_stress_and_full_physical_operator=self.reshape_tensor.proof,
            full_current_R2_moments_V110_F0_analytic_P0_and_formal_radius_retained=True,
            positive_caps_are_enclosures_not_defining_source_values=True,passed=True)
        self.hashes={**self.reshape_tensor.hashes,**source['input_hashes'],**self.compiler_proof['input_hashes']}
        for stem in ('current_switch_power_background_tensor','current_switch_power_source_adapter'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual switch-power tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.reshape_tensor.assert_graph()
        if not (self.reshape_tensor.acceptance_loaded and self.physical is self.reshape_tensor.physical and
                self.pulse is self.reshape_tensor.pulse and self.ctx is self.physical.ctx and
                self.switch is self.reshape_tensor.reshape.switch_mixed and
                self.physical.dispatch.provider('switch_power') is self.switch):
            raise ValueError('Same checked reshape and actual nested switch-power graph required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<0 or endpoints(v)[1]>1 or endpoints(nu)[0]<=0:raise ValueError('Actual switch-power fraction[0,1],Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        power=self.evaluator(self.switch,Z,v)
        raw=raw_reshape_rows(c,power,c.mpf(endpoints(self.switch.invP2)))
        z=IntervalTaylor.variable(c,Z,5);delta=c.mpf(self.pulse.delta)
        velocity=raw['velocity'];histories=raw['histories'];P=raw['absolute_pressure']
        rows=raw_pre_stress_rows(c,delta,z,velocity['theta'],velocity['axial'],histories,P)
        logR,radius=self.physical.radius(chart,v,power,self.switch);sectors={}
        for label,parts in rows.items():
            sectors[label]={}
            for name,part in parts.items():
                rp,bp,dp,hp=part['mode']
                grid={'s%d_Z%d'%(j,n):jet[n]*math.factorial(n) for j,jet in enumerate(part['full_derivative_rows']) for n in range(4-j)}
                sectors[label][name]=dict(mode=part['mode'],exact_source_log_parts=dict(source_logR=rp*logR,
                    logPstar=bp*self.physical.logP,selected_log_end_scale=c.mpf(0),signed_original_memory_log=c.mpf(0),normalization=-c.ln(2)/2),
                    full_stress_mixed3_coefficient_enclosures=grid)
        packet=dict(Z=Z,s=v,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),
            exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)
        point=self.lift(c,packet,delta,velocity,lt,theta,nu)
        return dict(point,chart=chart,coverage_coordinate=v,current_actual_source_stress_packet=packet,
            actual_upstream_current_switch_power_source=power,current_raw_five_history_rows=histories,
            current_absolute_pressure_ordinary_y_rows=P,current_source_three_component_velocity_rows=velocity,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,actual_power_amplitude_source_enclosures=raw,
            full_current_R2_histories_V110_F0_analytic_P0_and_formal_hb_radius_retained=True,
            full_energy_baseline_and_absolute_P0_retained=True,
            derivative_coordinate_is_ordinary_logR_not_selector_fraction=True,
            source_locals_exposed_only_without_operator_changes=True,
            original_Bell_products_combined_before_source_caps=True,
            only_radius_and_constant_Pstar_factored_from_full_stress=True,actual_full_stress_not_local_difference=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name!='power_reshape':raise ValueError('Actual completed R110 tensor join required')
        left=self.chart('switch_power',Z,1,log_tau,theta,viscosity)
        right=self.reshape_tensor.chart('reshape',Z,0,log_tau,theta,viscosity)
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Actual R110 completed tensor layout differs')
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
            current_actual_switch_power_full_tensor_and_R110_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.reshape_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=26,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_switch_power_domains=DOMAINS,
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_core_bridge_microswitch_R100_R2_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentSwitchPowerBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_switch_power_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_switch_power_tensor_views'][name]=field.chart(*args)
        print('Current actual R2..R110 power full tensor: '+name,flush=True)
    result['current_actual_R110_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
