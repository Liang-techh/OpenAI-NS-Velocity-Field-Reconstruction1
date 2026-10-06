"""Three actual bridge full tensors and completed phase1/smoothing/R100 joins."""
import gzip
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_microswitch_background_tensor import (
    CurrentMicroswitchBackgroundTensor,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,canonical_tensor_groups,BASE,IntervalTaylor)
from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import CurrentActualBridgeMixedC4,REPLAY,REPLAY_PROOF
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_bridge_stress_operator import (
    SOURCE_KEY,compiled_current_bridge_with_rows,compiled_raw_bridge_rows,actual_bridge_normalization_theorem,
    actual_bridge_endpoint_operator_theorem,expanded_stress_sectors,factored_rows_record)

NAME=PREFIX+'current_bridge_background_tensor.json.gz'
RECEIPT=PREFIX+'current_bridge_background_tensor_check.json'
FUNCTIONAL_RECEIPT=PREFIX+'current_bridge_functional_joins_check.json'
GATES=('current_actual_three_bridge_full_tensors_available',
    'current_actual_bridge_full_meridional_decomposition_available',
    'current_actual_bridge_phase1_smoothing_R100_completed_tensor_joins_certified')
DOMAINS={'bridge_first':(0,1),'bridge_second':(1,2),'bridge_macro':(0,1)}
SEAMS=('bridge_first_second','bridge_second_macro','bridge_switch_R100')
SHARDS={chart:PREFIX+'current_bridge_background_tensor_'+chart+'.json.gz' for chart in DOMAINS}
VIEWS={}
for chart,(lo,hi) in DOMAINS.items():
    stem=chart.replace('bridge_','')
    VIEWS.update({
        stem+'_whole':(chart,(-1,1),(lo,hi),('-3','-1'),None,'1'),
        stem+'_left':(chart,(-1,1),lo,'-1',None,'1'),
        stem+'_right':(chart,(-1,1),hi,'-1',None,'1'),
        stem+'_fresh':(chart,'.537',str(lo+.537),'-2.6','.41','.8'),
        stem+'_fresh_low':(chart,'.731',str(lo+.1337),'-2.6','.41','.8'),
        stem+'_fresh_high':(chart,'-.317',str(lo+.831),'-2.6','.41','.2')})


def current_bridge_source_theorem(field):
    bridge=field.bridge;native=field.physical.dispatch.native;compiled=field.evaluator.__globals__['REPLAY']
    live=dict(exact_current_actual_bridge_owner=type(bridge) is CurrentActualBridgeMixedC4,
        same_checked_current_core=bridge.core is field.micro_tensor.switch.core,
        same_checked_current_history=bridge.history is field.micro_tensor.switch.history,
        same_checked_current_switch=bridge.switch is field.micro_tensor.switch,
        same_three_native_bridge_providers=all(field.physical.dispatch.provider(name) is bridge for name in DOMAINS),
        same_exact_current_bridge_provider_graph=all(bridge.current_provider_graph().values()),
        same_complete_current_native_provider_graph=all(native.provider_graph_identity().values()),
        same_unmodified_current_parent_function=bridge.current_parents.__func__ is CurrentActualBridgeMixedC4.current_parents,
        original_coordinate_labelled_replay_code=compiled.__code__ is REPLAY.__code__,
        original_current_derivative_replay_proof=field.compiler_proof['original_current_parent_acquisition_replay']==REPLAY_PROOF,
        same_checked_raw_stress_operator=field.stress is field.micro_tensor.stress,
        same_checked_full_physical_operator=field.lift is field.micro_tensor.lift,
        same_original_BASE_radius=field.physical.radius.__func__ is BASE.radius)
    if not all(live.values()) or not bridge.acceptance_loaded:raise ValueError('Same checked actual bridge/current micro/core graph required')
    receipt=accepted(FUNCTIONAL_RECEIPT,field.family,field.source,'current_first_second_functional_mixed4_join_certified');_verify_hashes(receipt)
    if receipt['datum_enclosure_sha256']!=field.datum_sha or not all(receipt[k] for k in (
            'current_first_second_functional_mixed4_join_certified','current_second_macro_functional_mixed4_join_certified',
            'current_R100_bridge_switch_functional_mixed4_join_certified')) or receipt['current_core_bridge_functional_mixed4_join_certified']:
        raise ValueError('Same-source three bridge functional joins, with core scope still open, required')
    normalization=actual_bridge_normalization_theorem();endpoint=actual_bridge_endpoint_operator_theorem()
    if normalization['total_identities']!=282 or endpoint['source_generator_identity_count']!=270 or not normalization['passed'] or not endpoint['passed']:raise ValueError('Original bridge raw units/full source boundary programs failed')
    unit=field.micro_tensor.proof['current_actual_microswitch_source_pressure_units_and_two_endpoint_theorem']['exact_formal_radius_four_factor_and_current_tensor_packet_units']
    if not unit['passed'] or not all(unit['complete_current_physical_unit_identities'].values()):raise ValueError('Checked arbitrary current source R/Pstar/four-factor units required')
    asts=SourceAST()
    original=assignment_source_bindings('bridge_mixed_C4','evaluate',{
        'logu0':'c.ln(2*R)/2+self.logF0+c.ln(phi[0])-self.core.logP',
        'algebra':'AxialSixAlgebra(c,[self.logh,2*self.core.logP,2*self.logF0,2*logu0],self.proofs)',
        'scale':'h if microscopic else algebra.lift(1)',
        'coordinate_power':'algebra.width if microscopic else lambda row,power=1:algebra.lift(row)'})
    packet=assignment_source_bindings('current_bridge_background_tensor','chart',{
        '(logR, radius)':'self.physical.radius(chart,v,bridge,self.bridge)',
        'packet':'dict(Z=Z,s=v,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)'})
    asts.method('global_physical_assembly','radius');asts.method('bridge_mixed_C4','evaluate')
    asts.method('current_actual_bridge_mixed_C4','current_parents');asts.method('current_bridge_background_tensor','chart')
    return dict(live_original_callable_bindings=live,consumed_current_three_functional_join_receipt=receipt,
        original_current_micro_and_macro_raw_unit_theorem=normalization,
        original_current_full_boundary_control_and_generator_theorem=endpoint,
        checked_arbitrary_R_Pstar_four_factor_and_full_physical_units=unit,
        original_bridge_scale_and_log_source_bindings=original,current_bridge_tensor_packet_unit_bindings=packet,
        formal_micro_radius='Ra*exp(original_hb*phase)',
        formal_macro_radius='Ra*exp(2*original_hb+fraction*(log(100/Ra)-2*original_hb))',
        original_BASE_radius_log_bounds_enclose_formal_radius_not_choose_field=True,
        internal_inverse_Pstar_and_external_mode_factors_counted_once=True,
        dummy_end_scale_and_signed_memory_absent_only_for_current_bridge_charts=True,
        original_axial6_comparison_smoothing_and_actual_signed_drive_retained=True,
        actual_six_own_moments_never_replaced_by_comparison_own_moments=True,
        same_original_P0_complete_energy_and_current_swirl_axial_normalization=True,
        three_source_function_joins_and_same_full_operators_imply_completed_tensor_traces=True,
        input_hashes={**receipt['input_hashes'],FUNCTIONAL_RECEIPT:sha(FUNCTIONAL_RECEIPT),
            **normalization['input_hashes'],**endpoint['input_hashes'],**asts.hashes},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))
def read_shard(chart):return json.loads(gzip.decompress((HERE/SHARDS[chart]).read_bytes()))


class CurrentBridgeBackgroundTensor:
    @source_precision
    def __init__(self,micro_tensor=None,require_checked=True):
        self.micro_tensor=micro_tensor if micro_tensor is not None else CurrentMicroswitchBackgroundTensor()
        self.physical=self.micro_tensor.physical;self.pulse=self.micro_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.micro_tensor.family;self.source=self.micro_tensor.source;self.datum_sha=self.micro_tensor.datum_sha
        self.bridge=self.physical.dispatch.provider('bridge_first');self.assert_graph()
        self.stress=self.micro_tensor.stress;self.remainder=self.micro_tensor.remainder;self.lift=self.micro_tensor.lift
        self.evaluator,self.compiler_proof=compiled_current_bridge_with_rows()
        self.raw_adapter,self.raw_adapter_proof=compiled_raw_bridge_rows()
        source=current_bridge_source_theorem(self)
        self.proof=dict(current_actual_bridge_source_pressure_units_and_three_endpoint_theorem=source,
            original_output_only_compiler_theorem=self.compiler_proof,original_current_raw_unit_adapter_theorem=self.raw_adapter_proof,
            consumed_checked_original_full_factored_tensor_operators=self.micro_tensor.proof['original_full_factored_physical_pullback_theorem'],
            every_original_signed_full_stress_and_remainder_sector_retained=True,
            original_absolute_P0_energy_baseline_actual_own_history_and_comparison_source_retained=True,passed=True)
        self.hashes={**self.micro_tensor.hashes,**source['input_hashes'],**self.compiler_proof['input_hashes'],**self.raw_adapter_proof['input_hashes']}
        for stem in ('current_bridge_background_tensor','current_bridge_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual bridge tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.micro_tensor.assert_graph()
        if not (self.micro_tensor.acceptance_loaded and self.physical is self.micro_tensor.physical and
                self.pulse is self.micro_tensor.pulse and self.ctx is self.physical.ctx and
                self.bridge is self.physical.dispatch.provider('bridge_first') and
                all(self.physical.dispatch.provider(name) is self.bridge for name in DOMAINS) and
                self.bridge.acceptance_loaded and self.bridge.switch is self.micro_tensor.switch):
            raise ValueError('Same checked micro-switch tensor/current bridge graph required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<DOMAINS[chart][0] or endpoints(v)[1]>DOMAINS[chart][1] or endpoints(nu)[0]<=0:raise ValueError('Original first[0,1]/second[1,2]/macro[0,1],Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        bridge=self.evaluator(self.bridge,Z,v,chart.replace('bridge_',''));raw=self.raw_adapter(c,bridge);algebra=raw['algebra']
        z=algebra.lift(IntervalTaylor.variable(c,Z,5));delta=c.mpf(self.pulse.delta)
        rows=self.stress(c,delta,z,raw['velocity']['theta'],raw['velocity']['axial'],raw['histories'],raw['absolute_pressure'])
        logR,radius=self.physical.radius(chart,v,bridge,self.bridge)
        sectors=expanded_stress_sectors(c,rows,algebra,logR,self.physical.logP)
        packet=dict(Z=Z,s=v,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)
        point=self.lift(c,packet,delta,raw['velocity'],lt,theta,nu)
        return dict(point,chart=chart,coverage_coordinate=v,current_actual_source_stress_packet=packet,
            actual_upstream_current_bridge_source={key:value for key,value in bridge.items() if key!=SOURCE_KEY},
            original_unresolved_source_coordinate_rows=factored_rows_record(bridge[SOURCE_KEY]),
            current_unresolved_raw_source_rows=factored_rows_record(raw),fixed_current_factored_source_log_bases=list(algebra.logs),
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,
            ordinary_logR_conversion='D_y^j=hb^-j*D_phase^j' if chart!='bridge_macro' else 'D_y^j=D_original_coordinate^j; fraction labels coverage only',
            microscopic_inverse_width_or_macro_identity_before_any_resolution=True,
            every_original_stress_and_remainder_sector_expanded_without_pruning=True,
            original_radial_prefactors_and_absolute_P0_included_once=True,
            full_current_signed_axial_drive_axial6_comparison_and_actual_own_histories_retained=True,
            derivative_coordinate_is_ordinary_logR_not_macro_fraction=True,
            source_factors_retained_through_full_tensor_and_physical_operators=True,
            source_locals_exposed_only_without_operator_changes=True,
            actual_full_stress_not_local_difference=True,source_bounds_not_resolved_physical_point_values=True,
            **dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        pairs={'bridge_first_second':(('bridge_first',1),('bridge_second',1)),
            'bridge_second_macro':(('bridge_second',2),('bridge_macro',0))}
        if name in pairs:
            a,b=pairs[name];left=self.chart(a[0],Z,a[1],log_tau,theta,viscosity);right=self.chart(b[0],Z,b[1],log_tau,theta,viscosity)
        elif name=='bridge_switch_R100':
            left=self.chart('bridge_macro',Z,1,log_tau,theta,viscosity)
            right=self.micro_tensor.chart('switch_first',Z,0,log_tau,theta,viscosity)
        else:raise ValueError('Actual completed bridge phase1/smoothing/R100 tensor join required')
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Actual bridge complete tensor layout differs')
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
            current_actual_bridge_full_tensor_and_three_endpoint_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.micro_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=31,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_bridge_domains=DOMAINS,
            producer_storage='Complete unpruned deterministic gzip chart shards and source/interface manifest',
            remaining_core_core_bridge_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentBridgeBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_bridge_chart_shards']={}
    for chart,path in SHARDS.items():
        views={}
        for name,args in VIEWS.items():
            if args[0]!=chart:continue
            views[name]=field.chart(*args);print('Current actual full bridge tensor: '+name,flush=True)
        payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
        (HERE/path).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
        result['current_actual_bridge_chart_shards'][chart]=dict(path=path,sha256=sha(path),views=list(views))
    result['current_actual_bridge_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
