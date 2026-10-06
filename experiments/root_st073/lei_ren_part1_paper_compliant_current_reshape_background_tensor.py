"""Actual whole long-reshape tensor and its completed Rsh attachment."""
import gzip
import json
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_restore_background_tensor import (
    CurrentRestoreBackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,endpoints,
    source_precision,accepted,_verify_hashes,canonical_tensor_groups,BASE,PULSE)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_compliant_current_reshape_stress_operator import raw_reshape_rows,actual_reshape_normalization_theorem
from lei_ren_part1_paper_compliant_actual_long_reshape_mixed_C4 import (
    CompliantActualLongReshapeMixedC4,CompliantActualLongReshapeProfiles,profile_source_bindings)
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import (
    CompliantLongReshapeMixedC4,reshape_mixed,exponential_derivatives,scaled_positive_source)
from lei_ren_part1_paper_compliant_long_reshape_profiles import CompliantLongReshapeProfiles
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import current_E_source_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import source_bindings,boundary_identities
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

NAME=PREFIX+'current_reshape_background_tensor.json.gz'
RECEIPT=PREFIX+'current_reshape_background_tensor_check.json'
GATES=('current_actual_long_reshape_full_tensor_available',
    'current_actual_long_reshape_full_meridional_decomposition_available',
    'current_actual_Rsh_completed_tensor_join_certified')
DOMAINS={'reshape':(0,1)}
SEAMS=('reshape_reference',)
VIEWS={'reshape_whole':('reshape',(-1,1),(0,1),('-3','-1'),None,'1'),
    'reshape_left':('reshape',(-1,1),0,'-1',None,'1'),
    'reshape_right':('reshape',(-1,1),1,'-1',None,'1'),
    'reshape_fresh':('reshape','.537','.537','-2.6','.41','.8'),
    'reshape_fresh_low':('reshape','.731','.1337','-2.6','.41','.8'),
    'reshape_fresh_high':('reshape','-.317','.831','-2.6','.41','.2')}


def current_reshape_source_theorem(field):
    """Bind actual R110 histories, variable source jets and exact Rsh functions."""
    mixed=field.reshape;profile=mixed.reshape;restore=field.restore_tensor.restore
    native=field.physical.dispatch.native
    live=dict(exact_current_mixed_owner=type(mixed) is CompliantActualLongReshapeMixedC4,
        exact_current_profile_owner=type(profile) is CompliantActualLongReshapeProfiles,
        same_nested_long_source=mixed is restore.long_mixed,
        same_current_profile=profile is restore.reference.reshape,
        same_current_history=mixed.history is restore.reference.history,
        same_current_core=profile.core is restore.reference.core,
        same_original_inverse_Pstar_unit=endpoints(mixed.invP2)==endpoints(restore.invP2),
        same_current_V_E_namespace=mixed.shared_axial_source==restore.shared_axial_source,
        same_native_reshape_provider=field.physical.dispatch.provider('reshape') is mixed,
        same_exact_current_provider_graph=all(native.provider_graph_identity().values()),
        original_profile_inputs=profile.inputs.__func__ is CompliantLongReshapeProfiles.inputs,
        original_profile_evaluate=profile.evaluate.__func__ is CompliantLongReshapeProfiles.evaluate,
        actual_mixed_evaluate=mixed.evaluate.__func__ is CompliantActualLongReshapeMixedC4.evaluate,
        original_mixed_evaluate=mixed.original_evaluate.__func__ is CompliantLongReshapeMixedC4.evaluate,
        original_mixed_source=mixed.original_evaluate.__func__.__globals__['reshape_mixed'] is reshape_mixed,
        original_Bell=reshape_mixed.__globals__['exponential_derivatives'] is exponential_derivatives,
        original_positive_combination=raw_reshape_rows.__globals__['scaled_positive_source'] is scaled_positive_source,
        original_BASE_radius=field.physical.radius.__func__ is BASE.radius,
        actual_raw_adapter=field.chart.__func__.__wrapped__.__globals__['raw_reshape_rows'] is raw_reshape_rows,
        checked_raw_stress=field.chart.__func__.__wrapped__.__globals__['raw_pre_stress_rows'] is raw_pre_stress_rows)
    if not all(live.values()):raise ValueError('Same actual R110/reshape/reference source graph required')
    long=native.receipts['actual_long_reshape_mixed_C4_check']
    rsh=native.receipts['actual_Rsh_source_join_check']
    if not all(long[k] for k in ('current_actual_long_reshape_mixed4_available',
            'current_actual_R110_feedback_in_long_reshape','original_selected_A_T_logref_and_pressure_retained')):raise ValueError('Current original R110 histories and selected reshape source required')
    if not all(rsh[k] for k in ('current_Rsh_source_functional_join_certified',
            'current_Rsh_parent_P0_V_E_and_physical_normalization_verified',
            'source_equality_proof_not_enclosure_overlap','constant_power_boundary_extension_does_not_replace_left_field')):raise ValueError('Exact current Rsh source theorem required')
    if not (mixed.original_check['symbolic_checks']['passed'] and mixed.original_check['independent_physical_fixture']['passed']):raise ValueError('Unchanged mixed source operator and independent physical fixture required')
    if endpoints(profile.T)!=endpoints(400*profile.A):raise ValueError('Actual original T=400*Abar changed')
    asts=SourceAST()
    asts.expression('long_reshape_mixed_C4','evaluate','log_y',
        wanted="[B*(-cutoff[k]*math.factorial(k)/T**k)+(c.mpf('.1') if k==1 else 0) for k in range(1,5)]")
    asts.expression('long_reshape_profiles','evaluate','moment_shapes',
        wanted="dict(theta=inherited['theta']+kernels['theta'],theta_z=inherited['theta_z']+inp['v']*kernels['theta'],pressure=inherited['pressure']+kernels['pressure'],swirl=inherited['swirl']+kernels['swirl'],mean=mean,axial=axial)")
    asts.method('long_reshape_profiles','backward_kernel')
    asts.method('long_reshape_profiles','normalized_amplitude_decay')
    asts.method('long_reshape_mixed_C4','scaled_positive_source')
    asts.expression('current_reshape_background_tensor','chart','reshape',wanted='self.reshape.evaluate(Z,v)')
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    T,D=s.symbols('actual_reshape_T actual_loggap',real=True)
    provider=SimpleNamespace(reshape=SimpleNamespace(T=T),reference=SimpleNamespace(reshape=SimpleNamespace(T=T),loggap=D))
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify,ln=s.log),logRref=s.log(110)+T+D)
    if s.simplify(radius(owner,'reshape',1,{},provider)[0]-radius(owner,'inner_reference',0,{},provider)[0])!=0:raise ArithmeticError('Actual Rsh radius sources differ')
    flat={}
    for phase in (0,1):
        cutoff=sigma_jets(field.ctx,field.ctx.mpf(phase))
        if endpoints(cutoff[0])!=(mp.mpf(phase),mp.mpf(phase)) or any(endpoints(cutoff[j])!=(mp.mpf(0),mp.mpf(0)) for j in range(1,5)):raise ValueError('Original flat reshape endpoint jets changed')
        flat[str(phase)]=True
    boundary=boundary_identities();bindings=source_bindings()
    if not boundary['passed'] or boundary['total_physical_mixed4_rows_implied']!=135:raise ValueError('Exact Rsh boundary function proof failed')
    pressure=field.restore_tensor.proof['current_actual_restore_source_pressure_and_three_endpoint_theorem']['checked_actual_analytic_pressure_function']
    if not pressure['passed']:raise ValueError('Same current analytic pressure defining function required')
    return dict(live_original_callable_bindings=live,current_R110_profile_bindings=profile_source_bindings(mixed),
        current_correlated_E_source_bindings=current_E_source_bindings(restore),
        actual_original_Rsh_source_bindings=bindings,actual_original_Rsh_boundary_function_theorem=boundary,
        exact_actual_Rsh_radius_function_identity=True,original_flat_reshape_endpoint_jets=flat,
        checked_actual_analytic_pressure_function=pressure,
        actual_variable_B_T_cutoff_and_four_log_amplitude_jets_retained=True,
        full_original_positive_kernels_and_nonzero_R110_histories_retained=True,
        inverse_T_factors_and_Bell_products_not_replaced_by_constant_rates=True,
        constant_power_extension_used_only_for_boundary_not_finite_reshape=True,
        source_function_identity_precedes_completed_tensor_triangle_bounds=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes=asts.hashes,passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentReshapeBackgroundTensor:
    @source_precision
    def __init__(self,restore_tensor=None,require_checked=True):
        self.restore_tensor=restore_tensor if restore_tensor is not None else CurrentRestoreBackgroundTensor()
        self.physical=self.restore_tensor.physical;self.pulse=self.restore_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.restore_tensor.family;self.source=self.restore_tensor.source;self.datum_sha=self.restore_tensor.datum_sha
        self.reshape=self.restore_tensor.restore.long_mixed;self.assert_graph();self.lift=self.restore_tensor.lift
        source=current_reshape_source_theorem(self);normalization=actual_reshape_normalization_theorem()
        self.proof=dict(current_actual_variable_reshape_source_pressure_and_Rsh_theorem=source,
            original_actual_variable_reshape_current_radius_normalization_theorem=normalization,
            consumed_checked_arbitrary_raw_stress_and_full_physical_operator=self.restore_tensor.proof,
            full_current_R110_moments_V110_original_B_T_kernels_analytic_P0_retained=True,
            positive_caps_are_enclosures_not_defining_source_values=True,passed=True)
        self.hashes=dict(self.restore_tensor.hashes)
        for value in (source,normalization):self.hashes.update(value['input_hashes'])
        for stem in ('current_reshape_background_tensor','current_reshape_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual reshape tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.restore_tensor.assert_graph()
        if not (self.restore_tensor.acceptance_loaded and self.physical is self.restore_tensor.physical and
                self.pulse is self.restore_tensor.pulse and self.ctx is self.physical.ctx and
                self.reshape is self.restore_tensor.restore.long_mixed and
                self.physical.dispatch.provider('reshape') is self.reshape):
            raise ValueError('Same checked restore and actual nested long-reshape graph required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<0 or endpoints(v)[1]>1 or endpoints(nu)[0]<=0:raise ValueError('Actual reshape phase[0,1],Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        reshape=self.reshape.evaluate(Z,v)
        raw=raw_reshape_rows(c,reshape,c.mpf(endpoints(self.reshape.invP2)))
        z=IntervalTaylor.variable(c,Z,5);delta=c.mpf(self.pulse.delta)
        velocity=raw['velocity'];histories=raw['histories'];P=raw['absolute_pressure']
        rows=raw_pre_stress_rows(c,delta,z,velocity['theta'],velocity['axial'],histories,P)
        logR,radius=self.physical.radius(chart,v,reshape,self.reshape);sectors={}
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
            actual_upstream_current_long_reshape_source=reshape,current_raw_five_history_rows=histories,
            current_absolute_pressure_ordinary_y_rows=P,current_source_three_component_velocity_rows=velocity,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,actual_variable_amplitude_source_enclosures=raw,
            full_current_R110_moments_V110_variable_B_T_kernels_and_analytic_P0_retained=True,
            full_energy_baseline_and_absolute_P0_retained=True,
            derivative_coordinate_is_ordinary_logR_not_selector_phase=True,
            original_variable_Bell_products_combined_before_source_caps=True,
            only_radius_and_constant_Pstar_factored_from_full_stress=True,actual_full_stress_not_local_difference=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name!='reshape_reference':raise ValueError('Actual completed Rsh tensor join required')
        left=self.chart('reshape',Z,1,log_tau,theta,viscosity)
        right=self.restore_tensor.chart('inner_reference',Z,0,log_tau,theta,viscosity)
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Actual Rsh completed tensor layout differs')
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
            current_actual_variable_reshape_full_tensor_and_Rsh_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.restore_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=25,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_reshape_domains=DOMAINS,
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_core_bridge_switch_R110_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentReshapeBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_reshape_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_reshape_tensor_views'][name]=field.chart(*args)
        print('Current actual whole long reshape full tensor: '+name,flush=True)
    result['current_actual_Rsh_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
