"""Actual current reference/restoration tensors attached to the five-moment patch."""
import gzip
import json
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_actual_patch_background_tensor import (
    CurrentActualPatchBackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,endpoints,
    source_precision,accepted,_verify_hashes,canonical_tensor_groups,BASE,PULSE)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_compliant_current_restore_stress_operator import raw_restore_rows,actual_restore_normalization_theorem
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import (
    CompliantActualReferenceRestoreMixedC4,CompliantActualReferenceRestoreProfiles,current_E_source_bindings)
from lei_ren_part1_paper_compliant_reference_restore_mixed_C4 import CompliantReferenceRestoreMixedC4,source_mixed
from lei_ren_part1_paper_compliant_reference_restore_profiles import CompliantReferenceRestoreProfiles,restore_centered
from lei_ren_part1_paper_compliant_current_matched_source_dispatcher import Rm_functional_identity
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z as SYMBOLIC_Z
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

NAME=PREFIX+'current_restore_background_tensor.json.gz'
RECEIPT=PREFIX+'current_restore_background_tensor_check.json'
GATES=('current_actual_inner_reference_axial_restore_buffer_full_tensors_available',
    'current_actual_restore_full_meridional_decomposition_available',
    'current_actual_three_restore_completed_tensor_joins_certified')
DOMAINS={'inner_reference':(0,1),'axial_restore':(0,1),'restore_buffer':(-7,-6)}
SEAMS=('reference_restore','restore_buffer','buffer_patch')
VIEWS={}
for _chart,_domain in DOMAINS.items():
    VIEWS[_chart+'_whole']=(_chart,(-1,1),_domain,('-3','-1'),None,'1')
    VIEWS[_chart+'_left']=(_chart,(-1,1),_domain[0],'-1',None,'1')
    VIEWS[_chart+'_right']=(_chart,(-1,1),_domain[1],'-1',None,'1')
    VIEWS[_chart+'_fresh']=(_chart,'.537','-6.337' if _chart=='restore_buffer' else '.537','-2.6','.41','.8')


def current_restore_source_theorem(field):
    """Same live centered histories/E/P0, original cutoff and actual Rm join."""
    restore=field.restore;ref=restore.reference;patch=field.patch_tensor.patch;native=field.physical.dispatch.native
    live=dict(exact_current_mixed_owner=type(restore) is CompliantActualReferenceRestoreMixedC4,
        exact_current_profile_owner=type(ref) is CompliantActualReferenceRestoreProfiles,
        same_current_patch_nested_restore=restore is patch.patch.reference_mixed,
        same_actual_patch_reference=ref is patch.patch.reference,same_current_core=ref.core is patch.patch.core,
        same_all_three_native_providers=all(field.physical.dispatch.provider(name) is restore for name in DOMAINS),
        same_exact_current_provider_graph=all(native.provider_graph_identity().values()),
        original_packet=restore.packet.__func__ is CompliantReferenceRestoreMixedC4.packet,
        original_source_mixed=restore.packet.__func__.__globals__['source_mixed'] is source_mixed,
        original_reference_branch=restore.reference_branch.__func__ is CompliantReferenceRestoreMixedC4.reference_branch,
        original_restoration=restore.restoration.__func__ is CompliantReferenceRestoreMixedC4.restoration,
        original_postrestore=restore.postrestore.__func__ is CompliantReferenceRestoreMixedC4.postrestore,
        original_sigma=restore.restoration.__func__.__globals__['sigma_jets'] is sigma_jets,
        original_centered_restoration=ref.restoration.__func__.__globals__['restore_centered'] is restore_centered,
        current_correlated_inputs=ref.inputs.__func__ is CompliantActualReferenceRestoreProfiles.inputs,
        original_BASE_radius=field.physical.radius.__func__ is BASE.radius,
        actual_raw_adapter=field.chart.__func__.__wrapped__.__globals__['raw_restore_rows'] is raw_restore_rows,
        checked_raw_stress=field.chart.__func__.__wrapped__.__globals__['raw_pre_stress_rows'] is raw_pre_stress_rows)
    for name in ('reference_centered','packet','reference','kernels','restoration','terminal','defects'):
        live['original_profile_'+name]=getattr(ref,name).__func__ is getattr(CompliantReferenceRestoreProfiles,name)
    if not all(live.values()):raise ValueError('Same original current reference/restore/patch source graph required')
    current=native.receipts['actual_reference_restore_mixed_C4_check']
    if not all(current[k] for k in ('current_actual_reference_restore_mixed4_available','current_correlated_E_installed',
            'current_original_restore_end_exact_4Z','current_E_source_C2_theorem_verified_before_intersection')):raise ValueError('Current correlated E and original source theorem required')
    if not (restore.original_check['structural_checks']['passed'] and restore.original_check['independent_physical_fixture']['passed']
            and native.join_proof['passed'] and native.receipts['actual_feedback_patch_mixed_C4_check']['current_Rm_source_transport_and_P0_retained']):raise ValueError('Checked unchanged source operator/functional Rm transport required')
    asts=SourceAST();bindings=current_E_source_bindings(restore)
    wanted={
        ('reference_restore_profiles','reference','gap'):'(self.loggap-8)*phase',
        ('reference_restore_profiles','reference','centered'):'self.reference_centered(Z,gap)',
        ('reference_restore_profiles','restoration','initial'):'self.reference_centered(Z,self.loggap-8)',
        ('reference_restore_profiles','restoration','centered'):"restore_centered(c,initial,inp['E'],t,kernels)",
        ('reference_restore_profiles','restoration','logu'):"-logarithm(1+square(inp['z']))+(-8+t)/10",
        ('reference_restore_profiles','terminal','initial'):'self.reference_centered(Z,self.loggap-8)',
        ('reference_restore_profiles','terminal','restored'):"restore_centered(c,initial,inp['E'],c.mpf(1),self.kernels(1))",
        ('reference_restore_profiles','terminal','length'):'offset+7',
        ('reference_restore_profiles','terminal','rates'):"dict(mean_error=1,angular_error='1.6',mixed_error='1.6',axial_square=1,swirl_error='1.2',pressure_error='.2')",
        ('reference_restore_profiles','terminal','logu'):"-logarithm(1+square(inp['z']))+offset/10",
        ('reference_restore_mixed_C4','packet','p0'):"jet('pressure_axis_axial5_coefficients')",
        ('actual_moment_patch','evaluate','p0'):"self.reference.inputs(Z)['original_axis_pressure']",
        ('long_reshape_profiles','__init__','self.logref'):'10*(self.core.logC+self.core.logP)',
        ('global_physical_assembly','__init__','self.logRref'):'c.ln(110)+10*(self.logC+self.logP)',
    }
    for (stem,method,target),value in wanted.items():asts.expression(stem,method,target,wanted=value)
    # Original finite kernels are exactly neutral at t=0. This identifies
    # all six centered functions, not merely intersecting output boxes.
    z=SYMBOLIC_Z;c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp)
    J=lambda value:FunctionJet.function(c,value,5)
    initial={name:J(s.Function('same_'+name)(z)) for name in ('mean_error','angular_error','mixed_error','axial_square','swirl_error','pressure_error')}
    fn=asts.replay('reference_restore_profiles','restore_centered',dict(square=lambda value:value*value))
    start=fn(c,initial,J(s.Function('same_current_E')(z)),0,dict(mean=0,mixed=0,square=0))
    neutral={name:s.cancel(start[name].expr-initial[name].expr)==0 for name in initial}
    if not all(neutral.values()):raise ArithmeticError('Original restore inlet is not the same centered function')
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    T,D=s.symbols('exact_reshape_T exact_loggap',real=True);rr=s.log(110)+T+D
    provider=SimpleNamespace(reference=SimpleNamespace(reshape=SimpleNamespace(T=T),loggap=D))
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify,ln=s.log),logRref=rr)
    pairs={'reference_restore':(('inner_reference',1),('axial_restore',0)),
        'restore_buffer':(('axial_restore',1),('restore_buffer',-7)),
        'buffer_patch':(('restore_buffer',-6),('actual_patch',1))}
    radii={}
    for name,(a,b) in pairs.items():
        if s.simplify(radius(owner,*a,{},provider)[0]-radius(owner,*b,{},provider)[0])!=0:raise ArithmeticError('Actual restoration radius differs: '+name)
        radii[name]=True
    sigma={}
    for t in (0,1):
        jet=sigma_jets(field.ctx,field.ctx.mpf(t))
        if endpoints(jet[0])!=(mp.mpf(t),mp.mpf(t)) or any(endpoints(jet[j])!=(mp.mpf(0),mp.mpf(0)) for j in range(1,5)):raise ArithmeticError('Original sigma flat endpoint changed')
        sigma[str(t)]=True
    rm=Rm_functional_identity()
    pressure=field.patch_tensor.proof['current_actual_patch_source_pressure_support_and_Rh_theorem']['actual_current_pressure_defining_function_bridge']
    if not pressure['passed']:raise ValueError('Same actual analytic P0 defining function required')
    return dict(live_original_callable_bindings=live,current_correlated_E_source_bindings=bindings,
        original_restore_start_six_centered_function_identities=neutral,three_exact_source_radius_identities=radii,
        original_sigma_flat_endpoint_jets=sigma,actual_current_Rm_open_neighborhood_function_theorem=rm,
        checked_actual_analytic_pressure_function=pressure,
        reference_Rz_endpoint_uses_exact_negative_point8_log_amplitude=True,
        restore_end_and_buffer_start_share_same_full_original_kernels_and_six_histories=True,
        original_positive_swirl_log_and_caps_enclose_source_not_replace_it=True,
        full_current_centered_energy_baseline_and_cross_terms_retained=True,
        same_source_mixed4_and_checked_original_operators_imply_three_completed_tensor_traces=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes=asts.hashes,passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentRestoreBackgroundTensor:
    @source_precision
    def __init__(self,patch_tensor=None,require_checked=True):
        self.patch_tensor=patch_tensor if patch_tensor is not None else CurrentActualPatchBackgroundTensor()
        self.physical=self.patch_tensor.physical;self.pulse=self.patch_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.patch_tensor.family;self.source=self.patch_tensor.source;self.datum_sha=self.patch_tensor.datum_sha
        self.restore=self.patch_tensor.patch.patch.reference_mixed;self.assert_graph();self.lift=self.patch_tensor.lift
        source=current_restore_source_theorem(self);normalization=actual_restore_normalization_theorem()
        self.proof=dict(current_actual_restore_source_pressure_and_three_endpoint_theorem=source,
            original_actual_restore_current_radius_normalization_theorem=normalization,
            consumed_checked_arbitrary_raw_stress_and_full_physical_operator=self.patch_tensor.proof,
            full_current_correlated_E_centered_histories_cutoff_axial5_and_analytic_P0_retained=True,
            positive_amplitude_caps_are_enclosures_not_defining_source_values=True,passed=True)
        self.hashes=dict(self.patch_tensor.hashes)
        for value in (source,normalization):self.hashes.update(value['input_hashes'])
        for stem in ('current_restore_background_tensor','current_restore_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual restore tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.patch_tensor.assert_graph()
        if not (self.patch_tensor.acceptance_loaded and self.physical is self.patch_tensor.physical and
                self.pulse is self.patch_tensor.pulse and self.ctx is self.physical.ctx and
                self.restore is self.patch_tensor.patch.patch.reference_mixed and
                all(self.physical.dispatch.provider(name) is self.restore for name in DOMAINS)):
            raise ValueError('Same checked patch and current nested reference/restore graph required')

    def source_packet(self,chart,Z,v):
        if chart=='inner_reference':return self.restore.reference_branch(Z,v)
        if chart=='axial_restore':return self.restore.restoration(Z,v)
        if chart=='restore_buffer':return self.restore.postrestore(Z,v)
        raise ValueError('Actual inner reference/restoration/buffer chart required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<DOMAINS[chart][0] or endpoints(v)[1]>DOMAINS[chart][1] or endpoints(nu)[0]<=0:raise ValueError('Actual reference/restoration chart domain,Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        restore=self.source_packet(chart,Z,v)
        raw=raw_restore_rows(c,restore,c.mpf(endpoints(self.restore.invP2)))
        z=IntervalTaylor.variable(c,Z,5);delta=c.mpf(self.pulse.delta)
        velocity=raw['velocity'];histories=raw['histories'];P=raw['absolute_pressure']
        rows=raw_pre_stress_rows(c,delta,z,velocity['theta'],velocity['axial'],histories,P)
        logR,radius=self.physical.radius(chart,v,restore,self.restore);sectors={}
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
            actual_upstream_current_reference_restore_source=restore,current_raw_five_history_rows=histories,
            current_absolute_pressure_ordinary_y_rows=P,current_source_three_component_velocity_rows=velocity,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,actual_positive_amplitude_source_enclosures=raw,
            full_current_correlated_E_centered_histories_and_original_cutoff_jets_retained=True,
            full_energy_baseline_cross_terms_and_absolute_P0_retained=True,
            derivative_coordinate_is_ordinary_logR_not_selector_phase=True,
            positive_amplitude_caps_enclose_original_log_source_not_define_it=True,
            only_radius_and_constant_Pstar_factored_from_full_stress=True,actual_full_stress_not_local_difference=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        pairs={'reference_restore':(('inner_reference',1),('axial_restore',0)),
            'restore_buffer':(('axial_restore',1),('restore_buffer',-7)),
            'buffer_patch':(('restore_buffer',-6),('actual_patch',1))}
        if name not in pairs:raise ValueError('Three actual reference/restore/patch tensor joins required')
        a,b=pairs[name];left=self.chart(a[0],Z,a[1],log_tau,theta,viscosity)
        right=(self.patch_tensor.chart if b[0]=='actual_patch' else self.chart)(b[0],Z,b[1],log_tau,theta,viscosity)
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Actual restore completed tensor layout differs')
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
            current_actual_restore_full_tensor_and_three_join_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.patch_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=24,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_restore_domains=DOMAINS,
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_core_bridge_switch_reshape_Rsh_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentRestoreBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_restore_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_restore_tensor_views'][name]=field.chart(*args)
        print('Current actual reference/restore full tensor: '+name,flush=True)
    result['current_actual_three_restore_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
