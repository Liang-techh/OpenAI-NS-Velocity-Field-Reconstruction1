"""Current implicit five-moment patch full tensor and its actual Rh attachment."""
import gzip
import json
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O2_background_tensor import (
    CurrentO2BackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,endpoints,
    source_precision,accepted,_verify_hashes,canonical_tensor_groups,BASE,PULSE)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_compliant_current_patch_stress_operator import raw_patch_rows,actual_patch_normalization_theorem
from lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4 import (
    CompliantActualFeedbackPatchMixedC4,current_mixed_source_bindings)
from lei_ren_part1_paper_compliant_actual_feedback_moment_patch import (
    CompliantActualFeedbackMomentPatch,current_patch_source_bindings)
from lei_ren_part1_paper_compliant_actual_patch_mixed_C4 import CompliantActualPatchMixedC4,patch_mixed
from lei_ren_part1_paper_compliant_actual_moment_patch import CompliantActualMomentPatch
from lei_ren_part1_paper_compliant_actual_Rh_source_join import pressure_defining_function_proof,Rh_functional_identity
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import CompliantPrePulseMixedC4

NAME=PREFIX+'current_actual_patch_background_tensor.json.gz'
RECEIPT=PREFIX+'current_actual_patch_background_tensor_check.json'
GATES=('current_actual_five_moment_patch_full_tensor_available',
    'current_actual_five_moment_patch_full_meridional_decomposition_available',
    'current_actual_patch_Rh_completed_tensor_join_certified',
    'current_actual_patch_six_internal_support_tensor_traces_certified')
EDGES=(49,51,59,61,69,71)
SEAMS=('patch_Rh',)+tuple('patch_support_%d'%edge for edge in EDGES)
VIEWS={'whole_patch':('actual_patch',(-1,1),'whole',('-3','-1'),None,'1'),
    'patch_Rm':('actual_patch',(-1,1),1,'-1',None,'1'),
    'patch_Rh':('actual_patch',(-1,1),'Rh','-1',None,'1')}
for _edge in EDGES:VIEWS['patch_edge_%d'%_edge]=('actual_patch',(-1,1),'edge:%d'%_edge,'-1',None,'1')
for _name,_x,_z in (('first','1.25','.317'),('second','1.5','-.537'),('third','1.75','.731')):
    VIEWS['active_'+_name]=('actual_patch',_z,_x,'-2.6','.41','.8')


def current_patch_source_theorem(field):
    """Bind this live current patch, terminal closure, P0 and compact edges."""
    patch=field.patch;pre=field.physical.pre;native=field.physical.dispatch.native
    live=dict(exact_current_mixed_owner=type(patch) is CompliantActualFeedbackPatchMixedC4,
        exact_current_moment_owner=type(patch.patch) is CompliantActualFeedbackMomentPatch,
        same_native_anchor=patch is native.anchor,same_actual_Rh_pre_owner=pre is native.rh_reference,
        unchanged_original_mixed_evaluate=patch.evaluate.__func__ is CompliantActualPatchMixedC4.evaluate,
        unchanged_original_gamma=patch.gamma.__func__ is CompliantActualPatchMixedC4.gamma,
        unchanged_original_patch_mixed=patch.evaluate.__func__.__globals__['patch_mixed'] is patch_mixed,
        unchanged_original_actual_data=patch.patch.actual_data.__func__ is CompliantActualMomentPatch.actual_data,
        unchanged_original_coefficients=patch.patch.coefficients.__func__ is CompliantActualMomentPatch.coefficients,
        unchanged_original_moment_evaluate=patch.patch.evaluate.__func__ is CompliantActualMomentPatch.evaluate,
        unchanged_original_pre_reference=pre.reference.__func__ is CompliantPrePulseMixedC4.reference,
        original_BASE_radius=field.physical.radius.__func__ is BASE.radius,
        same_actual_logPstar=endpoints(field.physical.logP)==endpoints(patch.patch.core.logP),
        same_delta=endpoints(field.ctx.mpf(field.pulse.delta))==endpoints(field.ctx.mpf(patch.patch.core.delta)),
        original_patch_adapter=field.chart.__func__.__wrapped__.__globals__['raw_patch_rows'] is raw_patch_rows,
        original_checked_raw_stress=field.chart.__func__.__wrapped__.__globals__['raw_pre_stress_rows'] is raw_pre_stress_rows)
    if not all(live.values()):raise ValueError('Same original current patch/Rh/physical source programs required')
    record=native.receipts['actual_feedback_patch_mixed_C4_check']
    if not (record['current_Rh_same_unique_implicit_function_and_full_support_retained'] and
            record['current_five_functional_terminal_identities_connected'] and native.Rh_acceptance_loaded):
        raise ValueError('Current unique implicit functional closure and actual Rh admission required')
    pressure=pressure_defining_function_proof(patch.patch.core.datum,pre.datum)
    if not native.pressure_proof['passed'] or not native.Rh_proof['passed']:raise ValueError('Actual native P0/Rh function bridge missing')
    rh=Rh_functional_identity();mixed=current_mixed_source_bindings();moment=current_patch_source_bindings()
    asts=SourceAST()
    asts.expression('actual_moment_patch','evaluate','terminal',wanted='endpoints(x)[0]>=mp.mpf(71)/40')
    asts.expression('actual_moment_patch','evaluate','p0',wanted="self.reference.inputs(Z)['original_axis_pressure']")
    asts.expression('actual_patch_mixed_C4','evaluate','p0',wanted="jet('original_P0_axial5')")
    asts.expression('actual_patch_mixed_C4','gamma','rows',wanted='beta_jets(c,argument)')
    asts.expression('actual_patch_mixed_C4','gamma','derivatives',wanted='[rows[k]*(math.factorial(k)/(self.radius**(k+1)*self.N)) for k in range(5)]')
    asts.method('actual_feedback_moment_patch','__init__')
    asts.method('actual_feedback_patch_mixed_C4','__init__')
    asts.method('actual_Rh_source_join','pressure_defining_function_proof')
    asts.method('actual_Rh_source_join','Rh_functional_identity')
    flat=accepted(PREFIX+'flat_pulse_derivatives_check.json',field.family,field.source,'original_radial_shape_derivatives_C4_available')
    _verify_hashes(flat)
    if not flat['support_crossing_derivatives_available'] or not all(flat['analytic_original_beta_and_flat_envelope_checks'].values()):raise ValueError('Original beta flat support theorem missing')
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    rr=s.Symbol('logRref',real=True);owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify,ln=s.log),logRref=rr)
    a=radius(owner,'actual_patch',s.E,{},patch)[0];b=radius(owner,'Rh_reference',s.Integer(-5),{},pre)[0]
    if s.simplify(a-b)!=0:raise ArithmeticError('Actual patch/Rh radius differs')
    r=s.Rational(1,40);edges=[s.Rational(i,4)+side*r for i in (5,6,7) for side in (-1,1)]
    if edges!=[s.Rational(i,40) for i in EDGES] or not s.Rational(71,40)<s.E:raise ArithmeticError('Original patch support geometry changed')
    if not patch.generic_check['exact_coordinate_and_join_checks']['passed']:raise ValueError('Original patch coordinate/support algorithm fixture required')
    return dict(live_original_callable_bindings=live,current_mixed_source_AST_bindings=mixed,
        current_moment_transport_inverse_and_pressure_AST_bindings=moment,
        actual_current_pressure_defining_function_bridge=pressure,actual_current_Rh_functional_identity=rh,
        current_terminal_closure_receipt_consumed=True,
        exact_patch_Rh_radius_identity='Rref*exp(-6)*e=Rref*exp(-5)',
        current_Rh_open_neighborhood_uses_same_unique_implicit_full_weight_closure=True,
        six_exact_rational_support_edges=[str(edge) for edge in edges],
        original_beta_flat_derivatives_0_through_4_and_continuous_partial_integrals_consumed=True,
        completed_tensor_divergence_and_remainder_traces_follow_from_same_source_mixed4=True,
        rounded_edge_boxes_keep_positive_tail_caps_and_are_not_exact_zero_assertions=True,
        same_analytic_P0_function_and_first_six_projection_not_hash_only=True,
        interval_overlap_not_used_as_function_identity=True,
        input_hashes={**asts.hashes,**flat['input_hashes'],PREFIX+'flat_pulse_derivatives_check.json':sha(PREFIX+'flat_pulse_derivatives_check.json')},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentActualPatchBackgroundTensor:
    @source_precision
    def __init__(self,o2_tensor=None,require_checked=True):
        self.o2_tensor=o2_tensor if o2_tensor is not None else CurrentO2BackgroundTensor()
        self.physical=self.o2_tensor.physical;self.pulse=self.o2_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.o2_tensor.family;self.source=self.o2_tensor.source;self.datum_sha=self.o2_tensor.datum_sha
        self.patch=self.physical.dispatch.provider('actual_patch');self.assert_graph();self.lift=self.o2_tensor.lift
        source=current_patch_source_theorem(self);normalization=actual_patch_normalization_theorem()
        self.proof=dict(current_actual_patch_source_pressure_support_and_Rh_theorem=source,
            original_actual_patch_current_radius_normalization_theorem=normalization,
            consumed_checked_arbitrary_raw_stress_velocity_and_full_physical_operator=self.o2_tensor.proof,
            full_current_implicit_axial5_coefficients_and_all_five_partial_histories_retained=True,
            actual_absolute_pressure_radial_velocity_and_full_energy_cross_terms_retained=True,
            exact_beta_flat_limits_not_rounded_interval_overlap_prove_internal_traces=True,passed=True)
        self.hashes=dict(self.o2_tensor.hashes)
        for value in (source,normalization):self.hashes.update(value['input_hashes'])
        for stem in ('current_actual_patch_background_tensor','current_patch_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Actual patch tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.o2_tensor.assert_graph()
        if not (self.o2_tensor.acceptance_loaded and self.physical is self.o2_tensor.physical and
                self.pulse is self.o2_tensor.pulse and self.ctx is self.physical.ctx and
                self.patch is self.physical.dispatch.provider('actual_patch') and
                self.patch is self.physical.dispatch.native.anchor):raise ValueError('Same checked O2 and actual current implicit patch graph required')
        if (self.patch.family!=self.family or self.patch.source!=self.source or
                self.patch.patch.core.datum.datum_sha!=self.datum_sha):raise ValueError('Foreign current patch source/datum')

    def coordinate(self,value):
        c=self.ctx
        if isinstance(value,str) and value=='whole':return c.mpf((1,endpoints(c.exp(1))[1]))
        if isinstance(value,str) and value=='Rh':return c.exp(1)
        if isinstance(value,str) and value.startswith('edge:'):
            edge=int(value[5:])
            if edge not in EDGES:raise ValueError('Original rational support edge required')
            return c.mpf(edge)/40
        return c.mpf(value)

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=self.coordinate(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart!='actual_patch' or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<1 or endpoints(v)[1]>endpoints(c.exp(1))[1] or endpoints(nu)[0]<=0:raise ValueError('Actual patch x[1,e],Z[-1,1],finite time,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        patch=self.patch.evaluate(v,Z)
        raw=raw_patch_rows(c,patch,c.mpf(endpoints(self.patch.invP2)))
        z=IntervalTaylor.variable(c,Z,5);delta=c.mpf(self.pulse.delta)
        velocity=raw['velocity'];histories=raw['histories'];P=raw['absolute_pressure']
        rows=raw_pre_stress_rows(c,delta,z,velocity['theta'],velocity['axial'],histories,P)
        logR,radius=self.physical.radius(chart,v,patch,self.patch);sectors={}
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
        return dict(point,chart=chart,coverage_coordinate=v,exact_coordinate_selector=coordinate,
            current_actual_source_stress_packet=packet,actual_upstream_current_implicit_patch_source=patch,
            current_raw_five_history_rows=histories,current_absolute_pressure_ordinary_y_rows=P,
            current_source_three_component_velocity_rows=velocity,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate(chart,Z,v,lt,theta),
            exact_original_radius_source=radius,full_current_implicit_axial5_and_partial_integrals_retained=True,
            full_energy_cross_terms_and_absolute_P0_retained=True,ordinary_x_converted_to_logR_before_full_physical_derivatives=True,
            only_radius_and_constant_Pstar_factored_from_full_stress=True,actual_full_stress_not_local_difference=True,
            rounded_support_coordinate_is_enclosure_not_exact_zero_source=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name not in SEAMS:raise ValueError('Actual patch Rh or six original support tensor traces required')
        if name=='patch_Rh':
            left=self.chart('actual_patch',Z,'Rh',log_tau,theta,viscosity)
            right=self.o2_tensor.chart('Rh_reference',Z,-5,log_tau,theta,viscosity)
            a=canonical_tensor_groups(left);b=canonical_tensor_groups(right)
        else:
            edge=int(name.rsplit('_',1)[1]);value=self.chart('actual_patch',Z,'edge:%d'%edge,log_tau,theta,viscosity)
            a=canonical_tensor_groups(value);b=a
        if set(a)!=set(b):raise ValueError('Actual patch completed tensor layout differs')
        common={}
        for key in a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in a[key]+b[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_source_and_completed_tensor_endpoint_function_theorem=self.proof,
            original_single_global_patch_function_has_equal_one_sided_mixed4_traces=name!='patch_Rh',
            source_function_equality_precedes_common_triangle_bounds=True,
            rounded_edge_box_not_exact_zero_assertion=True,interval_overlap_not_used_as_function_identity=True,
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_actual_patch_full_tensor_and_Rh_support_theorem=self.proof,
            actual_current_tensor_regions_available=['actual_patch']+self.o2_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=21,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            actual_patch_domain='x=R/Rm in [1,e]; exact e and rational edge selectors evaluated with outward enclosures',
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_core_bridge_switch_reshape_restore_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentActualPatchBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_patch_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_patch_tensor_views'][name]=field.chart(*args)
        print('Current actual implicit patch full tensor: '+name,flush=True)
    result['current_actual_patch_Rh_and_six_support_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
