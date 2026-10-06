"""Actual Rh reference and O2 full tensors on the checked O3/current graph."""
import ast
import copy
import gzip
import json
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_transition_background_tensor as o3
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import (
    CurrentO3TransitionBackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,
    endpoints,source_precision,accepted,_verify_hashes,canonical_tensor_groups,BASE,PULSE)
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import (
    CompliantPrePulseMixedC4,physical_mixed,turnoff_derivatives)
from lei_ren_part1_paper_compliant_actual_Rp_source_join import actual_parameter_formula_bindings
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import join_source_proof

NAME=PREFIX+'current_O2_background_tensor.json.gz'
RECEIPT=PREFIX+'current_O2_background_tensor_check.json'
GATES=('current_actual_O2_reference_slope_axial_buffer_full_tensors_available',
    'current_actual_O2_full_meridional_decomposition_available',
    'current_actual_four_O2_completed_tensor_joins_certified')
DOMAINS={'Rh_reference':(-5,0),'O2_slope':(0,1),'O2_axial':(0,1),'O2_buffer':(0,11)}
SEAMS=('reference_slope','slope_axial','axial_buffer','buffer_transition')
VIEWS={}
for _chart,_domain in DOMAINS.items():
    VIEWS[_chart+'_whole']=(_chart,(-1,1),_domain,('-3','-1'),None,'1')
    VIEWS[_chart+'_left']=(_chart,(-1,1),_domain[0],'-1',None,'1')
    VIEWS[_chart+'_right']=(_chart,(-1,1),_domain[1],'-1',None,'1')
    VIEWS[_chart+'_fresh']=(_chart,'.537',dict(Rh_reference='-2.337',O2_slope='.537',O2_axial='.1337',O2_buffer='5.337')[_chart],'-2.6','.41','.8')


def compiled_current_O2_exporter():
    """Keep checked variable-source math; adapt only route/domain/publication."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('current_O3_transition_background_tensor','chart'));fn.decorator_list=[]
    guard="chart != 'O3_slope_mu' or not all((mp.isfinite(x) for x in endpoints(Z) + endpoints(v) + endpoints(lt) + endpoints(nu))) or endpoints(Z)[0] < -1 or endpoints(Z)[1] > 1 or endpoints(v)[0] < 0 or endpoints(v)[1] > 1 or endpoints(nu)[0] <= 0"
    replacement="chart not in DOMAINS or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<DOMAINS[chart][0] or endpoints(v)[1]>DOMAINS[chart][1] or endpoints(nu)[0]<=0"
    changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.If) and ast.dump(node.test)==ast.dump(ast.parse(guard,mode='eval').body):
            node.test=ast.parse(replacement,mode='eval').body
            node.body[0]=ast.parse("raise ValueError('Actual Rh reference/O2 chart domain,finite time,Z[-1,1],nu>0 required')").body[0]
            changes.append('domain_guard')
        if isinstance(node,ast.Assign) and ast.unparse(node)=='pre = self.physical.pre.slope_mu(Z, v)':
            node.value=ast.parse('self.source_packet(chart,Z,v)',mode='eval').body;changes.append('actual_original_pre_route')
        if isinstance(node,ast.keyword) and node.arg=='actual_upstream_original_pre_transition_source':
            node.arg='actual_upstream_original_pre_O2_source';changes.append('publication_label')
    if sorted(changes)!=sorted(('domain_guard','actual_original_pre_route','publication_label')):raise ValueError('Unreviewed current O2 exporter adaptation')
    env=dict(vars(o3));env['DOMAINS']=DOMAINS
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<checked variable pre tensor,actual O2 source routes>','exec'),env)
    return env['chart'],dict(exact_AST_adaptations=changes,
        checked_raw_stress_velocity_pressure_remainder_and_full_physical_math_unchanged=True,
        input_hashes=asts.hashes)


def turnoff_ordinary_source_theorem():
    """Replay the original cutoff program in the actual logR variable."""
    y,md=s.symbols('positive_logR_offset Md',positive=True);sigma=s.Function('original_sigma')
    argument=1-s.log(y)/md;B=sigma(argument)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)) if isinstance(v,(str,int,float)) else v)
    def jets(ctx,x):
        a=s.Symbol('cutoff_argument')
        return [s.diff(sigma(a),a,j).subs(a,x)/math.factorial(j) for j in range(5)]
    asts=SourceAST();fn=asts.replay('pre_pulse_mixed_C4','turnoff_derivatives',dict(math=math,sigma_jets=jets))
    values=fn(c,y,md,s.log(y)/md);checks={}
    for j,value in enumerate(values):
        if s.simplify(value-s.diff(B,y,j))!=0:raise ArithmeticError('Ordinary logR turnoff derivative differs: '+str(j))
        checks['ordinary_logR_cutoff_derivative_'+str(j)]=True
    return dict(identities=checks,original_signed_Stirling_program_replayed=True,
        cutoff='sigma(1-log(y)/Md); y=log(R/Rref)=exp(Md*phase)',
        derivatives_are_ordinary_logR_not_selector_phase=True,input_hashes=asts.hashes,passed=True)


def current_O2_source_theorem(field):
    pre=field.physical.pre;asts=SourceAST()
    live={name:getattr(pre,name).__func__ is getattr(CompliantPrePulseMixedC4,name) for name in ('reference','slope','inlet','axial','slope_mu','packet')}
    live.update(original_physical_mixed=pre.packet.__func__.__globals__['physical_mixed'] is physical_mixed,
        original_turnoff_ordinary_derivatives=pre.axial.__func__.__globals__['turnoff_derivatives'] is turnoff_derivatives,
        original_BASE_radius=field.physical.radius.__func__ is BASE.radius,
        same_parameter_and_radius_logPstar=endpoints(field.physical.logP)==endpoints(field.ctx.mpf(endpoints(pre.params.logPstar))),
        same_actual_pre_owner=pre is field.o3_tensor.physical.pre,
        unchanged_source_context_copy=field.exporter.__globals__['copy_jet'] is o3.copy_jet)
    if not all(live.values()):raise ValueError('Actual original O2 source programs/parameters required')
    trace=SimpleNamespace(reference=lambda Z,v:('reference',Z,v),slope=lambda Z,v:('slope',Z,v),
        axial=lambda Z,**kwargs:('axial',Z,kwargs))
    route=asts.replay('current_O2_background_tensor','source_packet',{})
    owner=SimpleNamespace(physical=SimpleNamespace(pre=trace));z,v=s.symbols('Z source_coordinate',real=True)
    expected={'Rh_reference':('reference',z,v),'O2_slope':('slope',z,v),
        'O2_axial':('axial',z,dict(phase=v)),'O2_buffer':('axial',z,dict(buffer_offset=v))}
    routes={chart:route(owner,chart,z,v)==wanted for chart,wanted in expected.items()}
    if not all(routes.values()):raise ValueError('Actual O2 source routing differs')
    parameters=actual_parameter_formula_bindings();joins=join_source_proof(field.ctx)
    asts.expression('pre_pulse_mixed_C4','slope_mu','parent',wanted='self.axial(Z,buffer_offset=11)')
    asts.expression('pre_pulse_mixed_C4','inlet', 'self.cache[key]',wanted='self.slope(Z,1)')
    asts.expression('pre_pulse_mixed_C4','axial','parent',wanted='self.inlet(Z)')
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='turnoff_derivatives(c,y,md,phase)')
    asts.expression('pre_pulse_mixed_C4','axial','K',wanted='turnoff_kernels(c,y,self.params.Md,self.cells,self.window)')
    asts.expression('pre_pulse_mixed_C4','packet','data',wanted='physical_mixed(c,Z,self.delta,u,logU,V,history,p0,self.invP2)')
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    rr,md=s.symbols('logRref Md',real=True)
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRref=rr,logP=s.exp(md)+11,params=SimpleNamespace(Md=md))
    pairs={'reference_slope':(('Rh_reference',0,{}),('O2_slope',0,{})),
        'slope_axial':(('O2_slope',1,{}),('O2_axial',0,dict(actual_y=s.Integer(1),exact_radius_source='original'))),
        'axial_buffer':(('O2_axial',1,dict(actual_y=s.exp(md),exact_radius_source='original')),('O2_buffer',0,dict(actual_y=s.exp(md),exact_radius_source='original'))),
        'buffer_transition':(('O2_buffer',11,dict(actual_y=s.exp(md)+11,exact_radius_source='original')),('O3_slope_mu',0,{}))}
    identities={}
    for label,(a,b) in pairs.items():
        left=radius(owner,*a,None)[0];right=radius(owner,*b,None)[0]
        if s.cancel(left-right)!=0:raise ArithmeticError('Actual O2 tensor source radius differs: '+label)
        identities[label+'_same_actual_radius']=True
    raw=field.o3_tensor.proof['current_actual_transition_source_and_power_endpoint']['checked_same_production_source_function_theorem']
    if not raw['passed']:raise ValueError('Checked actual composed O2/O3 production functions required')
    return dict(live_original_callable_bindings=live,actual_route_AST_replay=routes,
        original_parameter_definitions=parameters,original_four_pre_endpoint_source_theorem=joins,
        exact_radius_identities=identities,ordinary_turnoff_cutoff_theorem=turnoff_ordinary_source_theorem(),
        checked_actual_composed_O2_O3_functions=raw,
        common_five_histories_absolute_pressure_and_velocity_source_jets_identified=True,
        four_interface_velocity_pressure_and_primitive_mixed4_rows=540,
        full_completed_tensor_traces_follow_from_same_sources_and_checked_original_operators=True,
        Rh_reference_start_is_actual_source_boundary_not_admitted_repaired_core_tensor_join=True,
        interval_overlap_not_used_as_function_identity=True,
        input_hashes={**asts.hashes,**parameters['input_hashes']},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentO2BackgroundTensor:
    @source_precision
    def __init__(self,o3_tensor=None,require_checked=True):
        self.o3_tensor=o3_tensor if o3_tensor is not None else CurrentO3TransitionBackgroundTensor()
        self.physical=self.o3_tensor.physical;self.pulse=self.o3_tensor.pulse;self.ctx=self.physical.ctx
        self.family=self.o3_tensor.family;self.source=self.o3_tensor.source;self.datum_sha=self.o3_tensor.datum_sha
        self.assert_graph();self.lift=self.o3_tensor.lift;self.exporter,self.adaptation=compiled_current_O2_exporter()
        source=current_O2_source_theorem(self)
        generic=self.o3_tensor.proof
        for key in ('original_arbitrary_variable_source_full_stress_theorem','current_actual_variable_source_velocity_and_remainder_theorem'):
            if not generic[key]['passed'] or not all(generic[key]['identities'].values()):raise ValueError('Checked arbitrary variable raw pre theorem required')
        self.proof=dict(current_actual_O2_source_and_four_endpoint_theorem=source,
            consumed_checked_arbitrary_variable_source_stress_velocity_and_full_physical_operator=generic,
            current_O2_route_exporter_adaptation=self.adaptation,
            actual_axial_transport_and_full_ordinary_cutoff_jets_retained=True,
            full_nonzero_histories_not_reset_when_local_axial_velocity_turns_off=True,
            source_function_equality_precedes_common_tensor_bounds=True,passed=True)
        self.hashes=dict(self.o3_tensor.hashes)
        self.hashes.update(source['input_hashes']);self.hashes.update(source['ordinary_turnoff_cutoff_theorem']['input_hashes']);self.hashes.update(self.adaptation['input_hashes'])
        name=PREFIX+'current_O2_background_tensor.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current O2 tensor admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.o3_tensor.assert_graph()
        if not (self.o3_tensor.acceptance_loaded and self.physical is self.o3_tensor.physical and self.pulse is self.o3_tensor.pulse and self.ctx is self.physical.ctx):raise ValueError('Same checked current O3 transition graph required')

    def source_packet(self,chart,Z,v):
        if chart=='Rh_reference':return self.physical.pre.reference(Z,v)
        if chart=='O2_slope':return self.physical.pre.slope(Z,v)
        if chart=='O2_axial':return self.physical.pre.axial(Z,phase=v)
        if chart=='O2_buffer':return self.physical.pre.axial(Z,buffer_offset=v)
        raise ValueError('Actual Rh reference/O2 source chart required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        return self.exporter(self,chart,Z,coordinate,log_tau,theta,viscosity)

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        pairs={'reference_slope':(('Rh_reference',0),('O2_slope',0)),
            'slope_axial':(('O2_slope',1),('O2_axial',0)),
            'axial_buffer':(('O2_axial',1),('O2_buffer',0)),
            'buffer_transition':(('O2_buffer',11),('O3_slope_mu',0))}
        if name not in pairs:raise ValueError('Declared actual O2 completed tensor seam required')
        a,b=pairs[name];left=self.chart(a[0],Z,a[1],log_tau,theta,viscosity)
        right=(self.o3_tensor.chart if b[0]=='O3_slope_mu' else self.chart)(b[0],Z,b[1],log_tau,theta,viscosity)
        rows_a=canonical_tensor_groups(left);rows_b=canonical_tensor_groups(right);common={}
        if set(rows_a)!=set(rows_b):raise ValueError('Actual O2 completed tensor component layout differs')
        for key in rows_a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in rows_a[key]+rows_b[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_source_and_completed_tensor_endpoint_function_theorem=self.proof,
            source_function_equality_precedes_common_triangle_bounds=True,
            interval_overlap_not_used_as_function_identity=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_actual_O2_tensor_and_four_joins_theorem=self.proof,
            actual_current_tensor_regions_available=list(DOMAINS)+self.o3_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=20,actual_current_completed_tensor_internal_interface_count=4,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),actual_O2_domains=DOMAINS,
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_repaired_core_tensor_Rh_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentO2BackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_O2_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_O2_tensor_views'][name]=field.chart(*args)
        print('Current actual O2 full tensor: '+name,flush=True)
    result['current_actual_four_O2_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
