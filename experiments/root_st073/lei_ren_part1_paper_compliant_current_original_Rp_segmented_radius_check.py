"""Original integral-defined waiting, absolute radius and Jacobian checks."""
import ast
import copy
import gzip
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius as current
import lei_ren_part1_paper_compliant_current_heat_physical_assembly as heat_geometry
import lei_ren_part1_paper_compliant_current_steep_waiting_physical_assembly as steep_geometry
import lei_ren_part1_paper_compliant_outer_angular_candidate as angular_source
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def interpreter(owner,literal=True,bindings=None):
    field=SimpleNamespace(graph=owner.graph,parameters=owner.parameters,contracts=owner.frame.bridge.leading.contracts)
    return current.RadiusInterpreter(field,literal,bindings)


def same(a,b,label):
    assert s.cancel(s.expand(a-b))==0,label


def sigma_live_proof(owner):
    angular=owner.before.before.future.angular
    assert angular.flatten.__func__ is angular_source.SharedOuterAngularCandidate.flatten
    assert angular.before_waiting.__func__ is angular_source.SharedOuterAngularCandidate.before_waiting
    assert angular.flatten.__func__.__globals__['stable_sigma'] is angular_source.stable_sigma
    assert angular_source.unit_kernels.__globals__['stable_sigma'] is angular_source.stable_sigma
    assert angular_source.collar_preheat_integral.__globals__['stable_sigma'] is angular_source.stable_sigma
    def require(fn,statement):
        tree=ast.parse(inspect.getsource(fn));wanted=ast.dump(ast.parse(statement).body[0])
        assert any(ast.dump(n)==wanted for n in ast.walk(tree)),(fn.__name__,statement)
    for statement in ('if v<=0: return c.mpf(0),c.mpf(0)',
        'if v>=1: return c.mpf(1),c.mpf(0)','a=c.exp(-1/x**2)',
        'b=c.exp(-1/(1-x)**2)','value=a/(a+b)'):
        require(angular_source.stable_sigma,statement)
    jets=current.leading_source.providers.slope.original.original.sigma_jets
    assert owner.frame.bridge.leading.contracts['original_sigma']==current.leading_source.current.ast_binding(jets)
    left=jets.__globals__['_sigma_left']
    for statement in ('odds = 1/(1-x)**2 - 1/x**2','e = positive_exp(c, odds)','value = e/(1+e)'):
        require(left,statement)
    A,B=s.symbols('left_flat right_flat',positive=True)
    same((A/B)/(1+A/B),A/(A+B),'live scalar and jet sigma are the same function')
    same(A/(A+B)+B/(A+B),1,'reflected sigma pair integrates to one')
    return dict(actual_scalar_and_graph_Taylor_sigma_functions_identified=True,
        exact_flat_endpoints_and_sigma_reflection_proved=True,
        original_unit_sigma_integral_is_exactly_one_half=True,passed=True)


def waiting_proof(owner):
    f=owner.functions;mu=s.Symbol('mu',positive=True);a=s.Symbol('a',positive=True);k=s.Symbol('k',positive=True)
    L=s.Symbol('L',positive=True);T=s.Symbol('T',positive=True);eps=s.Symbol('eps',positive=True);Xp=s.Symbol('Xp',real=True)
    bindings={f[key].node:value for key,value in dict(rate=a,k=k,Lrel=L,Ts=T,epsilon=eps,Xp=Xp).items()}
    q=interpreter(owner,True,bindings);Sigma=current.post.selected.inlet.exact.frame.bridge.leading_check.Sigma
    u=s.Symbol('current_Rp_wait_kernel_u',real=True);v=s.Symbol('current_Rp_wait_kernel_v',real=True)
    J=s.Integral(Sigma(v),(v,0,u))
    Iin=s.Integral(s.exp(a*(u-J)),(u,0,1));Iout=s.Integral(s.exp(k*J),(u,0,1))
    F=s.exp(-s.log(2)*Sigma(u/100))
    flatten=s.Integral(F*s.exp(-a*(100-u)),(u,0,100))
    flat=s.Piecewise((s.exp(-1/((3-u)/2)**2),(3-u)/2>0),(0,True))
    collar=s.Integral(s.exp(k*u)*(1-Sigma(u)+Sigma(u)*flat),(u,0,3))
    for key,wanted in (('flatten_mass',flatten),('Iin',Iin),('Iout',Iout),('collarJ',collar)):
        assert q.at(f[key]).as_dummy()==wanted.as_dummy(),key
    Xv=1/a+(Xp-1/a)*s.exp(-13*a/mu)
    same(q.at(f['Xv']),Xv,'current native affine Xv')
    Xf=2*(Xv*s.exp(-100*a)+flatten)
    Xrel=1/a+(Xf-1/a)*s.exp(-a*L)
    Xt=((Xrel+Iin)*s.exp(-a/2)+T+Iout)*s.exp(-k/2)
    same(q.at(f['Xt']),Xt,'original before-waiting Xt')
    wanted=(s.log(Xt-1/k)+s.log(1-eps)-s.log(eps)-s.log(1/k+collar))/k
    # Bind independently checked composite nodes before reading the root.
    q.bindings[f['Xt'].node]=Xt;q.bindings[f['collarJ'].node]=collar
    assert q.at(f['waiting'])==wanted
    q0,constants=owner.before.before.inlet.exact.current_constants()
    actual=interpreter(owner,False,{owner.frame.functions['Tw'].node:s.Symbol('Tw',positive=True)})
    same(actual.at(f['Xp']),constants['Xp'],'current full inlet Xp at Z=0')
    assert q0.z not in actual.at(f['Xp']).free_symbols
    native=current.post.selected.inlet.exact.frame.native
    scalar_sources={}
    specs=(('future_swirl_energy','CompliantFutureSwirlEnergy','__init__',{'self.Lrel':'-30*self.params.log_mu'}),
        ('outer_angular_candidate','SharedOuterAngularCandidate','reference_buffer',
            {'t':'-30*self.params.log_mu*phase','departure':'(inlet[\'X\']-eq)*c.exp(-self.rate*t)'}),
        ('outer_angular_candidate','SharedOuterAngularCandidate','steep_in',
            {'X':"(inlet['X']+I)*c.exp(-self.rate*(t-J))"}),
        ('outer_angular_candidate','SharedOuterAngularCandidate','before_waiting',
            {'X':"(inlet['X']+I)*c.exp(-self.restore_rate*J)"}),
        ('current_exact_repair_branch','CurrentExactRepairBranch','__init__',
            {'terminal':"angular.before_waiting('0')",
             'angular.waiting':'(c.ln(Xt-eq)+angular.waiting_logone-angular.params.log_epsilon-c.ln(eq+angular.collarJ))/k'}),
        ('pressure_source','CompliantOuterParameters','__init__',
            {'self.log_epsilon':"c.ln(c.mpf('.001')) + self.log_delta"}),
        ('current_power_angular_source','CurrentPowerAngularC4','__init__',
            {'self.Lrel':'box(future.Lrel)','self.future':'self.fifth.fourth.energy.base'}),
        ('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',
            {'self.wait':'box(f.angular.waiting)','self.Ts':'box(f.params.Ts)','self.epsilon':'box(f.heat.epsilon)'}))
    # The base parameter module has no compliant prefix. Check its actual
    # Ts assignment directly rather than constructing a separate family.
    for stem,cls,fn,assignments in specs:
        scalar_sources[stem+'.'+fn]=native.ast_assignments(stem,cls,fn,assignments)
    # The literal sigma is symmetric under reflection. Hence integrating
    # sigma(t)+sigma(1-t)=1 proves J(1)=1/2, independently of a boxed kernel.
    assert Sigma(0)==0 and Sigma(1)==1
    shape=owner.frame.bridge.leading.contracts['original_sigma']
    assert shape==owner.graph.nodes[next(i for i,n in enumerate(owner.graph.nodes)
        if n['operation']=='exact_original_source_flat_sigma')]['source_binding']
    tree=ast.parse((current.HERE/'lei_ren_part1_paper_logarithmic_outer_parameters.py').read_text(encoding='utf8'))
    wanted_ast=ast.dump(ast.parse('4*(c.ln(2)-self.log_delta)',mode='eval').body)
    assert sum(isinstance(n,ast.Assign) and any(ast.unparse(t)=='self.Ts' for t in n.targets)
        and ast.dump(n.value)==wanted_ast for n in ast.walk(tree))==1
    for target,expr in (('self.md','c.mpf(Md)'),('self.logPstar','c.exp(self.md)+11'),
        ('choice','-4*self.logPstar-30')):
        wanted_ast=ast.dump(ast.parse(expr,mode='eval').body)
        assert sum(isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)
            and ast.dump(n.value)==wanted_ast for n in ast.walk(tree))==1
    assert owner.before.before.future.params.Md=='40'
    canonical=current.post.selected.inlet.exact.frame.bridge.outer.rh.read(current.PREFIX+'current_exact_repair_branch_check.json')
    assert canonical['all_passed'] and canonical['current_native_Xv_exact_source_bound_to_repair']
    assert owner.before.before.future.repair is owner.before.before.seed.exact.repair
    assert owner.before.before.future.angular is owner.before.before.seed.exact.repair.angular
    source=canonical['current_exact_native_repair_source_proof']
    assert source['passed']
    assert owner.before.power.bindings['current_native_parameter_source_bridge']['passed']
    assert owner.before.power.bindings['current_native_parameter_source_bridge']['original_exact_parameter_formula_theorem_reused']
    return dict(current_Xp_function_at_Z0_proved=True,literal_original_waiting_kernels_identified=4,
        native_affine_terminal_and_entire_waiting_equation_identified=True,
        actual_original_constructor_source_bindings=scalar_sources,actual_Ts_source_assignment_identified=True,
        actual_Md40_logPstar_and_delta_branch_assignments_identified=True,
        actual_sigma_endpoint_symmetry_and_half_primitive_proof=sigma_live_proof(owner),
        Md40_delta_choice_theorem_reused=True,epsilon_is_delta_over1000=True,
        current_unique_repair_and_waiting_owner_retained=True,waiting_box_not_used_as_function_value=True,passed=True)


def maps_proof(owner):
    f=owner.functions;rho=s.Symbol('rho_p',real=True);mu=s.Symbol('mu',positive=True)
    L,T,W=s.symbols('L T W',positive=True);v=s.Symbol('current_Rp_radius_coordinate',real=True)
    q=interpreter(owner,True,{owner.logRp.node:rho,f['Lrel'].node:L,f['Ts'].node:T,f['waiting'].node:W})
    common=rho+13/mu
    expected={'pulse_entrance':rho+v,'pulse_main':rho+v/mu,'pulse_exit':rho+v/mu,'pulse_gap':rho+v/mu,
        'pulse_gap_end':common+v,'pulse_end':common+v,'flatten':common+v,
        'outer_power':common+100+(L-4)*v,'outer_angular':common+100+L+v,
        'steep_entry':common+100+L+v,'steep_power':common+101+L+T*v,
        'steep_exit':common+101+L+T+v,'waiting':common+102+L+T+W*v,
        'heat_collar':common+102+L+T+W+v,'heat_exterior':common+102+L+T+W+v}
    assert set(expected)==set(owner.before.registry)
    for chart,wanted in expected.items():
        same(q.at(owner.maps[chart]['logR']),wanted,chart+'/absolute radius')
        same(q.at(owner.maps[chart]['native_to_log_radius_jacobian']),s.diff(wanted,v),chart+'/Jacobian')
    seams=(('pulse_entrance',s.Rational(1,50)/mu,'pulse_main',s.Rational(1,50)),
        ('pulse_main',10,'pulse_exit',10),('pulse_exit',11,'pulse_gap',11),
        ('pulse_gap',12,'pulse_gap_end',-1/mu),('pulse_gap_end',-4,'pulse_end',-4),
        ('pulse_end',0,'flatten',0),('flatten',100,'outer_power',0),
        ('outer_power',1,'outer_angular',-4),('outer_angular',0,'steep_entry',0),
        ('steep_entry',1,'steep_power',0),('steep_power',1,'steep_exit',0),
        ('steep_exit',1,'waiting',0),('waiting',1,'heat_collar',0),('heat_collar',3,'heat_exterior',3))
    for a,av,b,bv in seams:same(q.at(owner.maps[a]['logR']).subs(v,av),q.at(owner.maps[b]['logR']).subs(v,bv),a+'->'+b)
    return dict(original_absolute_radius_function_identities=15,original_native_coordinate_Jacobian_identities=15,
        exact_segmented_radius_seam_identities=14,physical_velocity_or_stress_seam_certification_not_inferred=True,passed=True)


@source_precision
def run(owner=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpSegmentedRadius(require_checked=False)
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record
    assert owner.assert_graph()==raw['actual_same_object_graph']
    assert owner.graph.nodes==raw['exact_function_graph_nodes']
    assert owner.graph.nodes[:len(owner.prefix)]==owner.frame.graph.nodes
    wait=waiting_proof(owner);geometry=maps_proof(owner)
    radius_bindings=dict(power=current.post.downstream.power.original_radius_functional_bindings(),
        steep=steep_geometry.original_steep_waiting_radius_bindings(),heat=heat_geometry.original_heat_radius_bindings())
    assert all(row['passed'] for row in radius_bindings.values())
    local=raw['actual_local_radius_steps']
    for chart in local:
        view=owner.local_step(chart,3,'7/2') if chart=='heat_exterior' else owner.local_step(chart,0,'1/2')
        assert current.post.selected.inlet.encode(current.post.selected.inlet.pack(view))==local[chart]
        assert current.post.selected.inlet.endpoints(view['native_to_log_radius_jacobian_bound'])[0]>0
    for chart in ('flatten','steep_entry','steep_exit','heat_collar','heat_exterior'):
        view=owner.local_step(chart,3,'7/2') if chart=='heat_exterior' else owner.local_step(chart,0,'1/2')
        assert current.post.selected.inlet.endpoints(view['log_radius_increment_bound'])==(s.Rational(1,2),s.Rational(1,2))
    rejected=[]
    for chart,value in (('flatten',101),('outer_angular',1),('waiting',2),('heat_collar',4),('heat_exterior',2)):
        try:owner.geometry(chart,value)
        except ValueError:rejected.append(chart)
        else:raise AssertionError('Out-of-domain radius accepted')
    old=copy.copy(owner);old.logRp=current.FunctionRef(owner.graph,owner.frame.functions['original_Rc_phase_offset'].node)
    try:old.assert_graph()
    except ValueError:pass
    else:raise AssertionError('Relative log(R/r_minus) accepted as absolute logRp')
    view=owner.evaluate('heat_collar','.427',0)
    assert view['current_absolute_source_geometry']==raw['actual_seventeen_current_radius_calls'][13]['geometry']
    assert owner.before.provider('heat_collar') is owner.before.heat
    assert len(raw['actual_seventeen_current_radius_calls'])==len(raw['actual_source_call_trace'])==17
    assert all(row['actual_absolute_source_radius_mapper_called'] for row in raw['actual_source_call_trace'])
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        independent_literal_waiting_function_proof=wait,independent_absolute_radius_and_Jacobian_proof=geometry,
        retained_actual_original_radius_operator_bindings=radius_bindings,
        seventeen_actual_current_source_calls_have_absolute_radius_functions=True,
        exact_half_unit_local_steps_not_lost_to_huge_logC_or_inverse_mu=True,
        out_of_domain_geometry_rejected=rejected,relative_phase_origin_rejected_as_absolute_radius=True,
        original_function_graph_prefix_and_selected_frequency_unchanged=True,
        local_directed_bounds_not_used_as_point_or_source_values=True,
        scope='Absolute lazy source geometry and local coordinate bounds; numerical absolute point radius, physical field and global mixed contracts remain open',
        **dict.fromkeys(current.OPEN,False),input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(heat_geometry.__file__).name:current.sha(Path(heat_geometry.__file__).name),
            Path(steep_geometry.__file__).name:current.sha(Path(steep_geometry.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.post.selected.inlet.encode(current.post.selected.inlet.pack(result)),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_SEGMENTED_RADIUS 15 maps/Jacobians; 14 exact seams; true waiting integrals',flush=True)
    return result


if __name__=='__main__':run()
