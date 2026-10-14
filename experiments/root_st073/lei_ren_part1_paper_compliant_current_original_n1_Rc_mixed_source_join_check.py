"""Actual flat O3 source, nonzero Rc inlet and finite own mixed join."""
from pathlib import Path
import json
import time
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_n1_Rc_mixed_source_join as current
import lei_ren_part1_paper_compliant_current_original_n1_combined_outer_mixed_source_functions_check as combined_check

compact_check=combined_check.compact_check
outer_check=combined_check.outer_check
zero,bind=compact_check.zero,compact_check.bind
Algebra,rows,encode=current.outer.Algebra,current.outer.rows,current.outer.encode


def bind_input(env,owner,ids,label,Z,Y=None):
    """An arbitrary smooth input consistent with exact constant rows."""
    constant=next((n for n,node in enumerate(ids) if node in (0,1)),None)
    if constant is None:body=sy.Function(label)(*(tuple([Y]) if Y is not None else ()),Z)
    else:
        assert all(node==0 for node in ids[constant+1:])
        body=sy.Integer(ids[constant])*Z**constant/sy.factorial(constant)
        for n in range(constant):
            coefficient=sy.Function(label+'_coefficient'+str(n))(Y) if Y is not None else sy.Symbol(label+'_coefficient'+str(n),real=True)
            body+=coefficient*Z**n
    bind(env,owner.jet(ids),body,Z);return body


class RcAlgebra(compact_check.CompactAlgebra):
    def __init__(self,owner,bindings,flat=False):
        super().__init__(owner.c,bindings);self.rc_owner=owner;self.flat=flat
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation'];owner=self.rc_owner
        if op=='original_lazy_flat_branch' and row['eta']==owner.leading['eta']:
            assert self.flat and row['Delta']==owner.leading['chart']['roots']['Delta'][0]
            assert row['active_body_not_evaluated_on_flat_branch'] and row['flat_value']==0
            assert self.g.nodes[row['flat_predicate']]==dict(operation='exact_real_comparison',operator='>=',left=row['Delta'],right=row['eta'])
            self.flat_nodes.add(i);result=sy.Integer(0)
        elif op=='analytic_unary' and row['name']=='expm1':result=sy.exp(self.at(row['argument']))-1
        elif op=='definite_integral':
            lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:result=sy.Integer(0)
            else:
                assert row['endpoint_independent_of_Z']
                t=sy.Symbol(row['variable'],real=True)
                body=sy.expand_mul(sy.expand_power_exp(self.at(row['integrand'])))
                terms=[]
                for term in sy.Add.make_args(body):
                    coefficient,integrand=term.as_independent(t,as_Add=False)
                    terms.append(coefficient*sy.Integral(integrand,(t,lo,hi)))
                result=sy.Add(*terms)
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            symbol=sy.Symbol(self.g.nodes[row['variable']]['name'],real=True);value=self.at(row['value'])
            env={key:body.subs(symbol,value) if hasattr(body,'subs') else body for key,body in self.bindings.items()}
            env[row['variable']]=value
            child=RcAlgebra(owner,env,self.flat);result=child.at(row['expression']);self.flat_nodes.update(child.flat_nodes)
        else:return super().at(i)
        self.cache[i]=result;return result


def source_contract(owner):
    owner.assert_graph()
    for module in (current.previous,current.leading_source,current.compact_source,current.outer.intermediate):
        combined_check.admitted(owner,module)
    c=owner.c.c;d=owner.leading;g=owner.g
    assert c['incoming_background']==d['background_endpoint']['history']
    assert c['incoming_correction']==d['endpoint']['correction']
    assert c['incoming_integrals']==d['endpoint']['integrals'] and c['initial']==d['initial']
    assert c['Rc']==d['endpoint']['radius'] and c['mu']==d['mu'] and c['N']==d['N']
    assert owner.intermediate['coordinate']==owner.coordinate.node and owner.intermediate['radius']==owner.radius.node
    assert owner.contract==current.previous.CONTRACT
    assert all(len(row)==5 for row in (owner.profiles,owner.histories))
    assert all(len(rows(jet))==6 for packet in owner.profiles+owner.histories for jet in packet.values())
    assert all(len(rows(jet))==(5 if key.startswith('Q') else 6) for key,jet in owner.physical.items())
    endpoint_rows=0
    def fixed_endpoint(jet,variable,value):
        nonlocal endpoint_rows
        for ref in rows(jet):
            row=g.nodes[ref.node]
            assert row['operation']=='function_substitution' and row['variable']==variable
            assert row['Z_independent_substitution']
            assert g.nodes[row['value']]==dict(operation='exact_rational',numerator=value,denominator=1)
            endpoint_rows+=1
    for side,packet,variable,value in (('left',owner.leading_left,owner.coordinate.node,2),('right',owner.leading_right,owner.c.x.node,1)):
        for group in ('profiles','histories'):
            for data in packet[group]:
                for jet in data.values():fixed_endpoint(jet,variable,value)
        for group in ('physical','density','higher_density'):
            for jet in packet[group].values():fixed_endpoint(jet,variable,value)
    for side,variable,value in (('left',owner.coordinate.node,2),('right',owner.c.phase.node,0)):
        for spec in owner.contract:
            jets=owner.join[side][spec['group']][spec['field']]
            assert len(jets)==spec['physical_y_order']+1
            for jet in jets:
                assert len(rows(jet))==spec['ordinary_Z_order']+1
                fixed_endpoint(jet,variable,value)
    integral_rows=0
    for data,coordinate,lower in ((owner.leading['chart']['local_integrals'],owner.coordinate.node,0),
        (c['pressure_local'],owner.c.x.node,1)):
        for key,ids in data.items():
            assert len(ids)==4
            for n,node in enumerate(ids):
                row=g.nodes[node]
                assert row['operation']=='definite_integral' and row['endpoint_independent_of_Z']
                assert row['upper']==coordinate and row['lower']==lower
                assert row['ordinary_slow_Z_derivative_order']==n
                assert row['physical_Jacobian_applied_once'] and row['physical_radius_power']==(1 if key=='S0' else 2)
                integral_rows+=1
    # Re-run only the fixed canonical parameter admission, not numerical
    # field oracles. Its actual eta log is restored from the pinned source.
    proof=g.nodes[d['eta']]['canonical_parameter_proof'];assert proof==d['mu_proof']['eta_proof']
    assert d['mu_proof']['source_family']==owner.family
    fn=current.leading_source.original_outer.o3.source_mu_admission
    binding=current.parent.active.current.ast_binding(fn)
    assert binding==d['mu_proof']['original_admission_binding'] and owner.hashes[binding['module']]==binding['module_sha256']
    ctx=MPIntervalContext();ctx.dps=500
    eta_log=current.leading_source.previous.previous.leading.restore_value(ctx,proof['original_canonical_eta_log'])
    logmu,Tw,margin=fn(ctx,eta_log)
    assert current.leading_source.original_outer.o3.ep(margin)[0]>=0
    assert current.leading_source.original_outer.o3.ep(Tw)[0]>current.leading_source.original_outer.o3.ep(2+ctx.ln(2))[1]
    mu_row=g.nodes[d['mu']];assert mu_row['operation']=='analytic_unary' and mu_row['name']=='exp'
    actual_log=outer_check.SourceAlgebra(g,{}).at(mu_row['argument'])
    zero(actual_log-(sy.log(sy.Rational(1,1000))-4*(sy.exp(40)+11)),'actual original logarithmic mu recipe')
    pins=compact_check.contract(owner.c)
    for fn in (current.leading_source.transport_geometry,current.leading_source.pressure_geometry,
        current.leading_source.engine.generic_roots,current.leading_source.engine.loop_functions,
        current.leading_source.engine.finite_N_functions):
        pin=current.parent.active.current.ast_binding(fn);assert owner.hashes[pin['module']]==pin['module_sha256'];pins[fn.__name__]=pin
    return dict(passed=True,actual_canonical_eta_log_and_mu_admission_bound=True,
        source_family=owner.family,original_admission_binding=binding,canonical_power_log_margin=str(margin),
        mu_and_eta_bounds_not_field_evaluations=True,actual_Rc_background_correction_integral_initial_node_maps_identified=True,
        actual_fixed_Z_independent_endpoint_substitution_rows=endpoint_rows,
        actual_pressure_integral_C3_commutation_source_rows=integral_rows,
        higher_Z_trace_rows_are_derivatives_of_directly_proved_whole_function_identities=True,
        ordinary_source_and_Duhamel_pressure_bindings=pins)


def power_context(owner):
    g=owner.g;d=owner.leading;Z=sy.Symbol('Z',real=True)
    Y=sy.Symbol(g.nodes[owner.coordinate.node]['name'],real=True);mu=sy.Symbol('mu',positive=True)
    P=sy.Symbol('actual_Pstar',positive=True);N=sy.Symbol('actual_N',positive=True,integer=True)
    delta=sy.Symbol('delta',positive=True);R0=sy.Symbol('actual_O3_power_R0',positive=True)
    env={g.symbol('Z').node:Z,owner.coordinate.node:Y,owner.mu.node:mu,owner.Pstar.node:P,
        owner.delta.node:delta,d['N']:N}
    bg=d['background'];bp=g.add(g.constant('1/2'),owner.mu)
    power=g.unary('exp',g.neg(g.mul(bp,owner.coordinate)))
    amplitude=[]
    for n,node in enumerate(bg['fields']['E']):
        row=g.nodes[node];assert row['operation']=='product' and power.node in row['arguments']
        factors=[i for i in row['arguments'] if i!=power.node]
        assert factors and all(not outer_check.free_query_dependency(owner.own,i,owner.coordinate.node) for i in factors)
        coefficient=g.mul(*(owner.ref(i) for i in factors));env[coefficient.node]=sy.diff(sy.Function('actual_power_inlet_E')(Z),Z,n)
        amplitude.append(coefficient.node)
    offset=outer_check.SourceAlgebra(g,env).at(bg['log_offset'])
    zero(sy.diff(offset,Y)-1,'actual O3 physical log radius Jacobian')
    env[owner.radius.node]=R0*sy.exp(Y)
    return dict(env=env,Z=Z,Y=Y,mu=mu,P=P,N=N,delta=delta,R0=R0,
        E=sy.Function('actual_power_inlet_E')(Z)*sy.exp(-(sy.Rational(1,2)+mu)*Y),
        power_node=power.node,amplitude_nodes=amplitude)


def actual_flat_power_source(owner):
    context=power_context(owner);env=context['env'];Z=context['Z'];Y=context['Y'];mu=context['mu'];E=context['E']
    g=owner.g;d=owner.leading['chart'];alg=RcAlgebra(owner,env);roots=[]
    for key,value in (('a',2+2*mu),('b',sy.Integer(0)),('Delta',2*mu)):
        for n,ref in enumerate(rows(owner.jet(d['roots'][key]))):
            zero(alg.at(ref)-(value if n==0 else 0),('actual O3 root function',key,n));roots.append(key+'_Z'+str(n))
    # The exact roots above and the canonical source admission establish
    # this actual predicate before any active square-root/inverse body.
    assert owner.leading['mu_proof']['positive_mu_and_power_2mu_ge_eta_and_Rp_reservation']
    alg=RcAlgebra(owner,env,True);flat=[];density=[];prefix=[]
    for key in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi'):
        for n,ref in enumerate(rows(owner.jet(d['loop'][key]))):zero(alg.at(ref),('actual O3 flat loop',key,n));flat.append(key+'_Z'+str(n))
    for key,ids in d['finite_N']['five_signed_density_increments'].items():
        for n,ref in enumerate(rows(owner.jet(ids))):zero(alg.at(ref),('actual O3 local density',key,n));density.append(key+'_Z'+str(n))
    for j,row in enumerate(owner.profiles):
        for key,wanted in (('E',sy.diff(E,Y,j)),('V',sy.Integer(0))):
            for n,ref in enumerate(rows(row[key])):
                zero(alg.at(ref)-sy.diff(wanted,Z,n),('actual power/first-native/higher profile',key,j,n));prefix.append(key+'_y'+str(j)+'_Z'+str(n))
    return dict(passed=True,actual_power_node=context['power_node'],actual_query_independent_amplitude_nodes=context['amplitude_nodes'],
        actual_constant_root_C5_rows=roots,actual_lazy_flat_loop_C5_rows=flat,actual_zero_local_density_C5_rows=density,
        actual_profile_y0_to_y4_C5_rows=prefix,canonical_source_guard_consumed_before_active_body=True,
        actual_first_radial_prefix_not_replaced_by_compact_rows=True,nonzero_incoming_memories_not_zeroed=True)


def ode_reduce(expr,H,density,Y,Z,rates):
    for _ in range(12):
        substitutions={}
        for term in expr.atoms(sy.Derivative):
            if term.expr not in H.values():continue
            orders=dict(term.variable_count);j=orders.get(Y,0);n=orders.get(Z,0)
            if j:
                key=next(k for k,value in H.items() if value==term.expr)
                substitutions[term]=sy.diff(density[key]-sy.Rational(rates[key])*H[key],Y,j-1,Z,n)
        if not substitutions:return expr
        expr=expr.xreplace(substitutions).doit()
    raise AssertionError('Unreduced actual history ODE derivatives')


def actual_history_and_recovery(owner):
    context=power_context(owner);env=dict(context['env']);Z,Y,E=context['Z'],context['Y'],context['E']
    g=owner.g;d=owner.leading;rates=current.compact_source.RATES
    for label,packet in (('background',d['background']['incoming']),('correction',d['chart']['incoming'])):
        for key,ids in packet.items():bind_input(env,owner,ids,'actual_'+label+'_incoming_'+key,Z)
    alg=RcAlgebra(owner,env,True);first=[];correction=[]
    densities=dict(m=sy.Integer(0),h=E,k=sy.Integer(0),e=-E**2/2,p=E**2/2)
    for key,rate in rates.items():
        base=alg.at(owner.histories[0][key].value)
        zero(sy.diff(base,Y)-(densities[key]-sy.Rational(rate)*base),('actual full O3 Duhamel ODE',key))
        first.append(key+'_whole_Duhamel_ODE')
        incoming=RcAlgebra(owner,env,True).at(d['chart']['incoming'][key][0])
        zero(alg.at(d['chart']['history'][key][0])-incoming*sy.exp(-sy.Rational(rate)*Y),('actual nonzero correction memory',key))
        correction.append(key+'_whole_nonzero_memory')
    zero(alg.at(owner.histories[1]['m'].value)+alg.at(owner.histories[0]['m'].value),'actual native mass_y prefix from Duhamel')
    H={key:sy.Function('actual_full_history_'+key)(Y,Z) for key in rates};env=dict(context['env'])
    for key,jet in owner.histories[0].items():bind(env,jet,H[key],Z)
    # The actual native first mass source has just been proved; this
    # input identification makes subsequent ordinary commutation explicit.
    bind(env,owner.histories[1]['m'],-H['m'],Z)
    alg=RcAlgebra(owner,env,True);whole=[];derived=[]
    for j,row in enumerate(owner.histories[1:],1):
        for key,jet in row.items():
            wanted=ode_reduce(sy.diff(H[key],Y,j),H,densities,Y,Z,rates)
            zero(alg.at(jet.value)-wanted,('actual full history higher source',key,j));whole.append(key+'_y'+str(j)+'_Z0')
            derived.extend(key+'_y'+str(j)+'_Z'+str(n) for n in range(1,6))
    # The exact mass ODE with zero local axial velocity identifies its
    # actual amplitude times exp(-y); no amplitude or datum is zeroed.
    env=dict(context['env']);M=sy.Function('actual_nonzero_O3_mass_amplitude')(Z)*sy.exp(-Y)
    bind(env,owner.histories[0]['m'],M/context['P'],Z);bind(env,owner.histories[1]['m'],-M/context['P'],Z)
    alg=RcAlgebra(owner,env,True);recovery=[]
    Q=(-(1-context['delta'])*Z*M-(1-Z**2)*sy.diff(M,Z))/(1-context['delta']*Z**2)
    for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
        for key,body in (('U',sy.Integer(0)),('m',sy.diff(M,Y,j)),('Q',sy.diff(Q,Y,j))):
            jet=owner.physical[key+suffix]
            for n,ref in enumerate(rows(jet)):
                zero(alg.at(ref)-sy.diff(body,Z,n),('actual physical O3 recovery',key,j,n));recovery.append(key+'_y'+str(j)+'_Z'+str(n))
    return dict(passed=True,actual_full_background_plus_correction_whole_ODE_rows=first,
        actual_nonzero_correction_whole_memory_rows=correction,actual_higher_history_whole_rows=whole,
        higher_history_Z_rows_by_pinned_ordinary_commutation=derived,actual_physical_U_m_Q_y0_to_y4_rows=recovery,
        original_first_mass_native_source_compared_before_input_identification=True,
        original_pressure_rate_zero_and_all_nonzero_memories_retained=True,leading_one_minus_delta_recovery=True)


def inlet_context(owner):
    context=power_context(owner);env=dict(context['env']);g=owner.g;d=owner.leading;Z,Y=context['Z'],context['Y']
    for label,packet in (('background',d['background']['fields']),('correction',d['chart']['history'])):
        for key,ids in packet.items():
            if key not in current.compact_source.RATES:continue
            bind_input(env,owner,ids,'actual_'+label+'_history_'+key,Z,Y)
    alg=RcAlgebra(owner,env,True)
    B=alg.at(d['background']['fields']['m'][0]);C=alg.at(d['chart']['history']['m'][0])
    bind(env,owner.histories[1]['m'],-B-C,Z)
    c=owner.c.c;env[owner.c.x.node]=sy.Symbol(g.nodes[owner.c.x.node]['name'],real=True)
    env[c['Rc']]=context['R0']*sy.exp(2);env[c['W']['normal']]=sy.Symbol('actual_compact_beta_normalization',positive=True)
    for j,ids in enumerate(c['h']):bind(env,owner.jet(ids),sy.Function('actual_selected_compact_control_'+str(j))(Z),Z)
    for key,ids in d['chart']['cumulative_integrals'].items():bind(env,owner.jet(ids),sy.Function('actual_O3_integral_'+key)(Y,Z),Z)
    for key,ids in d['initial'].items():bind(env,owner.jet(ids),sy.Function('actual_initial_'+key)(Z),Z)
    return dict(context,env=env)


def actual_Rc_leading_inlet(owner):
    context=inlet_context(owner);env=context['env'];Z=context['Z'];g=owner.g;c=owner.c.c;d=owner.leading
    alg=RcAlgebra(owner,env,True);profiles=[];histories=[];physical=[];pressure=[];empty=[];derived=[]
    for group in ('contributions','pressure_local'):
        for key,ids in c[group].items():
            row=g.nodes[ids[0]];assert row['operation']=='definite_integral' and row['lower']==1 and row['upper']==owner.c.x.node
            for n,ref in enumerate(rows(Algebra(g,len(ids)-1).at(owner.jet(ids),owner.c.x,g.one))):
                zero(alg.at(ref),('actual compact empty inlet primitive',group,key,n));empty.append(group+'_'+key+'_Z'+str(n))
    for j in range(5):
        for key in ('E','V'):
            left,right=owner.leading_left['profiles'][j][key],owner.leading_right['profiles'][j][key]
            assert len(rows(left))==len(rows(right))==6
            zero(alg.at(left.value)-alg.at(right.value),('actual Rc whole profile source join',key,j));profiles.append(key+'_y'+str(j)+'_Z0')
            derived.extend('profile_'+key+'_y'+str(j)+'_Z'+str(n) for n in range(1,6))
        for key in current.compact_source.RATES:
            left,right=owner.leading_left['histories'][j][key],owner.leading_right['histories'][j][key]
            assert len(rows(left))==len(rows(right))==6
            zero(alg.at(left.value)-alg.at(right.value),('actual Rc whole full history source join',key,j));histories.append(key+'_y'+str(j)+'_Z0')
            derived.extend('history_'+key+'_y'+str(j)+'_Z'+str(n) for n in range(1,6))
    for key,left in owner.leading_left['physical'].items():
        right=owner.leading_right['physical'][key];assert len(rows(left))==len(rows(right))
        zero(alg.at(left.value)-alg.at(right.value),('actual Rc whole physical source join',key));physical.append(key+'_Z0')
        derived.extend('physical_'+key+'_Z'+str(n) for n in range(1,len(rows(left))))
    for key,ids in d['endpoint']['pressure'].items():
        left,right=owner.jet(ids),owner.jet(c['inlet']['pressure'][key]);assert len(rows(left))==len(rows(right))==4
        zero(alg.at(left.value)-alg.at(right.value),('actual Rc whole nonzero pressure primitive identity',key));pressure.append(key+'_Z0')
        derived.extend('primitive_'+key+'_Z'+str(n) for n in range(1,4))
    for packet in ('density','higher_density'):
        for key,left in owner.leading_left[packet].items():
            right=owner.leading_right[packet][key];assert len(rows(left))==len(rows(right))==4
            zero(alg.at(left.value)-alg.at(right.value),('actual Rc whole pressure radial source identity',packet,key));pressure.append(packet+'_'+key+'_Z0')
            derived.extend(packet+'_'+key+'_Z'+str(n) for n in range(1,4))
    return dict(passed=True,actual_compact_empty_inlet_integral_rows=empty,actual_leading_Rc_profile_rows=profiles,
        actual_leading_Rc_complete_history_rows=histories,actual_leading_Rc_physical_rows=physical,
        actual_Rc_pressure_primitive_and_radial_source_rows=pressure,
        higher_Z_rows_from_proved_whole_function_identity_and_pinned_ordinary_commutation=derived,
        actual_pressure_outputs_not_bound=True,actual_background_and_correction_inlets_not_zeroed=True,
        actual_shared_initial_S0_S1_and_Q_Qy_finite_boundary_data=True)


def actual_pressure_FTC_and_higher_source(owner):
    context=power_context(owner);env=dict(context['env']);Z,Y=context['Z'],context['Y'];g=owner.g
    M=sy.Function('actual_O3_physical_mass_amplitude')(Z)*sy.exp(-Y)
    bind(env,owner.histories[0]['m'],M/context['P'],Z);bind(env,owner.histories[1]['m'],-M/context['P'],Z)
    alg=RcAlgebra(owner,env,True);oldS=owner.jet(owner.leading['chart']['S']);prefix=[];higher=[];derived=[]
    zero(alg.at(oldS.value)-alg.at(owner.density['S'].value),'actual whole old/new O3 reduced S source');prefix.append('S_Z0')
    derived.extend('S_Z'+str(n) for n in range(1,4))
    Py=alg.at(owner.density['P1_y'].value)
    for key,j in (('P1_yy',1),('P1_yyy',2)):
        zero(alg.at(owner.higher_density[key].value)-sy.diff(Py,Y,j),('actual whole O3 pressure high derivative source',key));higher.append(key+'_Z0')
        derived.extend(key+'_Z'+str(n) for n in range(1,4))
    # Actual S0/S1 integrals and finite boundary terms, with the common
    # original S input identified above. Pressure outputs stay unbound.
    S=sy.Function('actual_O3_reduced_S')(Y,Z);Q=sy.Function('actual_O3_Q_amplitude')(Z)*sy.exp(-Y)
    env=dict(context['env'])
    for jet in (oldS,owner.density['S']):bind(env,jet,S,Z)
    for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):bind(env,owner.physical['Q'+suffix],sy.diff(Q,Y,j),Z)
    d=owner.leading['chart']
    for group in ('incoming_integrals',):
        for key,ids in d[group].items():
            assert all(not outer_check.free_query_dependency(owner.own,node,owner.coordinate.node) for node in ids)
            bind_input(env,owner,ids,'actual_retained_O3_'+key,Z)
    for key,ids in owner.leading['initial'].items():bind_input(env,owner,ids,'actual_initial_boundary_'+key,Z)
    alg=RcAlgebra(owner,env,True);FTC=[];R=context['R0']*sy.exp(Y)
    for key,weight in (('m13',sy.Integer(1)),('m15',-R)):
        zero(sy.diff(alg.at(d['pressure'][key][0]),Y)-weight*alg.at(owner.density['P1_y'].value),('actual whole O3 signed pressure FTC',key))
        FTC.append(key+'_y_Z0');derived.extend(key+'_y_Z'+str(n) for n in range(1,4))
    for key in ('m13','m15'):
        chart=owner.intermediate;complete=owner.own.inputs[current.outer.intermediate.previous.previous.NAME]
        expected=Algebra(g,3).add(owner.jet(complete['segments']['inner'][key]),owner.jet(d['pressure'][key]),owner.jet(complete['segments']['lower'][key]))
        assert encode(expected)==chart['cumulative'][key]
        for segment in ('inner','lower'):
            for node in complete['segments'][segment][key]:assert not outer_check.free_query_dependency(owner.own,node,owner.coordinate.node)
    return dict(passed=True,actual_old_and_new_reduced_S_C3_rows=prefix,
        actual_pressure_yy_yyy_C3_derivative_rows=higher,actual_signed_S0_S1_pressure_FTC_rows=FTC,
        higher_Z_rows_by_proved_whole_source_identities_and_pinned_integral_commutation=derived,
        actual_own_primitive_is_same_inner_lower_plus_actual_O3_pressure=True,
        pressure_outputs_not_bound=True,actual_original_R_and_Rsquared_measures_and_signed_finite_boundaries_kept=True)


def actual_own_Rc_mixed_join(owner):
    g=owner.g;context=inlet_context(owner);env=dict(context['env']);Z=context['Z'];c=owner.c
    compact_rec=combined_check.admitted(owner,current.previous.previous)
    assert compact_rec['checks']['actual_total_pressure_source']['actual_whole_compact_reserved_radius_source_identity']
    own_context=compact_check.own_context(c)['bindings'];theta=compact_check.theta0_context(c)['bindings']
    env.update(own_context);env.update(theta)
    env.update({g.symbol('Z').node:Z,owner.mu.node:context['mu'],owner.Pstar.node:context['P'],owner.delta.node:context['delta'],
        owner.c.c['Rc']:context['R0']*sy.exp(2),owner.c.own_R.node:context['R0']*sy.exp(2)*sy.exp(sy.Symbol(g.nodes[c.phase.node]['name'],real=True))})
    complete=owner.own.inputs[current.outer.intermediate.previous.previous.NAME]
    for key in ('m11','m12','m14'):
        assert owner.intermediate['cumulative'][key]==c.repaired['incoming'][key]
        for query in (owner.coordinate.node,c.phase.node):
            assert all(not outer_check.free_query_dependency(owner.own,node,query) for node in c.repaired['incoming'][key])
        bind_input(env,owner,c.repaired['incoming'][key],'actual_shared_own_'+key,Z)
    for key in ('m13','m15'):
        for segment in ('inner','lower'):
            ids=complete['segments'][segment][key]
            for query in (owner.coordinate.node,c.phase.node):
                assert all(not outer_check.free_query_dependency(owner.own,node,query) for node in ids)
            bind_input(env,owner,ids,'actual_retained_'+segment+'_'+key,Z)
    alg=RcAlgebra(owner,env,True);evidence=[]
    # No own state, pressure, cumulative output or seam trace is bound.
    # Both complete recovery recipes consume the same identified inputs.
    for spec in owner.contract:
        group,key=spec['group'],spec['field']
        for j in range(spec['physical_y_order']+1):
            left,right=owner.join['left'][group][key][j],owner.join['right'][group][key][j]
            assert len(rows(left))==len(rows(right))==spec['ordinary_Z_order']+1
            zero(alg.at(left.value)-alg.at(right.value),('actual own Rc whole mixed source join',group,key,j))
            for n,(l,r) in enumerate(zip(rows(left),rows(right))):
                evidence.append(dict(group=group,field=key,physical_y_order=j,ordinary_Z_order=n,
                    actual_left_node=l.node,actual_right_node=r.node,
                    proof='direct whole function identity' if n==0 else 'ordinary Z derivative of proved whole function identity; pinned source algebra and fixed endpoints'))
    assert len(evidence)==246
    # Radius identification is an input theorem, not a new free common
    # radius assumption. The compact Rc node is the actual O3 endpoint.
    assert c.c['Rc']==owner.leading['endpoint']['radius']
    return dict(passed=True,actual_Rc_contracted_mixed_join_rows=evidence,
        actual_compact_reserved_whole_radius_identity_receipt=current.previous.previous.RECEIPT,
        original_Rc_selected_parameter_and_own_recovery_sign_preserved=True,
        original_nonzero_cumulative_and_pressure_data_retained=True,whole_outputs_not_bound=True)


CHECKS={'source_contract','actual_flat_power_source','actual_history_and_recovery','actual_Rc_leading_inlet',
    'actual_pressure_FTC_and_higher_source','actual_own_Rc_mixed_join'}


def save(owner,checks,execution_seconds):
    owner.assert_graph();assert set(checks)==CHECKS and all(proof['passed'] for proof in checks.values())
    evidence=checks['actual_own_Rc_mixed_join']['actual_Rc_contracted_mixed_join_rows']
    expected={(s['group'],s['field'],j,n) for s in owner.contract for j in range(s['physical_y_order']+1) for n in range(s['ordinary_Z_order']+1)}
    assert len(evidence)==len(expected)==246
    assert {(r['group'],r['field'],r['physical_y_order'],r['ordinary_Z_order']) for r in evidence}==expected
    raw=current.parent.linear.read(current.NAME)
    assert raw['Rc_mixed_source_join_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_Rc_contracted_mixed_join_rows=246,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=execution_seconds)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_OWN_RC_246_MIXED_SOURCE_ROWS_ACCEPTED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in (source_contract,actual_flat_power_source,actual_history_and_recovery,actual_Rc_leading_inlet,
        actual_pressure_FTC_and_higher_source,actual_own_Rc_mixed_join):
        checks[fn.__name__]=fn(owner);print('RC_'+fn.__name__.upper()+'_PASSED',flush=True)
    return save(owner,checks,time.monotonic()-began)
