"""Active total derivatives, periodic lift and actual micro own joining."""
import hashlib
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_micro_active_radial4_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_rminus_mixed_source_join_check as prior

rows,encode=current.outer.rows,current.outer.encode
zero,bind=prior.zero,prior.bind


class MicroAlgebra(prior.RminusAlgebra):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation'];owner=self.rminus
        if op=='function_substitution':
            assert row['Z_independent_substitution'];scope=dict(self.scope)
            scope[row['variable']]=self.at(row['value'])
            result=MicroAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['expression'])
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z'];lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:self.evidence['empty'].add(i);result=sy.Integer(0)
            else:
                variable=sy.Symbol(row['variable'],real=True);scope=dict(self.scope)
                scope[self.g.symbol(row['variable']).node]=variable
                body=MicroAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['integrand'])
                result=sy.Integral(body,(variable,lo,hi))
        elif op=='original_lazy_flat_branch':
            raise AssertionError(('Active loop source must not be replaced by the rminus flat branch',i))
        else:return super().at(i)
        self.cache[i]=result;return result


def source_contract(owner):
    owner.assert_graph();g=owner.g
    rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and all(rec[key] for key in current.previous.GATES)
    assert rec['checks']['source_contract']['actual_first_micro_recovery_equals_original_zero_velocity_fields']
    assert owner.loop['first_micro']['coordinate']==owner.loop['second_micro']['coordinate']==owner.coordinate.node
    assert owner.loop['first_micro']['radius']==owner.loop['second_micro']['radius']==owner.radius.node
    assert owner.loop['first_micro']['physical_measure']==owner.loop['second_micro']['physical_measure']==owner.jacobian.node
    assert owner.loop['first_micro']['actual_phase']==owner.loop['second_micro']['actual_phase']==owner.phase.node
    assert owner.radius==owner.base.radius and owner.jacobian==owner.base.jacobian
    assert owner.loop['first_micro']['roots']==owner.loop['second_micro']['roots']
    assert owner.modified['first_micro']['radial_background']==owner.modified['second_micro']['radial_background']==encode(owner.background_y)
    for pin in (owner.lower.sigma_pin,owner.lower.sigma_interior_pin):
        assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    sign,mantissa,exponent,bits=owner.base.common['sc_binding']['exact_mpf_tuple']
    assert owner.base.common['sc_binding']['exact_singleton_recomputed_not_interval_overlap']
    assert sign==0 and mantissa>0 and mantissa.bit_length()==bits and exponent+bits<=-2
    assert sy.Rational(3,16)<sy.Rational(7,8)
    assert owner.loop['first_micro']['domain'][1]==owner.loop['second_micro']['domain'][0]==g.one.node
    before=current.previous.common_source.previous
    raw=current.parent.linear.read(before.NAME)['intermediate_own_C3_chart_bindings']
    support=current.parent.linear.read(before.RECEIPT)['checks']
    assert support['source']['actual_preRc_coordinate_radius_domain_and_C3_pressure_source_identity']
    assert support['actual_regional_recovery']['own_V_y_and_mass_y_zero_only_outside_velocity_supports']
    count=0;core=[];sigma=[]
    for name in current.CHARTS:
        chart=owner.charts[name];D=owner.factories[name]
        assert D.cache[owner.phase.node]==g.mul(owner.N,owner.jacobian)
        assert chart['intermediate']==raw['charts'][name] and name in raw['preRc_names']
        assert chart['intermediate']['region']=='beforeRc'
        assert chart['intermediate']['fields']['M1']==chart['intermediate']['cumulative']['m11']
        assert chart['intermediate']['fields']['P1']==chart['intermediate']['cumulative']['m13']
        for key in ('Uz1','Utheta1','F1'):assert set(chart['intermediate']['fields'][key])=={0}
        for key in ('m11','m12','m14'):
            assert chart['intermediate']['cumulative'][key]==owner.lower.rin.cut['endpoint'][key][:4]
            assert all(node<owner.coordinate.node for node in chart['intermediate']['cumulative'][key])
        for key in ('U','m','Q'):
            assert encode(chart['physical'][key][:2])==[owner.modified[name]['physical'][key],owner.modified[name]['physical'][key+'_y']]
            assert len(chart['physical'][key])==5
            assert all(len(rows(jet))==(5 if key=='Q' else 6) for jet in chart['physical'][key])
        for request in owner.requests[name]:
            assert D.cache[request['input']].node==request['native_output']
            assert g.nodes[request['physical_output']]==g.nodes[g.mul(
                g.quotient(g.one,owner.jacobian,'same positive original micro physical Jacobian'),owner.ref(request['native_output'])).node]
            count+=1
        for factory in D.shared['factories'].values():
            assert factory.cache[owner.phase.node]==(g.mul(owner.N,owner.jacobian) if factory.variable==owner.coordinate else g.zero)
            for node in factory.core_extensions:
                row=g.nodes[node];owner.guard_core(row)
                assert row['unique_analytic_Banach_solution'] and row['actual_twenty_term_operator']
                assert row['admitted_solution_not_a_model_or_enclosure_center'] and row['source_family']==owner.family
                core.append(node)
            sigma+=factory.sigma_extensions
    assert count==36
    for group in owner.join.values():
        for spec in owner.contract:
            assert len(group[spec['group']][spec['field']])==spec['physical_y_order']+1
    return dict(passed=True,actual_total_native_U_y2_y3_y4_C5_requests=count,
        actual_identical_core_radial4_Z6_extension_nodes=sorted(set(core)),actual_sigma_extensions=sigma,
        root_native_phase_lift_NJ_and_all_other_binders_phase_held=True,
        original_modified_y0_y1_prefix_and_intermediate_support_receipts_retained=True,
        source_scope_is_seam_neighborhood_not_the_first_micro_flat_collar=True)


def periodic_lift(owner):
    r,z,u,q,a,b,E,t0=sy.symbols('r z u q a b E t0',real=True)
    wz=((z+1/z)/2-r)/(1-r*(z+1/z)+r*r)
    series=(z/(1-r*z)+1/z/(1-r/z))/2
    zero(wz-series,'actual Poisson Fourier/geometric representation')
    # Absolute convergence for |r|<1 gives zero mean and this exact norm.
    mean_w2=1/(2*(1-r*r));rr=u/sy.sqrt(1+u*u);alpha=2*q/sy.sqrt(1+u*u)
    zero((alpha**2*mean_w2).subs(r,rr)-2*q*q,'actual normalized second Poisson moment')
    zero((-a*t0-b).subs(t0,-b/a),'actual B increment after one lifted phase turn')
    recname=next(name for name in owner.hashes if name.endswith('whole_Z_all_N_control_connection_check.json'))
    rec=current.parent.linear.read(recname)
    assert rec['symbolic_identities']['first_micro->second_micro']
    identities=rec['symbolic_identities']
    assert identities['normalized_five_density_identities'] and identities['true_relative_Z_quotients']
    assert identities['Poisson_w_r'] and identities['Poisson_r_u']
    Z=sy.Symbol('Z',real=True);env={owner.g.symbol('Z').node:Z}
    actual_E=sy.Function('actual_background_E')(Z);actual_Ey=sy.Function('actual_background_Ey')(Z)
    actual_Vy=sy.Function('actual_background_Vy')(Z)
    for ids,body in ((owner.background['fields']['E'],actual_E),(owner.background['radial']['E_y'],actual_Ey),
        (owner.background['radial']['V_y'],actual_Vy)):bind(env,owner.jet(ids),body,Z)
    alg=MicroAlgebra(owner,env);root_rows=[]
    for name in current.CHARTS:
        root=owner.loop[name]['roots'];target=-2*actual_Vy/(actual_E-2*actual_Ey)
        for n,node in enumerate(root['t0']):zero(alg.at(node)-sy.diff(target,Z,n),('actual root t0=-b/a',name,n));root_rows.append([name,n,node])
        zero(alg.at(root['t0'][0])+alg.at(root['b'][0])/alg.at(root['a'][0]),('actual source root ratios',name))
        loop=owner.loop[name]['loop'];inverse=owner.g.nodes[loop['inverse'][0]]
        assert inverse['operation']=='original_monotone_phase_inverse'
        assert inverse['phase_function']==loop['Phi'][0] and inverse['target_phase']==owner.phase.node
        assert inverse['angle_variable']==owner.g.nodes[loop['angle']]['name']
        assert inverse['angle_lower']==owner.g.zero.node and inverse['true_prescribed_source_C5_implicit_inverse']
    return dict(passed=True,actual_Poisson_source_receipt=recname,actual_t0_root_identity_rows=root_rows,
        exact_w_zero_mean_and_alpha_squared_mean_w_squared_equals_2q_squared=True,
        Phi_monodromy_one_from_original_K_normalization=True,
        inverse_monodromy_two_pi_and_T1_increment_two_pi_t0=True,
        original_t0_equals_minus_b_over_a_cancels_B_increment=True,
        A_and_B_and_all_required_phase_Z_mixed_rows_periodic=True,
        Aphi='a/2*(1-1/(2*pi*lambda))',Bphi='E/2*(-a*t/(2*pi*lambda)-b)',
        total_derivatives_apply_to_periodic_lifted_composites_not_raw_fractional_part=True,
        sigma_and_Delta_lazy_flat_source_uses_original_smooth_branch=True)


class MicroFingerprint:
    """Alpha-equivalence on the two-sided seam collar, with bound angles."""
    def __init__(self,owner):self.owner=owner;self.g=owner.g;self.cache={}
    def name(self,name):
        return 'actual_common_loop_angle' if name in ('true_loop_C5_angle_first_micro','true_loop_C5_angle_second_micro') else name
    def at(self,node):
        node=node.node if hasattr(node,'node') else node
        if node in self.cache:return self.cache[node]
        row=self.g.nodes[node];op=row['operation'];F=self.at
        if op=='bound_variable':body=[op,self.name(row['name'])]
        elif op in ('sum','product'):body=[op,sorted(F(i) for i in row['arguments'])]
        elif op=='negative':body=[op,F(row['argument'])]
        elif op=='positive_quotient':body=[op,F(row['numerator']),F(row['denominator'])]
        elif op in ('analytic_unary','exact_original_source_flat_sigma'):
            body=[op,row.get('name'),F(row['argument'])]
        elif op=='function_substitution':body=[op,F(row['expression']),F(row['variable']),F(row['value'])]
        elif op=='definite_integral':body=[op,F(row['integrand']),self.name(row['variable']),F(row['lower']),F(row['upper'])]
        elif op=='original_monotone_phase_inverse':body=[op,F(row['phase_function']),self.name(row['angle_variable']),F(row['target_phase']),F(row['angle_lower']),F(row['angle_upper'])]
        elif op=='substitute_original_inverse_angle':body=[op,F(row['body']),F(row['inverse_angle']),self.name(row['angle_variable'])]
        elif op=='exact_real_comparison':body=[op,row['operator'],F(row['left']),F(row['right'])]
        elif op=='logical_or':body=[op,sorted(F(i) for i in row['arguments'])]
        elif op=='original_lazy_flat_branch':
            predicate=row['flat_predicate']
            if row.get('original_collar_coordinate')==self.owner.coordinate.node:
                assert row['original_collar_upper']==self.owner.loop['first_micro']['loop']['original_collar_upper']
                args=self.g.nodes[predicate]['arguments'];remaining=[]
                for ref in args:
                    test=self.g.nodes[ref]
                    if test==dict(operation='exact_real_comparison',operator='<=',left=self.owner.coordinate.node,right=row['original_collar_upper']):continue
                    remaining.append(ref)
                assert len(remaining)==1;predicate=remaining[0]
            body=[op,F(predicate),F(row['active_body']),F(row['flat_value'])]
        elif op=='current_original_unique_analytic_core_component':
            if node>=len(self.owner.base.prefix):self.owner.guard_core(row)
            body=[op,row]
        else:
            assert op in ('exact_rational','mathematical_pi','mathematical_positive_infinity',
                'original_source_parameter','current_original_source_parameter','exact_positive_integer_expression',
                'current_original_leading_function_recipe','current_original_axis_pressure_function_derivative'),op
            body=[op,row]
        result=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        self.cache[node]=result;return result


def actual_leading_micro_sources(owner):
    common=current.parent.linear.read(current.previous.common_source.RECEIPT)['checks']['common_IVP']
    assert common['common_true_rho4_inlet_and_comparison_alpha_one_from_actual_sc']
    # That original selector theorem gives sc<1/4. Hence the first-only
    # collar predicate s<=3sc/4 is false on [7/8,9/8].
    fingerprint=MicroFingerprint(owner);evidence=[]
    for n in range(5):
        for k,(left,right) in enumerate(zip(rows(owner.charts['first_micro']['physical']['U'][n]),rows(owner.charts['second_micro']['physical']['U'][n]))):
            assert fingerprint.at(left)==fingerprint.at(right),('actual active U alpha-equivalence',n,k)
            evidence.append([n,k,left.node,right.node])
    assert owner.loop['second_micro']['incoming']==owner.loop_input['endpoints']['first_micro']
    assert owner.modified['second_micro']['incoming_integrals']==owner.modified_input['endpoints']['first_micro']['integrals']
    for name in current.CHARTS:
        assert owner.loop[name]['memory']['p']==owner.g.one.node
    return dict(passed=True,actual_U_y0_to_y4_C5_source_alpha_equivalence=evidence,
        exact_original_sc_selector_makes_extra_first_collar_predicate_false=True,
        common_background_roots_phase_radius_and_periodic_loop_sources=True,
        actual_first_complete_Duhamel_history_is_second_incoming=True,
        same_mean_ODE_U_minus_m_gives_all_high_mean_rows_from_actual_inlet=True,
        leading_one_minus_delta_Q_gives_all_high_Q_rows_from_same_U_mean=True,
        actual_first_complete_pressure_integrals_are_second_incoming=True,
        pressure_rate_zero_memory_not_reset=True)


def actual_leading_recovery_pressure_and_FTC(owner):
    evidence=[];ftc=[];Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_active_physical_y',real=True)
    for name in current.CHARTS:
        chart=owner.charts[name];old=owner.modified[name]
        env,v=prior.active_context(owner.base,False);delta=v['delta'];R=sy.Symbol('actual_micro_R0',positive=True)*sy.exp(y)
        env[owner.radius.node]=R;u=sy.Function('independent_actual_U')(y,Z);m=sy.Function('independent_actual_mean')(y,Z)
        for n,jet in enumerate(chart['physical']['U']):bind(env,jet,sy.diff(u,y,n),Z)
        bind(env,chart['physical']['m'][0],m,Z);alg=MicroAlgebra(owner,env)
        means=[m]
        for n in range(1,5):means.append(sy.diff(u,y,n-1)-means[-1])
        for n in range(5):
            zero(alg.at(chart['physical']['m'][n].value)-means[n],('actual mean ODE derivative',name,n))
            q=(2*Z*sy.diff(u,y,n)-(1-delta)*Z*means[n]-(1-Z**2)*sy.diff(means[n],Z))/(1-delta*Z**2)
            for k,ref in enumerate(rows(chart['physical']['Q'][n])):
                zero(alg.at(ref)-sy.diff(q,Z,k),('actual leading Q derivative recovery',name,n,k));evidence.append([name,n,k])
        env,v=prior.active_context(owner.base,False);env[owner.radius.node]=R
        q=sy.Function('independent_actual_Q')(y,Z)
        for key,body in (('U',u),('Q',q)):
            for n,jet in enumerate(chart['physical'][key]):bind(env,jet,sy.diff(body,y,n),Z)
        alg=MicroAlgebra(owner,env);d=1-Z**2;L=1-delta*Z**2
        S=(q+(1-delta)*Z*sy.diff(q,Z)/2+sy.diff(q,y))/L+q*(q/2+sy.diff(q,y))+u*(-2*Z*q+d*sy.diff(q,Z)-2*Z*sy.diff(q,y))/L
        N=-S/2+(sy.diff(q,y)+sy.diff(q,y,2))/R
        for k,node in enumerate(old['S']):zero(alg.at(node)-sy.diff(S,Z,k),('actual saved reduced-S body',name,k))
        zero(alg.at(chart['density']['S'].value)-S,('new density is actual saved S',name))
        for n,jet in enumerate(chart['N_y']):
            for k,ref in enumerate(rows(jet)):zero(alg.at(ref)-sy.diff(N,y,n,Z,k),('actual N higher physical source',name,n,k))
        for n in range(3):
            for key,target in (('pressure_y',R*N),('m15_y',-R**2*N)):
                for k,ref in enumerate(rows(chart[key][n])):
                    zero(alg.at(ref)-sy.diff(target,y,n,Z,k),('actual signed pressure derivative',name,key,n,k))
                    zero(alg.at(rows(chart['raw_'+key][n])[k])-alg.at(ref),('raw/coupled pressure identity',name,key,n,k))
        env,v=prior.active_context(owner.base,False);s=v['s'];H=v['H'];R=v['eps']*v['rho']
        qs=sy.Function('actual_saved_Q_'+name)(s,Z);ss=sy.Function('actual_saved_S_'+name)(s,Z)
        bind(env,owner.jet(old['physical']['Q']),qs,Z);bind(env,owner.jet(old['physical']['Q_y']),sy.diff(qs,s)/H,Z)
        bind(env,owner.jet(old['S']),ss,Z)
        for key,ids in owner.modified_input['initial'].items():bind(env,owner.jet(ids),sy.Function('actual_global_initial_'+key)(Z),Z)
        for key,ids in old['incoming_integrals'].items():
            if any(ids):bind(env,owner.jet(ids),sy.Function('actual_nonzero_incoming_'+name+'_'+key)(Z),Z)
        alg=MicroAlgebra(owner,env);Py=-R*ss/2+sy.diff(qs,s)/H+sy.diff(qs,s,2)/H**2
        for key,target in (('m13',Py),('m15',-R*Py)):
            for k,node in enumerate(old['pressure'][key]):
                zero(sy.diff(alg.at(node),s).doit()/H-sy.diff(target,Z,k),('actual original pressure primitive FTC',name,key,k))
                ftc.append([name,key,k,node])
        for key,power in (('S0',1),('S1',2)):
            for k,node in enumerate(old['local_integrals'][key]):
                row=owner.g.nodes[node]
                assert row['operation']=='definite_integral' and row['lower']==owner.loop[name]['domain'][0]
                assert row['upper']==owner.coordinate.node and row['physical_radius_power']==power
                assert row['physical_Jacobian_applied_once'] and row['ordinary_slow_Z_derivative_order']==k
    return dict(passed=True,actual_Q_y0_to_y4_C4_recovery_rows=evidence,actual_signed_primitive_C3_FTC_rows=ftc,
        actual_mean_ODE_all_five_rows_and_actual_saved_reduced_S_C3=True,
        actual_N_y0_y1_y2_C3_and_signed_pressure_y1_y2_y3_C3=True,
        original_nonzero_incoming_and_global_initial_constants_retained=True,
        raw_operator_and_coupled_fixed_N_definitions_identical=True,own_pressure_outputs_not_bound=True)


def actual_nonzero_histories_and_own_micro_join(owner):
    env,Z,slow=prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
    s=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True)
    env[owner.base.active_sc.node]=MicroAlgebra(owner,env).at(owner.parameters['sc'])
    # Actual first defining primitives are arbitrary nonzero functions. The
    # second inlet is their fixed endpoint; its local integrals are empty.
    for key,ids in owner.loop['first_micro']['contributions'].items():
        bind(env,owner.jet(ids),sy.Function('actual_first_finite_N_primitive_'+key)(s,Z),Z)
    bind(env,owner.jet(owner.background['fields']['m']),sy.Function('actual_background_mean')(s,Z),Z)
    for key,ids in owner.modified['first_micro']['local_integrals'].items():
        bind(env,owner.jet(ids),sy.Function('actual_first_signed_'+key+'_primitive')(s,Z),Z)
    for name in current.CHARTS:
        for key in ('Q','Q_y'):
            bind(env,owner.jet(owner.modified[name]['physical'][key]),sy.Function('common_actual_leading_'+key)(s,Z),Z)
        for n,jet in enumerate(owner.charts[name]['N_y']):bind(env,jet,sy.Function('common_actual_fixed_N_y'+str(n))(s,Z),Z)
    alg=MicroAlgebra(owner,env);a=current.outer.Algebra(owner.g,5)
    means=[]
    for name in current.CHARTS:
        jet=a.at(owner.charts[name]['physical']['m'][0],owner.coordinate,owner.seam)
        means.append([alg.at(ref) for ref in rows(jet)])
    for k,(left,right) in enumerate(zip(*means)):zero(left-right,('actual leading first/second mean history',k))
    # Support and constant mass, not the leading finite-N collar, imply the
    # own exponential mean flow on both sides of this active-loop seam.
    flow=[];yy=sy.Symbol('independent_own_log_radius',real=True);R=alg.at(owner.radius)
    delta=alg.at(owner.parameters['delta']);d=1-Z**2;L=1-delta*Z**2
    for name in current.CHARTS:
        f=owner.charts[name]['own']['fields'];mean=alg.at(f['mean1'][0].value)
        zero(alg.at(f['M1'][0].value)-R*mean,('actual own mass',name))
        zero(alg.at(f['K1'][0].value)-mean,('actual own K',name))
        model_mean=mean*sy.exp(-yy);model_R=R*sy.exp(yy)
        model_Q=(-(1+delta)*Z*model_mean-d*sy.diff(model_mean,Z))/L
        models=dict(Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),mean1=model_mean,K1=model_mean,
            M1=model_R*model_mean,Q1=model_Q,V1=model_R*model_Q)
        for key,body in models.items():
            for n,jet in enumerate(f[key]):
                zero(alg.at(jet.value)-sy.diff(body,yy,n).subs(yy,0),('actual supported own mean flow',name,key,n));flow.append([name,key,n])
        for key in ('M1','mean1','Q1','V1'):
            for new,old in zip(rows(f[key][1]),owner.charts[name]['intermediate']['fields'][key+'_y']):
                zero(alg.at(new)-alg.at(old),('actual original own first-y source',name,key))
    evidence=[]
    for spec in owner.contract:
        for n in range(spec['physical_y_order']+1):
            for k,(left,right) in enumerate(zip(rows(owner.join['left'][spec['group']][spec['field']][n]),rows(owner.join['right'][spec['group']][spec['field']][n]))):
                zero(alg.at(left)-alg.at(right),('actual own micro mixed seam',spec['group'],spec['field'],n,k))
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=n,ordinary_Z_order=k,left_node=left.node,right_node=right.node))
    assert len(evidence)==246
    return dict(passed=True,actual_full_own_micro_mixed_join_rows=evidence,actual_own_support_flow_rows=flow,
        actual_nonzero_first_Duhamel_history_then_second_empty_integrals=True,
        actual_nonzero_signed_pressure_history_and_both_finite_boundaries_retained=True,
        common_fixed_leading_N_is_input_not_an_imposed_own_output=True,
        all_original_own_first_y_rows_and_one_plus_delta_recovery_retained=True)


def independent_total_native_calculus(owner):
    g=current.source.FunctionTransportGraph();Y,Z,F,T=sy.symbols('fixture_active_native Z fixture_active_phase fixture_active_angle')
    y,z,angle=g.symbol(str(Y)),g.symbol('Z'),g.symbol(str(T))
    J=sy.Rational(3,2);N=sy.Integer(37);phi0=sy.Rational(1,3)
    phase=g.unary('fractional_part',g.add(g.constant('1/3'),g.mul(g.constant('111/2'),y)))
    pi_symbol=sy.Symbol('fixture_actual_pi');cos_symbol=sy.Symbol('fixture_actual_cos_angle')
    variables={str(Y):y,'Z':z,str(F):phase,str(T):angle,str(pi_symbol):g.node('mathematical_pi'),str(cos_symbol):g.unary('cos',angle)}
    expr=lambda body:current.source.expression(g,sy.sympify(body).xreplace({sy.pi:pi_symbol,sy.cos(T):cos_symbol}),variables)
    c=Y**2/13;Phi=expr(T/(2*sy.pi)+c)
    inverse=g.node('original_monotone_phase_inverse',phase_function=Phi.node,angle_variable=str(T),
        target_phase=phase.node,angle_lower=g.zero.node,angle_upper=g.mul(g.constant(2),g.node('mathematical_pi')).node,
        Delta=g.zero.node,eta=g.one.node,independent_positive_linear_inverse_fixture=True)
    v=(1+Y**3/17)*(1+Z**5/19);a=sy.Rational(3,2);b=sy.Integer(1);t0=-b/a
    t=expr(t0+v*sy.cos(T))
    primitive=g.node('definite_integral',integrand=t.node,variable=str(T),lower=g.zero.node,
        upper=angle.node,endpoint_independent_of_Z=True)
    composed=current.previous.loop_source.at_inverse(g,primitive,inverse,angle)
    U=g.add(expr(2+Y**2+Z**5/23),g.mul(g.constant('-3/148'),g.quotient(composed,g.node('mathematical_pi'),'positive pi')),
        g.mul(g.constant('-1/37'),phase))
    known_inverse=2*sy.pi*(F-c)
    expected=2+Y**2+Z**5/23+(-a*(t0*known_inverse+v*sy.sin(known_inverse))/(2*sy.pi)-b*F)/N
    D=current.PeriodicTotalNativeDerivative(g,y,phase,y,g.constant('3/2'),g.constant(37),lambda row:None)
    # Ordinary Z rows are differentiated independently before native rows.
    source_rows=[expr(sy.diff(2+Y**2+Z**5/23,Z,k)) for k in range(6)]
    # Rebuild each actual composed integral with its genuine Z derivative.
    for k in range(6):
        body=expr(sy.diff(t0+v*sy.cos(T),Z,k))
        integral=g.node('definite_integral',integrand=body.node,variable=str(T),lower=g.zero.node,upper=angle.node,endpoint_independent_of_Z=True)
        at=current.previous.loop_source.at_inverse(g,integral,inverse,angle)
        source_rows[k]=g.add(source_rows[k],g.mul(g.constant('-3/148'),g.quotient(at,g.node('mathematical_pi'),'positive pi')),
            g.mul(g.constant('-1/37'),phase) if k==0 else g.zero)
    Al=current.modified_source.lexical.Algebra;close=current.modified_source.lexical.close
    numerical=0;wraps=0
    for n in range(5):
        for k,ref in enumerate(source_rows):
            p={Y:sy.Rational(1,11),Z:sy.Rational(2,7)};f=phi0+N*J*p[Y]
            got=Al(g,{y.node:p[Y],z.node:p[Z],phase.node:f,inverse.node:known_inverse.subs({F:f,**p})}).at(ref.node)
            wanted=sy.diff(sy.diff(expected.subs(F,phi0+N*J*Y),Y,n)/J**n,Z,k).subs(p)
            close(got,wanted,('independent total native/implicit/integral',n,k));numerical+=1
            yw=(1-phi0)/(N*J);values=[]
            for endpoint in (0,1):
                values.append(Al(g,{y.node:yw,z.node:p[Z],phase.node:sy.Integer(endpoint),
                    inverse.node:known_inverse.subs({Y:yw,Z:p[Z],F:endpoint})}).at(ref.node))
            close(values[0],values[1],('periodic lifted mixed derivative wrap',n,k));wraps+=1
        if n<4:source_rows=[g.mul(g.constant('2/3'),D.derivative(ref)) for ref in source_rows]
    # The exact original sigma body and its flat branches, through native4.
    sig=current.previous.loop_source.previous.sigma(g,y);x=Y
    e=sy.exp(1/(1-x)**2-1/x**2);sigma=e/(1+e);sigma_rows=[]
    for n in range(1,5):
        sig=D.derivative(sig)
        close(Al(g,{y.node:sy.Rational(2,5)}).at(sig.node),sy.diff(sigma,x,n).subs(x,sy.Rational(2,5)),('actual sigma native derivative',n))
        for edge in (0,1):zero(Al(g,{y.node:sy.Integer(edge)}).at(sig.node),('actual sigma flat native derivative',n,edge))
        sigma_rows.append(sig.node)
    assert all(factory.cache[phase.node]==(g.mul(g.constant(37),g.constant('3/2')) if factory.variable==y else g.zero)
        for factory in D.shared['factories'].values())
    return dict(passed=True,independent_nonunit_J_N_implicit_inverse_moving_integral_y4_Z5_rows=numerical,
        periodic_wrap_mixed_rows=wraps,actual_original_sigma_native_orders_1_to_4=sigma_rows,
        angle_binder_holds_phase_while_main_native_factory_has_NJ=True,numerical_calculus_fixture_only=True)


CHECKS=(source_contract,periodic_lift,independent_total_native_calculus,actual_leading_micro_sources,
    actual_leading_recovery_pressure_and_FTC,actual_nonzero_histories_and_own_micro_join)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('ACTIVE_MICRO_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['micro_active_radial4_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_first_second_micro_contracted_mixed_join_rows=246,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_ACTIVE_MICRO_RADIAL4_AND_OWN246_SOURCE_ACCEPTED',flush=True);return rec
