"""Scoped six-chart source, physical geometry and Duhamel/pressure checks."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as current
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions_check as prior

class Algebra(prior.Algebra):
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i]
        if row['operation']=='exact_original_source_flat_sigma':
            x=self.at(row['argument'])
            if x.free_symbols:
                self.cache[i]=sy.Function('original_sigma')(x);return self.cache[i]
        return super().at(i)


lexical=prior.lexical


def jet(g,variables,body,order):return prior.jet(g,variables,body,order)


def actual_source(owner):
    g=owner.g;r=current.rows;free={i:set(v) for i,v in owner.leading.free_source_seeds.items()};count=0
    bg_in,corr,pressure=owner.initial_background,owner.initial_correction,owner.initial_integrals
    assert owner.patch.acceptance_loaded and tuple(owner.background)==tuple(owner.charts)==current.CHARTS
    assert owner.initial_correction is owner.patch.Rh_correction and owner.initial_integrals is owner.patch.Rh['integrals']
    assert owner.initial is owner.patch.initial and owner.leading_Rh is owner.leading.endpoints['Rh']
    assert all(owner.initial_background[k] is owner.leading_Rh['fields'][k] for k in current.RATES)
    for name in ('N','eta','dstar','sc','P0','Pstar','phase_origin'):assert getattr(owner,name) is getattr(owner.patch,name)
    assert owner.bindings['original_leading']==current.original_leading.source_contracts()
    assert owner.bindings['original_outer_frontends']==current.original_outer.FRONTEND_BINDINGS
    assert owner.bindings['original_geometry']==current.parent.active.current.ast_binding(current.original_geometry.exact_geometry)
    assert owner.mu_proof['eta_proof']==g.nodes[owner.eta.node]['canonical_parameter_proof']
    assert owner.mu_proof['positive_mu_and_power_2mu_ge_eta_and_Rp_reservation']
    for name in current.CHARTS:
        bg=owner.background[name];chart=owner.charts[name];coordinate=bg['phase'];native=g.nodes[coordinate.node]['name']
        assert bg['incoming'] is bg_in and chart['incoming'] is corr and chart['incoming_integrals'] is pressure
        assert chart['source_chart'] is bg and chart['roots']['p2'] is bg['generic']['p2']
        assert chart['loop']['original_collar_upper'] is None
        assert chart['actual_phase']==g.unary('fractional_part',g.mul(owner.N,g.sub(g.unary('log',
            g.quotient(bg['radius'],owner.Ra,'same original positive actual R/Ra')),owner.phase_origin)))
        assert lexical.free_variables(g,chart['actual_phase'],free)=={native}
        for group,order in ((bg['fields'].values(),6),(bg['radial'].values(),6),(bg['generic'].values(),5),
            (chart['history'].values(),5),(chart['modified_fields'].values(),5),
            (chart['finite_N']['five_signed_density_increments'].values(),5),(chart['generic'].values(),4),
            (chart['pressure'].values(),3),(chart['cumulative_integrals'].values(),3)):
            for q in group:
                assert len(r(q))==order+1
                for ref in r(q):assert lexical.free_variables(g,ref,free)<={'Z',native},(name,ref.node)
                count+=order+1
        D=owner.derivatives[name];assert D.cache[chart['actual_phase'].node]==g.zero
        for key,q in bg['fields'].items():
            for ref,dr in zip(r(q),r(chart['radial_background'][key])):
                assert D.cache[ref.node]==g.mul(bg['physical_measure'],dr)
        for key in ('E_y','V_y','m_y','slow_A_y','slow_B_y'):
            assert len(r(chart['radial'][key]))==6
            for ref in r(chart['radial'][key]):assert lexical.free_variables(g,ref,free)<={'Z',native}
        for key,q in chart['physical'].items():assert len(r(q))==(5 if key in ('Q','Q_y') else 6)
        for group,leading in ((bg['history_contributions'],True),(chart['contributions'],False)):
            for key,q in group.items():
                for n,ref in enumerate(r(q)):
                    row=g.nodes[ref.node]
                    assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==n
                    assert row['lower']==bg['domain'][0].node and row['upper']==coordinate.node
                    assert row['physical_Jacobian']==bg['physical_measure'].node
                    assert row['physical_log_offset']==bg['log_offset'].node
                    assert row['leading_background']==leading and row['finite_N_correction']==(not leading)
                    assert row['predecessor_memory_never_reset'] and row['original_rate']==str(current.RATES[key])
        assert bg['memory']['p']==chart['memory']['p']==g.one
        assert r(bg['fields']['P0'])==r(owner.P0)[:7] and r(chart['modified_fields']['P0'])==r(owner.P0)[:6]
        for key,q in chart['local_integrals'].items():
            for n,ref in enumerate(r(q)):
                row=g.nodes[ref.node]
                assert row['lower']==bg['domain'][0].node and row['upper']==coordinate.node
                assert row['physical_radius_power']==(1 if key=='S0' else 2)
                assert row['physical_Jacobian_applied_once'] and row['nonzero_preceding_S0_S1_retained']
        bg_in=owner.background_endpoints[name]['history'];corr=owner.endpoints[name]['correction'];pressure=owner.endpoints[name]['integrals']
    assert owner.Rc is owner.endpoints['O3_power'] and owner.Rc_background is owner.background_endpoints['O3_power']
    return dict(passed=True,actual_six_chart_source_history_generic_pressure_rows=count,
        true_Rh_leading_and_selected_N_correction_and_pressure_sources=True,
        every_nonzero_background_correction_and_pressure_predecessor_retained=True,
        genuine_C6_background_C5_generic_and_finite_N_C5_rows=True,
        same_original_N_global_phase_P0_and_complete_native_derivative_seeds=True,
        physical_Jacobian_and_radius_pressure_weights_once=True,
        terminal_Rc_and_compact_then_quiet_Rp_steps_not_conflated=True,
        numeric_mixed_norms_compact_controls_preRp_and_global_goals_not_claimed=True)


def exact_geometry(owner):
    g=owner.g;Rm=sy.Rational(7,3);X=sy.Symbol('X',real=True);Y=sy.exp(40)
    offsets=dict(Rh_reference=6+X,O2_slope=6+X,O2_axial=6+sy.exp(40*X),O2_buffer=6+Y+X,
        O3_transition=17+Y+X,O3_power=18+Y+X)
    expected_domain=dict(Rh_reference=(-5,0),O2_slope=(0,1),O2_axial=(0,1),O2_buffer=(0,11),O3_transition=(0,1),O3_power=(0,2))
    ends=[]
    for name,bg in owner.background.items():
        alg=Algebra(g,{owner.Rm.node:Rm,bg['phase'].node:X})
        assert sy.simplify(alg.at(bg['radius'].node)-Rm*sy.exp(offsets[name]))==0,name
        assert sy.simplify(alg.at(bg['physical_measure'].node)-sy.diff(offsets[name],X))==0,name
        assert tuple(alg.at(q.node) for q in bg['domain'])==expected_domain[name]
        ends.append((offsets[name].subs(X,expected_domain[name][0]),offsets[name].subs(X,expected_domain[name][1])))
    assert ends[0][0]==1 and all(ends[i][1]==ends[i+1][0] for i in range(5))
    assert ends[-1][1]==Y+20
    alg=Algebra(g,{})
    assert sy.simplify(alg.at(owner.logmu.node)-(sy.log(sy.Rational(1,1000))-4*(Y+11)))==0
    assert sy.simplify(alg.at(owner.Tw.node)+60*alg.at(owner.logmu.node))==0
    return dict(passed=True,six_original_maps_domains_Jacobians_and_five_radius_seams=True,
        O2_axial_Jacobian_40_exp_40x_retained=True,
        original_Rc_Rm_exp_exp40_plus20_and_exact_mu_Tw_recipe=True,
        Rp_requires_original_compact_Rc_to_2Rc_then_quiet_Tw_minus2_minus_log2=True)


def original_seed_theorem(owner):
    """Bind function identities, rather than reinterpret a C3 range as C6."""
    identity=owner.seed_identity;leading=current.parent.linear.read(current.previous.previous.RECEIPT)
    Rh=current.parent.linear.read(current.original_Rh.RECEIPT)
    assert identity['actual_Rh_fields'] is owner.leading_Rh['fields']
    assert leading['all_passed'] and leading['same_exact_original_unique_branch_theorem_bound']
    assert leading['checks']['branch_theorem']['same_exact_original_unique_branch_theorem_bound']
    assert leading['checks']['patch_matching']['passed']
    assert Rh['all_passed'] and Rh[current.original_Rh.GATES[0]]
    assert Rh['source_family']==leading['source_family']==owner.family
    assert Rh['fixed_leading_input_N0']==owner.leading.fixed_N0
    assert Rh['current_full_weight_map_and_source_AST_checked'] and Rh['interval_overlap_not_used_as_function_identity']
    assert identity['original_full_weight_function_identity']==current.original_Rh.exact_join_proof()
    assert owner.bindings['original_Rh_leading_identity']==current.original_Rh.source_bindings()
    assert identity['original_normalized_history_ODE_binding']==current.parent.active.current.ast_binding(current.original_leading.own_first_y)
    # Verify the normalized reference initial functions satisfy the actual
    # defining moment ODE. Its identity holds for all Z, so differentiation
    # preserves it beyond the old packet's C3 representation.
    g=current.source.FunctionTransportGraph();Y,Z=sy.symbols('reference_y Z');P=sy.Symbol('P',positive=True)
    variables={'reference_y':g.symbol('reference_y'),'Z':g.symbol('Z'),'P':g.symbol('P')}
    E=sy.exp(Y/10)/(1+Z**2);V=4*Z/P
    H=dict(m=V,h=5*E/8,k=5*E*V/8,e=V**2-5*E**2/12,p=5*E**2/2)
    make=lambda body:current.original_leading.target.C3Function(*current.rows(jet(g,variables,sy.sympify(body),3)))
    got=current.original_leading.own_first_y(g,make(E),make(V),{k:make(v) for k,v in H.items()})
    alg=Algebra(g,{variables['reference_y'].node:Y,variables['Z'].node:Z,variables['P'].node:P})
    for key,q in got.items():
        for n,ref in enumerate(current.original_leading.target.rows(q)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(H[key],Y,Z,n))==0,(key,n)
    assert identity['propagated_original_histories_identified_by_linear_ODE_uniqueness']
    return dict(passed=True,actual_Rh_C6_endpoint_bound_to_same_exact_original_unique_branch=True,
        original_full_support_leading_function_identity_not_derivative_range_overlap=True,
        original_normalized_reference_moment_ODE_symbolic_rows=20,
        equal_E_V_and_inlets_identify_all_propagated_histories_by_linear_ODE_uniqueness=True,
        C6_derivatives_follow_equal_smooth_functions_and_predecessor_implicit_theorem=True,
        old_C1_correction_and_C3_ranges_not_used_as_selected_N_C5_or_C6_fields=True,
        quantitative_high_order_norm_admission_still_open=True)


def actual_profile_calculus(owner):
    g=owner.g;r=current.rows;Z=sy.Symbol('Z');X=sy.Symbol('X',real=True);point=sy.Rational(3,8);mu=sy.Rational(1,13)
    A=2+Z+Z**2/17+Z**6/19;B=3+Z**2/23+Z**6/29;count=0
    sigma=lambda x:sy.Function('original_sigma')(x)
    prime=lambda x:sy.Function('original_sigma_prime')(x)
    for name,bg in owner.background.items():
        bindings={bg['phase'].node:X,g.symbol('Z').node:point,owner.mu.node:mu}
        incoming=owner.leading_Rh if name=='Rh_reference' else owner.background_endpoints[current.CHARTS[current.CHARTS.index(name)-1]]
        for n,ref in enumerate(r(incoming['fields']['E'])):bindings[ref.node]=sy.diff(A,Z,n).subs(Z,point)
        for n,ref in enumerate(r(owner.leading_Rh['fields']['V'])):bindings[ref.node]=sy.diff(B,Z,n).subs(Z,point)
        J=sy.Symbol('J_'+name,real=True)
        for i,row in enumerate(g.nodes):
            if row['operation']=='definite_integral' and row.get('exact_original_sigma_integral') and row['upper']==bg['phase'].node:
                bindings[i]=J
        if name=='Rh_reference':E=A*sy.exp((X+5)/10);V=B;Ey=E/10;Vy=0
        elif name=='O2_slope':E=A*sy.exp(X/10-3*J/5);V=B;Ey=E*(sy.Rational(1,10)-3*sigma(X)/5);Vy=0
        elif name=='O2_axial':E=A*sy.exp(-(sy.exp(40*X)-1)/2);V=B*sigma(1-X);Ey=-E/2;Vy=-B*prime(1-X)/(40*sy.exp(40*X))
        elif name=='O2_buffer':E=A*sy.exp(-X/2);V=0;Ey=-E/2;Vy=0
        elif name=='O3_transition':E=A*sy.exp(-X/2-mu*J);V=0;Ey=-(sy.Rational(1,2)+mu*sigma(X))*E;Vy=0
        else:E=A*sy.exp(-(sy.Rational(1,2)+mu)*X);V=0;Ey=-(sy.Rational(1,2)+mu)*E;Vy=0
        alg=Algebra(g,bindings)
        for key,q,body in (('E',bg['fields']['E'],E),('V',bg['fields']['V'],V),('E_y',bg['radial']['E_y'],Ey),('V_y',bg['radial']['V_y'],Vy)):
            for n,ref in enumerate(r(q)):
                assert sy.simplify(alg.at(ref.node)-sy.diff(sy.sympify(body),Z,n).subs(Z,point))==0,(name,key,n);count+=1
    return dict(passed=True,actual_profile_and_physical_first_radial_C6_rows=count,
        original_reference_slope_turnoff_buffer_transition_power_recipes=True,
        nonzero_sixth_amplitude_and_axial_rows_and_O2_physical_chain_rule=True,
        independent_sigma_integral_symbols_are_recipe_fixture_not_field_values=True,
        actual_mu_not_replaced_by_fixture_parameter=True,mathematical_fixture_only=True)


def exponential_transport_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_exp_outer_x');X=sy.Symbol('fixture_exp_outer_x',real=True);Z=sy.Symbol('Z')
    variables={'fixture_exp_outer_x':x,'Z':g.symbol('Z')};offset=g.unary('exp',g.mul(g.constant(2),x));jac=g.mul(g.constant(2),offset)
    for order in (5,6):
        density={key:1+Z+Z**order/(17+i) for i,key in enumerate(current.RATES)}
        incoming={key:2+Z**2/(23+i)+Z**order/(31+i) for i,key in enumerate(current.RATES)}
        make=lambda body:jet(g,variables,body,order)
        out,_,memory=current.transport_geometry(g,{k:make(v) for k,v in density.items()},
            {k:make(v) for k,v in incoming.items()},x,g.zero,offset,jac,'fixture_exp',order)
        alg=Algebra(g,{x.node:X,g.symbol('Z').node:Z});width=sy.exp(2*X)-1
        for key,rate in current.RATES.items():
            r=sy.Rational(rate);decay=sy.exp(-r*width);mass=width if rate==0 else (1-decay)/r
            want=incoming[key]*decay+density[key]*mass
            for n,ref in enumerate(current.rows(out[key])):
                assert sy.simplify(alg.at(ref.node)-sy.diff(want,Z,n))==0,(order,key,n)
        assert memory['p']==g.one
    return dict(passed=True,independent_exponential_measure_background_C6_and_correction_C5_rows=65,
        exact_integrated_log_width_and_nonzero_inlet_and_pressure_rate_zero=True,
        constant_Jacobian_shortcut_would_fail=True,mathematical_fixture_only=True)


def actual_quiet_power(owner):
    g=owner.g;r=current.rows;chart=owner.charts['O3_power'];Z=sy.Symbol('Z');X=sy.Symbol('quiet_x',real=True)
    mu=sy.Rational(1,13);point=sy.Rational(3,8);E=(2+Z+Z**5/17)*sy.exp(-(sy.Rational(1,2)+mu)*X)
    bindings={chart['coordinate'].node:X,g.symbol('Z').node:point,owner.mu.node:mu,owner.eta.node:sy.Rational(1,100),
              owner.N.node:sy.Integer(37),chart['actual_phase'].node:sy.Rational(2,5)}
    for n,ref in enumerate(r(chart['source_chart']['fields']['E'])):bindings[ref.node]=sy.diff(E,Z,n).subs(Z,point)
    for j,q in enumerate(chart['incoming'].values()):
        for n,ref in enumerate(r(q)):bindings[ref.node]=sy.diff(5+j+Z**5/(23+j),Z,n).subs(Z,point)
    alg=Algebra(g,bindings)
    for key,value in (('a',2+2*mu),('b',0),('Delta',2*mu)):
        for n,ref in enumerate(r(chart['roots'][key])):
            assert sy.simplify(alg.at(ref.node)-(value if n==0 else 0))==0,('quiet root',key,n)
    for key in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi'):
        for ref in r(chart['loop'][key]):assert alg.at(ref.node)==0,('quiet loop',key)
    for key,q in chart['finite_N']['five_signed_density_increments'].items():
        for ref in r(q):assert alg.at(ref.node)==0,('quiet density',key)
    for key,rate in current.RATES.items():
        for n,ref in enumerate(r(chart['history'][key])):
            incoming=alg.at(r(chart['incoming'][key])[n].node)
            assert sy.simplify(alg.at(ref.node)-incoming*sy.exp(-sy.Rational(rate)*X))==0,('actual quiet memory',key,n)
    return dict(passed=True,actual_constant_a_b_Delta_C5_and_flat_loop_and_density_rows=True,
        actual_nonzero_correction_C5_quiet_semigroup_rows=30,
        rate_zero_pressure_inlet_retained=True,canonical_source_mu_eta_guard_separate_from_fixture=True,
        mathematical_fixture_only=True)


def exponential_pressure_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_outer_pressure_x');X=sy.Symbol('fixture_outer_pressure_x',real=True)
    Z,Y=sy.symbols('Z Y');variables={'fixture_outer_pressure_x':x,'Z':g.symbol('Z')}
    P,Rm,de=sy.Rational(3,2),sy.Rational(7,3),sy.Rational(1,5);native=sy.exp(2*X);jac=2*native
    m=2+Z+Z**5/11+(3+Z**2)*Y+(1+Z**3)*Y**2;U=m+sy.diff(m,Y);Uy=sy.diff(U,Y)
    make=lambda body,order=5:jet(g,variables,body.subs(Y,native),order)
    physical=current.radial.physical_pressure_fields(g,{'V':make(U/P),'m':make(m/P)},
        {'V_y':make(Uy/P)},g.constant(de),g.constant(P))
    S=current.radial.reduced_pressure_S(g,physical,g.constant(de))
    R=current.source.expression(g,Rm*sy.exp(native),variables);J=current.source.expression(g,jac,variables)
    incoming={k:make(v,3) for k,v in {'S0':7+Z**3/19,'S1':11+Z**3/23}.items()}
    initial={k:make(v,3) for k,v in {'Q_plus_Q_y':13+Z**3/29,'R_Q_y':17+Z**3/31}.items()}
    _,_,_,pressure=current.pressure_geometry(g,S,physical,R,x,g.zero,J,incoming,initial,'fixture_exp')
    L=1-de*Z**2;d=1-Z**2;q=(2*Z*U-(1-de)*Z*m-d*sy.diff(m,Z))/L;qy=sy.diff(q,Y);qyy=sy.diff(q,Y,2)
    ss=(q+(1-de)*Z*sy.diff(q,Z)/2+qy)/L+q*(q/2+qy)+U*(-2*Z*q+d*sy.diff(q,Z)-2*Z*qy)/L
    N1p=-ss/2+(qy+qyy)/(Rm*sy.exp(Y));point=sy.Rational(3,8);alg=Algebra(g,{x.node:X,g.symbol('Z').node:point});count=0
    for key,multiplier in (('m13',Rm*sy.exp(native)*jac),('m15',-Rm**2*sy.exp(2*native)*jac)):
        for n,ref in enumerate(current.rows(pressure[key])):
            value=alg.at(ref.node);want=sy.diff(N1p,Z,n).subs({Y:native,Z:point})*multiplier
            assert sy.simplify(sy.diff(value,X)-want)==0,(key,n,'exponential measure');count+=1
    return dict(passed=True,independent_direct_N1p_exponential_pressure_C3_rows=count,
        exponential_Jacobian_and_physical_R_Rsquared_and_m15_sign_once=True,
        inherited_integrals_and_original_global_boundary_present=True,mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['true_outer_C6_C5_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_source',actual_source),('geometry',exact_geometry),('original_seed_theorem',original_seed_theorem),('actual_profiles',actual_profile_calculus),
                    ('actual_quiet_power',actual_quiet_power)):
        checks[name]=fn(owner);print('TRUE_OUTER_'+name.upper()+'_CHECKED',flush=True)
    checks['exponential_transport']=exponential_transport_calculus();print('TRUE_OUTER_EXPONENTIAL_TRANSPORT_CHECKED',flush=True)
    checks['exponential_pressure']=exponential_pressure_calculus();print('TRUE_OUTER_EXPONENTIAL_PRESSURE_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_SIX_OUTER_C6_C5_AND_RC_PRESSURE_C3_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses checked actual Rh sources')
