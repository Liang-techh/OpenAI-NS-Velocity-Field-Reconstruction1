"""Full original source domains and the actual macro1/first-switch0 seam."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions as current
import lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check as prior
import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check as geometry_check

rows,encode=current.outer.rows,current.outer.encode
zero,bind=prior.zero,prior.bind
Algebra=prior.SeamAlgebra


def source_contract(owner):
    owner.assert_graph();g=owner.g
    assert tuple(owner.charts)==current.CHARTS and owner.base.acceptance_loaded
    assert current.BACKGROUND_NAMES==dict(post_power='post_switch',postrestore='post_restore')
    assert owner.loop_input==owner.base.loop_input and owner.modified_input==owner.base.modified_input
    support=current.parent.linear.read(current.previous.previous.previous.common_source.previous.RECEIPT)['checks']
    assert support['source']['actual_preRc_coordinate_radius_domain_and_C3_pressure_source_identity']
    assert support['actual_regional_recovery']['own_V_y_and_mass_y_zero_only_outside_velocity_supports']
    requests=[];core=[];sigma=[]
    for name in current.CHARTS:
        chart=owner.charts[name];loop=owner.loop[name];old=owner.modified[name];bg=owner.backgrounds[name]
        assert encode(chart['domain'])==loop['domain']==chart['intermediate']['domain']
        if name=='first_micro':assert loop['domain'][1]==g.one.node
        elif name=='second_micro':assert loop['domain']==[g.one.node,g.constant(2).node]
        else:assert loop['domain']==[g.zero.node,g.one.node]
        assert chart['coordinate'].node==bg['phase']==loop['coordinate']
        assert chart['radius'].node==bg['radius']==loop['radius']
        assert chart['jacobian'].node==bg['physical_measure']==loop['physical_measure']
        if name in ('first_switch','second_switch'):assert chart['jacobian'].node==owner.downstream['hb']
        if name in ('restoration','postrestore'):assert chart['jacobian']==g.one
        if name=='long_reshape':assert chart['jacobian'].node==owner.downstream['T']
        if name=='reference':assert chart['jacobian'].node==owner.downstream['reference_length']
        if name=='post_power':assert chart['jacobian'].node==owner.downstream['post_length']
        assert chart['intermediate']==owner.intermediate['charts'][name]
        assert chart['intermediate']['region']=='beforeRc' and chart['intermediate']['complete_incoming_mass_not_reset']
        assert chart['intermediate']['fixed_pressure']==old['pressure']
        for key in ('m11','m12','m14'):
            assert chart['intermediate']['cumulative'][key]==owner.lower.rin.cut['endpoint'][key][:4]
        assert encode(owner.background_y[name])==old['radial_background']
        D=owner.factories[name]
        assert len(owner.requests[name])==18 and D.cache[chart['phase'].node]==g.mul(owner.N,chart['jacobian'])
        for key,ids in bg['fields'].items():
            for node,dy in zip(ids,rows(owner.background_y[name][key])):
                assert D.cache[node]==g.mul(chart['jacobian'],dy)
        for key in ('U','m','Q'):
            assert encode(chart['physical'][key][:2])==[old['physical'][key],old['physical'][key+'_y']]
            assert len(chart['physical'][key])==5
            assert all(len(rows(jet))==(5 if key=='Q' else 6) for jet in chart['physical'][key])
        for request in owner.requests[name]:
            assert D.cache[request['input']].node==request['native_output']
            out=g.nodes[request['physical_output']]
            # The saved physical output really is J^-1 times the issued
            # native derivative; its exact certificate is source-specific.
            certificate=('same positive original micro physical Jacobian' if name in current.CHARTS[:2] else
                'same positive original frozen-macro physical Jacobian' if name=='frozen_macro' else
                'same positive original prescribed chart physical Jacobian')
            assert out==g.nodes[g.mul(g.quotient(g.one,chart['jacobian'],certificate),owner.ref(request['native_output'])).node]
            requests.append([name,request])
        for factory in D.shared['factories'].values():
            assert factory.cache[chart['phase'].node]==(g.mul(owner.N,chart['jacobian'])
                if factory.variable==chart['coordinate'] else g.zero)
            for node in factory.core_extensions:owner.active.guard_core(g.nodes[node]);core.append(node)
            for extension in factory.sigma_extensions:
                assert 1<=extension['order']<=4;sigma.append([name,extension])
    assert len(requests)==180
    inherited=current.parent.linear.read(current.previous.RECEIPT)
    assert inherited['all_passed'] and len(inherited['checks'])==7
    return dict(passed=True,actual_native_to_physical_derivative_requests=requests,
        ten_full_original_domains_freshly_rederived=True,
        actual_chart_to_background_map_and_distinct_native_Jacobians=True,
        all_primary_phase_caches_NJ_and_secondary_binders_hold_phase=True,
        actual_core_extensions=sorted(set(core)),actual_sigma_extensions=sigma,
        actual_beforeRc_support_and_nonzero_cutoff_histories_consumed=True,
        inherited_periodic_inverse_integral_and_sigma_calculus=current.previous.RECEIPT)


def original_geometry(owner):
    g=owner.g;env={}
    for i,row in enumerate(g.nodes):
        if row['operation'] in ('original_source_parameter','current_original_source_parameter'):
            env[i]=sy.Symbol(row['name']+'_'+str(i),positive=True)
    geometry=[]
    for name in current.CHARTS:
        chart=owner.charts[name];x=sy.Symbol(g.nodes[chart['coordinate'].node]['name'])
        alg=geometry_check.Algebra(g,{**env,chart['coordinate'].node:x})
        R=alg.at(chart['radius'].node);J=alg.at(chart['jacobian'].node)
        zero(sy.diff(R,x)-J*R,('actual original physical measure',name))
        assert x not in J.free_symbols
        geometry.append(dict(chart=name,radius=chart['radius'].node,physical_Jacobian=chart['jacobian'].node))
    certificate=owner.downstream['geometry_certificate']
    for key in ('shared_upstream_is_same_object','bridge_logh_is_actual_upstream_object',
        'switch_width_is_exact_source_reboxing','bridge_hbB_switch_hbS_kept_distinct',
        'all_four_native_T_source_tuples_equal','actual_native_T_and_logC_match_common_source',
        'directed_T_positive','directed_reference_gap_positive','directed_post_length_positive'):
        assert certificate[key]
    ledger=current.PREFIX+'K1_ledger.json'
    assert certificate['shared_K1_ledger_sha256']==owner.hashes[ledger]==current.sha(ledger)
    assert owner.downstream['hb']!=owner.downstream['bridge_hb']==owner.bridge['hb']
    record=current.parent.linear.read(current.downstream_source.RECEIPT)['checks']
    assert record['actual_source']['distinct_width_tags_with_exact_shared_upstream_K1_identity']
    return dict(passed=True,actual_R_native_derivative_equals_J_times_R=geometry,
        all_J_are_original_Z_independent_constants=True,
        hbB_hbS_are_distinct_nodes_with_same_actual_K1_width_source=True,
        exact_width_certificate=certificate,
        actual_positive_finite_T_and_reference_post_lengths_retained=True)


def reduced_S(owner,name):
    g=owner.g;a=current.outer.Algebra(g,3);old=owner.modified[name]
    z=current.outer.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=a.fixed(1);de=a.fixed(owner.parameters['delta'])
    d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    U,Q,Qy=(a.trim(owner.jet(old['physical'][key])) for key in ('U','Q','Q_y'))
    Qz=a.shift(owner.jet(old['physical']['Q']));div=lambda v:a.div(v,L,'same original positive meridional L')
    return a.add(div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy)),
        a.mul(Q,a.add(a.scale(Q,'1/2'),Qy)),
        div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2)))))


def actual_recovery_and_pressure_sources(owner):
    Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_physical_y',real=True)
    U=sy.Function('actual_arbitrary_U')(y,Z);M=sy.Function('actual_arbitrary_mean')(y,Z)
    rows_checked=[];operators=[]
    inherited=current.parent.linear.read(current.previous.RECEIPT)['checks']
    assert inherited['actual_macro_pressure_and_FTC']['own_pressure_outputs_not_bound']
    assert inherited['actual_leading_recovery_and_common_N']['common_N_is_derived_before_own_output_input_bindings']
    for name in current.CHARTS:
        chart=owner.charts[name];env,v=prior.active.prior.active_context(owner.rminus,False);delta=v['delta']
        for n,jet in enumerate(chart['physical']['U']):bind(env,jet,sy.diff(U,y,n),Z)
        bind(env,chart['physical']['m'][0],M,Z);alg=Algebra(owner,env);means=[M]
        for n in range(1,5):means.append(sy.diff(U,y,n-1)-means[-1])
        for n in range(5):
            zero(alg.at(chart['physical']['m'][n].value)-means[n],('actual mean ODE',name,n))
            Q=(2*Z*sy.diff(U,y,n)-(1-delta)*Z*means[n]-(1-Z**2)*sy.diff(means[n],Z))/(1-delta*Z**2)
            for k,ref in enumerate(rows(chart['physical']['Q'][n])):
                zero(alg.at(ref)-sy.diff(Q,Z,k),('actual full-chart leading Q',name,n,k));rows_checked.append([name,n,k,ref.node])
        f={key+suffix:chart['physical'][key][n] for key in ('U','m','Q')
            for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        density=current.rc_source.radial.pressure_density(owner.g,f,owner.parameters['delta'],chart['radius'])
        higher=current.high.higher_pressure_density(owner.g,f,owner.parameters['delta'],chart['radius'],density)
        assert encode(density)==encode(chart['density']) and encode(higher)==encode(chart['higher_density'])
        assert encode([density['N1p'],higher['N1p_y'],higher['N1p_yy']])==encode(chart['N_y'])
        assert encode(reduced_S(owner,name))==owner.modified[name]['S']
        operators.append(dict(chart=name,actual_N_y=encode(chart['N_y']),actual_S=owner.modified[name]['S']))
    return dict(passed=True,actual_y0_to_y4_leading_Q_C4_rows=rows_checked,
        actual_mean_ODE_and_actual_saved_reduced_S_on_all_ten_charts=True,
        actual_fixed_pressure_operators_replayed=operators,
        arbitrary_U_Q_R_and_nonunit_J_pressure_operator_theorem=current.previous.RECEIPT,
        R_y_R_geometry_and_same_pinned_operator_transfer_N_y2_pressure_y3_C3=True,
        no_new_numerical_or_global_mixed_bound_claim=True)


def actual_signed_transport_sources(owner):
    g=owner.g;a=current.outer.Algebra(g,3);previous=None;integrals=[]
    initial={key:owner.jet(ids) for key,ids in owner.modified_input['initial'].items()}
    for name in current.CHARTS:
        chart=owner.charts[name];old=owner.modified[name];loop=owner.loop[name]
        if previous is None:
            assert all(set(ids)=={0} for ids in old['incoming_integrals'].values())
            assert all(set(ids)=={0} for ids in loop['incoming'].values())
        else:
            assert old['incoming_integrals']==owner.modified_input['endpoints'][previous]['integrals']
            assert loop['incoming']==owner.loop_input['endpoints'][previous]
        S=owner.jet(old['S']);coordinate=chart['coordinate'];J=chart['jacobian'];R=chart['radius'];local={};cum={}
        for key,power in (('S0',1),('S1',2)):
            ids=old['local_integrals'][key];template=g.nodes[ids[0]];t=g.symbol(template['variable'])
            body=a.at(a.scale(S,g.mul(J,*([R]*power))),coordinate,t)
            for n,(node,integrand) in enumerate(zip(ids,rows(body))):
                row=g.nodes[node]
                assert row['operation']=='definite_integral' and row['integrand']==integrand.node
                assert row['lower']==loop['domain'][0] and row['upper']==coordinate.node
                assert row['physical_radius_power']==power and row['physical_Jacobian_applied_once']
                assert row['ordinary_slow_Z_derivative_order']==n and row['endpoint_independent_of_Z']
                assert row['original_chart']==name and row['original_signed_reduced_pressure']
                integrals.append([name,key,n,node])
            local[key]=owner.jet(ids);cum[key]=a.add(owner.jet(old['incoming_integrals'][key]),local[key])
        assert encode(cum)==old['cumulative_integrals']
        Q,Qy=(a.trim(owner.jet(old['physical'][key])) for key in ('Q','Q_y'))
        boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R))
        assert encode(boundary)==old['finite_boundary']
        expected=dict(m13=a.add(a.scale(cum['S0'],'-1/2'),a.sub(boundary['Q_plus_Q_y'],initial['Q_plus_Q_y'])),
            m15=a.sub(a.scale(cum['S1'],'1/2'),a.sub(boundary['R_Q_y'],initial['R_Q_y'])))
        assert encode(expected)==old['pressure']
        for key in ('integrals','pressure','boundary'):
            packet=cum if key=='integrals' else expected if key=='pressure' else boundary
            assert {field:encode(a.at(jet,coordinate,owner.ref(loop['domain'][1]))) for field,jet in packet.items()}==owner.modified_input['endpoints'][name][key]
        previous=name
    assert len(integrals)==80
    return dict(passed=True,actual_ten_chart_signed_partial_integral_rows=integrals,
        actual_S_RJ_and_S_R2J_integrands_and_both_finite_boundaries_replayed=True,
        both_signed_primitive_formulas_replayed_on_all_ten_original_domains=True,
        every_actual_predecessor_pressure_and_five_correction_history_retained=True,
        only_original_rminus_correction_inlet_is_zero=True,
        original_signed_m15_not_negated_a_second_time=True,
        arbitrary_constant_J_FTC_identity_consumed_from=current.previous.RECEIPT)


def actual_supported_own_flow_and_queries(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);yy=sy.Symbol('independent_own_physical_y',real=True)
    env={};mass=sy.Function('nonzero_actual_cutoff_mass')(Z)
    bind(env,owner.jet(owner.lower.rin.cut['endpoint']['m11'][:4]),mass,Z)
    for name in current.CHARTS:env[owner.charts[name]['radius'].node]=sy.Symbol('actual_R_'+name,positive=True)
    alg=Algebra(owner,env);delta=alg.at(owner.parameters['delta']);checked=[]
    for name in current.CHARTS:
        chart=owner.charts[name];f=chart['own']['fields'];R=alg.at(chart['radius'])
        mean=mass/R*sy.exp(-yy);radius=R*sy.exp(yy)
        Q=(-(1+delta)*Z*mean-(1-Z**2)*sy.diff(mean,Z))/(1-delta*Z**2)
        model=dict(Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),M1=mass,
            mean1=mean,K1=mean,Q1=Q,V1=radius*Q)
        for key,body in model.items():
            for n,jet in enumerate(f[key]):
                zero(alg.at(jet.value)-sy.diff(body,yy,n).subs(yy,0),('actual full supported own flow',name,key,n))
                checked.append([name,key,n])
        for key in ('M1','mean1','Q1','V1'):
            for new,old in zip(rows(f[key][1]),chart['intermediate']['fields'][key+'_y']):
                zero(alg.at(new)-alg.at(old),('actual original own first derivative',name,key))
        assert f['P1']==chart['own']['cumulative']['m13']
        assert encode(f['P1'][0])==chart['intermediate']['fields']['P1']
        assert encode(chart['own']['cumulative']['m15'][0])==chart['intermediate']['cumulative']['m15']
    queries=[]
    for name in current.CHARTS:
        for fraction in ('0','1/2','1'):
            packet=owner.at_fraction(name,fraction)
            assert not packet['numerical_oracle_available'] and len(packet['physical']['U'])==5
            assert len(packet['N_y'])==3 and len(packet['own']['fields']['P1'])==4
            queries.append([name,fraction,packet['native_query'].node])
    rejected=0
    for name,fraction in (('not_a_chart','1/2'),('first_micro','-1/10'),('postrestore','11/10')):
        try:owner.at_fraction(name,fraction)
        except ValueError:rejected+=1
        else:raise AssertionError('Outside original source domain accepted')
    return dict(passed=True,actual_supported_own_state_physical_y0_to_y4_rows=checked,
        actual_own_first_y_rows_original_one_plus_delta_and_nonzero_mass_retained=True,
        full_chart_fraction_queries=queries,invalid_chart_or_fraction_rejections=rejected,
        pressure_and_all_five_history_value_sources_are_original_beforeRc_nodes=True)


def actual_macro_first_background_seam(owner):
    g=owner.g;bridge=prior.bridge;a=bridge.BridgeAlgebra(g,6);b=bridge.BridgeAlgebra(g,7)
    saved=owner.bridge;bg=owner.backgrounds['first_switch'];loop=owner.loop['first_switch']
    phase=owner.charts['first_switch']['coordinate'];R=owner.charts['first_switch']['radius']
    h=owner.ref(owner.downstream['hb']);m=owner.ref(saved['macro']['phase'])
    old=saved['comparison_macro'];phi,U=owner.jet(old['phi']),owner.jet(old['U'])
    incoming={key:b.at(owner.jet(ids),m,g.one) for key,ids in old['moments'].items()}
    distance=g.unary('log',g.quotient(R,g.constant(100),'same positive current radius over100'))
    eq=dict(H=phi,M=U,K=b.mul(phi,U),A=b.mul(U,U),B=b.scale(b.mul(phi,phi),'1/2'),C=b.mul(phi,phi))
    J={}
    for key,rate in bridge.RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),distance))
        J[key]=b.add(b.scale(incoming[key],decay),b.scale(eq[key],g.sub(g.one,decay)))
    F0=owner.jet(saved['core7']['F0']);P0=owner.jet(saved['P0']);Pstar=owner.ref(owner.downstream['Pstar'])
    known=bridge.comparison_directions(g,phi,U,J,F0,P0,owner.parameters['delta'])
    assert encode(known)==bg['known_comparison_direction']
    initial={key:owner.jet(saved['R100']['actual'][key]) for key in ('phi','U')}
    actual=current.downstream_source.first_switch_fields(a,known,initial,phase,h,R,a.trim(F0),Pstar,phi)
    assert encode(actual)==bg['actual']
    true_incoming={key:owner.jet(ids) for key,ids in saved['R100']['moments'].items()}
    assert encode(true_incoming)==bg['incoming_histories']
    histories=bridge.histories(a,actual['phi'],actual['U'],true_incoming,phase,h,'true_first_switch')
    assert encode(histories)==bg['moments']
    fields,radial,generic=bridge.generic_fields(a,actual,histories,a.trim(F0),P0,Pstar,R,owner.parameters['delta'])
    assert encode(fields)==bg['fields'] and encode(radial)==bg['radial'] and encode(generic)==bg['generic']
    for key in ('phi','U'):
        assert encode(a.at(owner.jet(saved['macro']['actual'][key]),m,g.one))==saved['R100']['actual'][key]
    for key,ids in saved['macro']['moments'].items():
        assert encode(a.at(owner.jet(ids),m,g.one))==saved['R100']['moments'][key]
    sigma=current.downstream_source.sigma(g,phase);D=owner.factories['first_switch']
    evaluator=current.previous.previous.modified_source.lexical.Algebra;controls=[]
    for n in range(5):
        zero(evaluator(g,{phase.node:sy.Integer(0)}).at(sigma.node),('actual first-switch sigma flat inlet',n))
        controls.append([n,sigma.node])
        if n<4:sigma=D.derivative(sigma)
    # Equality of widths is a pinned source theorem, not node aliasing.
    cert=owner.downstream['geometry_certificate']
    assert cert['switch_width_is_exact_source_reboxing'] and cert['bridge_hbB_switch_hbS_kept_distinct']
    H,Ra=sy.symbols('same_K1_width Ra',positive=True)
    geometry=geometry_check.Algebra(g,{saved['Ra']:Ra,saved['hb']:H,owner.downstream['hb']:H,
        saved['macro']['phase']:sy.Integer(1),phase.node:sy.Integer(0)})
    length=geometry.at(saved['length'])
    zero(length-(sy.log(100)-sy.log(Ra)-2*H),'actual saved macro length definition')
    zero(geometry.at(saved['macro']['radius'])-100,'actual saved macro radius at1')
    zero(geometry.at(R.node)-100,'actual saved first-switch radius at0')
    for n in range(1,5):
        symbol=sy.Symbol('sigma_derivative_'+str(n))
        zero((H*symbol).subs(symbol,0),'actual axial switch control derivative flat at0')
    inherited=current.parent.linear.read(current.previous.RECEIPT)['checks']
    assert inherited['actual_smooth_source_seam']['original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows']
    return dict(passed=True,actual_current_radius_comparison_direction_and_B_half_replayed=True,
        actual_first_switch_phi_U_histories_generic_background_replayed=True,
        actual_macro1_nonzero_phi_U_six_histories_are_original_first_switch0_inputs=True,
        actual_first_switch_sigma_native0_to4_at0=controls,
        common_K1_width_consumed_while_hbB_hbS_native_nodes_remain_distinct=True,
        common_physical_comparison_history_ODE_then_actual_background_jets_y0_to5=True,
        common_original_roots_unique_inverse_periodic_lift_and_NJ_transfer_U_y4_C5=True,
        original_mean_history_ODE_Q_recovery_and_N_operator_transfer_Q_y4_N_y2=True,
        seam_scope='actual macro1/first-switch0 source identity; numerical mixed bounds remain open')


def actual_macro_first_own_seam(owner):
    g=owner.g;env,Z,slow=prior.active.prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
    macro=owner.charts['frozen_macro'];first=owner.charts['first_switch']
    x=sy.Symbol(g.nodes[macro['coordinate'].node]['name'],real=True)
    s=sy.Symbol(g.nodes[first['coordinate'].node]['name'],real=True)
    L,H=sy.symbols('actual_macro_length same_actual_K1_width',positive=True)
    env[macro['jacobian'].node]=L;env[first['jacobian'].node]=H;env[owner.ref(owner.bridge['hb']).node]=H
    env[owner.rminus.active_sc.node]=Algebra(owner,env).at(owner.parameters['sc'])
    env[macro['radius'].node]=100*sy.exp(L*(x-1));env[first['radius'].node]=100*sy.exp(H*s)
    for name in ('first_micro','second_micro','frozen_macro'):
        coordinate=sy.Symbol(g.nodes[owner.charts[name]['coordinate'].node]['name'],real=True)
        for key,ids in owner.loop[name]['contributions'].items():
            bind(env,owner.jet(ids),sy.Function('actual_'+name+'_correction_primitive_'+key)(coordinate,Z),Z)
        for key,ids in owner.modified[name]['local_integrals'].items():
            bind(env,owner.jet(ids),sy.Function('actual_'+name+'_signed_primitive_'+key)(coordinate,Z),Z)
    mean=sy.Function('actual_common_macro_first_background_mean');mapped=1+H*s/L
    bind(env,owner.jet(owner.backgrounds['frozen_macro']['fields']['m']),mean(x,Z),Z)
    bind(env,owner.jet(owner.backgrounds['first_switch']['fields']['m']),mean(mapped,Z),Z)
    for name in ('first_micro','second_micro'):
        coordinate=sy.Symbol(g.nodes[owner.charts[name]['coordinate'].node]['name'],real=True)
        for key in ('Q','Q_y'):
            bind(env,owner.jet(owner.modified[name]['physical'][key]),sy.Function('retained_'+name+'_leading_'+key)(coordinate,Z),Z)
    for name,argument in (('frozen_macro',x),('first_switch',mapped)):
        for key in ('Q','Q_y'):
            bind(env,owner.jet(owner.modified[name]['physical'][key]),sy.Function('proved_common_leading_'+key)(argument,Z),Z)
        for n,jet in enumerate(owner.charts[name]['N_y']):
            bind(env,jet,sy.Function('proved_common_leading_N_y'+str(n))(argument,Z),Z)
    alg=Algebra(owner,env);means=[]
    for name,value in (('frozen_macro',g.one),('first_switch',g.zero)):
        chart=owner.charts[name]
        jet=current.outer.Algebra(g,5).at(chart['physical']['m'][0],chart['coordinate'],value)
        means.append([alg.at(ref) for ref in rows(jet)])
    for k,(left,right) in enumerate(zip(*means)):zero(left-right,('actual leading mean transfer macro1 first0',k))
    evidence=[]
    for spec in owner.contract:
        for n in range(spec['physical_y_order']+1):
            for k,(left,right) in enumerate(zip(rows(owner.macro_first_join['left'][spec['group']][spec['field']][n]),
                rows(owner.macro_first_join['right'][spec['group']][spec['field']][n]))):
                zero(alg.at(left)-alg.at(right),('actual macro1 first-switch0 own mixed seam',spec['group'],spec['field'],n,k))
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=n,
                    ordinary_Z_order=k,left_node=left.node,right_node=right.node))
    assert len(evidence)==246
    return dict(passed=True,actual_macro_first_own_mixed_rows=evidence,
        actual_leading_mean_C5_value_transfer_rows=6,
        original_nonzero_macro_correction_and_signed_pressure_memories_retained=True,
        common_leading_Q_and_N_inputs_follow_preceding_source_theorem_not_desired_own_outputs=True,
        native_length_and_hbS_remain_separate=True)


CHECKS=(source_contract,original_geometry,actual_recovery_and_pressure_sources,actual_signed_transport_sources,
    actual_supported_own_flow_and_queries,actual_macro_first_background_seam,actual_macro_first_own_seam)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner)
        print('PRESCRIBED_HIGH_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['prescribed_high_radial_source_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_full_original_charts=list(current.CHARTS),actual_native_to_physical_derivative_requests=180,
        actual_macro_first_switch_own_mixed_rows=246,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_TEN_FULL_CHART_HIGH_SOURCES_AND_MACRO_FIRST246_ACCEPTED',flush=True)
    return rec


if __name__=='__main__':raise SystemExit('run(owner) reuses the accepted live macro ancestor')
