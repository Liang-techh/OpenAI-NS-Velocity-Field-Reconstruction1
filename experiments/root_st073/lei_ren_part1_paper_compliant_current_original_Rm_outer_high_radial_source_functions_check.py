"""Actual variable-J high sources, exact beta calculus and original memories."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_Rm_outer_high_radial_source_functions as current
import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check as high_check
import lei_ren_part1_paper_compliant_current_original_annular_mixed_source_joins_check as annular_check

rows,encode=current.outer.rows,current.outer.encode
zero,bind,Algebra=high_check.zero,high_check.bind,high_check.Algebra
lexical=high_check.current.previous.previous.modified_source.lexical


def source_contract(owner):
    owner.assert_graph();g=owner.g;requests=[];beta=[];core=[]
    assert tuple(owner.charts)==current.CHARTS and owner.base.acceptance_loaded
    inherited=current.parent.linear.read(current.previous.RECEIPT)
    assert inherited['all_passed'] and inherited['actual_own_mixed_rows']==984
    assert owner.patch_input['N']==owner.outer_input['N']==owner.N.node
    assert owner.leading_input['branch_transfer']['actual_source_Z_domain']==['-1','1']
    # Every opaque control has the actual source dependency Z only. Native
    # zero derivatives follow from dependency, not a zero control value.
    controls=owner.leading_input['controls'];control_values=[ids[0] for ids in controls]
    for node in control_values:
        row=g.nodes[node]
        assert row['operation']=='original_Rm_unique_leading_implicit_control'
        assert row['source_family']==owner.family and owner.free_seeds[node]=={'Z'}
    for name in current.CHARTS:
        chart=owner.charts[name];bg=owner.backgrounds[name];old=owner.modified[name];loop=owner.loop[name];D=owner.factories[name]
        assert encode(chart['domain'])==loop['domain']==chart['intermediate']['domain']
        assert chart['coordinate'].node==bg['phase']==loop['coordinate']
        assert chart['radius'].node==bg['radius']==loop['radius']
        assert chart['jacobian'].node==bg['physical_measure']==loop['physical_measure']
        assert encode(owner.background_y[name])==old['radial_background']
        assert chart['intermediate']['region']=='beforeRc' and chart['intermediate']['complete_incoming_mass_not_reset']
        assert chart['intermediate']['fixed_pressure']==old['pressure']
        for field in ('m11','m12','m14'):
            assert chart['intermediate']['cumulative'][field]==owner.lower.rin.cut['endpoint'][field][:4]
        for field,ids in bg['fields'].items():
            for node,dy in zip(ids,rows(owner.background_y[name][field])):
                assert D.cache[node]==g.mul(chart['jacobian'],dy)
        for node in control_values:
            assert D.shared['free'][node]=={'Z'} and D.derivative(owner.ref(node))==g.zero
        for field in ('U','m','Q'):
            assert encode(chart['physical'][field][:2])==[old['physical'][field],old['physical'][field+'_y']]
            assert len(chart['physical'][field])==5
            assert all(len(rows(jet))==(5 if field=='Q' else 6) for jet in chart['physical'][field])
        assert len(owner.requests[name])==18
        assert g.nodes[chart['coordinate'].node]['name'] in lexical.free_variables(g,chart['phase'],D.shared['free'])
        for request in owner.requests[name]:
            assert D.cache[request['input']].node==request['native_output']
            assert g.mul(owner.ref(request['Jacobian_inverse']),owner.ref(request['native_output'])).node==request['physical_output']
            requests.append([name,request])
        for factory in D.shared['factories'].values():
            assert factory.cache[chart['phase'].node]==(g.mul(owner.N,chart['jacobian']) if factory.variable==chart['coordinate'] else g.zero)
            for node in factory.core_extensions:owner.base.active.guard_core(g.nodes[node]);core.append(node)
            for item in factory.beta_extensions:
                row=g.nodes[item['input']]
                assert row==g.nodes[current.patch_leading.raw_beta(g,owner.ref(item['argument'])).node]
                beta.append([name,item])
    assert len(requests)==108 and beta
    return dict(passed=True,actual_native_to_physical_derivative_requests=requests,
        actual_six_original_domains_and_first_y_sources_bound=True,
        original_opaque_unique_Rm_controls_are_Z_only_not_zero_values=True,
        actual_beta_rule_extensions=beta,actual_guarded_core_extensions=sorted(set(core)),
        every_main_phase_cache_NJ_and_secondary_binder_holds_phase=True,
        existing_O3_power_Rc_compact_quiet_and_combined_outer_not_rebuilt=True)


def original_geometry(owner):
    g=owner.g;p=owner.parameters
    env={p['logC'].node:sy.Symbol('actual_logC',real=True),p['logP'].node:sy.Symbol('actual_logP',real=True)}
    # Bind Rm as the actual positive physical inlet. Geometry identities
    # below retain all native offsets; no interval value defines Rm.
    Rm=owner.leading_input['inlet']['radius'];env[Rm]=sy.Symbol('actual_positive_Rm',positive=True)
    evidence=[]
    for name in current.CHARTS:
        chart=owner.charts[name];x=sy.Symbol('native_'+name,positive=True)
        a=annular_check.SourceAlgebra(g,{**env,chart['coordinate'].node:x})
        R,J=a.at(chart['radius'].node),a.at(chart['jacobian'].node)
        zero(sy.diff(R,x)-J*R,('actual R_x equals J R',name))
        if name=='Rm_patch':zero(J-1/x,'actual Rm J=1/x')
        elif name=='O2_axial':zero(J-40*sy.exp(40*x),'actual O2 axial exponential J')
        else:zero(J-1,('actual constant outer J',name))
        evidence.append(dict(chart=name,radius=chart['radius'].node,jacobian=chart['jacobian'].node,nonconstant=name in ('Rm_patch','O2_axial')))
    return dict(passed=True,actual_R_native_derivative_equals_J_R=evidence,
        original_Rm_patch_and_O2_axial_nonconstant_J_retained=True,
        physical_Dy_repeated_instead_of_J_power_times_native_order=True)


def independent_variable_J_and_beta_calculus(owner):
    g=current.source.FunctionTransportGraph();x=g.symbol('independent_native');z=g.symbol('Z');phase=g.symbol('independent_lifted_phase')
    X,Z,Y=sy.symbols('independent_native Z independent_physical_y',real=True);N=sy.Symbol('fixed_positive_N',positive=True)
    n=g.symbol('fixed_positive_N');model=(1+Z**2+Z**5/19)*sy.sin(N*Y)+Z**3*sy.exp(2*Y)
    evidence=[]
    for name,yexpr,jexpr in (('Rm_patch',sy.log(X),1/X),('O2_axial',sy.exp(40*X),40*sy.exp(40*X))):
        variables={'independent_native':x,'Z':z};yref=current.source.expression(g,yexpr,variables)
        J=current.source.expression(g,jexpr,variables)
        D=current.RmOuterTotalNativeDerivative(g,x,phase,x,J,n,lambda row:(_ for _ in ()).throw(AssertionError('Fixture has no core')))
        D.shared['free'][phase.node]={'independent_native'}
        invJ=g.quotient(g.one,J,'positive fixture native Jacobian')
        jets=[]
        for k in range(6):
            slow=sy.diff((1+Z**2+Z**5/19),Z,k);exp_slow=sy.diff(Z**3,Z,k)
            jet=g.add(g.mul(current.source.expression(g,slow,variables),g.unary('sin',phase)),
                g.mul(current.source.expression(g,exp_slow,variables),g.unary('exp',g.mul(g.constant(2),yref))))
            jets.append(jet)
        env={x.node:X,z.node:Z,n.node:N,phase.node:N*yexpr}
        for order in range(1,5):
            jets=[g.mul(invJ,D.derivative(ref)) for ref in jets]
            a=annular_check.SourceAlgebra(g,env)
            for k,ref in enumerate(jets):
                expected=sy.diff(model,Y,order,Z,k).subs(Y,yexpr)
                zero(a.at(ref.node)-expected,('independent actual variable-J mixed derivative',name,order,k))
                evidence.append([name,order,k])
    # Exact interior derivative identity and exact lazy zero extension at
    # both support endpoints are checked separately on the actual beta rule.
    B=current.patch_leading.raw_beta(g,x);beta=sy.exp(-1/(1-X**2))
    D=current.RmOuterTotalNativeDerivative(g,x,phase,x,g.one,g.one,lambda row:None)
    D.shared['free'][B.node]={'independent_native'}
    class Interior(annular_check.SourceAlgebra):
        def at(self,index):
            if index not in self.cache and self.g.nodes[index]['operation']=='original_lazy_flat_branch':
                self.cache[index]=self.at(self.g.nodes[index]['active_body'])
            return super().at(index)
    beta_rows=[];current_ref=B
    for order in range(1,6):
        current_ref=D.derivative(current_ref)
        actual=Interior(g,{x.node:X,B.node:beta}).at(current_ref.node)
        assert sy.simplify(sy.together((actual-sy.diff(beta,X,order))/beta))==0,('original beta derivative',order)
        for edge in (-2,-1,1,2):
            zero(lexical.Algebra(g,{x.node:sy.Integer(edge)}).at(current_ref.node),('exact flat beta derivative',order,edge))
        beta_rows.append([order,current_ref.node])
    return dict(passed=True,independent_variable_J_y1_to4_Z0_to5_rows=evidence,
        exact_original_beta_interior_derivatives_and_four_flat_queries_per_order=beta_rows,
        repeated_physical_operator_includes_J_derivatives_and_NJ_phase=True,
        implicit_inverse_and_moving_integral_calculus_consumed_from=current.high_source.previous.previous.RECEIPT,
        beta_derivatives_not_numerical_or_zero_inside_support=True)


def actual_recovery_and_pressure_sources(owner):
    Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_physical_y',real=True)
    U=sy.Function('actual_arbitrary_U')(y,Z);M=sy.Function('actual_arbitrary_mean')(y,Z);evidence=[];operators=[]
    inherited=current.parent.linear.read(current.high_source.RECEIPT)['checks']['actual_recovery_and_pressure_sources']
    assert inherited['R_y_R_geometry_and_same_pinned_operator_transfer_N_y2_pressure_y3_C3']
    for name in current.CHARTS:
        chart=owner.charts[name];env,v=high_check.prior.active.prior.active_context(owner.rminus,False);delta=v['delta']
        for n,jet in enumerate(chart['physical']['U']):bind(env,jet,sy.diff(U,y,n),Z)
        bind(env,chart['physical']['m'][0],M,Z);a=Algebra(owner,env);means=[M]
        for n in range(1,5):means.append(sy.diff(U,y,n-1)-means[-1])
        for n in range(5):
            zero(a.at(chart['physical']['m'][n].value)-means[n],('actual Rm/outer mean ODE',name,n))
            Q=(2*Z*sy.diff(U,y,n)-(1-delta)*Z*means[n]-(1-Z**2)*sy.diff(means[n],Z))/(1-delta*Z**2)
            for k,ref in enumerate(rows(chart['physical']['Q'][n])):
                zero(a.at(ref)-sy.diff(Q,Z,k),('actual Rm/outer leading Q',name,n,k));evidence.append([name,n,k,ref.node])
        f={field+suffix:chart['physical'][field][n] for field in ('U','m','Q') for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        density=current.radial.pressure_density(owner.g,f,owner.parameters['delta'],chart['radius'])
        higher=current.high.higher_pressure_density(owner.g,f,owner.parameters['delta'],chart['radius'],density)
        assert encode(density)==encode(chart['density']) and encode(higher)==encode(chart['higher_density'])
        assert encode(high_check.reduced_S(owner,name))==owner.modified[name]['S']
        assert encode(chart['N_y'])==encode([density['N1p'],higher['N1p_y'],higher['N1p_yy']])
        operators.append([name,encode(chart['N_y'])])
    return dict(passed=True,actual_leading_Q_y0_to4_Z0_to4_rows=evidence,
        actual_mean_ODE_and_original_saved_S_replayed_on_six_charts=True,
        actual_N_y2_and_pressure_y3_fixed_operator_bindings=operators,
        arbitrary_physical_y_pressure_operator_theorem=current.high_source.RECEIPT,
        variable_native_J_not_inserted_into_physical_pressure_operator=True)


def actual_signed_history_sources(owner):
    g=owner.g
    initial={key:owner.jet(ids) for key,ids in owner.patch_input['initial'].items()};evidence=[]
    prior_name=None
    for name in current.CHARTS:
        old=owner.modified[name];loop=owner.loop[name];chart=owner.charts[name];x=chart['coordinate'];J=chart['jacobian'];R=chart['radius']
        incoming={field:owner.jet(ids) for field,ids in loop['incoming'].items()}
        densities={field:owner.jet(ids) for field,ids in loop['finite_N']['five_signed_density_increments'].items()}
        if name=='Rm_patch':
            assert loop['incoming']==owner.base.loop_input['endpoints']['postrestore']
            assert old['incoming_integrals']==owner.base.modified_input['endpoints']['postrestore']['integrals']
            history,parts,memory=current.patch_source.transport_log(g,densities,incoming,x,g.one,'Rm_patch')
        else:
            endpoint=owner.patch_input['endpoints']['Rh'] if prior_name=='Rm_patch' else owner.outer_input['endpoints'][prior_name]
            assert loop['incoming']==endpoint['correction'] and old['incoming_integrals']==endpoint['integrals']
            history,parts,memory=current.outer_source.transport_geometry(g,densities,incoming,x,chart['domain'][0],
                owner.ref(owner.backgrounds[name]['log_offset']),J,name,5)
        assert encode(history)==loop['history'] and encode(parts)==loop['contributions'] and encode(memory)==loop['memory']
        physical={field:owner.jet(ids) for field,ids in old['physical'].items()};S=owner.jet(old['S'])
        pressure_in={field:owner.jet(ids) for field,ids in old['incoming_integrals'].items()}
        if name=='Rm_patch':
            result=current.patch_source.pressure_continuation(g,S,physical,owner.ref(owner.leading_input['inlet']['radius']),x,pressure_in,initial)
        else:result=current.outer_source.pressure_geometry(g,S,physical,R,x,chart['domain'][0],J,pressure_in,initial,name)
        for packet,key in zip(result,('boundary','local_integrals','cumulative_integrals','pressure')):
            expected=old['finite_boundary'] if key=='boundary' and 'finite_boundary' in old else old[key]
            assert encode(packet)==expected,(name,key)
        for key,power in (('S0',1),('S1',2)):
            for n,node in enumerate(old['local_integrals'][key]):
                row=g.nodes[node]
                assert row['original_signed_reduced_pressure'] and row['physical_Jacobian_applied_once']
                assert row['physical_radius_power']==power and row['ordinary_slow_Z_derivative_order']==n
                assert row['lower']==loop['domain'][0] and row['upper']==x.node
                evidence.append([name,key,n,node])
        prior_name=name
    return dict(passed=True,actual_signed_partial_pressure_rows=evidence,
        actual_five_correction_transport_bodies_and_every_nonzero_predecessor_replayed=True,
        original_Rm_dx_over_x_and_outer_offset_J_measures_once=True,
        original_finite_boundaries_and_signed_m13_m15_formulas_replayed=True,
        source_coverage_does_not_claim_new_246_row_interface_closure=True)


def actual_own_sources_and_queries(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_own_physical_y',real=True)
    mass=sy.Function('nonzero_actual_cutoff_mass')(Z);env={}
    bind(env,owner.jet(owner.lower.rin.cut['endpoint']['m11'][:4]),mass,Z)
    for name in current.CHARTS:env[owner.charts[name]['radius'].node]=sy.Symbol('actual_R_'+name,positive=True)
    a=Algebra(owner,env);delta=a.at(owner.parameters['delta']);states=[];queries=[]
    for name in current.CHARTS:
        chart=owner.charts[name];f=chart['own']['fields'];R=a.at(chart['radius'])
        mean=mass/R*sy.exp(-y);radius=R*sy.exp(y)
        Q=(-(1+delta)*Z*mean-(1-Z**2)*sy.diff(mean,Z))/(1-delta*Z**2)
        model=dict(Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),M1=mass,mean1=mean,K1=mean,Q1=Q,V1=radius*Q)
        for field,body in model.items():
            for n,jet in enumerate(f[field]):
                zero(a.at(jet.value)-sy.diff(body,y,n).subs(y,0),('actual own Rm/outer physical flow',name,field,n));states.append([name,field,n])
        for field in ('M1','mean1','Q1','V1'):
            for ref,old in zip(rows(f[field][1]),chart['intermediate']['fields'][field+'_y']):
                zero(a.at(ref)-a.at(old),('actual own original first-y source',name,field))
        assert f['P1']==chart['own']['cumulative']['m13']
        assert encode(f['P1'][0])==chart['intermediate']['fields']['P1']
        assert encode(chart['own']['cumulative']['m15'][0])==chart['intermediate']['cumulative']['m15']
        for fraction in ('0','1/2','1'):
            packet=owner.at_fraction(name,fraction)
            assert not packet['numerical_oracle_available'] and len(packet['physical']['U'])==5
            assert len(packet['N_y'])==3 and len(packet['own']['fields']['P1'])==4
            queries.append([name,fraction,packet['native_query'].node])
    for name,fraction in (('not_a_chart','0'),('Rm_patch','-1/10'),('O2_axial','11/10')):
        try:owner.at_fraction(name,fraction)
        except ValueError:pass
        else:raise AssertionError('Outside accepted source chart')
    return dict(passed=True,actual_own_state_y0_to4_scalar_rows=states,actual_full_native_chart_fraction_queries=queries,
        original_own_first_y_and_one_plus_delta_recovery_retained=True,
        actual_fixed_pressure_and_all_five_cumulative_value_sources_retained=True,
        three_invalid_domain_queries_rejected=True,numerical_oracles_and_interface_joins_not_claimed=True)


CHECKS=(source_contract,original_geometry,independent_variable_J_and_beta_calculus,
        actual_recovery_and_pressure_sources,actual_signed_history_sources,actual_own_sources_and_queries)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('RM_OUTER_HIGH_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['Rm_outer_high_radial_source_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),actual_full_original_charts=list(current.CHARTS),
        actual_native_to_physical_derivative_requests=108,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,source_interfaces_closed=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_SIX_RM_OUTER_HIGH_RADIAL_SOURCE_CHARTS_ACCEPTED',flush=True)
    return rec
