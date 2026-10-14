"""Genuine analytic lower radial derivatives, pressure and full Rin traces."""
import copy
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_lower_higher_radial_Rin_join as current
import lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join_check as previous_check

rows,encode=current.outer.rows,current.outer.encode
zero,bind=previous_check.zero,previous_check.bind


class HigherAlgebra(previous_check.RinAlgebra):
    def __init__(self,owner,env,evidence=None):
        super().__init__(owner.rin,env,evidence);self.higher=owner
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation'];owner=self.higher
        if op=='current_original_unique_analytic_core_component':
            owner.guard_core(row);r=self.at(row['rho']);Z=sy.Symbol('Z',real=True)
            x=sy.Symbol('analytic_core_argument',real=True)
            expression=sy.Function('actual_core_'+row['quantity'])(x,Z)
            result=sy.diff(expression,x,row['radial_order'],Z,row['ordinary_slow_Z_derivative_order']).subs(x,r)
        elif op=='analytic_unary' and row['name']=='original_flat_sigma_order3':
            assert row['active_interval']==[0,1] and row['flat_value']==0
            assert row['active_body_not_evaluated_on_flat_branch'] and row['ordinary_from_Taylor_factor']==6
            assert row['derivative_source']==owner.sigma_interior_pin
            x=sy.cancel(self.at(row['argument']))
            if x.is_nonpositive or (x-1).is_nonnegative:result=sy.Integer(0)
            else:
                t=sy.Symbol('actual_sigma_argument',real=True)
                result=sy.diff(sy.Function('actual_original_sigma')(t),t,3).subs(t,x)
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            symbol=sy.Symbol(self.g.nodes[row['variable']]['name'],real=True);value=self.at(row['value'])
            env={key:body.subs(symbol,value) if hasattr(body,'subs') else body for key,body in self.bindings.items()}
            env[row['variable']]=value;result=HigherAlgebra(owner,env,self.evidence).at(row['expression'])
        else:return super().at(i)
        self.cache[i]=result;return result


def leading_environment(owner,transport=False):
    Z=sy.Symbol('Z',real=True);s=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True)
    r=sy.Symbol('independent_physical_rho',positive=True)
    H,SC,L=sy.symbols('original_hb original_sc original_logLambda',positive=True)
    p=owner.parameters;env={p['hb'].node:H,p['sc'].node:SC,p['logLambda'].node:L}
    ratio_s=sy.Function('retained_ratio')(s,Z);psi_s=sy.Function('retained_collar_Psi')(s,Z)
    mean_s=sy.Function('retained_collar_mean')(s,Z)
    if transport:
        bind(env,owner.ratio[0],ratio_s,Z)
        bind(env,owner.radial['Psi'][0],psi_s,Z)
        bind(env,owner.radial['mean_Psi'][0],mean_s,Z)
    phase=sy.log(r/4)/H
    ratio=sy.Function('retained_ratio_by_rho')(r,Z)
    psi=sy.Function('retained_Psi_by_rho')(r,Z);mean=sy.Function('retained_mean_by_rho')(r,Z)
    def normalize(expression):
        value=expression.subs(s,phase)
        return value.xreplace({ratio_s.subs(s,phase):ratio,psi_s.subs(s,phase):psi,mean_s.subs(s,phase):mean}).doit()
    return env,dict(Z=Z,s=s,r=r,H=H,SC=SC,L=L,phase=phase,ratio=ratio,psi=psi,mean=mean,normalize=normalize)


def source_contract(owner):
    owner.assert_graph();raw=current.parent.linear.read(current.previous.NAME)
    rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and all(rec[key] for key in current.previous.GATES)
    assert not any(rec[key] for key in current.previous.OPEN)
    assert rec['actual_Rin_contracted_join_rows']==138
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    for name,digest in rec['input_hashes'].items():assert owner.hashes[name]==digest==current.sha(name)
    assert owner.contract==owner.rin.rc.contract==current.CONTRACT
    for key,data in owner.old['radial_fields'].items():
        name=key.removesuffix('_rho_rows')
        assert encode(owner.radial[name][:3])==data
    for key,data in owner.old['core_radial_rows'].items():assert encode(owner.core_radial[key][:3])==data
    assert encode(owner.N[0])==owner.old['N1p']
    for key,jets in owner.radial.items():
        assert len(jets)==5 and all(len(rows(jet))==(5 if key=='Q' else 6) for jet in jets)
    assert len(owner.N)==3 and all(len(rows(jet))==4 for jet in owner.N)
    count=0
    for side in ('left','right'):
        for spec in owner.contract:
            jets=owner.join[side][spec['group']][spec['field']]
            assert len(jets)==spec['physical_y_order']+1
            for jet in jets:
                assert len(rows(jet))==spec['ordinary_Z_order']+1
                for ref in rows(jet):
                    row=owner.g.nodes[ref.node]
                    assert row['operation']=='function_substitution' and row['Z_independent_substitution']
                    assert row['variable']==owner.coordinate.node and row['value']==owner.rin.endpoint.node
                    count+=1
    assert count==492
    for side,old in (('left',owner.rin.left),('right',owner.rin.right)):
        data=getattr(owner,side)
        for group,packet in old.items():
            for key,jets in packet.items():assert encode(data[group][key][:len(jets)])==encode(jets)
    return dict(passed=True,actual_old_radial0_1_2_and_N0_preserved=True,
        actual_Rin138_prefix_retained=True,fixed_Z_independent_endpoint_rows=count,
        quantitative_radial4_Z5_bounds_and_numerical_oracles_not_claimed=True)


def actual_core_and_sigma_sources(owner):
    g=owner.g;evidence=[]
    for key,jets in owner.core_radial.items():
        for n,jet in enumerate(jets):
            for k,ref in enumerate(rows(jet)):
                row=g.nodes[ref.node];owner.guard_core(row)
                assert row['radial_order']==n and row['ordinary_slow_Z_derivative_order']==k
                assert row['actual_twenty_term_operator'] and row['unique_analytic_Banach_solution']
                assert row['admitted_solution_not_a_model_or_enclosure_center']
                for pin_key in ('radial_recurrence_provider','actual_source_provider','actual_green_inverse_provider'):
                    pin=row[pin_key];assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
                evidence.append([key,n,k])
    bad=copy.deepcopy(owner.core_templates['Phi']);bad['rho']=-1
    try:owner.guard_core(bad)
    except ValueError:pass
    else:raise AssertionError('Unadmitted analytic component accepted')
    env,v=leading_environment(owner);x=v['s'];S=sy.Function('sigma_for_independent_derivative')(x)
    D=2/x**3+2/(1-x)**3;first=S*(1-S)*D;expected=first
    for _ in range(2):expected=sy.diff(expected,x).subs(sy.diff(S,x),first)
    scalar=sy.Symbol('actual_sigma_value',real=True)
    env[owner.sigma[0].node]=scalar
    body=g.nodes[owner.sigma[3].node]['interior_body']
    zero(HigherAlgebra(owner,env).at(body)-expected.subs(S,scalar),'actual sigma3 interior source')
    env,v=leading_environment(owner);alg=HigherAlgebra(owner,env);r,H=v['r'],v['H'];norm=v['normalize']
    sigma=sy.Function('actual_original_sigma')(v['phase']);chi=1+(H-1)*sigma
    for n,jet in enumerate(owner.chi):zero(norm(alg.at(jet.value))-sy.diff(chi,r,n),('actual chi rho derivative',n))
    return dict(passed=True,actual_analytic_core_mixed_rows=evidence,changed_core_coordinate_rejected=True,
        original_sigma3_ordinary_factor_six_and_lazy_flat_source=True,
        actual_chi_rho0_rho1_rho2_rho3_chain_rules=True,
        analytic_function_derivatives_not_numerical_enclosure_centers=True)


def actual_original_transport_FTC(owner):
    env,v=leading_environment(owner,True);alg=HigherAlgebra(owner,env);g=owner.g
    f=owner.rin.lower['fields'];s=owner.coordinate
    def upper_body(node):
        row=g.nodes[node]
        assert row['operation']=='definite_integral' and row['lower']==0 and row['upper']==s.node
        assert row['endpoint_independent_of_Z'] and row['ordinary_slow_Z_derivative_order']==0
        var=g.symbol(row['variable'])
        return g.node('function_substitution',expression=row['integrand'],variable=var.node,value=s.node,Z_independent_substitution=True)
    jac=alg.at(owner.parameters['hb'])*alg.at(owner.rho)
    logbody=upper_body(f['logratio'][0])
    zero(alg.at(logbody)/jac-alg.at(owner.g_rows[0].value),'actual original logratio FTC')
    candidates=[i for i,row in enumerate(g.nodes[:len(owner.rin.prefix)]) if row['operation']=='definite_integral'
        and row.get('actual_Psi_radial_transport') and row.get('ordinary_slow_Z_derivative_order')==0
        and not row.get('ordinary_Z_derivative_order',0)]
    assert len(candidates)==1;psibody=upper_body(candidates[0])
    assert candidates[0] in g.nodes[f['Psi_rho_rows'][0][0]]['arguments']
    zero(alg.at(psibody)/jac-alg.at(owner.reconstructed['Psi'][1].value),'actual original Psi FTC')
    mean_candidates=[i for i,row in enumerate(g.nodes[:len(owner.rin.prefix)]) if row['operation']=='definite_integral'
        and row.get('actual_own_mean_FTC') and row.get('ordinary_slow_Z_derivative_order')==0
        and not row.get('ordinary_Z_derivative_order',0)]
    assert len(mean_candidates)==1
    mean_numerator=g.nodes[f['mean_Psi_rho_rows'][0][0]]['numerator']
    assert mean_candidates[0] in g.nodes[mean_numerator]['arguments']
    meanbody=upper_body(mean_candidates[0])
    zero(alg.at(meanbody)/jac-alg.at(owner.radial['Psi'][0].value),'actual nonzero mean history FTC')
    assert g.nodes[f['ratio'][0]]==dict(operation='analytic_unary',name='exp',argument=f['logratio'][0])
    return dict(passed=True,actual_logratio_Psi_and_mean_integral_FTC=True,
        defining_nonzero_rho4_inlets_and_full_histories_preserved=True,
        ratio_derivatives_follow_original_ODE_not_free_differentiation=True)


def actual_collar_radial4(owner):
    env,v=leading_environment(owner,True);alg=HigherAlgebra(owner,env)
    r,Z,H=v['r'],v['Z'],v['H'];norm=v['normalize'];ratio=v['ratio']
    chi=1+(H-1)*sy.Function('actual_original_sigma')(v['phase'])
    phi=sy.Function('actual_core_Phi')(r,Z);psi=sy.Function('actual_core_Psi')(r,Z)
    gg=(chi-1)*sy.diff(phi,r)/phi
    ratio_values=[ratio,ratio*gg,ratio*(gg**2+sy.diff(gg,r)),ratio*(gg**3+3*gg*sy.diff(gg,r)+sy.diff(gg,r,2))]
    ratio_replace={sy.diff(ratio,r,k):ratio_values[k] for k in (1,2,3)}
    for n,jet in enumerate(owner.ratio):zero(norm(alg.at(jet.value))-ratio_values[n],('actual ratio derivative',n))
    expected={}
    for key,core in (('Phi',phi),('Psi',psi)):
        expected[key]=[norm(alg.at(owner.radial[key][0].value))]
        for n in range(1,5):
            value=sy.diff(ratio*chi*sy.diff(core,r),r,n-1).xreplace(ratio_replace)
            expected[key].append(value)
            zero(norm(alg.at(owner.reconstructed[key][n].value))-value,('independent transport derivative',key,n))
            zero(norm(alg.at(owner.radial[key][n].value))-value,('retained/new actual radial prefix',key,n))
    M,Psi=v['mean'],v['psi'];mean_values=[M,(Psi-M)/r]
    replacements={sy.diff(M,r):mean_values[1]}
    for n in range(1,4):replacements[sy.diff(Psi,r,n)]=expected['Psi'][n]
    for n in range(2,5):
        value=sy.diff(mean_values[-1],r).xreplace(replacements).xreplace(ratio_replace)
        mean_values.append(value)
    for n,jet in enumerate(owner.radial['mean_Psi']):
        zero(norm(alg.at(jet.value))-mean_values[n],('independent mean FTC derivative',n))
    return dict(passed=True,actual_old_lower_radial0_1_2_reproduced=True,
        actual_Phi_Psi_radial3_radial4_from_independent_ODE_product_differentiation=True,
        actual_mean_radial0_to_4_from_nonzero_history_FTC=True,
        ordinary_Z0_to_5_follow_same_admitted_analytic_and_product_algebra=True)


def actual_Q_and_pressure_derivatives(owner):
    env,v=leading_environment(owner);Z,r=v['Z'],v['r'];s=v['s'];H=v['H'];norm=v['normalize']
    U=sy.Function('arbitrary_leading_U')(r,Z);M=sy.Function('arbitrary_leading_mean')(r,Z)
    de=sy.Symbol('fixed_delta',real=True);env[owner.parameters['delta'].node]=de
    for key,body in (('U',U),('m',M)):
        for n,jet in enumerate(owner.radial[key]):bind(env,jet,sy.diff(body,r,n),Z)
    alg=HigherAlgebra(owner,env);L=1-de*Z**2;d=1-Z**2;h=1-de
    q=(2*Z*U-h*Z*M-d*sy.diff(M,Z))/L
    for n,jet in enumerate(owner.radial['Q']):
        zero(alg.at(jet.value)-sy.diff(q,r,n),('actual leading Q radial recovery',n))
    # Use arbitrary correlated Q/U derivative inputs for an independent
    # differentiation of the physical pressure operator, including1/R.
    env,v=leading_environment(owner);env[owner.parameters['delta'].node]=de
    q=sy.Function('arbitrary_actual_Q')(r,Z);u=sy.Function('arbitrary_actual_U')(r,Z)
    for key,body in (('Q',q),('U',u)):
        for n,jet in enumerate(owner.radial[key]):bind(env,jet,sy.diff(body,r,n),Z)
    eps=sy.exp(-v['L'])
    # The same exact coordinate map also covers factors flattened into
    # products by the original graph, not just a top-level radius node.
    env[owner.coordinate.node]=v['phase'];env[owner.rho.node]=r
    env[owner.radius.node]=eps*r;alg=HigherAlgebra(owner,env)
    qy=r*sy.diff(q,r)
    S=(q+h*Z*sy.diff(q,Z)/2+qy)/L+q*(q/2+qy)+u*(-2*Z*q+d*sy.diff(q,Z)-2*Z*qy)/L
    N=-S/2+sy.exp(v['L'])*(2*sy.diff(q,r)+r*sy.diff(q,r,2))
    for n in range(3):
        zero(alg.at(owner.S[n].value)-sy.diff(S,r,n),('actual S rho derivative',n))
        zero(alg.at(owner.N_reconstructed[n].value)-sy.diff(N,r,n),('actual N rho derivative',n))
        zero(alg.at(owner.N[n].value)-sy.diff(N,r,n),('preserved original N and new derivatives',n))
    # Physical radius must agree with epsilon*rho in this independent frame.
    env[owner.radius.node]=eps*r;alg=HigherAlgebra(owner,env)
    Py=eps*r*N;m15y=-(eps*r)**2*N
    by13=eps*v['H']*r*N;by15=-eps**2*v['H']*r**2*N
    for n in range(3):
        zero(alg.at(owner.pressure_y[n].value)-Py,('actual P primitive y derivative',n+1))
        zero(alg.at(owner.m15_y[n].value)-m15y,('actual m15 primitive y derivative',n+1))
        zero(alg.at(owner.signed_bodies['m13'][n].value)-sy.diff(by13,r,n),('signed physical m13 rho derivative',n))
        zero(alg.at(owner.signed_bodies['m15'][n].value)-sy.diff(by15,r,n),('signed physical m15 rho derivative',n))
        Py=r*sy.diff(Py,r);m15y=r*sy.diff(m15y,r)
    return dict(passed=True,actual_Q_radial0_to_4_leading_one_minus_delta_recovery=True,
        actual_N_rho0_rho1_rho2_independent_pressure_operator_differentiation=True,
        physical_y_equals_rho_d_rho_and_inverse_R_derivatives_retained=True,
        actual_pressure_and_m15_primitive_y1_y2_y3_and_signed_rho0_rho1_rho2=True,
        ordinary_Z0_to_3_by_actual_C5_C4_C3_jet_composition=True)


def actual_flat_high_own_and_Rin_join(owner):
    env,Z,s=previous_check.environment(owner.rin,history=True,pressure_theorem=True)
    for n in (1,2):bind(env,owner.N[n],sy.Function('actual_common_N_rho'+str(n))(s,Z),Z)
    alg=HigherAlgebra(owner,env);a=current.outer.Algebra(owner.g,3)
    # The high state extension is valid on the actual whole flat collar,
    # where the original cutoff argument is1+t, t>=0 (Rin is t=1).
    t=sy.Symbol('flat_collar_fraction',nonnegative=True);flat_env=dict(env)
    q=alg.at(owner.rin.endpoint);flat_env[owner.coordinate.node]=q*(3+t)/4
    flat=HigherAlgebra(owner,flat_env)
    for key in ('F1','Uz1'):
        for jet in owner.rin.left['fields'][key]:zero(flat.at(jet.value),('whole flat cutoff source',key))
    body=owner.rin.cut['history_bodies']['first']['mean_integrand']
    candidates=[row for row in owner.g.nodes[:len(owner.rin.prefix)] if row['operation']=='definite_integral'
        and row.get('actual_cutoff_own_mean') and row['integrand']==body[0]]
    assert len(candidates)==1;row=candidates[0]
    upper=owner.g.node('function_substitution',expression=row['integrand'],variable=owner.g.symbol(row['variable']).node,
        value=owner.coordinate.node,Z_independent_substitution=True)
    zero(flat.at(upper),'actual mean-history upper density on entire flat collar')
    # Derive all high own rows from the whole-source relations. Endpoint
    # agreement alone cannot establish that an assigned high row is a
    # derivative of the original own field.
    fields=owner.left['fields'];radius=flat.at(owner.radius)
    mean=flat.at(fields['mean1'][0].value);Uz=flat.at(fields['Uz1'][0].value)
    de=flat.at(owner.parameters['delta']);L=1-de*Z**2;d=1-Z**2
    own_Q=(2*Z*Uz-(1+de)*Z*mean-d*flat.at(rows(fields['mean1'][0])[1]))/L
    zero(flat.at(fields['K1'][0].value)-mean+Uz,'actual own K definition on flat collar')
    zero(flat.at(fields['Q1'][0].value)-own_Q,'actual own one-plus-delta Q recovery on flat collar')
    zero(flat.at(fields['V1'][0].value)-radius*own_Q,'actual own V=R Q definition')
    zero(flat.at(fields['M1'][0].value)-radius*mean,'actual own M=R mean definition')
    # This independent y-flow is the unique solution of the actual mean
    # FTC ODE mean_y=-mean with its retained, arbitrary nonzero datum.
    y=sy.Symbol('independent_flat_log_radius',real=True)
    mean_flow=mean*sy.exp(-y);R_flow=radius*sy.exp(y)
    Q_flow=(-(1+de)*Z*mean_flow-d*sy.diff(mean_flow,Z))/L
    model=dict(Uz1=sy.Integer(0),F1=sy.Integer(0),Utheta1=sy.Integer(0),
        mean1=mean_flow,K1=mean_flow,Q1=Q_flow,V1=R_flow*Q_flow,M1=R_flow*mean_flow)
    flat_state_rows=[]
    for key,flow in model.items():
        for j,jet in enumerate(fields[key]):
            zero(flat.at(jet.value)-sy.diff(flow,y,j).subs(y,0),('actual own flat-source y derivative',key,j))
            flat_state_rows.append([key,j])
    zero_density_rows=[];zero_cumulative_rows=[];moment_recipes=[]
    for key in ('m11','m12','m14'):
        for n,node in enumerate(owner.rin.cut['integrands']['first'][key]):
            zero(flat.at(node),('actual whole-flat signed density',key,n))
            zero_density_rows.append([key,n,node])
        for j,jet in enumerate(owner.left['cumulative'][key][1:],1):
            zero(flat.at(jet.value),('actual flat cumulative FTC derivative',key,j))
            zero_cumulative_rows.append([key,j])
        if key in ('m12','m14'):
            # Anchor the constant-tail inference to the actual original
            # cumulative integrals and their signed density inputs.
            body_env=dict(env);body_function=sy.Function('actual_signed_'+key)(s,Z)
            bind(body_env,owner.rin.jet(owner.rin.cut['integrands']['first'][key]),body_function,Z)
            body_alg=HigherAlgebra(owner,body_env)
            for n,node in enumerate(owner.rin.cut['cumulative']['first'][key]):
                matches=[i for i,rec in enumerate(owner.g.nodes[:len(owner.rin.prefix)]) if rec['operation']=='definite_integral'
                    and rec.get('actual_signed_n1_inner_moment') and rec.get('paper_moment')==key
                    and rec.get('full_core_to_first_history_retained') and rec.get('ordinary_slow_Z_derivative_order')==n
                    and not rec.get('ordinary_Z_derivative_order',0)]
                assert len(matches)==1;integral=owner.g.nodes[matches[0]]
                assert matches[0] in owner.g.nodes[node]['arguments']
                assert integral['lower']==0 and integral['upper']==owner.coordinate.node
                assert integral['endpoint_independent_of_Z']
                variable=sy.Symbol(integral['variable'],real=True);hb=body_alg.at(owner.parameters['hb'])
                expected=hb*4*sy.exp(hb*variable)*sy.diff(body_function,Z,n).subs(s,variable)
                zero(body_alg.at(integral['integrand'])-expected,('actual signed cumulative FTC body',key,n))
                moment_recipes.append([key,n,matches[0],integral['integrand']])
    for j in range(1,4):
        zero(flat.at(fields['P1'][j].value)-flat.at(owner.pressure_y[j-1].value),('actual nonzero pressure flat source',j))
        zero(flat.at(owner.left['cumulative']['m13'][j].value)-flat.at(owner.pressure_y[j-1].value),('actual m13 flat source',j))
        zero(flat.at(owner.left['cumulative']['m15'][j].value)-flat.at(owner.m15_y[j-1].value),('actual signed m15 flat source',j))
    evidence=[]
    for spec in owner.contract:
        for j in range(spec['physical_y_order']+1):
            left=rows(owner.join['left'][spec['group']][spec['field']][j])
            right=rows(owner.join['right'][spec['group']][spec['field']][j])
            for n,(l,r) in enumerate(zip(left,right)):
                zero(alg.at(l)-alg.at(r),(spec['group'],spec['field'],j,n))
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=j,
                    ordinary_Z_order=n,left_node=l.node,right_node=r.node,
                    proof='actual shared high pressure source and flat-support/FTC derivative'))
    assert len(evidence)==246
    return dict(passed=True,actual_Rin_contracted_mixed_join_rows=evidence,
        whole_original_flat_collar_and_own_mean_ODE_supply_y3_y4=True,
        actual_own_K_Q_one_plus_delta_V_and_M_source_relations=True,
        independent_flat_own_state_y0_to_y4_rows=flat_state_rows,
        actual_zero_signed_density_C3_rows=zero_density_rows,
        actual_constant_m11_m12_m14_FTC_rows=zero_cumulative_rows,
        actual_m12_m14_retained_cumulative_integral_recipes=moment_recipes,
        actual_nonzero_pressure_m13_m15_y1_y2_y3_source_relations=True,
        own_one_plus_delta_and_nonzero_pressure_and_moment_histories_preserved=True,
        actual_m15_is_primitive_its_first_derivative_is_negative_R_squared_N=True,
        global_n1_other_joins_numerical_bounds_recursion_corrected_NS_open=True)


CHECKS=(source_contract,actual_core_and_sigma_sources,actual_original_transport_FTC,
    actual_collar_radial4,actual_Q_and_pressure_derivatives,actual_flat_high_own_and_Rin_join)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('LOWER_HIGH_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['lower_higher_radial_Rin_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_Rin_contracted_mixed_join_rows=246,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_LOWER_RADIAL4_PRESSURE_RHO2_RIN246_ACCEPTED',flush=True)
    return rec
