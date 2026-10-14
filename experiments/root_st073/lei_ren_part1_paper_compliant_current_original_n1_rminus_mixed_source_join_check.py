"""Original common IVP, active high radial source, own pressure FTC/rminus."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_rminus_mixed_source_join as current
import lei_ren_part1_paper_compliant_current_original_lower_higher_radial_Rin_join_check as lower_check

rin_check=lower_check.previous_check
rows,encode=current.outer.rows,current.outer.encode
zero,bind=lower_check.zero,lower_check.bind


class RminusAlgebra(lower_check.HigherAlgebra):
    def __init__(self,owner,env,evidence=None,scope=None):
        # Keep source inputs unevaluated across nested coordinate bindings.
        # A fixed cutoff endpoint inside a later collar query must shadow the
        # outer lower coordinate, rather than inherit its already-substituted
        # history. Both sides retain the very same nonzero primitive.
        self.source_bindings=env;self.scope={} if scope is None else scope
        substitutions={sy.Symbol(owner.g.nodes[key]['name'],real=True):value
            for key,value in self.scope.items()}
        context={key:body.subs(substitutions,simultaneous=True) if hasattr(body,'subs') else body
            for key,body in env.items()}
        context.update(self.scope)
        super().__init__(owner.lower,context,evidence);self.rminus=owner
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation'];owner=self.rminus
        if op=='original_lazy_flat_branch' and row.get('original_collar_coordinate')==owner.coordinate.node:
            assert row['original_collar_upper']==owner.loop['loop']['original_collar_upper']
            assert row['active_body_not_evaluated_on_flat_branch'] and row['flat_value']==0
            predicate=self.g.nodes[row['flat_predicate']];assert predicate['operation']=='logical_or'
            assert any(self.g.nodes[n]==dict(operation='exact_real_comparison',operator='<=',
                left=owner.coordinate.node,right=row['original_collar_upper']) for n in predicate['arguments'])
            self.evidence['flat'].add(i);result=sy.Integer(0)
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            value=self.at(row['value']);scope=dict(self.scope);scope[row['variable']]=value
            result=RminusAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['expression'])
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z']
            lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:self.evidence['empty'].add(i);result=sy.Integer(0)
            else:
                # The integral variable is bound inside its body even when
                # its name is also the collar's moving query coordinate.
                variable=sy.Symbol(row['variable'],real=True)
                scope=dict(self.scope);scope[self.g.symbol(row['variable']).node]=variable
                body=RminusAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['integrand'])
                if row.get('actual_finite_N_signed_correction'):
                    assert row['lower']==owner.right_endpoint.node
                    zero(body,'actual first-micro local finite-N density on flat collar');result=sy.Integer(0)
                else:
                    terms=[]
                    for term in sy.Add.make_args(sy.expand_mul(sy.expand_power_exp(body))):
                        coefficient,integrand=term.as_independent(variable,as_Add=False)
                        terms.append(coefficient*sy.Integral(integrand,(variable,lo,hi)))
                    result=sy.Add(*terms)
        elif op=='analytic_unary' and row['name']=='expm1':result=sy.exp(self.at(row['argument']))-1
        else:return super().at(i)
        self.cache[i]=result;return result


def active_context(owner,common=True):
    env,v=lower_check.leading_environment(owner.lower,False)
    s=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True);Z=v['Z'];H=v['H'];r=v['r']
    p=owner.parameters;eps=sy.exp(-v['L']);P=sy.Symbol('original_Pstar',positive=True)
    de,j=sy.symbols('fixed_delta fixed_j',real=True)
    env.update({owner.active_sc.node:v['SC'],p['delta'].node:de,p['j'].node:j,
        owner.common['inlet_binding']['Pstar']:P})
    bind(env,owner.jet(owner.bridge['core7']['U0']),4*Z+j,Z)
    bind(env,owner.jet(owner.bridge['core7']['F0']),sy.Function('actual_axis_F0')(Z),Z)
    ratio=sy.Function('retained_ratio')(s,Z);psi=sy.Function('retained_collar_Psi')(s,Z)
    mean=sy.Function('retained_collar_mean')(s,Z);rho=4*sy.exp(H*s)
    corePhi=sy.Function('actual_core_Phi');corePsi=sy.Function('actual_core_Psi')
    u=sy.Symbol('independent_core_rho',positive=True)
    if common:
        # These are the registered whole-source IVP identities consumed in
        # source_contract; they identify actual source inputs, not own outputs.
        bind(env,owner.ratio[0],ratio,Z)
        bind(env,owner.jet(owner.micro['actual']['phi'][:6]),ratio*corePhi(rho,Z),Z)
        bind(env,owner.jet(owner.bridge['comparison_micro']['phi'][:6]),corePhi(rho,Z),Z)
        bind(env,owner.jet(owner.micro['actual']['U'][:6]),4*Z+j+eps*psi,Z)
        bind(env,owner.jet(owner.micro['moments']['M'][:6]),4*Z+j+eps*mean,Z)
        phi_r=sy.diff(corePhi(u,Z),u).subs(u,rho);psi_r=sy.diff(corePsi(u,Z),u).subs(u,rho)
        bind(env,owner.jet(owner.bridge['comparison_micro']['directions']['D_over_R'][:6]),-2*phi_r/(eps*corePhi(rho,Z)),Z)
        bind(env,owner.jet(owner.micro['actual']['full_signed_axial_drive'][:6]),-eps*rho*psi_r,Z)
    phase=sy.log(r/4)/H
    def normalize(expression):
        return expression.subs(s,phase).xreplace({ratio.subs(s,phase):v['ratio'],psi.subs(s,phase):v['psi'],mean.subs(s,phase):v['mean']}).doit()
    return env,dict(v,s=s,Z=Z,eps=eps,P=P,delta=de,j=j,rho=rho,normalize=normalize)


def source_contract(owner):
    owner.assert_graph();g=owner.g;records={}
    for module in (current.previous,current.common_source,current.common_source.previous,
        current.bridge_source,current.loop_source,current.modified_source):
        rec=current.parent.linear.read(module.RECEIPT);raw=current.parent.linear.read(module.NAME)
        assert rec['all_passed'] and rec['source_family']==owner.family==raw['source_family']
        assert all(rec[key] for key in module.GATES) and not any(rec[key] for key in module.OPEN)
        assert raw['exact_graph_nodes']==g.nodes[:len(raw['exact_graph_nodes'])]
        for name,digest in rec['input_hashes'].items():assert owner.hashes[name]==digest==current.sha(name)
        records[module.NAME]=rec
    common=records[current.common_source.NAME]['checks'];ivp=common['common_IVP']
    for key in ('actual_lower_and_prescribed_axial_integral_function_identity',
        'common_true_rho4_inlet_and_comparison_alpha_one_from_actual_sc',
        'angular_positive_linear_IVP_uniqueness_with_nonconstant_chi',
        'all_six_same_inlet_Duhamel_histories_then_have_identical_solutions',
        'mean_and_absolute_pressure_follow_same_M_C_FTC_IVP'):assert ivp[key]
    assert common['flat_collar']['common_seam_sc_half_lies_strictly_inside_original_two_sided_flat_collar']
    assert common['modified_decomposition']['actual_common_IVP_background_fields_feed_all_eight_modified_outputs']
    intermediate_module=current.common_source.previous
    intermediate_receipt=records[intermediate_module.NAME]['checks']
    for key in ('positive_order_velocities_zero_without_zeroing_leading_velocity',
        'actual_preRc_coordinate_radius_domain_and_C3_pressure_source_identity'):
        assert intermediate_receipt['source'][key]
    for key in ('arbitrary_nonzero_incoming_mass_and_pressure_functions',
        'positive_order_one_plus_delta_and_true_mean_Z3',
        'own_mean_y_and_Q_y_are_minus_the_own_mean_and_Q',
        'own_V_y_and_mass_y_zero_only_outside_velocity_supports'):
        assert intermediate_receipt['actual_regional_recovery'][key]
    intermediate_raw=current.parent.linear.read(intermediate_module.NAME)['intermediate_own_C3_chart_bindings']
    chart=owner.intermediate
    assert chart==intermediate_raw['charts']['first_micro']
    assert 'first_micro' in intermediate_raw['preRc_names'] and chart['region']=='beforeRc'
    assert chart['complete_incoming_mass_not_reset']
    assert chart['fields']['M1']==chart['cumulative']['m11']
    assert chart['fields']['P1']==chart['cumulative']['m13']
    for key in ('Uz1','Utheta1','F1'):assert all(node==g.zero.node for node in chart['fields'][key])
    for key in ('m11','m12','m14'):
        assert chart['cumulative'][key]==owner.lower.rin.cut['endpoint'][key][:4]
        # These exact original cutoff endpoint functions precede the new
        # active coordinate in the checked append-only source DAG.
        assert all(node<owner.coordinate.node for node in chart['cumulative'][key])
    # The pinned intermediate receipt proves the original recovery for these
    # exact archived nodes. Retain its actual first-y rows as well; the
    # independent whole-flow check below derives the additional orders.
    assert chart['fields']['K1']==chart['fields']['mean1']
    for key in ('M1','V1'):
        assert encode(owner.right['fields'][key][1])==chart['fields'][key+'_y']
    sc=owner.common['sc_binding']
    assert sc['resolved_exact_sc_node']==owner.parameters['sc'].node and sc['active_sc_node']==owner.active_sc.node
    assert sc['exact_singleton_recomputed_not_interval_overlap'] and sc['rho4point1_common_analytic_continuation_margin_positive']
    assert owner.coordinate.node==owner.micro['phase']==owner.loop['coordinate']==owner.intermediate['coordinate']
    assert owner.radius.node==owner.micro['radius']==owner.intermediate['radius']
    assert owner.jacobian.node==owner.parameters['hb'].node==owner.micro['physical_measure']
    assert owner.right_endpoint.node==owner.loop['domain'][0]
    assert owner.left_endpoint.node==owner.parameters['sc_half'].node
    assert owner.loop['roots']['E']==owner.micro['fields']['E'][:6]
    assert owner.loop['roots']['V']==owner.micro['fields']['V'][:6]
    for key in ('U','m','Q'):
        assert encode(owner.physical[key][0])==owner.modified['physical'][key]
        assert encode(owner.physical[key][1])==owner.modified['physical'][key+'_y']
    count=0
    for side,variable,value in (('left',owner.lower.coordinate.node,owner.left_endpoint.node),('right',owner.coordinate.node,owner.right_endpoint.node)):
        for spec in owner.contract:
            jets=owner.join[side][spec['group']][spec['field']]
            assert len(jets)==spec['physical_y_order']+1
            for jet in jets:
                assert len(rows(jet))==spec['ordinary_Z_order']+1
                for ref in rows(jet):
                    row=g.nodes[ref.node];assert row['operation']=='function_substitution'
                    assert row['variable']==variable and row['value']==value and row['Z_independent_substitution']
                    count+=1
    assert count==492
    return dict(passed=True,actual_common_whole_function_IVP_receipt=current.common_source.RECEIPT,
        actual_intermediate_outside_support_receipt=intermediate_module.RECEIPT,
        actual_first_micro_recovery_equals_original_zero_velocity_fields=True,
        actual_constant_m11_m12_m14_are_original_nonzero_cutoff_endpoints=True,
        actual_original_source_selector_binding=sc,actual_fixed_endpoint_substitution_rows=count,
        rminus_uses_sc_half_not_Rin_sc_quarter=True,actual_physical_Jacobian_hb_applied_once=True)


def actual_flat_first_micro_and_native_prefix(owner):
    env,v=active_context(owner);alg=RminusAlgebra(owner,env);Z=v['Z'];g=owner.g
    checked=[]
    for key,ids in [(key,owner.loop['loop'][key]) for key in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi')]+[
        ('delta_E',owner.loop['finite_N']['delta_E']),('delta_V',owner.loop['finite_N']['delta_V'])]+[
        ('density_'+key,ids) for key,ids in owner.loop['finite_N']['five_signed_density_increments'].items()]:
        for n,node in enumerate(ids):zero(alg.at(node),('actual original flat loop source',key,n));checked.append([key,n,node])
    for key,ids in owner.loop['incoming'].items():assert all(node==0 for node in ids)
    for key,ids in owner.loop['history'].items():
        for node in ids:zero(alg.at(node),('actual empty/zero finite-N local history',key))
    for key in ('U','m'):
        wanted=4*Z+v['j']+v['eps']*sy.Function('retained_collar_'+('Psi' if key=='U' else 'mean'))(v['s'],Z)
        zero(alg.at(owner.physical[key][0].value)-wanted,('actual modified/background decomposition',key))
        zero(alg.at(owner.physical[key][1].value)-alg.at(owner.reconstructed_y[key][1].value),('actual native first-y prefix',key))
    zero(alg.at(owner.radius)-v['eps']*v['rho'],'actual whole active physical radius')
    assert alg.evidence['flat']
    return dict(passed=True,actual_flat_loop_and_five_finite_N_density_rows=checked,
        actual_zero_leading_correction_inlet_and_histories_not_own_incoming_data=True,
        actual_modified_U_m_and_first_native_y_rows_identified=True,
        actual_background_six_nonzero_histories_and_core_IVP_retained=True)


def actual_leading_high_radial_source_join(owner):
    env_left,vleft=lower_check.leading_environment(owner.lower,True)
    env_right,vright=active_context(owner)
    # Supply the same actual axis gauge in both views before comparing the
    # independently constructed collar derivatives under rho(s).
    bind(env_left,owner.jet(owner.bridge['core7']['U0']),4*vleft['Z']+vright['j'],vleft['Z'])
    env_left[owner.parameters['j'].node]=vright['j']
    env_left[owner.parameters['delta'].node]=vright['delta']
    left=RminusAlgebra(owner,env_left);right=RminusAlgebra(owner,env_right)
    a=current.outer.Algebra(owner.g,5);qalg=current.outer.Algebra(owner.g,4)
    evidence=[]
    for key in ('U','m','Q'):
        alg=qalg if key=='Q' else a
        for n in range(5):
            lower=current.y_row(alg,owner.g,owner.lower.rho,owner.lower.radial[key],n)
            active=owner.physical[key][n]
            zero(vleft['normalize'](left.at(lower.value))-vright['normalize'](right.at(active.value)),('actual distinct high radial sources',key,n))
            evidence.append(dict(field=key,physical_y_order=n,left_node=lower.value.node,right_node=active.value.node))
    # The same physical-Y linear leading recovery defines every higher Q;
    # the original Q/Q_y prefix stays unchanged.
    for n in range(5):zero(right.at(owner.physical['Q'][n].value)-right.at(owner.recovered_Q[n].value),('actual active leading Q recovery',n))
    return dict(passed=True,actual_distinct_leading_U_m_Q_y0_to_y4_whole_source_identities=evidence,
        new_active_high_rows_use_actual_micro_ratio_U_mean_and_core_not_lower_outputs=True,
        original_common_IVP_and_two_sided_flat_source_consumed=True,
        ordinary_Z_C5_Q_C4_by_actual_source_and_product_composition=True,
        leading_one_minus_delta_recovery_preserved=True)


def actual_pressure_operator_and_FTC(owner):
    env,v=active_context(owner,False);Z=v['Z'];y=sy.Symbol('independent_active_log_radius',real=True)
    R0=sy.Symbol('actual_R0',positive=True);R=R0*sy.exp(y);env[owner.radius.node]=R
    u=sy.Function('actual_leading_U')(y,Z);q=sy.Function('actual_leading_Q')(y,Z)
    for key,body in (('U',u),('Q',q)):
        for n,jet in enumerate(owner.physical[key]):bind(env,jet,sy.diff(body,y,n),Z)
    alg=RminusAlgebra(owner,env);de=v['delta'];d=1-Z**2;L=1-de*Z**2;h=1-de
    S=(q+h*Z*sy.diff(q,Z)/2+sy.diff(q,y))/L+q*(q/2+sy.diff(q,y))+u*(-2*Z*q+d*sy.diff(q,Z)-2*Z*sy.diff(q,y))/L
    for n,node in enumerate(owner.modified['S']):
        zero(alg.at(node)-sy.diff(S,Z,n),('actual saved modified reduced-S source',n))
    zero(alg.at(owner.density['S'].value)-S,'actual new density uses the same saved reduced-S source')
    N=-S/2+(sy.diff(q,y)+sy.diff(q,y,2))/R
    for n,jet in enumerate(owner.N_y):zero(alg.at(jet.value)-sy.diff(N,y,n),('actual active N physical-y derivative',n))
    for n in range(3):
        zero(alg.at(owner.pressure_y[n].value)-sy.diff(R*N,y,n),('actual own pressure-y source',n+1))
        zero(alg.at(owner.raw_pressure_y[n].value)-alg.at(owner.pressure_y[n].value),('actual original/coupled pressure source',n+1))
        zero(alg.at(owner.m15_y[n].value)-sy.diff(-R**2*N,y,n),('actual signed m15-y source',n+1))
        zero(alg.at(owner.raw_m15_y[n].value)-alg.at(owner.m15_y[n].value),('actual original/coupled m15 source',n+1))
    # Differentiate the saved original primitives as whole functions. No
    # pressure output or derivative is bound to a desired join value.
    env,v=active_context(owner,False);s=v['s'];Z=v['Z'];H=v['H'];R=v['eps']*v['rho']
    qs=sy.Function('original_active_Q')(s,Z);ss=sy.Function('original_reduced_S')(s,Z)
    bind(env,owner.jet(owner.modified['physical']['Q']),qs,Z)
    bind(env,owner.jet(owner.modified['physical']['Q_y']),sy.diff(qs,s)/H,Z)
    bind(env,owner.jet(owner.modified['S']),ss,Z)
    alg=RminusAlgebra(owner,env);original_pressure={key:alg.at(ids[0]) for key,ids in owner.modified['pressure'].items()}
    Py=-R*ss/2+sy.diff(qs,s)/H+sy.diff(qs,s,2)/H**2
    zero(sy.diff(original_pressure['m13'],s).doit()/H-Py,'actual saved first-micro m13 FTC')
    zero(sy.diff(original_pressure['m15'],s).doit()/H+R*Py,'actual saved first-micro signed m15 FTC')
    recipe=[]
    for key,ids in owner.modified['local_integrals'].items():
        for n,node in enumerate(ids):
            row=owner.g.nodes[node]
            assert row['operation']=='definite_integral' and row['lower']==owner.right_endpoint.node and row['upper']==owner.coordinate.node
            assert row['ordinary_slow_Z_derivative_order']==n and row['endpoint_independent_of_Z']
            assert row['physical_Jacobian_applied_once'] and row['physical_radius_power']==(1 if key=='S0' else 2)
            recipe.append([key,n,node,row['integrand']])
    return dict(passed=True,actual_pressure_N_y0_y1_y2_and_own_P_y1_y2_y3_source_derivatives=True,
        actual_original_reduced_S_and_signed_pressure_FTC=True,
        actual_signed_integral_C3_recipes=recipe,original_hb_divides_native_derivative_once=True,
        common_high_leading_derivative_theorem_identifies_lower_and_active_fixed_N_sources=True,
        whole_own_pressure_outputs_not_bound=True)


def actual_own_flow_and_rminus_join(owner):
    env,Z,slow=rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
    N0=env[owner.lower.N[0].value.node]
    N1=sy.Function('common_fixed_N_rho1')(slow,Z);N2=sy.Function('common_fixed_N_rho2')(slow,Z)
    bind(env,owner.lower.N[1],N1,Z);bind(env,owner.lower.N[2],N2,Z)
    env[owner.active_sc.node]=RminusAlgebra(owner,env).at(owner.parameters['sc'])
    active=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True)
    alg0=RminusAlgebra(owner,env);rho=alg0.at(owner.rho)
    for jet,body in zip(owner.N_y,(N0.subs(slow,active),rho*N1.subs(slow,active),rho*N1.subs(slow,active)+rho**2*N2.subs(slow,active))):bind(env,jet,body,Z)
    bind(env,owner.jet(owner.modified['physical']['Q']),sy.Function('retained_boundary_Q')(active,Z),Z)
    bind(env,owner.jet(owner.modified['physical']['Q_y']),sy.Function('retained_boundary_Qy')(active,Z),Z)
    alg=RminusAlgebra(owner,env);f=owner.right['fields'];mean=alg.at(f['mean1'][0].value);R=alg.at(owner.radius)
    for key in ('M1','mean1','Q1','V1'):
        for n,(new,old) in enumerate(zip(rows(f[key][1]),owner.intermediate['fields'][key+'_y'])):
            zero(alg.at(new)-alg.at(old),('retained actual intermediate first-y row',key,n))
    de=alg.at(owner.parameters['delta']);d=1-Z**2;L=1-de*Z**2
    ownQ=(-(1+de)*Z*mean-d*sy.diff(mean,Z))/L
    zero(alg.at(f['K1'][0].value)-mean,'actual first-micro own K')
    zero(alg.at(f['Q1'][0].value)-ownQ,'actual own one-plus-delta Q')
    zero(alg.at(f['M1'][0].value)-R*mean,'actual own mass primitive')
    zero(alg.at(f['V1'][0].value)-R*ownQ,'actual own V=RQ')
    y=sy.Symbol('independent_own_log_radius',real=True)
    model_mean=mean*sy.exp(-y);model_R=R*sy.exp(y)
    model_Q=(-(1+de)*Z*model_mean-d*sy.diff(model_mean,Z))/L
    model=dict(Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),mean1=model_mean,K1=model_mean,
        Q1=model_Q,M1=model_R*model_mean,V1=model_R*model_Q)
    flow=[]
    for key,body in model.items():
        for j,jet in enumerate(f[key]):zero(alg.at(jet.value)-sy.diff(body,y,j).subs(y,0),('actual first-micro own flow',key,j));flow.append([key,j])
    for key in ('m11','m12','m14'):
        for j,jet in enumerate(owner.right['cumulative'][key][1:],1):zero(alg.at(jet.value),('actual constant own moment',key,j))
    # Check the complete lower partial equals the retained full integral
    # at sc_half; preserve this nonzero input when the active increment is0.
    evidence=[]
    for spec in owner.contract:
        for j in range(spec['physical_y_order']+1):
            left=rows(owner.join['left'][spec['group']][spec['field']][j]);right=rows(owner.join['right'][spec['group']][spec['field']][j])
            for n,(l,r) in enumerate(zip(left,right)):
                zero(alg.at(l)-alg.at(r),(spec['group'],spec['field'],j,n))
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=j,ordinary_Z_order=n,
                    left_node=l.node,right_node=r.node,proof='actual own recovery, nonzero retained histories and fixed pressure FTC'))
    assert len(evidence)==246
    return dict(passed=True,actual_rminus_contracted_mixed_join_rows=evidence,actual_own_state_flow_rows=flow,
        original_own_nonzero_moment_mass_pressure_incoming_data_preserved=True,
        empty_first_micro_relative_pressure_integrals_and_boundary_increments=True,
        actual_common_fixed_N_is_input_not_an_imposed_own_pressure_output=True,
        leading_finite_N_zero_incoming_not_reused_as_own_n1_data=True)


CHECKS=(source_contract,actual_flat_first_micro_and_native_prefix,actual_leading_high_radial_source_join,
    actual_pressure_operator_and_FTC,actual_own_flow_and_rminus_join)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('RMINUS_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['rminus_mixed_source_join_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_rminus_contracted_mixed_join_rows=246,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_OWN_RMINUS246_SOURCE_ACCEPTED',flush=True)
    return rec
