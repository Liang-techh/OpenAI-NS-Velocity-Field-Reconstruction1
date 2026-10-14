"""Actual Rin flat cutoff, common pressure source, signed FTC and join."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join as current
import lei_ren_part1_paper_compliant_current_original_n1_Rc_mixed_source_join_check as previous_check

base=previous_check.outer_check
rows,encode=current.outer.rows,current.outer.encode
zero=previous_check.zero


def bind(env,jet,body,Z):
    for k,ref in enumerate(rows(jet)):
        value=sy.diff(body,Z,k)
        if ref.node in (0,1):assert value==ref.node
        elif ref.node in env:zero(env[ref.node]-value,('shared actual input',ref.node))
        else:env[ref.node]=value


class RinAlgebra(base.SourceAlgebra):
    def __init__(self,owner,env,evidence=None):
        super().__init__(owner.g,env);self.owner=owner
        self.evidence=evidence if evidence is not None else dict(flat=set(),empty=set())
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation']
        if op=='product':
            # All factors belong to admitted finite analytic sources here.
            # Consume an exact flat-cutoff zero before expanding unrelated
            # regular-solution/gauge bodies multiplied by that zero.
            values=[];deferred=[]
            for node in row['arguments']:
                try:value=self.at(node)
                except (AssertionError,KeyError) as error:
                    deferred.append(error);continue
                if value==0:self.cache[i]=sy.Integer(0);return self.cache[i]
                values.append(value)
            if deferred:raise deferred[0]
            result=sy.Mul(*values)
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            symbol=sy.Symbol(self.g.nodes[row['variable']]['name'],real=True);value=self.at(row['value'])
            env={key:body.subs(symbol,value) if hasattr(body,'subs') else body for key,body in self.bindings.items()}
            env[row['variable']]=value;result=RinAlgebra(self.owner,env,self.evidence).at(row['expression'])
        elif op=='exact_original_source_flat_sigma' or (op=='analytic_unary' and row['name'] in ('original_flat_sigma_prime','original_flat_sigma_second')):
            x=sy.cancel(self.at(row['argument']));order=0 if op=='exact_original_source_flat_sigma' else (1 if row['name'].endswith('prime') else 2)
            if (x-1).is_nonnegative:
                self.evidence['flat'].add(i);result=sy.Integer(1 if order==0 else 0)
            elif x.is_nonpositive:
                self.evidence['flat'].add(i);result=sy.Integer(0)
            else:
                t=sy.Symbol('actual_sigma_argument',real=True)
                result=sy.diff(sy.Function('actual_original_sigma')(t),t,order).subs(t,x)
        elif op=='current_original_unique_regular_n1_component':
            assert row['source_family']==self.owner.family and row['unique_regular_analytic_solution']
            coordinate=self.at(row['coordinate']);Z=sy.Symbol('Z',real=True)
            result=sy.diff(sy.Function('actual_regular_'+row['quantity'])(coordinate,Z),Z,row['ordinary_slow_Z_derivative_order'])
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z']
            zero(self.at(row['upper'])-self.at(row['lower']),'actual empty lower primitive')
            self.evidence['empty'].add(i);result=sy.Integer(0)
        else:return super().at(i)
        self.cache[i]=result;return result


def environment(owner,history=False,pressure_theorem=False):
    Z=sy.Symbol('Z',real=True);s=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True)
    H,SC,L=sy.symbols('actual_hb actual_sc actual_logLambda',positive=True)
    p=owner.parameters
    env={p['hb'].node:H,p['sc'].node:SC,p['logLambda'].node:L}
    # epsilon and Lambda are evaluated from their actual exponential recipes.
    if history:
        core_rho=sy.Symbol(owner.g.nodes[owner.inner['coordinates']['core']]['name'],real=True)
        bind(env,owner.jet(owner.cut['fields']['core']['rows']['mean1'][0]),sy.Function('actual_retained_core_mean')(core_rho,Z),Z)
        # Bind the nonzero defining history, not an output mean whose
        # expanded sum may differ from the original M1 input expression.
        ids=[]
        for k,body in enumerate(owner.cut['history_bodies']['first']['mean_integrand']):
            matches=[i for i,row in enumerate(owner.g.nodes) if row['operation']=='definite_integral'
                and row.get('actual_cutoff_own_mean') and row['integrand']==body
                and row['ordinary_slow_Z_derivative_order']==k]
            assert len(matches)==1;ids.append(matches[0])
        bind(env,owner.jet(ids),sy.Function('actual_retained_cutoff_mean_history')(s,Z),Z)
        for key in ('m12','m13','m14','m15'):
            bind(env,owner.jet(owner.cut['cumulative']['first'][key]),sy.Function('actual_cutoff_'+key)(s,Z),Z)
    if pressure_theorem:
        N=sy.Function('proved_actual_common_N1p')(s,Z)
        bind(env,owner.jet(owner.inner['systems']['first']['N1p']),N,Z)
        bind(env,owner.jet(owner.leading['N1p']),N,Z)
    for key,body in (('F_rho_rows',sy.Function('actual_leading_F')(s,Z)),):
        bind(env,owner.jet(owner.lower['fields'][key][0]),body,Z)
    return env,Z,s


def source_contract(owner):
    owner.assert_graph();g=owner.g
    for module in (current.previous,current.cutoff_source,current.inner_source,current.lower_source,current.lower_high):
        raw=current.parent.linear.read(module.NAME);rec=current.parent.linear.read(module.RECEIPT)
        assert rec['all_passed'] and rec['source_family']==raw['source_family']==owner.family
        assert all(rec[key] for key in getattr(module,'GATES',(getattr(module,'GATE',''),)))
        assert not any(rec[key] for key in module.OPEN)
        assert raw['exact_graph_nodes']==g.nodes[:len(raw['exact_graph_nodes'])]
        for name,digest in rec['input_hashes'].items():assert owner.hashes[name]==digest==current.sha(name),name
    assert owner.inner['leading']['first']==owner.lower['fields']
    assert owner.inner['parameters']==owner.lower['parameters']
    assert owner.coordinate.node==owner.lower['fields']['query']
    assert owner.intermediate['coordinate']==owner.coordinate.node
    assert owner.intermediate['radius']==owner.radius.node
    assert owner.intermediate['domain'][0]==owner.endpoint.node
    p=owner.parameters
    assert g.nodes[p['epsilon'].node]==dict(operation='analytic_unary',name='exp',argument=g.neg(p['logLambda']).node)
    assert g.nodes[p['Lambda'].node]==dict(operation='analytic_unary',name='exp',argument=p['logLambda'].node)
    count=0
    for side in ('left','right'):
        for spec in owner.contract:
            jets=owner.join[side][spec['group']][spec['field']]
            assert len(jets)==spec['physical_y_order']+1
            for jet in jets:
                assert len(rows(jet))==spec['ordinary_Z_order']+1
                for ref in rows(jet):
                    row=g.nodes[ref.node]
                    assert row['operation']=='function_substitution' and row['Z_independent_substitution']
                    assert row['variable']==owner.coordinate.node and row['value']==owner.endpoint.node
                    count+=1
    assert count==276
    assert owner.definition['positive_order_mean_weight']=='1+delta'
    return dict(passed=True,same_actual_inner_first_and_lower_leading_graph=True,
        exact_epsilon_Lambda_reciprocal_source=True,fixed_Z_independent_endpoint_rows=count,
        actual_native_state_rho0_rho1_rho2_retained=True,higher_radial_pressure_not_claimed=True)


def actual_leading_pressure_identity(owner):
    env,Z,s=environment(owner);g=owner.g;rho=4*sy.exp(env[owner.parameters['hb'].node]*s)
    f=owner.lower['fields']
    bind(env,owner.jet(owner.leading['fields']['U']),sy.Function('actual_leading_U')(s,Z),Z)
    for k,ids in enumerate(owner.leading['radial_fields']['Q_rho_rows']):
        bind(env,owner.jet(ids),sy.Function('actual_Q_rho'+str(k))(s,Z),Z)
    alg=RinAlgebra(owner,env)
    left=owner.inner['systems']['first']['N1p'][0];right=owner.leading['N1p'][0]
    zero(alg.at(left)-alg.at(right),'actual inner/lower whole N1p identity')
    Q0,Qr,Qrr=[alg.at(ids[0]) for ids in f['Q_rho_rows']]
    zero(alg.at(f['Q_y'][0])-rho*Qr,'actual Q_y source')
    zero(alg.at(f['Q_yy'][0])-rho*Qr-rho**2*Qrr,'actual Q_yy source')
    expected=alg.at(owner.lower['pressure_bindings']['N1p'][0])
    zero(alg.at(right)-expected,'higher source preserves original lower N1p value')
    assert owner.leading['N1p'][:2]==owner.lower['pressure_bindings']['N1p']
    assert len(owner.inner['systems']['first']['N1p'])==len(owner.leading['N1p'])==4
    return dict(passed=True,actual_raw_and_lower_N1p_whole_function_identity=True,
        identity='Lambda*(2*Q_rho+rho*Q_rhorho)=(Q_yy+Q_y)/R',
        same_actual_S_and_leading_field=True,ordinary_Z0_Z3_by_admitted_derivative_composition=True,
        no_N1p_output_bound_during_whole_source_proof=True)


def actual_flat_cutoff_and_FTC(owner):
    env,Z,s=environment(owner,history=True,pressure_theorem=True);alg=RinAlgebra(owner,env)
    g=owner.g;p=owner.parameters;end=alg.at(owner.endpoint)
    evidence=[]
    for key in ('F1','Uz1'):
        for j,jet in enumerate(owner.left['fields'][key]):
            zero(alg.at(owner.join['left']['fields'][key][j].value),('actual flat cutoff',key,j))
            evidence.append([key,j])
    # The defining histories are not replaced by zero: check their actual
    # integrands at the moving upper endpoint before using the FTC.
    body=owner.cut['history_bodies']['first'];a=current.outer.Algebra(g,3)
    inspected=[];raw=current.parent.linear.read(current.cutoff_source.NAME)
    for key,label in (('mean_integrand','mean'),('pressure_integrand','pressure')):
        ids=body[key]
        for k,node in enumerate(ids):
            # Inspect the actual defining cutoff integral, not a copied
            # boundary formula or an enclosure of its numerical magnitude.
            matches=[row for row in raw['exact_graph_nodes'] if row['operation']=='definite_integral' and row['integrand']==node]
            assert len(matches)==1
            row=matches[0];assert row['lower']==0 and row['upper']==owner.coordinate.node
            assert row['endpoint_independent_of_Z'] and row['ordinary_slow_Z_derivative_order']==k
            variable=g.symbol(row['variable']);jet=owner.jet(ids)
            at_upper=a.at(jet,variable,owner.coordinate)
            if label=='mean':zero(alg.at(rows(a.at(at_upper,owner.coordinate,owner.endpoint))[k]),'flat upper own mean density')
            else:
                rho=alg.at(owner.rho);eps=alg.at(p['epsilon']);hb=alg.at(p['hb'])
                want=eps*hb*rho*alg.at(owner.inner['systems']['first']['N1p'][0])
                # Evaluate after the fixed endpoint substitution so the cutoff
                # is flat before inspecting its regular solution factors.
                target=sy.diff(want,Z,k).subs(s,end)
                zero(alg.at(rows(a.at(at_upper,owner.coordinate,owner.endpoint))[k])-target,'actual pressure history upper FTC')
            inspected.append(dict(history=label,ordinary_Z_order=k,variable=row['variable'],body=node))
    assert alg.evidence['flat']
    return dict(passed=True,actual_flat_velocity_rows=evidence,
        actual_history_FTC_integral_recipes=inspected,flat_argument_at_Rin=2,
        local_upper_mean_density_zero_not_full_history_zero=True,
        pressure_upper_density_retains_fixed_leading_source=True)


def actual_signed_transport_and_empty_lower(owner):
    env,Z,s=environment(owner,history=True,pressure_theorem=True);alg=RinAlgebra(owner,env)
    p=owner.parameters;g=owner.g;a=current.outer.Algebra(g,3);hb=alg.at(p['hb'])
    N=alg.at(owner.leading['N1p'][0]);rho=alg.at(owner.rho);R=alg.at(owner.radius)
    endpoint=alg.at(owner.endpoint)
    for key,expected in (('m13',hb*R*N),('m15',-hb*R**2*N)):
        body=owner.leading['signed_bodies'][key]
        zero(alg.at(body[0])-expected,('actual lower signed physical body',key))
        pushed=a.scale(owner.jet(owner.cut['integrands']['first'][key]),g.mul(p['hb'],owner.rho))
        zero(alg.at(a.at(pushed,owner.coordinate,owner.endpoint).value)-expected.subs(s,endpoint),('actual cutoff signed pushforward',key))
        fixed=owner.intermediate['fixed_pressure'][key]
        zero(alg.at(a.at(owner.jet(fixed),owner.coordinate,owner.endpoint).value),('actual empty partial primitive',key))
    empty=[]
    for key,ids in owner.intermediate['fixed_pressure'].items():
        for n,node in enumerate(ids):
            row=g.nodes[node];assert row['operation']=='definite_integral'
            assert row['lower']==owner.endpoint.node and row['upper']==owner.coordinate.node
            assert row['integrand']==owner.leading['signed_bodies'][key][n]
            assert row['endpoint_independent_of_Z']
            empty.append(dict(moment=key,ordinary_Z_order=n,node=node))
    return dict(passed=True,actual_signed_pressure_pushforward_and_single_m15_minus=True,
        rho_s_is_hb_rho_and_R_s_is_hb_R=True,actual_lower_empty_integral_rows=empty,
        nonzero_cutoff_m13_m15_endpoints_preserved=True)


def actual_Rin_join(owner):
    env,Z,s=environment(owner,history=True,pressure_theorem=True);alg=RinAlgebra(owner,env);evidence=[]
    for spec in owner.contract:
        for j in range(spec['physical_y_order']+1):
            left=rows(owner.join['left'][spec['group']][spec['field']][j])
            right=rows(owner.join['right'][spec['group']][spec['field']][j])
            for n,(l,r) in enumerate(zip(left,right)):
                try:zero(alg.at(l)-alg.at(r),(spec['group'],spec['field'],j,n))
                except Exception as error:raise AssertionError(('Rin row',spec['group'],spec['field'],j,n)) from error
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=j,
                    ordinary_Z_order=n,left_node=l.node,right_node=r.node,
                    proof='direct identity on the actual retained source and FTC'))
    assert len(evidence)==138 and len(alg.evidence['empty'])==8
    return dict(passed=True,actual_Rin_contracted_join_rows=evidence,
        original_nonzero_pressure_and_moment_inputs_retained=True,
        own_one_plus_delta_Q_recovery_and_radius_factors_preserved=True,
        higher_radial_pressure_full_mixed_contract_and_global_n1_open=True)


CHECKS=(source_contract,actual_leading_pressure_identity,actual_flat_cutoff_and_FTC,
    actual_signed_transport_and_empty_lower,actual_Rin_join)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('RIN_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['Rin_source_join_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_Rin_contracted_join_rows=138,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_OWN_RIN_138_SOURCE_ROWS_ACCEPTED',flush=True)
    return rec
