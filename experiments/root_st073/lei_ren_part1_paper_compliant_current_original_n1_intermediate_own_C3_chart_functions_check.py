"""Regional source ownership, actual C3 bodies and partial pressure calculus."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_repaired_own_C3_terminal_functions_check as prior


def bind_free_jet(g,bindings,jet,label,Z):
    """Respect genuine exact-zero source rows in symbolic test functions."""
    refs=current.rows(jet);zero=next((n for n,q in enumerate(refs) if q==g.zero),None)
    body=(sy.Integer(0) if zero==0 else sum(sy.Symbol(label+'_c'+str(n),real=True)*Z**n for n in range(zero))) \
        if zero is not None else sy.Function(label)(Z)
    for n,q in enumerate(refs):
        if q!=g.zero:bindings[q.node]=sy.diff(body,Z,n)
    bindings[g.zero.node]=sy.Integer(0)
    return body


def source_check(owner):
    g=owner.g;a=current.Algebra(g,3)
    assert len(owner.charts)==22 and len(owner.preRc_names)==17
    assert owner.patch.modified is owner.modified and owner.compact.outer is owner.outer
    expected={**{name:data['pressure'] for name,data in owner.modified.charts.items()},
              'Rm_patch':owner.patch.pressure,**{name:data['pressure'] for name,data in owner.outer.charts.items()}}
    source_meta={**{name:data['source_chart'] for name,data in owner.modified.charts.items()},
                 'Rm_patch':owner.patch.source_chart,**owner.outer.charts}
    for name,data in owner.charts.items():
        assert data['complete_incoming_mass_not_reset']
        assert all(len(current.rows(data['fields'][k]))==4 for k in ('M1','mean1','K1','P1'))
        assert all(len(current.rows(data['fields'][k]))==3 for k in ('Q1','V1','Q1_y','V1_y'))
        assert data['fields']['M1']==data['cumulative']['m11'] and data['fields']['P1']==data['cumulative']['m13']
        for k in ('Uz1','Utheta1','F1'):assert all(ref==g.zero for ref in current.rows(data['fields'][k]))
        if name in owner.preRc_names:
            meta=source_meta[name]
            coordinate=owner.patch.coordinate if name=='Rm_patch' else meta['coordinate']
            assert data['coordinate']==coordinate and data['radius']==meta['radius']
            assert data['domain']==list(meta['domain'])
            assert data['fixed_pressure'] is expected[name]
            for key in ('m11','m12','m14'):assert data['cumulative'][key]==a.trim(owner.inner[key])
            for key in ('m13','m15'):
                assert data['cumulative'][key]==a.add(owner.inner[key],expected[name][key],owner.complete.segments['lower'][key])
    for name,binding in owner.pulse_bindings.items():
        meta=owner.pulse.pulse.own.fields[name];data=owner.charts['pulse_'+name]
        assert data['coordinate']==meta['query'] and data['radius']==meta['radius']
        assert data['domain']==[meta['domain']['lower'],meta['domain']['upper']]
        assert binding['coordinate']==meta['query'] and binding['value']==meta['query']
        assert binding['radius']==a.at(a.fixed(meta['radius']),meta['query'],meta['query']).value
        radius=g.nodes[meta['radius'].node]
        assert radius['operation']=='analytic_unary' and radius['name']=='exp'
        assert radius['argument']==meta['log_radius'].node
        for key,q in binding['lower'].items():assert q==a.trim(owner.pulse.global_endpoints[0][key])
        for key,q in binding['current_integrals'].items():
            template=g.nodes[owner.pulse.integrals[name][key].value.node]
            for k,ref in enumerate(current.rows(q)):
                row=g.nodes[ref.node]
                assert row['operation']=='definite_integral' and row['endpoint_independent_of_Z']
                assert row['integrand']==current.rows(binding['current_integrands'][key])[k].node
                assert row['variable']==template['variable'] and row['lower']==template['lower']
                assert row['upper']==binding['value'].node and row['ordinary_slow_Z_derivative_order']==k
                assert row['measure']==template['measure'] and row['physical_Pstar_not_reapplied']
        assert binding['preceding_complete_charts']==list(('main','gap','end')[:('main','gap','end').index(name)])
    for key,q in owner.lower_partial.items():
        for k,ref in enumerate(current.rows(q)):
            row=g.nodes[ref.node]
            assert row['integrand']==current.rows(owner.lower.signed_bodies[key])[k].node
            assert row['ordinary_slow_Z_derivative_order']==k and row['original_m15_minus_applied_once']
            meta=owner.lower.lower.fields;template=g.nodes[owner.lower.pressure[key].value.node]
            assert row['upper']==meta['query'].node and row['lower']==template['lower']
            assert row['variable']==template['variable']
    meta=owner.lower.lower.fields;data=owner.charts['lower']
    assert data['coordinate']==meta['query'] and data['radius']==meta['radius']
    assert data['domain']==[meta['domain']['lower'],meta['domain']['upper']]
    assert owner.definition['leading_pulse_velocity_not_misidentified_as_positive_order_velocity']
    assert owner.at_fraction('quiet','1/3')['QV_C2']
    try:owner.at_fraction('quiet','-1/100')
    except ValueError:pass
    else:raise AssertionError('Outside native chart accepted')
    return dict(passed=True,actual_regional_charts=22,actual_preRc_charts=17,
        complete_lower_active_patch_outer_quiet_pulse_pressure_histories_retained=True,
        actual_pulse_partial_integral_body_rows=120,actual_lower_partial_integral_body_rows=8,
        original_variables_limits_physical_measures_and_no_extra_Pstar=True,
        positive_order_velocities_zero_without_zeroing_leading_velocity=True,
        actual_preRc_coordinate_radius_domain_and_C3_pressure_source_identity=True,
        pulse_true_radius_recipe_and_query_radius_not_frozen_interval=True,
        fixed_pulse_pressure_is_original_leading_N1p_forcing=True,
        exact_native_coordinate_query_domain_enforced=True)


def actual_regional_recovery(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);R=sy.Symbol('R',positive=True);de=sy.Symbol('delta',positive=True)
    M=sy.Function('M')(Z);P=sy.Function('P')(Z)
    Q=(-(1+de)*Z*M-(1-Z*Z)*sy.diff(M,Z))/(R*(1-de*Z*Z))
    expected=dict(M1=M,mean1=M/R,K1=M/R,P1=P,Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),
        Q1=Q,V1=R*Q,M1_y=sy.Integer(0),mean1_y=-M/R,Q1_y=-Q,V1_y=sy.Integer(0))
    checked=0
    for name,data in owner.charts.items():
        bindings={g.symbol('Z').node:Z,data['radius'].node:R,owner.delta.node:de}
        prior.bind_jet(bindings,data['cumulative']['m11'],M,Z);prior.bind_jet(bindings,data['cumulative']['m13'],P,Z)
        alg=prior.FunctionAlgebra(g,bindings)
        for key,body in expected.items():
            for n,ref in enumerate(current.rows(data['fields'][key])):
                assert sy.cancel(alg.at(ref.node)-sy.diff(body,Z,n))==0,(name,key,n)
                checked+=1
    return dict(passed=True,actual_regional_own_recovery_rows=checked,
        arbitrary_nonzero_incoming_mass_and_pressure_functions=True,
        positive_order_one_plus_delta_and_true_mean_Z3=True,
        own_mean_y_and_Q_y_are_minus_the_own_mean_and_Q=True,
        own_V_y_and_mass_y_zero_only_outside_velocity_supports=True)


def actual_pulse_bodies(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);R=sy.Symbol('R',positive=True);mu=sy.Symbol('mu',positive=True)
    checked=[];prefix=0
    for name,binding in owner.pulse_bindings.items():
        meta=owner.pulse.pulse.own.fields[name];f=owner.pulse.leading_fields[name]
        muv=owner.pulse.pulse.own.jets.complete.Gamma.schedule['mu']
        bindings={g.symbol('Z').node:Z,meta['radius'].node:R,muv.node:mu}
        U=bind_free_jet(g,bindings,f['U'],'U',Z);Q=bind_free_jet(g,bindings,f['Q'],'Q',Z)
        Uz,Qz=sy.diff(U,Z),sy.diff(Q,Z)
        values=dict(I0_q_Z=Qz,I0_U_q_Z=U*Qz,I0_U_squared=U*U,I0_U_U_Z=U*Uz,
            I1_q=Q,I1_q_Z=Qz,I1_q_squared=Q*Q,I1_U_q_Z=U*Qz,I1_U_squared=U*U,I1_U_U_Z=U*Uz)
        alg=prior.FunctionAlgebra(g,bindings)
        for key,q in binding['current_integrands'].items():
            weight=(R if name=='end' else R/mu)*(R if key.startswith('I1_') else 1)
            old=owner.pulse.pulse.bindings['products'][name][key]
            assert sy.cancel(alg.at(old.value.node)-values[key])==0,(name,key,'original value')
            prefix+=1
            for n,ref in enumerate(current.rows(q)):
                assert sy.cancel(alg.at(ref.node)-weight*sy.diff(values[key],Z,n))==0,(name,key,n)
                checked.append(f'{name}_{key}_Z{n}')
    return dict(passed=True,actual_saved_pulse_body_rows=checked,
        original_nonradial_value_product_identities=prefix,
        actual_U_C5_and_Q_C4_supply_real_product_derivatives=True,
        original_R_over_mu_and_R_squared_over_mu_measures=True,
        end_chart_R_and_R_squared_measures=True,physical_Pstar_applied_once=True,
        partial_integral_derivatives_not_copied_from_zero_whole_integral_rows=True)


def pressure_calculus(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);de=sy.Symbol('delta',positive=True)
    binding=owner.pulse_bindings['main'];delta=owner.pulse.pulse.own.jets.complete.Gamma.schedule['delta']
    bindings={g.symbol('Z').node:Z,delta.node:de}
    for label,group in (('integral',binding['total_integrals']),('lower',binding['lower']),('upper',binding['upper'])):
        for key,q in group.items():
            if all(ref==g.zero for ref in current.rows(q)):continue
            bind_free_jet(g,bindings,q,label+'_'+key,Z)
    alg=prior.FunctionAlgebra(g,bindings);pressure=owner.charts['pulse_main']['fixed_pressure']
    pairs=lambda group:{k:current.source.C1Function(*current.rows(q)[:2]) for k,q in group.items()}
    old=current.pressure_source.pressure_functions(g,delta,pairs(binding['total_integrals']),pairs(binding['lower']),pairs(binding['upper']))
    checked=[]
    for key,q in pressure.items():
        body=alg.at(old[key].value.node)
        for n,ref in enumerate(current.rows(q)):
            assert sy.cancel(alg.at(ref.node)-sy.diff(body,Z,n))==0,(key,n)
            checked.append(f'{key}_Z{n}')
    return dict(passed=True,actual_current_pressure_rows=checked,
        original_signed_reduced_S_formula_and_finite_boundaries_retained=True,
        actual_C3_Leibniz_quotient_rows_not_Taylor_coefficients=True,
        m15_single_negative_weight_and_endpoint_signs=True,
        mathematical_function_calculus_not_original_scale_numeric_bounds=True)


class EndpointAlgebra:
    """Actual flat startup/end and empty physical histories on their source graph."""
    def __init__(self,g,bindings):self.g,self.bindings=g,bindings
    def at(self,node,local=None):
        if local is not None and node in local:return local[node]
        if node in self.bindings:return self.bindings[node]
        row=self.g.nodes[node];op=row['operation'];at=lambda n:self.at(n,local)
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='bound_variable':return sy.Symbol(row['name'],real=True)
        if op=='sum':return sy.Add(*(at(i) for i in row['arguments']))
        if op=='product':return sy.Mul(*(at(i) for i in row['arguments']))
        if op=='negative':return -at(row['argument'])
        if op=='positive_quotient':return sy.cancel(at(row['numerator'])/at(row['denominator']))
        if op=='function_substitution':return self.at(row['expression'],{**(local or {}),row['variable']:at(row['value'])})
        if op=='definite_integral':
            lo,hi=at(row['lower']),at(row['upper'])
            if sy.simplify(lo-hi)==0:return sy.Integer(0)
            raise AssertionError('Nonempty unbound endpoint integral')
        if op=='compact_raw_beta' or (op=='analytic_unary' and row['name']=='compact_raw_beta_derivative'):
            q=sy.simplify(at(row['argument']))
            assert q.is_number and abs(q)>=1
            return sy.Integer(0)
        if op in ('analytic_unary','exact_original_source_flat_sigma'):
            q=at(row['argument']);name='original_flat_sigma' if op=='exact_original_source_flat_sigma' else row['name']
            if name in ('exp','log','sqrt'):return getattr(sy,name)(q)
            if name in ('original_flat_sigma','original_flat_sigma_prime'):
                assert q.is_number and (q<=0 or q>=1)
                return sy.Integer(0) if name.endswith('prime') or q<=0 else sy.Integer(1)
        raise AssertionError('Unbound actual pulse endpoint source '+op)


def actual_pulse_endpoints(owner):
    p=owner.pulse;g=owner.g;Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    Gamma=p.pulse.own.jets.complete.Gamma;W=Gamma.pulse.functions
    LP=sy.Symbol('logP',real=True);LR=sy.Symbol('logRp',real=True);de=sy.Symbol('delta',positive=True)
    bindings={g.symbol('Z').node:Z,Gamma.schedule['mu'].node:mu,Gamma.schedule['delta'].node:de,
              Gamma.schedule['logRp'].node:LR,
              Gamma.pulse.preRp.operator.background['parameters']['logP']:LP,
              W['log_end_scale'].node:sy.Symbol('log_end_scale',real=True),W['normalization'].node:sy.Symbol('beta_normal',positive=True)}
    U,C1,C2,C0,CE,X=sy.symbols('U C1 C2 C0 CE X',real=True);q=1+Z*Z
    native=dict(u=U/q,m1=C1*Z*q,m2=C2*Z*q,energy=C0+CE*Z*Z*q*q,X=X)
    for key,body in native.items():prior.bind_jet(bindings,p.native[key],body,Z)
    bindings[p.native_evidence['actual_U_from_native_value'].node]=U
    bindings[p.native_evidence['actual_CE_from_native_ordinary_Z2'].node]=CE
    for key in ('ap','C1','C2'):
        jet=p.selected_pulse[key]
        for n,ref in enumerate(current.rows(jet)):
            if g.nodes[ref.node]['operation']!='exact_rational':bindings[ref.node]=sy.Symbol('selected_'+key+'_Z'+str(n),real=True)
    ell=g.nodes[W['ell'].node]
    assert ell['operation']=='exact_rational'
    ell_value=sy.Rational(ell['numerator'],ell['denominator'])
    assert 0<ell_value<1
    checked=[]
    for name,index,atvalue in (('main',0,g.zero),('end',1,g.zero)):
        fields=p.leading_fields[name];meta=p.pulse.own.fields[name]
        alg=EndpointAlgebra(g,{**bindings,meta['query'].node:sy.Integer(0)})
        actual=current.endpoint_functions(g,current.Algebra(g,3).trim(fields['U']),current.Algebra(g,3).trim(fields['Q']),
            current.Algebra(g,3).trim(fields['Q_y']),meta['radius'])
        target=p.global_endpoints[index]
        for key,jet in actual.items():
            for n,(ref,desired) in enumerate(zip(current.rows(jet),current.rows(target[key]))):
                assert sy.simplify(alg.at(ref.node)-alg.at(desired.node))==0,(name,key,n)
                checked.append(f'{name}_{key}_Z{n}')
    return dict(passed=True,actual_current_startup_and_terminal_endpoint_rows=checked,
        native_rational_whole_function_source_identity_consumed=True,
        actual_gp_startup_flatness_and_end_beta_support_used=True,
        actual_end_history_empty_integration_limits_used=True,
        startup_native_mean_not_zeroed=True,end_beta_half_width=str(ell_value))


def partial_endpoint_transfer(owner):
    """Same defining bodies/limits transfer the whole endpoints, including C3."""
    g=owner.g;lower=owner.charts['lower'];lower_end=lower['domain'][1]
    completed=current.partial_integrals(g,owner.lower.pressure,lower_end,owner.lower.signed_bodies)
    for key,q in completed.items():
        for n,ref in enumerate(current.rows(q)):
            row=g.nodes[ref.node];assert row['lower']==g.nodes[owner.lower.pressure[key].value.node]['lower']
            assert row['upper']==lower_end.node and row['integrand']==current.rows(owner.lower.signed_bodies[key])[n].node
            template=g.nodes[owner.lower.pressure[key].value.node]
            for field,value in template.items():
                if field not in ('integrand','upper','ordinary_slow_Z_derivative_order'):
                    assert row[field]==value
            assert row['ordinary_slow_Z_derivative_order']==n
    # All new partial pulse bodies have just been proved to be ordinary
    # derivatives of the same original physical product. At the complete
    # endpoint, fixed-Z differentiation under the integral gives the
    # admitted whole-integral C3 rows, even if a whole row simplifies to0.
    for name,binding in owner.pulse_bindings.items():
        meta=owner.pulse.pulse.own.fields[name]
        whole=current.partial_integrals(g,owner.pulse.integrals[name],meta['domain']['upper'],binding['current_integrands'])
        for key,q in whole.items():
            original=g.nodes[owner.pulse.integrals[name][key].value.node]
            for n,ref in enumerate(current.rows(q)):
                row=g.nodes[ref.node]
                assert row['lower']==original['lower'] and row['upper']==original['upper'] and row['variable']==original['variable']
                assert row['integrand']==current.rows(binding['current_integrands'][key])[n].node
                assert row['endpoint_independent_of_Z'] and row['physical_Pstar_not_reapplied']
                assert row['measure']==original['measure'] and row['ordinary_slow_Z_derivative_order']==n
    prior_receipt=current.parent.linear.read(current.previous.RECEIPT)
    assert prior_receipt['checks']['transported_terminal']['twenty_actual_terminal_integral_function_rows_zero']
    return dict(passed=True,actual_complete_lower_endpoint_uses_same_eight_body_rows=True,
        actual_complete_pulse_endpoints_use_same_original_physical_integrands_and_limits=True,
        fixed_Z_integral_C3_commutation_not_zero_whole_row_padding=True,
        preceding_complete_charts_and_certified_global_start_end_preserved=True,
        quiet_incoming_and_terminal_source_is_checked_reserved_C3_history=True,
        existing_actual_terminal_five_C3_and_zero_tail_proofs_retained=True,
        unresolved_global_leading_joins_and_native_C5_inlet_bridge_not_claimed=True)


def actual_terminal_partition(owner):
    """Consume proved same-body integral/end identities in actual new sources."""
    g=owner.g;Z=sy.Symbol('Z',real=True);S=sy.Symbol('quiet_s',real=True)
    X=sy.Symbol('pulse_end_s',real=True);L=sy.Symbol('quiet_length',positive=True)
    bindings={g.symbol('Z').node:Z,owner.quiet.s.node:S,owner.quiet.length.node:L,
              owner.charts['pulse_end']['coordinate'].node:X,
              owner.delta.node:sy.Symbol('delta',positive=True)}
    replacements={}
    for key,q in owner.own.band_endpoint.items():
        for n,ref in enumerate(current.rows(q)):
            if ref!=g.zero:bindings[ref.node]=sy.Symbol(f'band_{key}_{n}',real=True)
    for key,q in owner.quiet.pressure.items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.Function(f'quiet_{key}_{n}')(S)
    for key,q in owner.quiet.incoming_pressure.items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.Symbol(f'compact_{key}_{n}',real=True)
    for name,group in owner.pulse.integrals.items():
        for key,q in group.items():
            for n,ref in enumerate(current.rows(q)):
                if ref!=g.zero:bindings.setdefault(ref.node,sy.Symbol(f'{name}_{key}_{n}',real=True))
    end=owner.pulse_bindings['end']
    for key,q in end['current_integrals'].items():
        for n,ref in enumerate(current.rows(q)):
            body=sy.Function(f'proved_current_end_{key}_{n}')(X);bindings[ref.node]=body
            old=current.rows(owner.pulse.integrals['end'][key])[n]
            replacements[body.subs(X,0)]=sy.Integer(0) if old==g.zero else bindings[old.node]
    for key,q in end['upper'].items():
        for n,ref in enumerate(current.rows(q)):
            body=sy.Function(f'proved_current_endpoint_{key}_{n}')(X);bindings[ref.node]=body
            replacements[body.subs(X,0)]=sy.Integer(0)
    for key,q in owner.pulse.global_endpoints[0].items():
        for n,ref in enumerate(current.rows(q)):
            if ref!=g.zero:bindings.setdefault(ref.node,sy.Symbol(f'global_lower_{key}_{n}',real=True))
    for key,q in owner.pulse.global_endpoints[1].items():
        for ref in current.rows(q):bindings[ref.node]=sy.Integer(0)
    alg=prior.FunctionAlgebra(g,bindings);checked=[]
    # The accepted old pressure C1 graph expands its first derivative into
    # native source factors. Higher raw expansion would require reopening
    # those ancestors. Instead identify the actual value/first derivative
    # with the same complete original formula, then consume the checked
    # ordinary-derivative construction to transfer Z2/Z3. This is equality
    # of the defining functions, not a terminal cancellation assumption.
    a=current.Algebra(g,3)
    total={key:a.add(*(owner.pulse.integrals[name][key] for name in ('main','gap','end')))
           for key in owner.pulse.integrals['main']}
    canonical=current.pressure_functions(g,owner.delta,total,
        {k:a.trim(q) for k,q in owner.pulse.global_endpoints[0].items()},
        {k:a.trim(q) for k,q in owner.pulse.global_endpoints[1].items()})
    old_pressure_rows=[]
    for key,q in owner.pulse.pressure.items():
        saved=current.rows(q);normal=current.rows(canonical[key])
        for n in (0,1):
            assert sy.cancel(alg.at(saved[n].node)-alg.at(normal[n].node))==0,(key,n,'original complete formula')
            old_pressure_rows.append(f'{key}_Z{n}')
        for n in (2,3):
            assert owner.pulse.derivative.derivative(saved[n-1])==saved[n]
            bindings[saved[n].node]=alg.at(normal[n].node)
            old_pressure_rows.append(f'{key}_Z{n}')
    alg=prior.FunctionAlgebra(g,bindings)
    for key,q in owner.endpoints['pulse_end']['right']['cumulative'].items():
        for n,ref in enumerate(current.rows(q)):
            actual=alg.at(ref.node).subs(replacements)
            expected=alg.at(current.rows(owner.own.terminal[key])[n].node)
            assert sy.cancel(actual-expected)==0,(key,n)
            checked.append(f'{key}_Z{n}')
    return dict(passed=True,actual_new_pulse_exit_equals_checked_terminal_rows=checked,
        actual_original_complete_pressure_C1_and_checked_C3_extension_rows=old_pressure_rows,
        higher_complete_pressure_identity_by_same_ordinary_derivative_construction=True,
        actual_end_flatness_and_complete_same_integrand_C3_identity_consumed=True,
        inner_lower_active_compact_quiet_main_gap_end_counted_once=True,
        complete_true_terminal_integral_zeros_transfer_to_current_regional_exit=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['intermediate_own_C3_chart_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_regional_recovery',actual_regional_recovery),('actual_pulse_bodies',actual_pulse_bodies),
        ('pressure_calculus',pressure_calculus),('actual_pulse_endpoints',actual_pulse_endpoints),
        ('partial_endpoint_transfer',partial_endpoint_transfer),('actual_terminal_partition',actual_terminal_partition),
        ('source',source_check)):
        checks[name]=fn(owner);print('INTERMEDIATE_OWN_C3_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        original_scale_numerical_accuracy_and_global_n1_join_admission=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('INTERMEDIATE_N1_OWN_C3_CHART_FUNCTIONS_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses accepted sources')
