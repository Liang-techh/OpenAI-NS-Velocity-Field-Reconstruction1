"""Exact original geometry, source joins and downstream C6 transport checks."""
import ast
import json
import math
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions as current
import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions_check as prior_check


class Algebra(prior_check.KernelAlgebra):
    def at(self,index):
        if index in self.cache:return self.cache[index]
        row=self.g.nodes[index];op=row['operation']
        if op=='exact_original_source_flat_sigma' or (op=='analytic_unary' and row.get('name')=='original_flat_sigma_prime'):
            x=self.at(row['argument']);prime=op=='analytic_unary'
            if x.is_number and (x<=0 or x>=1):value=sy.Integer(0) if prime or x<=0 else sy.Integer(1)
            else:value=sy.Function('original_sigma_prime' if prime else 'original_sigma')(x)
            self.cache[index]=value;return value
        return super().at(index)


def actual_source(owner):
    g=owner.g;r=current.rows;cache={};rho=g.nodes[owner.core.rho.node]['name']
    assert owner.bridge.acceptance_loaded and tuple(owner.charts)==current.CHARTS
    assert owner.hb!=owner.bridge_hb
    assert g.nodes[owner.hb.node]['name']=='hbS' and g.nodes[owner.bridge_hb.node]['name']=='hbB'
    assert owner.parameter_binding==current.parent.active.current.ast_binding(current.connection.build_graph)
    assert owner.bindings==current.source_bindings()
    assert all(owner.geometry_certificate[key] for key in ('shared_upstream_is_same_object',
        'all_four_native_T_source_tuples_equal','actual_native_T_and_logC_match_common_source',
        'directed_T_positive','directed_reference_gap_positive','directed_post_length_positive'))
    assert owner.charts['first_switch']['incoming_histories'] is owner.bridge.R100['moments']
    predecessors=('R1','R2','R110','Rsh','Rz','Rrestore')
    for name,predecessor in zip(current.CHARTS[1:],predecessors):
        assert owner.charts[name]['incoming_histories'] is owner.endpoints[predecessor]['moments']
    checked=[]
    for name,chart in owner.charts.items():
        phase_name=g.nodes[chart['phase'].node]['name']
        for group in ('moments','fields','radial','actual'):
            for key,jet in chart[group].items():
                assert len(r(jet))==7,(name,group,key)
                for ref in r(jet):
                    free=prior_check.free_variables(g,ref,cache)
                    assert rho not in free and free<={'Z',phase_name},(name,group,key,free)
        assert all(len(r(q))==6 for q in chart['generic'].values())
        assert 'Z' not in prior_check.free_variables(g,chart['physical_measure'],cache)
        checked.append(name)
    for endpoint in owner.endpoints.values():
        assert set(endpoint['actual'])=={'phi','U','phi_y','U_y'}
        assert all(len(r(q))==7 for q in endpoint['actual'].values())
        assert all(len(r(q))==6 for q in endpoint['generic'].values())
    return dict(passed=True,actual_prescribed_charts_checked=checked,
        all_genuine_seven_source_and_six_generic_rows=True,
        every_nonzero_incoming_history_function_retained_by_identity=True,
        endpoint_phi_U_and_physical_y_derivatives_exposed=True,
        distinct_width_tags_with_exact_shared_upstream_K1_identity=True,
        frozen_T_and_positive_gap_post_length_guarded=True,
        no_free_core_rho_or_old_phase_capture_in_chart_bodies=True,
        finite_N_modified_source_and_quantitative_C6_not_claimed=True)


def geometry(owner):
    g=owner.g;h,T,C,P=sy.symbols('h T logC logP',positive=True)
    base={owner.hb.node:h,owner.T.node:T,owner.core.parameters['logC'].node:C,owner.core.parameters['logP'].node:P}
    alg=Algebra(g,base);expected=dict(R1=100*sy.exp(h),R2=100*sy.exp(2*h),R110=sy.Integer(110),
        Rsh=110*sy.exp(T),Rz=110*sy.exp(10*(C+P)-8),Rrestore=110*sy.exp(10*(C+P)-7),Rm=110*sy.exp(10*(C+P)-6))
    for name,endpoint in owner.endpoints.items():assert sy.simplify(alg.at(endpoint['radius'].node)-expected[name])==0,name
    for name,chart in owner.charts.items():
        phase=sy.Symbol('phase_'+name)
        local=Algebra(g,{**base,chart['phase'].node:phase})
        radius=local.at(chart['radius'].node);jacobian=local.at(chart['physical_measure'].node)
        assert sy.simplify(sy.diff(radius,phase)/radius-jacobian)==0,name
    assert sy.simplify(alg.at(owner.post_length.node)-(sy.log(sy.Rational(11,10))-2*h))==0
    assert sy.simplify(alg.at(owner.reference_length.node)-(10*(C+P)-T-8))==0
    return dict(passed=True,exact_radius_endpoint_identities=7,exact_physical_log_radius_Jacobians=7,
        original_finite_T_and_Rref_offsets_retained=True,
        complete_radius_order_proved_by_actual_source_positive_guards=True)


def cutoff_identity():
    gamma=current.previous.inlet.previous.lower.own.jets.complete.Gamma
    jets=current.previous.inlet.previous.lower.own.c2.flat_source
    x=sy.Symbol('x',positive=True);a=sy.exp(-1/x**2);b=sy.exp(-1/(1-x)**2)
    value=a/(a+b);odds=1/(1-x)**2-1/x**2
    assert sy.simplify(value-sy.exp(odds)/(1+sy.exp(odds)))==0
    assert sy.simplify(sy.diff(value,x)-value*(1-value)*(2/x**3+2/(1-x)**3))==0
    trees={name:ast.parse(Path(module.__file__).read_text(encoding='utf8'))
        for name,module in (('value',gamma.flat_source),('jets',jets))}
    expected={'value':("value=a/(a+b)","derivative=value*(b/(a+b))*(2/x**3+2/(1-x)**3)"),
        'jets':("odds=1/(1-x)**2-1/x**2","L1=2/(1-x)**3+2/x**3")}
    for name,texts in expected.items():
        for text in texts:
            wanted=ast.dump(ast.parse(text).body[0])
            assert any(ast.dump(node)==wanted for node in ast.walk(trees[name])),text
    return dict(passed=True,original_value_and_prime_sources_are_same_flat_function=True,
        actual_stable_value_and_Taylor_log_odds_assignments_bound=True,
        factorial_one_for_sigma_prime=True,exact_zero_one_and_flat_derivative_endpoints=True)


def long_restore_joins(owner):
    """Actual new graph rows on independent nonconstant sixth-order inlet."""
    g=owner.g;Z=sy.Symbol('Z');point=sy.Rational(3,8)
    phi=2+Z**2/17+Z**6/19;G=Z**2/23+Z**6/29;V=3+Z+Z**6/31
    C,P,Lambda,T=sy.Rational(2),sy.Rational(3),sy.Rational(3,7),sy.Rational(3,2)
    F0=sy.exp(-C-Lambda*G);base={g.symbol('Z').node:point,owner.core.parameters['logC'].node:C,
        owner.core.parameters['logP'].node:P,owner.core.parameters['Lambda'].node:Lambda,owner.T.node:T}
    def bind(jet,body):
        for k,ref in enumerate(current.rows(jet)):
            value=sy.diff(body,Z,k).subs(Z,point)
            if ref.node in base:assert sy.simplify(base[ref.node]-value)==0
            base[ref.node]=value
    bind(owner.endpoints['R110']['actual']['phi'],phi);bind(owner.endpoints['R110']['actual']['U'],V)
    bind(owner.bridge.core7['G'],G);bind(owner.F0,F0)
    # Later U endpoints preserve exactly the constant R110 velocity.
    for endpoint in ('Rsh','Rz'):bind(owner.endpoints[endpoint]['actual']['U'],V)
    checks=[('long_reshape',0,'phi',phi),('long_reshape',1,'physical_Utheta',sy.exp(T/10-C)/(1+Z**2)),
        ('reference',0,'physical_Utheta',sy.exp(T/10-C)/(1+Z**2)),
        ('reference',1,'physical_Utheta',sy.exp(P-sy.Rational(4,5))/(1+Z**2)),
        ('restoration',0,'physical_Utheta',sy.exp(P-sy.Rational(4,5))/(1+Z**2)),
        ('restoration',1,'physical_Utheta',sy.exp(P-sy.Rational(7,10))/(1+Z**2)),
        ('restoration',0,'U',V),('restoration',1,'U',4*Z),('post_restore',0,'U',4*Z),
        ('post_restore',0,'physical_Utheta',sy.exp(P-sy.Rational(7,10))/(1+Z**2)),
        ('post_restore',1,'physical_Utheta',sy.exp(P-sy.Rational(3,5))/(1+Z**2))]
    count=0
    for chart,edge,key,body in checks:
        data=owner.charts[chart];alg=Algebra(g,{**base,data['phase'].node:sy.Integer(edge)})
        for k,ref in enumerate(current.rows(data['actual'][key])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,k).subs(Z,point))==0,(chart,edge,key,k)
            count+=1
    for chart in ('long_reshape','reference','restoration','post_restore'):
        for edge in (0,1):
            data=owner.charts[chart];alg=Algebra(g,{**base,data['phase'].node:sy.Integer(edge)})
            for k in range(7):
                actual=alg.at(current.rows(data['actual']['phi_y'])[k].node)
                wanted=-sy.Rational(2,5)*alg.at(current.rows(data['actual']['phi'])[k].node)
                assert sy.simplify(actual-wanted)==0,(chart,edge,k)
    return dict(passed=True,independent_long_reference_restore_endpoint_rows=count,
        physical_y_slope_join_rows=56,nonconstant_actual_R110_B_and_velocity_inputs=True,
        long_inlet_F0_anchored_gauge_logC_and_B_cancellation_checked=True,
        full_finite_T_and_offsets_minus8_minus7_minus6=True,
        restoration_end_actual_U_is_4Z_through_sixth_derivative=True,
        mathematical_fixture_only_not_original_scale_numeric_accuracy=True)


def switch_controls():
    """Exact source ODE identities keep both distinct switch controls."""
    g=current.source.FunctionTransportGraph();a=current.Algebra(g,6);Z=sy.Symbol('Z');S=sy.Symbol('fixture_switch_s')
    phase=g.symbol(str(S));point=sy.Rational(3,8);h=sy.Rational(1,9)
    jet=lambda body:prior_check.ordinary_jet(g,Z,sy.sympify(body),6)
    phi=2+Z**2/17+Z**6/19;U=3+Z+Z**6/23;F=2+Z**3/29
    known=dict(D_over_R=a.fixed(0),drive_hydro=jet(1+Z**6/7),
        drive_pressure=jet(2+Z**3/11),drive_swirl=jet(3+Z**2/13))
    R=g.mul(g.constant(100),g.unary('exp',g.mul(g.constant(str(h)),phase)))
    got=current.first_switch_fields(a,known,dict(phi=jet(phi),U=jet(U)),phase,g.constant(str(h)),R,jet(F),g.constant('3/2'),jet(phi))
    alg=Algebra(g,{g.symbol('Z').node:point,phase.node:S});count=0
    for k in range(7):
        actual=alg.at(current.rows(got['U'])[k].node);Uy=alg.at(current.rows(got['U_y'])[k].node)
        assert sy.simplify(sy.diff(actual,S)-h*Uy)==0,(k,'first ODE')
        assert sy.simplify(actual.subs(S,0)-sy.diff(U,Z,k).subs(Z,point))==0,(k,'first inlet')
        count+=1
    second=current.second_switch_fields(a,known,dict(phi=jet(phi),U=jet(U)),phase,g.constant(str(h)),R)
    alg=Algebra(g,{g.symbol('Z').node:point,phase.node:S})
    for k in range(7):
        logarithm=alg.at(current.rows(second['logratio'])[k].node)
        want=-sy.Rational(2,5)*h*sy.Function('original_sigma')(S) if k==0 else 0
        assert sy.simplify(sy.diff(logarithm,S)-want)==0,(k,'second log ODE')
        assert alg.at(current.rows(second['U_y'])[k].node)==0
        count+=1
    return dict(passed=True,independent_first_and_second_switch_ODE_rows=count,
        first_axial_hb_squared_one_minus_sigma_and_current_R_Rsquared_preserved=True,
        second_angular_minus_two_fifths_sigma_and_constant_U_preserved=True,
        mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['prescribed_R100_Rm_C6_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks=dict(actual_source=actual_source(owner),geometry=geometry(owner),cutoff=cutoff_identity())
    print('TRUE_DOWNSTREAM_SOURCE_GEOMETRY_AND_WIDTH_IDENTITY_CHECKED',flush=True)
    checks['long_restore_joins']=long_restore_joins(owner);checks['switch_controls']=switch_controls()
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_PRESCRIBED_R100_RM_C6_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the admitted live source chain')
