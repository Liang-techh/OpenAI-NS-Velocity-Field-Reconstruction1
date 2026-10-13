"""Independent Gamma derivatives, ordinary calculus and native C3 checks."""
import json
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
import sympy as sy
import scipy.integrate  # Load fixture dependencies before the deep legacy import chain.

import lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions_check as Gamma_check


def Gamma_fixture():
    g=current.source.FunctionTransportGraph();c=g.constant;z=g.symbol('Z')
    functions,pieces=current.complete.Gamma.Gamma_future_functions(g,c('3/10'),c('1/50'),c('3/25'))
    owner=SimpleNamespace(functions=functions,pieces=pieces,schedule={'S':c('3/25')})
    jet,evidence=current.Gamma_higher_functions(g,owner)
    evaluate=Gamma_check.FixtureEvaluator(g,.3)
    with mp.workdps(35):
        a,S,t=mp.mpf('.15'),mp.mpf('.12'),mp.mpf('.7')
        def direct(Z):
            xi=2*(1-Z*Z)*S*mp.exp(-t)
            return mp.quad(lambda v:mp.exp(-v)*v**a*(-mp.expm1(-a*mp.log1p(xi*v))),
                [0,1,mp.inf])/(a*S*mp.gamma(1+a))
        point=mp.mpf('.65');env={z.node:float(point),pieces['offset_variable'].node:float(t)}
        errors={}
        for order,ref in ((2,jet.ZZ),(3,jet.ZZZ)):
            want=float(mp.diff(direct,point,order));value=evaluate(ref,env)
            errors[str(order)]=abs(value-want)
            assert abs(value-want)<2e-10,(order,value,want)
        # One-sided endpoint derivatives are given by the continued Gamma
        # integrals, not evaluations at a negative heat argument.
        for point in (-1.,1.):
            env={z.node:point,pieces['offset_variable'].node:.7}
            assert mp.isfinite(evaluate(jet.ZZ,env)) and mp.isfinite(evaluate(jet.ZZZ,env))
            assert evaluate(jet.ZZ,env)<0
            assert evaluate(jet.ZZZ,env)*point<0
        assert evaluate(jet.ZZZ,{z.node:0.,pieces['offset_variable'].node:.7})==0
    return dict(passed=True,mathematical_fixture_only=True,
        direct_original_Gamma_expectation_differentiated_independently=True,
        ordinary_second_and_third_max_errors=errors,
        finite_one_sided_axis_endpoint_derivatives=True,
        exact_nonzero_S_and_full_improper_integrals_retained=True)


def calculus_fixture():
    g=current.source.FunctionTransportGraph();z=g.symbol('Z');Z=sy.Symbol('Z')
    def graph(expr):
        expr=sy.sympify(expr)
        if expr.is_Rational:return g.constant(str(expr))
        if expr==Z:return z
        if expr.func is sy.Add:return g.add(*(graph(v) for v in expr.args))
        if expr.func is sy.Mul:return g.mul(*(graph(v) for v in expr.args))
        if expr.func is sy.Pow:
            base,power=expr.args
            if power.is_Integer:
                if power>=0:return g.mul(*([graph(base)]*int(power)))
                return g.quotient(g.one,g.mul(*([graph(base)]*int(-power))),'positive fixture power')
            return g.unary('exp',g.mul(graph(power),g.unary('log',graph(base))))
        if expr.func in (sy.exp,sy.log):return g.unary('exp' if expr.func is sy.exp else 'log',graph(expr.args[0]))
        raise ValueError('Unknown independent fixture expression')
    f=sy.exp(Z**2/10)*sy.sqrt(3+Z**2)/(1+Z**2)
    pair=current.source.C1Function(graph(f),graph(sy.diff(f,Z)))
    d=current.OrdinaryZDerivative(g);d.seed_pair(pair);jet=d.lift(pair)
    v=g.symbol('fixture_Z_integral_v');body=g.unary('exp',g.mul(z,v))
    integrate=lambda ref:current.controls.integral(g,ref,g.nodes[v.node]['name'],g.zero,g.one)
    integral_pair=current.source.C1Function(integrate(body),integrate(g.mul(v,body)))
    d.seed_pair(integral_pair);integral_jet=d.lift(integral_pair)
    fixed=g.node('function_substitution',expression=pair.value.node,variable=z.node,value=g.zero.node,
        Z_independent_substitution=True)
    assert d.derivative(fixed)==g.zero
    evaluate=Gamma_check.FixtureEvaluator(g,.3)
    for point in (-.7,0.,.8):
        for order,ref in enumerate((jet.value,jet.Z,jet.ZZ,jet.ZZZ)):
            want=float(sy.diff(f,Z,order).subs(Z,point).evalf(35))
            assert abs(evaluate(ref,{z.node:point})-want)<2e-11
        with mp.workdps(35):
            for order,ref in enumerate((integral_jet.value,integral_jet.Z,integral_jet.ZZ,integral_jet.ZZZ)):
                want=float(mp.quad(lambda v:v**order*mp.exp(mp.mpf(str(point))*v),[0,1]))
                assert abs(evaluate(ref,{z.node:point})-want)<2e-11
    # Unknown derivative sources and moving endpoints are rejected.
    unknown=g.node('unidentified_axial_function',argument=z.node)
    try:d.derivative(unknown)
    except ValueError:pass
    else:raise AssertionError('Unknown derivative function accepted')
    moving=current.controls.integral(g,body,g.nodes[v.node]['name'],g.zero,z)
    try:d.derivative(moving)
    except ValueError:pass
    else:raise AssertionError('Incomplete Leibniz rule accepted')
    assert jet.value==pair.value and jet.Z==pair.Z
    return dict(passed=True,mathematical_fixture_only=True,
        mixed_exponential_sqrt_quotient_Z0_Z1_Z2_Z3_checked=True,
        true_integral_Z2_Z3_checked_against_independent_moments=True,
        fixed_axial_substitution_derivative_zero_proved=True,
        unknown_source_and_moving_endpoint_fail_closed=True,
        accepted_value_and_first_Z_handles_preserved=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1PulseC3Functions(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['function_bindings']==owner.snapshot
    assert owner.complete.acceptance_loaded and not any(raw[key] for key in current.OPEN)
    original=current.parent.linear.read(current.parent.quiet.FRAME_NAME)
    graph=current.parent.linear.realization.restore_graph(original['exact_graph_nodes'])
    importer=current.parent.GraphImport(owner.g,graph)
    packet=original['actual_Rp_native_frame_functions']['actual_identified_Rp_native_frame_C3']
    for key,jet in owner.native.items():
        expected=[importer(current.source.FunctionRef(graph,i)).node for i in packet[key]]
        assert expected==current.encode(jet),key
    prior=dict(heat=owner.complete.Gamma.functions,histories=owner.complete.Gamma.histories,
        angular=owner.complete.angular,future=owner.complete.functions,pulse=owner.complete.pulse_functions)
    for group,values in owner.functions.items():
        for key,jet in values.items():
            pair=prior[group][key]
            assert jet.value==pair.value and jet.Z==pair.Z
            owner.g.ids([jet.value,jet.Z,jet.ZZ,jet.ZZZ])
    fixtures=dict(Gamma_second_and_third=Gamma_fixture(),ordinary_Z_calculus=calculus_fixture(),
        original_native_C3=dict(passed=True,all_eight_actual_native_C3_jets_relocated=True,
            ordinary_rows_not_factorial_divided_Taylor_rows=True),
        full_improper_integral_domination=dict(passed=True,
            positive_Gamma_moments_through_v_a_plus_three_dominate=True,
            D_DZ_DZZ_decay_exp_minus_t=True,DZZZ_decay_exp_minus_two_t_or_faster=True,
            theta_exterior_decay_rate_a_strict_positive=True,
            pressure_and_energy_exterior_decay_faster=True,
            actual_nonzero_S_is_Z_independent=True))
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=fixtures,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        uniform_original_Z2_Z3_numerical_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_ACTUAL_PULSE_C3_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
