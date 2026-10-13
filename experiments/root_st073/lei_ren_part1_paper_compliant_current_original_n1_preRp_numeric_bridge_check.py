"""Actual function replay, reduced-S coefficient recipes and pre-Rp bounds."""
import ast
import inspect
import json
import textwrap
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_preRp_numeric_bridge as current


def differentiated_coefficient_recipes():
    """Compare all 7/9 production cap coefficients with exact reduced-S terms.

    IntervalTaylor differentiates these coefficients before taking absolute
    values. The factors 2 cover both finite endpoints; integral factors 1
    cover the full positive physical radius measure. This check concerns
    mathematical recipes, not numerical fixture values.
    """
    z,de=sy.symbols('Z delta',real=True)
    L,d=1-de*z*z,1-z*z
    env=dict(z=z,delta=de,a=1-de,d=d,invL=1/L,half=sy.Rational(1,2))
    tree=ast.parse(textwrap.dedent(inspect.getsource(current.caps.active_pressure_C1)))
    observed={}
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('S0','S1') for t in node.targets):
            name=node.targets[0].id;rows=[]
            for term in node.value.args[0].elts:
                assert isinstance(term.func,ast.Name) and term.func.id=='term'
                coefficient=eval(compile(ast.Expression(term.args[0]),'<original-preRp-coefficient>','eval'),{'__builtins__':{}},env)
                pair=term.args[1];assert isinstance(pair,ast.Subscript) and pair.value.id=='pairs'
                field=ast.literal_eval(pair.slice)
                factor=ast.literal_eval(term.args[2]) if len(term.args)>2 else 1
                rows.append((sy.sympify(coefficient),field,factor))
            observed[name]=rows
    expected=dict(S0=[(1/L,'q',2),((1-de)*z/(2*L),'qZ',1),(sy.Rational(1,2),'q2',2),
        (d/L,'UqZ',1),(-2*z*z/L**2,'U2',2),(-2*de*z*z/L**2,'U2',1),(2*z*d/L**2,'UUZ',1)],
        S1=[(1/L,'q',2),(-1/L,'q',1),((1-de)*z/(2*L),'qZ',1),(sy.Rational(1,2),'q2',2),
        (-sy.Rational(1,2),'q2',1),(d/L,'UqZ',1),(-2*z*z/L**2,'U2',2),
        (2*(1-de)*z*z/L**2,'U2',1),(2*z*d/L**2,'UUZ',1)])
    assert set(observed)==set(expected)
    for name in expected:
        assert len(observed[name])==len(expected[name])
        for got,want in zip(observed[name],expected[name]):
            assert got[1:]==want[1:] and sy.cancel(got[0]-want[0])==0
            assert sy.cancel(sy.diff(got[0],z)-sy.diff(want[0],z))==0
    return dict(passed=True,exact_S0_coefficients_and_ordinary_Z_derivatives=7,
        exact_S1_coefficients_and_ordinary_Z_derivatives=9,
        finite_endpoint_two_sided_radius_majorants_retained=True,
        all_pressure_coefficient_derivatives_retained=True,
        generic_q_U_product_rules_and_native_derivative_bounds_inherited=True)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported true pre-Rp numerical query admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();owner.assert_graph();proof=differentiated_coefficient_recipes()
    assert owner.source_identity['passed'] and raw['source_identity']==owner.source_identity
    assert len(packets)==len(raw['samples'])==3
    assert owner.bindings==raw['true_function_bindings']
    checks=[]
    for name,packet in packets.items():
        view=owner.report(packet);assert current.op.snapshot(view)==current.op.snapshot(raw['samples'][name])
        assert all(view[k] is False for k in current.OPEN)
        beta=view['original_positive_cover_per_nonzero_piece'];assert current.ep(beta)[0]>0
        assert view['ordinary_Z0_Z1_enclosures_not_defining_function_values']
        assert view['quiet_m13_positive_m15_negative_and_Z1_exact_zero']
        assert view['no_original_radius_or_hyper_exponential_cap_materialized']
        for cell in view['cell_piece_enclosures']:
            for key,pieces in cell['piece_enclosures'].items():
                for label in ('active','compact'):
                    for value in pieces[label]:
                        lo,hi=current.ep(value);assert lo<0<hi
                value,derivative=pieces['quiet']
                assert current.ep(derivative)==(0,0)
                assert current.ep(value)[0]>=0 if key=='m13' else current.ep(value)[1]<=0
        for key,pair in view['normalized_true_preRp_C1_enclosures'].items():
            assert len(pair)==2
            # A directed sequence of sums can exceed the independently
            # rounded 3*b/2*b expression by an ulp. Retain this rounding
            # width in the reported interval and use a positive cover.
            for value,maximum in zip(pair,(4*beta,3*beta)):
                assert current.ep(abs(value))[1]<=current.ep(maximum)[1]
        checks.append(dict(name=name,active_cells=view['active_cells'],actual_query_executed=True,
            true_pressure_value_Z1_enclosure_rows=4,quiet_Z1_exact_zero=True))
    assert len(owner.report(packets['boundary'])['active_cells'])==2
    assert len(owner.report(packets['whole_Z'])['active_cells'])==4
    assert rejected(lambda:owner.evaluate('1.01'))
    assert rejected(lambda:owner.report(current.lower.OriginalTrueLowerNumericPacket()))
    original=owner._bindings
    try:
        owner._bindings={};assert rejected(owner.assert_graph)
    finally:owner._bindings=original
    owner.assert_graph()
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=owner.source_identity,
        differentiated_coefficient_proof=proof,query_results=checks,
        actual_source_function_replay_binds_numeric_magnitude_providers=True,
        signed_field_centers_and_uniform_relative_accuracy_not_claimed=True,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.op.snapshot(receipt),indent=2)+'\n').encode('utf8'))
    print('TRUE_PRERP_NUMERIC_C1_BRIDGE_CHECKED',len(checks),flush=True)
    return receipt


if __name__=='__main__':print(json.dumps(differentiated_coefficient_recipes(),indent=2))
