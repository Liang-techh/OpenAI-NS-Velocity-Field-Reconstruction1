"""Actual source replay and independent second-derivative majorant checks.

Finite scalar arithmetic checks validate the bound calculus only. The four
current Z-cell replay separately binds the astronomical source functions.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_target_ranges as current

phase=current.phase


def independent_majorant_calculus():
    z=s.Symbol('Z');eta=s.Symbol('eta',positive=True)
    a=s.Function('a')(z);delta=s.Function('Delta')(z)
    gamma=2*eta-delta;root=s.sqrt(gamma/(2*a))
    L=s.log(gamma)/2-s.log(a)/2
    expected=-s.diff(delta,z,2)/(2*gamma)-s.diff(delta,z)**2/(2*gamma**2) \
        -s.diff(a,z,2)/(2*a)+s.diff(a,z)**2/(2*a**2)
    assert s.simplify(s.diff(L,z,2)-expected)==0
    assert s.simplify(s.diff(root,z,2)-root*(s.diff(L,z,2)+s.diff(L,z)**2))==0
    # Bounds.div accepts log(lower), hence 2*log(gamma_min) gives
    # gamma_min**2, with no extra half factor in either squared term.
    m=s.Symbol('m',positive=True)
    for n,C,maximum in ((1,4,32),(2,28,1792)):
        envelope=C*s.exp(4-m**-2)*m**(-3*n)
        assert s.simplify(s.diff(s.log(envelope),m)-(2-3*n*m*m)/m**3)==0
        assert 2-s.Rational(3*n,4)>0
        assert envelope.subs(m,s.Rational(1,2))==maximum
    assert phase.flat_source.SIGMA_CONSTANTS[1:3]==(4,28)
    u=s.Function('u')(z)
    for f,fu,fuu in (
        (u/s.sqrt(1+u*u),(1+u*u)**s.Rational(-3,2),-3*u*(1+u*u)**s.Rational(-5,2)),
        ((1+u*u)**s.Rational(-1,2),-u*(1+u*u)**s.Rational(-3,2),
            (2*u*u-1)*(1+u*u)**s.Rational(-5,2))):
        assert s.simplify(s.diff(f,z,2)-fu*s.diff(u,z,2)-fuu*s.diff(u,z)**2)==0
    # The implicit inverse follows twice differentiating Phi(psi(Z),Z)=phi.
    x=s.Function('psi')(z);X,Y=s.symbols('X Y')
    c00,c10,c01,c20,c11,c02=s.symbols('c00 c10 c01 c20 c11 c02')
    # A completely arbitrary local two-jet checks all independent chain terms.
    f=c00+c10*X+c01*Y+c20*X**2/2+c11*X*Y+c02*Y**2/2
    composed=s.diff(f.subs({X:x,Y:z}),z,2)
    formula=s.diff(f,Y,2).subs({X:x,Y:z})+2*s.diff(f,X,Y).subs({X:x,Y:z})*s.diff(x,z) \
        +s.diff(f,X,2).subs({X:x,Y:z})*s.diff(x,z)**2+s.diff(f,X).subs({X:x,Y:z})*s.diff(x,z,2)
    assert s.simplify(composed-formula)==0
    n=s.Function('numerator')(z);d=s.Function('positive_denominator')(z);q=n/d
    assert s.simplify(d*s.diff(q,z,2)+2*s.diff(d,z)*s.diff(q,z)+s.diff(d,z,2)*q-s.diff(n,z,2))==0
    return dict(root_second_log_derivative_full_squared_terms=True,
        denominator_power_convention='div(v, k*log_lower) divides by lower**k',
        sigma_global_bounds_derived_from_original_constants=[32,1792],
        signed_u_chain_rules_and_inverse_chain_rule_checked=True,
        positive_quotient_second_row_binomial_factor_two_checked=True,
        scope='independent calculus identities; not current numerical point values')


def directed_arithmetic_check():
    c=MPIntervalContext();c.dps=90;bd=current.Bounds(c)
    def upper(q):return mp.mpf(0) if q.log is None else current.ep(c.exp(q.log))[1]
    loglower=c.ln(c.mpf(3))
    cap=bd.div(bd.constant(18),2*loglower)
    assert upper(cap)>=2 and upper(cap)<mp.mpf('2.0000001')
    a=current.JetBound(*(bd.constant(v) for v in (2,3,5)))
    b=current.JetBound(*(bd.constant(v) for v in (7,11,13)))
    product=bd.product(a,b)
    for q,want in zip((product.value,product.Z,product.ZZ),(14,43,127)):
        assert upper(q)>=want
    n=current.JetBound(*(bd.constant(v) for v in (19,23,29)))
    q=bd.quotient(n,b,c.ln(c.mpf(7)))
    q0=mp.mpf(19)/7;q1=(23+11*q0)/7;q2=(29+22*q1+13*q0)/7
    for got,want in zip((q.value,q.Z,q.ZZ),(q0,q1,q2)):assert upper(got)>=want
    assert bd.div(bd.zero,2*loglower).log is None
    return dict(passed=True,log_lower_square_division_checked=True,
        product_and_quotient_second_row_directed_magnitude_bounds_checked=True,
        scalar_fixture_is_not_current_source=True)


def encoded(value):return json.loads(json.dumps(phase.encoded(current.record(value))))


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C2_target_ranges_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2TargetRanges(require_checked=False)
        assert not field.acceptance_loaded and field.target.acceptance_loaded
        assert raw['source_family']==field.identity
        parent=json.loads((current.HERE/current.target.RECEIPT).read_bytes())
        assert parent['all_passed'] and all(parent[k] for k in current.target.GATES)
        assert raw['actual_current_regular_majorant_theorem']==encoded(field.proof)
        calculus=independent_majorant_calculus();arithmetic=directed_arithmetic_check()
        bindings=dict(primitive_bounds=current.current.ast_binding(current.primitive_bounds),
            density_bounds=current.current.ast_binding(current.density_bounds),
            actual_source_projection=current.current.ast_binding(phase.CurrentC2PhasePrimitives.source_packet),
            original_log_bound_arithmetic=current.current.ast_binding(current.LogUpper),
            original_phase_density_transport_source=current.target.RECEIPT,
            original_phase_density_transport_sha256=current.sha(current.target.RECEIPT))
        assert raw['actual_C2_function_source_bindings']==bindings
        counts=dict(actual_Z_cells=0,actual_chart_transports=0,actual_native_cells=0,
            own_rate_chart_steps=0,nonzero_quiet_second_row_memories=0,actual_second_target_bounds=0,
            exact_flat_cells=0,body_cells=0,active_flat_union_cells=0)
        for ends,saved in zip(current.outer.CELLS,raw['actual_four_Z_C2_target_transports']):
            replay=encoded(field.transport(ends));assert replay==saved,ends
            assert tuple(saved['exact_Z_cell'])==tuple(ends)
            assert saved['actual_five_target_C2_function_handles']==encoded(field.target.target_functions())
            for chart,window in zip(phase.CHARTS,saved['actual_17_chart_C2_transports']):
                assert window['chart']==chart and window['second_row_predecessor_never_reset']
                assert window['global_phase_cover']==[0,1]
                partitions=current.current.native_partitions()[chart]
                assert len(window['cells'])==len(partitions)-1
                for (left,right),cell in zip(zip(partitions,partitions[1:]),window['cells']):
                    assert cell['exact_left']==encoded(left) and cell['exact_right']==encoded(right)
                    assert cell['physical_Jacobian_applied_once']
                    query=field.cell(ends,chart,left,right)
                    proof=query['primitive']['proof'];assert proof['exact_same_source_P0']
                    assert proof['source_context_basis_and_ledger_retained']
                    assert proof['actual_C2_function_handles']==phase.encoded(field.phase.functions[chart])
                    branch=proof['branch'];counts[dict(exact_original_flat='exact_flat_cells',
                        body_sigma_one='body_cells',active_flat_union='active_flat_union_cells')[branch]]+=1
                    for key,rate in current.current.RATES.items():
                        weights=field.weights(ends,chart,left,right,rate)
                        for name in ('original_mass','cell_decay','suffix_decay'):
                            assert current.ep(weights[name].coefficient)[0]>=0
                        if branch=='exact_original_flat':
                            assert all(q.log is None for q in (
                                query['density']['N_scaled'][key].value,
                                query['density']['N_scaled'][key].Z,
                                query['density']['N_scaled'][key].ZZ))
                    counts['actual_native_cells']+=1
                for key in current.current.RATES:
                    if chart=='O3_power':
                        assert not window['actual_N_scaled_incoming_C2'][key][2]['exact_zero']
                        assert window['actual_N_scaled_local_C2'][key][2]['exact_zero']
                        assert not window['actual_N_scaled_outgoing_C2'][key][2]['exact_zero']
                        counts['nonzero_quiet_second_row_memories']+=1
                    counts['own_rate_chart_steps']+=1
                counts['actual_chart_transports']+=1
            for row in saved['actual_N_scaled_five_target_C2_ranges'].values():
                assert len(row)==3 and not row[2]['exact_zero']
                assert row[2]['source_value_or_exponential_not_materialized']
                counts['actual_second_target_bounds']+=1
            assert saved['uniform_for_all_integer_N_ge_N0'] and saved['range_caps_not_point_function_values']
            assert saved['original_independent_P0_in_source_and_not_added_twice']
            counts['actual_Z_cells']+=1
        assert len(raw['actual_four_Z_C2_target_transports'])==len(current.outer.CELLS)==4
        assert counts['actual_chart_transports']==68 and counts['actual_native_cells']==228
        assert counts['own_rate_chart_steps']==340 and counts['nonzero_quiet_second_row_memories']==20
        assert counts['actual_second_target_bounds']==20
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            quantitative_current_C2_target_ranges_installed=True,
            independent_majorant_calculus=calculus,independent_directed_arithmetic=arithmetic,
            actual_replay_counts=counts,actual_C2_repaired_limit_controls_installed=False,
            current_numeric_point_field_oracle_installed=False,range_bounds_are_not_signed_function_values=True,
            **dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
                Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(phase.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual four-Z 17-chart C2 magnitude ranges and five-target bounds checked',flush=True)
    return result


if __name__=='__main__':run()
