"""Independent cubic magnitude calculus and actual four-cell source replay."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_target_ranges as current

phase=current.phase


def encoded(value):return json.loads(json.dumps(phase.encoded(current.record(value))))


def independent_majorant_calculus():
    z=s.Symbol('Z');eta=s.Symbol('eta',positive=True);a=s.Function('a')(z);D=s.Function('Delta')(z)
    gamma=2*eta-D;L=s.log(gamma)/2-s.log(a)/2;root=s.sqrt(gamma/(2*a))
    want=-s.diff(D,z,3)/(2*gamma)-3*s.diff(D,z)*s.diff(D,z,2)/(2*gamma**2)-s.diff(D,z)**3/gamma**3 \
        -s.diff(a,z,3)/(2*a)+3*s.diff(a,z)*s.diff(a,z,2)/(2*a**2)-s.diff(a,z)**3/a**3
    assert s.simplify(s.diff(L,z,3)-want)==0
    assert s.simplify(s.diff(root,z,3)-root*(s.diff(L,z,3)+3*s.diff(L,z)*s.diff(L,z,2)+s.diff(L,z)**3))==0
    m=s.Symbol('m',positive=True);envelope=256*s.exp(4-m**-2)*m**-9;peak=s.sqrt(2)/3
    assert s.simplify(s.diff(s.log(envelope),m)-(2-9*m*m)/m**3)==0
    assert s.simplify(envelope.subs(m,peak)-256*s.exp(-s.Rational(1,2))*(3/s.sqrt(2))**9)==0
    assert peak<s.Rational(1,2) and phase.flat_source.SIGMA_CONSTANTS[3]==256
    v=s.Symbol('u',real=True);u=s.Function('u')(z)
    for f in (v/s.sqrt(1+v*v),1/s.sqrt(1+v*v)):
        assert s.simplify(s.diff(f.subs(v,u),z,3)-s.diff(f,v,3).subs(v,u)*s.diff(u,z)**3
            -3*s.diff(f,v,2).subs(v,u)*s.diff(u,z)*s.diff(u,z,2)-s.diff(f,v).subs(v,u)*s.diff(u,z,3))==0
    # Uniform r''' bound follows 12u^2+3 <=15(1+u^2).
    assert s.expand(15*(1+v*v)-(12*v*v+3))==3*v*v+12
    # For h''', |u|/sqrt(1+u^2)<=1 bounds each 9u and 6u^3 term by 9 and 6.
    n=s.Function('numerator')(z);d=s.Function('denominator')(z);q=n/d
    assert s.simplify(d*s.diff(q,z,3)+3*s.diff(d,z)*s.diff(q,z,2)
        +3*s.diff(d,z,2)*s.diff(q,z)+s.diff(d,z,3)*q-s.diff(n,z,3))==0
    psi=s.Symbol('psi');r=s.Function('r')(z);den=1-2*r*s.cos(psi)+r*r
    assert s.simplify(s.diff(den,z,3)-6*s.diff(r,z)*s.diff(r,z,2)-2*(r-s.cos(psi))*s.diff(r,z,3))==0
    w=(s.cos(psi)-r)/den;dp=s.diff(den,psi)
    assert s.simplify(s.diff(w,psi)-(-s.sin(psi)-w*dp)/den)==0
    assert s.simplify(s.diff(w,psi,2)-(-s.cos(psi)-2*s.diff(w,psi)*dp-w*s.diff(den,psi,2))/den)==0
    return dict(root_third_log_derivative_all_cubic_terms_checked=True,
        denominator_power_convention='k*log(lower) means lower**k',
        sigma_third_source_envelope_global_peak_checked=True,signed_u_third_chain_rule_checked=True,
        uniform_signed_u_third_bound=15,Poisson_third_and_angle_cross_partials_checked=True,
        third_positive_quotient_binomial_coefficients_checked=True)


def directed_arithmetic_check():
    c=MPIntervalContext();c.dps=120;bd=current.Bounds(c)
    def upper(q):return mp.mpf(0) if q.log is None else current.ep(c.exp(q.log))[1]
    with mp.workdps(180):
        loglower=c.ln(c.mpf(3));got=bd.div(bd.constant(54),3*loglower)
        assert upper(got)>=2 and upper(got)<mp.mpf('2.0000001')
        a=current.JetBound(*(bd.constant(v) for v in (2,3,5,17)))
        b=current.JetBound(*(bd.constant(v) for v in (7,11,13,19)))
        product=bd.product(a,b)
        for got,want in zip(current.rows(product),(14,43,127,439)):assert upper(got)>=want
        n=current.JetBound(*(bd.constant(v) for v in (23,29,31,37)))
        q=bd.quotient(n,b,c.ln(c.mpf(7)))
        q0=mp.mpf(23)/7;q1=(29+11*q0)/7;q2=(31+22*q1+13*q0)/7;q3=(37+33*q2+39*q1+19*q0)/7
        for got,want in zip(current.rows(q),(q0,q1,q2,q3)):assert upper(got)>=want
        exact=256*mp.exp(-mp.mpf('.5'))*(3/mp.sqrt(2))**9
        lo,hi=current.ep(current.sigma_third_cap(c));assert hi>=exact and lo<=exact
        assert bd.div(bd.zero,3*loglower).log is None
    return dict(passed=True,log_lower_cubic_division_checked=True,
        directed_cubic_product_and_quotient_checked=True,source_sigma_third_peak_enclosed=True,
        scope='finite arithmetic fixture only, not current point-field values')


def actual_formula_checks(field):
    """Reconstruct majorant terms from stored source-bound intermediate caps."""
    c=field.c;bd=current.Bounds(c);count=0
    def value(record):return current.LogUpper(c,record['log_absolute_upper'])
    for key,query in field.cache.items():
        p=query['primitive'];proof=p['proof']
        if proof['branch']=='exact_original_flat':continue
        t,tp,tpp=[current.JetBound(*(value(v) for v in proof[k])) for k in
            ('fixed_angle_t_C3_upper','fixed_angle_tpsi_C3_upper','fixed_angle_tpsipsi_C3_upper')]
        Phi,PhiPsi,PhiPsiPsi=[current.JetBound(*(value(v) for v in proof[k])) for k in
            ('fixed_angle_Phi_C3_upper','fixed_angle_PhiPsi_C3_upper','fixed_angle_PhiPsiPsi_C3_upper')]
        pp3=value(proof['fixed_angle_PhiPsiPsiPsi_upper']);inv=value(proof['actual_inverse_Jacobian_reciprocal_upper'])
        p1,p2,p3=[value(v) for v in proof['same_original_inverse_psi_C3_upper']]
        want=inv*bd.sum(Phi.ZZZ,bd.scale(PhiPsi.ZZ*p1,3),bd.scale(PhiPsiPsi.Z*bd.power(p1,2),3),
            pp3*bd.power(p1,3),bd.scale(bd.sum(PhiPsi.Z,PhiPsiPsi.value*p1)*p2,3))
        assert encoded(want)==encoded(p3)
        want=bd.sum(bd.constant(2*c.pi)*t.ZZZ,bd.scale(t.ZZ*p1,3),bd.scale(tp.Z*bd.power(p1,2),3),
            tpp.value*bd.power(p1,3),bd.scale(t.Z*p2,3),bd.scale(tp.value*p1*p2,3),t.value*p3)
        assert encoded(want)==encoded(proof['composed_T1_C3_upper'][3])
        A,B,E,V=[p[k] for k in ('A','B','E','V')];N0=field.phase.outer.N0
        eps=bd.constant(c.mpf(1)/N0);fac=bd.constant(c.exp(c.mpf('1.25')/N0))
        want=bd.sum(E.ZZZ*A.value,bd.scale(E.ZZ*A.Z,3),bd.scale(E.Z*bd.sum(A.ZZ,bd.power(A.Z,2)*eps),3),
            E.value*bd.sum(A.ZZZ,bd.scale(A.Z*A.ZZ*eps,3),bd.power(A.Z,3)*bd.power(eps,2)))*fac
        assert encoded(want)==encoded(query['density']['F'].ZZZ)
        count+=1
    return dict(actual_nonflat_inverse_T1_full_F_third_bound_checks=count)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C3_target_ranges_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC3TargetRanges(require_checked=False)
        assert not field.acceptance_loaded and field.target.acceptance_loaded
        assert raw['source_family']==field.identity
        signed=json.loads((current.HERE/current.target.RECEIPT).read_bytes())
        assert signed['all_passed'] and all(signed[k] for k in current.target.GATES)
        assert raw['actual_current_regular_majorant_theorem']==encoded(field.proof)
        assert raw['actual_C3_function_source_bindings']==current.source_bindings()
        calculus=independent_majorant_calculus();arithmetic=directed_arithmetic_check()
        counts=dict(actual_Z_cells=0,actual_chart_transports=0,actual_native_cells=0,own_rate_chart_steps=0,
            nonzero_quiet_third_memories=0,actual_third_target_bounds=0,exact_flat_cells=0,body_cells=0,active_flat_union_cells=0)
        assert len(raw['actual_four_Z_C3_target_transports'])==len(current.outer.CELLS)==4
        for ends,saved in zip(current.outer.CELLS,raw['actual_four_Z_C3_target_transports']):
            assert encoded(field.transport(ends))==saved and tuple(saved['exact_Z_cell'])==tuple(ends)
            assert saved['actual_five_target_C3_function_handles']==current.target.encoded(field.target.target_functions())
            for chart,window in zip(phase.CHARTS,saved['actual_17_chart_C3_transports']):
                assert window['chart']==chart and window['third_row_predecessor_never_reset']
                assert window['global_phase_cover']==[0,1];parts=current.current.native_partitions()[chart]
                assert len(window['cells'])==len(parts)-1
                for (left,right),cell in zip(zip(parts,parts[1:]),window['cells']):
                    assert cell['exact_left']==encoded(left) and cell['exact_right']==encoded(right)
                    assert cell['physical_Jacobian_applied_once'];query=field.cell(ends,chart,left,right)
                    proof=query['primitive']['proof'];assert proof['exact_same_source_P0']
                    assert proof['source_context_basis_and_ledger_retained'] and proof['live_source_phase_Z_exact_zero']
                    assert proof['actual_C3_function_handles']==current.target.encoded(field.target.primitives[chart])
                    branch=proof['branch'];counts[dict(exact_original_flat='exact_flat_cells',body_sigma_one='body_cells',
                        active_flat_union='active_flat_union_cells')[branch]]+=1
                    for key,rate in current.current.RATES.items():
                        weight=field.weights(ends,chart,left,right,rate)
                        for k in ('original_mass','cell_decay','suffix_decay'):assert current.ep(weight[k].coefficient)[0]>=0
                        if branch=='exact_original_flat':assert all(q.log is None for q in current.rows(query['density']['N_scaled'][key]))
                    counts['actual_native_cells']+=1
                for key in current.current.RATES:
                    if chart=='O3_power':
                        assert not window['actual_N_scaled_incoming_C3'][key][3]['exact_zero']
                        assert window['actual_N_scaled_local_C3'][key][3]['exact_zero']
                        assert not window['actual_N_scaled_outgoing_C3'][key][3]['exact_zero']
                        counts['nonzero_quiet_third_memories']+=1
                    counts['own_rate_chart_steps']+=1
                counts['actual_chart_transports']+=1
            for row in saved['actual_N_scaled_five_target_C3_ranges'].values():
                assert len(row)==4 and not row[3]['exact_zero'] and row[3]['source_value_or_exponential_not_materialized']
                counts['actual_third_target_bounds']+=1
            assert saved['uniform_for_all_integer_N_ge_N0'] and saved['range_caps_not_point_function_values']
            assert saved['original_independent_P0_in_source_and_not_added_twice'];counts['actual_Z_cells']+=1
        assert counts['actual_chart_transports']==68 and counts['actual_native_cells']==228
        assert counts['own_rate_chart_steps']==340 and counts['nonzero_quiet_third_memories']==20
        assert counts['actual_third_target_bounds']==20
        formula=actual_formula_checks(field)
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            quantitative_current_C3_target_ranges_installed=True,independent_majorant_calculus=calculus,
            independent_directed_arithmetic=arithmetic,actual_bound_formula_checks=formula,actual_replay_counts=counts,
            actual_C3_repaired_limit_controls_installed=False,current_numeric_point_field_oracle_installed=False,
            range_bounds_are_not_signed_function_values=True,**dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(phase.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual four-Z 17-chart C3 magnitude ranges and five-target bounds checked',flush=True)
    return result


if __name__=='__main__':run()
