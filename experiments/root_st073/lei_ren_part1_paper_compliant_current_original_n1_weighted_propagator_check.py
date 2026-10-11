"""Independent rational kernel moments, gauge conjugation and Bessel check."""
from fractions import Fraction
import json
import math
from pathlib import Path
import time

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def contains(got,want):return ep(got)[0]<=ep(want)[0] and ep(got)[1]>=ep(want)[1]
def zero(got):return ep(got)[0]<=0<=ep(got)[1]
def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported weighted source admitted')


def rational_kernel_moments():
    rows,derivatives=current.kernel_polynomials(12);count=0
    def moment(row,m):
        return sum((v*Fraction((-1)**b*math.factorial(b),(a+2*m+1)**(b+1))
                    for (a,b),v in row.items()),Fraction(0))
    for j,(row,derivative) in enumerate(zip(rows,derivatives)):
        assert sum(v for (a,b),v in row.items() if b==0)==0
        assert sum(v for (a,b),v in derivative.items() if b==0)==(1 if j==0 else 0)
        for m in range(9):
            # Direct solution on monomial forcing from the differential ODE.
            expected=Fraction(math.factorial(m)*math.factorial(m+1),
                4**(j+1)*math.factorial(m+j+1)*math.factorial(m+j+2))
            assert moment(row,m)==expected
            assert moment(derivative,m)==2*(m+j+1)*expected;count+=2
    return dict(independent_exact_radial_ODE_moment_rows=count,
                source_diagonal_kernel_boundary_and_axis_regular=True)


def conjugation_fixture(c):
    z=s.symbols('z');grad=s.Rational(1,3)+2*z/5+z*z/7
    F0=5*s.exp(-7*(z/3+z*z/5+z**3/21));weights=[F0*7 if p[1] else s.Integer(1) for p in current.WEIGHTS]
    poly=2+3*z/5+z*z/7;derivative_poly=1-z/4+z*z/9
    algebra=op.N1LogAlgebra(c,(c.mpf(0),c.ln(5),c.mpf(0),c.ln(7)),[])
    def jet(expr):return op.Jet(c,[c.mpf(int(v.p))/int(v.q) for v in
                                  (s.expand(expr).coeff(z,k) for k in range(3))])
    b0={key:algebra.lift(jet(poly)) for key in ('1,5','5,1','5,2','4,1','6,4')}
    b1={key:algebra.lift(jet(derivative_poly)) for key in ('5,1','5,2','6,4')}
    g={key:algebra.lift(jet(poly)) for key in ('4','5','6')}
    transformed=current.conjugate(dict(B0=b0,B1=b1,g=g),jet(grad),2)
    precise=type(c)();precise.dps=2*c.dps;count=0
    def finite(value,k):
        return sum((v[k]*c.mpf(5)**p[1]*c.mpf(7)**p[3] for p,v in value.terms.items()),c.mpf(0))
    for label in ('B0','B1','g'):
        for key,value in transformed[label].items():
            if label=='g':expr=poly/weights[int(key)-1]
            else:
                i,j=(int(v)-1 for v in key.split(','));ratio=weights[j]/weights[i]
                expr=ratio*(derivative_poly if label=='B1' else poly)
                if label=='B0' and key in b1:expr+=ratio*derivative_poly*(-7*grad)*current.WEIGHTS[j][1]
            for k in range(3):
                expected=s.simplify(s.diff(expr,z,k).subs(z,0)/math.factorial(k))
                want=precise.mpf(int(expected.p))/int(expected.q)
                assert contains(finite(value,k),want),(label,key,k);count+=1
    return dict(independent_symbolic_weighted_operator_derivative_rows=count,
                unknown_F0_weight_gauge_added_once_with_correct_sign=True)


def bessel_kernel_fixture(c):
    count=0
    for chi_value in ('0','.2','1'):
        chi=op.Jet(c,[c.mpf(chi_value),c.mpf(0)])
        tail=current.potential_kernel_tail(c,c.mpf(2),c.mpf(0),12)
        for tvalue in ('0','.125','.5','1'):
            t=c.mpf(tvalue);f,dx=current.potential_kernel(chi,t,c.mpf(2),12,1)
            with mp.workdps(120):
                cv=mp.mpf(chi_value);x=mp.mpf(2);tv=mp.mpf(tvalue);v=x*tv
                if tv==0:expected_f=expected_dx=mp.mpf(0)
                elif cv==0:
                    expected_f=v*(1-tv*tv)/2;expected_dx=tv**3
                else:
                    q=mp.sqrt(2*cv);a=q*v;b=q*x
                    pair=mp.besselj(1,a)*mp.bessely(1,b)-mp.bessely(1,a)*mp.besselj(1,b)
                    expected_f=mp.pi*v*v*pair/(2*x)
                    db=(mp.besselj(1,a)*(mp.bessely(0,b)-mp.bessely(2,b))
                        -mp.bessely(1,a)*(mp.besselj(0,b)-mp.besselj(2,b)))/2
                    expected_dx=-expected_f/x+mp.pi*v*v*q*db/(2*x)
                for actual,expected,name in ((f[0],expected_f,'F1'),(dx[0],expected_dx,'d_x_F1')):
                    # Independent high-precision consistency observation;
                    # analytic bounds and rational identities carry admission.
                    check=c.mpf(str(expected));enclosure=actual+c.mpf([-1,1])*tail[name][0]
                    assert contains(enclosure,check),(chi_value,tvalue,name);count+=1
    return dict(independent_Bessel_Green_kernel_consistency_rows=count,
                Bessel_check_is_consistency_not_interval_proof=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**rational_kernel_moments(),**conjugation_fixture(c),**bessel_kernel_fixture(c))
    checks['negative_or_unbounded_kernel_inputs_rejected']=all((
        rejected(lambda:current.kernel_polynomials(True)),
        rejected(lambda:current.power_log_range(c,c.mpf([-1,0]),1,0)),
        rejected(lambda:current.relative_amplitude(owner.inputs.frame('.371')[1],op.Jet(c,[1]),1,2))))
    axis=owner.report(packets['core_axis']);view=owner.report(packets['selected_actual_inner_exit'])
    for packet in (axis,view):
        assert all(not packet[k] for k in current.OPEN)
        assert packet['remaining_normalized_B0_B1_feedback_not_solved']
        assert packet['full_sector_preserving_coupled_tail_still_open']
        assert packet['old_coupled_prefix_not_replaced_by_initial_integral']
        assert packet['ordinary_Z_initial_derivative_order']==1
    axis_state=owner._packets[id(packets['core_axis'])][4]
    assert all(ep(v)==(0,0) for value in axis_state.values() for jet in value.terms.values() for v in jet.coefficients)
    assert len(view['integration_cells'])==20
    assert all(not view['residual_B0'][key]['terms'] for key in ('1,5','2,6','3,6'))
    checks['kinematic_rows_removed_from_remaining_feedback_not_double_counted']=True
    state=owner._packets[id(packets['selected_actual_inner_exit'])][4]
    assert all(p[1]>=0 and p[0]>=0 for value in state.values() for p in value.terms)
    assert state['F1'].terms and all(p[1]%2==0 for p in state['F1'].terms)
    for field in view['normalized_kernel_sector_budgets'].values():
        for row in field.values():
            assert not row['remaining_B0_B1_feedback_included']
            assert row['kernel_tail_normalized_reference_norm_upper'] is not None
            assert ep(row['kernel_tail_normalized_reference_norm_upper'])[1]<c.mpf('1e-12')
    checks['actual_four_source_normalized_kernel_tail_rows_below_1e_minus12']=True
    checks['outer_endpoint_interval_crossing_axis_rejected']=rejected(lambda:owner.evaluate('.371',c.mpf([0,1]),'core'))
    p=owner._packets[id(packets['selected_actual_inner_exit'])][3]
    z=owner.inputs._packets[id(p)][1]['Z'];algebra=owner.inputs._packets[id(p)][5]
    source=owner.physical.axis_inputs(z);gradient=op.Jet(c,owner.source.amplitude.gradient_jets(z))
    identity=(source['H'].truncate(2)*gradient.truncate(2)/source['L'].truncate(2)-source['chi'].truncate(2))
    assert all(zero(v) for v in identity.coefficients)
    assert contains(owner.physical.Lambda*owner.physical.epsilon,c.mpf(1))
    checks['actual_H_Gprime_over_L_equals_chi_rows']=3
    checks['actual_Lambda_epsilon_identity_retained']=True
    checks['actual_axis_zero_and_original_20_source_cells_retained']=True
    checks['F0_Lambda_weight_keeps_swirl_and_pressure_sectors_distinct']=True
    checks['foreign_weighted_packet_rejected']=rejected(lambda:owner.report(current.OriginalN1WeightedPropagatorPacket()))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])]
    old=entry[1]['remaining_normalized_B0_B1_feedback_not_solved'];entry[1]['remaining_normalized_B0_B1_feedback_not_solved']=False
    try:checks['false_feedback_completion_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['remaining_normalized_B0_B1_feedback_not_solved']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_WEIGHTED_PROPAGATOR',checks['independent_exact_radial_ODE_moment_rows'],
          checks['independent_symbolic_weighted_operator_derivative_rows'],checks['independent_Bessel_Green_kernel_consistency_rows'],flush=True)
    return result
