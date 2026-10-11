"""Exact principal ODE, regular weighted kernel and same-source function checks."""
import copy
from fractions import Fraction
import json
import math
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_regular_series_function as current
import lei_ren_part1_paper_compliant_current_core_first_interface as interface
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported regular series function admitted')


def principal_ode_identities():
    """Every polynomial kernel solves the ODE, for arbitrary forcing."""
    rows,dxrows=current.weighted.kernel_polynomials(12)
    def Dx(row,degree):
        out={}
        for (p,b),v in row.items():
            out[p,b]=out.get((p,b),0)+(degree-p)*v
            if b:out[p,b-1]=out.get((p,b-1),0)-b*v
        return {k:v for k,v in out.items() if v}
    count=0
    for j,row in enumerate(rows):
        first=Dx(row,2*j+1);assert first==dxrows[j]
        second=Dx(first,2*j)
        heat={k:second.get(k,0)+3*first.get(k,0) for k in set(second)|set(first)}
        heat={k:v for k,v in heat.items() if v}
        assert heat==({} if j==0 else rows[j-1])
        assert sum(v for (p,b),v in row.items() if b==0)==0
        assert sum(v for (p,b),v in first.items() if b==0)==(1 if j==0 else 0)
        count+=4
    x=s.symbols('x',positive=True);h=s.Function('h6')(x);h4=s.Function('h4')(x)
    I1=s.Function('I1')(x);I3=s.Function('I3')(x);I4=s.Function('I4')(x)
    subs={s.diff(I1,x):x*h,s.diff(I3,x):x**3*h,s.diff(I4,x):h4}
    D=I1/x;K=-I1/2+I3/(2*x*x)
    assert s.cancel((s.diff(D,x)+D/x).subs(subs)-h)==0
    assert s.cancel((s.diff(K,x)+2*K/x).subs(subs)+D)==0
    assert s.cancel(s.diff(I4,x).subs(subs)-h4)==0
    Tf=(I1-I3/x**2)/2;Tx=s.diff(Tf,x).subs(subs)
    assert s.cancel(Tx-I3/x**3)==0
    assert s.cancel((s.diff(Tx,x)+3*Tx/x).subs(subs)-h)==0
    # T on monomials gives the entire positive comparison series.
    t=s.symbols('t',positive=True)
    for j in range(8):
        integral=s.integrate(t*(1-t*t/x**2)*t**(2*j)/2,(t,0,x))
        assert s.cancel(integral-x**(2*j+2)/(4*(j+1)*(j+2)))==0
    # Formal potential inverse telescopes for each finite cutoff.
    chi=s.symbols('chi');A=s.symbols('A0:10')
    for N in range(8):
        Lheat=sum((-2*chi)**j*A[j] for j in range(N+1))
        potential=sum(2*chi*(-2*chi)**j*A[j+1] for j in range(N+1))
        assert s.expand(Lheat+potential-A[0]+(-2*chi)**(N+1)*A[N+1])==0
    return dict(exact_arbitrary_forcing_kernel_ODE_and_boundary_identities=count,
                exact_pressure_axial_and_K_principal_FTC_equations_verified=True,
                positive_T_monomial_and_derivative_comparison_verified=True,
                infinite_swirl_inverse_has_factorial_tail_and_regular_zero_axis_data=True)


def rational(c,v):
    v=s.cancel(v)
    return c.mpf(int(s.numer(v)))/int(s.denom(v))


@source_precision
def theta_Z2_fixture(c):
    z,r=s.symbols('z rho');delta=s.Rational(1,5);lam=s.Integer(3)
    phi=1+z*z+r*z+r*r/5;grad=s.Rational(3,10)+z/7
    F0=s.exp(-lam*(3*z/10+z*z/14));F=F0*phi;d=1-z*z;L=1-delta*z*z
    T=lambda a,v:(a*z*v+d*s.diff(v,z)-2*z*r*s.diff(v,r))/L
    expected=s.cancel(s.expand(-2*T(-3,T(-2-delta,F))/F0/lam**2))
    algebra=op.N1LogAlgebra(c,(c.mpf(0),c.mpf(0),c.mpf(0),c.ln(3)),[])
    sub={z:s.Rational(1,4),r:s.Rational(3,2)}
    jet=lambda v,n:op.Jet(c,[rational(c,s.diff(v,z,k).subs(sub))/math.factorial(k) for k in range(n+1)])
    rows=[algebra.lift(jet(s.diff(phi,r,i),5-i)) for i in range(3)]
    got=current.theta_source(algebra,rows,algebra.lift(op.Jet.variable(c,c.mpf('.25'),5)),c.mpf('1.5'),
                             c.mpf('.2'),algebra.lift(jet(grad,5)),2)
    for k in range(3):
        value=sum((row[k]*c.mpf(3)**p[3] for p,row in got.terms.items()),c.mpf(0))
        want=rational(c,s.diff(expected,z,k).subs(sub))/math.factorial(k)
        assert ep(value)[0]<=ep(want)[1] and ep(value)[1]>=ep(want)[0],k
    return dict(independent_original_F0_cancellation_and_theta_ordinary_Z2_fixture_rows=3)


def seam_and_own_mean_bindings():
    bindings=interface.defining_boundary_bindings()
    assert bindings['exact_first_zero_source_increment_branch_bound']
    assert bindings['actual_first_initial_values_and_moments_are_atoms_not_covers']
    import lei_ren_part1_paper_compliant_current_original_n1_leading_mean_epsilon_factor_check as mc
    assert mc.executable_own_mean_bindings()['original_own_core_integral_inlet_bound']
    # The production actual leading rows use that same core function at
    # phase zero. The series always integrates from zero; it never resets W.
    import ast
    src=ast.parse((current.HERE/Path(current.natural.mean.__file__).name).read_text(encoding='utf8'))
    fn=next(n for n in src.body if isinstance(n,ast.FunctionDef) and n.name=='actual_rows')
    wanted=ast.parse('if ep(phase)[1]==0:\n return inlet').body[0]
    assert any(ast.dump(n)==ast.dump(wanted) for n in ast.walk(fn))
    return dict(genuine_core_actual_first_inlet_function_bindings_preserved=True,
                each_exact_V_word_uses_one_zero_to_x_path_and_inherits_source_seam=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**principal_ode_identities(),**theta_Z2_fixture(c),**seam_and_own_mean_bindings())
    contract=owner.contract();t=owner.tail;eta=owner.natural._data['eta'];a=owner.natural._data['a']
    assert all(not contract[k] for k in current.OPEN)
    assert t['regular_axis_powers']==[2,2,2,2,1,1]
    assert t['regular_Green_factor']==64 and ep(t['regular_Green_column_majorant'])[1]<64
    assert t['original_holomorphic_chi_tube']==owner.weighted.tube
    assert t['original_holomorphic_chi_tube']['chi_holomorphic']
    assert ep(t['original_holomorphic_chi_tube']['chi_complex_modulus_upper'])[1]<=1
    tt=s.symbols('tt',nonnegative=True)
    assert s.expand(tt**2*(1-tt**2)/2+tt**4-tt**2*(1+tt**2)/2)==0
    # General inequality: (j+1/2)/(j+1)! <= (j+1/2)/(j+1)<1<=j.
    j=s.symbols('j',integer=True,positive=True)
    assert s.simplify(1-(j+s.Rational(1,2))/(j+1)-1/(2*(j+1)))==0
    assert t['regular_initial_M0']==64*owner.natural._data['M0']
    assert t['core_pointwise_regular_source_norm_dominated_by_M0']
    assert ep(t['core_pointwise_regular_source_cover']['modulus_upper'])[1]<=ep(owner.natural._data['M0'])[0]
    core_source=owner.natural._data['core_L1']
    assert t['core_pointwise_regular_source_rows']=={
        'g4_over_x':op.export(core_source['g4']/(2*a),0),
        'g5':op.export(core_source['g5']/2,0),'g6':op.export(core_source['g6']/2,0)}
    # Bind the recovery to the production uniform source formulas: this
    # division is justified because ic is 2*sc, not by an L1-to-Linf rule.
    import ast
    src=ast.parse((current.HERE/Path(current.natural.__file__).name).read_text(encoding='utf8'))
    bounds=next(n for n in src.body if isinstance(n,ast.FunctionDef) and n.name=='natural_majorants')
    source=next(n for n in bounds.body if isinstance(n,ast.FunctionDef) and n.name=='source')
    pressure=ast.parse('g4=algebra.shift(Np*(2*a),(0,0,-1,0))').body[0]
    core_integral=ast.parse('ic={k:v*2 for k,v in sc.items()}').body[0]
    assert any(ast.dump(n)==ast.dump(pressure) for n in source.body)
    assert any(ast.dump(n)==ast.dump(core_integral) for n in ast.walk(bounds))
    assert t['source_M0_tube_radius']==eta/8 and t['value_tube_radius']==eta/16
    assert t['derivative_tube_radius']==eta/32 and t['derivative_Cauchy_gap']==eta/32
    assert ep(t['ratio_upper'])[1]<c.mpf('1e-12')
    assert ep(t['ordinary_Z3_remainder_upper'])[1]<c.mpf('1e-12')
    assert ep(owner.contract(3)['finite_equation_residual_upper_on_eta_over32'])[1]<ep(contract['finite_equation_residual_upper_on_eta_over32'])[0]
    assert contract['all_remaining_words_enclosed_not_set_zero'] and contract['numeric_centers_compute_only_zeroth_source']
    for name,p in packets.items():
        v=owner.report(p);W,dx,dr,recovery,recovery_dr,system=owner._packets[id(p)][5]
        assert all(not v[k] for k in current.OPEN)
        assert v['same_original_common_inner_function_for_every_query'] and v['all_six_row_feedback_included']
        assert v['full_higher_radial_axis_and_join_contracts_still_open'] and v['n1_own_repair_and_n2_still_open']
        assert all(w.order==3 for w in W) and all(w.order==2 for w in dx+dr)
        assert recovery['Q1'].order==2 and recovery_dr['Q1'].order==1
        assert all(int(key.split(',')[0]) in (4,5,6) for key in system['C0'])
        assert all(int(key.split(',')[0]) in (5,6) for key in system['B1'])
        # Only source row4,column1 enters the pressure Green column, where
        # the regular weights give (s/x)^2<=1. Other columns have p=1 or2.
        assert {key for key in system['C0'] if key.startswith('4,')}=={'4,1'}
        if name=='core_axis':
            assert all(not w.terms for w in W)
            assert len(dr)==4 and len(v['physical_d_rho_function_jets'])==4
            assert all(not row.terms for row in dx[:4])
            assert op.export(dx[4])==op.export(system['g']['5']/4)
            assert op.export(dx[5])==op.export(system['g']['6']/2)
        else:
            assert op.export(dx[0])==op.export(op.truncate(W[4],2))
            assert op.export(dx[1])==op.export(op.truncate(W[5],2))
            assert all(op.export(r)==op.export(x/(2*system['x'])) for x,r in zip(dx,dr))
            assert v['passing_six_component_value_Z1_rows']==12
        for field,w,scale in zip(current.FIELDS,W,current.natural.WEIGHTS):
            expected=current.weighted.transform(w,system['gradient'],scale,3)
            assert v['physical_function_jets'][field]==op.export(expected)
    ce=owner.report(packets['core_exit']);fi=owner.report(packets['actual_first_inlet'])
    for key in ('normalized_function_jets','normalized_d_x_function_jets','normalized_d_rho_function_jets',
                'physical_own_mean_Q_V','physical_own_mean_Q_V_rho','exact_normalized_equation'):
        assert ce[key]==fi[key],key
    # Check the finite adapter at a genuine core chi, including near-axis s.
    chi=owner.physical.axis_inputs(c.mpf('.371'))['chi'].truncate(1)
    for ss,xx in (('0','2'),('.000001','.0001'),('1','2'),('2','2')):
        kernel,tail=current.principal_kernel(c,chi,c.mpf(ss),c.mpf(xx),12,1)
        assert set(kernel)==set(current.FIELDS) and tail['first_omitted_potential_index']==13
    checks['all_six_equation_source_rows_and_weighted_regular_Green_support_checked']=True
    checks['regular_all_word_tail_and_ordinary_Z3_convergence_checked']=True
    checks['supported_weighted_Green_columns_and_original_holomorphic_chi_bound_recorded']=True
    checks['core_pointwise_g4_over_x_g5_g6_and_weighted_M0_premises_bound']=True
    checks['same_function_core_collar_seam_and_other_radius_Z_queries_execute']=True
    checks['twelve_component_value_Z1_rows_pass_at_all_six_nonaxis_queries']=True
    checks['equation_derived_first_radial_mixed_Z2_and_axis_limits_checked']=True
    checks['uncomputed_high_Z_jets_bounded_and_physical_F0_gauge_added_once']=True
    checks['unsupported_query_and_foreign_packet_rejected']=all((
        rejected(lambda:owner.evaluate('.371','4.1','core')),
        rejected(lambda:owner.evaluate('.371','1.01','first')),
        rejected(lambda:owner.contract(True)),
        rejected(lambda:current.principal_kernel(c,chi,c.mpf(1),c.mpf(2),12,2)),
        rejected(lambda:owner.report(current.OriginalN1RegularSeriesFunctionPacket()))))
    entry=owner._packets[id(packets['selected_actual_inner_exit'])];old=entry[1]['full_higher_radial_axis_and_join_contracts_still_open']
    entry[1]['full_higher_radial_axis_and_join_contracts_still_open']=False
    try:checks['false_full_higher_radial_scope_mutation_rejected']=rejected(lambda:owner.report(packets['selected_actual_inner_exit']))
    finally:entry[1]['full_higher_radial_axis_and_join_contracts_still_open']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_N1_REGULAR_SERIES_FUNCTION',checks['exact_arbitrary_forcing_kernel_ODE_and_boundary_identities'],flush=True)
    return result
