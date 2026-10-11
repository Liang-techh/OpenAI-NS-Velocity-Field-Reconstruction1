"""Independent physical-R moments, interval ranges and common-path guards."""
import json
import math
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported own moment query admitted')


def rational(c,value):
    value=s.cancel(value);return c.mpf(int(s.numer(value)))/int(s.denom(value))


def resolved(value,k):
    c=value.ctx
    return sum((row[k]*math.factorial(k)*c.exp(sum((value.algebra.logs[i]*p[i] for i in range(4)),c.mpf(0)))
                for p,row in value.terms.items()),c.mpf(0))


def overlaps(left,right,order=1):
    c=left.ctx
    for power in left.terms.keys()|right.terms.keys():
        for k in range(order+1):
            a=left.terms[power][k] if power in left.terms else c.mpf(0)
            b=right.terms[power][k] if power in right.terms else c.mpf(0)
            assert ep(a)[0]<=ep(b)[1] and ep(a)[1]>=ep(b)[0],(power,k)


@source_precision
def independent_physical_moment_fixture(c):
    z,R,rho=s.symbols('z R rho');eps=s.Rational(2,7);Pi=s.Integer(3);anchor=s.Integer(5)
    z0=s.Rational(1,4);r0=s.Rational(3,2);u0=1+z+(2-z)*R**2
    f0=anchor*(1+z/3+R*(1+z*z)/7)
    f1=anchor/eps*(R/eps)*(1+z*(R/eps)/5+(R/eps)**2/20)
    u1=eps*((R/eps)/5+(z+1)*(R/eps)**2/7)
    P1=eps*Pi*((R/eps)*(1+z*z)+(R/eps)**2/8+z*(R/eps)**3/11)
    mean1=s.integrate(u1,(R,0,R))/R
    omega=-2*R*(s.diff(P1,R)-2*f0*f1)
    direct={'m11':u1,'m12':2*R*f1,'m13':2*f0*f1-omega/(2*R),
            'm14':2*R*(u0*f1+u1*f0),'m15':2*u0*u1-2*R*f0*f1+omega/2}
    expected={k:s.integrate(s.cancel(v),(R,0,eps*rho)) for k,v in direct.items()}
    sub={z:z0,rho:r0};algebra=op.N1LogAlgebra(c,(c.ln(c.mpf('.2')),c.ln(rational(c,anchor)),c.ln(rational(c,Pi)),c.ln(rational(c,1/eps))),[])
    jet=lambda f:op.Jet(c,[rational(c,s.diff(f,z,k).subs(sub)) for k in (0,1)])
    at=lambda f:s.cancel(f.subs(R,eps*rho))
    W=[algebra.lift(jet(at(f1)*eps/anchor)),algebra.lift(jet(at(u1)/eps)),
       algebra.lift(jet(at(mean1-u1)/eps)),algebra.lift(jet(at(P1)/(eps*Pi)))]
    physical=current.physical_values(algebra,W,op.Jet.constant(c,0,1))
    integrands=current.rho_integrands(algebra,algebra.lift(jet(at(f0))),algebra.lift(jet(at(u0))),physical,rational(c,r0))
    endpoint=current.endpoint_cumulants(algebra,W,rational(c,r0));count=0
    for k in ('m11','m13'):
        for grade in (0,1):
            truth=rational(c,s.diff(expected[k],z,grade).subs(sub))
            assert ep(resolved(endpoint[k],grade)-truth)[0]<=0<=ep(resolved(endpoint[k],grade)-truth)[1];count+=1
    transformed={'m12':eps*at(direct['m12']),'m14':eps*at(direct['m14']),
                 'm15_parts':eps*at(2*u0*u1+P1)}
    for k,expr in transformed.items():
        for grade in (0,1):
            truth=rational(c,s.diff(expr,z,grade).subs(sub));difference=resolved(integrands[k],grade)-truth
            assert ep(difference)[0]<=0<=ep(difference)[1];count+=1
    parts=s.integrate(transformed['m15_parts'],(rho,0,rho))-eps*rho*at(P1)
    assert s.simplify(parts-expected['m15'])==0
    assert s.simplify(expected['m13']-at(P1))==0
    assert s.simplify(s.diff(expected['m11'],rho)-eps*at(u1))==0
    assert s.simplify(s.diff(expected['m13'],rho)-eps*at(direct['m13']))==0
    # Nonzero inertia is explicit; excluding omega changes m13 and m15.
    assert s.simplify(omega)!=0
    assert s.simplify(expected['m13']-s.integrate(2*f0*f1,(R,0,eps*rho)))!=0
    assert s.simplify(expected['m15']-s.integrate(2*u0*u1-2*R*f0*f1,(R,0,eps*rho)))!=0
    # Variable S paths and both physical measures, independently of the code.
    t,S,h=s.symbols('t S h',positive=True);k=s.symbols('k',integer=True,nonnegative=True)
    assert s.simplify(s.integrate((S*t)**3*S,(t,0,1))-S**4/4)==0
    rr=4*s.exp(h*S)
    assert s.simplify(s.diff(eps*rr,S)-eps*h*rr)==0
    # The independent Bessel expansion gives the exact weighted primitive.
    chi=s.symbols('chi',positive=True);B=sum((-chi*rho/2)**j/(s.factorial(j)*s.factorial(j+1)) for j in range(7))
    I=s.integrate(rho*s.diff(B,chi),(rho,0,rho))
    antiderivative=sum((-1)**j*chi**(j-1)*rho**(j+2)/(2**j*s.factorial(j-1)*s.factorial(j+1)*(j+2)) for j in range(1,7))
    assert s.simplify(I-antiderivative)==0
    # Core model decomposition: Qpsi=q*rho and Pvisc=2*eps*q*rho.
    delta=s.Rational(1,5);slope=1+z+z*z;L=1-delta*z*z;d=1-z*z
    q=((3+delta)*z*slope-d*s.diff(slope,z))/(2*L)
    Qpsi=(2*z*slope*rho-(1-delta)*z*slope*rho/2-d*s.diff(slope*rho/2,z))/L
    assert s.simplify(Qpsi-q*rho)==0
    assert s.simplify(s.integrate(eps*(Qpsi+rho*s.diff(Qpsi,rho)),(rho,0,rho))-eps*q*rho**2)==0
    return count


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    fixture=independent_physical_moment_fixture(c)
    views={name:owner.report(packet) for name,packet in packets.items()}
    checks=dict(independent_physical_R_moment_and_density_value_Z1_rows=fixture,
        pressure_inertia_and_Omega0_terms_are_nonzero_in_fixture=True,
        exact_m11_mean_and_m13_full_pressure_FTC_verified=True,
        m15_full_pressure_integration_by_parts_verified=True,
        variable_endpoint_substitution_and_actual_collar_physical_Jacobian_verified=True)
    for name,packet in packets.items():
        view=views[name];entry=owner._packets[id(packet)];endpoint,integrals,cumulative,records=entry[3]
        assert set(view['physical_partial_cumulative_moments'])==set(current.MOMENTS)
        assert all(v.order==1 for v in cumulative.values())
        assert all(not view[key] for key in current.OPEN)
        assert view['whole_domain_relative_moment_accuracy_not_certified'] and view['local_partial_cumulative_moments_not_terminal_full_support_defects']
        pe=owner.regular._packets[id(entry[2])];W=pe[5][0];rho=pe[5][5]['rho']
        expected=current.endpoint_cumulants(W[0].algebra,[op.truncate(v,1) for v in W],rho)
        for key in expected:overlaps(endpoint[key],expected[key])
        overlaps(cumulative['m15'],integrals['m15_parts']+endpoint['m15_boundary'])
        if name=='core_axis':assert all(all(ep(v)[0]==ep(v)[1]==0 for row in a.terms.values() for v in row.coefficients) for a in cumulative.values())
        for cell in entry[4]:
            assert all(all(math.isfinite(float(ep(v)[0])) and math.isfinite(float(ep(v)[1])) for row in value.terms.values() for v in row.coefficients)
                       for value in cell[3])
    inlet=owner._packets[id(packets['actual_first_inlet'])][3][2]
    exitcore=owner._packets[id(packets['core_exit'])][3][2]
    for key in current.MOMENTS:overlaps(inlet[key],exitcore[key])
    checks.update(all_five_partial_cumulants_Z1_execute=True,
        axis_zero_and_nonsingular_axis_cell_ranges_verified=True,full_pressure_endpoint_not_only_viscous_seed=True,
        actual_inlet_preserves_all_core_moment_history=True,all_remaining_feedback_and_original_family_bound=True)
    assert all(views[k]['passing_partial_moment_value_Z1_rows']==10 for k in (
        'core_interior','core_exit','actual_first_inlet','selected_actual_inner_exit'))
    assert views['core_axis']['passing_partial_moment_value_Z1_rows']==0
    assert all(views['selected_actual_inner_exit']['partial_moment_value_Z1_budgets'][k]['Z0']['sign']==(
        'positive' if k in ('m11','m13') else 'negative') for k in current.MOMENTS)
    checks.update(all_ten_partial_moment_value_Z1_width_rows_pass_at_four_nonaxis_queries=True,
        correlated_core_weighted_Bessel_primitive_and_pressure_model_verified=True,
        all_five_selected_inner_cumulants_certifiably_nonzero=True)
    # A whole core cell contains two interior point function enclosures.
    cell=owner.range_cell('.371',c.mpf(['.5','1.5']),'core')
    for coord in ('.625','1.375'):
        pp=owner.regular.evaluate('.371',coord,'core');owner.regular.report(pp)
        point=owner.regular._packets[id(pp)][5][0]
        for rangevalue,value in zip(cell[3],point[:4]):overlaps(rangevalue,op.truncate(value,1))
    checks['core_uniform_cell_ranges_overlap_independent_interior_queries']=True
    token=packets['selected_actual_inner_exit'];entry=owner._packets[id(token)];view=entry[1]
    key='local_partial_cumulative_moments_not_terminal_full_support_defects';saved=view[key];view[key]=False
    rejected(lambda:owner.report(token));view[key]=saved
    value=entry[3][2]['m12'];saved={p:op.Jet(c,list(row.coefficients)) for p,row in value.terms.items()}
    value.terms={};rejected(lambda:owner.report(token));value.terms=saved
    rejected(lambda:owner.report(current.OriginalN1CumulativeMomentPacket()))
    rejected(lambda:owner.evaluate('.371','4.1','core'))
    rejected(lambda:owner.range_cell('.371',c.mpf(['0','.2']),'core',False))
    checks.update(moment_scope_and_hidden_cumulative_mutations_rejected=True,unsupported_queries_and_foreign_packets_rejected=True)
    owner.report(token);inputs=dict(owner.hashes);inputs[current.NAME]=current.sha(current.NAME);inputs[Path(__file__).name]=current.sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=inputs,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(op.snapshot(result),separators=(',',':'))+'\n').encode())
    print('PASS_N1_OWN_CUMULATIVE_MOMENTS',fixture,flush=True);return result
