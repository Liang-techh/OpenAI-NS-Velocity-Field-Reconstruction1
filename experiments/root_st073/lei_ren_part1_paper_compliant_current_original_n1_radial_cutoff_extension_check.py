"""Independent physical-R cutoff/recovery identities and live source guards."""
import json
import math
from pathlib import Path
import time

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported cutoff query admitted')


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


def exact_zero(value):
    return all(ep(v)[0]==ep(v)[1]==0 for row in value.terms.values() for v in row.coefficients)


def physical_overlap(left,right,order=1):
    # Different cumulative partitions need not share identical hb-sector
    # decompositions. Compare the actual sum at the fixed original scales.
    difference=left-right
    for k in range(order+1):
        budget=current.parent.moment_budget(difference,k)
        if budget['exact_zero']:continue
        interval=budget['signed_physical_partial_moment_normalized_enclosure']
        assert interval is not None and ep(interval)[0]<=0<=ep(interval)[1],k


@source_precision
def independent_physical_fixture(c):
    z,R,rho=s.symbols('z R rho');eps=s.Rational(2,7);hb=s.Rational(2,5)
    Pi=s.Integer(3);anchor=s.Integer(5);z0=s.Rational(1,4);r0=s.Rational(3,2)
    u0=1+z+(2-z)*R**2;f0=anchor*(1+z/3+R*(1+z*z)/7)
    f1=anchor/eps*(R/eps)*(1+z*(R/eps)/5+(R/eps)**2/20)
    u1=eps*((R/eps)/5+(z+1)*(R/eps)**2/7)
    sub={z:z0,rho:r0};algebra=op.N1LogAlgebra(c,(c.ln(rational(c,hb)),c.ln(rational(c,anchor)),
        c.ln(rational(c,Pi)),c.ln(rational(c,1/eps))),[])
    jet=lambda f,n:op.Jet(c,[rational(c,s.diff(f,z,k).subs(sub))/math.factorial(k) for k in range(n+1)])
    at=lambda f:s.cancel(f.subs(R,eps*rho))
    change=s.Rational(-2,5)
    physical={'F1':algebra.lift(jet(at(f1),1)),'Uz1':algebra.lift(jet(at(u1),1))}
    densities=current.change_densities(algebra,physical,algebra.lift(jet(at(f0),1)),
        algebra.lift(jet(at(u0),1)),rational(c,r0),rational(c,change))
    direct={'m11':change*u1,'m12':change*2*R*f1,'m13':change*2*f0*f1,
        'm14':change*2*R*(u0*f1+u1*f0),'m15':change*(2*u0*u1-2*R*f0*f1)}
    count=0
    for name,expr in direct.items():
        for k in (0,1):
            truth=rational(c,s.diff(eps*at(expr),z,k).subs(sub));diff=resolved(densities[name],k)-truth
            assert ep(diff)[0]<=0<=ep(diff)[1],(name,k);count+=1
    # Independent polynomial cutoff fixture: own M determines mean and V.
    chi=(1-R)**3;cutU=chi*u1;M=s.integrate(cutU,(R,0,R));delta=s.Rational(1,5)
    L=1-delta*z*z;d=1-z*z
    V=(2*z*R*cutU-(1+delta)*z*M-d*s.diff(M,z))/L
    assert s.simplify(L*s.diff(V,R)-(1-delta)*z*cutU-2*z*R*s.diff(cutU,R)+d*s.diff(cutU,z))==0
    recovered=current.recover_velocity(algebra,algebra.lift(jet(at(cutU),2)),algebra.lift(jet(at(M),2)),
        rational(c,z0),rational(c,r0),rational(c,delta))
    for name,truth in {'mean1':M/R,'K1':M/R-cutU,'V1':V}.items():
        for k in range(recovered[name].order+1):
            diff=resolved(recovered[name],k)-rational(c,s.diff(at(truth),z,k).subs(sub))
            assert ep(diff)[0]<=0<=ep(diff)[1],(name,k);count+=1
    # The fixed leading pressure forcing survives the velocity cutoff.
    N=z+R;P=s.integrate(2*f0*f1+N,(R,0,R));deltaP=s.integrate((chi-1)*2*f0*f1,(R,0,R))
    assert s.simplify(s.diff(P+deltaP,R)-2*f0*chi*f1-N)==0
    assert s.simplify((P+deltaP-chi*P).subs(R,1))!=0
    # Compact velocities are supported in [0,1], Rin=2, leading force in [0,3].
    pressure_inner=s.integrate(2*f0*chi*f1,(R,0,1))+s.integrate(N,(R,0,2))
    pressure_full=s.integrate(2*f0*chi*f1,(R,0,1))+s.integrate(N,(R,0,3))
    fifth_inner=s.integrate(2*u0*chi*u1-2*R*f0*chi*f1,(R,0,1))-s.integrate(R*N,(R,0,2))
    fifth_full=s.integrate(2*u0*chi*u1-2*R*f0*chi*f1,(R,0,1))-s.integrate(R*N,(R,0,3))
    assert s.simplify(pressure_full-pressure_inner-s.integrate(N,(R,2,3)))==0
    assert s.simplify(fifth_full-fifth_inner+s.integrate(R*N,(R,2,3)))==0
    assert s.simplify(pressure_full-pressure_inner)!=0 and s.simplify(fifth_full-fifth_inner)!=0
    return count


@source_precision
def independent_cutoff_derivatives(c):
    h=mp.mpf(2)/5;S=mp.mpf(3)/7;q=mp.mpf(61)/100;r=4*mp.exp(h*S*q);z=mp.mpf(1)/4
    algebra=op.N1LogAlgebra(c,(c.ln(c.mpf(2)/5),c.mpf(0),c.mpf(0),c.mpf(0)),[])
    chi=current.cutoff_radial_jets(algebra,c.mpf(61)/100,c.mpf(str(r)),c.mpf(3)/7)
    def chi_function(rr):
        x=4*mp.log(rr/4)/(h*S)-2
        return 1/(1+mp.exp(1/(1-x)**2-1/x**2))
    tol=mp.mpf('1e-250');count=0
    for k,v in enumerate(chi):
        a,b=ep(resolved(v,0));truth=mp.diff(chi_function,r,k)
        assert a-tol<=truth<=b+tol,(k,truth);count+=1
        assert all(ep(row[j])[0]==ep(row[j])[1]==0 for row in v.terms.values() for j in (1,2))
    assert (-1,0,0,0) in chi[1].terms
    assert {(-1,0,0,0),(-2,0,0,0)}.issubset(chi[2].terms)
    def f(rr,zz):return rr**2*(1+zz)+rr**3*zz**2
    jets=[algebra.lift(op.Jet(c,[c.mpf(str(mp.diff(lambda rr:f(rr,z),r,k))),
        c.mpf(str(mp.diff(lambda rr:rr**2+2*rr**3*z,r,k)))])) for k in range(3)]
    pairs=current.cut_pair(*jets,chi)
    for k,name in enumerate(('value','rho1','rho2')):
        for j in (0,1):
            truth=mp.diff(lambda rr:(f(rr,z) if j==0 else rr**2+2*rr**3*z)*chi_function(rr),r,k)
            a,b=ep(resolved(pairs[name],j));assert a-tol<=truth<=b+tol,(name,j);count+=1
    # The geometric midpoint is no later than the requested arithmetic midpoint.
    a=mp.exp(h*S/2);b=mp.exp(h*S)
    assert mp.exp(h*S*mp.mpf('.75'))<=(a+b)/2
    return count


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    fixture=independent_physical_fixture(c);derivatives=independent_cutoff_derivatives(c)
    views={k:owner.report(p) for k,p in packets.items()}
    checks=dict(independent_physical_R_density_and_own_recovery_Z_rows=fixture,
        independent_smooth_cutoff_rho0_rho1_rho2_and_product_rule_Z1_rows=derivatives,
        explicit_hb_inverse_inverse_squared_and_sigma_Taylor_second_coefficient_verified=True,
        divergence_recovery_from_own_cutoff_axial_primitive_verified=True,
        pressure_recovery_and_nonzero_fixed_leading_forcing_tail_signs_verified=True,
        cutoff_zero_radius_before_paper_arithmetic_midpoint_verified=True)
    for name,packet in packets.items():
        view=views[name];entry=owner._packets[id(packet)];chi,pairs,physical,M,own,P,totals=entry[2]
        assert set(view['full_velocity_support_moments'])=={'m11','m12','m14'}
        assert view['unresolved_leading_pressure_forcing_tails']['both_tails_uncomputed_and_not_set_zero']
        assert all(not view[k] for k in current.OPEN)
        assert all(v.order==1 for rows in physical.values() for v in rows.values())
        assert M.order==2 and own['V1'].order==1 and P.order==1
        if name in ('zero_endpoint','inner_exit'):
            assert all(exact_zero(v) for rows in physical.values() for v in rows.values())
        if name in ('core_axis','core_exit','keep_endpoint'):
            assert ep(resolved(chi[0],0))==(1,1) and exact_zero(chi[1]) and exact_zero(chi[2])
            assert all(exact_zero(v) for v in entry[3][4][0].values())
        if name=='core_axis':assert exact_zero(M) and all(exact_zero(v) for v in own.values())
        assert all(row['relative_partial_moment_width_certified'] for budgets in view['full_velocity_support_moment_budgets'].values() for row in budgets.values())
    zero=owner._packets[id(packets['zero_endpoint'])][2];exit=owner._packets[id(packets['inner_exit'])][2]
    overlaps(zero[3],exit[3],2);overlaps(zero[4]['V1'],exit[4]['V1'])
    for name in ('m11','m12','m14'):overlaps(zero[6][name],exit[6][name])
    # The full m11 equals the same terminal primitive, with an independently integrated enclosure.
    physical_overlap(exit[6]['m11'],op.truncate(exit[3],1))
    checks.update(six_original_queries_with_own_recovery_execute=True,
        untouched_keep_and_exact_zero_velocity_radial_jets_after_cutoff=True,
        own_M_and_V_constant_after_cutoff_verified=True,
        three_full_velocity_support_moments_have_all_six_value_Z1_width_rows=True,
        full_first_moment_matches_terminal_own_M=True,
        remaining_two_leading_force_tails_and_repair_recursion_global_flags_not_admitted=True)
    token=packets['cutoff_interior'];entry=owner._packets[id(token)]
    key='mean_V_pressure_recovered_not_multiplied_by_cutoff';saved=entry[1][key];entry[1][key]=False
    rejected(lambda:owner.report(token));entry[1][key]=saved
    value=entry[2][3];saved=dict(value.terms);value.terms={};rejected(lambda:owner.report(token));value.terms=saved
    change=entry[3][4];saved=change[2][0]['measure'];change[2][0]['measure']='wrong'
    rejected(lambda:owner.report(token));change[2][0]['measure']=saved
    cell=change[3][0];value=cell[3][0];saved=dict(value.terms);value.terms={}
    rejected(lambda:owner.report(token));value.terms=saved
    rejected(lambda:owner.report(current.OriginalN1RadialCutoffPacket()))
    rejected(lambda:owner.evaluate('.371','1.01','first'))
    rejected(lambda:owner.changes('.371','.625',False))
    rejected(lambda:owner.changes('.371',c.mpf(['.6','.7'])))
    owner.report(token)
    checks.update(cutoff_scope_primitive_ledger_and_retained_range_mutations_rejected=True,
        foreign_packets_invalid_partitions_and_interval_endpoints_rejected=True)
    inputs=dict(owner.hashes);inputs[current.NAME]=current.sha(current.NAME);inputs[Path(__file__).name]=current.sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=inputs,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(op.snapshot(result),separators=(',',':'))+'\n').encode())
    print('PASS_N1_RADIAL_CUTOFF_EXTENSION',fixture,derivatives,flush=True);return result
