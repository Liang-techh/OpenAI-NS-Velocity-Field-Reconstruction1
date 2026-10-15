"""Evaluate the original full n1 core with principal/feedback error budgets."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_n1_core_volterra_numeric import (
    OriginalN1CoreVolterraNumeric,DifferenceJet,forcing_difference)

DIRECTORY = ROOT/'experiments/root_st073'
NAME = 'lei_ren_part1_paper_compliant_current_original_n1_core_volterra_numeric.json'


def exact_forcing_difference(solver):
    """Check the error map against physical, unnormalized weighted sources."""
    r,z = s.symbols('r z')
    de,e,P2 = s.Rational(1,17),s.Rational(1,3),s.Integer(5)
    A = 1+z*r/7+r*r/11
    U = 2+z/3+r*z*z/5
    da = r*(1+z+z*z)/13+r**3*z/17
    du = r*r*(1+z)/19
    mean = lambda v:s.cancel(s.integrate(v,r)/r)
    R = solver.R
    def local(value):
        rows = [[s.diff(value,r,i,z,j).subs({r:s.Rational(1,2),z:0})/
            (math.factorial(i)*math.factorial(j)) for j in range(4)] for i in range(3)]
        return R(rows,2,3,s.sympify)
    pair = lambda value,error=0:DifferenceJet(local(value),local(error))
    g = e*(2+2*z)
    result = forcing_difference(pair(A,da),pair(U,du),pair(mean(U),mean(du)),
        pair(z),pair(r),pair(de),pair(e),pair(P2),pair(g))
    F0 = s.exp(-z*z-2*z)
    d,L = 1-z*z,1-de*z*z
    T = lambda weight,value:(weight*z*value+d*s.diff(value,z)-2*z*r*s.diff(value,r))/L
    def original(a,u):
        q = s.cancel((2*z*u-(1-de)*z*mean(u)-d*s.diff(mean(u),z))/L)
        S = (q+(1-de)*z*s.diff(q,z)/2+r*s.diff(q,r))/L
        S += q*(q/2+r*s.diff(q,r))+u*T(-2,q)
        np = -S/2+(2*s.diff(q,r)+r*s.diff(q,r,2))/e
        return dict(f=-2*e*e*T(-3,T(-2-de,F0*a))/F0,
            u=-2*T(-2,T(-1-de,u))-4*z*r*e*np/L,p=np/P2)
    old,new = original(A,U),original(A+da,U+du)
    count = 0
    for field,value in result.items():
        for k in (0,1):
            at = lambda expr:s.simplify(s.diff(expr,z,k).subs({r:s.Rational(1,2),z:0}))
            assert s.simplify(value.base.rows[0][k]-at(old[field])) == 0,(field,k,'base')
            assert s.simplify(value.error.rows[0][k]-at(new[field]-old[field])) == 0,(field,k,'difference')
            count += 2
    return count


def exact_principal_equations(solver):
    """Check integrated polynomials by their differential equations."""
    R = solver.R
    r,z = s.symbols('r z')
    source = {key:R(rows,2,1,s.sympify) for key,rows in dict(
        f=[[s.Rational(2,3),s.Rational(1,7)],[s.Rational(-1,5),s.Rational(2,11)],[s.Rational(1,13),0]],
        u=[[s.Rational(1,3),s.Rational(-2,7)],[s.Rational(2,5),0],[s.Rational(-1,11),s.Rational(1,17)]],
        p=[[s.Rational(3,5),0],[s.Rational(-1,7),s.Rational(1,13)],[s.Rational(2,17),0]]).items()}
    chi = R([[s.Rational(1,3),s.Rational(1,7)]],0,1,s.sympify)
    J = 8
    result = solver._principal(source,chi,J)
    expr = lambda v:sum(v.rows[n][k]*r**n*z**k for n in range(v.radial+1) for k in range(2))
    f,u,k,p = (expr(result[key]) for key in ('f','u','k','p'))
    hf,hu,hp = (expr(source[key]) for key in ('f','u','p'))
    residual = s.expand(4*(r*s.diff(f,r,2)+2*s.diff(f,r))+2*expr(chi)*f-hf)
    count = 0
    for j in (0,1):
        for i in range(J+1):
            assert residual.coeff(r,i).coeff(z,j) == 0
            count += 1
        for body in (4*(r*s.diff(u,r,2)+s.diff(u,r))-hu,s.diff(p,r)-hp,
                     r*s.diff(k,r)+k+r*s.diff(u,r)):
            assert s.expand(body).coeff(z,j) == 0
            count += 1
    assert all(expr(value).subs(r,0) == 0 for value in result.values())
    # This verifies actual Volterra integration, not a fitted point table.
    return count


def solve(solver,degree,kernel):
    values = []
    for Z in ('-0.3','0','0.3'):
        for rho in ('0','1','2','4'):
            values.append(solver.evaluate(Z,rho,degree=degree,kernel_degree=kernel))
        print(f"Original full n1 core: dps={solver.precision}, Z={Z}, source degree={degree}, Kchi degree={kernel}",flush=True)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args = parser.parse_args()
    if args.inspect_report:
        report = json.loads((DIRECTORY/NAME).read_text(encoding='utf-8'))
        for path,digest in report['implementation_hashes'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest,path
        for name,digest in report['input_hashes'].items():
            assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest,name
        assert report['checks']['passed']
        print(json.dumps(dict(checks=report['checks'],archive_only=True),indent=2))
        return
    started = time.monotonic()
    coarse = OriginalN1CoreVolterraNumeric(precision=50,quadrature_nodes=96)
    exact = dict(exact_source_difference_and_derivative_comparisons=exact_forcing_difference(coarse),
        exact_principal_equation_and_kinematic_comparisons=exact_principal_equations(coarse))
    print('Independent exact source-difference and principal-equation checks pass',flush=True)
    first = solve(coarse,12,16)
    fine = OriginalN1CoreVolterraNumeric(precision=70,quadrature_nodes=192)
    second = solve(fine,18,24)
    c = mp.mp.clone();c.dps = 120
    comparisons = []
    maximum_relative = {key:c.mpf(0) for key in ('f','u','k','p','fx','ux','mean_u','q')}
    nonzero_increments = 0
    for old,new in zip(first,second):
        assert (old['Z'],old['rho']) == (new['Z'],new['rho'])
        for field,view in new['normalized_original_n1'].items():
            previous = old['normalized_original_n1'][field]
            rows = enumerate(zip(view,previous)) if isinstance(view,list) else [(0,(view,previous))]
            for derivative,(current,before) in rows:
                point = c.mpf(current['approximation'])
                difference = abs(point-c.mpf(before['approximation']))
                combined = c.mpf(current['absolute_error_upper'])+c.mpf(before['absolute_error_upper'])
                assert difference <= combined,(new['Z'],new['rho'],field,derivative,difference,combined)
                if derivative == 0 and point:
                    relative = c.mpf(current['relative_error_upper'])
                    assert relative < c.mpf('1e-16'),(new['Z'],new['rho'],field,relative)
                    maximum_relative[field] = max(maximum_relative[field],relative)
                comparisons.append(dict(Z=new['Z'],rho=new['rho'],field=field,Z_derivative=derivative,
                    difference=c.nstr(difference,35),combined_original_full_n1_error_upper=c.nstr(combined,35)))
        for increment in new.get('nominal_feedback_increments',[]):
            nonzero_increments += sum(c.mpf(v) != 0 for row in increment.values() for v in row)
    assert nonzero_increments > 0
    # Whole-core APIs retain their exact axis and reject unsupported domains.
    rejected = 0
    for Z,rho,kwargs in (('1','2',{}),('0.3','4.01',{}),('0.3','-1',{}),
                         ('0.3','1',{'degree':3}),('0.3','1',{'feedback_iterations':2})):
        try:fine.evaluate(Z,rho,**kwargs)
        except ValueError:rejected += 1
        else:raise AssertionError((Z,rho,kwargs))
    implementation = dict(fine.code_hashes)
    for path in (Path(__file__),ROOT/'src/openai_ns_reconstruction/original_n1_core_volterra_numeric.py'):
        implementation[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    report = dict(source_family=fine.family,checks=dict(passed=True,**exact,
        cross_precision_original_full_n1_comparisons=len(comparisons),
        original_core_points_per_level=len(second),nonzero_nominal_feedback_components=nonzero_increments,
        rejected_unsupported_requests=rejected,
        maximum_refined_value_relative_error_upper={k:c.nstr(v,30) for k,v in maximum_relative.items()},
        accepted_Kchi_decomposition_and_natural_component_frame=True,
        feedback_word_tail_applies_to_principal_origin=True,
        leading_source_tail_potential_kernel_tail_and_feedback_word_tail_separate=True,
        no_source_midpoint_selected=True,n1_radial_Taylor_tail_substituted=False),
        levels=[dict(precision=50,pressure_quadrature_nodes=96,source_degree=12,kernel_degree=16,samples=first),
                dict(precision=70,pressure_quadrature_nodes=192,source_degree=18,kernel_degree=24,samples=second)],
        comparisons=comparisons,accepted_regular_word_tail=fine.tail,input_hashes=fine.hashes,
        implementation_hashes=implementation,
        limitations=dict(original_full_core_n1_numerical_approximation_with_absolute_error=True,
            exact_original_point_oracle=False,full_inner_collar_not_evaluated=True,
            nominal_feedback_iteration_certified_separately=False,
            nominal_feedback_charged_against_index_zero_tail=True,
            common_global_N=False,global_moment_repair_evaluated=False,
            n0_plus_n1_physical_velocity_assembled=False,heat_stress_layers_completed=False,
            genuine_n_ge_2_temporal_recursion=False,full_corrected_NS_solution=False),
        execution_seconds=round(time.monotonic()-started,3))
    if args.write_report:(DIRECTORY/NAME).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)


if __name__ == '__main__':main()
