"""Compute original n1 finite radial prefixes and independent equation checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_n1_core_numeric import (
    OriginalN1CoreNumeric, RAW_SOURCE, linear_rhs, pure_definitions)

DIRECTORY = ROOT/"experiments/root_st073"
NAME = "lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json"


def raw_conjugation():
    """Independently contract the authoritative six-row raw matrices."""
    raw = pure_definitions(DIRECTORY/RAW_SOURCE,("raw_matrix_expressions",),
        {"sy":s})["raw_matrix_expressions"]()
    z,r,x,de,e,P2,F0,A,Ar,Az,U,Ur,Uz,Q,g = s.symbols(
        "Z rho x delta epsilon P2 F0 A Ar Az U Ur Uz Q g", nonzero=True)
    f,u,k,p,fr,ur,fz,uz,kz,pz = s.symbols("f u k p fr ur fz uz kz pz")
    d,L = 1-z*z,1-de*z*z
    source = dict(rho=r,x=x,Z=z,delta=de,epsilon=e,F=F0*A,Fr=F0*Ar,
        FZ=F0*(Az-g*A/e),U=U,Ur=Ur,UZ=Uz,Q=Q)
    state = (F0*f/e,e*u,e*k,e*P2*p,2*x*F0*fr/e,2*x*e*ur)
    zstate = (F0*(fz-g*f/e)/e,e*uz,e*kz,e*P2*pz,0,0)
    action = {i:s.Integer(0) for i in (4,5,6)}
    for name,values in (("B0",state),("B1",zstate)):
        for key,expression in raw[name].items():
            i,j = map(int,key.split(","))
            if i in action:
                actual = expression.xreplace({symbol:source[str(symbol)] for symbol in expression.free_symbols})
                action[i] += actual*values[j-1]
    normalized = (s.expand(action[5]*e/F0).subs(x*x,r),
        s.expand(action[6]/e).subs(x*x,r),s.cancel(action[4]/(2*x*e*P2)))
    values = dict(z=z,rho=r,delta=de,epsilon=e,P2=P2,F2=F0**2,L=L,d=d,A=A,C=A+r*Ar,g=g,
        E=z*U-s.Rational(1,2),H=z*(1-de)/2+d*U,B=Q+(1-2*z*U)/L,Q=Q,
        Zf_without_gauge=((-2-de)*z*A+d*Az-2*z*r*Ar)/L,
        Zw=((-1-de)*z*U+d*Uz-2*z*r*Ur)/L,
        f=f,u=u,k=k,p=p,fZ=fz,uZ=uz,kZ=kz,pZ=pz,Euler_f=r*fr,Euler_u=r*ur,Euler_U=r*Ur)
    got = linear_rhs(values)
    for a,b in zip(got,normalized):
        assert s.cancel(a-b) == 0
    # Each coefficient is checked separately, so pressure, gauge and radial
    # transport cannot hide behind one hand-picked numerical state.
    state_symbols = (f,u,k,p,fr,ur,fz,uz,kz,pz)
    pairs = 0
    for a,b in zip(got,normalized):
        for variable in state_symbols:
            assert s.cancel(s.diff(a-b,variable)) == 0
            pairs += 1
    symbols = dict(z=z,r=r,de=de,e=e,P2=P2,F0=F0,A=A,Ar=Ar,Az=Az,U=U,Ur=Ur,Uz=Uz,Q=Q,g=g,
        f=f,u=u,k=k,p=p,fr=fr,ur=ur,fz=fz,uz=uz,kz=kz,pz=pz)
    return normalized,symbols,pairs


def exact_mixed_fixture(solver,normalized,symbols):
    """Compare the executable recurrence to direct raw coefficient extraction.

    A nonconstant nonzero F0 tests its gauge derivative twice and the
    F0-squared pressure/axial feedback. This fixture is independent of the
    enormous physical parameter hierarchy; actual sources are solved below.
    """
    z,r = symbols["z"],symbols["r"]
    e,de,P2 = s.Rational(1,3),s.Rational(1,17),s.Integer(5)
    F0 = s.exp(-2*z-z*z)
    g = e*(2+2*z)
    A = 1+(1+z*z)*r/7+z*r*r/11+r**3/19+r**4/23
    U = 2+z/3+z*r/5+(1+z)*r*r/13+r**3/29+r**4/31
    d,L = 1-z*z,1-de*z*z
    mean = s.cancel(s.integrate(U,r)/r)
    Q = s.cancel((2*z*U-(1-de)*z*mean-d*s.diff(mean,z))/L)
    T = lambda weight,v:(weight*z*v+d*s.diff(v,z)-2*z*r*s.diff(v,r))/L
    Nt = -T(-3,T(-2-de,F0*A))
    Nz = -T(-2,T(-1-de,U))
    S = (Q+(1-de)*z*s.diff(Q,z)/2+r*s.diff(Q,r))/L
    S += Q*(Q/2+r*s.diff(Q,r))+U*T(-2,Q)
    Np = -S/2+(2*s.diff(Q,r)+r*s.diff(Q,r,2))/e
    forcing = (2*e*e*Nt/F0,2*Nz-4*z*r*e*Np/L,Np/P2)
    constants = {symbols[key]:value for key,value in dict(e=e,de=de,P2=P2,F0=F0,A=A,
        Ar=s.diff(A,r),Az=s.diff(A,z),U=U,Ur=s.diff(U,r),Uz=s.diff(U,z),Q=Q,g=g).items()}
    degree,depth = 3,1
    count = degree+depth+3
    coefficient = lambda expr,i,j:s.cancel(s.diff(expr,r,i,z,j).subs({r:0,z:0})/(math.factorial(i)*math.factorial(j)))
    leading = {key:[[coefficient(expression,i,j) for j in range(count-i)]
        for i in range(degree+2)] for key,expression in (("A",A),("Uz",U))}
    fixed = dict(ell_Z_taylor=[s.diff(-g,z,j).subs(z,0)/math.factorial(j) for j in range(count)],
        S_Z_taylor=[s.diff(e*e*F0**2,z,j).subs(z,0)/math.factorial(j) for j in range(count)])
    state = {key:[[s.Integer(0)]*(degree+depth+1)] for key in ("f","u","k","p")}
    c = SimpleNamespace(mpf=s.sympify)
    for n in range(degree):
        solver._advance(fixed,leading,state,n,degree,depth,s.Integer(0),de,e,P2,c)
    polynomials = {key:s.Integer(0) for key in state}
    checked = 0
    feedback_nonzero = False
    for order in range(1,degree+1):
        substitutions = dict(constants)
        substitutions.update({symbols[key]:value for key,value in polynomials.items()})
        substitutions.update({symbols[key+"z"]:s.diff(polynomials[key],z) for key in polynomials})
        substitutions.update({symbols["fr"]:s.diff(polynomials["f"],r),
                              symbols["ur"]:s.diff(polynomials["u"],r)})
        for field,index,denominator in (("f",0,4*order*(order+1)),("u",1,4*order*order),("p",2,order)):
            body = normalized[index].subs(substitutions,simultaneous=True)+forcing[index]
            row = s.cancel(s.diff(body,r,order-1).subs(r,0)/(math.factorial(order-1)*denominator))
            polynomials[field] += row*r**order
            for j,value in enumerate(state[field][order]):
                expected = s.simplify(s.diff(row,z,j).subs(z,0)/math.factorial(j))
                assert s.simplify(value-expected) == 0,(field,order,j,value,expected)
                checked += 1
        polynomials["k"] += -s.Rational(order,order+1)*s.expand(polynomials["u"]).coeff(r,order)*r**order
        for j,value in enumerate(state["k"][order]):
            expected = -s.Rational(order,order+1)*state["u"][order][j]
            assert value == expected
            checked += 1
        if order == 3:
            # Raw B6,1 is indispensable at rho^3, even though its value
            # vanishes at Z=0. Its first Z derivative is nonzero there.
            term = -8*z*r*F0**2*A*polynomials["f"]/L
            test = s.diff(term,r,2,z).subs({r:0,z:0})/2
            assert test != 0
            feedback_nonzero = True
    return dict(exact_raw_equation_mixed_coefficient_comparisons=checked,
        gauge_and_axial_second_source_checked=True,nonzero_F0_squared_feedback_checked=feedback_nonzero)


def solve(solver,degree):
    results,coefficients = [],[]
    for Z in ("-0.3","0","0.3"):
        packet = solver.coefficients(Z,degree=degree)
        for rho in ("0","1","2","4"):
            results.append(solver.evaluate_prefix(packet,rho))
        coefficients.append(packet)
        print(f"Computed coupled original n1 coefficients: dps={solver.precision}, Z={Z}, degree={degree}",flush=True)
    return results,coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report",action="store_true")
    parser.add_argument("--inspect-report",action="store_true")
    args = parser.parse_args()
    if args.inspect_report:
        report = json.loads((DIRECTORY/NAME).read_text(encoding="utf-8"))
        for relative,digest in report["implementation_hashes"].items():
            assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() == digest,relative
        for name,digest in report["input_hashes"].items():
            assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest,name
        assert report["checks"]["passed"]
        print(json.dumps(dict(checks=report["checks"],archive_only=True),indent=2))
        return
    started = time.monotonic()
    normalized,symbols,pairs = raw_conjugation()
    print("Exact conjugation agrees with all original raw matrix entries",flush=True)
    coarse = OriginalN1CoreNumeric(precision=50,quadrature_nodes=96)
    exact = exact_mixed_fixture(coarse,normalized,symbols)
    print("Exact mixed fixture agrees with raw equation coefficient extraction",flush=True)
    first,old_packets = solve(coarse,16)
    fine = OriginalN1CoreNumeric(precision=70,quadrature_nodes=192)
    second,packets = solve(fine,16)
    c = mp.mp.clone()
    c.dps = 120
    comparisons = []
    for old,new in zip(first,second):
        assert (old["Z"],old["rho"]) == (new["Z"],new["rho"])
        for field,value in new["normalized_prefix"].items():
            previous = old["normalized_prefix"][field]
            difference = abs(c.mpf(value["approximation"])-c.mpf(previous["approximation"]))
            joint = c.mpf(value["finite_prefix_absolute_error_upper"])+c.mpf(previous["finite_prefix_absolute_error_upper"])
            assert difference <= joint,(new["Z"],new["rho"],field)
            comparisons.append(dict(Z=new["Z"],rho=new["rho"],field=field,
                difference=c.nstr(difference,35),combined_finite_prefix_error_upper=c.nstr(joint,35)))
    # Refinement differences are observed diagnostics, not a radial tail proof.
    refinements = []
    for packet in packets:
        for rho in ("1","2","4"):
            views = [fine.evaluate_prefix(packet,rho,degree=N) for N in (8,12,16)]
            for field in views[-1]["normalized_prefix"]:
                points = [c.mpf(v["normalized_prefix"][field]["approximation"]) for v in views]
                differences = [abs(points[i+1]-points[i]) for i in range(2)]
                refinements.append(dict(Z=packet["Z"],rho=rho,field=field,
                    degree8_to12_difference=c.nstr(differences[0],30),
                    degree12_to16_difference=c.nstr(differences[1],30),
                    observed_refinement_only=True,radial_remainder_proved=False))
    coefficient_rows = []
    nonzero_rows = 0
    for packet in packets:
        for field,rows in packet["point"].items():
            for n,row in enumerate(rows):
                for j,value in enumerate(row):
                    if value != 0:
                        nonzero_rows += 1
                    coefficient_rows.append(dict(Z=packet["Z"],field=field,rho_degree=n,Z_degree=j,
                        approximation=fine.ctx.nstr(value,fine.precision),
                        source_enclosure_exact_mpf_tuples=[list(v) for v in packet["directed"][field][n][j]._mpi_]))
    implementation = dict(fine.code_hashes)
    for path in (Path(__file__),ROOT/"src/openai_ns_reconstruction/original_n1_core_numeric.py"):
        implementation[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    report = dict(source_family=fine.family,checks=dict(passed=True,
        raw_matrix_normalization_coefficient_comparisons=pairs,**exact,
        cross_precision_finite_prefix_comparisons=len(comparisons),
        directed_coefficient_slots=len(coefficient_rows),nonzero_point_coefficient_slots=nonzero_rows,
        radial_refinement_diagnostics=len(refinements),
        same_original_leading_n1_source_family=True,no_source_midpoint_selected=True,
        all_B0_B1_pressure_feedback_and_Z_derivatives_retained=True),
        levels=[dict(precision=50,quadrature_nodes=96,degree=16,samples=first),
                dict(precision=70,quadrature_nodes=192,degree=16,samples=second)],
        coefficients=coefficient_rows,comparisons=comparisons,radial_refinements=refinements,
        input_hashes=fine.hashes,implementation_hashes=implementation,
        limitations=dict(finite_radial_prefix_only=True,original_n1_coefficient_enclosures=True,
            nonzero_implicit_swirl_retained_in_directed_branch=True,
            numerical_point_proxy_F0_squared_zero_is_not_an_original_source_choice=True,
            pressure_atom_proxy_error_covered_by_original_directed_source=True,
            radial_Taylor_remainder_certified=False,Volterra_word_tail_reused_as_radial_tail=False,
            original_n1_full_point_oracle=False,original_n1_full_point_error=False,
            global_physical_velocity_evaluated=False,temporal_n_ge_2_recursion=False,
            full_corrected_NS_solution=False),execution_seconds=round(time.monotonic()-started,3))
    if args.write_report:
        (DIRECTORY/NAME).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(dict(checks=report["checks"],execution_seconds=report["execution_seconds"]),indent=2),flush=True)


if __name__ == "__main__":
    main()
