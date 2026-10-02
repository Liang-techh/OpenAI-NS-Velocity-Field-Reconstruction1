"""Independent weighted integrals, differentiated product tails and source replay."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_compliant_core_integral_atoms import (
    CompliantCoreIntegralAtoms,finite_atom_coefficients,product_tail)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_core_integral_atoms.json'


def independent_integral_fixture():
    with mp.workdps(240):
        c=MPIntervalContext();c.dps=180;rho,z,t=s.symbols('rho Z t');center=s.Rational(3,10);order=6
        phi=1-rho/12+rho*rho/150+z*rho/90+z*z*rho*rho/700
        uz=4*z+s.Rational(1,10)+s.Rational(7,100)*rho+z*rho*rho/250+z**3*rho/1000
        rows=[]
        for expression in (phi,uz):
            polynomial=s.Poly(expression.subs(z,center+t),rho,t)
            rows.append([[c.mpf(int(polynomial.coeff_monomial(rho**n*t**k).p))
                          /int(polynomial.coeff_monomial(rho**n*t**k).q) for k in range(order+1)] for n in range(3)])
        values=finite_atom_coefficients(c,*rows,order)
        actual_phi=phi+rho**3*(1+z+z*z)/10**8
        actual_uz=uz-2*rho**4*(1-z+z**3)/10**8
        integrands=dict(H=rho*phi/8,M=uz/4,K=rho*phi*uz/8,A=uz**2/4,B=rho*phi**2/16,C=phi**2/4)
        actual=dict(H=rho*actual_phi/8,M=actual_uz/4,K=rho*actual_phi*actual_uz/8,
                    A=actual_uz**2/4,B=rho*actual_phi**2/16,C=actual_phi**2/4)
        pnorm=[sum((abs(row[k])*c.mpf(4)**n for n,row in enumerate(rows[0])),c.mpf(0)) for k in range(order+1)]
        unorm=[sum((abs(row[k])*c.mpf(4)**n for n,row in enumerate(rows[1])),c.mpf(0)) for k in range(order+1)]
        tails=[]
        for expression in (actual_phi-phi,actual_uz-uz):
            coeffs=[]
            for k in range(order+1):
                v=s.diff(expression,z,k).subs({rho:4,z:center})/math.factorial(k)
                coeffs.append(abs(c.mpf(int(v.p))/int(v.q)))
            tails.append(coeffs)
        ep,eu=tails
        errors=dict(H=ep,M=eu,K=[product_tail(c,pnorm,unorm,ep,eu,k) for k in range(order+1)],
                    A=[product_tail(c,unorm,unorm,eu,eu,k) for k in range(order+1)],
                    B=[product_tail(c,pnorm,pnorm,ep,ep,k)/2 for k in range(order+1)],
                    C=[product_tail(c,pnorm,pnorm,ep,ep,k) for k in range(order+1)])
        count=tail_count=0
        for name,expression in integrands.items():
            integrated=s.integrate(expression,(rho,0,4));difference=s.integrate(actual[name]-expression,(rho,0,4))
            for k in range(order+1):
                exact=s.diff(integrated,z,k).subs(z,center)/math.factorial(k)
                target=mp.mpf(int(exact.p))/int(exact.q)
                if not endpoints(values[name][k])[0]<=target<=endpoints(values[name][k])[1]:
                    raise ArithmeticError('Independent weighted atom integral failed: '+str((name,k)))
                change=s.diff(difference,z,k).subs(z,center)/math.factorial(k)
                if abs(mp.mpf(int(change.p))/int(change.q))>endpoints(errors[name][k])[1]:
                    raise ArithmeticError('Differentiated product integral tail failed: '+str((name,k)))
                count+=1;tail_count+=1
        return dict(independent_exact_weighted_integral_coefficients=count,
            independent_signed_nonlinear_product_tail_coefficients=tail_count,
            factorial_and_affine_U0_units_checked=True,synthetic_fixture_not_paper_source=True,passed=True)


def independent_factorial_tails():
    with mp.workdps(160):
        c=MPIntervalContext();c.dps=180
        f=object.__new__(CompliantCoreIntegralAtoms);f.ctx=c;f.axial_weight=c.mpf('.00000001')
        f.core=SimpleNamespace(delta=c.mpf('.02'),j=c.mpf('.12'),sigma=c.mpf('.12')/500)
        j=mp.mpf('.12');delta=mp.mpf('.02');sigma=j/500;center=mp.mpf('.3')
        H=lambda z:-4*z**3-j*z*z+(9-delta)*z/2+j
        chi=lambda z:H(z)**2/(H(z)**2+sigma*sigma)
        coefficients=mp.taylor(chi,center,6);counts=0
        # Powers are independently convolved scalar Taylor rows. Check signed
        # tail coefficients against the uniform modulus bound, through order6.
        powers=[mp.mpf(1)]+[mp.mpf(0)]*6
        sums={N:[mp.mpf(0)]*7 for N in (6,24)}
        for n in range(1,145):
            powers=[sum((powers[ell]*coefficients[k-ell] for ell in range(k+1)),mp.mpf(0)) for k in range(7)]
            factor=mp.mpf(4)**n/(mp.mpf(2)**n*math.factorial(n)*math.factorial(n+1))
            for N in sums:
                if n>N:
                    for k in range(7):sums[N][k]+=abs(powers[k])*factor
        for N in sums:
            bounds=f.model_tail_coefficients(c.mpf('.3'),N,6)
            for k in range(7):
                if sums[N][k]>endpoints(bounds[k])[1]:raise ArithmeticError('Axial6 factorial model tail too small')
                counts+=1
        return dict(independent_factorial_radial_tail_sums_through_axial6=counts,passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Atom source changed: '+name)
    fixture=independent_integral_fixture();tails=independent_factorial_tails()
    with mp.workdps(400):
        f=CompliantCoreIntegralAtoms();c=f.ctx;atom_count=endpoint_count=0
        for center,packed in receipt['core_atom_packets'].items():
            fresh=f.root_atoms() if center=='exact_shared_root' else f.core_atoms(center)
            if packed['cover_endpoints_used_as_atom_values'] or packed['point_parameter_representatives_selected'] or packed['old_finite_rows_read']:
                raise ValueError('Actual coefficientwise atoms replaced by cover/old rows/point choices')
            for name,coefficients in packed['actual_core_atom_axial_coefficients'].items():
                for k,coefficient in enumerate(coefficients):
                    old=read_interval(c,coefficient);new=fresh['actual_core_atom_axial_coefficients'][name][k]
                    if old._mpi_!=new._mpi_:raise ArithmeticError('Actual atom source replay failed: '+str((center,name,k)))
                    if read_interval(c,packed['directed_integrated_tail_coefficient_bounds'][name][k])._mpi_!=fresh['directed_integrated_tail_coefficient_bounds'][name][k]._mpi_:
                        raise ArithmeticError('Atom tail replay failed')
                    if read_interval(c,packed['actual_core_atom_ordinary_axial_derivatives'][name][k])._mpi_!=(new*math.factorial(k))._mpi_:
                        raise ArithmeticError('Axial derivative/Taylor factorial mismatch')
                    atom_count+=1
            for name in ('H','B','C'):
                if endpoints(fresh['actual_core_atom_axial_coefficients'][name][0])[0]<=0:raise ArithmeticError('Positive atom lost')
            for k in range(7):
                C=fresh['actual_core_atom_axial_coefficients']['C'][k]
                if read_interval(c,packed['pressure_primitive_at_exit_axial_coefficients'][k])._mpi_!=(4*C)._mpi_:
                    raise ArithmeticError('Original pressure primitive normalization failed')
                M=fresh['actual_core_atom_axial_coefficients']['M'][k]
                if read_interval(c,packed['core_exit_profile_axial_coefficients']['mean'][k])._mpi_!=M._mpi_:
                    raise ArithmeticError('Inlet mean is not same primitive')
                endpoint_count+=1
            if center=='exact_shared_root':
                if not packed['exact_shared_H_root_used_before_enclosure']:raise ValueError('Lost exact root identity')
                if any(endpoints(read_interval(c,v))!=(mp.mpf(0),mp.mpf(0)) for v in packed['radial_tail_coefficient_bounds']['Phi_model']):
                    raise ArithmeticError('Root model radial tail must vanish through axial6')
        for flag in ('original_signed_bridge_integrals_resolved','all_annular_source_values_resolved',
                     'full_point_physical_field_evaluation','measured_blowup_dynamics','physical_energy_integral_certified',
                     'admissible_stress_lift_constructed','temporal_recursion'):
            if receipt[flag]:raise ValueError('Core atom scope promoted: '+flag)
    for name in (NAME,Path(__file__).name):hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        independent_exact_integral_fixture=fixture,independent_factorial_tail_fixture=tails,
        actual_integrated_atom_axial_coefficients_checked=atom_count,
        shared_mean_and_pressure_primitive_coefficients_checked=endpoint_count,
        actual_core_atoms_coefficientwise_integrated=True,actual_core_atoms_axial6_available=True,
        original_signed_bridge_integrals_resolved=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual core atoms PASS:210 source coefficients,42 independent integrals,42 signed product tails,axial6',flush=True)
    return result


if __name__=='__main__':run()
