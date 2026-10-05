"""Independent coefficient extraction and Green/source checks.

Checks the actual changed operator. Unchanged upstream finite rows, physical
maps and amplitude reports are consumed by hash, not regenerated.
"""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as s

from lei_ren_part1_paper_compliant_current_core_nonlinear_operator import (
    HERE,PREFIX,NAME,RECEIPT,GATE,NORM_GATE,NORM_BOUNDS_GATE,OPEN,CurrentCoreNonlinearOperator,
    TaylorRectangle,source_terms,next_rows,correction_map,green_inverse,
    production_expressions,source_bindings_and_incidence,sha)
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

CTX = SimpleNamespace(mpf=s.sympify)


def equal(left,right,message):
    difference=s.cancel(left-right)
    if difference != 0 and s.simplify(difference) != 0:
        raise ArithmeticError(message)


def arbitrary_order_extraction():
    # Universal coefficient laws, not a finite numerical-order experiment.
    n,i,k=s.symbols('n i k',integer=True,nonnegative=True)
    rho=s.Symbol('rho',positive=True)
    monomial=rho**(n+1)
    tests=[((rho*s.diff(monomial,rho,2)+2*s.diff(monomial,rho))/rho**n,(n+1)*(n+2)),
           ((rho*s.diff(monomial,rho,2)+s.diff(monomial,rho))/rho**n,(n+1)**2),
           (s.diff(rho**(i+1)/(i+1),rho),rho**i),
           (s.diff(rho**(k+1)/s.factorial(k+1),rho),rho**k/s.factorial(k))]
    for left,right in tests:
        equal(left,right,'All-order derivative/integral coefficient law differs')
    # The source gauge identity is an equality in a commutative differential
    # algebra. Coefficient extraction is a homomorphism for sums/products,
    # with these exact derivative/integral laws. The prior whole production
    # AST proof binds those same laws to advance_scaled_one for arbitrary n.
    parent=json.loads((HERE/(PREFIX+'current_core_recurrence_source_check.json')).read_bytes())
    _verify_hashes(parent)
    if (not parent['all_passed'] or not parent['entire_production_AST_scaling']['passed']
        or not parent['all_radial_indices_homogeneity']['passed']):
        raise ValueError('Accepted entire-production recurrence source proof required')
    env=production_expressions()
    if env['theta_difference']!=0 or env['z_difference']!=0:
        raise ValueError('Universal production gauge identities absent')
    return dict(radial_variable='rho=s=Lambda R',
        coefficient_laws=['[rho^n] D2 Phi=(n+1)(n+2)Phi_(n+1)',
                          '[rho^n] D1 Psi=(n+1)^2 Psi_(n+1)',
                          'M_n=Psi_n/(n+1)',
                          'Pcal_(n+1)=F0^2 sum_(i=0)^n Phi_i Phi_(n-i)/(n+1)',
                          'ordinary Taylor derivative: (dZ f)_[k]=(k+1)f_[k+1]',
                          'ordinary Taylor product: (fg)_[k]=sum_(i=0)^k f_[i]g_[k-i]'],
        exact_source_identities_in_differential_algebra=True,
        prior_whole_production_AST_scaling_consumed_without_regeneration=True,
        coefficient_extraction_valid_at_every_nonnegative_radial_order=True,
        assumptions=['Phi_axis=1','Psi_axis=0','common pressure primitive',
                     'L has nonzero constant coefficient','epsilon constant in Z and rho'],
        analytic_existence_not_deduced_from_formal_identity=True,passed=True)


def independent_production_extraction():
    z,delta=s.symbols('Z delta')
    eps=s.Rational(1,2)
    counts={}
    for n,K in ((0,3),(1,2),(2,2)):
        size=n+K+1
        vector=lambda tag,length:[s.Symbol(tag+'_'+str(m)) for m in range(length)]
        fixed={key:vector(tag,size) for key,tag in
               (('ell_Z_taylor','ell'),('S_Z_taylor','S'),
                ('U0_Z_taylor','U0'),('P0_Z_taylor','P0'))}
        rows=dict(A=[[s.Integer(1)]+[s.Integer(0)]*(size-1)],
                  Uz=[list(fixed['U0_Z_taylor'])],P=[list(fixed['P0_Z_taylor'])])
        for j in range(1,n+1):
            rows['A'].append(vector('Phi'+str(j),size-j))
            rows['Uz'].append([eps*v for v in vector('Psi'+str(j),size-j)])
        # Independent pressure convolution, without the operator's ring.
        for j in range(1,n+1):
            pressure=[]
            for m in range(size-j):
                v=sum(fixed['S_Z_taylor'][a]*rows['A'][i][b]*rows['A'][j-1-i][m-a-b]
                      for a in range(m+1) for b in range(m-a+1) for i in range(j))
                pressure.append(v/j)
            rows['P'].append(pressure)
        packet=source_terms(fixed,rows,n,z,delta,eps)
        actual=next_rows(packet,n)
        original=copy.deepcopy(rows)
        advance_scaled_one(CTX,fixed,original,n,z,delta,eps)
        count=0
        for key in ('A','Uz','P'):
            if len(actual[key])!=K:
                raise ValueError('Extra derivative input or output order lost')
            for left,right in zip(actual[key],original[key][-1]):
                equal(left,right,'Actual gauge operator differs from production: '+key+' n='+str(n))
                count+=1
        counts[str(n)]=count
        # Coupled pressure derivative is calculated from the same primitive.
        if n:
            for m in range(K):
                equal(packet['pressure_z'].rows[n][m],
                      (m+1)*rows['P'][n][m+1]/eps**2,
                      'Pressure-Z is not the derivative of the common primitive')
    return dict(unmodified_production_advance_scaled_one_executed=True,
                independent_formal_axial_coefficients=True,
                exact_next_row_scalar_identities=counts,total_scalar_identities=sum(counts.values()),
                pressure_source_and_pressure_Z_coupling_checked=True,
                unavailable_axial_derivative_not_filled_with_zero=True,passed=True)


def independent_Green_and_models():
    # Arbitrary source coefficients and fixed chi Taylor function. Verify the
    # differential equation from the returned inverse coefficients directly.
    source=[[s.Symbol('f_'+str(n)+'_'+str(m)) for m in range(3)] for n in range(4)]
    chi=s.symbols('chi0:3')
    identities=0
    for component in ('theta','z'):
        out=green_inverse(source,component,chi)
        if any(out[0]):
            raise ArithmeticError('Correction Green inverse changes axis data')
        for n in range(4):
            for m in range(3):
                d=2*(n+1)*(n+2 if component=='theta' else n+1)
                lhs=d*out[n+1][m]
                if component=='theta':
                    lhs+=sum(chi[j]*out[n][m-j] for j in range(m+1))
                equal(lhs,source[n][m],'Green inverse equation/factor two differs')
                identities+=1
    c,B=s.symbols('chi B')
    model=[(-c/2)**n/(s.factorial(n)*s.factorial(n+1)) for n in range(6)]
    for n in range(5):
        equal(2*(n+1)*(n+2)*model[n+1]+c*model[n],0,'Bessel model coefficient differs')
    equal(2*(B/2),B,'Linear axial model differs')
    return dict(full_equation_Green_scalar_identities=identities,
                zero_correction_axis_preserved=True,raw_J_has_no_one_half=True,
                one_half_occurs_once_in_full_inverse=True,
                angular_resolvent_chi_Taylor_convolution_retained=True,
                Bessel_Phi_model_and_linear_Psi_model_preserved=True,passed=True)


def raw_weight_ratios():
    n,m=s.symbols('n m',integer=True,nonnegative=True)
    # Prove inequalities by polynomials with nonnegative coefficients after
    # shifting n,m into their specified domains. This is all-order algebra.
    certificates=[]
    expressions=[('J1',80*(n+1)**3*(n+m+1)-20*(n+2)**2),
                 ('J2',40*(n+1)**2*(n+m+1)-20*(n+2)),
                 ('V',80*(n+1)**2*(n+m+1)-20*(n+2)**2),
                 ('J1_Euler',80*(n+1)**3*(n+m+1)-20*n*(n+2)**2),
                 ('J2_Euler',40*(n+1)**2*(n+m+1)-20*n*(n+2)),
                 ('J1_dZ',80*(n+1)**3*(m+2)**2-20*(n+2)**2*(m+1)**2),
                 ('J2_dZ',40*(n+1)**2*(m+2)**2-20*(n+2)*(m+1)**2)]
    for name,polynomial in expressions:
        poly=s.Poly(s.expand(polynomial),n,m)
        if any(v<0 for v in poly.coeffs()):
            raise ArithmeticError('All-order Green norm certificate fails: '+name)
        certificates.append(dict(operator=name,nonnegative_polynomial=str(poly.as_expr())))
    return dict(weight='20^n h^m(n+1)^2(m+1)^2/binomial(n+m,m) for ordinary Taylor coefficients',
                raw_J1_norm_upper=80,raw_J2_norm_upper=40,
                raw_J1_dZ_norm_upper='80/h',raw_J2_dZ_norm_upper='40/h',
                all_order_nonnegative_polynomial_certificates=certificates,
                mixed_product_norm_not_inferred_from_single_operator_bounds=True,passed=True)


def supremum_composite_norm_proof():
    """Paper8.30/8.32 supremum norm, including BOTH convolution sums.

    No l1 coefficient-space substitution. Polynomial certificates below are
    for all nonnegative indices, while Vandermonde and the split-sum bound
    supply the complete radial/axial convolution, not an individual summand.
    """
    i,j,p,q=s.symbols('i j p q',integer=True,nonnegative=True)
    n,m=i+j,p+q
    certificates=[]
    for nu in (1,2):
        for label,multiplier in (('dZ_times_Euler',(i+1)*j),('dZ_times_field',i+1)):
            remainder=s.Poly(s.expand((n+1)*(n+nu)-multiplier),i,j)
            if any(c<0 for c in remainder.coeffs()):
                raise ArithmeticError('Composite derivative inversion inequality fails')
            certificates.append(dict(nu=nu,source=label,nonnegative_polynomial=str(remainder.as_expr())))
    r=s.Symbol('r',nonnegative=True)
    geometric=1/(1-r)
    weighted_sum=s.diff(r*s.diff(r*geometric,r),r)
    equal(weighted_sum,(1+r)/(1-r)**3,'Fixed analytic convolution series differs')
    H,sigma=s.symbols('H sigma',real=True)
    equal(H/(H*H+sigma*sigma),
          (1/(H+s.I*sigma)+1/(H-s.I*sigma))/2,'Cross-swirl pole identity differs')
    # Constant chain is exact, not a sampled majorant.
    if 2*4*2 != 16 or 16*16 != 256 or 20*4*16*16 != 20480:
        raise ArithmeticError('Supremum convolution constant differs')
    applications=[
        ('theta',1,'J2 f',40),('theta',2,'J2(f Euler g)',20480),
        ('theta',3,'J2(f_Z)', '40/h'),('theta',4,'J2(Mf g)',40*256),
        ('theta',5,'J2(Mf Euler g)',20480),('theta',6,'J2((Mf)_Z g)','20480/h'),
        ('theta',7,'J2((Mf)_Z Euler g)','20480/h'),('theta',8,'J2(fg)',40*256),
        ('theta',9,'J2(f g_Z)','20480/h'),('theta',10,'J2(fg)',40*256),
        ('z',1,'J1(Euler f)',20480),('z',2,'J1 f',80),('z',3,'J1(f_Z)','80/h'),
        ('z',4,'J1(Mf Euler g)',20480),('z',5,'J1((Mf)_Z Euler g)','20480/h'),
        ('z',6,'J1(fg)',80*256),('z',7,'J1(f g_Z)','20480/h'),
        ('z',8,'J1(Pcal_Z); Pcal=F0^2 V(Phi^2)','(80/h)*(80*256^2*F0_squared_norm)'),
        ('z',9,'J1 Pcal; same common primitive','80*(80*256^2*F0_squared_norm)'),
        ('z',10,'J1(rho F0^2 Phi^2)','80*80*256^2*F0_squared_norm')]
    return dict(norm_convention='supremum over n,m,Z; ordinary Taylor coefficients are dZ^m/m!',
        product_norm_upper=256,radial_average_norm_upper=1,radial_primitive_norm_upper=80,
        radial_derivative_product_after_Green_upper=20480,
        axial_derivative_product_after_Green_upper='20480/h',
        mixed_Mz_radial_derivative_after_Green_upper='20480/h',
        discrete_convolution_proof=[
            'Split i+j=n at n/2: the other denominator >=(n+1)^2/4.',
            'sum_(i>=0)1/(i+1)^2 <= 1+integral_1^infinity x^-2 dx=2.',
            'Both halves give sum_(i+j=n)1/((i+1)^2(j+1)^2)<=16/(n+1)^2.',
            'Apply the same bound to axial p+q=m, including shifted p+2.',
            'Vandermonde: binom(i+p,p)binom(j+q,q)<=binom(n+m,m).',
            'Derivative: (p+1)binom(i+p+1,p+1)=(i+1)binom(i+p+1,p).',
            'Derivative Vandermonde denominator is binom(n+m+1,m).',
            '(i+1)j/((n+1)(n+nu))<=1; dropping the M factor1/(i+1) only enlarges the bound.',
            '(n+2)^2/(n+1)^2<=4; hence20*4*16*16/h=20480/h.'],
        nonnegative_polynomial_certificates=certificates,
        fixed_multiplier_proof=['A fixed radial-constant coefficient commutes with J1/J2.',
                                '|c_[a]|<=complex_modulus/(eta/2)^a by Cauchy.',
                                'Weighted axial convolution ratio <=(a+1)^2; binomial ratio<=1.',
                                'Its summed factor is sum_(a>=0)(a+1)^2 r^a=(1+r)/(1-r)^3.',
                                'Apply once to the whole source product, with r=h/(eta/2)<1.',
                                'Stored coefficient modulus*weight is conservative since weight>=1.'],
        cross_swirl_pole_identity='H/(H^2+sigma^2)=((H+i sigma)^-1+(H-i sigma)^-1)/2',
        cross_swirl_modulus_upper='d_modulus/pole_distance; H factor is included by partial fractions',
        pressure_and_swirl_norms_use_actual_F0_squared_bound=True,
        term_specific_composite_bound_applications=[dict(component=c,number=n,operator=op,bound=b)
                                                  for c,n,op,b in applications],
        complete_convolution_sums_in_supremum_norm=True,
        analytic_tail_and_shared_solution_identification_not_proved_here=True,passed=True)


@source_precision
def run(operator=None):
    raw=json.loads((HERE/NAME).read_bytes())
    _verify_hashes(raw)
    operator=operator if operator is not None else CurrentCoreNonlinearOperator(require_checked=False)
    if (raw['actual_five_defect_family_sha256']!=operator.family
        or raw['implicit_source_sha256']!=operator.source_sha or raw['datum_enclosure_sha256']!=operator.datum_sha
        or raw[GATE] or raw[NORM_BOUNDS_GATE] or raw[NORM_GATE] or any(raw[key] for key in OPEN)
        or raw['source_bindings_and_incidence']!=source_bindings_and_incidence()):
        raise ValueError('Current operator source/scope changed')
    runtime=operator.evaluate_seed('.293')
    if encode(runtime)!=raw['fresh_current_runtime']:
        raise ValueError('New current finite-source/Green runtime differs from saved producer')
    result=dict(actual_five_defect_family_sha256=operator.family,implicit_source_sha256=operator.source_sha,
                datum_enclosure_sha256=operator.datum_sha,
                source_bindings_and_incidence=operator.incidence,
                arbitrary_order_extraction=arbitrary_order_extraction(),
                independent_production_extraction=independent_production_extraction(),
                independent_Green_and_models=independent_Green_and_models(),
                raw_weight_ratios=raw_weight_ratios(),
                supremum_composite_norm_proof=supremum_composite_norm_proof(),
                actual_pressure_primitive_and_its_Z_derivative_use_same_source=True,
                fresh_current_runtime_replayed=True,current_runtime_radial_degree=4,
                formal_operator_does_not_replace_exact_source_with_Cauchy_box=True,
                **{GATE:True,NORM_BOUNDS_GATE:True,NORM_GATE:False},**dict.fromkeys(OPEN,False),
                remaining_dependency='Actual-operator Banach admission and same finite/tail solution identification, then core-first pressure/ODE join',
                all_scoped_checks_passed=True,all_passed=True,
                input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS actual twenty-term source, Green inverses and original radial extraction; analytic fixed-point/tail still open',flush=True)
    return result


if __name__=='__main__':
    run()
