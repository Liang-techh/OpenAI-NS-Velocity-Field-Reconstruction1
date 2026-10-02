"""Independent identities, real scaled integrals, peak levels and units."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_compliant_root_centered_peak import CompliantRootCenteredPeak
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_root_centered_peak.json'


def contains(bound,point,label):
    lo,hi=endpoints(bound)
    if not lo<=point<=hi:raise ArithmeticError('Independent value outside source enclosure: '+label)


def identities():
    a,b,xi,delta,h,A,epsilon=s.symbols('a b xi delta h A epsilon',real=True)
    La=1-delta*a*a;q=h+A*b*xi-4*b*b*xi*xi;L=1-delta*(a+b*xi)**2
    exact=s.Poly(s.expand(L*q),xi)
    expected=[La*h,b*(La*A-2*delta*a*h),b*b*(-4*La-2*delta*a*A-delta*h),
              b**3*(8*delta*a-delta*A),4*delta*b**4]
    for k,row in enumerate(expected):
        if s.expand(exact.nth(k)-row)!=0:raise ArithmeticError('Scaled exact numerator identity failed')
    antiderivative=sum(row*xi**(k+2)/(k+2) for k,row in enumerate(expected))
    if s.expand(s.diff(antiderivative,xi)-xi*L*q)!=0:
        raise ArithmeticError('Signed polynomial primitive failed')
    v=s.symbols('v',real=True)
    # H(a+x)=x*q(x), sigma^2=Lambda*b^2 and epsilon=1/Lambda.
    scaled=xi*L*q/(1+epsilon*xi*xi*q*q)
    original=1/epsilon*L*(b*xi*q)/(b*b/epsilon+b*b*xi*xi*q*q)*b
    if s.cancel(original-scaled)!=0:raise ArithmeticError('Original to root-scaled gradient identity failed')
    omitted=xi*L*q-scaled
    if s.cancel(omitted-epsilon*xi**3*L*q**3/(1+epsilon*xi*xi*q*q))!=0:
        raise ArithmeticError('Positive denominator remainder identity failed')
    Z=s.symbols('Z',real=True);tau,D=s.symbols('tau D',positive=True)
    p=(1-delta)/2
    # D=1-Z^2>0 on the finite physical chart. Use its exact chain rule
    # rather than asking a CAS to split fractional powers without a domain.
    physical_z=Z*tau**p*D**(-p)
    derivative=s.diff(physical_z,Z)-2*Z*s.diff(physical_z,D)
    coefficient=s.simplify(derivative/(tau**p*D**(-p-1)))
    if s.simplify(coefficient.subs(D,1-Z**2)-(1-delta*Z**2))!=0:
        raise ArithmeticError('Original physical axial map Jacobian failed')
    lam,rho,F,Phi,C,K=s.symbols('lambda rho F Phi C K',positive=True)
    swirl=lam**(-1-delta)*s.sqrt(2*rho*epsilon)*s.exp(-C)*s.exp(-K)*Phi
    factors=Phi*s.exp(-K)*s.exp(-C)*s.sqrt(2*rho*epsilon)/lam*lam**(-delta)
    if s.simplify(swirl-factors)!=0:raise ArithmeticError('Split logarithmic physical swirl units failed')
    if s.expand((1-(a+b*xi)**2)-(1-a*a)-(-2*a*b*xi-b*b*xi*xi))!=0:
        raise ArithmeticError('Physical lambda ratio increment identity failed')
    phi0,change=s.symbols('phi0 change',nonzero=True)
    if s.cancel((phi0+change)/phi0-1-change/phi0)!=0:
        raise ArithmeticError('Shared-source Phi relative normalization identity failed')
    return dict(five_exact_numerator_coefficients=True,signed_polynomial_primitive=True,
        original_scaled_gradient=True,positive_denominator_remainder_identity=True,
        physical_axial_map_Jacobian=True,split_logarithmic_swirl_units=True,
        physical_lambda_ratio_increment_identity=True,shared_source_Phi_normalization_identity=True,
        physical_ratio_mean_value_bound='q in[1/2,1], t>=1/2: q*t^(q-1)<=sqrt(2)<2',passed=True)


def original_source_diagnostics(f,receipt):
    nominal=lambda value:sum(endpoints(value))/2
    j=nominal(f.core.j);delta=nominal(f.delta);epsilon=nominal(f.epsilon);sigma=j/500
    H=lambda z:-4*z**3-j*z*z+(9-delta)*z/2+j
    anchor=mp.findroot(H,-j/((9-delta)/2),tol=mp.eps*16)
    hp=(9-delta)/2-12*anchor**2-2*j*anchor;A=-12*anchor-j;b=sigma*mp.sqrt(epsilon)
    # Evaluate factored polynomials; no subprecision b*xi is added to a.
    q=lambda x:hp+A*b*x-4*b*b*x*x
    L=lambda x:1-delta*anchor**2-2*delta*anchor*b*x-delta*b*b*x*x
    grad=lambda x:x*L(x)*q(x)/(1+epsilon*x*x*q(x)**2)
    X,aa,bb,dd,hh,AA,ee=s.symbols('X aa bb dd hh AA ee',real=True)
    qq=hh+AA*bb*X-4*bb*bb*X*X
    LL=1-dd*aa*aa-2*dd*aa*bb*X-dd*bb*bb*X*X
    exact_grad=X*LL*qq/(1+ee*X*X*qq*qq)
    differentiated=[s.lambdify((X,aa,bb,dd,hh,AA,ee),s.diff(exact_grad,X,k),'mpmath') for k in range(5)]
    integrals=derivatives=amplitude_derivatives=0;rows=[]
    for text,packet in receipt['peak_packets'].items():
        x=mp.mpf(text);K=mp.quad(grad,[0,x])
        bound=read_interval(f.ctx,packet['actual_scaled_G']);contains(bound,K,'actual source K '+text)
        contains(read_interval(f.ctx,packet['normalized_F0_relative_to_shared_anchor']),mp.exp(-K),'actual source exp(-K) '+text)
        for k,row in enumerate(packet['K_derivatives']):
            # Scalar finite differences cannot resolve derivatives as tiny
            # as b or epsilon. Symbolic differentiation retains these atoms.
            contains(read_interval(f.ctx,row),differentiated[k](x,anchor,b,delta,hp,A,epsilon),'actual source gradient '+str((text,k)))
            derivatives+=1
        # Independent derivative via scalar Taylor ODE A'=-K' A.
        g=[differentiated[k](x,anchor,b,delta,hp,A,epsilon)/math.factorial(k) for k in range(5)];coef=[mp.exp(-K)]
        for n in range(1,6):coef.append(-sum(g[k]*coef[n-1-k] for k in range(n))/n)
        for n,row in enumerate(packet['normalized_F0_derivatives']):
            contains(read_interval(f.ctx,row),coef[n]*math.factorial(n),'normalized amplitude derivative '+str((text,n)))
            amplitude_derivatives+=1
        if text!='0' and endpoints(bound)[1]-endpoints(bound)[0]>mp.mpf('1e-180'):
            raise ArithmeticError('Scaled K peak is not resolved independently of huge Lambda')
        rows.append(dict(xi=text,K=mp.nstr(K,45),normalized_F0=mp.nstr(mp.exp(-K),45)))
        integrals+=1
    return dict(independent_original_scaled_real_integrals=integrals,independent_K_derivatives=derivatives,
        independent_normalized_amplitude_derivatives=amplitude_derivatives,rows=rows,
        quadrature_is_diagnostic_not_certificate=True,provisional_parameters_not_selected_source=True,
        tiny_denominator_effect_proved_algebraically_not_resolved_by_scalar_quadrature=True,
        independent_gradient_differentiation_keeps_b_and_epsilon_atoms=True,passed=True)


def moderate_scale_fixture():
    """Diagnostic with a visible denominator correction, never a paper datum."""
    c=MPIntervalContext();c.dps=120
    j=mp.mpf('.02');delta=mp.mpf('.001');sigma=j/500;Lam=mp.mpf(25)
    H=lambda z:-4*z**3-j*z*z+(9-delta)*z/2+j
    a=mp.findroot(H,-j/((9-delta)/2),tol=mp.eps*16)
    f=object.__new__(CompliantRootCenteredPeak);f.ctx=c
    f.core=SimpleNamespace(j=c.mpf('.02'),Lambda=c.mpf(25),logLambda=c.ln(25),logC=c.mpf(1000))
    f.anchor=c.mpf([a-mp.mpf('1e-110'),a+mp.mpf('1e-110')]);f.delta=c.mpf('.001');f.sigma=c.mpf('.02')/500
    f.initialize_chart();b=sigma/mp.sqrt(Lam)
    # Independent source H evaluated in ordinary coordinates is safe here.
    gradient=lambda x:Lam*(1-delta*(a+b*x)**2)*H(a+b*x)/(H(a+b*x)**2+sigma*sigma)*b
    count=0;max_visible=mp.mpf(0)
    for text in ('-4','-2','-.5','.5','2','4'):
        x=mp.mpf(text);K=mp.quad(gradient,[0,x]);packet=f.evaluate(text)
        contains(packet['actual_scaled_G'],K,'moderate fixture integral '+text)
        polynomial=packet['integrated_numerator_polynomial'];tail=packet['denominator_remainder_abs_upper']
        if K>endpoints(polynomial)[1] or endpoints(polynomial)[0]-K>endpoints(tail)[1]:
            raise ArithmeticError('One-sided denominator remainder failed at visible correction scale')
        max_visible=max(max_visible,abs(sum(endpoints(polynomial))/2-K))
        for k,row in enumerate(packet['K_derivatives']):
            contains(row,mp.diff(gradient,x,k),'moderate fixture derivative '+str((text,k)));count+=1
    if max_visible<1:raise ArithmeticError('Fixture did not expose omitted denominator correction')
    return dict(synthetic_fixture_not_paper_source=True,Lambda='25',j='.02',delta='.001',
        independent_visible_denominator_integrals=6,independent_derivatives=count,
        largest_visible_denominator_effect=mp.nstr(max_visible,25),passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Peak prerequisite changed: '+name)
    with mp.workdps(400):
        f=CompliantRootCenteredPeak();c=f.ctx;algebra=identities()
        if receipt['implicit_source_sha256']!=f.core.source or receipt['actual_five_defect_family_sha256']!=f.core.family:
            raise ValueError('Peak chart must use the accepted compliant pressure/core source')
        if 'lei_ren_part1_paper_compliant_core_uniform_bounds.py' not in hashes:
            raise ValueError('Full analytic Phi derivative embedding must be source-hash bound')
        chart=receipt['scaled_chart']
        if endpoints(read_interval(c,chart['positive_b']))[0]<=0:
            raise ArithmeticError('Original microscopic width was zeroed')
        for name in ('whole_chart_positive_q_lower','whole_chart_positive_L_lower'):
            if endpoints(read_interval(c,chart[name]))[0]<=0:raise ArithmeticError('Lost source-uniform monotone peak')
        original=original_source_diagnostics(f,receipt);fixture=moderate_scale_fixture()
        widths=[]
        for name,row in receipt['level_widths'].items():
            drop=read_interval(c,row['log_drop'])
            for side,packet in row['sides'].items():
                inner=read_interval(c,packet['inner_endpoint_K']);outer=read_interval(c,packet['outer_endpoint_K'])
                if not endpoints(inner)[1]<endpoints(drop)[0] or not endpoints(outer)[0]>endpoints(drop)[1]:
                    raise ArithmeticError('Stored peak level signs are not strict')
                root=read_interval(c,packet['absolute_xi_root_interval'])
                # Independent nominal level root of exact scaled real integral.
                j=sum(endpoints(f.core.j))/2;d=sum(endpoints(f.delta))/2;e=sum(endpoints(f.epsilon))/2
                H=lambda z:-4*z**3-j*z*z+(9-d)*z/2+j
                a=mp.findroot(H,-j/((9-d)/2));b=j/500*mp.sqrt(e);h=(9-d)/2-12*a*a-2*j*a;A=-12*a-j
                g=lambda x:x*(1-d*a*a-2*d*a*b*x-d*b*b*x*x)*(h+A*b*x-4*b*b*x*x)/(1+e*x*x*(h+A*b*x-4*b*b*x*x)**2)
                target=sum(endpoints(drop))/2;sign=packet['sign']
                nominal=mp.findroot(lambda w:mp.quad(g,[0,sign*w])-target,mp.sqrt(2*target/((1-d*a*a)*h)))
                contains(root,nominal,'peak fractional level '+name+side)
            if row['whole_vortex_aspect_ratio_measured']:raise ValueError('Core cutoff mislabeled as measured vortex width')
            widths.append(name)
        for geometry in receipt['physical_geometry_packets']:
            t=read_interval(c,geometry['log_tau']);relative=-f.delta*t/2
            if read_interval(c,geometry['relative_elongation_log_factor'])._mpi_!=relative._mpi_:
                raise ArithmeticError('Original positive delta/time correlation was lost')
            if not geometry['coordinate_scaling_observation_not_fitted_dynamics']:
                raise ValueError('Imposed time exponent promoted to measured dynamics')
        contains(read_interval(c,receipt['physical_geometry_packets'][-1]['relative_elongation_log_factor']),mp.mpf(1),'slow time log factor')
        swirls=0
        for xi,packets in receipt['root_centered_physical_swirl_packets'].items():
            for row in packets:
                fresh=f.swirl_value(xi,read_interval(c,row['rho']))
                for key in ('anchor_Phi','transported_Phi_enclosure','physical_utheta_relative_to_shared_anchor','physical_utheta_coefficient'):
                    if read_interval(c,row[key])._mpi_!=fresh[key]._mpi_:raise ArithmeticError('Shared-root swirl source replay failed: '+key)
                for stored,current in zip(row['physical_utheta_positive_scale_log_terms'],fresh['physical_utheta_positive_scale_log_terms']):
                    if read_interval(c,stored)._mpi_!=current._mpi_:raise ArithmeticError('Physical swirl logarithmic scale changed')
                if xi=='0' and endpoints(read_interval(c,row['physical_utheta_relative_to_shared_anchor']))!=(mp.mpf(1),mp.mpf(1)):
                    raise ArithmeticError('Shared-source Phi normalization does not cancel exactly')
                if not row['same_nonlinear_Phi_from_fresh_radial_rows'] or not row['shared_root_transport_uses_full_analytic_norm']:
                    raise ValueError('Nonlinear Phi replaced with an independent model')
                swirls+=1
        for packet in receipt['peak_packets'].values():
            if len(packet['absolute_log_F0_as_separate_terms'])!=2 or not packet['root_offset_not_added_to_anchor']:
                raise ValueError('Finite peak shape was absorbed into the enormous base amplitude')
        for flag in ('full_swirl_peak_including_Phi_measured','whole_vortex_aspect_ratio_measured',
                     'full_point_physical_field_evaluation','all_annular_source_values_resolved',
                     'measured_blowup_dynamics','physical_energy_integral_certified','admissible_stress_lift_constructed','temporal_recursion'):
            if receipt[flag]:raise ValueError('Root-centered peak scope promoted: '+flag)
    for name in (NAME,Path(__file__).name):hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        source_identities=algebra,original_source_diagnostics=original,visible_denominator_fixture=fixture,
        certified_fractional_levels_checked=widths,source_uniform_monotone_left_right_level_roots=6,
        reproduced_same_source_physical_swirl_packets=swirls,original_delta_slow_time_factor_preserved=True,
        actual_scaled_G_peak_resolved=True,absolute_F0_not_materialized=True,
        full_swirl_peak_including_Phi_measured=False,whole_vortex_aspect_ratio_measured=False,
        measured_blowup_dynamics=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Root-centered peak PASS:9 source integrals,6 level roots,6 Phi-weighted swirl packets,visible denominator fixture',flush=True)
    return result


if __name__=='__main__':run()
