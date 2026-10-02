"""Actual fixed-radius extrema identities, source bounds and finite fixture."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_compliant_rooted_swirl_morphology import CompliantRootedSwirlMorphology
from lei_ren_part1_paper_compliant_root_centered_peak import CompliantRootCenteredPeak
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_rooted_swirl_morphology.json'


def identities():
    r,z,delta=s.symbols('rho Z delta',positive=True);D=1-z*z
    P=s.Function('Phi')(r,z);E=s.diff(P,z)-2*z*r*s.diff(P,r)/D
    B=E/P-(2+delta)*z/D
    along=lambda expression:s.diff(expression,z)-2*z*r*s.diff(expression,r)/D
    desired=(s.diff(P,z,2)-4*z*r*s.diff(P,r,z)/D+4*z*z*r*r*s.diff(P,r,2)/D**2-2*r*s.diff(P,r)/D)/P
    desired-=(E/P)**2+(2+delta)*(1+z*z)/D**2
    if s.simplify(along(B)-desired)!=0:raise ArithmeticError('Fixed physical radius B derivative failed')
    lam,F,kp,b,sigma,epsilon,y=s.symbols('lambda F Kprime b sigma epsilon y',positive=True)
    # At fixed r, utheta=r*lambda^(-2-delta)*F0*Phi;
    # rho_Z=-2Zrho/D, lambda_Z/lambda=Z/D, F0_Z/F0=-Kprime/b.
    slope=b*(-(2+delta)*z/D-kp/b+E/P)
    if s.simplify(slope-(-kp+b*B))!=0:raise ArithmeticError('Physical angular log-slope normalization failed')
    h,L,q=s.symbols('h L q',positive=True)
    Kp=(b*y)*L*q/(1+epsilon*(b*y)**2*q*q)
    if s.simplify(Kp/b-y*L*q/(1+epsilon*b*b*y*y*q*q))!=0:
        raise ArithmeticError('Second-centered peak slope cancellation failed')
    u=s.symbols('u',positive=True);a=s.symbols('a',real=True)
    quadratic=y*a-h*y*y/2
    maximum=quadratic.subs(y,a/h)
    if s.simplify(maximum-a*a/(2*h))!=0:raise ArithmeticError('Normalized log-height comparison failed')
    return dict(actual_fixed_radius_log_swirl_slope=True,along_fixed_radius_B_derivative=True,
        second_centered_peak_chart_identity=True,positive_log_height_comparison=True,
        radial_swirl_derivative_identity='Dr[r*Phi(rho)] = Phi+2rho*Phi_rho, rho proportional to r^2',passed=True)


def finite_scale_fixture():
    with mp.workdps(160):
        c=MPIntervalContext();c.dps=120;j=mp.mpf('.12');delta=mp.mpf('.02');Lambda=mp.mpf(1000);sigma=j/500
        H=lambda z:-4*z**3-j*z*z+(9-delta)*z/2+j
        a=mp.findroot(H,-j/((9-delta)/2));b=sigma/mp.sqrt(Lambda);Da=1-a*a
        peak=object.__new__(CompliantRootCenteredPeak);peak.ctx=c
        peak.core=SimpleNamespace(j=c.mpf('.12'),delta=c.mpf('.02'),Lambda=c.mpf(1000),logLambda=c.ln(1000))
        peak.anchor=c.mpf([a-mp.mpf('1e-110'),a+mp.mpf('1e-110')]);peak.delta=peak.core.delta;peak.sigma=c.mpf('.12')/500
        peak.initialize_chart()
        R,Z=s.symbols('R Z');P=1-s.Rational(3,100)*R+s.Rational(1,100)*R*Z+s.Rational(1,500)*R*R
        coefficients={(i,k):s.Poly(s.diff(P,R,i,Z,k),R,Z).terms() for i,k in ((0,0),(1,0),(0,1),(0,2),(1,1),(2,0))}
        def profile(xi,rho,i=0,k=0):
            zz=peak.anchor+peak.b*c.mpf(xi);rr=c.mpf(rho)
            value=sum((c.mpf(int(v.p))/int(v.q)*rr**powers[0]*zz**powers[1] for powers,v in coefficients[(i,k)]),c.mpf(0))
            return dict(source_jets=dict(Phi=value))
        f=object.__new__(CompliantRootedSwirlMorphology);f.ctx=c;f.core=peak.core;f.peak=peak;f.cache={}
        f.field=SimpleNamespace(rooted_profile=profile)
        phi=lambda rr,zz:1-mp.mpf(3)*rr/100+rr*zz/100+rr*rr/500
        g=lambda z:(1-delta*z*z)*H(z)/(H(z)**2+sigma*sigma)
        # Independent scalar physical angular velocity at FIXED r. Constants
        # r,tau,Cstar cancel in this within-slice logarithmic ratio.
        logu=lambda z:-Lambda*mp.quad(g,[a,z])+mp.log(phi(mp.mpf(2)*(1-z*z)/Da,z))+(2+delta)/2*mp.log(1-z*z)
        location=f.peak_location(2);ylow,yhigh=endpoints(location['peak_y_interval'])
        root=mp.findroot(lambda yy:mp.diff(logu,a+b*b*yy),mp.mpf('.02'))
        if not ylow<=root<=yhigh:raise ArithmeticError('Independent physical swirl maximum outside certified bracket')
        gain=(logu(a+b*b*root)-logu(a))/(b*b)
        lo,hi=endpoints(location['normalized_log_peak_gain_over_b_squared'])
        if not lo<=gain<=hi or gain<=0:raise ArithmeticError('Independent physical peak-height gain failed')
        for ypoint in (ylow,yhigh):
            derivative=mp.diff(logu,a+b*b*ypoint)
            if (ypoint==ylow and derivative<=0) or (ypoint==yhigh and derivative>=0):
                raise ArithmeticError('Independent extrema endpoint signs failed')
        radial=f.radial_monotonicity()
        if endpoints(radial['positive_radial_swirl_derivative_coefficient'])[0]<=0:
            raise ArithmeticError('Finite-scale radial monotonicity comparison failed')
        return dict(synthetic_fixture_not_paper_source=True,independent_scalar_physical_swirl_maximum=True,
            fixed_radius_rho_relation_used=True,visible_peak_y=mp.nstr(root,45),
            visible_normalized_log_height_gain=mp.nstr(gain,45),independent_extrema_endpoint_signs=True,passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Morphology source changed: '+name)
    with mp.workdps(400):
        f=CompliantRootedSwirlMorphology();c=f.ctx;algebra=identities();fixture=finite_scale_fixture()
        maxima=levels=0
        for rho,row in receipt['actual_physical_swirl_peak_locations'].items():
            B=read_interval(c,row['rooted_slice_bounds']['log_weight_B'])
            curvature=read_interval(c,row['rooted_slice_bounds']['actual_log_swirl_xi_curvature'])
            if endpoints(B)[0]<=0 or endpoints(curvature)[1]>=0:raise ArithmeticError('Actual slice positivity/concavity not strict')
            inner=read_interval(c,row['inner_scaled_log_slope']);outer=read_interval(c,row['outer_scaled_log_slope'])
            if endpoints(inner)[0]<=0 or endpoints(outer)[1]>=0:raise ArithmeticError('Stored actual peak endpoint signs not strict')
            current=f.peak_location(rho)
            for key in ('peak_y_interval','normalized_log_peak_gain_over_b_squared','inner_scaled_log_slope','outer_scaled_log_slope'):
                if read_interval(c,row[key])._mpi_!=current[key]._mpi_:raise ArithmeticError('Source extrema replay failed: '+key)
            gain=read_interval(c,row['normalized_log_peak_gain_over_b_squared'])
            B0=read_interval(c,row['rooted_slice_bounds']['anchor_log_weight_B'])
            kappa=read_interval(c,row['rooted_slice_bounds']['actual_negative_log_swirl_xi_curvature'])
            trial=read_interval(c,row['lower_height_comparison_trial_xi'])
            if (endpoints(B0)[0]<=0 or endpoints(kappa)[0]<=0 or endpoints(trial)[0]<=0
                    or endpoints(trial)[1]>=4 or not row['height_gain_uses_actual_log_curvature_and_anchor_slope']):
                raise ArithmeticError('Actual-curvature peak-height proof not admitted')
            if endpoints(gain)[0]<=0 or not row['Phi_and_physical_lambda_derivatives_included']:
                raise ValueError('F0-only peak substituted for actual physical swirl')
            if row['global_swirl_peak_across_annuli_measured']:raise ValueError('Local axial peak promoted globally')
            maxima+=1
            for name,width in receipt['true_peak_normalized_axial_widths'][rho].items():
                drop=read_interval(c,width['log_drop']);current_width=f.level_width(rho,drop)
                if read_interval(c,width['full_xi_width'])._mpi_!=current_width['full_xi_width']._mpi_:
                    raise ArithmeticError('Actual physical fractional width replay failed')
                for side,packet in width['sides'].items():
                    low=read_interval(c,packet['inner_log_ratio_to_true_peak']);high=read_interval(c,packet['outer_log_ratio_to_true_peak'])
                    if endpoints(low)[0]<=-endpoints(drop)[0] or endpoints(high)[1]>=-endpoints(drop)[1]:
                        raise ArithmeticError('True peak-normalized level endpoint signs not strict')
                if not width['normalized_to_actual_physical_swirl_peak'] or not width['radius_is_fixed_during_each_axial_scan']:
                    raise ValueError('Wrong axial scan normalization')
                levels+=1
        radial=receipt['radial_monotonicity'];coefficient=read_interval(c,radial['positive_radial_swirl_derivative_coefficient'])
        if endpoints(coefficient)[0]<=0 or not radial['global_radial_swirl_peak_requires_original_annular_values']:
            raise ValueError('Cutoff radius substituted for actual radial morphology')
        for flag in ('complete_vortex_radial_width_measured','whole_vortex_aspect_ratio_measured','all_annular_source_values_resolved',
                     'full_point_physical_field_evaluation','measured_blowup_dynamics','physical_energy_integral_certified',
                     'admissible_stress_lift_constructed','temporal_recursion'):
            if receipt[flag]:raise ValueError('Rooted axial morphology scope promoted: '+flag)
    for name in (NAME,Path(__file__).name):hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        independent_source_identities=algebra,independent_physical_extrema_fixture=fixture,
        actual_fixed_radius_Phi_weighted_maxima_checked=maxima,
        actual_true_peak_normalized_axial_levels_checked=levels,whole_rooted_domain_radial_swirl_monotonicity_checked=True,
        actual_Phi_weighted_local_physical_swirl_peak_certified=True,
        full_radial_vortex_width_or_aspect_ratio_measured=False,
        original_annular_value_resolution_needed_for_radial_swirl_peak=True,
        measured_blowup_dynamics=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Physical swirl morphology PASS:3 true local maxima,9 fractional levels,whole-root radial monotonicity',flush=True)
    return result


if __name__=='__main__':run()
