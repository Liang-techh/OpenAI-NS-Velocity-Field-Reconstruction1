"""Fresh analytic core tube for the Section 9/10-selected inlet tolerance.

Pressure data is unchanged, but j, sigma, analytic radii, G bound and the
core-family fingerprint are fresh. No old finite core is admitted here.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent


def run():
    name = 'lei_ren_part1_paper_shared_bump_constants.json'
    raw = json.loads((HERE / name).read_bytes())
    for path, digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE / path).read_bytes()).hexdigest() != digest:
            raise ValueError('Constants dependency changed: ' + path)
    if not raw['fixed_CA_CQ_directed_bounds_certified'] or not raw['fixed_CS_analytic_bound_certified_under_stated_inputs']:
        raise ValueError('Fixed constants required')
    c = MPIntervalContext()
    c.dps = 160
    with mp.workdps(c.dps + 40):
        j = restore_value(c, raw['required_j'])
        dt = c.mpf('1e-200')
        sigma = j / 500
        eta = sigma / 40
        h = eta / 8
        a = 1 + eta
        Hder = (9+dt)/2 + 2*j*a + 12*a*a
        displacement = eta*Hder
        pole = sigma-displacement
        L = 1-dt*a*a
        q = 1-eta*(2+eta)
        if not (endpoints(Hder)[1] < 20 and min(endpoints(pole)[0],endpoints(L)[0],endpoints(q)[0]) > 0):
            raise ArithmeticError('Fresh common tube exclusion failed')
        chi = 1+sigma*sigma/pole**2
        beta = (3+dt/2+j*a)/L
        H = (9+dt)*a/2+j+4*a**3+j*a*a
        gradient = (1+dt*a*a)*H/pole**2
        Greal = (1+j)*(1+dt)/(2*sigma)
        G = Greal+eta*gradient
        Gbar = c.mpf(endpoints(G)[1])
        definition = dict(j_definition='min(epsilon0/4,e_star/100)/8',
            CA=raw['CA'],CQ=raw['CQ'],CS=raw['CS'],KN=raw['KN'],
            sigma_definition='j/500',eta_definition='sigma/40',h_definition='eta/8',
            delta_upper='1e-200',core_U0='4Z+j',core_H0='(1-delta)Z/2+(1-Z^2)(4Z+j)',
            G_anchor='unique H0 root in [-j,0]',
            logC='Lambda*Gbar+2logLambda+1000',
            constants_receipt_sha256=hashlib.sha256((HERE/name).read_bytes()).hexdigest())
        family_sha = hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        report = dict(j=j,delta=dt,sigma=sigma,complex_tube_radius=eta,Xh_parameter=h,
            H0_derivative_modulus_upper=Hder,H0_variation_upper=displacement,
            denominator_factor_modulus_lower=pole,L_modulus_lower=L,
            pressure_q_modulus_lower=q,chi_modulus_upper=chi,beta_modulus_upper=beta,
            H0_modulus_upper=H,complex_gradient_modulus_upper=gradient,
            G_modulus_upper=G,A_Omega_upper=G,Gbar_for_symbolic_logC=Gbar,
            common_complex_axis_poles_excluded=True,
            tube_definition='union of disks |Z-x|<=eta over real x in [-1,1]',
            norm_definition='paper8.30:20^n h^m(n+1)^2(m+1)^2|d_Z^m f_n|/[m! binomial(n+m,m)]',
            core_family_definition=definition,analytic_core_family_sha256=family_sha,
            source_j_eta_tol_relation_verified=True,
            Cstar_guard_by_symbolic_definition=True,Cstar_margin=1000,
            actual_Cstar_or_F0_not_materialized=True,
            analytic_P0_modulus_certified=False,nonlinear_contraction_certified=False,
            fresh_finite_core_generated=False,old_finite_coefficients_reused=False,
            full_Section9_parameter_admission=False,temporal_recursion=False,
            input_hashes={**raw['input_hashes'],name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('Fresh shared-tolerance tube: j',mp.nstr(endpoints(j)[1],14),
              'eta',mp.nstr(endpoints(eta)[1],14),'poles excluded')
        return report


if __name__ == '__main__':
    run()
