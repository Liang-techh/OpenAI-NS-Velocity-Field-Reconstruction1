"""Independent finite-parameter integral fixture and admitted-view checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_absolute_moment_closure import source_identities
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_absolute_moment_closure.json'


def finite_parameter_fixture():
    """True Gamma/collar and squared bump densities, without interval caps.

    Moderate parameters exercise the physical unit conversion. This fixture
    does not replace admission of the actual Md40 source. Infinite pressure
    tails are bounded after an explicit finite endpoint; angular tails use
    the exact integrated Gamma expectation.
    """
    with mp.workdps(70):
        mu, a, eps, Ts = map(mp.mpf, ('.025', '.1', '.0002', '.7'))
        rate, k, bp, bh, prate = 1-mu, 1-a, mp.mpf('.5')+mu, mp.mpf('.5')+a, 1+2*mu
        phrate = 1+2*a

        def sigma(t):
            if t <= 0: return mp.mpf(0)
            if t >= 1: return mp.mpf(1)
            A, B = mp.exp(-1/t**2), mp.exp(-1/(1-t)**2)
            return A/(A+B)

        def shape(t):
            phi = mp.exp(-4/(3-t)**2) if t < 3 else mp.mpf(0)
            sig = sigma(t)
            return 1-eps*(1-sig+sig*phi), sig*(1-eps*phi)

        J = mp.quad(lambda t: mp.exp(k*t)*(1-shape(t)[0])/eps, [0, '.5', 1, 2, 3])
        Xt0 = mp.mpf(5)
        tau = (mp.log(Xt0-1/k)+mp.log1p(-eps)-mp.log(eps)-mp.log(1/k+J))/k
        length = 2+Ts+tau
        ell = -mp.mpf('1.5')*(2+Ts)+(rate+k)/2-bh*tau
        Rrel = mp.mpf(3); Rtail = Rrel*mp.exp(length)
        Etail = mp.exp(ell); cinf = Etail*Rtail**bh/(1-eps)
        transfer = (Rrel/Rtail)**mp.mpf('1.5')/Etail
        tail_scale = cinf**2*Rtail**(-phrate)
        rawbeta = lambda x: mp.exp(-1/(1-x*x)) if abs(x) < 1 else mp.mpf(0)
        norm = mp.quad(rawbeta, [-1, 0, 1]); width = mp.mpf('.15')
        beta = lambda v: rawbeta(v/width)/(width*norm)
        A = mp.quad(lambda v: mp.exp(rate*v)*beta(v), [-width, 0, width])
        B = mp.quad(lambda v: mp.exp(-prate*v)*beta(v), [-width, 0, width])
        D = mp.quad(lambda v: mp.exp(-prate*v)*beta(v)**2, [-width, 0, width])
        breaks = list(map(mp.mpf, ('-4', '-3.15', '-3', '-2.85', '-1.15', '-1', '-.85', '0')))
        rawpast = (mp.exp(4*prate)-1)/(2*prate)
        rows = []; max_error = mp.mpf(0)
        for Z in map(mp.mpf, ('0', '.5', '1')):
            xi = 2*(1-Z**2)/Rtail

            def H(t):
                x = xi*mp.exp(-t)
                return x**(-1-a)*mp.hyperu(1+a, 2, 1/x) if x else mp.mpf(1)

            def actual(t):
                pre, Csig = shape(t)
                return pre-Csig*(1-H(t))

            collar_A = mp.quad(lambda t: mp.exp(k*t)*(shape(t)[0]-actual(t)), [0, '.5', 1, 2, 3])
            xi3 = xi*mp.exp(-3)
            gamma_A3 = xi3**(-1-a)*mp.hyperu(1+a, 3, 1/xi3)/k if xi3 else 1/k
            delta_A = collar_A+mp.exp(3*k)*(gamma_A3-1/k)
            preP = mp.quad(lambda t: mp.exp(-phrate*t)*shape(t)[0]**2/2, [0, '.5', 1, 2, 3])+mp.exp(-3*phrate)/(2*phrate)
            limit = mp.mpf(45)
            actualP = mp.quad(lambda t: mp.exp(-phrate*t)*actual(t)**2/2, [0, '.5', 1, 2, 3])
            actualP += mp.quad(lambda t: mp.exp(-phrate*t)*H(t)**2/2, [3, 8, 20, limit])
            actualP += mp.exp(-phrate*limit)/(2*phrate)
            # 0 <= 1-H <= a*(1+a)*xi*exp(-t), so the replacement
            # of the remaining pure-heat integral has this directed error.
            pressure_remainder = a*(1+a)*xi*mp.exp(-(phrate+1)*limit)/(phrate+1)
            rpre = mp.mpf('.004')*Z**2
            # Convert by actual radii/amplitudes, independently of the
            # production logarithmic multiplier expressions.
            rH = cinf*Rtail**k*delta_A/Rrel**mp.mpf('1.5')
            sH = tail_scale*(preP-actualP)
            r = rpre+rH
            angular = lambda d1, d2: A*(mp.exp(-3*rate)*d1+mp.exp(-rate)*d2)
            pressure = lambda d1, d2: B*(mp.exp(3*prate)*d1+mp.exp(prate)*d2)+D*(mp.exp(3*prate)*d1**2+mp.exp(prate)*d2**2)/2
            lin = mp.matrix([[A*mp.exp(-3*rate), A*mp.exp(-rate)], [B*mp.exp(3*prate), B*mp.exp(prate)]])
            guess = mp.lu_solve(lin, mp.matrix([r, sH]))
            d1, d2 = mp.findroot(lambda x, y: (angular(x, y)-r, pressure(x, y)-sH), tuple(guess), tol=mp.mpf('1e-60'))
            h = lambda v: d1*beta(v+3)+d2*beta(v+1)
            measured_A = mp.quad(lambda v: mp.exp(rate*v)*h(v), breaks)
            actualpast = mp.quad(lambda v: mp.exp(-prate*v)*(1+h(v))**2/2, breaks)
            datumP0 = -(mp.mpf('2.5')+rawpast+1+tail_scale*preP)
            actualMp = mp.mpf('2.5')+actualpast+1+tail_scale*actualP
            physical_pressure_constant = datumP0+actualMp
            XtZ = Xt0-mp.exp(-(rate+k)/2)*rpre
            corrected_Xtail = 1/k+(XtZ+mp.exp(-(rate+k)/2)*measured_A-1/k)*mp.exp(-k*tau)
            angular_constant = (1-eps)*corrected_Xtail-(1/k+eps*J+delta_A)
            for error in (abs(physical_pressure_constant), abs(angular_constant), abs(measured_A-r)):
                max_error = max(max_error, error)
                if error > mp.mpf('1e-40')+tail_scale*pressure_remainder:
                    raise ArithmeticError('Independent physical moment fixture failed at Z='+str(Z))
            if Z == 0 and (delta_A <= 0 or abs(d1)+abs(d2) == 0):
                raise ArithmeticError('Positive Z0 heat deficit incorrectly erased')
            rows.append(dict(Z=str(Z), d1=str(d1), d2=str(d2),
                             physical_pressure_constant=str(physical_pressure_constant),
                             renormalized_angular_constant=str(angular_constant),
                             pressure_infinite_tail_error_bound=str(pressure_remainder)))
        return dict(rows=rows, max_absolute_error=str(max_error),
                    exact_Gamma_used=True, positive_Z0_heat_correction_retained=True,
                    true_physical_radii_used=True, actual_Md40_source_admission=False, passed=True)


def run():
    r = json.loads((HERE/NAME).read_bytes()); hashes = dict(r['input_hashes'])
    reference = r['paper_reference']
    if hashlib.sha256((HERE.parents[1]/reference['repository_path']).read_bytes()).hexdigest() != reference['sha256']:
        raise ValueError('Authoritative v2 paper extraction changed')
    for name, digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Closure dependency changed: '+name)
    proof = source_identities()
    if proof != r['proof']:
        raise ValueError('Recorded symbolic source proof changed')
    c = MPIntervalContext(); c.dps = 230
    read = lambda v: read_interval(c, v)
    for point in r['samples']:
        if not point['all_five_terminal_moment_identities_certified'] or not point['full_global_pressure_terminal_identity_verified']:
            raise ValueError('Admitted view lost absolute closure')
        if point['whole_outer_cone_certified'] or point['temporal_recursion'] or point['C4_outer_bounds_certified']:
            raise ValueError('Leading moment closure overclaimed downstream completion')
        if endpoints(read(point['P_over_Utheta_squared']['coefficients'][0]))[1] >= 0:
            raise ArithmeticError('Certified equivalent pressure tail must be negative')
        for key in ('Mp_over_Pstar_squared', 'P0_over_Pstar_squared'):
            if len(point[key]['coefficients']) != 2:
                raise ValueError('Original C1 history lost')
        for new, old in zip(point['P_over_Pstar_squared']['coefficients'], point['original_forward_pressure_enclosure']['coefficients']):
            nl, nh = endpoints(read(new)); ol, oh = endpoints(read(old))
            if max(nl, ol) > min(nh, oh):
                raise ArithmeticError('Equivalent pressure and original history are disjoint')
        if 'original_forward_angular_enclosure' in point:
            for new, old in zip(point['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'], point['original_forward_angular_enclosure']['coefficients']):
                nl, nh = endpoints(read(new)); ol, oh = endpoints(read(old))
                if max(nl, ol) > min(nh, oh):
                    raise ArithmeticError('Equivalent heat angular target and original history are disjoint')
    fixture = finite_parameter_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(source_identities=proof['identities'], independent_physical_integral_fixture=fixture,
                  source_sha256=r['source_sha256'], actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
                  original_histories_retained_and_equivalent_views_consistent=True,
                  all_five_terminal_moment_identities_certified=True, full_outer_five_moment_match=True,
                  C4_outer_bounds_certified=False, whole_outer_cone_certified=False,
                  global_admissible_stress_lift_constructed=False, temporal_recursion=False,
                  all_passed=True, input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('Absolute leading five moments:', len(proof['identities']), 'source identities and independent physical integral fixture PASS; higher regularity/stress/recursion pending', flush=True)
    return result


if __name__ == '__main__':
    run()
