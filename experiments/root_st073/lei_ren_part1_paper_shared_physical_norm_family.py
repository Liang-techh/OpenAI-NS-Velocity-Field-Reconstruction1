"""Uniform physical Section 9 norms for an increased-Cstar analytic family.

All mixed norms are in paper coordinates s=R/Ra or y=log(R/Ra), Z.
Positive majorants are represented by their logarithms; exp(logC) is never
formed. Existing finite/exit receipts are not relabelled as a new solution.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_shared_core_uniform_bounds import embedding
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent


class LogBounds:
    """Nonnegative scalar norm algebra. None is the exact zero bound.

sum uses max(log terms)+log(count), avoiding exponential subtraction of
enormous logarithms. This is a conservative upper bound, not log-sum-exp.
"""
    def __init__(self, c):
        self.c = c

    def hi(self, value):
        return self.c.mpf(endpoints(self.c.mpf(value))[1])

    def number(self, value):
        value = self.c.mpf(value)
        if endpoints(value)[0] < 0:
            raise ValueError('Positive norm bound required')
        return None if endpoints(value)[1] == 0 else self.hi(self.c.ln(value))

    def add(self, *logs):
        logs = [self.hi(x) for x in logs if x is not None]
        if not logs:
            return None
        top = max(endpoints(x)[1] for x in logs)
        return self.hi(self.c.mpf(top) + self.c.ln(len(logs)))

    def mul(self, *logs):
        if any(x is None for x in logs):
            return None
        return self.hi(sum(logs, self.c.mpf(0)))

    def scale(self, log, factor):
        return self.mul(log, self.number(factor))

    def reciprocal(self, norm_log, positive_floor, order):
        # At each base point, invert its Taylor jet by the finite Neumann
        # series. The nonconstant part is nilpotent through total order m.
        if endpoints(positive_floor)[0] <= 0:
            raise ValueError('Reciprocal requires a strictly positive value floor')
        floor_log = self.c.ln(positive_floor)
        if endpoints(norm_log-floor_log)[0] < 0:
            raise ValueError('Norm upper bound is below its asserted value floor')
        terms = [self.c.mpf(0)] + [self.hi(n*(norm_log-floor_log))
                                  for n in range(1, order+1)]
        return self.hi(-floor_log + self.add(*terms))


def run():
    names = ['lei_ren_part1_paper_shared_core_majorant.json',
             'lei_ren_part1_paper_shared_analytic_tube.json',
             'lei_ren_part1_paper_shared_linear_resolvent.json',
             'lei_ren_part1_paper_shared_core_uniform_bounds.json',
             'lei_ren_part1_paper_shared_frozen_angular_bounds.json',
             'lei_ren_part1_paper_shared_frozen_H_bounds.json']
    major, tube, linear, uniform, angular, frozen = [
        json.loads((HERE/n).read_bytes()) for n in names]
    hashes = {}
    for record in (major, tube, uniform, angular, frozen):
        if record['analytic_core_family_sha256'] != major['analytic_core_family_sha256']:
            raise ValueError('Analytic norm family mismatch')
        for name, digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Norm dependency changed: '+name)
            hashes[name] = digest
    if not (major['contraction_proved'] and uniform['all_supplied_core_exit_inputs_certified']
            and frozen['full_frozen_profile_test_9_14_certified']):
        raise ValueError('Analytic core/frozen input gates required')
    if (major['logC_definition'] != 'Lambda*Gupper + 2*logLambda + 1000'
            or len(major['terms']) != 20 or not major['all_paper_8_50_terms_included']):
        raise ValueError('Uniform amplitude majorant hypotheses changed')
    c = MPIntervalContext(); c.dps = 160
    with mp.workdps(200):
        b = LogBounds(c)
        get = lambda record, key: restore_value(c, record[key])
        lo = lambda x: c.mpf(endpoints(x)[0])
        logLam = get(major, 'logLambda'); lam = get(major, 'Lambda')
        logP = get(major, 'logPstar'); eps = get(major, 'epsilon')
        eta = get(tube, 'complex_tube_radius'); rho = lo(eta/2)
        h = get(tube, 'Xh_parameter'); G = b.hi(get(major, 'Gupper_in_logC_definition'))
        amp_log = b.hi(lam*G)  # |log(Cstar*F0)|<=Lambda*Gbar
        min_logC = b.hi(lam*G + 2*logLam + 1000)
        dt = b.hi(get(tube, 'delta')); j = get(major, 'required_j')
        ell = c.mpf('.01'); S = c.exp(ell); r = 4*eps
        if endpoints(4*S)[1] >= mp.mpf('4.1'):
            raise ValueError('Core continuation outside uniform radial domain')
        phi_floor = lo(get(uniform, 'actual_Phi_lower'))
        phi_norm = get(linear, 'Phi_model_Xh_norm_upper') + get(major, 'scaled_map_size_upper')
        psi_norm = get(major, 'fresh_Psi_model_Xh_norm_upper') + get(major, 'scaled_map_size_upper')

        # Weighted mixed norm sum_{i+k<=m} sup|d_s^i d_Z^k|/(i!k!).
        # It is submultiplicative. Its unweighted C3 sum is <=6 times it.
        ph = {}; vel = {}; exp_normalized = {}; N = {}; invN = {}
        for m in (3, 4):
            ph_terms = []; vel_terms = []
            for total in range(m+1):
                for i in range(total+1):
                    k = total-i
                    factor = 4**i*embedding(c, h, c.mpf('4.1'), i, k)
                    ph_terms.append(phi_norm*factor/(math.factorial(i)*math.factorial(k)))
                    base = (4+j if i == k == 0 else c.mpf(4) if i == 0 and k == 1 else c.mpf(0))
                    vel_terms.append((base+eps*psi_norm*factor)/(math.factorial(i)*math.factorial(k)))
            ph[m] = b.number(sum(ph_terms, c.mpf(0)))
            vel[m] = b.number(sum(vel_terms, c.mpf(0)))
            # Cauchy bounds for Lambda*G: factorial-weighted Z terms.
            gnorm = b.number(lam*G*sum((rho**(-k) for k in range(m+1)), c.mpf(0)))
            # Taylor jet of exp(g) is exp(g0)*sum_{n<=m} (g-g0)^n/n!.
            exp_normalized[m] = b.hi(amp_log + m*b.add(c.mpf(0), gnorm) + c.ln(m+1))
            N[m] = b.mul(exp_normalized[m], ph[m])
            invPhi = b.reciprocal(ph[m], phi_floor, m)
            invN[m] = b.mul(exp_normalized[m], invPhi)

        logphi_value = max(abs(x) for v in (phi_floor, b.hi(get(uniform, 'actual_Phi_upper')))
                           for x in endpoints(c.ln(v)))
        logphi_norm_log = b.add(b.number(logphi_value), *[
            b.hi(n*(ph[3]-c.ln(phi_floor))-c.ln(n)) for n in range(1, 4)])
        # This scalar norm is modest enough to materialize. A itself is
        # huge but has a compact mp exponent; exp(A) must stay implicit.
        logphi_norm = c.exp(logphi_norm_log)
        gnorm3 = lam*G*sum((rho**(-k) for k in range(4)), c.mpf(0))
        Vnorm3 = c.exp(vel[3])
        A = b.hi(10+6*(gnorm3+logphi_norm+Vnorm3))

        # Physical P0 C4 from the same complex pressure majorant; needed
        # because stress already contains one Z derivative.
        P = {m: b.mul(b.number(get(major, 'fresh_pressure_norm_upper')),
                       b.number(sum((rho**(-k) for k in range(m+1)), c.mpf(0))))
             for m in (3, 4)}

        # Core moment coefficient bounds: factor Cstar^{-1} from theta,
        # theta_z; Cstar^{-2} from pressure and the swirl energy term.
        # At Cstar>=1 dropping negative powers safely bounds all moments.
        core = {}
        for m in (3, 4):
            primitive = b.number(S+1); weight = primitive
            core[m] = {
                'theta': b.mul(b.number(2*r*r), primitive, weight, N[m]),
                'z': b.mul(b.number(r), primitive, vel[m]),
                'theta_z': b.mul(b.number(2*r*r), primitive, weight, N[m], vel[m]),
                'ztheta': b.add(b.mul(b.number(r), primitive, vel[m], vel[m]),
                               b.mul(b.number(r*r), primitive, weight, N[m], N[m])),
                'p': b.mul(b.number(r), primitive, N[m], N[m]),
            }

        # Frozen radius-function norms in y; constant subtraction is
        # bounded by addition. All five moments retain their entry values.
        def Rpower(power, m):
            value = (c.mpf(110)**power if power >= 0 else r**power)
            return b.number(value*sum((c.mpf(abs(power))**i/math.factorial(i)
                                       for i in range(m+1)), c.mpf(0)))

        moments = {}
        for m in (3, 4):
            Rminus = b.add(Rpower(1, m), b.number(r))
            R2minus = b.add(Rpower(2, m), b.number(r*r))
            moments[m] = {
                'theta': b.add(core[m]['theta'], b.mul(N[m], R2minus)),
                'z': b.add(core[m]['z'], b.mul(vel[m], Rminus)),
                'theta_z': b.add(core[m]['theta_z'], b.mul(N[m], vel[m], R2minus)),
                'ztheta': b.add(core[m]['ztheta'], b.mul(vel[m], vel[m], Rminus),
                               b.scale(b.mul(N[m], N[m], R2minus), c.mpf('.5'))),
                'p': b.add(core[m]['p'], b.mul(N[m], N[m], Rminus)),
            }
        zNorm = b.number(2); dNorm = b.number(4)
        deriv = lambda log: b.scale(log, 4)  # ||d_Z q||_3<=4||q||_4
        Lmin = lo(1-dt)
        Lnorm = b.number(1+3*dt)
        Linv = b.reciprocal(Lnorm, Lmin, 3)
        AopMz = b.add(b.mul(zNorm, moments[3]['z']),
                     b.mul(dNorm, deriv(moments[4]['z'])))
        Bnorm = b.add(Rpower(1, 3), AopMz)
        angular_numerator = b.add(moments[3]['theta'],
            b.scale(b.mul(zNorm, deriv(moments[4]['theta'])), c.mpf('.5')),
            b.mul(dNorm, deriv(moments[4]['theta_z'])),
            b.mul(zNorm, moments[3]['theta_z']))
        D = b.mul(Linv, b.add(Bnorm,
            b.scale(b.mul(Rpower(-1, 3), invN[3], angular_numerator), c.mpf('.5'))))
        pressure = {m: b.add(P[m], moments[m]['p']) for m in (3, 4)}
        Poperator = b.add(b.scale(b.mul(zNorm, pressure[3]), 2*(1+dt)),
                         b.mul(dNorm, deriv(pressure[4])))
        axial_numerator = b.add(b.mul(Bnorm, vel[3]),
            b.scale(moments[3]['z'], c.mpf('.5')),
            b.scale(b.mul(zNorm, deriv(moments[4]['z'])), c.mpf('.5')),
            b.scale(b.mul(zNorm, moments[3]['ztheta']), 2*dt),
            b.mul(dNorm, deriv(moments[4]['ztheta'])),
            b.mul(Rpower(1, 3), Poperator))
        # E has Cstar times the C^0 terms and Cstar^{-1} swirl terms.
        # Both are bounded by Cstar times the displayed coefficient norm.
        E = b.scale(b.mul(Linv, Rpower(c.mpf('-.5'), 3), invN[3], axial_numerator),
                    1/c.sqrt(2))
        Dfloor = lo(get(angular, 'frozen_D_whole_domain_lower'))
        invD = b.reciprocal(D, Dfloor, 3)

        K_terms = {
            'one_million': b.number(1000000), 'inverse_ell_c': b.number(1/ell),
            'inverse_Ra': b.number(1/r), 'Cstar_coefficient': c.mpf(0),
            'Pstar': b.hi(logP), 'A': b.number(A),
            'Fcore_coefficient': b.scale(N[3], 6),
            'inverse_Fcore_coefficient': b.scale(invN[3], 6),
            'P0_C3': b.scale(P[3], 6),
            'five_core_moments_coefficients': b.scale(b.add(*core[3].values()), 6),
            'Df_C3': b.scale(D, 6), 'Ef_C3_coefficient': b.scale(E, 6),
            'inverse_Df_C3': b.scale(invD, 6),
        }
        logKbar = b.add(*K_terms.values())
        log_one_plus_Kbar = b.add(c.mpf(0), logKbar)

        # Lemma10.3 (10.24) turns simultaneous9.17 into lower bounds
        # on Cstar. Use K(Cstar)<=Kbar*Cstar, with no circular choice.
        restrictions = {
            'complex_amplitude_guard': min_logC,
            'Cstar_at_least_exp_4A': b.hi(4*A),
            'long_reshape_before_reference': b.hi(40*A+1+c.ln(1+A)-logP),
            'axial_inherited_stress_radius': b.hi((8+2*log_one_plus_Kbar-12*logP)/8),
            'Rm_at_least_16': b.hi((c.ln(c.mpf(16)/110)+6)/10-logP),
        }
        selected_logC = b.hi(2*max(endpoints(v)[1] for v in restrictions.values()))
        logRref = b.hi(c.ln(110)+10*(selected_logC+logP))
        logKupper = b.hi(selected_logC+logKbar)
        margins = {
            'Cstar_minus_exp4A_log_margin': lo(selected_logC-4*A),
            'long_reshape_radius_log_margin': lo(10*(selected_logC+logP)-400*A-10-10*c.ln(1+A)),
            'axial_radius_log_margin': lo(8*selected_logC+12*logP-8-2*log_one_plus_Kbar),
            'Rm16_log_margin': lo(logRref-6-c.ln(16)),
        }
        if any(endpoints(v)[0] <= 0 for v in margins.values()):
            raise ArithmeticError('Selected family radius compatibility not proved')
        definition = dict(base_analytic_core_family_sha256=major['analytic_core_family_sha256'],
            pressure_source_sha256=major['implicit_source_sha256'],
            pressure_datum_sha256=major['datum_enclosure_sha256'],
            fixed_data='j,Lambda,delta,Pstar,P0,G and analytic tube held fixed',
            Cstar_family='logCstar>=Lambda*Gbar+2logLambda+1000',
            amplitude='F0=exp(-logCstar-Lambda*G)',
            selection='logCstar=2*max(explicit logarithmic lower bounds)',
            selected_logC=encode(selected_logC))
        family_sha = hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest()
        physical_units = dict(
            normalized_swirl='Fhat=Cstar*Fcore=exp(-Lambda*G)*Phi',
            inverse_normalized_swirl='1/Fhat=exp(Lambda*G)/Phi',
            angular_moments='Cstar*Mtheta and Cstar*Mtheta_z have the displayed coefficient bounds',
            pressure_moment='Cstar^2*Mp has the displayed coefficient bound',
            axial_moment='Mz has Cstar degree0',
            energy_moment='Mztheta has degrees0 and -2; bound the sum of both positive coefficients',
            frozen_D='angular Cstar^-1 factors cancel against f; remaining coefficient is degree0',
            frozen_E='E=Cstar*E0+Cstar^-1*E2 for each normalized profile; |E|/Cstar is bounded by both coefficients',
            variable_profiles='Phi and V can depend on Cstar through the nonlinear equations; their bounds are uniform, not identical point solutions',
            K='all displayed terms have degree at most1 after factoring physical amplitudes; Cstar>=1 gives K<=Kbar*Cstar')
        jet_proofs = dict(
            product='factorial-weighted Taylor-jet l1 norm is submultiplicative by coefficient convolution',
            exp='write g=g0+t in the degree-m truncated Taylor ring; t^(m+1)=0; exp(g)=exp(g0)*sum_{n=0}^m t^n/n!; bounded by exp(sup g)*(m+1)*(1+||g||_m)^m',
            log_Phi='write Phi=Phi0+t; bound constant log Phi0 separately and use log(1+t/Phi0)=sum_{n=1}^m (-1)^(n+1)(t/Phi0)^n/n; ||t||<=||Phi||, Phi0>=phi_floor',
            reciprocal='finite inverse jet (1/q0)*sum_{n=0}^m (-t/q0)^n; positive zeroth-order floors for Phi,L,Df come from the analytic positivity receipts',
            primitive='radial value plus shifted radial derivative coefficients are bounded by (S+1)*||integrand||_m; this deliberately overcounts positive terms',
            stress_derivative='weighted C3 of a Z derivative is <=4*weighted C4; all frozen stress numerators are charged four derivatives',
            physical_C3='for i+k<=3, i!k!<=6, so the unweighted derivative sum is <=6 times the weighted norm')
        frozen_extension = dict(
            angular='same uniform Phi floor/ceiling, Psi norms and correction norms hold for the full amplitude envelope; all forcing/barrier/Qentry/S0/moment-error bounds are unchanged as bounds',
            high_chi='chi and Lambda are fixed; polynomial slope majorant and Phi correction remain uniform; the signed-logF decomposition uses g_axis and normalized derivatives, not logCstar',
            low_chi='pressure atoms,Pstar,j and root localization stay fixed; velocity and cumulative-average bounds are uniform; F0<=exp(-2logLambda-1000) and its derivative envelope hold for all larger Cstar',
            H='low-chi N pressure floor survives the same physical F^2/F^2_Z upper bounds; f stays positive; the upper f bound gives the same lower |E|; D>0 permits H>=2|E|',
            Dfloor='the angular proof therefore keeps Df>=4*epsilon*Qglobal>0 for every family member',
            no_pointwise_monotonicity_claim='normalized Phi,V need not be monotone in Cstar; the unchanged rectangular analytic envelopes, not a monotone point profile, transfer the proof')
        result = dict(
            base_analytic_core_family_sha256=major['analytic_core_family_sha256'],
            uniform_Cstar_family_sha256=family_sha, definition=definition,
            implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            core_coordinates='s=R/Ra in [0,exp(.01)], Z in [-1,1]',
            frozen_coordinates='y=log(R/Ra) in [0,log(110/Ra)], Z in [-1,1]',
            weighted_norm='sum_{i+k<=m} sup|d_radial^i d_Z^k|/(i!k!); physical unweighted C3<=6*weighted',
            physical_amplitude_units=physical_units, finite_jet_proofs=jet_proofs,
            core_C4_used_for_frozen_stress_C3=True, A_upper=A,
            normalized_core_log_weighted_norms={str(m):dict(Phi=ph[m],V=vel[m],
                Cstar_Fcore=N[m],inverse_Cstar_Fcore=invN[m]) for m in (3,4)},
            core_moment_log_coefficient_norms={str(m):core[m] for m in (3,4)},
            frozen_moment_log_coefficient_norms={str(m):moments[m] for m in (3,4)},
            frozen_log_C3_weighted_bounds=dict(Df=D,Ef_divided_by_Cstar=E,inverse_Df=invD),
            K_term_log_coefficient_upper_bounds=K_terms, log_Kbar_upper=logKbar,
            K_growth_bound='K(Cstar)<=Kbar*Cstar for all Cstar>=complex minimum',
            Cstar_log_lower_bounds=restrictions,selected_logCstar=selected_logC,
            selected_logRref=logRref,selected_log_K_upper=logKupper,
            radius_log_margins=margins,
            uniform_analytic_fixed_point_admission_extended=True,
            extension_proof='the normalized fixed-point majorants depend on Cstar only through the F0^2 upper; the physical amplitude powers are separately recorded in physical_amplitude_units',
            frozen_branchwise_uniform_extension_proof=frozen_extension,
            frozen_input_proofs_extend_uniformly=True,
            full_physical_C3_K_norms_certified_for_uniform_analytic_family=True,
            radius_9_17_compatibility_certified_for_selected_analytic_family=True,
            actual_Cstar_F0_or_h_materialized=False, pressure_datum_held_fixed=True,
            finite_core_or_exit_receipts_relabelled=False,
            selected_family_finite_core_and_exit_bound=False,
            corrected_outer_profile_at_selected_Rref_built=False,
            K1_numeric_bound_certified=False,full_Section9_parameter_admission=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,
            input_hashes={**hashes,**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                for n in names+[Path(__file__).name]}})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result), indent=2)+'\n', encoding='utf-8')
        print('Uniform physical C3 K(Cstar)<=Kbar*Cstar bound generated; radius9.17 PASS', flush=True)
        print('log(A upper)=', mp.nstr(mp.log(endpoints(A)[1]),14),
              'log(logCstar)=',mp.nstr(mp.log(endpoints(selected_logC)[1]),14),flush=True)
        return result


if __name__ == '__main__':
    run()
