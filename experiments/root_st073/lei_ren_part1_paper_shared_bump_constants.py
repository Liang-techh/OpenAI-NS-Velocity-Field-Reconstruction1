"""Directed fixed-bump constants and Section 9/10 inlet tolerance.

Uses fixed integral enclosures, not any older candidate's defect or inverse
solution. Integer majorants are independent of physical input parameters.
The C_S proof and its inherited input hypotheses are recorded separately
from their admission by the current field.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_five_bump_inverse import weights, matmul
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent


def upper(c, value):
    lo, hi = endpoints(value)
    return c.mpf(max(abs(lo), abs(hi)))


def maximum(values):
    return max(values, key=lambda x: endpoints(x)[1])


def norm1(c, matrix):
    return maximum([sum((upper(c, row[j]) for row in matrix), c.mpf(0))
                    for j in range(len(matrix[0]))])


def minimum(c, a, b):
    al, ah = endpoints(a)
    bl, bh = endpoints(b)
    return c.mpf([min(al, bl), min(ah, bh)])


def ceiling(value):
    return max(1, int(mp.ceil(endpoints(value)[1])))


def run():
    c = MPIntervalContext()
    c.dps = 160
    names = ['lei_ren_part1_paper_bump_integral_enclosures_check.json',
             'lei_ren_part1_paper_shared_pressure_Kp.json',
             'lei_ren_part1_paper_analytic_radial_tail.json']
    raw, pressure, tube = [json.loads((HERE / n).read_bytes()) for n in names]
    if not raw['directed_integrals_certified']:
        raise ValueError('Directed fixed bump weights required')
    for name, key in [('lei_ren_part1_paper_bump_integral_enclosures.py', 'module_sha256'),
                      ('lei_ren_part1_paper_bump_integral_enclosures_check.py', 'check_source_sha256')]:
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != raw[key]:
            raise ValueError('Fixed bump dependency changed: ' + name)
    for name, digest in pressure['input_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pressure dependency changed: ' + name)
    with mp.workdps(c.dps + 40):
        W = weights(c, raw)
        L = W['L']
        midpoint = mp.matrix([[(endpoints(v)[0] + endpoints(v)[1]) / 2
                               for v in row] for row in L])
        inverse = mp.inverse(midpoint)
        R = [[c.mpf(inverse[i, j]) for j in range(5)] for i in range(5)]
        RL = matmul(R, L)
        error = [[c.mpf(int(i == j)) - RL[i][j] for j in range(5)] for i in range(5)]
        q = norm1(c, error)
        if endpoints(q)[1] >= 1:
            raise ArithmeticError('Neumann inverse gate failed')
        inverse_bound = norm1(c, R) / (1 - q)
        CA = ceiling(inverse_bound)

        # The norm is sum_i(sup|h_i|+sup|h_i_Z|), a Banach algebra.
        # Off-support products vanish identically. Each ordered pair sums
        # its absolute output coefficients in the symmetric polarization.
        pair = [[c.mpf(0) for _ in range(5)] for _ in range(5)]
        for i, angular in enumerate((2, 4)):
            pair[i][angular] = pair[angular][i] = W['fg'][i] / 2
            pair[i][i] = 40 * W['gg'][i]
        for i in range(3):
            pair[i+2][i+2] = (W['ff'][i] + W['ff_over_x'][i]) / 2
        quadratic_bound = maximum([upper(c, v) for row in pair for v in row])
        CQ = ceiling(quadratic_bound)

        # beta(x)=exp(-1/(1-(x/r)^2))/(r*N), r=1/40.
        # b<=1; |b'|<=2 q^2 exp(-q)<=8 (q>=1), including endpoint limits.
        # Conservative simple analytic bounds avoid a derivative sample grid.
        N = restore_value(c, raw['normalization'])
        radius = c.mpf(1) / 40
        B = 1 / (radius * N)
        D = 8 / (radius * radius * N)
        da = 4 * (2 * D + B / 10)
        dzeta = 2 * B
        velocity_norm = c.mpf('2.1') * B + 4 * D
        sq_error = 1 + 2 * da + 12 * B
        CS = ceiling(1000 * (1 + B + 2 * D))
        if any(endpoints(v)[1] > CS for v in (da, dzeta, B+2*D, velocity_norm, sq_error, c.mpf(3))):
            raise ArithmeticError('Selected C_S misses a displayed coefficient')
        t_cap = c.mpf(1) / (100 * CS)
        denominator = 1 - B * t_cap
        if endpoints(denominator)[0] <= c.mpf('.99'):
            raise ArithmeticError('Corrected velocity denominator not protected')
        # Extremal consequences at C_S*t<=.01, delta<=.001, epsilon0<=1e-6.
        amplitude_lower = c.exp(c.mpf('-.6')) * denominator / 2
        amplitude_upper = c.exp(c.mpf('-.6')) * (c.exp(c.mpf('.1')) + B*t_cap)
        a_error = da * t_cap
        sq_lower = c.mpf('1.8') - c.mpf('.001')/2 - 10*c.mpf('1e-6') - CS*t_cap
        if not (endpoints(amplitude_lower)[0] > mp.mpf(1)/8 and
                endpoints(amplitude_upper)[1] < 1 and
                endpoints(a_error)[1] < mp.mpf('.1') and
                endpoints(sq_lower)[0] > mp.mpf('1.4')):
            raise ArithmeticError('Elementary velocity/source consequences failed')

        KN = pressure['KN']
        t_star = c.mpf(1) / (10000 * KN * CS)
        e_contraction = c.mpf(1) / (8 * CA**2 * CQ)
        e_source = t_star / (2 * CA)
        e_star = minimum(c, e_contraction, e_source)
        epsilon0 = restore_value(c, pressure['epsilon0'])
        eta_tol = minimum(c, epsilon0/4, e_star/100)
        required_j = eta_tol/8
        current_j = restore_value(c, tube['j'])
        disjoint = endpoints(current_j)[0] > endpoints(required_j)[1]
        current_delta = restore_value(c, tube['delta'])
        elementary_inputs = (endpoints(epsilon0)[1] <= mp.mpf('1e-6') and
                             endpoints(current_delta)[1] <= mp.mpf('.001'))
        result = dict(
            paper_equations=['10.1','10.3','10.6','10.8','10.9','10.12','9.4'],
            norm='sum_i(sup_Z |h_i|+sup_Z |h_i_Z|), Z in [-1,1]',
            fixed_bump_radius='1/40', centers=['5/4','3/2','7/4'],
            fixed_matrix=L, fixed_preconditioner=R, neumann_l1_error=q,
            inverse_l1_bound=inverse_bound, CA=CA,
            quadratic_ordered_pair_bounds=pair, quadratic_C1_bound=quadratic_bound, CQ=CQ,
            amplitude_inverse_square_C1_bound=40,
            bump_sup_bound=B, bump_first_derivative_sup_bound=D, CS=CS,
            elementary_coefficients=dict(a_error=da,zeta_error=dzeta,
                velocity_C1_radial_norm=velocity_norm,SQ_source_error=sq_error),
            consequences_at_smallness_cap=dict(u_over_Pstar_lower=amplitude_lower,
                u_over_Pstar_upper=amplitude_upper,a_deviation=a_error,SQ_lower=sq_lower),
            CS_proof_hypotheses=dict(delta_upper='.001',epsilon0_upper='1e-6',
                Pstar_lower=1,coefficient_smallness='CS*||h||_1<=.01',
                entry_average='||Mz/R-4Z||_1<=2epsilon0',
                entry_pressure='||P||_1<=(Kp+100)Pstar^2',
                reference_velocity='u=Am*x^.1,V=4Z before correction'),
            CS_proof=dict(
                disjoint_unit_bumps='||f||_1,||g||_1<=B*t; ||x*f_x||_1,||x*g_x||_1<=2D*t',
                a='|a-.8|<=4(2D+.1B)t, since x^.1+f>=.99',
                zeta='|zeta+2Z/(1+Z^2)|<=2B*t',
                average='partial unit bump mass<=1, so |W-W0|<=2epsilon0+t',
                pressure='||Am^2||_1<=2Pstar^2; integral beta^2<=B; change<=3t*Pstar^2',
                source='|SQ-SQ0|<=10epsilon0+(1+2*a_coefficient+12B)t',
                relative_velocity='norm in10.6 <=(2.1B+4D)t, then t<=2CA*e'),
            KN=KN,epsilon0=epsilon0,t_star=t_star,
            e_star_contraction_branch=e_contraction,e_star_source_branch=e_source,
            e_star=e_star,eta_tol=eta_tol,required_j=required_j,current_source_j=current_j,
            current_j_over_required_j=current_j/required_j,
            current_source_j_exceeds_selected_tolerance=disjoint,
            fixed_CA_CQ_directed_bounds_certified=True,
            fixed_CS_analytic_bound_certified_under_stated_inputs=True,
            current_delta_epsilon0_satisfy_elementary_limits=elementary_inputs,
            current_source_j_eta_tol_relation_verified=False,
            tighter_constants_may_change_selected_tolerance=True,
            current_source_rejected_for_this_selected_hierarchy=disjoint,
            full_Section9_parameter_admission=False,
            full_input_average_pressure_or_defect_tests_certified=False,
            old_source_inverse_or_controls_used=False,temporal_recursion=False,
            input_hashes={**pressure['input_hashes'],**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                for n in names+[Path(__file__).name,'lei_ren_part1_paper_interval_five_bump_inverse.py',
                'lei_ren_part1_paper_interval_exit_continuation_enclosure.py',
                'lei_ren_part1_paper_bump_integral_enclosures.py',
                'lei_ren_part1_paper_bump_integral_enclosures_check.py']}}
        )
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Fixed bump constants: CA=%d CQ=%d CS=%d' % (CA,CQ,CS))
        print('Required j:',mp.nstr(endpoints(required_j)[1],14),
              'current source exceeds this tolerance:',disjoint)
        return result


if __name__ == '__main__':
    run()
