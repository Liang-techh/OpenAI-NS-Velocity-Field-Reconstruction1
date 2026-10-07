"""Whole-source modulation state errors and active-loop cone tolerance.

Reuse the accepted live all-N route. The bounds apply inside each continuous
cell, not just at its right endpoint. Original full inertial/pressure terms
are retained. The active-loop stability condition does not admit the q-flat
input cone, repaired band, global joins or a selected whole construction N.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls as allN
import lei_ren_part1_paper_compliant_current_native_C1_history_transfer as history_transfer
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as functions

HERE, PREFIX, sha = allN.HERE, allN.PREFIX, allN.sha
packets, repair, target = allN.packets, allN.repair, allN.target_module
LogUpper, Poly = repair.LogUpper, functions.FrequencyLogBound
ep, ZERO = allN.ep, allN.ZERO
NAME = PREFIX + 'current_native_generic_cone_state_errors.json'
RECEIPT = PREFIX + 'current_native_generic_cone_state_errors_check.json'
GATE = 'current_original_whole_modulation_C0_state_errors_and_active_loop_frequency_bound_certified'


def stress_increment(E, V, m, mZ, dE, dV, D, DZ, Z, delta):
    """Exact differences of the ORIGINAL full common-unit recovery operator.

    D/DZ are the five actual history differences and their ordinary Z rows.
    The original P0 cancels; original transport and all meridional cross
    terms remain. This callable accepts functions, symbolic or numeric data.
    """
    L, d = 1-delta*Z*Z, 1-Z*Z
    T = (1-delta)*Z*m+d*mZ
    DT = (1-delta)*Z*D['m']+d*DZ['m']
    return dict(
        transport=DT, pressure=D['p'], radial=(2*Z*dV-DT)/L,
        theta_linear=(-dE+(1-delta/2)*D['h']-(1-delta)*Z*DZ['h']/2)/L,
        theta_quadratic=((2*delta-1)*Z*D['k']-d*DZ['k']+E*DT+dE*T+dE*DT)/L,
        axial_linear=(-dV+(1-delta)*(D['m']-Z*DZ['m'])/2)/L,
        axial_quadratic=(V*DT+dV*T+dV*DT+2*delta*Z*D['e']-d*DZ['e']
                         +2*(1+delta)*Z*D['p']-d*DZ['p'])/L)


def inertial_error_envelopes(c, E, V, T, dE, dV, D, DZ):
    """Absolute C0 covers for full Z[-1,1], 0<=delta<=1/2.

    L^-1<=2; all ordinary Z history rows were integrated independently.
    Every polynomial has only negative N powers. No total-y jet is used.
    """
    C = lambda v: LogUpper.constant(c, v)
    constant = lambda v: Poly(c, {0: v})
    DT = D['m']+DZ['m']
    return dict(
        transport=DT, pressure=D['p'], radial=(dV.scale(C(2))+DT).scale(C(2)),
        theta_linear=(dE+D['h']+DZ['h'].scale(C('.5'))).scale(C(2)),
        theta_quadratic=(D['k']+DZ['k']+DT.scale(E)+dE.scale(T)+dE*DT).scale(C(2)),
        axial_linear=(dV+(D['m']+DZ['m']).scale(C('.5'))).scale(C(2)),
        axial_quadratic=(DT.scale(V)+dV.scale(T)+dV*DT+D['e']+DZ['e']
                         +D['p'].scale(C(3))+DZ['p']).scale(C(2)))


def prefix_history_polynomials(c, cell):
    """For 0<=s<=w: |D_p(s)|<=|D_p(in)|+density_cap*mass_rate(w).

    The accepted cell contribution is the same uniform density cover times
    its own full true-width mass. Its magnitude bounds every partial source
    integral. This use of absolute bounds does not assert signed cancellation.
    """
    def rows(incoming, contributions):
        return {key: Poly(c, {p: LogUpper.add(c, [target.magnitude(incoming[key][p]),
                            target.magnitude(contributions[key][p])]) for p in allN.ORDERS})
                for key in allN.RATES}
    return rows(cell['incoming'], cell['values']), rows(cell['incomingZ'], cell['Z_derivatives'])


def active_loop_margin(c, scales):
    """Explicit G1..G4 lower bounds on q!=0, from the original loop.

    |t|<=B*, aL>=2/(1+B*^2), v-2>=3eta/8,
    H-v>=m/2, 2(H-v)^2-(v-2)J^2>=m^2/4.
    These last two statements are ACTIVE-branch estimates only.
    """
    read = lambda row: packets.interval(c, row)
    upper = scales['logarithmic_conservative_upper_constants']
    lower = scales['logarithmic_selected_positive_lower_constants']
    add = lambda *logs: LogUpper.add(c, [LogUpper(c, v) for v in logs]).log
    logB = read(upper['B_star'])
    loga = c.ln(2)-add(c.mpf(0), 2*logB)
    logm, logeta = read(lower['margin_min']), read(scales['selected_positive_eta_log'])
    margins = dict(G1=loga, G2=loga+logeta+c.ln(c.mpf(3)/8),
                   G3=loga+logm-c.ln(2), G4=3*loga+2*logm-c.ln(4))
    logg = c.mpf(min(ep(v)[0] for v in margins.values()))
    state_logs = [c.ln(3)] + [read(upper[key]) for key in ('p1_abs_max', 'p2_abs_max') if upper[key] is not None]
    max_state = c.mpf(max(ep(v)[1] for v in state_logs))
    logH = add(c.ln(2), max_state)
    logrho = c.mpf(min(0, ep(logg-c.ln(512)-5*logH)[0]))
    return dict(active_q_nonzero_only=True, log_a_L_positive_lower=loga,
                G_log_positive_lowers=margins, log_g_positive_lower=logg,
                log_H_state_upper=logH, log_rho_stability_positive=logrho,
                stability='G_j(X+error)>=g/2 if max|error|<=rho=min(1,g/(512*H^5))',
                inactive_q_zero_and_post_modulation_input_margins_not_supplied=True)


def negative_power_threshold(c, rows, logrho):
    """For N>=1, sum C_p*N^p <= (sum C_p)/N, p<=-1."""
    result = {}
    for key, poly in rows.items():
        if any(p>=0 for p in poly.terms):
            raise ValueError('Cone C0 state error must have only negative N powers')
        total = LogUpper.add(c, list(poly.terms.values()))
        threshold = c.mpf(0) if total.log is None else c.mpf(max(0, ep(total.log-logrho)[1]))
        result[key] = dict(coefficient_sum=total.record(), sufficient_log_N_lower=threshold,
                           higher_negative_powers_bounded_by_N_minus1_only_for_N_ge1=True)
    return result


class NativeGenericConeStateErrors:
    def __init__(self, owner):
        if type(owner) is not allN.NativeRcAllNFunctionControls:
            raise TypeError('Existing accepted original all-N function owner required')
        checked = json.loads((HERE/allN.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[allN.GATE] or checked['source_family']!=owner.family:
            raise ValueError('Checked same-original all-N function/range family required')
        self.owner, self.ctx, self.family = owner, owner.ctx, owner.family
        self.hashes = {**owner.hashes, **checked['input_hashes'], allN.NAME: sha(allN.NAME),
                       allN.RECEIPT: sha(allN.RECEIPT), Path(__file__).name: sha(Path(__file__).name)}
        owner.service.bind_hashes(self.hashes)

    @allN.paired.native.inlet.source_precision
    def compute(self, live):
        c = self.ctx
        if len(live['cells'])!=24 or live['Z_box']._mpi_!=c.mpf((-1,1))._mpi_:
            raise ValueError('Accepted whole original24-cell/full-Z all-N route required')
        delta = c.mpf(self.owner.service.data['delta'])
        if ep(delta)[0]<0 or ep(delta)[1]>ep(c.mpf('.5'))[0]:
            raise ValueError('Original 0<=delta<=1/2 required')
        root = self.owner.target.q_owner.owner.owner
        signed = self.owner.target.transfer.owner.signed_owner
        coords = self.owner.target.coordinates
        C = lambda v: LogUpper.constant(c, v)
        magnitude = target.magnitude
        decode = lambda row: LogUpper(c, None if row['exact_zero'] else packets.interval(c, row['log_absolute_upper']))
        poly = lambda cap, p=-1: Poly(c, {p: cap})
        scale_exp = LogUpper(c, c.mpf(1)/128)
        relative_delta = poly(scale_exp*C(c.mpf(5)/4))
        margin = active_loop_margin(c, root.scales)
        cells = []; branch_count = 0; max_logN = c.mpf(0)
        for cell in live['cells']:
            label, chart = cell['record']['label'], cell['record']['chart']
            D, DZ = prefix_history_polynomials(c, cell)
            branches = []
            for item in cell['branches']:
                conditional = item['conditional']; roots = conditional['source']['roots']
                packet = conditional['source']['packet']; coords.require_family(packet.source_family)
                original = history_transfer.packet_history_functions(packet, coords, signed)
                E = magnitude(roots['E'][ZERO])
                ordinaryV = target.prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0], 0)
                Vleaf = target.prior.signed.expressions.RadiusPolynomial(packet.algebra, {0: ordinaryV})
                V = magnitude(signed.leaf(Vleaf, roots['E'][ZERO].scale.bases, roots['E'][ZERO].ledger))
                T = LogUpper.add(c, [magnitude(original['originals']['m']), magnitude(original['Z_derivatives']['m'])])
                dE = poly(magnitude(item['density']['values']['h'][-1]))
                dV = poly(magnitude(item['density']['values']['m'][-1]))
                errors = inertial_error_envelopes(c, E, V, T, dE, dV, D, DZ)
                positive = root.decode(root.inventory[chart]['actual_positive_denominator_theorem'])
                logE = positive['log_E_positive_lower']
                logR = packet.provenance['logR_cover']
                R_over_E = LogUpper(c, logR-logE)
                S = LogUpper(c, c.mpf(self.owner.service.data['logP']))
                pstate = {}
                for name, axis in (('p1', 'theta'), ('p2', 'axial')):
                    numerator = errors[axis+'_linear']+errors[axis+'_quadratic'].scale(S)
                    pstate[name] = (numerator.scale(R_over_E)+relative_delta.scale(magnitude(roots[name][ZERO]))).scale(scale_exp)
                loops = root.inventory[chart]['slow_phase_held_log_bounds']
                Ay = decode(loops['A']['y1_Z0']); By = decode(loops['B']['y1_Z0'])
                bL = LogUpper.add(c, [C(3), magnitude(roots['b'][ZERO])])
                pstate['a'] = poly(Ay*C(2))
                pstate['b'] = relative_delta.scale(bL)+poly(By*C(2)*LogUpper(c, -logE)*scale_exp)
                requirements = negative_power_threshold(c, pstate, margin['log_rho_stability_positive'])
                branch_logN = c.mpf(max(ep(v['sufficient_log_N_lower'])[1] for v in requirements.values()))
                max_logN = c.mpf(max(ep(max_logN)[1], ep(branch_logN)[1]))
                record = dict(original_source_provenance=packet.provenance,
                    original_cutoff_branch=item['record']['original_cutoff_branch'],
                    original_log_E_positive_lower=logE, original_log_R_cover=logR,
                    original_Pstar_log= S.log,
                    original_transport_C0_absolute_cap=T.record(),
                    full_inertial_and_radial_error_log_polynomials={k:v.record() for k,v in errors.items()},
                    normalized_state_error_log_polynomials={k:v.record() for k,v in pstate.items()},
                    active_loop_state_frequency_requirements=requirements,
                    active_loop_sufficient_log_N_lower=branch_logN,
                    original_P0_cancels_in_pressure_difference=True,
                    nonzero_original_V_transport_and_meridional_cross_terms_retained=True,
                    source_caps_are_bounds_only_not_defining_function_values=True)
                branches.append(dict(record=record, original=original, errors=errors, state=pstate, requirements=requirements))
            branch_count += len(branches)
            cells.append(dict(record=dict(label=label, chart=chart,
                original_geometry=cell['geometry']['record'],
                inside_continuous_cell_history_C0_coefficient_caps={k:v.record() for k,v in D.items()},
                inside_continuous_cell_history_Z_coefficient_caps={k:v.record() for k,v in DZ.items()},
                partial_history_recipe='exp(-rate*s)*incoming+integral_0^s exp(-rate*(s-t))*density(t)dt, 0<=s<=true_width',
                own_rate_mass_applied_once_and_incoming_memory_retained=True,
                original_branch_state_error_records=[v['record'] for v in branches],
                exact_original_zero_initial_collar=label=='initial_flat_collar'), D=D, DZ=DZ, branches=branches))
        if branch_count!=live['source_branch_count']:
            raise ValueError('All original all-N conditional source branches required')
        source_repair = allN.combined_log_conditions(self.owner, live)
        combined = c.mpf(max(ep(max_logN)[1], ep(source_repair['source_and_repair_sufficient_log_N_lower'])[1]))
        return dict(cells=cells, margin=margin, active_state_logN=max_logN,
                    source_repair_conditions=source_repair, combined_logN=combined, branch_count=branch_count)


@allN.paired.native.inlet.source_precision
def run(owner, live, *, return_live=False):
    began = time.monotonic(); adapter = NativeGenericConeStateErrors(owner); got = adapter.compute(live)
    result = dict(source_family=adapter.family, **{GATE: True},
        active_frozen_loop_polynomial_margin_and_stability=got['margin'],
        original_continuous_cell_state_error_records=[row['record'] for row in got['cells']],
        actual_source_branches=got['branch_count'], original_continuous_cells=24, full_Z_domain=[-1,1],
        active_loop_sufficient_log_N_lower=got['active_state_logN'],
        original_source_repair_and_active_loop_log_N_lower=got['combined_logN'],
        existing_source_repair_conditions=got['source_repair_conditions'],
        local_domain_integer_N_lower=160,
        exact_recovery_operator='stress_increment with actual full C0/Z original and defect functions',
        independent_of_total_y_higher_jet_smallness=True,
        existing_live_all_N_route_reused=True, original_source_ancestor_constructors_called=False,
        active_branch_cone_persistence_conditional_on_same_functions_and_recorded_N=True,
        inactive_input_and_post_modulation_cone_margins_certified=False,
        correction_band_partial_moments_and_cone_errors_certified=False,
        actual_five_controls_installed=False, certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False, current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN, False), input_hashes=adapter.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Whole original24-cell/full-Z C0 modulation state error envelopes, including within-cell histories, full inertial/pressure/radial cross terms and normalized shears. Quantitative frozen q-active polynomial cone tolerance and source+repair+active-state lower threshold. No q-flat/quiet/corrected-band/global cone or actual solved field admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result), indent=2)+'\n', encoding='utf8')
    print('Original continuous modulation state errors and active-loop log-N tolerance produced; global cone remains open', flush=True)
    return (result, adapter, got) if return_live else result
