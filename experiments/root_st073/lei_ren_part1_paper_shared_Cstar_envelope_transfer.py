"""Bind completed interval envelopes to the enlarged-Cstar analytic family.

This is an inclusion proof. It does not rename historical solutions or
pretend a fresh point-amplitude run took place. The existing S jets and h
boxes enclose every selected family member; all other scaled inputs stay
fixed. A corrected outer construction at the new Rref is still required.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_gauge_core import _unpack_vector
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def run():
    stems = ['shared_physical_norm_family', 'shared_core_majorant',
             'shared_analytic_tube', 'shared_interval_core_Z049_Z051',
             'shared_core_seed_check', 'shared_core_tail_admission',
             'shared_tolerance_core_exit', 'shared_tolerance_comparison',
             'shared_tolerance_exit_bridge', 'shared_tolerance_exit_continuation',
             'shared_tolerance_exit_switch', 'shared_tolerance_R110_cone',
             'shared_tolerance_phase_check']
    records = {n: json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in stems}
    hashes = {}
    for record in records.values():
        for name, digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Transfer dependency changed: '+name)
            hashes[name] = digest
    family = records['shared_physical_norm_family']
    major = records['shared_core_majorant']; seed = records['shared_core_seed_check']
    tube = records['shared_analytic_tube']; core = records['shared_interval_core_Z049_Z051']
    state_path = HERE/core['state_file']
    if hashlib.sha256(state_path.read_bytes()).hexdigest() != core['state_sha256']:
        raise ValueError('Terminal core state hash changed')
    state = json.loads(state_path.read_bytes())
    if not (core['completed_target'] and seed['finite_core_completed']
            and seed['seed_acceptance_passed'] and records['shared_core_tail_admission']['target_met']):
        raise ValueError('Completed core/seed/tail gates required')
    if hashlib.sha256(json.dumps(state['fixed'], sort_keys=True).encode()).hexdigest() != seed['fixed_seed_sha256']:
        raise ValueError('Accepted scaled input jets changed')
    base_sha = family['base_analytic_core_family_sha256']
    if any(r['analytic_core_family_sha256'] != base_sha for r in (major, seed, tube)):
        raise ValueError('Base analytical identity mismatch')
    if state['identity']['analytic_core_family_sha256'] != base_sha:
        raise ValueError('Finite envelope base identity mismatch')
    phase = records['shared_tolerance_phase_check']
    parameter_sha = phase['production_family_sha256']
    downstream = stems[7:12]
    if not phase['all_five_receipts_share_same_h_parameter']:
        raise ValueError('Shared h identity missing')
    if any(records[n]['parameter_family_sha256'] != parameter_sha for n in downstream):
        raise ValueError('Exit envelopes use different h families')
    if not records['shared_tolerance_R110_cone']['relaxed_3_23_cone_certified_at_R110_local_axis']:
        raise ValueError('Local endpoint direction gate required')

    c = MPIntervalContext(); c.dps = state['identity']['precision']
    with mp.workdps(c.dps+40):
        get = lambda r, n: read_interval(c, r[n])
        logC = get(family, 'selected_logCstar')
        logLambda = get(major, 'logLambda'); lam = get(major, 'Lambda')
        Gbar = get(major, 'Gupper_in_logC_definition')
        minLogC = lam*Gbar+2*logLambda+1000
        if endpoints(logC-minLogC)[0] <= 0:
            raise ValueError('Selected Cstar outside amplitude envelope')
        # epsilon^2 F0^2: two epsilon powers and two physical amplitudes.
        new_logS_upper = -2*logLambda-2*logC+2*lam*Gbar
        old_logS_upper = -6*logLambda-2000
        if endpoints(new_logS_upper-old_logS_upper)[1] >= 0:
            raise ValueError('New amplitude jets not included')
        rho = get(tube, 'complex_tube_radius')/2
        S = _unpack_vector(c, state['fixed']['S_Z_taylor'])
        old_bound = c.exp(old_logS_upper)
        for k, stored in enumerate(S):
            upper = endpoints(old_bound/rho**k)[1]
            expected = (mp.mpf(0), upper) if k == 0 else (-upper, upper)
            if endpoints(stored) != expected:
                raise ValueError('Stored Cauchy envelope mismatch at order '+str(k))
        hBox = get(records['shared_tolerance_R110_cone']['parameters'], 'shared_h_epsilon_box')
        log_h_upper = get(records['shared_tolerance_R110_cone']['parameters'], 'log_h_upper')
        new_log_h_upper = -100*logC-c.ln(2)
        if endpoints(new_log_h_upper-log_h_upper)[1] >= 0 or endpoints(hBox)[0] != 0:
            raise ValueError('Selected h outside completed local exit envelope')
        proof = dict(
            axis='ell_scaled=-g_axis, U0=4Z+j, epsilon*P0 unchanged; Cstar enters scaled recurrence only via S=epsilon^2 F0^2',
            S='selected logC increases: log|S|<=-2logLambda-2logC+2Lambda Gbar<=old complex supremum',
            derivatives='Cauchy on the same eta/2 disks includes all 148 stored S Taylor jets',
            finite_rows='directed polynomial/rational recurrence is inclusion-preserving on stored denominator guards; no point amplitude was used',
            tail='same uniform Xh fixed-point majorant and degree144 tail hold for every smaller amplitude',
            width='K>=Cstar and cstar<=1/2 imply log h<=-100logC-log2, inside the original shared positive-width box',
            exit='scaled transfer, cutoff phase, comparison and actual exit use only included S/h boxes and unchanged fixed inputs',
            limitation='local exit C1 family is [0.49,0.51]; no global actual cone follows from this inclusion')
        result = dict(
            base_analytic_core_family_sha256=base_sha,
            selected_uniform_Cstar_family_sha256=family['uniform_Cstar_family_sha256'],
            included_exit_envelope_parameter_family_sha256=parameter_sha,
            implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            terminal_state_sha256=core['state_sha256'], finite_radial_degree=144,
            scaled_S_Cauchy_jet_inclusion_checks=len(S),
            selected_logCstar=logC,selected_logS_upper=new_logS_upper,
            original_logS_upper=old_logS_upper,selected_log_h_upper=new_log_h_upper,
            original_log_h_upper=log_h_upper,proof=proof,
            larger_Cstar_core_and_local_exit_envelopes_bound=True,
            existing_interval_rows_reused_by_explicit_inclusion=True,
            stored_receipts_or_point_solutions_relabelled=False,
            fresh_point_amplitude_or_finite_run_claimed=False,
            local_R110_relaxed_direction_bound_transferred=True,
            whole_axis_actual_exit_cone_certified=False,
            corrected_outer_at_selected_radius_built=False,
            five_terminal_moment_identities_repaired=False,
            K1_numeric_bound_certified=False,full_Section9_parameter_admission=False,
            temporal_recursion=False,
            input_hashes={**hashes,core['state_file']:core['state_sha256'],
                **{PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in stems},
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result), indent=2)+'\n',encoding='utf-8')
        print('Enlarged Cstar inclusion: 148 S jets and completed local R110 envelopes bound; no recurrence rerun',flush=True)
        return result


if __name__ == '__main__':
    run()
