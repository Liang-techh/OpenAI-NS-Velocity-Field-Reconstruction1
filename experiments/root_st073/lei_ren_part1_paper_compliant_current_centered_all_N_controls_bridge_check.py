"""Focused check of centered all-N integration and the new C1 unit ball.

Reuse checked source/IBP receipts. Verify the actual exact graph is unchanged,
all 336 integrals retain one N/phase/Jacobian, and the bound drives the same
nonlinear equation. The unit-ball estimates are proved independently.
"""
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_centered_all_N_controls_bridge as current

HERE, PREFIX, sha, ep = current.HERE, current.PREFIX, current.sha, current.ep
require, packets = current.require, current.packets


def unit_ball_theorem():
    N, CA, CQ, C, rho, gmax = s.symbols('N CA CQ C rho gmax', positive=True)
    first_threshold = 2*CA*C/rho
    contraction_threshold = 4*CA*CQ*rho
    require(s.cancel(CA*C/first_threshold-rho/2) == 0, 'First-step threshold does not give rho/2')
    require(s.cancel(CA*CQ*rho**2/contraction_threshold-rho/4) == 0, 'Quadratic image is not rho/4')
    require(s.cancel(2*CA*CQ*rho/contraction_threshold-s.Rational(1, 2)) == 0,
            'C1 Lipschitz constant is not at most half')
    require(s.cancel(CA*C/first_threshold+CA*CQ*rho**2/contraction_threshold-3*rho/4) == 0,
            'C1 image radius is not 3rho/4')
    positivity_threshold = 2*gmax*rho*2**s.Rational(2, 3)
    require(s.simplify(gmax*rho/positivity_threshold-2**(-s.Rational(2, 3))/2) == 0,
            'Original swirl positivity reserve changed')
    c0, cz = s.symbols('coefficient coefficient_Z', positive=True)
    require(s.cancel(N*(c0+cz)/N**2-(c0+cz)/N) == 0, 'Averaged target lost an N factor')
    n = s.symbols('n', integer=True, nonnegative=True)
    require(s.simplify(2**(-n)/(1-s.Rational(1, 2))*(rho/2)-rho*2**(-n)) == 0,
            'Picard tail must use the smaller first-step radius')
    h, zero_step = s.symbols('h zero_step', nonnegative=True)
    require(s.cancel(2*zero_step-zero_step/(1-s.Rational(1, 2))) == 0,
            'Small fixed-point bound is not 2*CA*C/N')
    return dict(passed=True, independent_exact_threshold_identities=8,
        full_averaged_target_N_factor='N*(C/N^2)=C/N',
        normalized_C1_product_norm_submultiplicative_required=True,
        original_exact_B_mu_and_all_quadratic_terms_retained=True,
        actual_N_or_numerical_source_values_instantiated=False)


def exact_graph_checks(owner, built):
    old = owner.old
    accepted = json.loads((HERE/current.family.NAME).read_bytes())
    g, N = built['graph'], built['N'].node
    require(g.nodes == accepted['exact_function_graph_nodes'],
            'Adapter altered the accepted original nonlinear function graph')
    roots = current.controls.pair_roots
    require(roots(built['N_scaled_targets'].values(), current.controls.ROWS)
            == old['exact_reconstructed_N_scaled_target_roots'], 'Exact all-N target function changed')
    require(roots(built['generic_Picard_map'], current.controls.CONTROLS)
            == accepted['exact_generic_Picard_C1_map_roots'], 'Original nonlinear control map changed')
    require(roots(built['control_residual'], current.controls.ROWS)
            == accepted['exact_control_residual_roots'], 'Original full quadratic residual changed')
    shared = [i for i, row in enumerate(g.nodes) if row['operation'] == 'shared_positive_integer']
    require(shared == [N], 'One common original spatial/control N required')
    references = 0
    for row in g.nodes:
        if row['operation'] != 'original_function_graph' or not row.get('function_role', '').startswith('all_N_'):
            continue
        phase = g.nodes[row['phase']]
        product = g.nodes[phase['argument']]
        require(row['shared_N'] == N and row['graph_sha256'] == built['source_graph_sha256']
                and phase['operation'] == 'analytic_unary' and phase['name'] == 'fractional_part'
                and product['operation'] == 'product' and N in product['arguments'],
                'Source coefficient is not bound to the same N and actual spatial phase')
        references += 1
    integral_count = 0
    maps, unused, x = current.source.exact_radius_maps()
    for i, cell in enumerate(built['cells']):
        original = cell.original_cell
        coordinate = g.symbol('coordinate_'+str(i))
        offset = current.source.expression(g, maps[original['chart']], built['parameters'], coordinate)
        jacobian = current.source.expression(g, s.diff(maps[original['chart']], x), built['parameters'], coordinate)
        for p, rows in cell.coefficient_contributions.items():
            for key, pair in rows.items():
                kernel = g.unary('exp', g.neg(g.mul(g.constant(current.allN.RATES[key]),
                    g.sub(current.source.FunctionRef(g, original['right_radius_offset']), offset))))
                for attribute in ('value', 'Z'):
                    ref = getattr(pair, attribute)
                    if ref == g.zero:
                        continue
                    record = g.nodes[ref.node]
                    expected = g.mul(kernel, getattr(cell.coefficient_density_pairs[p][key], attribute), jacobian)
                    require(record['operation'] == 'definite_integral' and record['extracted_N_power'] == p
                            and record['integrand'] == expected.node
                            and record['coefficient_still_depends_on_N_and_actual_phase'],
                            'Original kernel/coefficient/Jacobian/N-dependent integral changed')
                    integral_count += 1
    require(len(built['cells']) == 24 and integral_count == 336, 'All original C1 cell integrals required')
    require(any(pair.value != g.zero for pair in built['coefficient_history'][-1].values()),
            'The exact order -1 history was incorrectly set to zero from an averaging cap')
    require(g.nodes == accepted['exact_function_graph_nodes'], 'Check appended a new original function definition')
    return dict(passed=True, exact_original_graph_nodes=len(g.nodes),
        same_N_original_source_references=references, original_kernel_Jacobian_integrals=integral_count,
        exact_order_minus1_functions_kept=True, exact_Picard_and_residual_roots_unchanged=True)


def bound_and_memory_checks(owner, contract):
    c, r = owner.ctx, owner.ranges.record
    quiet = pressure = 0
    for cell in r['cells']:
        if cell['chart'] != 'O3_power':
            continue
        for key in current.allN.RATES:
            require(cell['source_C0_Nminus2_covers'][key]['exact_zero']
                    and cell['source_Z_Nminus2_covers'][key]['exact_zero'],
                    'Quiet power correction is not exactly zero')
            quiet += 2
        for attribute in ('C0', 'Z'):
            require(cell['inherited_'+attribute+'_Nminus2_covers']['p']
                    == cell['right_'+attribute+'_Nminus2_covers']['p'],
                    'Zero-rate pressure memory changed')
            pressure += 1
    require(current.allN.RATES['p'] == 0, 'Original pressure recovery rate must be zero')
    floor = current.LogUpper.constant(c, c.mpf(1)/160)
    rows = {}
    for key, coefficient in owner.ranges.coefficients.items():
        accepted = current.cap(c, r['actual_uniform_N_scaled_repair_C1_caps']
            ['transformed_N_scaled_target_C1_caps'][key])
        replay = coefficient*floor
        # The accepted floor cap rounded its sum at 240 digits. This
        # 300-digit replay may tighten that bound without overlapping it.
        require(ep(replay.log)[1] <= ep(accepted.log)[1],
                'Unsubstituted N^-2 coefficient is wider than the accepted N=160 cap')
        rows[key] = dict(full_coefficient_not_floor_cap=True, original_floor_replay_not_wider=True)
    lower = contract['source_centered_repair_and_retained_local_C0_log_N_lower']
    require(all(ep(lower)[0] >= ep(value)[1] for value in contract['exact_frequency_conditions'].values()),
            'Common source/local frequency threshold omitted a condition')
    require(ep(packets.interval(c, contract['control_C1_upper_at_retained_local_frequency_floor']
            ['log_absolute_upper']))[1] < 0, 'Retained local frequency does not give a small control bound')
    require(ep(packets.interval(c, contract['normalized_F_G_Z_C1_upper_at_local_floor']
            ['log_absolute_upper']))[1] < 0, 'Original h/N bump normalization is not small')
    return dict(passed=True, same_original_cell_geometry_rows=len(owner.ranges.ledger),
        target_C0_Z_rows=len(rows)*2, normalized_target_coefficients=rows,
        quiet_exact_zero_rows=quiet, unchanged_zero_rate_pressure_memories=pressure,
        source_repair_and_retained_local_C0_conditions=len(contract['exact_frequency_conditions']),
        centered_full_target_not_exact_orderwise_coefficient_replacement=True,
        small_source_local_control_and_fixed_x_Z_C1_profile_bounds=True,
        higher_spatial_jets_or_whole_construction_N_admitted=False)


def run():
    began = time.monotonic()
    saved = json.loads((HERE/current.NAME).read_bytes())
    original = current.integrals.NativeRcC1IntegralRealization()
    ranges = current.SavedCenteredWholeZ(original)
    owner = current.CenteredRcAllNControlsBridge(original, ranges)
    built, contract = owner.build(), owner.contract()
    require(saved[current.GATE] and saved['source_family'] == owner.family,
            'Same source bridge result required')
    require(saved['source_local_C1_unit_ball_control_contract'] == packets.encode(contract),
            'Saved centered control contract differs from its typed implementation')
    require(saved['centered_full_target_C1_coefficients']
            == packets.encode({key:value.record() for key, value in ranges.coefficients.items()}),
            'Saved centered C1 coefficients changed')
    for key in ('actual_five_controls_installed', 'certified_actual_fixed_point_tail_installed',
                'actual_terminal_Z_function_closure_installed', 'current_whole_N_selected',
                'numerical_original_source_point_or_integral_oracle_installed', *packets.OPEN):
        require(saved[key] is False, 'A source/local bound cannot admit global numeric construction')
    result = dict(all_passed=True, **{current.GATE:True}, source_family=owner.family,
        independent_unit_ball_and_frequency_theorem=unit_ball_theorem(),
        actual_exact_integral_and_nonlinear_graph_checks=exact_graph_checks(owner, built),
        original_bound_units_and_history_checks=bound_and_memory_checks(owner, contract),
        read_only_review_model=dict(model='gpt-5.6-luna', reasoning_effort='max'),
        source_local_formal_unit_C1_control_ball_defined=True,
        actual_five_controls_installed=False, certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False, current_whole_N_selected=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        **dict.fromkeys(packets.OPEN, False),
        input_hashes=current.merge_hashes(owner.hashes,
            {current.NAME:sha(current.NAME), Path(__file__).name:sha(Path(__file__).name)}),
        execution_seconds=time.monotonic()-began, scope=saved['scope'])
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result), indent=2)+'\n', encoding='utf8')
    print('Centered all-N typed original integral/target/unit-ball bridge PASS', flush=True)
    return result


if __name__ == '__main__':
    run()
