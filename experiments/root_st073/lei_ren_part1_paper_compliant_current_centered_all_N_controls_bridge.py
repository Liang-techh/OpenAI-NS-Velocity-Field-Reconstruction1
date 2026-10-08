"""Bind centered all-N target bounds to the original realized C1 control graph.

The N^-2 cover bounds the full target, not the individual N^-1 coefficient.
Exact integrals, the common spatial N and P0 remain unchanged. No point
values are selected from caps.
"""
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_centered_phase_conditioning as centered
import lei_ren_part1_paper_compliant_current_native_Rc_convergent_control_family as family

integrals, allN, controls, source, packets = family.integrals, family.allN, family.controls, family.source, family.packets
HERE, PREFIX, sha, ep = family.HERE, family.PREFIX, family.sha, family.ep
LogUpper, require = allN.repair.LogUpper, integrals.require
NAME = PREFIX+'current_centered_all_N_controls_bridge.json'
RECEIPT = PREFIX+'current_centered_all_N_controls_bridge_check.json'
GATE = 'current_original_centered_whole_Z_all_N_target_bound_and_unit_ball_control_bridge_bound'


def merge_hashes(*groups):
    result = {}
    for group in groups:
        for name, digest in group.items():
            require(name not in result or result[name] == digest, 'Different original source bytes: '+name)
            result[name] = digest
    return result


def digest_record(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def cap(c, record):
    return LogUpper(c, None if record['exact_zero'] else packets.interval(c, record['log_absolute_upper']))


class SavedCenteredWholeZ:
    """Typed checked full-target range sidecar; not a source evaluator."""
    def __init__(self, original):
        if type(original) is not integrals.NativeRcC1IntegralRealization:
            raise TypeError('Accepted original C1 integral owner required')
        self.original, self.ctx, self.family = original, original.ctx, original.family
        checked = json.loads((HERE/centered.RECEIPT).read_bytes())
        data = json.loads((HERE/centered.NAME).read_bytes())
        require(checked['all_passed'] and checked[centered.GATE] and data[centered.GATE]
                and checked['source_family'] == data['source_family'] == self.family,
                'Checked centered same-family whole-source certificate required')
        require(checked['input_hashes'].get(centered.NAME) == sha(centered.NAME),
                'Centered data bytes differ from their receipt')
        self.hashes = merge_hashes(original.hashes, checked['input_hashes'], data['input_hashes'],
            {centered.NAME:sha(centered.NAME), centered.RECEIPT:sha(centered.RECEIPT)})
        self.record = data['actual_native_centered_phase_conditioning_records']['whole_Z']
        r = self.record
        require(packets.interval(self.ctx, r['Z_box'])._mpi_ == self.ctx.mpf((-1, 1))._mpi_
                and r['all_integer_N_lower'] == 160 and r['original_true_cells'] == 24
                and r['original_charts'] == 17 and len(r['cells']) == 24,
                'Whole-Z all-N/24-cell/17-chart range sidecar required')
        require(r['original_all_N_zero_inlet_and_phase_and_P0_retained']
                and r['native_normalization_and_ordinary_width_applied_once'],
                'Original inlet, phase, P0 and source units must be retained')
        old = original.data['current_native_Rc_all_N_function_controls']
        graph_name = source.sources.VIEWS
        require(self.hashes[graph_name] == old['input_hashes'][graph_name] == sha(graph_name),
                'Centered and exact C1 functions must use the same original graph')
        require(old['original_conditional_coefficient_branch_count'] == 35,
                'Original branchwise nonlinear density certificate required')
        self.ledger = []
        for cell, direct in zip(r['cells'], old['actual_all_N_continuous_cell_range_records']):
            require((cell['label'], cell['chart']) == (direct['label'], direct['chart'])
                    and cell['original_geometry'] == direct['original_geometry'],
                    'Centered range and original integral have different cell geometry')
            require(cell['true_mass_once_and_all_endpoints_preserved']
                    and cell['source_increment_functions_unchanged'],
                    'Endpoint, width or source function changed')
            if cell['label'] != 'initial_flat_collar':
                proof = cell['original_source']['native_periodic_mixed_proof']
                if not proof['original_whole_support_flat']:
                    require(proof['native_active_kappa_correlation']['active_support_only']
                            and proof['joint_weighted_kernel_proof']['no_q_or_chi_division_used']
                            and proof['centered_original_phase_proof']['a_nu_equals_v_used_before_bounding']
                            and proof['centered_original_phase_proof']['complete_mixed_phase_and_inverse_cross_terms_retained'],
                            'Centered active-support or complete inverse proof lost')
                require(cell['original_source']['formal_mixed_source_definitions_unchanged']
                        and cell['original_source']['source_point_values_not_selected_from_caps'],
                        'Range sidecar cannot define new function values')
            self.ledger.append(dict(label=cell['label'], chart=cell['chart'],
                original_geometry_sha256=digest_record(cell['original_geometry']),
                same_original_true_width_endpoints_and_normalization=True))
        for key in controls.ROWS:
            for rows in (r['target_C0_orders'], r['target_Z_orders']):
                require(rows[key]['-1']['exact_zero'], 'Averaged full-target cover must have only N^-2 order')
                require(rows[key]['-2']['encloses_original_source_function']
                        and not rows[key]['-2']['point_value_selected'], 'Source cap used as a point value')
        self.coefficients = {
            key:LogUpper.add(self.ctx, [cap(self.ctx, rows[key]['-2'])
                for rows in (r['target_C0_orders'], r['target_Z_orders'])])
            for key in controls.ROWS}
        nonzero = [value.log for value in self.coefficients.values() if value.log is not None]
        self.whole_coefficient = LogUpper(self.ctx, self.ctx.mpf(max(ep(v)[1] for v in nonzero)))
        self.theorem = data['exact_correlated_periodic_kernel_theorem']


def unit_ball_conditions(c, coefficient, old, source_requirements, local_logN):
    """Use ||N*r||C1<=C/N, rather than its uniform value at N=160."""
    read = lambda key:cap(c, old[key])
    CA, CQ, gmax = [read(key) for key in ('exact_integral_matrix_inverse_log_cap',
        'general_transformed_quadratic_C1_log_cap', 'compact_bump_profile_log_upper')]
    require(old['B_mu_Z_independent_and_C1_product_norm_submultiplicative'],
            'Original Z-independent exact inverse and quadratic C1 norm required')
    constant = lambda v:LogUpper.constant(c, v)
    first = CA*coefficient*constant(2)
    contraction = CA*CQ*constant(4)
    positivity = gmax*constant(2)*LogUpper(c, c.mpf(2)/3*c.ln(2))
    entries = dict(exponential_safety=c.ln(160), unit_ball_first_step=first.log,
        unit_ball_half_contraction=contraction.log, original_swirl_positivity=positivity.log)
    entries.update({'original_source_'+key:value for key, value in source_requirements.items()})
    repair_lower = c.mpf(max(ep(value)[1] for value in entries.values()))
    entries['accepted_original_q_flat_active_quiet_and_band_C0'] = local_logN
    local_lower = c.mpf(max(ep(value)[1] for value in entries.values()))
    small_control = first.divide_positive(local_lower)
    small_profile = gmax*small_control.divide_positive(local_lower)
    return dict(whole_Z_domain=[-1, 1], all_integer_N_floor=160,
        original_exact_inverse_C1_cap=CA.record(), original_quadratic_C1_cap=CQ.record(),
        original_bump_profile_cap=gmax.record(), averaged_full_target_C1_coefficient=coefficient.record(),
        target_bound='norm_C1(N*r(N,.))<=C/N, by the original whole-source N^-2 averaging theorem',
        control_ball_radius=1, first_Picard_C1_norm_upper='1/2',
        image_C1_radius_upper='3/4', Lipschitz_upper='1/2',
        exact_frequency_conditions=entries, source_and_centered_repair_log_N_lower=repair_lower,
        source_centered_repair_and_retained_local_C0_log_N_lower=local_lower,
        finite_integer_definition='N=ceil(exp(L)); L is the recorded directed local log-N threshold',
        integer_or_exponential_not_materialized=True,
        unit_ball_contained_in_prior_local_cone_control_ball=True,
        retained_local_C0_budget_reoptimized=False,
        exact_same_source_C1_limit='h_star=lim_C1 h_n; h0=0; h_next=-B(mu)^-1*(N*r+Q(mu,h_n)/N)',
        frequency_dependent_control_C1_upper='min(1,2*CA*C/N)',
        control_C1_upper_at_retained_local_frequency_floor=small_control.record(),
        normalized_F_G_Z_C1_upper_at_local_floor=small_profile.record(),
        small_profile_formula='g_max*min(1,2*CA*C/N)/N; C1 in Z at fixed repair x, using h/N in the original compact bumps',
        higher_spatial_profile_derivatives_not_bounded=True,
        implicit_Z_identity='(B+DQ(h_star)/N)*h_star_Z=-N*r_Z',
        mathematical_iteration_tail='norm_C1(h_star-h_n)<=2^-n',
        numerical_quadrature_or_oracle_errors_included=False,
        source_local_conditional_C1_unit_ball_family_defined=True,
        actual_global_integer_N_selected=False)


class CenteredRcAllNControlsBridge:
    """Original C1 integral owner plus a checked bound on its complete target."""
    def __init__(self, original, ranges):
        if type(original) is not integrals.NativeRcC1IntegralRealization or type(ranges) is not SavedCenteredWholeZ:
            raise TypeError('Typed original C1 integral owner and checked centered sidecar required')
        if ranges.original is not original or ranges.ctx is not original.ctx or ranges.family != original.family:
            raise ValueError('Same exact original owner, context and source family required')
        self.original, self.ranges, self.ctx, self.family = original, ranges, original.ctx, original.family
        self.hashes = merge_hashes(original.hashes, ranges.hashes,
            {Path(__file__).name:sha(Path(__file__).name)})
        for module in (integrals, family, integrals.previous):
            checked = json.loads((HERE/module.RECEIPT).read_bytes())
            require(checked['all_passed'] and checked[module.GATE] and checked['source_family'] == self.family,
                    'Checked original C1 realization/control/local budget required')
            self.hashes = merge_hashes(self.hashes, checked['input_hashes'],
                {module.NAME:sha(module.NAME), module.RECEIPT:sha(module.RECEIPT)})
        self.old = original.data['current_native_Rc_all_N_function_controls']
        self.local = json.loads((HERE/integrals.previous.NAME).read_bytes())
        old_conditions = self.old['actual_uniform_repair_C1_log_N_conditions']
        require(ep(cap(self.ctx, old_conditions['formal_control_C1_ball_radius_log']).log)[0] >= 0
                and old_conditions['B_mu_Z_independent_and_C1_product_norm_submultiplicative'],
                'Unit C1 ball must lie in the normalized ball used by retained local cone bounds')
        self.hashes[Path(family.evaluator.__file__).name] = sha(Path(family.evaluator.__file__).name)

    def build(self):
        built = self.original.build()
        original_nodes = list(built['graph'].nodes)
        built = family.generic_control_operator(built)
        require(built['graph'].nodes[:len(original_nodes)] == original_nodes,
                'Centering changed original function/integral recipes')
        return built

    def contract(self):
        c = self.ctx
        requirements = {chart:packets.interval(c, value) for chart, value in self.old[
            'original_source_and_repair_combined_log_conditions']['original_generic_source_log_N_requirements'].items()}
        require(len(requirements) == 17, 'All original source-frequency requirements required')
        return unit_ball_conditions(c, self.ranges.whole_coefficient,
            self.ranges.record['conditional_repair_log_conditions'], requirements,
            packets.interval(c, self.local['source_repair_active_flat_quiet_band_C0_sufficient_log_N_lower']))

    def evaluator(self, built, *, oracle, Z, N, ctx):
        """Actual supplied oracle required; no range-to-point fallback."""
        return family.CachedC1ControlEvaluator(built, oracle=oracle, Z=Z, N=N, ctx=ctx,
            original_log_N_required=self.contract()['source_centered_repair_and_retained_local_C0_log_N_lower'])

    def factored_tail(self, depth):
        if type(depth) is not int or not 0 <= depth <= 4096:
            raise ValueError('Exact finite iteration depth in[0,4096] required')
        return dict(depth=depth, control_C1_tail=dict(numerator=1, denominator=2**depth),
            preconditioned_residual_C1_tail=dict(numerator=3, denominator=2**(depth+1)),
            mathematical_tail_only=True, original_numeric_fixed_point_certified=False)


def run():
    began = time.monotonic()
    original = integrals.NativeRcC1IntegralRealization()
    ranges = SavedCenteredWholeZ(original)
    owner = CenteredRcAllNControlsBridge(original, ranges)
    built, contract = owner.build(), owner.contract()
    row = controls.pair_roots
    records = dict(source_family=owner.family, **{GATE:True},
        exact_original_all_N_graph=allN.NAME, exact_original_integral_receipt=integrals.RECEIPT,
        checked_centered_whole_Z_receipt=centered.RECEIPT,
        same_original_source_graph_sha256=built['source_graph_sha256'],
        original_integral_count=len(built['exact_integral_nodes']), original_radial_cells=len(built['cells']),
        same_original_24_cell_geometry_bindings=ranges.ledger,
        exact_N_scaled_target_roots=row(built['N_scaled_targets'].values(), controls.ROWS),
        exact_generic_Picard_C1_map_roots=row(built['generic_Picard_map'], controls.CONTROLS),
        exact_control_residual_roots=row(built['control_residual'], controls.ROWS),
        exact_shared_N_node=built['N'].node,
        exact_function_graph_nodes_unchanged_from_accepted_control_family=True,
        centered_full_target_C1_coefficients={key:value.record() for key, value in ranges.coefficients.items()},
        source_local_C1_unit_ball_control_contract=contract,
        factored_unit_ball_tail_examples=[owner.factored_tail(n) for n in (0, 3, 16, 64)],
        prior_all_N_repair_log_N_lower=owner.old['actual_uniform_repair_C1_log_N_conditions']['repair_sufficient_common_log_N_lower'],
        original_N_dependent_order_minus1_coefficient_functions_retained=True,
        averaging_cover_not_a_zero_coefficient_function_identity=True,
        exact_actual_phase='fractional_part(N*original_log_radius_minus_inlet)',
        original_P0_P0_Z_unchanged_and_not_added_to_correction_histories=True,
        original_source_ancestor_constructors_called=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False, certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False, current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN, False), input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Typed same-source centered whole-Z/all-N averaged target sidecar now drives original realized C1 integrals, exact Picard/residual graph, source/local frequency inequalities and a formal unit C1 ball. Original N-dependent coefficient functions unchanged. No cap-valued oracle, global N, numeric control/terminal field, heat/stress recursion or corrected uvw.')
    (HERE/NAME).write_text(json.dumps(packets.encode(records), indent=2)+'\n', encoding='utf8')
    print('Centered all-N original C1 graph/target bridge and unit control ball bound', flush=True)
    return records


if __name__ == '__main__':
    run()
