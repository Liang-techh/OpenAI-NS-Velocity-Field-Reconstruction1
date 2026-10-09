"""Live signed bridge functions and source-owned zero-inlet C1 transport.

The three original micro/macro charts retain their actual source, positive
width, global phase origin, P0 and background. No ancestor support producer
is called and no interval cap is selected as a field value.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from types import FunctionType, SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_signed_functions as previous
import lei_ren_part1_paper_compliant_current_original_whole_Z_micro_macro as source_module
import lei_ren_part1_paper_compliant_current_original_whole_Z_R100_finite_N as r100

backend = previous.backend
switch, original, macro = source_module.switch, r100.original, source_module.macro
HERE, PREFIX, sha, bind, ep = previous.HERE, previous.PREFIX, previous.sha, previous.bind, previous.ep
CELLS, RATES, C0, Z, OPEN = previous.CELLS, previous.RATES, previous.C0, previous.Z, previous.OPEN
CHARTS = ('first_micro', 'second_micro', 'frozen_macro')
INLET = ('sc_half',)
NAME = PREFIX+'current_original_whole_Z_bridge_signed_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_bridge_signed_functions_check.json'
GATE = 'current_original_whole_Z_all_three_bridge_signed_functions_and_zero_inlet_C1_R100_transport_installed'
serialized, encode = original.serialized, previous.encode


def bounded_signed_exponent(value):
    """Materialize only the proved bounded A/N, with directed tiny tails.

    A broad negative log interval is legal. Its extreme negative exponential
    is enclosed by the original bounded_exp service, never expanded or set
    to zero. The original formal row and its sign remain unchanged.
    """
    if value.zero:
        return value.ctx.mpf(0)
    if ep(value.record()['log_absolute_upper'])[1] > 0:
        raise ArithmeticError('Actual signed source exponent must be proved within [-1,1]')
    result = value.coefficient*value.bounded_exp(value.scale.evaluate())
    if max(abs(x) for x in ep(result)) > 1:
        raise ArithmeticError('Directed bounded exponent extension exceeds admitted [-1,1]')
    return result


def original_density_callback():
    """Same original density/expm1 bodies; only bounded scalar evaluation."""
    module = backend.current.reference.phase.densities
    original_density, original_increment = module.density_Z_kernels, module.density.factored_expm1
    scalar_phase = SimpleNamespace(**{**vars(module.phase), 'bounded_value': bounded_signed_exponent})
    increment_phase = SimpleNamespace(**{**vars(module.density.phase), 'bounded_value': bounded_signed_exponent})
    increment = FunctionType(original_increment.__code__,
        {**original_increment.__globals__, 'phase': increment_phase},
        original_increment.__name__, original_increment.__defaults__)
    density = FunctionType(original_density.__code__,
        {**original_density.__globals__, 'phase': scalar_phase,
            'density': SimpleNamespace(**{**vars(module.density), 'factored_expm1': increment})},
        original_density.__name__, original_density.__defaults__)
    reference_phase = SimpleNamespace(**{**vars(backend.current.reference.phase),
        'densities': SimpleNamespace(**{**vars(module), 'density_Z_kernels': density})})
    reference = SimpleNamespace(**{**vars(backend.current.reference), 'phase': reference_phase})
    query = previous.previous._identity_query
    bound_query = FunctionType(query.__code__,
        {**query.__globals__, 'current': SimpleNamespace(**{**vars(backend.current), 'reference': reference})},
        query.__name__, query.__defaults__)
    proof = dict(original_density_and_factored_expm1_code_bodies_unchanged=True,
        only_bounded_scalar_exponent_callback_supplied=True,
        actual_signed_primitives_not_replaced_or_capped=True,
        actual_A_over_N_budget_checked_before_scalar_evaluation=True,
        directed_original_small_exponential_tails_retained=True,
        original_density_module=Path(module.__file__).name,
        original_density_module_sha256=sha(Path(module.__file__).name),
        original_expm1_module=Path(module.density.__file__).name,
        original_expm1_module_sha256=sha(Path(module.density.__file__).name))
    return bound_query, density, increment, proof


BRIDGE_QUERY, SIGNED_DENSITY, FACTORED_INCREMENT, DENSITY_BINDING = original_density_callback()


def function_digest(value):
    """Retain the original micro width-polynomial source evidence as data."""
    return hashlib.sha256(json.dumps(encode(serialized(value)), sort_keys=True,
        separators=(',', ':')).encode()).hexdigest()


def binary_fraction(value):
    """Read an exact singleton MPF without decimal rounding."""
    lo, hi = ep(value)
    if lo != hi or not mp.isfinite(lo):
        raise ValueError('Exact finite selected source singleton required')
    sign, mantissa, exponent, _ = lo._mpf_
    if abs(exponent) > 100000:
        raise ValueError('Compact source coordinate required')
    return Fraction((-1 if sign else 1)*mantissa)*Fraction(2)**exponent


def rational(value):
    return (value.numerator, value.denominator)


class CurrentBridgeRadiusPhase:
    def __init__(self, op, identity, N, sc):
        self.op, self.identity, self.N, self.sc = op, identity, N, sc
        self.base = backend.phase.RmRadiusPhase(op, identity['actual_five_defect_family_sha256'], identity)
        if self.base.sc._mpi_ != sc._mpi_ or op.flow.Y._mpi_ != (op.c.ln(100)-op.flow.logs[3])._mpi_:
            raise ValueError('Same original selected sc and exact macro Y required')
        self.start = sc/2
        self.cache = {}

    def query(self, chart, left, right=None):
        right = left if right is None else right
        def native(c, value):
            if value == INLET:
                if chart != 'first_micro':
                    raise ValueError('Original symbolic inlet belongs to first_micro')
                return c.mpf(self.sc)/2
            q = backend.exact(value)
            return c.mpf(q.numerator)/q.denominator
        l, r = native(self.op.c, left), native(self.op.c, right)
        lower, upper = (self.start, self.op.c.mpf(1)) if chart == 'first_micro' else \
            (self.op.c.mpf(1), self.op.c.mpf(2)) if chart == 'second_micro' else (self.op.c.mpf(0), self.op.c.mpf(1))
        if chart not in CHARTS or ep(l)[0] < ep(lower)[0] or ep(l)[0] > ep(r)[1] or ep(r)[1] > ep(upper)[1]:
            raise ValueError('Actual bridge coordinate in the original correction domain required')
        key = (chart, left, right)
        if key in self.cache:
            return self.cache[key]
        c = MPIntervalContext()
        c.dps = max(260, 80+(self.N.bit_length()*30103+99999)//100000+60)
        with mp.workdps(c.dps+40):
            cv = lambda q: c.mpf(q.numerator)/q.denominator
            l, r = native(c, left), native(c, right)
            t = c.mpf((ep(l)[0], ep(r)[1]))
            sc = c.mpf(self.sc)
            Y = c.ln(25)+1000+4*(c.exp(40)+11)
            regular = Y*t if chart == 'frozen_macro' else c.mpf(0)
            coefficient = 2*(1-t)-sc/2 if chart == 'frozen_macro' else t-sc/2
            projection = backend.signed.first.spatial.ordinary_mod_one(c, self.N*regular)
            epsilon = c.mpf(10)**-100
            lo, hi = ep(coefficient)
            maximum = max(abs(lo), abs(hi))
            logupper = None if not maximum else c.ln(self.N)+c.mpf(self.op.flow.logs[0])+c.ln(c.mpf(maximum))
            if logupper is not None and ep(logupper)[1] >= ep(c.ln(epsilon))[0]:
                raise ArithmeticError('Actual positive bridge width exceeds phase error budget')
            error = c.mpf(0) if not maximum else c.mpf((0, ep(epsilon)[1])) if lo >= 0 \
                else c.mpf((-ep(epsilon)[1], 0)) if hi <= 0 else c.mpf((-ep(epsilon)[1], ep(epsilon)[1]))
            projected = dict(full_period=True, boxes=[c.mpf((0,1))]) if projection['full_period'] else \
                backend.radius.periodic_add(c, projection['boxes'], error)
            full_image = chart == 'frozen_macro' and ep(r-l)[0] > 0 and ep(self.N*Y*(r-l)-2*epsilon)[0] >= 1
            if projection['full_period'] and not full_image:
                raise ArithmeticError('Full bridge period requires a proved actual radius image')
            result = dict(source_identity=self.identity, candidate_N=self.N, actual_chart=chart,
                exact_left=list(left), exact_right=list(right),
                actual_phase_boxes=[self.op.c.mpf(row) for row in projected['boxes']],
                full_period=projected['full_period'], full_period_proved_from_actual_N_physical_width=full_image,
                actual_phase_Z_exact_zero=True, actual_same_source_Rm_radius=self.op.Rm_factor,
                exact_same_positive_bridge_width=self.op.flow.h, actual_log_bridge_width=self.op.flow.logs[0],
                exact_source_macro_Y=self.op.flow.Y, selected_source_sc=self.sc,
                original_parameter_binding=self.base.parameter_binding,
                actual_source_logRa_and_analytic_logP_bound_by_original_Rm_binder=True,
                actual_N_bridge_width_phase_log_upper=logupper,
                signed_actual_positive_bridge_width_error=self.op.c.mpf(error),
                actual_positive_width_retained_not_zeroed=True, exact_inlet_phase_zero=left == right == INLET and chart == 'first_micro',
                adaptive_analytic_interval_digits=c.dps, global_phase_not_restarted=True,
                phase_not_supplied_as_free_angle=True,
                exact_original_global_phase='frac(N*(logR-logRa-hb*s_c/2))',
                exact_original_bridge_radius_maps=dict(first_micro='logRa+hb*s', second_micro='logRa+hb*s',
                    frozen_macro='logRa+2hb+(Y-2hb)*rho'))
        self.cache[key] = result
        return result


class WholeZBridgeSignedFunctions(previous.WholeZLongRmSignedFunctions):
    def __init__(self, dps=500):
        super().__init__(dps)
        checked = json.loads((HERE/previous.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(previous.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked actual current signed long/Rm backend required')
        for name, digest in checked['input_hashes'].items():
            bind(self.hashes, name, digest)
        bind(self.hashes, previous.RECEIPT, sha(previous.RECEIPT))
        self.bridge_background = self.switch_source.background
        if type(self.bridge_background) is not source_module.WholeZMicroMacro or self.bridge_background.identity != self.identity:
            raise ValueError('Same live original whole-Z bridge background required')
        self.inlet_proof = self.switch_source.saved['source_owned_left_inlet_proof']
        self.sc = source_module.read(self.c, self.inlet_proof['selected_positive_s_c'])
        self.start = self.sc/2
        if not 0 < ep(self.start)[0] < mp.mpf('.125'):
            raise ValueError('Original positive source sc below 1/4 required')
        if not all(self.inlet_proof[key] for key in (
                'exact_zero_correction_follows_from_defining_Duhamel_lower_bound',
                'actual_background_incoming_histories_remain_nonzero',
                'q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar',
                'no_saved_frame_or_N1024_initial_correction_used')):
            raise ValueError('Actual source-owned unchanged-left inlet proof required')
        for key, value in (('selected_eta_log', self.switch_source.eta_log), ('selected_dstar_log', self.switch_source.dstar_log)):
            if source_module.read(self.c, self.inlet_proof[key])._mpi_ != value._mpi_:
                raise ValueError('Exact source-owned eta/dstar inlet tuple required')
        excess = source_module.read(self.c, self.inlet_proof['original_left_kappa_excess'])
        if ep(excess)[0] <= 0 or ep(self.switch_source.eta_log)[1] > ep(self.c.ln(excess)-self.c.ln(2))[0]:
            raise ValueError('Strict original flat collar eta margin required')
        self.bridge_radius, self.bridge_source_cache = {}, {}
        for module in (previous, source_module, r100, switch, original, macro):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def coordinate(self, chart, value):
        if (value == 'inlet' or value == INLET) and chart == 'first_micro':
            return INLET
        q = backend.exact(value)
        lower, upper = (Fraction(0), Fraction(1)) if chart == 'first_micro' else \
            (Fraction(1), Fraction(2)) if chart == 'second_micro' else (Fraction(0), Fraction(1))
        if chart not in CHARTS or not lower <= q <= upper \
                or chart == 'first_micro' and ep(self.c.mpf(q.numerator)/q.denominator)[0] < ep(self.start)[0]:
            raise ValueError('Original bridge correction coordinate required')
        return rational(q)

    def bridge_owner(self, ends):
        op, source = self.owner(ends), self.bridge_background.owner(ends)
        if source.flow is not op.flow or source.c is not op.c or source.flow.logs[0]._mpi_ != self.bridge_background.logh._mpi_:
            raise ValueError('Identical bridge flow/context/positive width required')
        p0 = source.reference.P0
        if len(p0) != len(op.P0) or any(a.coefficient._mpi_ != b.coefficient._mpi_
                or a.scale.powers != b.scale.powers or a.scale.offset._mpi_ != b.scale.offset._mpi_
                or a.zero != b.zero or a.ctx is not b.ctx or a.scale.bases is not b.scale.bases
                or a.ledger is not b.ledger for a, b in zip(p0, op.P0)):
            raise ValueError('Exact source-equivalent original bridge P0 rows required')
        return source

    def native(self, chart, coordinate):
        if chart == 'frozen_macro':
            return coordinate
        if coordinate == INLET and chart == 'first_micro':
            return self.start
        q = backend.exact(coordinate)
        value = self.c.mpf(q.numerator)/q.denominator
        if ep(value)[0] != ep(value)[1] or binary_fraction(value) != q:
            raise ValueError('Exactly representable singleton original micro source coordinate required')
        return value

    def query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().query(ends, chart, left, right)
        l, r = self.coordinate(chart, left), self.coordinate(chart, left if right is None else right)
        result = BRIDGE_QUERY(self, ends, chart, l, r)
        result['actual_source_bounded_exponent_callback_binding'] = DENSITY_BINDING
        return result

    def source_query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().source_query(ends, chart, left, right)
        l, r = self.coordinate(chart, left), self.coordinate(chart, left if right is None else right)
        if chart == 'frozen_macro' and backend.exact(l) > backend.exact(r) \
                or chart != 'frozen_macro' and ep(self.native(chart, l))[0] > ep(self.native(chart, r))[0]:
            raise ValueError('Ordered original bridge cell required')
        key = (tuple(ends), chart, l, r)
        if key in self.bridge_source_cache:
            return self.bridge_source_cache[key]
        source, op = self.bridge_owner(ends), self.owner(ends)
        nl, nr = self.native(chart, l), self.native(chart, r)
        value = self.bridge_background.generic(ends, chart, nl, nr)
        recovered = value['full_signed_generic_source']
        roots, q, proof = switch.general_quotients(value['proxy'], recovered, value['a'], self.switch_source.eta_log)
        support = r100.WholeZR100FiniteN.support(SimpleNamespace(c=self.c, sc=self.sc), chart, nl, nr)
        diagnostic_q = q
        if support in ('unchanged_left', 'original_flat_left_collar'):
            q = {C0: op.flow.scalar(0), Z: op.flow.scalar(0)}
        packet = dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            chart=chart, exact_left=list(l), exact_right=list(r), exact_common_P0_axial5=op.P0,
            original_bridge_independent_P0_axial5=source.reference.P0,
            identical_source_P0_tuples_rebound_to_live_descendant_object=True,
            actual_background_source=value['background'], original_raw_source=value['original_raw_source'],
            original_generic_source=recovered, original_roots=roots['roots'], original_q_C0_Z=q,
            original_full_source_quotients=proof, original_generic_q_diagnostic_before_source_collar_selector=diagnostic_q,
            source_owned_modification_support=support, actual_source_owned_initial_condition=self.inlet_proof,
            source_owned_collar_q_C0_Z_exact_zero=support in ('unchanged_left', 'original_flat_left_collar'),
            original_active_inverse_not_evaluated_on_source_owned_flat_collar=True,
            original_generic_roots_preserved_without_all_u_support_evaluation=True,
            actual_phase_Z_exact_zero=True)
        self.bridge_source_cache[key] = packet
        return packet

    def radius_query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().radius_query(ends, chart, left, right)
        l, r = self.coordinate(chart, left), self.coordinate(chart, left if right is None else right)
        self.bridge_owner(ends)
        key = tuple(ends)
        if key not in self.bridge_radius:
            self.bridge_radius[key] = CurrentBridgeRadiusPhase(self.owner(ends), self.identity, self.N, self.sc)
        return self.bridge_radius[key].query(chart, l, r)

    def partitions(self):
        return (('first_micro', (INLET, (1,4), (1,2), (3,4), (1,1)), (1,1)),
                ('second_micro', ((1,1), (5,4), (3,2), (7,4), (2,1)), (2,1)),
                ('frozen_macro', ((0,1), (1,2), (1,1)), None))

    def weights(self, ends, chart, left, right, rate):
        source = self.bridge_owner(ends)
        l, r = self.coordinate(chart, left), self.coordinate(chart, right)
        if chart == 'frozen_macro' and backend.exact(l) >= backend.exact(r) \
                or chart != 'frozen_macro' and ep(self.native(chart, l))[0] >= ep(self.native(chart, r))[0]:
            raise ValueError('Positive original bridge integral cell required')
        if chart == 'frozen_macro':
            return macro.weights(source.series, l, r, rate)
        return original.weights(source.series, self.native(chart, l), self.native(chart, r),
            self.c.mpf(1 if chart == 'first_micro' else 2), rate)

    def local_integral(self, ends, chart, left, right):
        if chart not in CHARTS:
            return super().local_integral(ends, chart, left, right)
        l, r = self.coordinate(chart, left), self.coordinate(chart, right)
        if chart == 'frozen_macro' and backend.exact(l) >= backend.exact(r) \
                or chart != 'frozen_macro' and ep(self.native(chart, l))[0] >= ep(self.native(chart, r))[0]:
            raise ValueError('Positive original bridge integral cell required')
        source, op = self.query(ends, chart, l, r), self.owner(ends)
        f, c = op.flow, op.c
        length = f.scalar(f.Y)-f.h*2 if chart == 'frozen_macro' else f.h
        if chart == 'frozen_macro':
            d = backend.exact(r)-backend.exact(l)
            ds = c.mpf(d.numerator)/d.denominator
        else:
            ds = self.native(chart, r)-self.native(chart, l)
        density, integrals, masses = source['actual_signed_nonlinear_density_C0_Z'], {}, {}
        for name, rate in RATES.items():
            mass, _, _ = self.weights(ends, chart, l, r, rate)
            integrals[name] = [density[part][name]*mass for part in ('kernels', 'Z_derivatives')]
            masses[name] = mass
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            actual_chart=chart, exact_left=list(l), exact_right=list(r), exact_common_P0_axial5=op.P0,
            actual_source_function=source, actual_original_window_length=length,
            actual_physical_log_width=length*ds,
            actual_positive_own_rate_masses=masses, actual_signed_local_integral_C0_Z=integrals,
            integral_recipe='int_left^right exp(-rate*L*(right-s))*signed_density(s,Z,N)*L ds',
            original_hb_hb_Y_minus_2hb_measure_and_cutoff_mass_preserved=True,
            physical_Jacobian_applied_once=True, actual_positive_micro_width_not_zeroed=True,
            live_inverse_function_extension_under_Z_independent_integral=True,
            directed_rectangle_function_enclosure_not_selected_integral_value=True,
            incoming_correction_not_supplied_or_reset=True, **dict.fromkeys(OPEN, False),
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False)

    def bridge_transport(self, ends):
        op = self.owner(ends)
        f = op.flow
        inlet = self.query(ends, 'first_micro', 'inlet')
        if not inlet['actual_original_source_packet']['source_owned_collar_q_C0_Z_exact_zero'] \
                or not all(row.zero for part in ('kernels', 'Z_derivatives')
                    for row in inlet['actual_signed_nonlinear_density_C0_Z'][part].values()):
            raise ValueError('Defining unchanged-left source must prove exact zero correction inlet')
        initial = {name: [f.scalar(0), f.scalar(0)] for name in RATES}
        incoming = {name: list(rows) for name, rows in initial.items()}
        windows = []
        for chart, partition, _ in self.partitions():
            local, cells = {name: [f.scalar(0), f.scalar(0)] for name in RATES}, []
            for left, right in zip(partition, partition[1:]):
                integral = self.local_integral(ends, chart, left, right)
                contributions, weights = {}, {}
                for name, rate in RATES.items():
                    mass, decay, suffix = self.weights(ends, chart, left, right, rate)
                    contributions[name] = [row*suffix for row in integral['actual_signed_local_integral_C0_Z'][name]]
                    for i in range(2):
                        local[name][i] += contributions[name][i]
                    weights[name] = dict(original_full_mass=mass, original_cell_decay=decay, original_suffix_decay=suffix)
                query = integral['actual_source_function']
                cells.append(dict(exact_left=list(left), exact_right=list(right),
                    current_signed_function_and_integral_sha256=function_digest(integral),
                    actual_signed_primitive_C0_Z_phi=query['actual_signed_primitive_C0_Z_phi'],
                    actual_signed_nonlinear_density_C0_Z=query['actual_signed_nonlinear_density_C0_Z'],
                    actual_signed_local_integral_C0_Z=integral['actual_signed_local_integral_C0_Z'],
                    actual_signed_contribution_at_chart_exit_C0_Z=contributions,
                    original_own_rate_weights=weights, actual_physical_log_width=integral['actual_physical_log_width'],
                    actual_global_phase_full_period_cover=query['actual_source_geometry']['full_period']))
            memory = {name: self.weights(ends, chart, partition[0], partition[-1], rate)[1] for name, rate in RATES.items()}
            inherited = {name: [row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing = {name: [inherited[name][i]+local[name][i] for i in range(2)] for name in RATES}
            windows.append(dict(actual_chart=chart, actual_signed_cells=cells, actual_incoming_correction_C0_Z=incoming,
                original_incoming_own_rate_memory=memory, actual_signed_local_driver_C0_Z=local,
                actual_retained_incoming_C0_Z=inherited, actual_exit_signed_correction_C0_Z=outgoing))
            incoming = outgoing
            print('Current signed bridge C1 transport: '+str(ends)+' '+chart, flush=True)
        terminal = self.query(ends, 'frozen_macro', (1,1))
        background = {name: list(rows[:2]) for name, rows in terminal['actual_original_source_packet']
            ['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete = {name: [background[name][i]+incoming[name][i] for i in range(2)] for name in RATES}
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=op.P0, actual_source_owned_initial_condition=self.inlet_proof,
            actual_source_owned_zero_inlet_signed_function=inlet,
            actual_current_inlet_correction_C0_Z=initial, actual_signed_bridge_windows=windows,
            actual_R100_signed_correction_C0_Z=incoming, actual_original_R100_background_C0_Z=background,
            actual_R100_signed_complete_own_history_C0_Z=complete, actual_original_R100_signed_source_function=terminal,
            exact_R100_radius=f.scalar(100), exact_micro_macro_source_join=self.bridge_background.seam(ends),
            zero_inlet_follows_from_source_definition_not_arbitrary_incoming_reset=True,
            original_background_nonzero_histories_remain_separate=True,
            interval_caps_not_substituted_for_defining_functions=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZBridgeSignedFunctions()
        points, transports = [], []
        for ends in CELLS:
            for chart, point in (('first_micro', (1,2)), ('second_micro', (3,2)), ('frozen_macro', (1,3))):
                points.append(serialized(owner.query(ends, chart, point)))
            transports.append(serialized(owner.bridge_transport(ends)))
        report = dict(**{GATE: True}, source_family=owner.identity, candidate_N=owner.N,
            exact_Z_partition=CELLS, current_bridge_charts=CHARTS, exact_original_native_partitions=owner.partitions(),
            actual_signed_bridge_point_functions=points, actual_signed_bridge_transports=transports,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='All three original whole-Z bridge signed inverse C0/Z/phase source functions and '
                '40 true-measure C1 integral cells from the source-defined zero correction inlet through R100. '
                'All 17 chart functions are now available; full inlet-to-Rc integral wiring, actual Rc '
                'defects/controls/terminal identities, global N/heat/cone and temporal recursion remain unfinished.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report), separators=(',',':'))+'\n').encode(),
            compresslevel=6, mtime=0))
    print('Current whole-Z bridge signed functions and source-owned zero-inlet C1 transport reach R100', flush=True)
    return report


if __name__ == '__main__':
    run()
