"""Current R100/R110 signed inverse functions and genuine C1 transport.

All three native coordinates are in [0,1]. Original positive microscopic
widths and the common phase origin remain; no width or incoming is reset.
"""
import gzip
import json
from pathlib import Path
from types import FunctionType
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_Rh_signed_functions as previous

backend = previous.previous
source_module = previous.patch.prefix
switch = source_module.upstream
original = switch.switch
HERE, PREFIX, sha, bind, ep = previous.HERE, previous.PREFIX, previous.sha, previous.bind, previous.ep
CELLS, RATES, C0, Z, OPEN = previous.CELLS, previous.RATES, previous.C0, previous.Z, previous.OPEN
CHARTS = switch.CHARTS
PARTITION = switch.PARTITION
NAME = PREFIX+'current_original_whole_Z_switch_signed_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_switch_signed_functions_check.json'
GATE = 'current_original_whole_Z_R100_R110_live_signed_functions_and_C1_transport_installed'
serialized, encode, function_digest = previous.serialized, previous.encode, previous.function_digest


def inverse_identity_branch(source, qrows, dstar_log, phi, branch):
    """Evaluate the same A through its exact, conditioned inverse identity."""
    got = backend.branch_Z_functions(source, qrows, dstar_log, phi, branch)
    if got['values'] is None or got['record']['geometry'] == 'flat':
        return got
    inverse = got['record']['original_C0_inverse']
    psi = inverse['psi_fraction_interval']
    loop = got['loop']
    old = got['values']['A']
    value = source['roots']['a'][C0]*loop.scalar(loop.c.mpf(phi)-psi)*loop.c.mpf('.5')
    got['values']['A'] = value
    got['record'].update(original_A_C0_before_conditioned_inverse_identity=old.record(),
        actual_inverse_identity_A_C0=value.record(),
        exact_original_inverse_identity='A=a/2*(phi-psi_fraction)',
        actual_selected_inverse_psi_used_not_support_cap=True,
        original_Z_phi_derivative_enclosures_retained=True,
        original_defining_inverse_primitive_function_unchanged=True)
    return got


def directed_dominating_union(rows):
    """Union in a proved dominating arithmetic coordinate, never a field unit."""
    if not rows:
        raise ValueError('Nonempty signed function alternatives required')
    first, c = rows[0], rows[0].ctx
    for row in rows:
        first.coerce(row)
    if len(rows) == 1:
        return first
    prior = backend.signed.density.local.prior
    reference = max(rows, key=lambda row: ep(row.scale.evaluate())[1]).scale
    if any(ep((row.scale-reference).evaluate())[1] > 0 for row in rows if not row.zero):
        upper = max(ep(row.scale.evaluate())[1] for row in rows if not row.zero)
        padding = c.mpf(max(mp.mpf(1), abs(upper)))*c.mpf(2)**(-c.prec//2)
        for _ in range(8):
            offset = c.mpf(ep(c.mpf(upper)+padding)[1])
            reference = prior.FormalScale(first.scale.bases, offset=offset)
            if all(ep((row.scale-reference).evaluate())[1] <= 0 for row in rows if not row.zero):
                break
            padding *= 2
        else:
            raise ArithmeticError('Actual signed function union requires scale refinement')
        first.ledger['directed_independent_log_rescalings'] += 1
    coefficients = [c.mpf(0) if row.zero else
        row.coefficient*row.bounded_exp((row.scale-reference).evaluate()) for row in rows]
    return prior.ScaledEnclosure(reference, c.mpf((min(ep(row)[0] for row in coefficients),
        max(ep(row)[1] for row in coefficients))), first.ledger)


_identity_query = FunctionType(backend.WholeZSignedLoopFunctions.query.__code__,
    {**backend.WholeZSignedLoopFunctions.query.__globals__, 'branch_Z_functions': inverse_identity_branch,
        'union': directed_dominating_union},
    'current_signed_inverse_identity_query', backend.WholeZSignedLoopFunctions.query.__defaults__)


class ZeroLengthCutoffMass:
    """The same source provider; its empty partial integral is exactly zero."""
    def __init__(self, first):
        self.first = first

    def __getattr__(self, name):
        return getattr(self.first, name)

    def mass(self, left, right):
        return self.first.c.mpf(0) if left == right else self.first.mass(left, right)


class CurrentSwitchRadiusPhase:
    def __init__(self, op, identity, N):
        self.op, self.identity, self.N = op, identity, N
        self.base = backend.phase.RmRadiusPhase(op, identity['actual_five_defect_family_sha256'], identity)
        self.cache = {}

    def query(self, chart, left, right=None):
        right = left if right is None else right
        l, r = backend.exact(left), backend.exact(right)
        if chart not in CHARTS or not 0 <= l <= r <= 1:
            raise ValueError('Actual switch native coordinate in [0,1] required')
        key = (chart, l, r)
        if key in self.cache:
            return self.cache[key]
        c = MPIntervalContext()
        c.dps = max(260, 80+(self.N.bit_length()*30103+99999)//100000+60)
        with mp.workdps(c.dps+40):
            cv = lambda q: c.mpf(q.numerator)/q.denominator
            t = c.mpf((ep(cv(l))[0], ep(cv(r))[1]))
            sc = c.mpf(self.base.sc)
            regular = c.ln(25)+1000+4*(c.exp(40)+11)
            coefficient = t-sc/2 if chart == 'first_switch' else 1+t-sc/2
            if chart == 'post_power':
                regular += c.ln(c.mpf(11)/10)*t
                coefficient = 2*(1-t)-sc/2
            projection = backend.signed.first.spatial.ordinary_mod_one(c, self.N*regular)
            epsilon = c.mpf(10)**-100
            lo, hi = ep(coefficient)
            maximum = max(abs(lo), abs(hi))
            logupper = None if not maximum else c.ln(self.N)+c.mpf(self.op.flow.logs[0])+c.ln(c.mpf(maximum))
            if logupper is not None and ep(logupper)[1] >= ep(c.ln(epsilon))[0]:
                raise ArithmeticError('Actual positive switch width exceeds phase error budget')
            error = c.mpf(0) if not maximum else c.mpf((0, ep(epsilon)[1])) if lo >= 0 \
                else c.mpf((-ep(epsilon)[1], 0)) if hi <= 0 else c.mpf((-ep(epsilon)[1], ep(epsilon)[1]))
            projected = dict(full_period=True, boxes=[c.mpf((0,1))]) if projection['full_period'] else \
                backend.radius.periodic_add(c, projection['boxes'], error)
            result = dict(source_identity=self.identity, candidate_N=self.N, actual_chart=chart,
                exact_left=list(left), exact_right=list(right),
                actual_phase_boxes=[self.op.c.mpf(row) for row in projected['boxes']],
                full_period=projected['full_period'], actual_phase_Z_exact_zero=True,
                actual_same_source_Rm_radius=self.op.Rm_factor,
                exact_same_positive_switch_width=self.op.flow.h,
                actual_log_switch_width=self.op.flow.logs[0], selected_source_sc=self.base.sc,
                original_parameter_binding=self.base.parameter_binding,
                actual_source_logRa_and_analytic_logP_collection_bound_by_original_Rm_binder=True,
                adaptive_analytic_interval_digits=c.dps,
                original_regular_phase_component=regular, original_regular_periodic_projection=projection,
                original_microscopic_width_coefficient=coefficient,
                actual_N_width_coefficient_log_upper=logupper,
                signed_actual_microscopic_phase_error=self.op.c.mpf(error),
                actual_positive_width_retained_not_zeroed=True,
                global_phase_not_restarted=True, phase_not_supplied_as_free_angle=True,
                exact_original_global_phase='frac(N*(log100+physical_log_radius_over100-logRa-hb*s_c/2))',
                original_radius_maps=dict(first_switch='log100+h*s',
                    second_switch='log100+h*(1+s)', post_power='log100+2*h+(log(11/10)-2*h)*t'),
                physical_log_radius_derivative='h' if chart != 'post_power' else 'log(11/10)-2*h')
        self.cache[key] = result
        return result


class WholeZSwitchSignedFunctions(previous.WholeZRmRhSignedFunctions):
    def __init__(self, dps=500):
        super().__init__(dps)
        checked = json.loads((HERE/previous.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(previous.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked current signed source/patch C1 backend required')
        for name, digest in checked['input_hashes'].items():
            bind(self.hashes, name, digest)
        bind(self.hashes, previous.RECEIPT, sha(previous.RECEIPT))
        self.switch_source = self.patch_source.source.source.upstream.source
        if type(self.switch_source) is not switch.WholeZSwitchFiniteN:
            raise ValueError('Genuine current whole-Z switch owner chain required')
        if self.switch_source.N != self.N or self.switch_source.identity != self.identity:
            raise ValueError('Actual switch source family and unchanged N required')
        native_prefix = self.switch_source.prefix
        self.switch_source.prefix = lambda first, second, chart, left, right: native_prefix(
            ZeroLengthCutoffMass(first), second, chart, left, right)
        self.switch_radius = {}
        for module in (previous, switch, original):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        primitive_module = switch.primitives
        name = Path(primitive_module.__file__).name
        bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def switch_owner(self, ends):
        op = self.owner(ends)
        source = self.switch_source.owner(ends)
        if source.flow is not op.flow or source.c is not op.c:
            raise ValueError('Identical live switch flow/context required; flow='+str(source.flow is op.flow)
                +'; context='+str(source.c is op.c))
        source_P0 = source.background.reference.P0
        if len(source_P0) != len(op.P0) or any(a.coefficient._mpi_ != b.coefficient._mpi_
                or a.scale.powers != b.scale.powers or a.scale.offset._mpi_ != b.scale.offset._mpi_
                or a.zero != b.zero or a.ctx is not b.ctx or a.scale.bases is not b.scale.bases
                or a.ledger is not b.ledger for a, b in zip(source_P0, op.P0)):
            raise ValueError('Exact source-equivalent independent P0 rows required')
        if source.flow.logs[0]._mpi_ != self.switch_source.background.logh._mpi_:
            raise ValueError('Same original positive microscopic switch width required')
        return source

    def query(self, ends, chart, left, right=None):
        if chart in CHARTS:
            return _identity_query(self, ends, chart, left, right)
        return super().query(ends, chart, left, right)

    def source_query(self, ends, chart, left, right=None):
        if chart in CHARTS:
            self.switch_owner(ends)
            packet = dict(self.switch_source.query(ends, chart, left, right))
            packet['original_switch_independent_P0_axial5'] = packet['exact_common_P0_axial5']
            packet['exact_common_P0_axial5'] = self.owner(ends).P0
            packet['identical_source_P0_tuples_rebound_to_live_descendant_object'] = True
            packet['empty_partial_cutoff_integral_exact_zero_nonempty_recipe_unchanged'] = True
            return packet
        return super().source_query(ends, chart, left, right)

    def radius_query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().radius_query(ends, chart, left, right)
        self.switch_owner(ends)
        key = tuple(ends)
        if key not in self.switch_radius:
            self.switch_radius[key] = CurrentSwitchRadiusPhase(self.owner(ends), self.identity, self.N)
        return self.switch_radius[key].query(chart, left, right)

    def local_integral(self, ends, chart, left, right):
        if chart not in CHARTS:
            return super().local_integral(ends, chart, left, right)
        l, r = backend.exact(left), backend.exact(right)
        if not 0 <= l < r <= 1:
            raise ValueError('Positive actual switch native integral interval required')
        source = self.query(ends, chart, left, right)
        op, original_owner = self.owner(ends), self.switch_owner(ends)
        f, c = op.flow, op.c
        length = f.h if chart != 'post_power' else f.scalar(c.ln(c.mpf(11)/10))-f.h*2
        width = length*(c.mpf((r-l).numerator)/(r-l).denominator)
        density = source['actual_signed_nonlinear_density_C0_Z']
        integrals, masses = {}, {}
        for name, rate in RATES.items():
            mass, _, _ = original.own_weights(original_owner.first, chart, left, right, rate)
            integrals[name] = [density[part][name]*mass for part in ('kernels', 'Z_derivatives')]
            masses[name] = mass
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            actual_chart=chart, exact_left=list(left), exact_right=list(right),
            exact_common_P0_axial5=op.P0, actual_source_function=source,
            actual_physical_log_width=width, original_positive_chart_log_length=length,
            actual_positive_own_rate_masses=masses, actual_signed_local_integral_C0_Z=integrals,
            integral_recipe='int_left^right exp(-rate*L*(right-s))*signed_density(s,Z,N)*L ds',
            physical_Jacobian_applied_once=True, actual_positive_micro_width_not_zeroed=True,
            live_inverse_function_extension_under_Z_independent_integral=True,
            directed_rectangle_function_enclosure_not_selected_integral_value=True,
            incoming_correction_not_supplied_or_reset=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN, False))

    def switch_transport(self, ends):
        op, original_owner = self.owner(ends), self.switch_owner(ends)
        f = op.flow
        initial = original_owner.incoming
        incoming = {name: list(rows) for name, rows in initial.items()}
        windows = []
        for chart in CHARTS:
            local = {name: [f.scalar(0), f.scalar(0)] for name in RATES}
            cells = []
            for left, right in zip(PARTITION, PARTITION[1:]):
                integral = self.local_integral(ends, chart, left, right)
                contribution, weights = {}, {}
                for name, rate in RATES.items():
                    mass, decay, suffix = original.own_weights(original_owner.first, chart, left, right, rate)
                    contribution[name] = [row*suffix for row in integral['actual_signed_local_integral_C0_Z'][name]]
                    for i in range(2):
                        local[name][i] += contribution[name][i]
                    weights[name] = dict(original_full_mass=mass, original_cell_decay=decay,
                        original_suffix_decay=suffix, own_rate=str(rate))
                query = integral['actual_source_function']
                cells.append(dict(exact_left=list(left), exact_right=list(right),
                    current_signed_function_and_integral_sha256=function_digest(integral),
                    actual_signed_primitive_C0_Z_phi=query['actual_signed_primitive_C0_Z_phi'],
                    actual_signed_nonlinear_density_C0_Z=query['actual_signed_nonlinear_density_C0_Z'],
                    actual_signed_local_integral_C0_Z=integral['actual_signed_local_integral_C0_Z'],
                    actual_signed_contribution_at_chart_exit_C0_Z=contribution,
                    actual_physical_log_width=integral['actual_physical_log_width'],
                    original_own_rate_weights=weights,
                    actual_global_phase_full_period_cover=query['actual_source_geometry']['full_period'],
                    full_function_packet_available_from_live_oracle_and_replay_digest=True))
            memory = {name: original.own_weights(original_owner.first, chart, (0,1), (1,1), rate)[1]
                for name, rate in RATES.items()}
            inherited = {name: [row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing = {name: [inherited[name][i]+local[name][i] for i in range(2)] for name in RATES}
            windows.append(dict(actual_chart=chart, actual_signed_cells=cells,
                actual_incoming_correction_C0_Z=incoming, original_incoming_own_rate_memory=memory,
                actual_signed_local_driver_C0_Z=local, actual_retained_incoming_C0_Z=inherited,
                actual_exit_signed_correction_C0_Z=outgoing))
            incoming = outgoing
            print('Current signed switch C1 transport: '+str(ends)+' '+chart, flush=True)
        terminal = self.query(ends, 'post_power', (1,1))
        background = {name: list(rows[:2]) for name, rows in terminal['actual_original_source_packet']
            ['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete = {name: [background[name][i]+incoming[name][i] for i in range(2)] for name in RATES}
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=op.P0, actual_R100_correction_incoming_C0_Z=initial,
            actual_signed_switch_windows=windows, actual_R110_signed_correction_C0_Z=incoming,
            actual_original_R110_background_C0_Z=background, actual_R110_signed_complete_own_history_C0_Z=complete,
            actual_original_R110_signed_source_function=terminal, exact_R110_radius=f.scalar(110),
            exact_total_log_width_identity='h+h+(log(11/10)-2*h)=log(11/10)',
            actual_R100_incoming_not_zeroed_or_replaced_by_background=True,
            actual_all_three_local_coordinates_in_zero_one=True,
            interval_caps_not_substituted_for_defining_functions=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZSwitchSignedFunctions()
        points, transports = [], []
        for ends in CELLS:
            for chart in CHARTS:
                points.append(serialized(owner.query(ends, chart, (1,3))))
            transports.append(serialized(owner.switch_transport(ends)))
        report = dict(**{GATE: True}, source_family=owner.identity, candidate_N=owner.N,
            exact_Z_partition=CELLS, current_switch_charts=CHARTS, exact_native_partition=PARTITION,
            actual_signed_switch_point_functions=points, actual_signed_switch_transports=transports,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Three current whole-Z R100/R110 switch source charts expose actual global-phase '
                'signed inverse C0/Z/phase functions and true local C1 integrals. Original positive '
                'microscopic widths, source P0/flow and genuine same-N R100 incoming are retained. '
                'Earlier bridge/long sources, complete Rc defect functions/controls and global '
                'N/heat/cone/temporal recursion remain unfinished.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report), separators=(',',':'))+'\n')
            .encode(), compresslevel=6, mtime=0))
    print('Current whole-Z switch signed functions and genuine C1 transport reach R110', flush=True)
    return report


if __name__ == '__main__':
    run()
