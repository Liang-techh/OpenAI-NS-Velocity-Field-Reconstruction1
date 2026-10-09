"""Actual whole-Z long/reference signed functions and R110-to-Rm C1 transport.

The current source frontend is reused before its all-u support evaluation.
Selected singleton T/logC use exact modulus; analytic logP is refined.
"""
import ast
from fractions import Fraction
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_switch_signed_functions as previous

backend = previous.backend
source_module = previous.source_module
HERE, PREFIX, sha, bind, ep = previous.HERE, previous.PREFIX, previous.sha, previous.bind, previous.ep
CELLS, RATES, C0, Z, OPEN = previous.CELLS, previous.RATES, previous.C0, previous.Z, previous.OPEN
CHARTS, PARTITION = source_module.CHARTS, source_module.PARTITION
NAME = PREFIX+'current_original_whole_Z_long_Rm_signed_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_long_Rm_signed_functions_check.json'
GATE = 'current_original_whole_Z_long_reference_live_signed_functions_and_C1_transport_through_Rm_installed'
serialized, encode, function_digest = previous.serialized, previous.encode, previous.function_digest


def original_frontend():
    """Keep every original source statement before the all-u support call."""
    fn = source_module.WholeZLongRmFiniteN.query
    node = ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    stops = [i for i, statement in enumerate(node.body) if isinstance(statement, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == 'got' for target in statement.targets)]
    if len(stops) != 2:
        raise ValueError('Original long source/support boundary changed')
    stop = stops[0]
    call = node.body[stop].value
    if ast.unparse(call.func) != 'primitives.all_u_primitive_bounds':
        raise ValueError('Actual original all-u support boundary required')
    original = node.body[:stop]
    proof = dict(original_source_frontend_before_all_u_AST_unchanged=True,
        exact_original_frontend_AST_sha256=hashlib.sha256(ast.dump(ast.Module(body=original,
            type_ignores=[])).encode()).hexdigest(),
        original_source_query_module=Path(source_module.__file__).name,
        original_source_query_sha256=sha(Path(source_module.__file__).name),
        support_caps_and_ancestor_contribution_not_evaluated=True)
    tail = ast.parse('''result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
        chart=chart,exact_left=list(left),exact_right=list(right),original_window_length=length,
        exact_common_P0_axial5=owner.reference.P0,actual_background_source=source,
        original_raw_source=raw,original_generic_source=recovered,original_full_source_quotients=proof,
        original_roots=roots['roots'],original_q_C0_Z=qr,actual_phase_Z_exact_zero=True,
        actual_original_source_frontend_binding=FRONTEND_BINDING)
self.cache[key]=result
return result
''').body
    node.name = 'current_original_long_source_frontend'
    node.body = original+tail
    module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
    scope = {**fn.__globals__, 'FRONTEND_BINDING': proof}
    exec(compile(module, '<unchanged current whole-Z long source frontend>', 'exec'), scope)
    return scope[node.name], proof


SOURCE_FRONTEND, FRONTEND_BINDING = original_frontend()


class CurrentLongRadiusPhase:
    def __init__(self, op, long_source, long_owner, identity, N):
        self.op, self.identity, self.N = op, identity, N
        self.base = backend.phase.RmRadiusPhase(op, identity['actual_five_defect_family_sha256'], identity)
        self.T, self.logC, self.gap = long_source.T, long_owner.reference.logC, long_owner.reference.gap
        if long_owner.flow is not op.flow or long_owner.reference.flow is not op.flow:
            raise ValueError('Same actual long/reference flow required')
        if self.T._mpi_ != (400*long_source.A)._mpi_ or self.T._mpi_ != long_owner.long.T._mpi_ \
                or self.T._mpi_ != long_owner.reference.T._mpi_:
            raise ValueError('Original defining T=400*A and actual live source tuples required')
        if self.logC._mpi_ != self.base.logC._mpi_:
            raise ValueError('Exact same source-selected logCstar required')
        if any(ep(row)[0] != ep(row)[1] for row in (self.T, self.logC)):
            raise ValueError('Actual source-selected singleton T/logC required for exact modulus')
        if self.gap._mpi_ != (10*(self.logC+op.flow.logs[1]/2)-self.T-8)._mpi_:
            raise ValueError('Exact original reference gap source required')
        self.cache = {}

    def length(self, chart):
        if chart not in CHARTS:
            raise ValueError('Current original long/reference chart required')
        return self.T if chart == 'long_reshape' else self.gap if chart == 'reference' else self.op.c.mpf(1)

    def point(self, chart, coordinate):
        q = backend.exact(coordinate)
        if chart not in CHARTS or not 0 <= q <= 1:
            raise ValueError('Actual long/reference native point in [0,1] required')
        c = MPIntervalContext()
        c.dps = max(260, 80+(self.N.bit_length()*30103+99999)//100000+60)
        with mp.workdps(c.dps+40):
            cv = lambda value: c.mpf(value.numerator)/value.denominator
            constant, p, tc, cc = Fraction(1000), Fraction(4), q, Fraction(0)
            if chart == 'reference':
                constant, p, tc, cc = 1000-8*q, 4+10*q, 1-q, 10*q
            elif chart in ('restoration', 'postrestore'):
                constant, p, tc, cc = 1000+(-8 if chart == 'restoration' else -7)+q, Fraction(14), Fraction(0), Fraction(10)
            projection = dict(full_period=False, boxes=[c.mpf(0)])
            components = []
            for name, value in (('constant', c.ln(c.mpf(110)/4)+cv(Fraction(constant))),
                    ('logP', cv(Fraction(p))*(c.exp(40)+11))):
                wrapped = backend.signed.first.spatial.ordinary_mod_one(c, self.N*value)
                if wrapped['full_period']:
                    raise ArithmeticError('Actual long analytic phase needs higher precision')
                pieces = []
                for part in wrapped['boxes']:
                    added = backend.radius.periodic_add(c, projection['boxes'], part)
                    if added['full_period']:
                        raise ArithmeticError('Actual long analytic phase sum needs refinement')
                    pieces.extend(added['boxes'])
                projection = backend.radius.periodic_add(c, pieces, c.mpf(0))
                components.append(dict(component=name, original_component=value, directed_projection=wrapped))
            exact = []
            for name, selected, coefficient in (('T', self.T, tc), ('logC', self.logC, cc)):
                if not coefficient:
                    continue
                value, proof = backend.signed.first.spatial.binary_mod_one(c, c.mpf(selected), Fraction(coefficient)*self.N)
                projection = backend.radius.periodic_add(c, projection['boxes'], value)
                exact.append(dict(component=name, **proof))
            epsilon = c.mpf(10)**-100
            logupper = c.ln(self.N)+c.mpf(self.op.flow.logs[0])+c.ln(c.mpf(self.base.sc)/2)
            if ep(logupper)[1] >= ep(c.ln(epsilon))[0]:
                raise ArithmeticError('Actual positive inlet width exceeds long phase budget')
            error = c.mpf((-ep(epsilon)[1], 0))
            projection = backend.radius.periodic_add(c, projection['boxes'], error)
            if projection['full_period']:
                raise ArithmeticError('Actual long point phase requires refinement')
            result = dict(source_identity=self.identity, candidate_N=self.N, actual_chart=chart,
                exact_left=list(coordinate), exact_right=list(coordinate),
                actual_phase_boxes=[self.op.c.mpf(row) for row in projection['boxes']], full_period=False,
                actual_phase_Z_exact_zero=True, actual_same_source_Rm_radius=self.op.Rm_factor,
                actual_original_window_length=self.length(chart), exact_live_source_T=self.T,
                exact_source_logCstar=self.logC, exact_live_source_reference_gap=self.gap,
                original_parameter_binding=self.base.parameter_binding,
                exact_defining_T_400_A_and_live_long_reference_tuples_bound=True,
                original_regular_component_proofs=components, exact_selected_dyadic_modulus_proofs=exact,
                actual_N_inlet_width_phase_log_upper=logupper, signed_actual_positive_inlet_width_error=self.op.c.mpf(error),
                actual_positive_width_retained_not_zeroed=True, selected_source_sc=self.base.sc,
                adaptive_analytic_interval_digits=c.dps, global_phase_not_restarted=True,
                phase_not_supplied_as_free_angle=True,
                exact_original_global_phase='frac(N*(logR-logRa-hb*s_c/2))',
                exact_original_long_radius_maps=dict(long_reshape='log110+T*t',
                    reference='log110+T+gap*t', restoration='log110+10*(C+P)-8+t',
                    postrestore='log110+10*(C+P)-7+t'))
        return result

    def query(self, chart, left, right=None):
        right = left if right is None else right
        l, r = backend.exact(left), backend.exact(right)
        if chart not in CHARTS or not 0 <= l <= r <= 1:
            raise ValueError('Ordered actual long native coordinates in [0,1] required')
        key = (chart, l, r)
        if key in self.cache:
            return self.cache[key]
        result = self.point(chart, left)
        result['exact_right'] = list(right)
        if r > l:
            c = self.op.c
            width = self.length(chart)*(c.mpf((r-l).numerator)/(r-l).denominator)
            variation = self.N*width
            if ep(variation)[0] >= 1:
                projection = dict(full_period=True, boxes=[c.mpf((0,1))])
            else:
                projection = backend.radius.periodic_add(c, result['actual_phase_boxes'], c.mpf((0, ep(variation)[1])))
            result.update(actual_phase_boxes=projection['boxes'], full_period=projection['full_period'],
                original_positive_physical_log_width=width, actual_monotone_phase_variation=variation,
                full_period_proved_from_actual_N_physical_width=ep(variation)[0] >= 1)
        self.cache[key] = result
        return result


class WholeZLongRmSignedFunctions(previous.WholeZSwitchSignedFunctions):
    def __init__(self, dps=500):
        super().__init__(dps)
        checked = json.loads((HERE/previous.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(previous.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked actual current signed switch backend required')
        for name, digest in checked['input_hashes'].items():
            bind(self.hashes, name, digest)
        bind(self.hashes, previous.RECEIPT, sha(previous.RECEIPT))
        self.long_source = self.patch_source.source.source.upstream
        if type(self.long_source) is not source_module.WholeZLongRmFiniteN \
                or self.long_source.source is not self.switch_source:
            raise ValueError('Same live whole-Z long/Rm/switch chain required')
        self.frontend_owner = SimpleNamespace(**{**vars(self.long_source), 'cache': {}, 'owner': self.long_source.owner})
        self.long_radius = {}
        for module in (previous, source_module, source_module.long, source_module.reference):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def long_owner(self, ends):
        op, source = self.owner(ends), self.long_source.owner(ends)
        if source.flow is not op.flow or source.c is not op.c:
            raise ValueError('Identical actual long flow/context required')
        p0 = source.reference.P0
        if len(p0) != len(op.P0) or any(a.coefficient._mpi_ != b.coefficient._mpi_
                or a.scale.powers != b.scale.powers or a.scale.offset._mpi_ != b.scale.offset._mpi_
                or a.zero != b.zero or a.ctx is not b.ctx or a.scale.bases is not b.scale.bases
                or a.ledger is not b.ledger for a, b in zip(p0, op.P0)):
            raise ValueError('Exact source-equivalent long P0 rows required')
        return source

    def query(self, ends, chart, left, right=None):
        if chart in CHARTS:
            return previous._identity_query(self, ends, chart, left, right)
        return super().query(ends, chart, left, right)

    def source_query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().source_query(ends, chart, left, right)
        l, r = backend.exact(left), backend.exact(left if right is None else right)
        if not 0 <= l <= r <= 1:
            raise ValueError('Current original long source native domain required')
        self.long_owner(ends)
        packet = dict(SOURCE_FRONTEND(self.frontend_owner, ends, chart, left, right))
        packet['original_long_independent_P0_axial5'] = packet['exact_common_P0_axial5']
        packet['exact_common_P0_axial5'] = self.owner(ends).P0
        packet['identical_source_P0_tuples_rebound_to_live_descendant_object'] = True
        return packet

    def radius_query(self, ends, chart, left, right=None):
        if chart not in CHARTS:
            return super().radius_query(ends, chart, left, right)
        key = tuple(ends)
        if key not in self.long_radius:
            self.long_radius[key] = CurrentLongRadiusPhase(self.owner(ends), self.long_source,
                self.long_owner(ends), self.identity, self.N)
        return self.long_radius[key].query(chart, left, right)

    def local_integral(self, ends, chart, left, right):
        if chart not in CHARTS:
            return super().local_integral(ends, chart, left, right)
        l, r = backend.exact(left), backend.exact(right)
        if not 0 <= l < r <= 1:
            raise ValueError('Positive actual long/reference integral interval required')
        source = self.query(ends, chart, left, right)
        op, original = self.owner(ends), self.long_owner(ends)
        f, c = op.flow, op.c
        length = self.long_source.T if chart == 'long_reshape' else original.reference.gap if chart == 'reference' else c.mpf(1)
        width = length*(c.mpf((r-l).numerator)/(r-l).denominator)
        density, integrals, masses = source['actual_signed_nonlinear_density_C0_Z'], {}, {}
        for name, rate in RATES.items():
            mass, _, _ = source_module.long.kernel_weight(f, width, c.mpf(0), rate)
            integrals[name] = [density[part][name]*mass for part in ('kernels', 'Z_derivatives')]
            masses[name] = mass
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            actual_chart=chart, exact_left=list(left), exact_right=list(right), exact_common_P0_axial5=op.P0,
            actual_source_function=source, actual_original_window_length=length, actual_physical_log_width=width,
            actual_positive_own_rate_masses=masses, actual_signed_local_integral_C0_Z=integrals,
            integral_recipe='int_left^right exp(-rate*L*(right-s))*signed_density(s,Z,N)*L ds',
            physical_Jacobian_applied_once=True, original_T_gap_one_one_windows_not_shortened=True,
            live_inverse_function_extension_under_Z_independent_integral=True,
            directed_rectangle_function_enclosure_not_selected_integral_value=True,
            incoming_correction_not_supplied_or_reset=True, **dict.fromkeys(OPEN, False),
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False)

    def long_transport(self, ends):
        op, original = self.owner(ends), self.long_owner(ends)
        f, c = op.flow, op.c
        inlet = self.switch_transport(ends)
        initial = inlet['actual_R110_signed_correction_C0_Z']
        incoming = {name: list(rows) for name, rows in initial.items()}
        windows = []
        for chart in CHARTS:
            length = self.long_source.T if chart == 'long_reshape' else original.reference.gap if chart == 'reference' else c.mpf(1)
            local, cells = {name: [f.scalar(0), f.scalar(0)] for name in RATES}, []
            for left, right in zip(PARTITION, PARTITION[1:]):
                integral = self.local_integral(ends, chart, left, right)
                r = backend.exact(right)
                suffix = length*(c.mpf((1-r).numerator)/(1-r).denominator)
                contributions, weights = {}, {}
                for name, rate in RATES.items():
                    mass, decay, tail = source_module.long.kernel_weight(f, integral['actual_physical_log_width'], suffix, rate)
                    contributions[name] = [row*tail for row in integral['actual_signed_local_integral_C0_Z'][name]]
                    for i in range(2):
                        local[name][i] += contributions[name][i]
                    weights[name] = dict(original_full_mass=mass, original_cell_decay=decay, original_suffix_decay=tail)
                query = integral['actual_source_function']
                cells.append(dict(exact_left=list(left), exact_right=list(right),
                    current_signed_function_and_integral_sha256=function_digest(integral),
                    actual_signed_primitive_C0_Z_phi=query['actual_signed_primitive_C0_Z_phi'],
                    actual_signed_nonlinear_density_C0_Z=query['actual_signed_nonlinear_density_C0_Z'],
                    actual_signed_local_integral_C0_Z=integral['actual_signed_local_integral_C0_Z'],
                    actual_signed_contribution_at_chart_exit_C0_Z=contributions,
                    original_own_rate_weights=weights, actual_physical_log_width=integral['actual_physical_log_width'],
                    actual_global_phase_full_period_cover=query['actual_source_geometry']['full_period']))
            memory = {name: source_module.long.kernel_weight(f, length, c.mpf(0), rate)[1] for name, rate in RATES.items()}
            inherited = {name: [row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing = {name: [inherited[name][i]+local[name][i] for i in range(2)] for name in RATES}
            windows.append(dict(actual_chart=chart, actual_original_window_length=length, actual_signed_cells=cells,
                actual_incoming_correction_C0_Z=incoming, original_incoming_own_rate_memory=memory,
                actual_signed_local_driver_C0_Z=local, actual_retained_incoming_C0_Z=inherited,
                actual_exit_signed_correction_C0_Z=outgoing))
            incoming = outgoing
            print('Current signed long/reference C1 transport: '+str(ends)+' '+chart, flush=True)
        terminal = self.query(ends, 'postrestore', (1,1))
        background = {name: list(rows[:2]) for name, rows in terminal['actual_original_source_packet']
            ['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete = {name: [background[name][i]+incoming[name][i] for i in range(2)] for name in RATES}
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=op.P0, actual_signed_R100_R110_inlet_transport_sha256=function_digest(inlet),
            actual_R110_current_signed_correction_incoming_C0_Z=initial, actual_signed_long_windows=windows,
            actual_Rm_signed_correction_C0_Z=incoming, actual_original_Rm_background_C0_Z=background,
            actual_Rm_signed_complete_own_history_C0_Z=complete, actual_original_Rm_signed_source_function=terminal,
            exact_Rm_radius=op.Rm_factor, exact_total_log_length_identity='T+gap+1+1=10*(logC+logP)-6',
            current_signed_switch_output_not_ancestor_support_incoming=True,
            original_R110_background_remains_separate_from_signed_correction=True,
            interval_caps_not_substituted_for_defining_functions=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZLongRmSignedFunctions()
        points, transports = [], []
        for ends in CELLS:
            for chart in CHARTS:
                points.append(serialized(owner.query(ends, chart, (1,3))))
            transports.append(serialized(owner.long_transport(ends)))
        report = dict(**{GATE: True}, source_family=owner.identity, candidate_N=owner.N,
            exact_Z_partition=CELLS, current_long_charts=CHARTS, exact_native_partition=PARTITION,
            original_live_source_frontend_binding=FRONTEND_BINDING,
            actual_signed_long_point_functions=points, actual_signed_long_transports=transports,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Four actual whole-Z long/reference/restoration/postrestore signed inverse C0/Z/phase '
                'source functions and genuine C1 transport from the current signed R110 inlet through Rm. '
                'Original T/gap/1/1 windows, exact source-selected dyadic parameters, P0 and background remain. '
                'Three bridge charts, full inner-to-Rc C1 defects/controls/global N/heat/cone and temporal recursion remain unfinished.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report), separators=(',',':'))+'\n').encode(),
            compresslevel=6, mtime=0))
    print('Current whole-Z long/reference signed functions and genuine C1 transport reach Rm', flush=True)
    return report


if __name__ == '__main__':
    run()
