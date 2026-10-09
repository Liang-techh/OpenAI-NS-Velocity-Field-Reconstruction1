"""Current signed correction-only graph from the genuine inlet through Rc.

Source-owned functions remain independent of correction histories. Only
the two stale incoming callbacks are replaced; original transport bodies,
physical measures, chart partitions and pressure memory are retained.
"""
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_signed_functions as previous

backend = previous.backend
switch_backend = previous.previous.previous
patch_backend = switch_backend.previous
o3 = backend.current
axial, slope, rh = o3.upstream, o3.upstream.upstream, o3.upstream.upstream.upstream
HERE, PREFIX, sha, bind, ep = previous.HERE, previous.PREFIX, previous.sha, previous.bind, previous.ep
CELLS, RATES, C0, Z, OPEN = previous.CELLS, previous.RATES, previous.C0, previous.Z, previous.OPEN
serialized, encode, function_digest = previous.serialized, previous.encode, previous.function_digest
OUTER_CHARTS = backend.CHARTS
OUTER_PARTITIONS = dict(Rh_reference=rh.PARTITION, O2_slope=slope.PARTITION,
    O2_axial=axial.PARTITIONS['axial'], O2_buffer=axial.PARTITIONS['buffer'],
    O3_transition=o3.PARTITIONS['transition'], O3_power=o3.PARTITIONS['power'])
NAME = PREFIX+'current_original_whole_Z_full_signed_transport.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_full_signed_transport_check.json'
GATE = 'current_original_whole_Z_actual_signed_zero_inlet_to_Rc_integral_graph_installed'


def incoming_callback(fn, role):
    """Replace exactly one original incoming-source expression."""
    node = ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    original = ast.dump(node)
    expected = 'initial = original_owner.incoming' if role == 'switch' else 'inlet = self.patch_source.inlet(ends)'
    replacement = "initial = self.bridge_transport(ends)['actual_R100_signed_correction_C0_Z']" \
        if role == 'switch' else 'inlet = self.current_patch_inlet(ends)'
    wanted = ast.dump(ast.parse(expected).body[0])
    hits = [i for i, statement in enumerate(node.body) if ast.dump(statement) == wanted]
    if len(hits) != 1:
        raise ValueError('Exactly one original correction-only incoming source required: '+role)
    index = hits[0]
    old = node.body[index]
    node.body[index] = ast.parse(replacement).body[0]
    changed = ast.dump(node)
    restored = node.body[index]
    node.body[index] = old
    if ast.dump(node) != original:
        raise ValueError('Original transport body changed outside incoming callback')
    node.body[index] = restored
    node.name = 'current_'+role+'_transport'
    module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
    scope = dict(fn.__globals__)
    exec(compile(module, '<same original transport with current '+role+' incoming>', 'exec'), scope)
    proof = dict(original_transport_module=Path(inspect.getfile(fn)).name,
        original_transport_module_sha256=sha(Path(inspect.getfile(fn)).name),
        original_transport_AST_sha256=hashlib.sha256(original.encode()).hexdigest(),
        callback_transport_AST_sha256=hashlib.sha256(changed.encode()).hexdigest(),
        exact_replaced_statement=expected, exact_current_statement=replacement,
        exact_replaced_statement_index=index,
        original_body_restored_exactly_by_reversing_single_incoming_callback=True,
        original_local_integrals_suffixes_memories_background_and_output_body_unchanged=True)
    return scope[node.name], proof


SWITCH_TRANSPORT, SWITCH_BINDING = incoming_callback(switch_backend.WholeZSwitchSignedFunctions.switch_transport, 'switch')
PATCH_TRANSPORT, PATCH_BINDING = incoming_callback(patch_backend.WholeZRmRhSignedFunctions.transport, 'patch')


class WholeZFullSignedTransport(previous.WholeZBridgeSignedFunctions):
    def __init__(self, dps=500):
        super().__init__(dps)
        checked = json.loads((HERE/previous.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(previous.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked actual signed whole-Z bridge source and inlet required')
        for name, digest in checked['input_hashes'].items():
            bind(self.hashes, name, digest)
        bind(self.hashes, previous.RECEIPT, sha(previous.RECEIPT))
        self.transports = {}
        for module in (previous, switch_backend, patch_backend, o3, axial, slope, rh, axial.original):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def query(self, ends, chart, left, right=None):
        if chart in previous.CHARTS:
            return super().query(ends, chart, left, right)
        # Same inverse/source body on every chart; only the already checked
        # bounded A/N scalar callback permits directed broad negative tails.
        result = previous.BRIDGE_QUERY(self, ends, chart, left, right)
        if chart == 'actual_patch':
            result['exact_left'] = patch_backend.coordinate_record(left)
            result['exact_right'] = patch_backend.coordinate_record(left if right is None else right)
        result['actual_source_bounded_exponent_callback_binding'] = previous.DENSITY_BINDING
        return result

    def bridge_transport(self, ends):
        key = ('bridge', tuple(ends))
        if key not in self.transports:
            self.transports[key] = super().bridge_transport(ends)
        return self.transports[key]

    def switch_transport(self, ends):
        key = ('switch', tuple(ends))
        if key not in self.transports:
            inlet = self.bridge_transport(ends)
            result = SWITCH_TRANSPORT(self, ends)
            result.update(actual_current_bridge_inlet_transport_sha256=function_digest(inlet),
                current_signed_bridge_output_not_ancestor_support_incoming=True,
                exact_current_incoming_callback_binding=SWITCH_BINDING)
            self.transports[key] = result
        return self.transports[key]

    def long_transport(self, ends):
        key = ('long', tuple(ends))
        if key not in self.transports:
            self.transports[key] = super().long_transport(ends)
        return self.transports[key]

    def current_patch_inlet(self, ends):
        op = self.owner(ends)
        inlet = self.long_transport(ends)
        return dict(source_identity=self.identity, candidate_N=self.N, exact_Z_cell=list(ends),
            exact_common_P0_axial5=op.P0, exact_Rm_radius=op.Rm_factor,
            actual_Rm_correction_C0_Z=inlet['actual_Rm_signed_correction_C0_Z'],
            actual_current_long_inlet_transport_sha256=function_digest(inlet),
            current_signed_long_output_not_ancestor_support_incoming=True,
            correction_only_not_complete_background_history=True)

    def transport(self, ends):
        key = ('patch', tuple(ends))
        if key not in self.transports:
            result = PATCH_TRANSPORT(self, ends)
            result['exact_current_incoming_callback_binding'] = PATCH_BINDING
            self.transports[key] = result
        return self.transports[key]

    def outer_weights(self, ends, chart, left, right, rate):
        f = self.owner(ends).flow
        if chart == 'Rh_reference':
            return rh.reference_weights(f, left, right, rate)
        if chart == 'O2_slope':
            return slope.slope_weights(f, left, right, rate)
        if chart in ('O2_axial', 'O2_buffer'):
            width, suffix, mass, tail = axial.original.physical_weights(f, chart.removeprefix('O2_'), left, right, rate)
        elif chart in ('O3_transition', 'O3_power'):
            width, suffix, mass, tail = o3.physical_weights(f, chart.removeprefix('O3_'), left, right, rate)
        else:
            raise ValueError('Actual original outer transport chart required')
        decay = f.scalar(1) if not rate else f.factor((0,0,0,0,0), -width*f.c.mpf(rate.numerator)/rate.denominator)
        return width, suffix, f.scalar(mass), decay, tail

    def outer_transport(self, ends):
        key = ('outer', tuple(ends))
        if key in self.transports:
            return self.transports[key]
        op = self.owner(ends)
        f = op.flow
        inlet = self.transport(ends)
        incoming = {name: list(rows) for name, rows in inlet['actual_Rh_signed_correction_C0_Z'].items()}
        windows = []
        for chart in OUTER_CHARTS:
            partition = OUTER_PARTITIONS[chart]
            local, cells = {name: [f.scalar(0), f.scalar(0)] for name in RATES}, []
            for left, right in zip(partition, partition[1:]):
                integral = self.local_integral(ends, chart, left, right)
                contributions, weights = {}, {}
                for name, rate in RATES.items():
                    width, suffix, mass, decay, tail = self.outer_weights(ends, chart, left, right, rate)
                    contributions[name] = [row*tail for row in integral['actual_signed_local_integral_C0_Z'][name]]
                    for i in range(2):
                        local[name][i] += contributions[name][i]
                    weights[name] = dict(physical_log_width=width, suffix_log_width=suffix,
                        original_full_mass=mass, original_cell_decay=decay, original_suffix_decay=tail)
                query = integral['actual_source_function']
                cells.append(dict(exact_left=list(left), exact_right=list(right),
                    current_signed_function_and_integral_sha256=function_digest(integral),
                    actual_signed_primitive_C0_Z_phi=query['actual_signed_primitive_C0_Z_phi'],
                    actual_signed_nonlinear_density_C0_Z=query['actual_signed_nonlinear_density_C0_Z'],
                    actual_signed_local_integral_C0_Z=integral['actual_signed_local_integral_C0_Z'],
                    actual_signed_contribution_at_chart_exit_C0_Z=contributions,
                    original_own_rate_weights=weights,
                    actual_physical_log_width=integral['actual_physical_log_width'],
                    actual_global_phase_full_period_cover=query['actual_source_geometry']['full_period']))
            memory = {name: self.outer_weights(ends, chart, partition[0], partition[-1], rate)[3]
                for name, rate in RATES.items()}
            inherited = {name: [row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing = {name: [inherited[name][i]+local[name][i] for i in range(2)] for name in RATES}
            terminal = self.query(ends, chart, partition[-1])
            background = {name: list(rows[:2]) for name, rows in terminal['actual_original_source_packet']
                ['original_generic_source']['common_own_five_histories_axial5'].items()}
            complete = {name: [background[name][i]+outgoing[name][i] for i in range(2)] for name in RATES}
            windows.append(dict(actual_chart=chart, actual_signed_cells=cells,
                actual_incoming_correction_C0_Z=incoming, original_incoming_own_rate_memory=memory,
                actual_retained_incoming_C0_Z=inherited, actual_signed_local_driver_C0_Z=local,
                actual_exit_signed_correction_C0_Z=outgoing, actual_original_exit_background_C0_Z=background,
                actual_exit_signed_complete_own_history_C0_Z=complete,
                current_terminal_source_function_sha256=function_digest(terminal)))
            incoming = outgoing
            print('Current full signed outer C1 transport: '+str(ends)+' '+chart, flush=True)
        terminal = self.query(ends, 'O3_power', OUTER_PARTITIONS['O3_power'][-1])
        result = dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=op.P0, actual_current_patch_inlet_transport_sha256=function_digest(inlet),
            actual_Rh_current_signed_correction_incoming_C0_Z=inlet['actual_Rh_signed_correction_C0_Z'],
            actual_signed_outer_windows=windows, actual_Rc_signed_correction_C0_Z=incoming,
            actual_original_Rc_background_C0_Z=windows[-1]['actual_original_exit_background_C0_Z'],
            actual_Rc_signed_complete_own_history_C0_Z=windows[-1]['actual_exit_signed_complete_own_history_C0_Z'],
            actual_original_Rc_signed_source_function=terminal,
            actual_positive_Rc_normalized_amplitude_C0_Z=[terminal['actual_original_source_packet']['original_roots']['E'][part]
                for part in (C0, Z)],
            exact_Rc_radius=op.Rm_factor*f.factor((0,.5,0,0,0), 9),
            exact_Rc_over_Rm_log_length='exp40+20=logPstar+9',
            current_signed_patch_output_not_ancestor_support_incoming=True,
            original_outer_partitions_and_physical_Jacobians_preserved=True,
            quiet_regions_retain_predecessor_correction_and_pressure_memory=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN, False))
        self.transports[key] = result
        return result

    def full_transport(self, ends):
        outer = self.outer_transport(ends)
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=self.owner(ends).P0,
            actual_source_owned_zero_inlet_to_Rc_signed_stages=dict(
                bridge=self.bridge_transport(ends), switch=self.switch_transport(ends),
                long=self.long_transport(ends), patch=self.transport(ends), outer=outer),
            actual_Rc_signed_correction_C0_Z=outer['actual_Rc_signed_correction_C0_Z'],
            actual_original_Rc_background_C0_Z=outer['actual_original_Rc_background_C0_Z'],
            actual_Rc_signed_complete_own_history_C0_Z=outer['actual_Rc_signed_complete_own_history_C0_Z'],
            exact_Rc_radius=outer['exact_Rc_radius'],
            actual_current_predecessor_output_used_at_every_radial_boundary=True,
            actual_same_source_signed_C0_Z_integral_graph_through_Rc=True,
            terminal_defects_controls_and_closure_not_inferred_from_history_transport=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZFullSignedTransport()
        transports = [serialized(owner.full_transport(ends)) for ends in CELLS]
        report = dict(**{GATE: True}, source_family=owner.identity, candidate_N=owner.N,
            exact_Z_partition=CELLS, exact_outer_native_partitions=OUTER_PARTITIONS,
            original_single_incoming_callback_bindings=dict(switch=SWITCH_BINDING, patch=PATCH_BINDING),
            actual_signed_full_transports=transports, all_seventeen_actual_signed_source_charts_connected=True,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Actual source-defined zero correction inlet -> bridge R100 -> switch R110 -> long Rm '
                '-> actual patch Rh -> six original outer charts -> Rc, with genuine current signed '
                'C0/Z predecessor correction and original background separate at every boundary. '
                'Rc defect normalization/control/terminal zero, global N/heat/cone and temporal recursion remain unfinished.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report), separators=(',',':'))+'\n').encode(),
            compresslevel=6, mtime=0))
    print('Current whole-Z signed zero-inlet C1 correction transport reaches Rc', flush=True)
    return report


if __name__ == '__main__':
    run()
