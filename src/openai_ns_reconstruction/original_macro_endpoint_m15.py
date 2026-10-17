"""Actual micro/macro pressure seams and the signed finite-N endpoint ledger.

Hydrate saved packets and prove only the two new interfaces.  No accepted
ancestor constructor, producer, checker, fixture or calculation is invoked.
"""
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT, DIRECTORY, PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15

PARENT = PREFIX + 'compliant_current_original_switch_m15_numeric.json'
PARENT_SHA = 'b4575dd765e43f9f02aac18b5fa11e3825380bd42f6a30aefd5b420bc86ba0e5'


class OriginalMacroEndpointM15(OriginalActiveViscousM15):
    def __init__(self, *, directory=DIRECTORY, precision=150):
        self.directory = directory
        self.ctx = mp.mp.clone(); self.ctx.dps = precision
        self.iv = MPIntervalContext(); self.iv.dps = precision + 30
        self.hashes = {}; self.code_hashes = {}
        raw = (directory / PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest() != PARENT_SHA:
            raise ValueError('Accepted switch report required')
        self.parent = json.loads(raw); self.family = self.parent['source_family']
        if not self.parent['checks']['passed']:
            raise ValueError('Checked switch parent required')
        self.merge(self.parent['input_hashes'], directory, self.hashes)
        self.merge(self.parent['implementation_hashes'], ROOT, self.code_hashes)
        self.hashes[PARENT] = PARENT_SHA
        self.loop = self.archive('prescribed_loop_C5_functions')
        self.modified = self.archive('modified_Q4_pressure_C3_functions')
        self.high_archive = self.archive('prescribed_high_radial_source_functions')
        self.switch_archive = self.archive('switch_mixed_source_joins')
        self.nodes = self.switch_archive['exact_graph_nodes']
        for archive in (self.loop, self.modified, self.high_archive):
            archived = archive['exact_graph_nodes']
            if self.nodes[:len(archived)] != archived:
                raise ValueError('Identical actual macro graph prefixes required')
        self.l = self.loop['prescribed_loop_C5_bindings']
        self.m = self.modified['modified_Q4_pressure_C3_bindings']
        self.high = self.high_archive['prescribed_high_radial_source_bindings']
        self.highbindings = self.high
        self.pstar = self.high['downstream']['Pstar']
        self.chart_names = ('second_micro', 'frozen_macro', 'first_switch')
        self.seam_pairs = {'micro_macro': ('second_micro', 'frozen_macro'),
                           'macro_first': ('frozen_macro', 'first_switch')}
        # These are archived query-node handles.  Node48 is the number2;
        # node2 is the number40 and must never be used as the micro query.
        self.seam_queries = {'micro_macro': (48, 0), 'macro_first': (1, 0)}
        self.chart_upper_queries = {'second_micro': 48, 'frozen_macro': 1, 'first_switch': 1}
        for node, value in ((0, 0), (1, 1), (48, 2)):
            if self.nodes[node] != dict(operation='exact_rational', numerator=value, denominator=1):
                raise ValueError('Actual native endpoint query node changed')
        self.background_charts = {name: self.high['backgrounds'][name] for name in self.chart_names}
        for name in self.chart_names:
            high_chart = self.high['charts'][name]; loop_chart = self.l['charts'][name]
            for high_key, loop_key in (('coordinate', 'coordinate'), ('radius', 'radius'),
                                       ('jacobian', 'physical_measure'), ('phase', 'actual_phase')):
                if high_chart[high_key] != loop_chart[loop_key]:
                    raise ValueError('Actual high-radial geometry differs from original loop chart')
            if high_chart['domain'][-1] != self.chart_upper_queries[name]:
                raise ValueError('Foreign actual upper native endpoint query')
            physical = high_chart['physical']
            modified = self.m['charts'][name]['physical']
            if (physical['U'][0] != modified['U'] or physical['U'][1] != modified['U_y']
                    or physical['m'][0] != modified['m'] or physical['Q'][0] != modified['Q']
                    or physical['Q'][1] != modified['Q_y']):
                raise ValueError('Actual high-radial physical ports differ from pressure archive')
        self.pressure_input_guards = self.parent['pressure_input_identity_guard_controls']
        if not self.pressure_input_guards['passed']:
            raise ValueError('Saved pressure input identity controls required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
        from .original_macro_native import build_macro_native
        self.derivative = build_macro_native(self.nodes, self.loop, self.modified,
                                             self.high_archive, self.parent, root=ROOT)
        self.native = self.derivative; self.native_guards = self.native
        if not self.native['passed']:
            raise ValueError('Actual three-chart native derivative source proof required')
        from .original_macro_source import build_macro_source
        self.source = build_macro_source(self.nodes, self.loop, self.modified,
                                         self.high_archive, self.parent, self.derivative, root=ROOT)
        if not self.source['passed']:
            raise ValueError('Actual micro/macro and macro/first source identities required')
        self.loop_identity = self.source['actual_loop_function_identity']
        if not self.loop_identity['passed']:
            raise ValueError('New actual scoped root/loop identities required')
        for packet in (self.source, self.native, self.loop_identity):
            self.bind_new_pins(packet)
        from .original_macro_pressure import build_macro_pressure_traces
        self.seam_proofs = build_macro_pressure_traces(self)
        for name in ('original_macro_endpoint_m15.py', 'original_macro_source.py', 'original_macro_native.py',
                     'original_macro_numeric_m15.py', 'original_macro_pressure.py'):
            path = 'src/openai_ns_reconstruction/' + name
            self.code_hashes[path] = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()

    def bind_new_pins(self, packet):
        for row in packet.get('source_pins', {}).values():
            if not isinstance(row, dict):
                continue
            if row.get('passed') is False:
                raise ValueError('Changed macro source pin')
            name = row.get('path') or row.get('relative')
            if not name:
                continue
            actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            if actual != row['sha256']:
                raise ValueError('Macro source pin mismatch ' + name)
            if name.startswith('experiments/root_st073/'):
                key = name.rsplit('/', 1)[-1]; target = self.hashes
            else:
                key = name; target = self.code_hashes
            if key in target and target[key] != actual:
                raise ValueError('Macro source pin conflicts with parent')
            target[key] = actual

    def evaluate(self, Z):
        if not hasattr(self, 'numeric'):
            raise ValueError('Build the new signed ledger first')
        row = next((sample for sample in self.numeric['samples'] if sample['Z'] == Z), None)
        if row is None:
            raise ValueError('Only archived axial queries admitted')
        return row

    def build(self):
        from .original_macro_numeric_m15 import build_macro_numeric
        self.numeric = build_macro_numeric(self, self.seam_proofs)
        return dict(source_family=self.family, parent_report=PARENT, parent_report_sha256=PARENT_SHA,
                    source_derivative_theorem=self.source, native_source_guards=self.native,
                    source_seam_proofs=self.seam_proofs,
                    pressure_input_identity_guard_controls=self.pressure_input_guards,
                    pressure_identity_guard_evidence_inherited_unchanged=True,
                    numerical_source_binding=self.numeric['source_binding'], samples=self.numeric['samples'],
                    contract=dict(accepted_ancestor_constructors_not_run=True, source_functions_not_cap_values=True,
                                  raw_213_and_independent_target_retained=True),
                    limits=dict(original_signed_finite_N_inlet_values=False, numeric_Rp_closure=False,
                                global_all_five_common_N_closure=False, exact_heat_energy_stress_flat=False,
                                true_temporal_n_recursion_summation=False, mean_oscillatory_correction=False,
                                independent_corrected_Cartesian_NS=False, full_goal_complete=False))
