"""First-switch/second-switch actual source and signed m15 continuation.

Only immutable packets are hydrated. No accepted ancestor constructor,
producer, checker, fixture or original calculation is rerun.
"""
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15
from .original_upstream_m15 import digest

PARENT=PREFIX+'compliant_current_original_upstream_m15_numeric.json'
PARENT_SHA='6028d49ac9ec1b94ae967b15d698366061b4a15b9e0ddef9b45be2d3971e355e'


class OriginalSwitchM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted upstream state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.leading=self.archive('prescribed_R100_Rm_C6_functions')
        self.switch=self.archive('switch_mixed_source_joins')
        self.nodes=self.switch['exact_graph_nodes']
        for archive in (self.loop,self.modified,self.leading):
            if self.nodes[:len(archive['exact_graph_nodes'])]!=archive['exact_graph_nodes']:
                raise ValueError('Identical original switch graph prefixes required')
        self.l=self.loop['prescribed_loop_C5_bindings'];self.m=self.modified['modified_Q4_pressure_C3_bindings']
        self.b=self.leading['prescribed_R100_Rm_C6_bindings'];self.pstar=self.b['Pstar']
        self.chart_names=('first_switch','second_switch');self.seam_pairs={'first_second':self.chart_names}
        self.background_charts={name:self.b['charts'][name] for name in self.chart_names}
        self.native=self.parent['native_source_guards'];self.native_guards=self.native
        if not self.native['passed']:raise ValueError('Accepted native source wiring required')
        for name in self.chart_names:
            guard=self.native['chart_guards'][name];chart=self.l['charts'][name]
            if not guard['passed'] or guard['coordinate']!=chart['coordinate'] or guard['radius']!=chart['radius'] or guard['native_J']!=chart['physical_measure']:
                raise ValueError('Native certificate bound to a foreign actual chart')
            pressure=self.native['pressure_ports'][name]
            if not pressure['passed'] or pressure['physical_rows']!=self.m['charts'][name]['physical']:
                raise ValueError('Saved actual physical pressure ports required')
        # These unchanged controls are saved evidence, not reexecuted tests.
        self.pressure_input_guards=self.parent['pressure_input_identity_guard_controls']
        if not self.pressure_input_guards['passed']:raise ValueError('Accepted input-identity guard contract required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),20000))
        from .original_switch_source import build_switch_source
        self.source=build_switch_source(self.nodes,self.loop,self.modified,self.leading,self.parent,root=ROOT)
        if not self.source['passed']:raise ValueError('Actual first/second switch source identity required')
        self.derivative=self.source['native_derivative_congruence']
        self.loop_identity=self.source['actual_loop_function_identity']
        if not self.derivative['passed'] or not self.loop_identity['passed']:
            raise ValueError('Actual native and all-root/loop certificates required')
        self.bind_new_pins(self.source);self.bind_new_pins(self.derivative);self.bind_new_pins(self.loop_identity)
        from .original_scoped_pressure_traces import build_scoped_pressure_traces
        self.seam_proofs=build_scoped_pressure_traces(self)
        for name in ('original_switch_m15.py','original_switch_source.py','original_switch_numeric_m15.py','original_scoped_pressure_traces.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

    def bind_new_pins(self,packet):
        for row in packet.get('source_pins',{}).values():
            if not isinstance(row,dict):continue
            if row.get('passed') is False:raise ValueError('Changed switch source pin')
            name=row.get('path') or row.get('relative')
            if not name:continue
            actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            if actual!=row['sha256']:raise ValueError('Switch source pin mismatch '+name)
            if name.startswith('experiments/root_st073/'):
                key=name.rsplit('/',1)[-1];target=self.hashes
            else:key=name;target=self.code_hashes
            if key in target and target[key]!=actual:raise ValueError('Switch pin conflicts with parent')
            target[key]=actual

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the correlated view first')
        row=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if row is None:raise ValueError('Only archived axial queries admitted')
        return row

    def build(self):
        from .original_switch_numeric_m15 import build_switch_numeric
        self.numeric=build_switch_numeric(self,self.seam_proofs)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            source_derivative_theorem=self.source,native_source_guards=self.native,
            source_seam_proofs=self.seam_proofs,pressure_input_identity_guard_controls=self.pressure_input_guards,
            pressure_identity_guard_evidence_inherited_unchanged=True,
            numerical_source_binding=self.numeric['source_binding'],samples=self.numeric['samples'],
            contract=dict(accepted_ancestor_constructors_not_run=True,source_functions_not_cap_values=True,
                raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_inlet_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
