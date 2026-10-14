"""Actual main/gap/end radial source traces and own second radial recovery.

All derivatives are physical y derivatives and ordinary Z rows. Higher
pressure radial derivatives, quantitative bounds and global n1 stay open.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_n1_outer_radial_pressure_source_joins as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer;radial=previous.previous;rp=previous.rp
NAME=PREFIX+'current_original_pulse_internal_radial_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_pulse_internal_radial_source_joins_check.json'
GATES=('current_original_pulse_main_gap_end_actual_radial_pressure_C3_source_joins_installed',
    'current_original_pulse_own_second_radial_C3_QV_C2_source_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(key for key in previous.OPEN if key not in GATES)+(
    'current_original_pulse_own_pressure_second_radial_C3_source_functions_installed',)))


class CurrentOriginalPulseInternalRadialSourceJoins:
    def __init__(self,joins=None,require_checked=True):
        self.joins=joins if joins is not None else previous.CurrentOriginalN1OuterRadialPressureSourceJoins()
        if not self.joins.acceptance_loaded:raise ValueError('Checked actual Rp/Rv radial source required')
        self.joins.assert_graph();self.pulse=self.joins.pulse;self.inlet=self.joins.inlet
        self.g=g=self.joins.g;self.family=copy.deepcopy(self.joins.family);self.hashes=dict(self.joins.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        kernel=parent.linear.read(radial.kernels.NAME)['exact_pulse_kernel_functions']
        self.W={key:self.decode(value) for key,value in kernel.items()}
        old_own=parent.linear.read(rp.pulse5.true.own.NAME)['own_main_gap_end_fields']
        self.history_kernels={}
        for chart,variable in (('main','original_own_main_history_v'),('end','original_own_end_history_v')):
            packets=[]
            for jet in old_own[chart]['normalized_linear_histories']:
                closure=self.joins.outer.operand_closure(jet[0])
                ids=sorted(i for i in closure if g.nodes[i]['operation']=='definite_integral'
                    and g.nodes[i]['variable']==variable)
                if len(ids)!=(1 if chart=='main' else 2):raise ValueError('Exact original normalized history kernels required')
                packets.append([self.ref(i) for i in ids])
            self.history_kernels[chart]=packets
        self.source_history_producer=parent.active.current.ast_binding(rp.pulse5.true.own.pulse_fields)
        self.source_selector=parent.active.current.ast_binding(radial.kernels.select_functions)
        self.leading_traces={
            'main_gap':dict(left=self.pulse.endpoints['main']['right'],right=self.pulse.endpoints['gap']['left']),
            'gap_end':dict(left=self.pulse.endpoints['gap']['right'],right=self.pulse.endpoints['end']['left'])}
        a=outer.Algebra(g,3);b=outer.Algebra(g,2);self.own_second_radial={}
        for name,packet in self.joins.current_own.items():
            fields=packet['fields']
            self.own_second_radial[name]=dict(Uz1_yy=a.fixed(0),Utheta1_yy=a.fixed(0),F1_yy=a.fixed(0),
                M1_yy=a.fixed(0),mean1_yy=fields['mean1'],K1_yy=fields['K1'],
                Q1_yy=fields['Q1'],V1_yy=b.fixed(0))
        self.definition=dict(actual_two_internal_leading_radial_source_traces=True,
            actual_finite_mu_history_kernels_and_original_divided_selector_retained=True,
            end_scale_is_exp_logpref_over_mu=True,physical_Jacobian_main_gap_mu_end_one=True,
            inactive_gap_mean_and_all_incoming_histories_retained=True,
            own_second_radial_recovery_one_plus_delta=True,
            leading_one_minus_delta_recovery_and_pressure_density_separate=True,
            pressure_second_radial_and_full_outer_mixed_numeric_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.joins,self.pulse,self.inlet,g);self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def decode(self,value):
        if isinstance(value,list):return [self.decode(i) for i in value]
        return self.ref(value)
    def _snapshot(self):
        return outer.encode(dict(W=self.W,history_kernels=self.history_kernels,
            source_history_producer=self.source_history_producer,source_selector=self.source_selector,
            leading_traces=self.leading_traces,own_second_radial=self.own_second_radial,definition=self.definition))
    def assert_graph(self):
        self.joins.assert_graph()
        if (self._owners!=(self.joins,self.pulse,self.inlet,self.g) or self._hashes!=self.hashes
            or self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.numerical_oracle_available):raise ValueError('Actual internal radial source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed internal radial source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or archive!=sha(NAME) or pins.pop(checker)!=sha(checker) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['pulse_internal_radial_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual internal radial source joins required')
        self.acceptance_loaded=True


def run(joins=None):
    began=time.monotonic();owner=CurrentOriginalPulseInternalRadialSourceJoins(joins,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        pulse_internal_radial_source_join_bindings=owner.snapshot,
        **dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_INTERNAL_PULSE_RADIAL_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
