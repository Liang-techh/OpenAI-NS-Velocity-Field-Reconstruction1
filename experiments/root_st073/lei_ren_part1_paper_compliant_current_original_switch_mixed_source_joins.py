"""Two actual switch interfaces on the accepted full prescribed source.

The first/second switch and second-switch/power interfaces retain actual
source nodes, complete incoming histories and distinct native measures.
These finite source identities do not supply numerical or global n1 bounds.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,downstream_source=previous.outer,previous.downstream_source
NAME=PREFIX+'current_original_switch_mixed_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_switch_mixed_source_joins_check.json'
JOINS=(('first_second','first_switch','second_switch','R1'),
       ('second_power','second_switch','post_power','R2'))
GATES=('current_original_n1_first_second_switch_full_mixed_source_join_installed',
       'current_original_n1_second_switch_post_power_full_mixed_source_join_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)


class CurrentOriginalSwitchMixedSourceJoins:
    def __init__(self,base=None,require_checked=True):
        self.base=base if base is not None else previous.CurrentOriginalPrescribedHighRadialSourceFunctions()
        if not self.base.acceptance_loaded:raise ValueError('Accepted complete prescribed high sources required')
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        self.joins={key:dict(left_chart=left,right_chart=right,original_endpoint=endpoint,
            left=self.base.at_own(left,g.one),right=self.base.at_own(right,g.zero),
            left_radius=outer.Algebra(g,0).at(outer.FunctionJet((self.charts[left]['radius'],)),
                self.charts[left]['coordinate'],g.one),
            right_radius=outer.Algebra(g,0).at(outer.FunctionJet((self.charts[right]['radius'],)),
                self.charts[right]['coordinate'],g.zero)) for key,left,right,endpoint in JOINS}
        self.controls={}
        for name in ('first_switch','second_switch'):
            sigma=downstream_source.sigma(g,self.charts[name]['coordinate']);D=self.factories[name]
            refs=[sigma]
            for n in range(1,5):refs.append(D.derivative(refs[-1]))
            self.controls[name]=refs
        self.definition=dict(actual_first_switch1_second_switch0_and_second_switch1_post_power0=True,
            original_sigma_native_derivatives_zero_through_order_four=True,
            actual_source_field_history_pressure_and_correction_inlets_not_reset=True,
            physical_background_y5_then_leading_y4_Q4_N2_pressure_y3=True,
            own_outputs_not_assigned_to_seam_targets=True,
            source_functions_not_numeric_evaluation_or_temporal_coefficients=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()

    def __getattr__(self,name):return getattr(self.base,name)
    def _snapshot(self):
        return outer.encode(dict(joins=self.joins,controls=self.controls,contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self._snapshot()!=self.snapshot or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):
            raise ValueError('Actual switch source joins changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed switch input '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=5
            or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['switch_mixed_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual switch mixed source joins required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalSwitchMixedSourceJoins(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        switch_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_TWO_SWITCH_MIXED_SOURCE_JOINS_DEFINED',len(owner.prefix),flush=True)
    return owner,raw
