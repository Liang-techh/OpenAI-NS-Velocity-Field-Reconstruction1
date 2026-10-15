"""Four actual R110/reshape/reference/restoration source interfaces."""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,downstream_source=previous.outer,previous.downstream_source
NAME=PREFIX+'current_original_annular_mixed_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_annular_mixed_source_joins_check.json'
JOINS=(('power_reshape','post_power','long_reshape','R110'),
       ('reshape_reference','long_reshape','reference','Rsh'),
       ('reference_restore','reference','restoration','Rz'),
       ('restore_postrestore','restoration','postrestore','Rrestore'))
CHARTS=('post_power','long_reshape','reference','restoration','postrestore')
GATES=tuple('current_original_n1_'+key+'_full_mixed_source_join_installed' for key,*_ in JOINS)
OPEN=tuple(key for key in previous.OPEN if key not in GATES)


def _restore_checked_layout(base,name,receipt):
    """Restore an accepted node layout on a freshly constructed owner only.

    Earlier symbolic queries may leave unused nodes after an owner's pinned
    prefix. Replaying those archived nodes preserves source IDs without
    replaying the queries. Subsequent constructors still derive and compare
    their own source snapshots; no output value is assigned here.
    """
    raw,rec=parent.linear.read(name),parent.linear.read(receipt)
    if not rec['all_passed'] or raw['source_family']!=base.family:
        raise ValueError('Same-family accepted source layout required')
    for filename,digest in rec['input_hashes'].items():
        if sha(filename)!=digest:raise ValueError('Changed layout input '+filename)
    n=len(base.prefix)
    if raw['exact_graph_nodes'][:n]!=base.g.nodes[:n]:raise ValueError('Exact ancestor prefix required')
    base.g.nodes=raw['exact_graph_nodes']
    base.g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(base.g.nodes)}
    if len(base.g.keys)!=len(base.g.nodes):raise ValueError('Unique exact source node layout required')
    base.assert_graph()


class CurrentOriginalAnnularMixedSourceJoins:
    def __init__(self,base=None,require_checked=True):
        if base is None:
            high=previous.previous.CurrentOriginalPrescribedHighRadialSourceFunctions()
            _restore_checked_layout(high,previous.NAME,previous.RECEIPT)
            base=previous.CurrentOriginalSwitchMixedSourceJoins(high)
            if require_checked:_restore_checked_layout(base,NAME,RECEIPT)
        self.base=base
        if not self.base.acceptance_loaded:raise ValueError('Accepted actual switch source interfaces required')
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        self.joins={};self.background_endpoints={}
        for key,left,right,endpoint in JOINS:
            packet=dict(left_chart=left,right_chart=right,original_endpoint=endpoint)
            bg={}
            for side,name,value in (('left',left,g.one),('right',right,g.zero)):
                chart=self.charts[name];x=chart['coordinate'];a=downstream_source.Algebra(g,6)
                packet[side]=self.base.at_own(name,value)
                packet[side+'_radius']=outer.Algebra(g,0).at(outer.FunctionJet((chart['radius'],)),x,value)
                bg[side]={group:{field:a.at(self.jet(ids),x,value) for field,ids in self.backgrounds[name][group].items()}
                    for group in ('actual','moments')}
            self.joins[key]=packet;self.background_endpoints[key]=bg
        self.annular_controls={}
        for name in ('long_reshape','restoration'):
            x=self.charts[name]['coordinate'];D=self.factories[name]
            sigma=[downstream_source.sigma(g,x)];prime=[downstream_source.sigma_prime(g,x)]
            for _ in range(5):sigma.append(D.derivative(sigma[-1]))
            for _ in range(4):prime.append(D.derivative(prime[-1]))
            self.annular_controls[name]=dict(sigma=sigma,prime=prime)
        self.definition=dict(actual_four_original_interfaces_to_postrestore=True,
            complete_ten_chart_high_sources_reused=True,
            actual_sigma_y5_and_sigma_prime_y4_from_original_flat_logistic_formula=True,
            actual_R110_B_anchored_G_and_F0_gauge_cancellation=True,
            exact_finite_T_reference_length_offsets_and_original_E_retained=True,
            nonzero_incoming_background_correction_and_signed_pressure_memory=True,
            source_background_y5_leading_y4_and_984_own_mixed_rows_only=True,
            numerical_oracles_global_n1_and_temporal_recursion_not_installed=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()

    def __getattr__(self,name):return getattr(self.base,name)
    def _snapshot(self):
        return outer.encode(dict(joins=self.joins,background_endpoints=self.background_endpoints,
            controls=self.annular_controls,contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self._snapshot()!=self.snapshot or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):
            raise ValueError('Actual annular source interfaces changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed annular source input '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=6
            or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['annular_mixed_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual annular mixed source interfaces required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalAnnularMixedSourceJoins(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        annular_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_FOUR_ANNULAR_MIXED_SOURCE_INTERFACES_DEFINED',len(owner.prefix),flush=True)
    return owner,raw
