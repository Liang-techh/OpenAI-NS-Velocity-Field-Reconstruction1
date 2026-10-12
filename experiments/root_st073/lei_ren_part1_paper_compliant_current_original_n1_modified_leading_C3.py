"""Same original finite-N velocity Z3 and own-mean endpoint increments.

Consume the accepted signed C3 loop/transport functions and their directed
ranges. This exposes corrected fields, not derivatives of magnitude caps.
Radial/phase derivatives and the complete n1 pressure tails remain open.
"""
import copy
from dataclasses import dataclass
from functools import wraps
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C3_target_ranges as parent
import lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension as cutoff

c3=parent.target;current=parent.current;source=current.source
HERE,PREFIX,sha,ep=parent.HERE,parent.PREFIX,parent.sha,parent.ep
NAME=PREFIX+'current_original_n1_modified_leading_C3.json.gz'
RECEIPT=PREFIX+'current_original_n1_modified_leading_C3_check.json'
GATE='current_original_modified_finite_N_velocity_Z3_and_own_mean_endpoint_increment_installed'
OPEN=tuple(dict.fromkeys((*parent.OPEN,*cutoff.OPEN)))


def source_precision(fn):
    @wraps(fn)
    def call(*args,**kwargs):
        with mp.workdps(540):return fn(*args,**kwargs)
    return call


def fingerprint(value):
    return json.dumps(c3.encoded(value),sort_keys=True,separators=(',',':'))


def corrected_functions(graph,roots,F,B,outgoing,incoming,N):
    alg=c3.C3Algebra(graph);eps=alg.fixed(graph.quotient(graph.one,N,'same original positive integer N'))
    dE=alg.mul(eps,F);dV=alg.mul(eps,B)
    return dict(E_N=alg.add(roots['E'],dE),V_N=alg.add(roots['V'],dV),
        delta_E=dE,delta_V=dV,own_mean_increment_at_right=alg.mul(eps,outgoing),
        own_mean_increment_at_left=alg.mul(eps,incoming))


def bounded_corrected_rows(c,primitive,density,mean,N0):
    bd=parent.Bounds(c);eps=bd.constant(c.mpf(1)/N0)
    dE=bd.scaled(density['F'],eps);dV=bd.scaled(primitive['B'],eps)
    return dict(E_N=bd.add(primitive['E'],dE),V_N=bd.add(primitive['V'],dV),
        delta_E=dE,delta_V=dV,own_mean_increment_at_right=bd.scaled(mean,eps))


@dataclass(frozen=True,eq=False)
class OriginalModifiedLeadingC3Packet:
    leading_hierarchy_order:int=0


class CurrentOriginalN1ModifiedLeadingC3:
    @source_precision
    def __init__(self,accepted_ranges,require_checked=True):
        if type(accepted_ranges) is not parent.CurrentC3TargetRanges or not accepted_ranges.acceptance_loaded:
            raise ValueError('Accepted signed original C3 range owner required')
        self.ranges=accepted_ranges;self.target=accepted_ranges.target;self.c=accepted_ranges.c
        self.family=copy.deepcopy(accepted_ranges.identity);self.hashes=dict(accepted_ranges.hashes)
        self.N0=accepted_ranges.phase.outer.N0;self._N0=self.N0
        self._ranges_hashes=copy.deepcopy(accepted_ranges.hashes)
        receipt=json.loads((HERE/cutoff.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[cutoff.GATE] or receipt['source_family']!=self.family:
            raise ValueError('Same accepted original n1 cutoff family required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Original cutoff source changed '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original source lineages differ '+name)
            self.hashes[name]=digest
        self.saved=json.loads(gzip.decompress((HERE/parent.NAME).read_bytes()))
        if self.saved['source_family']!=self.family:raise ValueError('Same original C3 transport family required')
        for name in (parent.NAME,parent.RECEIPT,cutoff.RECEIPT,Path(cutoff.__file__).name,Path(__file__).name):self.hashes[name]=sha(name)
        if current.RATES['m']!=1:raise ValueError('Own normalized axial mean must have original rate one')
        g=self.target.phase.built['graph'];self.graph=g;self.N=self.target.phase.built['N']
        self.prefix=copy.deepcopy(g.nodes);self.functions={};alg=c3.C3Algebra(g)
        for window in self.target.functions['windows']:
            chart=window['chart'];primitive=self.target.primitives[chart]
            F=window['F_N'] if window['F_N'] is not None else alg.fixed(0)
            self.functions[chart]=corrected_functions(g,primitive['source_roots'],F,primitive['B'],
                window['outgoing']['m'],window['incoming']['m'],self.N)
        self._identities=(self.ranges,self.target,self.c,self.graph,self.N)
        self._graph_length=len(g.nodes);self._graph=fingerprint(g.nodes)
        self._functions=fingerprint(self.functions);self._saved=fingerprint(self.saved)
        self._family=copy.deepcopy(self.family);self._packets={};self.acceptance_loaded=False
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed'] or not checked[GATE] or checked['source_family']!=self.family or any(checked[k] for k in OPEN):
                raise ValueError('Modified leading C3 scope/receipt differs')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Modified C3 source changed '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if ((self.ranges,self.target,self.c,self.graph,self.N)!=self._identities or
            self.target is not self.ranges.target or not self.ranges.acceptance_loaded or not self.target.acceptance_loaded or
            self.graph is not self.target.phase.built['graph'] or self.N is not self.target.phase.built['N'] or
            self.N0!=self._N0 or self.ranges.phase.outer.N0!=self._N0 or self.ranges.hashes!=self._ranges_hashes or
            self.family!=self._family or self.family!=self.target.identity or self.hashes!=self._hashes or
            len(self.graph.nodes)<self._graph_length or fingerprint(self.graph.nodes[:self._graph_length])!=self._graph or
            fingerprint(self.functions)!=self._functions or
            fingerprint(self.saved)!=self._saved):raise ValueError('Modified leading C3 source graph changed')

    @source_precision
    def endpoint(self,ends,chart):
        self.assert_graph();began=time.monotonic();ends=tuple(ends)
        if chart not in self.functions:raise ValueError('Original admitted chart required')
        transport=next((v for v in self.saved['actual_four_Z_C3_target_transports'] if tuple(v['exact_Z_cell'])==ends),None)
        if transport is None:raise ValueError('Original admitted whole-Z cell required')
        stored=next(v for v in transport['actual_17_chart_C3_transports'] if v['chart']==chart)
        right=current.native_partitions()[chart][-1]
        cell=self.ranges.cell(ends,chart,right,right)
        mean=parent.JetBound(*(parent.LogUpper(self.c,current.packets.interval(self.c,row['log_absolute_upper']))
            if not row.get('exact_zero',False) else parent.LogUpper.constant(self.c,0)
            for row in stored['actual_N_scaled_outgoing_C3']['m']))
        bounds=bounded_corrected_rows(self.c,cell['primitive'],cell['density'],mean,self.N0)
        view=dict(source_family=self.family,exact_Z_cell=ends,chart=chart,exact_native_right=right,
            corrected_signed_velocity_and_own_mean_C3_functions=c3.encoded(self.functions[chart]),
            corrected_source_C3_magnitude_ranges=parent.record(bounds),
            frequency=dict(original_N_function_node=self.N.node,uniform_for_integer_N_ge_N0=True,N0=self.N0,
                N_not_reselected_or_materialized=True),
            identities=dict(E_N='E+F_N/N=E*exp(A/N)',V_N='V+B/N',
                own_mean_increment='outgoing N-scaled m history divided by N',
                own_mean_ODE='d_y delta_m=delta_V-delta_m; zero at r_minus'),
            source_units='same common velocity unit S as original E,V; physical mean increment=S*delta_m',
            incoming_history_and_global_phase_not_reset=True,
            unchanged_original_axis_pressure_retained=True,
            complete_modified_Q_or_pressure_tail_not_installed=True,
            phase_radial_and_mixed_derivatives_not_exported=True,
            magnitude_ranges_are_not_point_values_or_derivatives_of_caps=True,
            endpoint_mean_increment_not_whole_cell_mean_oracle=True,
            physical_global_field_or_n1_repair_not_implied=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalModifiedLeadingC3Packet();self._packets[id(packet)]=(packet,view,cell,bounds,
            fingerprint((view,parent.record(cell),parent.record(bounds))))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalModifiedLeadingC3Packet or entry is None or entry[0] is not packet:
            raise ValueError('Live modified leading C3 packet required')
        if fingerprint((entry[1],parent.record(entry[2]),parent.record(entry[3])))!=entry[4]:
            raise ValueError('Modified leading C3 packet/range changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


def run(accepted_ranges=None):
    with mp.workdps(540):
        ranges=accepted_ranges if accepted_ranges is not None else parent.CurrentC3TargetRanges()
        owner=CurrentOriginalN1ModifiedLeadingC3(ranges,require_checked=False)
        packets={chart:owner.endpoint(('0','.5'),chart) for chart in ('first_micro','second_micro')}
        raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
            corrected_functions=c3.encoded(owner.functions),new_graph_nodes=owner.graph.nodes[len(owner.prefix):],
            original_C3_graph_prefix_length=len(owner.prefix),input_hashes=owner.hashes,
            **{GATE:False},**dict.fromkeys(OPEN,False))
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(c3.encoded(raw),separators=(',',':'))+'\n').encode(),mtime=0))
        print('MODIFIED_LEADING_C3_COMPUTED',len(packets),flush=True);return owner,packets
