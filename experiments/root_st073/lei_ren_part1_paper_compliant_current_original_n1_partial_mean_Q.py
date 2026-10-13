"""Actual finite-N partial own mean Z3 and radial-Q increment Z2.

The signed Duhamel function has a free query endpoint distinct from the
integration coordinate. Whole-cell bounds use the original positive logR
measure, incoming history and exact selected N, without materializing N.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_modified_leading_C3 as parent

c3,current,source=parent.c3,parent.current,parent.source
rg=parent.parent;HERE,PREFIX,sha,ep=parent.HERE,parent.PREFIX,parent.sha,parent.ep
source_precision=parent.source_precision;fingerprint=parent.fingerprint;OPEN=parent.OPEN
NAME=PREFIX+'current_original_n1_partial_mean_Q.json.gz'
RECEIPT=PREFIX+'current_original_n1_partial_mean_Q_check.json'
GATE='current_original_finite_N_partial_own_mean_Z3_mean_y_Z3_and_Q_increment_Z2_installed'


def at_native(g,body,chart,query):
    return g.node('substitute_original_native_coordinate',body=body.node,
        native_variable='native_'+chart,query_coordinate=query.node,
        definition='evaluate original body at native=query while holding Z and original N fixed',
        same_original_source_and_global_phase=True,
        defining_source=current.ast_binding(at_native))


def partial_mean_functions(g,base,window,velocity,N):
    """Original own-rate integral, with a separate free upper coordinate."""
    chart=base['chart'];variable='native_'+chart;query=g.symbol('partial_mean_query_'+chart)
    ref=lambda n:source.FunctionRef(g,n);alg=c3.C3Algebra(g)
    offset=at_native(g,ref(base['offset']),chart,query)
    memory=g.unary('exp',g.neg(g.sub(offset,ref(base['left_offset']))))
    kernel=g.unary('exp',g.neg(g.sub(offset,ref(base['offset']))))
    density=window['density']['m'];incoming=window['incoming']['m']
    scaled=[]
    for k,(local,before) in enumerate(zip(c3.rows(density),c3.rows(incoming))):
        integral=(g.zero if local==g.zero else current.controls.integral(g,
            g.mul(kernel,local,ref(base['Jacobian'])),variable,ref(base['domain'][0]),query,
            measure='original dy=Jacobian*dnative exactly once',N_scaled_density=True,
            ordinary_slow_Z_derivative_order=k,free_query_endpoint=query.node,
            incoming_mean_not_reset=True,query_kernel_endpoints_phase_Z_independent=True,
            original_native_source_partitions=base['actual_current_source_cells']))
        scaled.append(g.add(g.mul(memory,before),integral))
    eps=alg.fixed(g.quotient(g.one,N,'same selected exact positive N'))
    mean=alg.mul(eps,c3.C3Function(*scaled))
    V=c3.C3Function(*(at_native(g,v,chart,query) for v in c3.rows(velocity['delta_V'])))
    return dict(query=query,mean=mean,mean_y=alg.add(V,alg.neg(mean)),delta_V=V,
        original_N_scaled_mean=c3.C3Function(*scaled),logR_offset_at_query=offset)


def recover_Q_increment(g,V,m,delta):
    alg=c3.phase.C2Algebra(g);low=lambda q:c3.phase.C2Function(*c3.rows(q)[:3])
    z=c3.phase.C2Function(g.symbol('Z'),g.one,g.zero);d=alg.add(alg.fixed(1),alg.neg(alg.mul(z,z)))
    L=alg.add(alg.fixed(1),alg.neg(alg.mul(alg.fixed(delta),alg.mul(z,z))))
    mz=c3.phase.C2Function(m.Z,m.ZZ,m.ZZZ)
    numerator=alg.add(alg.scale(alg.mul(z,low(V)),2),
        alg.neg(alg.mul(alg.mul(alg.add(alg.fixed(1),alg.neg(alg.fixed(delta))),z),low(m))),
        alg.neg(alg.mul(d,mz)))
    return alg.div(numerator,L)


def coefficient_Q_bounds(c,ends,delta,V,m):
    z=IntervalTaylor.variable(c,c.mpf(ends),2);L=1-z*z*delta;d=1-z*z
    coefficients=(2*z/L,-z*(1-delta)/L,-d/L);bd=rg.Bounds(c)
    vr,mr=rg.rows(V),rg.rows(m);out=[]
    import math
    for k in range(3):
        pieces=[]
        for j in range(k+1):
            for coefficient,row in zip(coefficients,(vr,mr,mr[1:])):
                magnitude=c.mpf(ep(abs(coefficient[j])*math.factorial(j))[1])
                pieces.append(bd.constant(math.comb(k,j)*magnitude)*row[k-j])
        out.append(bd.sum(*pieces))
    return out


def rational(value):
    if (not isinstance(value,tuple) or len(value)!=2 or
        any(type(v) is not int for v in value) or value[1]<=0):
        raise ValueError('Exact rational native tuple required')
    return Fraction(*value)


def clipped_prefix(chart,right):
    if chart not in ('first_micro','second_micro'):raise ValueError('Numerical partial-mean micro charts only')
    end=rational(right);parts=current.native_partitions()[chart];cells=[]
    for a,b in zip(parts,parts[1:]):
        upper=rational(tuple(b))
        if a==('sc_half',):
            lower=None
        else:lower=rational(tuple(a))
        if lower is not None and lower>=end:break
        clipped=min(upper,end)
        if lower is None or clipped>lower:
            cells.append((tuple(a),(clipped.numerator,clipped.denominator)))
        if upper>=end:break
    return cells


@dataclass(frozen=True,eq=False)
class OriginalPartialMeanQPacket:
    leading_hierarchy_order:int=0


class CurrentOriginalN1PartialMeanQ:
    @source_precision
    def __init__(self,accepted_modified,require_checked=True):
        if type(accepted_modified) is not parent.CurrentOriginalN1ModifiedLeadingC3 or not accepted_modified.acceptance_loaded:
            raise ValueError('Accepted same-family finite-N C3 velocity owner required')
        accepted_modified.assert_graph();self.parent=accepted_modified;self.g=accepted_modified.graph
        self.c=accepted_modified.c;self.N=accepted_modified.N;self.family=copy.deepcopy(accepted_modified.family)
        self.hashes=dict(accepted_modified.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        logP=accepted_modified.target.phase.built['parameters']['logP']
        self.delta_node=self.g.unary('exp',self.g.neg(self.g.add(self.g.mul(self.g.constant(4),logP),self.g.constant(30))))
        self.delta=self.c.exp(-4*(self.c.exp(40)+11)-30)
        if not 0<ep(self.delta)[0]<=ep(self.delta)[1]<self.c.mpf('.5'):
            raise ValueError('Same original delta recipe must lie in(0,1/2)')
        node=self.g.nodes[self.N.node];selection=current.DoubleDyadicInteger(node['J'])
        self.logN=selection.log_interval(self.c)
        if ep(self.logN)[0]<ep(self.c.ln(accepted_modified.N0))[1]:raise ValueError('Original selected N must exceed input N0')
        self.epsilon_N=rg.LogUpper(self.c,-self.c.mpf(ep(self.logN)[0]));self.functions={}
        for base,window in zip(accepted_modified.target.phase.windows.values(),accepted_modified.target.functions['windows']):
            if base['chart']!=window['chart']:raise ValueError('Original chart ordering required')
            chart=base['chart'];got=partial_mean_functions(self.g,base,window,accepted_modified.functions[chart],self.N)
            got['delta_Q_Z2']=recover_Q_increment(self.g,got['delta_V'],got['mean'],self.delta_node)
            self.functions[chart]=got
        self._identities=(self.parent,self.g,self.c,self.N);self._length=len(self.g.nodes)
        self._graph=fingerprint(self.g.nodes);self._functions=fingerprint(self.functions)
        self._delta=self.delta._mpi_;self._logN=self.logN._mpi_;self._epsilon=fingerprint(self.epsilon_N.record())
        self._family=copy.deepcopy(self.family);self._packets={};self.acceptance_loaded=False
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed'] or not checked[GATE] or checked['source_family']!=self.family or any(checked[k] for k in OPEN):
                raise ValueError('Partial mean/Q scope or receipt differs')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Partial mean source changed '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if ((self.parent,self.g,self.c,self.N)!=self._identities or self.g is not self.parent.graph or self.N is not self.parent.N or
            not self.parent.acceptance_loaded or self.family!=self._family or self.family!=self.parent.family or self.hashes!=self._hashes or
            self.delta._mpi_!=self._delta or self.logN._mpi_!=self._logN or fingerprint(self.epsilon_N.record())!=self._epsilon or
            len(self.g.nodes)<self._length or fingerprint(self.g.nodes[:self._length])!=self._graph or fingerprint(self.functions)!=self._functions):
            raise ValueError('Partial mean/Q original graph changed')
        self.parent.assert_graph()

    @source_precision
    def range_cell(self,ends,chart,left,right):
        self.assert_graph();began=time.monotonic();ends=tuple(ends);l,r=rational(left),rational(right)
        limits={'first_micro':(Fraction(1,4),Fraction(1)), 'second_micro':(Fraction(1),Fraction(2))}
        if chart not in limits or not limits[chart][0]<=l<r<=limits[chart][1]:raise ValueError('Supported ordered actual micro cell required')
        transport=next((v for v in self.parent.saved['actual_four_Z_C3_target_transports'] if tuple(v['exact_Z_cell'])==ends),None)
        if transport is None:raise ValueError('Original whole-Z cell required')
        stored=next(v for v in transport['actual_17_chart_C3_transports'] if v['chart']==chart)
        bd=rg.Bounds(self.c);read=lambda row:(bd.zero if row.get('exact_zero',False) else rg.LogUpper(self.c,current.packets.interval(self.c,row['log_absolute_upper'])))
        scaled_mean=rg.JetBound(*(read(v) for v in stored['actual_N_scaled_incoming_C3']['m']))
        ledgers=[];sources=[]
        for a,b in clipped_prefix(chart,right):
            cell=self.parent.ranges.cell(ends,chart,a,b);weight=self.parent.ranges.weights(ends,chart,a,b,Fraction(0))
            width=bd.row(weight['original_mass']);scaled_mean=bd.add(scaled_mean,bd.scaled(cell['primitive']['B'],width))
            sources.append((cell,weight));ledgers.append(dict(exact_left=a,exact_right=b,original_positive_logR_width=width.record(),
                original_B_Z3_magnitude=rg.record(cell['primitive']['B']),rate_one_kernel_upper='<=1 for source before its query',
                incoming_decay_upper='<=1; predecessor not reset'))
        mean=bd.scaled(scaled_mean,self.epsilon_N)
        local=self.parent.ranges.cell(ends,chart,left,right);deltaV=bd.scaled(local['primitive']['B'],self.epsilon_N)
        mean_y=bd.add(deltaV,mean);Q=coefficient_Q_bounds(self.c,ends,self.delta,deltaV,mean)
        view=dict(source_family=self.family,exact_Z_cell=ends,chart=chart,exact_native_cell=[left,right],
            signed_partial_mean_query_functions=c3.encoded(self.functions[chart]),
            whole_cell_partial_own_mean_increment_Z3_magnitude=rg.record(mean),
            whole_cell_mean_y_increment_Z3_magnitude=rg.record(mean_y),
            whole_cell_Q_increment_Z2_magnitude=[v.record() for v in Q],
            selected_N_log=self.logN,epsilon_N_log_bound=self.epsilon_N.record(),original_delta_recipe='exp(-4*(exp(40)+11)-30)',
            original_incoming_N_scaled_mean_Z3=stored['actual_N_scaled_incoming_C3']['m'],partial_history_cells=ledgers,
            common_unit='physical increments are S*delta_m and S*delta_Q, original common S retained',
            same_selected_N_division_once_no_integer_materialization=True,
            partial_mean_free_query_is_distinct_from_bound_integration_coordinate=True,
            whole_cell_range_not_outgoing_endpoint_value=True,
            modified_Q_increment_only_not_complete_Q=True,
            internal_radial_phase_derivative_of_velocity_not_exported=True,
            complete_pressure_source_tail_support_and_repair_not_installed=True,
            magnitude_bounds_not_point_values_and_no_relative_width_claim=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalPartialMeanQPacket();state=(mean,mean_y,Q,scaled_mean,deltaV)
        self._packets[id(packet)]=(packet,view,state,sources,local,
            fingerprint((view,rg.record(state),[(rg.record(c),rg.record(w)) for c,w in sources],rg.record(local))))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalPartialMeanQPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live original partial mean/Q packet required')
        if fingerprint((entry[1],rg.record(entry[2]),[(rg.record(c),rg.record(w)) for c,w in entry[3]],rg.record(entry[4])))!=entry[5]:
            raise ValueError('Partial mean/Q packet or original ranges changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_modified=None):
    if accepted_modified is None:accepted_modified=parent.CurrentOriginalN1ModifiedLeadingC3(rg.CurrentC3TargetRanges())
    owner=CurrentOriginalN1PartialMeanQ(accepted_modified,require_checked=False)
    packets={chart:owner.range_cell(('0','.5'),chart,a,b) for chart,a,b in (
        ('first_micro',(1,4),(1,2)),('second_micro',(5,4),(3,2)))}
    raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
        signed_functions=c3.encoded(owner.functions),input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(c3.encoded(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('PARTIAL_MEAN_Q_COMPUTED',len(packets),flush=True);return owner,packets
