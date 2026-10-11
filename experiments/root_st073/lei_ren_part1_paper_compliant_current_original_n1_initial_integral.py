"""Computed regular source-integral enclosures on the ACTUAL inner domain.

Wbase=(I+J)Gg integrates the known n1 source and the three kinematic
rows. This is a computed initial iterate, not a solution of the coupled
equation. Interval cells enclose entire integrands; no midpoint values,
representative widths or finite-axis polynomials are substituted.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_operator as operator
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=operator.HERE,operator.PREFIX,operator.sha
NAME=PREFIX+'current_original_n1_initial_integral.json.gz'
RECEIPT=PREFIX+'current_original_n1_initial_integral_check.json'
GATE='current_original_n1_initial_regular_source_integrals_computed'
OPEN=operator.OPEN
ep=operator.ep
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')


def positive(c,value):
    lo,hi=ep(value)
    if hi<0:raise ValueError('Positive exact integration kernel required')
    return c.mpf([max(lo,0),hi])


def log_kernel(c,cell,endpoint):
    """Enclose s log(x/s), including its exact zero-axis continuation.

    Concavity locates its only maximum at x/e. Endpoint intervals are
    numerical enclosures of exact geometry, not selected source values.
    """
    a,b=ep(cell);l,u=ep(endpoint)
    if a<0 or l<=0 or a>u:raise ValueError('Nonnegative inner integration cell required')
    def at(v,X):return c.mpf(0) if v==0 else c.mpf(v)*c.ln(c.mpf(X)/v)
    low=0 if a==0 or b>=l else min(ep(at(a,l))[0],ep(at(b,l))[0])
    high=max(ep(at(a,u))[1],ep(at(min(b,u),u))[1])
    vertex=c.mpf(u)/c.e
    if a<=ep(vertex)[1] and b>=ep(vertex)[0]:high=max(high,ep(vertex)[1])
    return c.mpf([max(0,low),max(0,high)])


def kernels(c,cell,endpoint):
    ratio=cell/endpoint
    tail=positive(c,1-ratio*ratio)
    return dict(F1=cell*tail/2,Uz1=log_kernel(c,cell,endpoint),
                K1=-cell*tail/2,P1=c.mpf(1),
                d_x_F1=ratio**3,d_x_Uz1=ratio)


def integrate_cell(algebra,system,cell,endpoint,measure,collar=False):
    """Integrate source + nilpotent kinematic feedback by exact kernels."""
    weights=kernels(algebra.ctx,cell,endpoint);g=system['g'];out={}
    for field,row in (('F1','5'),('Uz1','6'),('K1','6'),('P1','4'),
                      ('d_x_F1','5'),('d_x_Uz1','6')):
        value=g[row]*weights[field]*measure
        if collar:value=algebra.width(value)
        if any(key[0]<0 for key in value.terms):
            raise ValueError('Inverse microscopic width did not cancel before integration')
        out[field]=value
    return out


@dataclass(frozen=True,eq=False)
class OriginalN1InitialIntegralPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1InitialIntegral:
    @source_precision
    def __init__(self,inputs,require_checked=True):
        if type(inputs) is not operator.CurrentOriginalN1InnerOperator or not inputs.acceptance_loaded:
            raise ValueError('Accepted actual n1 input-enclosure owner required')
        inputs.assert_graph();self.inputs=inputs;self.source=inputs.source;self.c=inputs.c
        self.family=copy.deepcopy(inputs.family);self.hashes=dict(inputs.hashes)
        for name in (operator.NAME,operator.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._inputs=inputs;self._source=self.source;self._c=self.c;self._family=copy.deepcopy(self.family)
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Initial regular integral receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Initial n1 integral input changed '+name)
                operator.domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.inputs is not self._inputs or self.source is not self._source or self.c is not self._c
            or self.family!=self._family or self.family!=self.inputs.family or self.c is not self.inputs.c
            or self.source is not self.inputs.source or not self.inputs.acceptance_loaded or self.hashes!=self._hashes):
            raise ValueError('Initial n1 integral source graph changed')
        return self.inputs.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',core_cells=16,collar_cells=4):
        self.assert_graph();c=self.c;z,algebra,_=self.inputs.frame(Z);requested_coordinate=coordinate;coordinate=c.mpf(coordinate)
        if isinstance(core_cells,bool) or not isinstance(core_cells,int) or not 1<=core_cells<=4096:
            raise ValueError('A bounded positive integer core partition required')
        if isinstance(collar_cells,bool) or not isinstance(collar_cells,int) or not 1<=collar_cells<=4096:
            raise ValueError('A bounded positive integer collar partition required')
        if not isinstance(requested_coordinate,(str,int,float)) and ep(coordinate)[0]!=ep(coordinate)[1]:
            raise ValueError('One exact endpoint coordinate required; integrands use interval cells')
        if chart=='core':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>4:raise ValueError('Actual core endpoint in[0,4] required')
            endpoint=c.sqrt(coordinate);core_end=endpoint;phase=c.mpf(0)
        elif chart=='first':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>.5:raise ValueError('First collar endpoint phase in[0,1/2] required')
            end=self.inputs.evaluate(z,coordinate,'first');view=self.inputs.report(end)
            endpoint=c.sqrt(view['geometry']['rho_enclosure_only']);core_end=c.mpf(2);phase=coordinate
        else:raise ValueError('Only actual core and first collar supported')
        totals={name:algebra.lift(0,1) for name in FIELDS};ledger=[];began=time.monotonic()
        if ep(core_end)[1]>0:
            for j in range(core_cells):
                cell=core_end*c.mpf([c.mpf(j).a/core_cells,c.mpf(j+1).b/core_cells])
                packet=self.inputs.evaluate(z,cell*cell,'core');self.inputs.report(packet)
                system=self.inputs._packets[id(packet)][4]
                values=integrate_cell(algebra,system,cell,endpoint,core_end/core_cells)
                totals={name:totals[name]+values[name] for name in FIELDS}
                ledger.append(dict(chart='core',x_cell=cell,rho_cell=cell*cell,measure=core_end/core_cells))
        if ep(phase)[1]>0:
            for j in range(collar_cells):
                cell=phase*c.mpf([c.mpf(j).a/collar_cells,c.mpf(j+1).b/collar_cells])
                packet=self.inputs.evaluate(z,cell,'first');view=self.inputs.report(packet)
                xcell=c.sqrt(view['geometry']['rho_enclosure_only']);system=self.inputs._packets[id(packet)][4]
                # dx=(x/2)*hb*ds. Multiply by hb BEFORE enclosing.
                values=integrate_cell(algebra,system,xcell,endpoint,xcell*phase/(2*collar_cells),True)
                totals={name:totals[name]+values[name] for name in FIELDS}
                ledger.append(dict(chart='first',phase_cell=cell,x_enclosure_only=xcell,
                    exact_dx='(sqrt(4*exp(hb*s))/2)*hb*ds',width_power_added_before_enclosure=1))
        result=dict(source_family=self.family,Z=z,endpoint=dict(chart=chart,coordinate=coordinate,
            x_enclosure_only=endpoint,exact_x='sqrt(rho)' if chart=='core' else '2*exp(hb*s/2)'),
            temporal_hierarchy_order=1,computed_iterate='Wbase=(I+J)Gg',
            G='component i: integral_0^x (s/x)^D_i g_i(s,Z) ds',
            J='G times constant B0 kinematic rows (1,5)=1,(2,6)=1,(3,6)=-1; J^2=0',
            initial_regular_field_enclosures={name:operator.export(value,1) for name,value in totals.items()},
            ordinary_Z_derivative_order=1,source_log_basis_names=['hb','F0_anchor','Pstar_squared','Lambda'],
            source_log_bases=list(algebra.logs),integration_cells=ledger,
            rigorous_entire_cell_interval_integration=True,midpoint_quadrature_used=False,
            inverse_width_cancelled_before_enclosure=True,actual_own_moment_enclosures_used=True,
            exact_point_moment_history_recovered=False,core_finite_axis_polynomial_used=False,
            coupled_B0_B1_feedback_still_open=True,full_field_error_not_yet_certified=True,
            computed_initial_integral_is_not_full_n1_solution=True,
            execution_seconds=time.monotonic()-began,**dict.fromkeys(OPEN,False))
        packet=OriginalN1InitialIntegralPacket();self._packets[id(packet)]=(packet,result,totals,operator.snapshot(result),operator.snapshot(totals))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1InitialIntegralPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live initial n1 integral packet from this owner required')
        if operator.snapshot(entry[1])!=entry[3] or operator.snapshot(entry[2])!=entry[4]:
            raise ValueError('Initial n1 integral packet changed')
        return dict(**copy.deepcopy(entry[1]),**{GATE:self.acceptance_loaded})


@source_precision
def run(inputs):
    began=time.monotonic();owner=CurrentOriginalN1InitialIntegral(inputs,require_checked=False)
    packets={name:owner.evaluate('.371',value,chart) for name,chart,value in
             (('core_exit','core','4'),('actual_first_keep','first','.5'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(v) for name,v in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    data=json.dumps(operator.snapshot(raw),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('INITIAL_N1_REGULAR_INTEGRALS_COMPUTED',len(packets),flush=True)
    return owner,packets
