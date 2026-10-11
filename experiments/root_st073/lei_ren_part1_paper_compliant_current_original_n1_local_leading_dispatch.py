"""Executable same-source unchanged leading branch on selected n1 inner.

Replay the accepted original Section11 support selector with its original
sc and selected eta. Preserve total background fields/pressure/histories;
only later modification increments vanish. Global completion stays open.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace

import lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection as selection
import lei_ren_part1_paper_compliant_current_original_micro_finite_N as original
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=selection.op;domain=selection.domain
HERE,PREFIX,sha=selection.HERE,selection.PREFIX,selection.sha
NAME=PREFIX+'current_original_n1_local_leading_dispatch.json.gz'
RECEIPT=PREFIX+'current_original_n1_local_leading_dispatch_check.json'
GATE='current_original_n1_selected_inner_local_leading_branch_identified'
OPEN=selection.OPEN;ep=selection.ep


@dataclass(frozen=True,eq=False)
class OriginalN1LocalLeadingPacket:
    temporal_hierarchy_order:int=0


class CurrentOriginalN1LocalLeadingDispatch:
    @source_precision
    def __init__(self,inner,require_checked=True):
        if type(inner) is not selection.CurrentOriginalN1CommonInnerSelection or not inner.acceptance_loaded:
            raise ValueError('Accepted corrected same-source n1 inner selection required')
        inner.assert_graph();self.inner=inner;self.inputs=inner.base.inputs;self.source=inner.source
        self.c=c=inner.c;self.family=copy.deepcopy(inner.family);self.hashes=dict(inner.hashes);self.records={}
        bind=domain.parent.parent.original.bind
        for stem,gate in (
            ('current_original_micro_finite_N',original.GATE),
            ('current_generic_shear_uniform_inputs','current_original_whole_generic_input_margins_and_log_scales_certified')):
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            expected=self.family['actual_five_defect_family_sha256'] if isinstance(receipt['source_family'],str) else self.family
            if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family']!=expected:
                raise ValueError('Same accepted original later-modification source required')
            for path,digest in {**receipt['input_hashes'],name:sha(name)}.items():
                if sha(path)!=digest:raise ValueError('Original leading support source changed '+path)
                bind(self.hashes,path,digest)
        name=PREFIX+'current_generic_shear_uniform_inputs.json'
        self.records['uniform']=json.loads((HERE/name).read_bytes());bind(self.hashes,name,sha(name))
        for name in (selection.NAME,selection.RECEIPT,Path(__file__).name):bind(self.hashes,name,sha(name))
        if self.records['uniform']['source_family']!=self.family:
            raise ValueError('Selected later eta source family differs')
        self.eta_log=selection.read(c,self.records['uniform']['current_actual_logarithmic_loop_scales']['selected_positive_eta_log'])
        loop=inner.loop;self.left_excess=selection.read(c,loop['left_collar_kappa_minus2_lower'])
        self.allowed_eta_log=selection.read(c,loop['allowed_eta_log_upper_for_automatic_constant_edges'])
        self.eta_gap=domain.lower(c,self.allowed_eta_log-self.eta_log)
        if ep(self.eta_log)[1]>ep(self.allowed_eta_log)[0] or ep(self.left_excess)[0]<=0:
            raise ValueError('Actual selected eta must satisfy original left flat collar')
        # Only the accepted selector is replayed. No macro/patch/heat
        # ancestor is constructed and no active inverse is evaluated.
        self.selector=SimpleNamespace(c=c,macro=SimpleNamespace(sc=inner.sc))
        self.support=original.OriginalMicroFiniteN.modification_support(self.selector,'first_micro',0,inner.phase_in)
        if not self.support['increments_exact_zero'] or self.support['region']!='unchanged_left':
            raise ValueError('Whole selected interval must use original unchanged-left branch')
        self.Rc_gap=domain.lower(c,10*self.source.core.logC+15*self.source.core.logP+c.ln(c.mpf(110)/4)+1003
                                -c.ln(c.mpf('4.1')/4)/4)
        if ep(inner.rm_log_gap_lower)[0]<=0 or ep(self.Rc_gap)[0]<=0:
            raise ValueError('Original leading and later repair radii must be outside selected Rin')
        self._inner=inner;self._inputs=self.inputs;self._source=self.source;self._c=c;self._family=copy.deepcopy(self.family)
        self._state=op.snapshot((self.records,self.eta_log,self.allowed_eta_log,self.left_excess,self.eta_gap,self.Rc_gap,self.support))
        self._selector_sc=inner.sc._mpi_;self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Local leading dispatch receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Local leading dispatch dependency changed '+name)
                bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if (self.inner is not self._inner or self.inputs is not self._inputs or self.source is not self._source
            or self.c is not self._c or self.inputs is not self.inner.base.inputs or self.source is not self.inner.source
            or self.c is not self.inner.c or self.family!=self._family or self.family!=self.inner.family
            or not self.inner.acceptance_loaded or self.hashes!=self._hashes
            or self.selector.c is not self.c or self.selector.macro.sc._mpi_!=self._selector_sc
            or self._state!=op.snapshot((self.records,self.eta_log,self.allowed_eta_log,self.left_excess,self.eta_gap,self.Rc_gap,self.support))):
            raise ValueError('Local leading original source/selector/eta state changed')
        return self.inner.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='core'):
        self.assert_graph();c=self.c;value=c.mpf(coordinate)
        if chart=='core':
            if ep(value)[0]<0 or ep(value)[1]>4:raise ValueError('Selected actual core subset[0,4] required')
            actual=value;support=dict(region='original_stress_free_core',increments_exact_zero=True,
                core_is_outside_Section11_loop_domain=True,core_cone_argument_not_used=True)
        elif chart=='first':
            if ep(value)[0]<0 or ep(value)[1]>1:raise ValueError('Selected first fraction subset[0,1] required')
            actual=value*self.inner.phase_in
            support=original.OriginalMicroFiniteN.modification_support(self.selector,'first_micro',0,actual)
            if not support['increments_exact_zero']:raise ValueError('Original changed source is outside this local owner')
        else:raise ValueError('Local leading dispatcher owns only selected core/first charts')
        packet=self.inputs.evaluate(Z,actual,chart);leading=self.inputs.report(packet)
        result=dict(source_family=self.family,Z=leading['Z'],query=dict(chart=chart,coordinate=value,actual_source_coordinate=actual),
            actual_leading_source=leading,executed_original_support_selector=support,
            later_velocity_pressure_increment={key:0 for key in ('F0','Uz0','Q0','P0')},
            later_five_moment_increment=dict.fromkeys(('m','h','k','e','p'),0),
            background_velocity_pressure_and_own_histories_retained=True,only_modification_increments_zeroed=True,
            selected_eta_log=self.eta_log,original_left_excess_lower=self.left_excess,selected_eta_admission_log_gap=self.eta_gap,
            flat_source_identity='Delta>=left_excess>=2eta => q=A=B=0 => E_N=E,V_N=V',
            any_positive_N_has_same_exact_local_branch=True,
            original_Rm_over_Rin_log_gap_lower=self.inner.rm_log_gap_lower,original_Rc_over_Rin_log_gap_lower=self.Rc_gap,
            no_additional_spatial_cutoff_introduced=True,source_graph_ancestor_constructors_called=False,
            active_Section11_inverse_evaluated=False,local_canonical_holomorphic_bounds_transfer_unchanged=True,
            source_point_histories_or_global_completed_dispatch_not_promoted=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1LocalLeadingPacket();self._packets[id(token)]=(token,result,packet,op.snapshot(result))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1LocalLeadingPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live local leading packet from this owner required')
        self.inputs.report(entry[2])
        if op.snapshot(entry[1])!=entry[3]:raise ValueError('Local leading packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(inner):
    began=time.monotonic();owner=CurrentOriginalN1LocalLeadingDispatch(inner,require_checked=False)
    packets={name:owner.evaluate(z,x,chart) for name,z,x,chart in (
        ('whole_core_axis',[-1,1],'0','core'),('core_exit','.371','4','core'),
        ('actual_core_first_inlet','.371','0','first'),('selected_actual_inner_exit','.371','1','first'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(packet) for name,packet in packets.items()},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_ORIGINAL_LOCAL_LEADING_DISPATCH_COMPUTED',len(packets),flush=True)
    return owner,packets
