"""Same-source n1 common-inner geometry below the later modification starts.

The selected sc is restored from its exact original binary endpoint.
hb remains the original formal width. This selects an admissible inner
interval; it does not install the incomplete completed-leading dispatch.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step as coupled
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=coupled.op;domain=coupled.domain
HERE,PREFIX,sha=coupled.HERE,coupled.PREFIX,coupled.sha
NAME=PREFIX+'current_original_n1_common_inner_selection.json.gz'
RECEIPT=PREFIX+'current_original_n1_common_inner_selection_check.json'
GATE='current_original_n1_common_inner_geometry_below_known_modification_starts_selected'
OPEN=coupled.OPEN;ep=coupled.ep


def read(c,row):
    if 'lower_exact_mpf_tuple' in row:
        return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                      mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))
    return c.mpf((row['lower'],row['upper']))


@dataclass(frozen=True,eq=False)
class OriginalN1CommonInnerPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1CommonInnerSelection:
    @source_precision
    def __init__(self,base,require_checked=True):
        if type(base) is not coupled.CurrentOriginalN1CoupledFirstStep or not base.acceptance_loaded:
            raise ValueError('Accepted genuine n1 coupled owner required')
        base.assert_graph();self.base=base;self.c=c=base.c;self.source=base.source
        self.family=copy.deepcopy(base.family);self.hashes=dict(base.hashes);self.records={}
        bind=domain.parent.parent.original.bind
        for stem,gate in (
            ('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified'),
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified')):
            receipt_name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/receipt_name).read_bytes())
            family=receipt.get('source_family') or {k:receipt[k] for k in self.family}
            if not receipt.get('all_passed') or not receipt.get(gate) or family!=self.family:
                raise ValueError('Same checked original inner/modification geometry required')
            for name,digest in {**receipt['input_hashes'],receipt_name:sha(receipt_name)}.items():
                if sha(name)!=digest:raise ValueError('Original geometry dependency changed '+name)
                bind(self.hashes,name,digest)
            name=PREFIX+stem+'.json';self.records[stem]=json.loads((HERE/name).read_bytes())
            bind(self.hashes,name,sha(name))
        for name in (coupled.NAME,coupled.RECEIPT,Path(__file__).name):bind(self.hashes,name,sha(name))
        # Its accepted radius/parameter frame binds basis slot 1 to
        # 2*logPstar, and Rm_factor to powers (0,5,0,0,0).
        name=PREFIX+'current_original_Rm_conditioned_phase_density_check.json'
        rm=json.loads((HERE/name).read_bytes())
        if not (rm.get('all_passed') and rm.get('original_actual_Rm_conditioned_Z_phase_and_five_signed_density_interface_installed')
                and rm['source_family']==self.family['actual_five_defect_family_sha256']):
            raise ValueError('Same accepted original Rm parameter/radius frame required')
        for path,digest in {**rm['input_hashes'],name:sha(name)}.items():
            if sha(path)!=digest:raise ValueError('Original Rm geometry dependency changed '+path)
            bind(self.hashes,path,digest)
        self.proof=self.records['current_inner_exit_strict_collar']['explicit_current_inner_exit_strict_collar']
        self.loop=self.records['current_generic_shear_loop_domain']['current_original_generic_loop_domain']
        bits=self.proof['selected_first_phase_endpoint']
        if bits['lower_exact_mpf_tuple']!=bits['upper_exact_mpf_tuple']:
            raise ValueError('Exact original selected sc, not a new point from a cover, required')
        self.sc=read(c,bits)
        if not 0<ep(self.sc)[0]==ep(self.sc)[1]<c.mpf('.25'):
            raise ValueError('Original exact positive sc below one quarter required')
        if self.loop['formal_radii']['r_minus']!='Ra*exp(hb*s_c/2)' or self.loop['source_family']!=self.family:
            raise ValueError('Original shear modification origin differs')
        if ep(read(c,self.proof['positive_margins']['subcollar_inside_original_rho4point1_continuation']))[0]<=0:
            raise ValueError('Original hb sc continuation bound required')
        if self.records['current_inner_exit_strict_collar']['exact_current_exit_source_attachment']['exact_width']!=\
                'hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum':
            raise ValueError('Original formal width definition differs')
        self.phase_in=self.sc/4;self.phase_keep=self.sc/8;self.phase_minus=self.sc/2
        # The already checked hb*sc<log(4.1/4) bounds only the retained
        # positive offset. It never supplies a replacement hb or radius.
        core=self.source.core
        self.rm_log_gap_lower=domain.lower(c,10*core.logC+14*core.logP+c.ln(c.mpf(110)/4)+994
                                          -c.ln(c.mpf('4.1')/4)/4)
        if ep(self.rm_log_gap_lower)[0]<=0:raise ValueError('Original Rm is not strictly outside selected Rin')
        self._base=base;self._c=c;self._source=self.source;self._family=copy.deepcopy(self.family)
        self._sc=self.sc._mpi_;self._rm_gap=self.rm_log_gap_lower._mpi_
        self._records=copy.deepcopy(self.records);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Common-inner selection receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Common-inner dependency changed '+name)
                bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.base is not self._base or self.c is not self._c or self.source is not self._source
            or self.source is not self.base.source or self.c is not self.base.c or not self.base.acceptance_loaded
            or self.family!=self._family or self.family!=self.base.family or self.sc._mpi_!=self._sc
            or self.rm_log_gap_lower._mpi_!=self._rm_gap
            or self.proof is not self.records['current_inner_exit_strict_collar']['explicit_current_inner_exit_strict_collar']
            or self.loop is not self.records['current_generic_shear_loop_domain']['current_original_generic_loop_domain']
            or self.records!=self._records or self.hashes!=self._hashes
            or self.phase_in!=self.sc/4 or self.phase_keep!=self.sc/8 or self.phase_minus!=self.sc/2):
            raise ValueError('Common-inner original source/geometry changed')
        return self.base.assert_graph()

    @source_precision
    def geometry(self):
        self.assert_graph();c=self.c
        return dict(source_family=self.family,
            selected_sc=self.sc,selected_sc_exact_mpf_tuple=self.proof['selected_first_phase_endpoint']['lower_exact_mpf_tuple'],
            formal_radii=dict(Ra='4*epsilon_core',Rkeep='Ra*exp(hb*s_c/8)',
                              Rin='Ra*exp(hb*s_c/4)',r_minus='Ra*exp(hb*s_c/2)'),
            actual_first_phases=dict(Rkeep=self.phase_keep,Rin=self.phase_in,r_minus=self.phase_minus),
            formal_log_radius_offsets={name:dict(exact_source='hb*'+fraction+'*s_c',
                log_positive_offset=self.source.logh+c.ln(self.sc)-c.ln(denominator))
                for name,fraction,denominator in (('Rkeep','1/8',8),('Rin','1/4',4),('r_minus','1/2',2))},
            exact_positive_gap_relations=dict(log_Rkeep_over_Ra='hb*s_c/8',
                log_Rin_over_Rkeep='hb*s_c/8',log_r_minus_over_Rin='hb*s_c/4'),
            original_Rm_log_radius='10*logCstar+10*logPstar+log110-6',
            exact_log_Rm_over_Rin='10*logCstar+14*logPstar+log(110/4)+994-hb*s_c/4',
            log_Rm_over_Rin_strict_lower=self.rm_log_gap_lower,
            original_Rm_basis_uses_5_times_2_logPstar=True,
            strict_Ra_Rkeep_Rin_r_minus_order=True,original_hb_preserved_not_materialized=True,
            old_phase_half_inner_radius_exceeded_shear_start=True,
            exact_original_sc_restored_not_selected_from_numeric_enclosure=True,
            canonical_actual_first_collar_bounds_restrict_to_selected_interval=True,
            completed_leading_dispatch_still_open=True,
            next_support_binding='Bind every assembled leading modification and its zero extension to this same source below Rin',
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',reuse_packet=None):
        """First chart uses fraction of the corrected common-inner phase."""
        self.assert_graph();c=self.c;value=c.mpf(coordinate)
        if chart=='core':
            if ep(value)[0]<0 or ep(value)[1]>4:raise ValueError('Actual core coordinate subset[0,4] required')
            actual=value
        elif chart=='first':
            if ep(value)[0]<0 or ep(value)[1]>1:raise ValueError('Common-inner fraction subset[0,1] required')
            actual=value*self.phase_in
        else:raise ValueError('Only common-inner core and actual first charts supported')
        upstream=self.base.evaluate(Z,actual,chart) if reuse_packet is None else reuse_packet
        view=self.base.report(upstream)
        if view['Z']!=c.mpf(Z) or view['endpoint']['chart']!=chart or view['endpoint']['coordinate']!=actual:
            raise ValueError('Reused coupled packet must have this exact common-inner query')
        result=dict(common_inner_geometry=self.geometry(),query=dict(chart=chart,coordinate=value,
            actual_source_coordinate=actual),canonical_n1_field=view,
            field_scope='canonical unrepaired n1 field restricted to the corrected common-inner interval; completed-leading dispatch unbound',
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalN1CommonInnerPacket();self._packets[id(packet)]=(packet,result,upstream,op.snapshot(result))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1CommonInnerPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live common-inner packet from this owner required')
        self.base.report(entry[2])
        if op.snapshot(entry[1])!=entry[3]:raise ValueError('Common-inner field packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(base,reuse_packet=None):
    began=time.monotonic();owner=CurrentOriginalN1CommonInnerSelection(base,require_checked=False)
    packets={'corrected_common_inner_exit':owner.evaluate('.371','1','first',reuse_packet=reuse_packet)}
    raw=dict(source_family=owner.family,geometry=owner.geometry(),
        samples={name:owner.report(packet) for name,packet in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_CORRECTED_COMMON_INNER_SAMPLE_COMPUTED',len(packets),flush=True)
    return owner,packets
