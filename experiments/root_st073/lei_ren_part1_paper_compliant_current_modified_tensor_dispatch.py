"""Actual modified33-chart tensor dispatcher and checked original power tail."""
import gzip
import hashlib
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor import (
    CurrentModifiedPrePhysicalTensor,HERE,PREFIX,NAME as LOCAL_NAME,RECEIPT as LOCAL_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,GATES as LOCAL_GATES,OPEN as LOCAL_OPEN)
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import modified_physical_packet
from lei_ren_part1_paper_compliant_current_modified_tensor_dispatch_operator import (
    original_signed_power,exact_modified_dispatch_theorem)

NAME=PREFIX+'current_modified_tensor_dispatch.json'
RECEIPT=PREFIX+'current_modified_tensor_dispatch_check.json'
VIEWS_NAME=PREFIX+'current_modified_tensor_dispatch_views.json.gz'
GATES=('current_modified_open_O3_power_continuation_source_certified',
       'current_modified_full_O3_power_normalized_phase_dispatch_available',
       'current_modified_power_q2_straddling_boxes_partitioned',
       'current_modified_full_33_chart_physical_tensor_dispatch_installed')
OPEN=tuple(key for key in LOCAL_OPEN if key not in GATES)

class CurrentModifiedTensorDispatch:
    @source_precision
    def __init__(self,local=None,require_checked=True):
        self.local=local if local is not None else CurrentModifiedPrePhysicalTensor()
        if type(self.local) is not CurrentModifiedPrePhysicalTensor or not self.local.acceptance_loaded:
            raise ValueError('Checked same-source modified local physical tensor required')
        self.ctx=c=self.local.ctx;self.delta=self.local.delta;self.N=self.local.N;self.mu=self.local.mu
        self.registry=self.local.registry;self.geometry=self.local.geometry
        self.pre=self.local.source.source.histories.pre;self.Tw=c.mpf(self.pre.params.Tw)
        self.phase_cut=endpoints(c.mpf('1.95')/c.mpf(endpoints(self.Tw)[1]))[1]
        self.theorem=exact_modified_dispatch_theorem(self)
        self.hashes={**self.local.hashes,**self.theorem['input_hashes'],LOCAL_NAME:sha(LOCAL_NAME),LOCAL_RECEIPT:sha(LOCAL_RECEIPT)}
        for stem in ('current_modified_tensor_dispatch_operator','current_modified_tensor_dispatch'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(local_modified_physical_definition=self.local.modified_physical_tensor_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_tensor_dispatch_operator.py'),
          dispatcher_sha256=sha(PREFIX+'current_modified_tensor_dispatch.py'))
        self.modified_dispatch_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Modified native dispatch receipt exceeds callable-source scope')
            if receipt['modified_dispatch_definition_sha256']!=self.modified_dispatch_definition_sha256:
                raise ValueError('Foreign modified source or continuation operator')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual modified dispatch manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def _arguments(self,Z,coordinate,log_tau,theta,viscosity):
        self.registry.assert_graph()
        if self.pre is not self.local.source.source.histories.pre or endpoints(self.Tw)!=endpoints(self.ctx.mpf(self.pre.params.Tw)):
            raise ValueError('Current source pre/Tw changed after continuation binding')
        c=self.ctx;z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(x) for x in endpoints(z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(z)[0]<-1 or endpoints(z)[1]>1 or endpoints(nu)[0]<=0:
            raise ValueError('Finite native coordinates,Z[-1,1],log(tau),nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        return z,v,lt,nu

    @source_precision
    def continuation(self,Z,coordinate,coordinate_kind='log_radius_offset',log_tau='-1',theta=None,viscosity='1'):
        c=self.ctx;Z,v,lt,nu=self._arguments(Z,coordinate,log_tau,theta,viscosity);lo,hi=endpoints(v)
        if coordinate_kind=='log_radius_offset':
            if lo<2 or hi>endpoints(self.Tw)[0]:raise ValueError('Independent q box requires2<=q<=guaranteed lower actual Tw; use original phase for true endpoint')
        elif coordinate_kind=='original_power_phase':
            if lo<0 or hi>1 or endpoints(self.Tw*v)[0]<=mp.mpf('1.9'):
                raise ValueError('Original phase[0,1] must be certified beyond all fixed repair supports')
        else:raise ValueError('Explicit original power coordinate kind required')
        signed=original_signed_power(self,Z,v,coordinate_kind)
        packet,radius=modified_physical_packet(self,signed)
        point=self.local.lift(c,packet,self.delta,None,lt,theta,nu)
        return dict(point,region='O3_power',coverage_coordinate=v,coverage_coordinate_kind=coordinate_kind,
          actual_original_continuation_signed_source=signed,current_actual_continuation_physical_packet=packet,
          exact_original_geometry_radius_source=radius,
          same_implicit_repair_source_defects_vanish_after_all_supports=True,
          actual_q_equals_same_original_Tw_times_phase=coordinate_kind=='original_power_phase',
          original_incoming_histories_and_same_axis_absolute_pressure_not_reset=True,
          source_enclosures_not_resolved_physical_point_values=True,
          **{key:self.acceptance_loaded for key in GATES},**dict.fromkeys(OPEN,False))

    def _piece(self,region,coordinate,kind,branch,view):
        return dict(source_coordinate_box=self.ctx.mpf(coordinate),source_coordinate_kind=kind,
          source_branch=branch,complete_tensor_source=self.registry._view(region,view,'cylindrical'))

    def _result(self,region,coordinate,kind,pieces):
        return dict(region=region,requested_coordinate=self.ctx.mpf(coordinate),requested_coordinate_kind=kind,
          source_pieces=pieces,pieces_are_alternative_source_charts_not_summed_fields=True,
          modified_dispatch_definition_sha256=self.modified_dispatch_definition_sha256,
          parent_local_physical_definition_sha256=self.local.modified_physical_tensor_definition_sha256,
          actual_original_33_native_region_inventory=list(self.registry.routes),
          modified_affected_native_regions=['O2_buffer','O3_slope_mu','O3_power'],
          full_native_dispatch_not_yet_global_modified_interfaces_or_admissibility=True,
          **{key:self.acceptance_loaded for key in LOCAL_GATES+GATES},**dict.fromkeys(OPEN,False))

    @source_precision
    def power_offset(self,Z,q,log_tau='-1',theta=None,viscosity='1'):
        c=self.ctx;Z,v,_,_=self._arguments(Z,q,log_tau,theta,viscosity);lo,hi=endpoints(v)
        if lo<0 or hi>endpoints(self.Tw)[0]:raise ValueError('Independent q must lie in[0,guaranteed lower actual Tw]; normalized phase covers actual endpoint')
        pieces=[]
        if lo<=2:
            box=c.mpf([lo,min(hi,mp.mpf(2))]);view=self.local.tensor('quiet_O3_power',Z,box,log_tau,theta,viscosity)
            pieces.append(self._piece('O3_power',box,'log_radius_offset','actual_modified_local',view))
        if hi>2:
            box=c.mpf([max(lo,mp.mpf(2)),hi]);view=self.continuation(Z,box,'log_radius_offset',log_tau,theta,viscosity)
            pieces.append(self._piece('O3_power',box,'log_radius_offset','exact_original_after_repair',view))
        return self._result('O3_power',v,'log_radius_offset',pieces)

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        c=self.ctx;Z,v,_,_=self._arguments(Z,coordinate,log_tau,theta,viscosity)
        if region not in self.registry.routes:raise ValueError('One of the original33 native regions required')
        if region in ('O2_buffer','O3_slope_mu'):
            local=self.local.tensor(region,Z,v,log_tau,theta,viscosity)
            return self._result(region,v,'original_native_coordinate',[self._piece(region,v,'original_native_coordinate','actual_modified_local',local)])
        if region=='O3_power':
            lo,hi=endpoints(v)
            if lo<0 or hi>1:raise ValueError('Original O3 power phase[0,1] required')
            pieces=[];cut=self.phase_cut
            if lo<=cut:
                phase=c.mpf([lo,min(hi,cut)]);q=self.Tw*phase
                view=self.local.tensor('quiet_O3_power',Z,q,log_tau,theta,viscosity)
                piece=self._piece(region,phase,'original_power_phase','actual_modified_local',view)
                piece['same_actual_q_offset_enclosure']=q;pieces.append(piece)
            if hi>=cut:
                phase=c.mpf([max(lo,cut),hi]);view=self.continuation(Z,phase,'original_power_phase',log_tau,theta,viscosity)
                pieces.append(self._piece(region,phase,'original_power_phase','exact_original_after_repair',view))
            return self._result(region,v,'original_power_phase',pieces)
        original=self.registry.native(region,Z,v,log_tau,theta,viscosity)
        piece=dict(source_coordinate_box=v,source_coordinate_kind='original_native_coordinate',
          source_branch='unchanged_same_checked_original_route',complete_tensor_source=original)
        return self._result(region,v,'original_native_coordinate',[piece])

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.registry.family,implicit_source_sha256=self.registry.source,
          datum_enclosure_sha256=self.registry.datum_sha,finite_integer_N=self.N,
          modified_dispatch_definition_sha256=self.modified_dispatch_definition_sha256,
          parent_local_modified_physical_definition_sha256=self.local.modified_physical_tensor_definition_sha256,
          exact_actual_modified_dispatch_source_theorem=self.theorem,
          actual_Tw_enclosure=self.Tw,original_phase_cut_for_source_coverage=self.phase_cut,
          actual_modified_33_chart_native_dispatch=list(self.registry.routes),
          affected_native_routes=['O2_buffer','O3_slope_mu','O3_power'],
          complete_modified_dispatch_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Callable full33 native physical tensor dispatch on same checked source graph; actual modified O2 buffer/O3 slope/quiet repair, same implicit-source exact original power continuation through true phase1 endpoint, directed q2/phase chart partitions. All original histories/pressure/radial terms retained, ordinary derivative convention unchanged. Source/physical/heat interface orders, own kinetic energy, common cone N, modified cones/global admissibility, recursion/waves/resolved full corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(LOCAL_GATES+GATES+OPEN,False))

QUERIES=(('whole_modified_native_power','native','O3_power',(0,1)),
  ('whole_modified_buffer','native','O2_buffer',(0,11)),
  ('whole_modified_transition','native','O3_slope_mu',(0,1)),
  ('exact_original_power_endpoint','native','O3_power',1),
  ('straddling_q2_with_nonzero_repair','offset','O3_power',('1.8','2.1')),
  ('fresh_open_power_offset','offset','O3_power','2.337'),
  ('exact_q2','offset','O3_power',2))

@source_precision
def query(field,kind,region,coordinate):
    return field.native(region,(-1,1),coordinate,('-3','-1'),None,'1') if kind=='native' else field.power_offset((-1,1),coordinate,('-3','-1'),None,'1')

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedTensorDispatch(require_checked=False)
    result=field.manifest();views={key:query(field,kind,region,value) for key,kind,region,value in QUERIES}
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual modified33 native tensor/source power continuation constructed; interfaces/common N/cones/full NS open',flush=True)
    return result

if __name__=='__main__':run()
