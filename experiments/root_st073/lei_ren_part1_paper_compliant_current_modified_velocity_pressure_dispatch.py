"""Same-source full33 native velocity/absolute-pressure factored bounds."""
import gzip
import hashlib
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces import (
    CurrentModifiedVelocityPressureInterfaces,HERE,PREFIX,NAME as PARENT_NAME,RECEIPT as PARENT_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,INHERITED_GATES as EARLIER_GATES,
    GATES as INTERFACE_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch_operator import (
    canonical_original_velocity,original_power_velocity,gap_source_coordinates,gap_main_cover_velocity,
    exact_modified_velocity_dispatch_theorem)

NAME=PREFIX+'current_modified_velocity_pressure_dispatch.json'
RECEIPT=PREFIX+'current_modified_velocity_pressure_dispatch_check.json'
VIEWS_NAME=PREFIX+'current_modified_velocity_pressure_dispatch_views.json.gz'
GATES=('current_modified_full_33_chart_velocity_pressure_dispatch_installed',
  'current_modified_power_q_phase_velocity_continuation_available',
  'current_modified_original_core_axis_velocity_pressure_route_available',
  'current_modified_full_unbounded_heat_velocity_pressure_bounds_available',
  'current_modified_native_velocity_pressure_coordinate_bridges_certified')
OPEN=tuple(key for key in EARLIER_OPEN if key not in GATES)+(
  'current_modified_all_global_velocity_pressure_interfaces_certified',
  'current_modified_velocity_pressure_physical_locator_integrated')
INHERITED_GATES=EARLIER_GATES+INTERFACE_GATES

class CurrentModifiedVelocityPressureDispatch:
    @source_precision
    def __init__(self,interfaces=None,require_checked=True):
        self.interfaces=interfaces if interfaces is not None else CurrentModifiedVelocityPressureInterfaces()
        if type(self.interfaces) is not CurrentModifiedVelocityPressureInterfaces or not self.interfaces.acceptance_loaded:
            raise ValueError('Checked actual twelve affected velocity/pressure interfaces required')
        self.velocity=self.interfaces.velocity;self.geometry=self.velocity.geometry
        self.registry=self.velocity.heat.registry;self.pre=self.velocity.heat.dispatch.pre
        self.ctx=self.velocity.ctx;self.delta=self.velocity.delta;self.N=self.velocity.N;self.mu=self.velocity.mu
        self.controls=self.velocity.source.repair.controls
        self.Tw=self.velocity.heat.dispatch.Tw;self.phase_cut=self.velocity.heat.dispatch.phase_cut
        self.gap_phase_cut=self.ctx.mpf(1)/20000
        self.theorem=exact_modified_velocity_dispatch_theorem(self)
        self.hashes={**self.interfaces.hashes,**self.theorem['input_hashes'],PARENT_NAME:sha(PARENT_NAME),PARENT_RECEIPT:sha(PARENT_RECEIPT)}
        for stem in ('current_modified_velocity_pressure_dispatch_operator','current_modified_velocity_pressure_dispatch'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_interfaces_definition=self.interfaces.modified_velocity_pressure_interfaces_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_velocity_pressure_dispatch_operator.py'),
          dispatcher_sha256=sha(PREFIX+'current_modified_velocity_pressure_dispatch.py'))
        self.modified_velocity_pressure_dispatch_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Full33 native velocity/pressure receipt exceeds source-query scope')
            if receipt['modified_velocity_pressure_dispatch_definition_sha256']!=self.modified_velocity_pressure_dispatch_definition_sha256:
                raise ValueError('Foreign current full33 velocity/pressure dispatcher')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual full33 velocity/pressure manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def _arguments(self,Z,coordinate,log_tau,theta):
        self.velocity.heat.assert_graph();self.geometry.assert_graph()
        if self.geometry is not self.registry.physical or self.pre is not self.geometry.pre or self.controls is not self.velocity.source.repair.controls:
            raise ValueError('Same actual canonical geometry, pre and implicit controls required')
        Z,v,lt,_=self.velocity.heat.dispatch._arguments(Z,coordinate,log_tau,theta,1)
        return Z,v,lt

    def _piece(self,coordinate,kind,branch,view):
        return dict(source_coordinate_box=self.ctx.mpf(coordinate),source_coordinate_kind=kind,
          source_branch=branch,complete_velocity_pressure_source=view)

    def _result(self,region,coordinate,kind,pieces):
        return dict(region=region,requested_coordinate=self.ctx.mpf(coordinate),requested_coordinate_kind=kind,
          source_pieces=pieces,pieces_are_alternative_source_charts_not_summed_fields=True,
          modified_velocity_pressure_dispatch_definition_sha256=self.modified_velocity_pressure_dispatch_definition_sha256,
          same_unique_implicit_repair_definition_sha256=self.velocity.source.repair.repair_definition_sha256,
          actual_original_33_native_region_inventory=list(self.registry.routes),
          modified_affected_native_regions=['O2_buffer','O3_slope_mu','O3_power'],
          resolved_physical_point_values_available=False,
          **{key:self.acceptance_loaded for key in INHERITED_GATES+GATES},**dict.fromkeys(OPEN,False))

    @source_precision
    def power_offset(self,Z,q,log_tau='-1',theta=None):
        c=self.ctx;Z,v,lt=self._arguments(Z,q,log_tau,theta);lo,hi=endpoints(v)
        if lo<0 or hi>endpoints(self.Tw)[0]:raise ValueError('Independent q requires[0,guaranteed lower Tw]; phase API covers true endpoint')
        pieces=[]
        if lo<=2:
            box=c.mpf([lo,min(hi,mp.mpf(2))]);view=self.velocity.physical('quiet_O3_power',Z,box,lt,theta)
            pieces.append(self._piece(box,'log_radius_offset','actual_modified_local',dict(view,physical_layout='four_label_Cartesian')))
        if hi>2:
            box=c.mpf([max(lo,mp.mpf(2)),hi]);view=original_power_velocity(self,Z,box,'log_radius_offset',lt,theta)
            pieces.append(self._piece(box,'log_radius_offset','exact_original_after_repair',view))
        return self._result('O3_power',v,'log_radius_offset',pieces)

    @source_precision
    def power_phase(self,Z,phase,log_tau='-1',theta=None):
        c=self.ctx;Z,v,lt=self._arguments(Z,phase,log_tau,theta);lo,hi=endpoints(v)
        if lo<0 or hi>1:raise ValueError('Original incoming power phase[0,1] required')
        pieces=[];cut=self.phase_cut
        if lo<=cut:
            box=c.mpf([lo,min(hi,cut)]);q=self.Tw*box
            view=self.velocity.physical('quiet_O3_power',Z,q,lt,theta)
            piece=self._piece(box,'original_power_phase','actual_modified_local',dict(view,physical_layout='four_label_Cartesian'))
            piece['same_actual_q_offset_enclosure']=q;pieces.append(piece)
        if hi>=cut:
            box=c.mpf([max(lo,cut),hi]);view=original_power_velocity(self,Z,box,'original_power_phase',lt,theta)
            pieces.append(self._piece(box,'original_power_phase','exact_original_after_repair',view))
        return self._result('O3_power',v,'original_power_phase',pieces)

    @source_precision
    def gap_phase(self,Z,phase,log_tau='-1',theta=None):
        c=self.ctx;Z,v,lt=self._arguments(Z,phase,log_tau,theta);lo,hi=endpoints(v)
        if lo<0 or hi>1:raise ValueError('Original reciprocal gap coverage phase[0,1] required')
        pieces=[];cl,ch=endpoints(self.gap_phase_cut)
        if lo<=ch:
            box=c.mpf([lo,min(hi,ch)]);view=gap_main_cover_velocity(self,Z,box,lt,theta)
            pieces.append(self._piece(box,'reciprocal_gap_phase','same_current_supplemental_gap_function',view))
        if hi>=cl:
            box=c.mpf([max(lo,cl),hi]);_,ss=gap_source_coordinates(self.geometry.pulse.mu,box)
            view=canonical_original_velocity(self,'pulse_gap_end',Z,ss,lt,theta)
            piece=self._piece(box,'reciprocal_gap_phase','same_current_gap_end_function',view)
            piece.update(exact_gap_end_source_function='s=-4-(1/mu-4)*(1-phase)',actual_same_source_s_enclosure=ss)
            pieces.append(piece)
        return self._result('pulse_gap_end',v,'reciprocal_gap_phase',pieces)

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None):
        Z,v,lt=self._arguments(Z,coordinate,log_tau,theta)
        if region not in self.registry.routes:raise ValueError('One of original33 native regions required')
        if region=='O3_power':return self.power_phase(Z,v,lt,theta)
        if region=='pulse_gap_end':return self.gap_phase(Z,v,lt,theta)
        if region in ('O2_buffer','O3_slope_mu'):
            view=dict(self.velocity.physical(region,Z,v,lt,theta),physical_layout='four_label_Cartesian')
            branch='actual_modified_local'
        else:
            view=canonical_original_velocity(self,region,Z,v,lt,theta,axis=region=='core_positive_radius' and endpoints(v)==(0,0))
            branch='unchanged_same_checked_canonical_route'
        return self._result(region,v,'original_native_coordinate',[self._piece(v,'original_native_coordinate',branch,view)])

    @source_precision
    def axis(self,Z,log_tau='-1',theta=None):
        return self.native('core_positive_radius',Z,0,log_tau,theta)

    @source_precision
    def unbounded_exterior(self,Z=(-1,1),log_tau=('-3','-1'),theta=None):
        Z,_,lt=self._arguments(Z,0,log_tau,theta)
        box=self.ctx.mpf([3,mp.inf]);view=canonical_original_velocity(self,'heat_exterior',Z,box,lt,theta)
        view.update(unbounded_coordinate_is_whole_domain_box_not_point_at_infinity=True,
          same_checked_full_Gamma_exterior_source_and_absolute_pressure_retained=True)
        return self._result('heat_exterior',box,'unbounded_heat_offset_box',[
          self._piece(box,'unbounded_heat_offset_box','unchanged_same_checked_full_Gamma_source',view)])

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.registry.family,implicit_source_sha256=self.registry.source,
          datum_enclosure_sha256=self.registry.datum_sha,finite_integer_N=self.N,
          modified_velocity_pressure_dispatch_definition_sha256=self.modified_velocity_pressure_dispatch_definition_sha256,
          parent_modified_velocity_pressure_interfaces_definition_sha256=self.interfaces.modified_velocity_pressure_interfaces_definition_sha256,
          exact_actual_modified_velocity_pressure_dispatch_theorem=self.theorem,
          actual_Tw_enclosure=self.Tw,original_power_phase_cut=self.phase_cut,reciprocal_gap_phase_cut=self.gap_phase_cut,
          actual_modified_33_chart_native_velocity_pressure_dispatch=list(self.registry.routes),
          affected_native_routes=['O2_buffer','O3_slope_mu','O3_power'],complete_velocity_pressure_dispatch_views=VIEWS_NAME,
          original_native_regions=33,original_adjacent_interfaces=32,original_internal_tensor_traces=14,
          spatial_cartesian_multiindex_count=35,
          scope='Full33 native velocity/absolute-pressure spatial4/fixed-x time1 function enclosures on the same actual checked canonical graph. Own O2/O3/compact repair and true original power continuation; explicit entrance xi/mu bridge, directed reciprocal gap phase cover, native nonsingular core/axis and full unbounded Gamma heat box. Twelve affected velocity/pressure interfaces inherited. Global velocity interfaces, physical-coordinate locator integration, total physical energy, common N/modified cones, true recursion/waves/resolved corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(INHERITED_GATES+GATES+OPEN,False))

def native_domains(field):
    c=field.ctx
    domains=dict(core_positive_radius=(0,4),bridge_first=(0,1),bridge_second=(1,2),bridge_macro=(0,1),
      switch_first=(0,1),switch_second=(1,2),switch_power=(0,1),reshape=(0,1),inner_reference=(0,1),
      axial_restore=(0,1),restore_buffer=(-7,-6),actual_patch=(1,endpoints(c.exp(1))[1]),Rh_reference=(-5,0),
      O2_slope=(0,1),O2_axial=(0,1),O2_buffer=(0,11),O3_slope_mu=(0,1),O3_power=(0,1),
      pulse_entrance=(0,'.02'),pulse_main=('.02',10),pulse_exit=(10,11),pulse_gap=(11,12),pulse_gap_end=(0,1),
      pulse_end=(-4,0),flatten=(0,100),outer_power=(0,1),outer_angular=(-4,0),steep_entry=(0,1),
      steep_power=(0,1),steep_exit=(0,1),waiting=(0,1),heat_collar=(0,3),heat_exterior=(3,7))
    if set(domains)!=set(field.registry.routes):raise ValueError('Exact original33 domain inventory required')
    return domains

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressureDispatch(require_checked=False)
    result=field.manifest();views={}
    for region,coordinate in native_domains(field).items():
        views[region]=field.native(region,(-1,1),coordinate,('-3','-1'),None)
        print('Current full33 native velocity/absolute pressure: '+region,flush=True)
    views['axis']=field.axis((-1,1),('-3','-1'))
    views['unbounded_exterior']=field.unbounded_exterior()
    views['straddling_q2']=field.power_offset((-1,1),('1.8','2.1'),('-3','-1'))
    views['exact_power_endpoint']=field.power_phase((-1,1),1,('-3','-1'))
    views['exact_gap_start']=field.gap_phase((-1,1),0,('-3','-1'))
    views['exact_gap_end']=field.gap_phase((-1,1),1,('-3','-1'))
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
