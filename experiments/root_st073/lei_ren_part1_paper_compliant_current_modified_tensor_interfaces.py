"""Callable same-source full tensor traces at twelve affected interfaces."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_modified_tensor_dispatch import (
    CurrentModifiedTensorDispatch,HERE,PREFIX,NAME as DISPATCH_NAME,RECEIPT as DISPATCH_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,LOCAL_GATES, GATES as DISPATCH_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import canonical_tensor_groups
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces_operator import (
    REPAIR_EDGES,compile_named_endpoint_germs,physical_history_germ,exact_modified_interfaces_theorem)

NAME=PREFIX+'current_modified_tensor_interfaces.json'
RECEIPT=PREFIX+'current_modified_tensor_interfaces_check.json'
VIEWS_NAME=PREFIX+'current_modified_tensor_interfaces_views.json.gz'
GATES=('current_modified_O2_buffer_inlet_source_tensor_join_certified',
       'current_modified_O2_O3_shared_source_tensor_join_certified',
       'current_modified_O3_transition_power_source_tensor_join_certified',
       'current_modified_modulation_flat_support_tensor_joins_certified',
       'current_modified_six_repair_flat_support_tensor_joins_certified',
       'current_modified_twelve_affected_tensor_traces_available')
OPEN=EARLIER_OPEN+('current_modified_full_spatial4_temporal1_interfaces_certified',)
SEAMS=('buffer_inlet','buffer_transition','transition_power','buffer_modulation_start',
       'transition_modulation_end')+tuple(REPAIR_EDGES)+('repair_exit',)

class CurrentModifiedTensorInterfaces:
    @source_precision
    def __init__(self,dispatch=None,require_checked=True):
        self.dispatch=dispatch if dispatch is not None else CurrentModifiedTensorDispatch()
        if type(self.dispatch) is not CurrentModifiedTensorDispatch or not self.dispatch.acceptance_loaded:
            raise ValueError('Checked complete modified tensor dispatcher required')
        self.ctx=self.dispatch.ctx;self.delta=self.dispatch.delta;self.N=self.dispatch.N;self.mu=self.dispatch.mu
        self.repaired=self.dispatch.local.source.source;self.controls=self.repaired.repair.controls
        self.flat_modulation,self.flat_repair,self.adaptation=compile_named_endpoint_germs()
        self.theorem=exact_modified_interfaces_theorem(self)
        self.hashes={**self.dispatch.hashes,**self.theorem['input_hashes'],**self.adaptation['input_hashes'],
          DISPATCH_NAME:sha(DISPATCH_NAME),DISPATCH_RECEIPT:sha(DISPATCH_RECEIPT)}
        for stem in ('current_modified_tensor_interfaces_operator','current_modified_tensor_interfaces'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_dispatch_definition=self.dispatch.modified_dispatch_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_tensor_interfaces_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_tensor_interfaces.py'))
        self.modified_interfaces_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Affected tensor interface receipt exceeds retained-order scope')
            if receipt['modified_interfaces_definition_sha256']!=self.modified_interfaces_definition_sha256:
                raise ValueError('Foreign affected interface source/adaptation')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual interface manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def endpoint_source(self,name,Z=(-1,1),side='left'):
        """Named exact germs; general crossing intervals stay with dispatcher."""
        if side not in ('left','right'):raise ValueError('Named left or right endpoint germ required')
        c=self.ctx;z,_,_,_=self.dispatch._arguments(Z,0,-1,None,1)
        if name in REPAIR_EDGES:
            a,b,_,_=REPAIR_EDGES[name];q=c.mpf(a)/b
            source=self.flat_repair(self.repaired,'quiet_O3_power',z,q,name)
            return dict(source,named_source_edge=name,named_source_side=side,
              exact_edge_rational=(a,b),both_germs_use_same_continuous_partial_integrals=True)
        if name=='transition_modulation_end':
            q=c.mpf(1)/2
            source=(self.repaired.history('O3_slope_mu',z,q) if side=='left' else
              self.flat_modulation(self.repaired,'O3_slope_mu',z,q))
            return dict(source,named_source_edge=name,named_source_side=side,
              flat_local_profiles_do_not_zero_own_cumulative_histories=True)
        raise ValueError('One of the six repair edges or modulation end required')

    def _original(self,region,Z,coordinate,lt,theta,nu):
        return self.dispatch.registry.native(region,Z,coordinate,lt,theta,nu)['original_complete_view']

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name not in SEAMS:raise ValueError('One of the twelve named affected interfaces required')
        c=self.ctx;Z,_,lt,nu=self.dispatch._arguments(Z,0,log_tau,theta,viscosity)
        local=self.dispatch.local
        if name=='buffer_inlet':
            left=self._original('O2_axial',Z,1,lt,theta,nu)
            right=local.tensor('O2_buffer',Z,0,lt,theta,nu)
        elif name=='buffer_transition':
            left=local.tensor('O2_buffer',Z,11,lt,theta,nu)
            right=local.tensor('O3_slope_mu',Z,0,lt,theta,nu)
        elif name=='transition_power':
            left=local.tensor('O3_slope_mu',Z,1,lt,theta,nu)
            right=local.tensor('quiet_O3_power',Z,0,lt,theta,nu)
        elif name=='buffer_modulation_start':
            left=self._original('O2_buffer',Z,9,lt,theta,nu)
            right=local.tensor('O2_buffer',Z,9,lt,theta,nu)
        elif name in REPAIR_EDGES or name=='transition_modulation_end':
            left=physical_history_germ(self,self.endpoint_source(name,Z,'left'),lt,theta,nu)
            right=physical_history_germ(self,self.endpoint_source(name,Z,'right'),lt,theta,nu)
        else:
            left=local.tensor('quiet_O3_power',Z,2,lt,theta,nu)
            right=self.dispatch.continuation(Z,2,'log_radius_offset',lt,theta,nu)
        a=canonical_tensor_groups(left);b=canonical_tensor_groups(right);common={}
        if set(a)!=set(b):raise ValueError('Complete tensor trace layout changed')
        for key in a:
            values=[endpoints(part['log_absolute_upper'])[1] for part in a[key]+b[key] if not part['exact_zero']]
            common[key]=dict(exact_zero=not values,
              log_absolute_upper=c.mpf(max(values))+c.ln(len(values)) if values else None)
        return dict(seam=name,left_complete_tensor_source=left,right_complete_tensor_source=right,
          common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
          modified_interfaces_definition_sha256=self.modified_interfaces_definition_sha256,
          source_function_identity_precedes_triangle_enclosures=True,
          shared_original_axis_P0_and_all_own_five_histories_retained=True,
          interval_overlap_not_used_as_source_identity=True,
          source_orders=dict(ordinary_logR=4,axial_Z=5,stress=3,divergence_and_remainder=2),
          **{key:self.acceptance_loaded for key in LOCAL_GATES+DISPATCH_GATES+GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.dispatch.registry.family,
          implicit_source_sha256=self.dispatch.registry.source,datum_enclosure_sha256=self.dispatch.registry.datum_sha,
          finite_integer_N=self.N,modified_interfaces_definition_sha256=self.modified_interfaces_definition_sha256,
          parent_modified_dispatch_definition_sha256=self.dispatch.modified_dispatch_definition_sha256,
          exact_actual_affected_interface_source_theorem=self.theorem,
          named_endpoint_source_adaptations=self.adaptation,affected_named_interfaces=list(SEAMS),
          complete_affected_interface_views=VIEWS_NAME,
          original_native_regions=33,original_adjacent_interfaces=32,original_internal_tensor_traces=14,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Twelve callable affected source-function tensor traces: three actual O2/O3 native joins, flat modulation start/end, six exact fixed repair edges and q2 exact implicit-source exit. Own logR4/Z5 profiles, five moments, absolute P0+Cp pressure and radial source yield full signed physical stress3/divergence2/remainder2 and Cartesian traces with71 common groups. Named source edge correlations are reduced before enclosure; intervals are never reset. Higher spatial4/time1/global physical and heat interfaces, own total energy, common cone N/modified cones, recursion/waves and full corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(LOCAL_GATES+DISPATCH_GATES+GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedTensorInterfaces(require_checked=False)
    result=field.manifest();views={}
    for name in SEAMS:
        views[name]=field.interface(name)
        print('Affected same-source tensor interface: '+name,flush=True)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
