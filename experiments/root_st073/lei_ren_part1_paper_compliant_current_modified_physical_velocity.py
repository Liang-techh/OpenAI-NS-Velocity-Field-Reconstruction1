"""Physical Cartesian/log-radius queries of the actual modified velocity/P."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch import (
    CurrentModifiedVelocityPressureDispatch,HERE,PREFIX,NAME as DISPATCH_NAME,RECEIPT as DISPATCH_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,
    INHERITED_GATES as EARLIER_GATES,GATES as DISPATCH_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_global_tensor_cover import (
    CurrentGlobalTensorCover,NAME as COVER_NAME,RECEIPT as COVER_RECEIPT)
from lei_ren_part1_paper_compliant_current_modified_physical_velocity_operator import (
    actual_candidate_velocity,existing_original_cover,exact_modified_physical_velocity_theorem)

NAME=PREFIX+'current_modified_physical_velocity.json'
RECEIPT=PREFIX+'current_modified_physical_velocity_check.json'
VIEWS_NAME=PREFIX+'current_modified_physical_velocity_views.json.gz'
GATES=('current_modified_velocity_pressure_physical_locator_integrated',
  'current_modified_physical_cartesian_velocity_pressure_queries_available',
  'current_modified_velocity_pressure_actual_lambda_bounds_available',
  'current_modified_velocity_pressure_constant_viscosity_bridge_certified')
OPEN=tuple(key for key in EARLIER_OPEN if key not in GATES)
INHERITED_GATES=EARLIER_GATES+DISPATCH_GATES

class CurrentModifiedPhysicalVelocity:
    @source_precision
    def __init__(self,dispatch=None,cover=None,require_checked=True):
        self.dispatch=dispatch if dispatch is not None else CurrentModifiedVelocityPressureDispatch()
        self.cover=cover if cover is not None else existing_original_cover(self.dispatch)
        self.ctx=self.dispatch.ctx;self.registry=self.dispatch.registry
        self.N=self.dispatch.N;self.delta=self.dispatch.delta;self.mu=self.dispatch.mu
        self.theorem=exact_modified_physical_velocity_theorem(self)
        self.hashes={**self.dispatch.hashes,**self.cover.hashes,**self.theorem['input_hashes'],
          DISPATCH_NAME:sha(DISPATCH_NAME),DISPATCH_RECEIPT:sha(DISPATCH_RECEIPT),
          COVER_NAME:sha(COVER_NAME),COVER_RECEIPT:sha(COVER_RECEIPT)}
        for stem in ('current_modified_physical_velocity_operator','current_modified_physical_velocity'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_velocity_dispatch_definition=self.dispatch.modified_velocity_pressure_dispatch_definition_sha256,
          original_global_physical_cover_sha256=sha(COVER_RECEIPT),
          operator_sha256=sha(PREFIX+'current_modified_physical_velocity_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_physical_velocity.py'))
        self.modified_physical_velocity_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Physical velocity query receipt exceeds actual source scope')
            if receipt['modified_physical_velocity_definition_sha256']!=self.modified_physical_velocity_definition_sha256:
                raise ValueError('Foreign physical velocity/pressure function or coordinate source')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual physical velocity manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def _consume(self,covered):
        self.cover.assert_graph();self.dispatch.velocity.heat.assert_graph()
        if self.cover.registry is not self.dispatch.registry or self.cover.locator.physical is not self.dispatch.geometry:
            raise ValueError('One actual physical candidate/source graph required')
        location=covered['original_physical_source_result']
        if not covered['independent_analytic_cover_applies'] or not covered['finite_directed_root_bracket_retained']:
            raise ValueError('Actual admitted all-point directed candidate cover required')
        views=[actual_candidate_velocity(self,candidate,location) for candidate in location['candidates']]
        return dict(location=location,current_candidate_velocity_pressure_views=views,
          actual_all_point_source_cover=covered,
          modified_physical_velocity_definition_sha256=self.modified_physical_velocity_definition_sha256,
          same_actual_velocity_dispatch_definition_sha256=self.dispatch.modified_velocity_pressure_dispatch_definition_sha256,
          candidate_union_and_native_pieces_are_alternative_charts_not_summed_fields=True,
          actual_implicit_lambda_Z_and_source_parameter_uncertainty_retained=True,
          resolved_physical_point_values_available=False,
          **{key:self.acceptance_loaded for key in INHERITED_GATES+GATES},**dict.fromkeys(OPEN,False))

    @source_precision
    def cartesian(self,x,y,z,time,viscosity=1,terminal_time=1):
        covered=self.cover.cartesian(x,y,z,time,viscosity,terminal_time,tensors=False)
        return self._consume(covered)

    @source_precision
    def log_radius(self,log_r,z,log_tau,theta=0,viscosity=1):
        covered=self.cover.locate_log_radius(log_r,z,log_tau,theta,viscosity)
        return self._consume(covered)

    @source_precision
    def correlated(self,request,z,log_tau,theta=0,viscosity=1):
        location=self.cover.field.locate(request,z,log_tau,theta,viscosity)
        covered=self.cover._attach(location,request=request)
        return self._consume(covered)

    def radius(self,anchor,offset=0):return self.cover.field.radius(anchor,offset)

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.registry.family,implicit_source_sha256=self.registry.source,
          datum_enclosure_sha256=self.registry.datum_sha,finite_integer_N=self.N,
          modified_physical_velocity_definition_sha256=self.modified_physical_velocity_definition_sha256,
          parent_modified_velocity_dispatch_definition_sha256=self.dispatch.modified_velocity_pressure_dispatch_definition_sha256,
          exact_actual_modified_physical_velocity_theorem=self.theorem,
          original_native_regions=33,original_adjacent_interfaces=32,original_internal_tensor_traces=14,
          complete_physical_velocity_pressure_query_views=VIEWS_NAME,
          scope='Finite physical Cartesian/log-radius and exact source-correlated velocity/absolute-pressure function queries on the same current33 native graph. Every directed locator candidate and native partition retained, actual lambda refined separately from tau, arbitrary positive constant viscosity source dilation, native core/axis and finite actual Gamma heat coordinate. Outputs are signed factored enclosures, not resolved point coefficients. Global smooth velocity interfaces, physical total energy, common N/modified cones, true coefficient recursion/waves/resolved corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(INHERITED_GATES+GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPhysicalVelocity(require_checked=False)
    result=field.manifest();views={}
    queries=(('physical_axis',lambda:field.cartesian(0,0,'.237','.6',viscosity='.8')),
      ('physical_cartesian',lambda:field.cartesian('1','.337','.237','.6',viscosity='.8')),
      ('source_transition',lambda:field.correlated(field.radius('Rd','.237337'),'.237','-2.337','.337','.8')),
      ('source_repair',lambda:field.correlated(field.radius('Rw','1.201337'),'.237','-2.337','.337','1.3')),
      ('source_heat',lambda:field.correlated(field.radius('Rt','4.337'),'.237','-2.337','.337','.8')))
    for name,query in queries:
        views[name]=query();print('Actual physical velocity/absolute-pressure query: '+name,flush=True)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
