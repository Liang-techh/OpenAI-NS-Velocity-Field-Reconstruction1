"""Modified compact-repair source carried through original analytic heat tail."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces import (
    CurrentModifiedTensorInterfaces,HERE,PREFIX,NAME as INTERFACES_NAME,RECEIPT as INTERFACES_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,LOCAL_GATES,DISPATCH_GATES,
    GATES as INTERFACE_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_modified_heat_inheritance_operator import (
    DOWNSTREAM,DOWNSTREAM_SEAMS,assert_datum_function_link,exact_modified_heat_inheritance_theorem)

NAME=PREFIX+'current_modified_heat_inheritance.json'
RECEIPT=PREFIX+'current_modified_heat_inheritance_check.json'
VIEWS_NAME=PREFIX+'current_modified_heat_inheritance_views.json.gz'
GATES=('current_modified_downstream_history_pressure_same_function_inherited',
       'current_modified_analytic_preheat_datum_preserved_after_repair',
       'current_modified_analytic_preheat_exterior_compatible_after_repair',
       'current_modified_full_unbounded_Gamma_exterior_source_identity_inherited',
       'current_modified_downstream_heat_source_dispatch_available')
OPEN=tuple(key for key in EARLIER_OPEN if key not in GATES)
ANCESTOR_GATES=LOCAL_GATES+DISPATCH_GATES+INTERFACE_GATES

class CurrentModifiedHeatInheritance:
    @source_precision
    def __init__(self,interfaces=None,require_checked=True):
        self.interfaces=interfaces if interfaces is not None else CurrentModifiedTensorInterfaces()
        if type(self.interfaces) is not CurrentModifiedTensorInterfaces or not self.interfaces.acceptance_loaded:
            raise ValueError('Checked affected modified source/tensor interfaces required')
        self.dispatch=self.interfaces.dispatch;self.registry=self.dispatch.registry;self.ctx=self.dispatch.ctx
        self.N=self.dispatch.N;self.mu=self.dispatch.mu;self.delta=self.dispatch.delta
        self.heat_owner=self.registry.owners['heat'];self.datum=self.dispatch.pre.datum
        self.heat_datum=self.heat_owner.history.flatten.inlet.datum
        self.datum_source_sha=self.datum.source_sha;self.datum_enclosure_sha=self.datum.datum_sha
        self.raw_pressure_witness=self.heat_owner.exterior.pressure.original_pressure_function
        self.controls=self.dispatch.local.source.source.repair.controls
        self.theorem=exact_modified_heat_inheritance_theorem(self)
        self.hashes={**self.interfaces.hashes,**self.theorem['input_hashes'],
          INTERFACES_NAME:sha(INTERFACES_NAME),INTERFACES_RECEIPT:sha(INTERFACES_RECEIPT)}
        for stem in ('current_modified_heat_inheritance_operator','current_modified_heat_inheritance'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_interfaces_definition=self.interfaces.modified_interfaces_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_heat_inheritance_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_heat_inheritance.py'))
        self.modified_heat_inheritance_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Modified heat source inheritance receipt exceeds regional scope')
            if receipt['modified_heat_inheritance_definition_sha256']!=self.modified_heat_inheritance_definition_sha256:
                raise ValueError('Foreign compact repair, datum or heat source graph')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual inherited heat source manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.registry.assert_graph();self.heat_owner.assert_graph();self.heat_owner.history.assert_graph()
        assert_datum_function_link(self)
        if self.registry.owners['heat'] is not self.heat_owner or self.dispatch.pre.datum is not self.datum:
            raise ValueError('Original pressure/heat provider replaced')
        if self.dispatch.local.source.source.repair.controls is not self.controls:
            raise ValueError('Actual implicit controls changed after inheritance binding')
        if self.heat_owner.exterior.pressure.original_pressure_function is not self.raw_pressure_witness:
            raise ValueError('Original analytic raw preheat integral witness changed')

    def _result(self,kind,view):
        return dict(inheritance_kind=kind,complete_same_source_view=view,
          modified_heat_inheritance_definition_sha256=self.modified_heat_inheritance_definition_sha256,
          actual_implicit_source_closure_and_whole_Z_datum_function_link_consumed=True,
          original_five_histories_P0_and_physical_source_units_retained=True,
          regional_heat_identity_not_global_modified_NS_or_energy=True,
          **{key:self.acceptance_loaded for key in ANCESTOR_GATES+GATES},**dict.fromkeys(OPEN,False))

    @source_precision
    def downstream(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph()
        if region not in DOWNSTREAM:raise ValueError('One of fifteen original downstream source regions required')
        view=self.dispatch.native(region,Z,coordinate,log_tau,theta,viscosity)
        return self._result(region,view)

    @source_precision
    def trace(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph()
        if name not in DOWNSTREAM_SEAMS:raise ValueError('One of fifteen original downstream source interfaces required')
        self.dispatch._arguments(Z,0,log_tau,theta,viscosity)
        view=self.registry.trace(name,Z,log_tau,theta,viscosity)
        return self._result(name,view)

    @source_precision
    def unbounded_exterior(self,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();self.dispatch._arguments(Z,0,log_tau,theta,viscosity)
        view=self.registry.exterior(Z,log_tau,theta,viscosity)
        return self._result('whole_unbounded_Gamma_exterior',view)

    @source_precision
    def pressure_datum(self,Z=(-1,1)):
        self.assert_graph();z,_,_,_=self.dispatch._arguments(Z,0,-1,None,1)
        view=self.datum.normalized_jets(endpoints(z),5)
        return self._result('original_analytic_preheat_P0_axial5_enclosure',view)

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.registry.family,implicit_source_sha256=self.registry.source,
          datum_enclosure_sha256=self.registry.datum_sha,finite_integer_N=self.N,
          modified_heat_inheritance_definition_sha256=self.modified_heat_inheritance_definition_sha256,
          parent_modified_interfaces_definition_sha256=self.interfaces.modified_interfaces_definition_sha256,
          original_analytic_pressure_source_sha256=self.datum_source_sha,
          exact_actual_modified_downstream_heat_inheritance_theorem=self.theorem,
          inherited_downstream_regions=list(DOWNSTREAM),inherited_downstream_interfaces=list(DOWNSTREAM_SEAMS),
          complete_inherited_downstream_heat_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Same exact implicit compact repair closure and original whole-Z power phase1 inlet compose with15 unchanged downstream native providers and15 original source tensor interfaces. Actual analytic P0/native14-stage raw preheat pressure function, original waiting root, full five terminal histories and complete unbounded Gamma heat exterior are inherited in callable pressure/downstream/trace/exterior APIs. No fitted pressure, reset histories or radial truncation. Regional exact heat tensor/remainder/momentum identity is inherited; higher modified physical interfaces, actual total kinetic energy, common cone N/modified/global cones, recursion/waves/resolved global corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(ANCESTOR_GATES+GATES+OPEN,False))

QUERIES=(('entrance','pulse_entrance',(0,'.02')),('main','pulse_main',('.02',10)),
  ('exit','pulse_exit',(10,11)),('gap','pulse_gap',(11,12)),('gap_end','pulse_gap_end',(0,1)),
  ('end','pulse_end',(-4,0)),('flatten','flatten',(0,100)),('outer_power','outer_power',(0,1)),
  ('outer_angular','outer_angular',(-4,0)),('steep_entry','steep_entry',(0,1)),
  ('steep_power','steep_power',(0,1)),('steep_exit','steep_exit',(0,1)),('waiting','waiting',(0,1)),
  ('collar','heat_collar',(0,3)),('exterior_far','heat_exterior',(10,1000)))
TRACE_QUERIES=('incoming_entrance','waiting_collar','collar_exterior')

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedHeatInheritance(require_checked=False)
    result=field.manifest();views={}
    for key,region,value in QUERIES:
        views[key]=field.downstream(region,(-1,1),value,('-3','-1'),None,'1')
        print('Inherited same-source downstream: '+region,flush=True)
    views['pressure_datum']=field.pressure_datum((-1,1))
    for name in TRACE_QUERIES:views['trace_'+name]=field.trace(name)
    views['unbounded_Gamma_exterior']=field.unbounded_exterior()
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
