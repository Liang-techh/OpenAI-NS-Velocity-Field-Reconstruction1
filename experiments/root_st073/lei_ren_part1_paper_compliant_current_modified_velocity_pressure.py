"""Callable own modified velocity and absolute-pressure spatial4/time1 bounds."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_modified_heat_inheritance import (
    CurrentModifiedHeatInheritance,HERE,PREFIX,NAME as HEAT_NAME,RECEIPT as HEAT_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,ANCESTOR_GATES,
    GATES as HEAT_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_operator import (
    modified_velocity_pressure_packet,modified_physical_velocity_pressure,exact_modified_velocity_pressure_theorem)

NAME=PREFIX+'current_modified_velocity_pressure.json'
RECEIPT=PREFIX+'current_modified_velocity_pressure_check.json'
VIEWS_NAME=PREFIX+'current_modified_velocity_pressure_views.json.gz'
GATES=('current_modified_own_velocity_absolute_pressure_mixed4_packet_available',
  'current_modified_absolute_P0_own_Cp_mixed4_source_certified',
  'current_modified_velocity_pressure_actual_signed_units_certified',
  'current_modified_local_velocity_pressure_cartesian_spatial4_time1_available')
OPEN=EARLIER_OPEN
INHERITED_GATES=ANCESTOR_GATES+HEAT_GATES

class CurrentModifiedVelocityPressure:
    @source_precision
    def __init__(self,heat=None,require_checked=True):
        self.heat=heat if heat is not None else CurrentModifiedHeatInheritance()
        if type(self.heat) is not CurrentModifiedHeatInheritance or not self.heat.acceptance_loaded:
            raise ValueError('Checked actual analytic pressure/repair lineage required')
        self.source=self.heat.interfaces.repaired;self.geometry=self.heat.dispatch.local.geometry
        self.ctx=self.heat.ctx;self.delta=self.heat.delta;self.N=self.heat.N;self.mu=self.heat.mu
        self.theorem=exact_modified_velocity_pressure_theorem(self)
        self.hashes={**self.heat.hashes,**self.theorem['input_hashes'],HEAT_NAME:sha(HEAT_NAME),HEAT_RECEIPT:sha(HEAT_RECEIPT)}
        for stem in ('current_modified_velocity_pressure_operator','current_modified_velocity_pressure'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_heat_definition=self.heat.modified_heat_inheritance_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_velocity_pressure_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_velocity_pressure.py'))
        self.modified_velocity_pressure_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Modified velocity/pressure packet receipt exceeds local source scope')
            if receipt['modified_velocity_pressure_definition_sha256']!=self.modified_velocity_pressure_definition_sha256:
                raise ValueError('Foreign own velocity/pressure source or physical operator')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual own velocity/pressure source manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def physical(self,region,Z,coordinate,log_tau='-1',theta=None):
        self.heat.assert_graph()
        Z,coordinate,lt,_=self.heat.dispatch._arguments(Z,coordinate,log_tau,theta,1)
        view=self.source.history(region,Z,coordinate)
        packet=modified_velocity_pressure_packet(self,view)
        physical=modified_physical_velocity_pressure(self,packet,view['Z'],lt,theta)
        public_packet=dict(packet,grids={label:{'y%d_Z%d'%index:terms for index,terms in grid.items()}
          for label,grid in packet['grids'].items()},
          amplitudes={label:{physical['source_log_base_names'][index]:power for index,power in parts.items()}
            for label,parts in packet['amplitudes'].items()})
        return dict(physical,region=region,coordinate=coordinate,Z=Z,
          modified_velocity_pressure_definition_sha256=self.modified_velocity_pressure_definition_sha256,
          parent_modified_source_definition_sha256=self.source.modified_source_definition_sha256,
          same_unique_implicit_repair_definition_sha256=self.source.repair.repair_definition_sha256,
          actual_own_velocity_pressure_mixed4_packet=public_packet,
          actual_own_history_pressure_source=view,
          source_function_identities_precede_physical_enclosures=True,
          complete_physical_velocity_and_absolute_pressure_source_bounds_available=True,
          quantified_affected_interfaces_require_separate_two_germ_construction=True,
          **{k:self.acceptance_loaded for k in INHERITED_GATES+GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.heat.registry.family,
          implicit_source_sha256=self.heat.registry.source,datum_enclosure_sha256=self.heat.registry.datum_sha,
          finite_integer_N=self.N,modified_velocity_pressure_definition_sha256=self.modified_velocity_pressure_definition_sha256,
          parent_modified_heat_inheritance_definition_sha256=self.heat.modified_heat_inheritance_definition_sha256,
          actual_modified_velocity_pressure_source_theorem=self.theorem,
          complete_modified_velocity_pressure_views=VIEWS_NAME,
          scoped_native_regions=('O2_buffer','O3_slope_mu','quiet_O3_power'),
          spatial_cartesian_multiindex_count=35,mixed_source_indices_per_label=15,
          physical_contribution_groups_per_query=216,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Own modified O2 buffer/O3 slope/quiet repaired power UT/UZ/UR/absolute P ordinary mixed4 packet. Actual signed radial Pstar0+1, theta/axial Pstar1 and pressure Pstar2 units, original P0+own Cp FTC, exact original radius, unchanged Cartesian spatial4 and fixed-x time1 source operators retained. Outputs are signed factored enclosures of actual functions. Twelve quantified physical interfaces, continuation/full33 velocity dispatch, total physical energy, common cone N, modified/global cones, true recursion/waves/resolved full corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(INHERITED_GATES+GATES+OPEN,False))

QUERIES=(('whole_buffer','O2_buffer',(0,11)),('whole_transition','O3_slope_mu',(0,1)),
  ('whole_repair','quiet_O3_power',(1,2)),('modulation_interior','O3_slope_mu','.237337'),
  ('repair_interior','quiet_O3_power','1.201337'),('restored_exit','quiet_O3_power',2))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressure(require_checked=False)
    result=field.manifest();views={}
    for name,region,value in QUERIES:
        views[name]=field.physical(region,(-1,1),value,('-3','-1'),None)
        print('Own mixed4 velocity/absolute pressure physical source: '+name,flush=True)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
