"""Twelve actual affected velocity/absolute-pressure physical interface traces."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure import (
    CurrentModifiedVelocityPressure,HERE,PREFIX,NAME as VELOCITY_NAME,RECEIPT as VELOCITY_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,INHERITED_GATES as EARLIER_GATES,
    GATES as VELOCITY_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces_operator import (
    SEAMS,physical_boundary_germ,physical_contribution_groups,exact_modified_velocity_pressure_interfaces_theorem)

NAME=PREFIX+'current_modified_velocity_pressure_interfaces.json'
RECEIPT=PREFIX+'current_modified_velocity_pressure_interfaces_check.json'
VIEWS_NAME=PREFIX+'current_modified_velocity_pressure_interfaces_views.json.gz'
GATES=('current_modified_full_spatial4_temporal1_interfaces_certified',
  'current_modified_twelve_affected_velocity_pressure_spatial4_time1_traces_available',
  'current_modified_original_and_own_signed_boundary_modes_preserved',
  'current_modified_common_216_velocity_pressure_interface_bounds_available')
OPEN=tuple(key for key in EARLIER_OPEN if key not in GATES)+(
  'current_modified_full_33_chart_velocity_pressure_dispatch_installed',)
INHERITED_GATES=EARLIER_GATES+VELOCITY_GATES

class CurrentModifiedVelocityPressureInterfaces:
    @source_precision
    def __init__(self,velocity=None,require_checked=True):
        self.velocity=velocity if velocity is not None else CurrentModifiedVelocityPressure()
        if type(self.velocity) is not CurrentModifiedVelocityPressure or not self.velocity.acceptance_loaded:
            raise ValueError('Checked actual own velocity/absolute-pressure mixed4 packet required')
        self.interfaces=self.velocity.heat.interfaces;self.controls=self.velocity.source.repair.controls
        self.ctx=self.velocity.ctx;self.delta=self.velocity.delta;self.N=self.velocity.N;self.mu=self.velocity.mu
        self.theorem=exact_modified_velocity_pressure_interfaces_theorem(self)
        self.hashes={**self.velocity.hashes,**self.theorem['input_hashes'],VELOCITY_NAME:sha(VELOCITY_NAME),VELOCITY_RECEIPT:sha(VELOCITY_RECEIPT)}
        for stem in ('current_modified_velocity_pressure_interfaces_operator','current_modified_velocity_pressure_interfaces'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_velocity_pressure_definition=self.velocity.modified_velocity_pressure_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_velocity_pressure_interfaces_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_velocity_pressure_interfaces.py'))
        self.modified_velocity_pressure_interfaces_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Affected physical interface receipt exceeds twelve-germ scope')
            if receipt['modified_velocity_pressure_interfaces_definition_sha256']!=self.modified_velocity_pressure_interfaces_definition_sha256:
                raise ValueError('Foreign two-germ source adapter, controls or physical source')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual affected physical interface manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None):
        if name not in SEAMS:raise ValueError('One of twelve affected physical interfaces required')
        self.velocity.heat.assert_graph()
        if self.interfaces.repaired.repair.controls is not self.controls:
            raise ValueError('Actual implicit boundary source controls replaced')
        Z,_,lt,_=self.velocity.heat.dispatch._arguments(Z,0,log_tau,theta,1)
        left=physical_boundary_germ(self,name,Z,'left',lt,theta)
        right=physical_boundary_germ(self,name,Z,'right',lt,theta)
        a=physical_contribution_groups(left);b=physical_contribution_groups(right);common={}
        if set(a)!=set(b):raise ValueError('Physical two-sided contribution layouts differ')
        for key in a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in (a[key],b[key]) if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,
              log_absolute_upper=self.ctx.mpf(max(values)) if values else None,
              both_actual_side_bounds_retained=True)
        return dict(seam=name,Z=Z,requested_log_tau=lt,
          left_complete_velocity_pressure_source=left,right_complete_velocity_pressure_source=right,
          common_physical_velocity_pressure_rows=common,current_common_physical_contribution_count=len(common),
          modified_velocity_pressure_interfaces_definition_sha256=self.modified_velocity_pressure_interfaces_definition_sha256,
          same_unique_implicit_repair_definition_sha256=self.velocity.source.repair.repair_definition_sha256,
          source_function_and_FTC_identities_precede_physical_bounds=True,
          interval_overlap_or_midpoint_not_used_as_function_identity=True,
          actual_P0_cumulative_moments_and_all_signed_boundary_units_retained=True,
          affected_interface_count=12,source_total_mixed_order=4,physical_spatial_order=4,physical_fixed_x_time_order=1,
          full33_physical_velocity_dispatch_and_global_interfaces_not_claimed=True,
          resolved_physical_point_values_available=False,
          **{key:self.acceptance_loaded for key in INHERITED_GATES+GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.velocity.heat.registry.family,
          implicit_source_sha256=self.velocity.heat.registry.source,datum_enclosure_sha256=self.velocity.heat.registry.datum_sha,
          finite_integer_N=self.N,modified_velocity_pressure_interfaces_definition_sha256=self.modified_velocity_pressure_interfaces_definition_sha256,
          parent_modified_velocity_pressure_definition_sha256=self.velocity.modified_velocity_pressure_definition_sha256,
          exact_actual_modified_velocity_pressure_interface_theorem=self.theorem,
          affected_physical_interface_names=list(SEAMS),actual_two_germ_physical_views=VIEWS_NAME,
          complete_physical_contribution_groups_per_seam=216,
          original_native_regions=33,original_adjacent_interfaces=32,original_internal_tensor_traces=14,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='All twelve actual affected velocity/absolute-pressure two-germ interfaces: original O2 turnoff/buffer, buffer/transition, transition/power, two modulation flat edges, six exact repair support edges and q2 exact implicit-source exit. Whole-Z ordinary mixed4 functions and original Cartesian spatial4/fixed-x time1 traces, with216 common contribution groups each. Actual original radial/axial Pstar0 and own radial Pstar0+1/theta/axial Pstar1/pressure Pstar2, original P0/moments/partial FTC and exact radius retained. Full33 velocity dispatch/all global physical interfaces, total physical energy, common N/modified cones, recursion/waves/resolved global corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(INHERITED_GATES+GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedVelocityPressureInterfaces(require_checked=False)
    result=field.manifest();views={}
    for name in SEAMS:
        views[name]=field.interface(name)
        print('Actual velocity/pressure physical spatial4/time1 interface: '+name,flush=True)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result

if __name__=='__main__':run()
