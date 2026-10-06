"""Local modified physical completed tensor/divergence/remainder source bounds."""
import gzip
import hashlib
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_modified_pre_stress import (
    CurrentModifiedPreStress,HERE,PREFIX,NAME as SIGNED_NAME,RECEIPT as SIGNED_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,GATES as SIGNED_GATES,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import (
    compiled_modified_pre_physical_lift,modified_physical_packet,exact_modified_pre_physical_theorem)

NAME=PREFIX+'current_modified_pre_physical_tensor.json'
RECEIPT=PREFIX+'current_modified_pre_physical_tensor_check.json'
VIEWS_NAME=PREFIX+'current_modified_pre_physical_tensor_views.json.gz'
GATES=('current_modified_local_physical_completed_tensor_source_available',
       'current_modified_local_diagonal_radial_divergence_completion_installed',
       'current_modified_diagonal_radial_tensor_and_divergence_completion_installed',
       'current_modified_local_full_physical_three_component_remainder_available',
       'current_modified_own_moment_physical_incompressibility_certified')
OPEN=tuple(key for key in EARLIER_OPEN if key!='current_modified_diagonal_radial_tensor_and_divergence_completion_installed')+(
  'current_modified_full_33_chart_physical_tensor_dispatch_installed',)

class CurrentModifiedPrePhysicalTensor:
    @source_precision
    def __init__(self,source=None,require_checked=True):
        self.source=source if source is not None else CurrentModifiedPreStress()
        if type(self.source) is not CurrentModifiedPreStress or not self.source.acceptance_loaded:
            raise ValueError('Checked actual signed modified pre stress/remainder required')
        self.ctx=self.source.ctx;self.delta=self.source.delta;self.N=self.source.N;self.mu=self.source.mu
        self.registry=self.source.source.histories.candidate.direction.registry
        self.geometry=self.registry.owners['o3'].physical
        self.lift,self.adaptation=compiled_modified_pre_physical_lift();self.theorem=exact_modified_pre_physical_theorem(self)
        self.hashes={**self.source.hashes,**self.theorem['input_hashes'],**self.adaptation['input_hashes'],
          SIGNED_NAME:sha(SIGNED_NAME),SIGNED_RECEIPT:sha(SIGNED_RECEIPT)}
        for stem in ('current_modified_pre_physical_tensor_operator','current_modified_pre_physical_tensor'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_stress_definition=self.source.modified_stress_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_pre_physical_tensor_operator.py'),
          source_sha256=sha(PREFIX+'current_modified_pre_physical_tensor.py'))
        self.modified_physical_tensor_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Modified local physical tensor receipt exceeds scope')
            if receipt['modified_physical_tensor_definition_sha256']!=self.modified_physical_tensor_definition_sha256:
                raise ValueError('Foreign signed modified source or physical operator')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Actual modified physical tensor manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def tensor(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        c=self.ctx;lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(x) for x in endpoints(lt)+endpoints(nu)) or endpoints(nu)[0]<=0:
            raise ValueError('Finite log(tau), finite nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        signed=self.source.stress(region,Z,coordinate)
        if not all(signed[key] for key in SIGNED_GATES):raise ValueError('Accepted same-source signed columns/remainder required')
        packet,radius=modified_physical_packet(self,signed)
        point=self.lift(c,packet,self.delta,None,lt,theta,nu)
        return dict(point,region=region,coordinate=c.mpf(coordinate),
          modified_source_definition_sha256=self.source.source.modified_source_definition_sha256,
          modified_stress_definition_sha256=self.source.modified_stress_definition_sha256,
          modified_physical_tensor_definition_sha256=self.modified_physical_tensor_definition_sha256,
          current_actual_modified_signed_source=signed,current_actual_modified_physical_packet=packet,
          exact_original_geometry_radius_source=radius,
          local_physical_source_bounds_not_resolved_physical_point_values=True,
          exact_zero_velocity_divergence_uses_own_moment_FTC_source_identity=True,
          own_physical_total_kinetic_energy_is_separate_required_construction=True,
          completed_local_tensor_not_yet_admitted_to_stress_cone=True,
          **{key:self.acceptance_loaded for key in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.source.source.histories.family,
          implicit_source_sha256=self.source.source.histories.source,datum_enclosure_sha256=self.source.source.histories.datum_sha,
          parent_modified_source_definition_sha256=self.source.source.modified_source_definition_sha256,
          parent_modified_stress_definition_sha256=self.source.modified_stress_definition_sha256,
          modified_physical_tensor_definition_sha256=self.modified_physical_tensor_definition_sha256,
          finite_integer_N=self.N,exact_actual_modified_physical_source_theorem=self.theorem,
          source_remainder_only_original_physical_lift_adaptation=self.adaptation,
          actual_modified_completed_physical_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Local actual modified O2 buffer/O3 slope/quiet q<=2 completed physical cylindrical stress3, diagonal/divergence2, Cartesian tensor, full signed physical remainder2 and momentum -div(T)+E source bounds. Original geometry/logPstar, own pressure/moments and all15 stress/11 remainder sectors retained. Own moment FTC proves velocity incompressibility; exact r*d_z(Trz) completion cancels radial tensor divergence. Resolved point values, open power continuation/full33 modified dispatcher, source/physical/heat joins, kinetic energy, common cone N, global admissibility, recursion/waves/full corrected NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

QUERIES=(('whole_O2_buffer','O2_buffer',(0,11),('-3','-1'),None,'1'),
  ('whole_O3_transition','O3_slope_mu',(0,1),('-3','-1'),None,'1'),
  ('whole_quiet_repair','quiet_O3_power',(1,2),('-3','-1'),None,'1'),
  ('first_overlap','quiet_O3_power','1.201337','-2.6','.41','.8'),
  ('middle_swirl','quiet_O3_power','1.501337','-2.6','.41','.8'),
  ('last_overlap','quiet_O3_power','1.801337','-2.6','.41','.8'),
  ('exact_exit','quiet_O3_power',2,'-1',None,'1'))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPrePhysicalTensor(require_checked=False)
    result=field.manifest();views={key:field.tensor(region,(-1,1),q,lt,theta,nu) for key,region,q,lt,theta,nu in QUERIES}
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual modified local completed physical tensor/divergence/remainder source constructed; cone/energy/global NS open',flush=True)
    return result

if __name__=='__main__':run()
