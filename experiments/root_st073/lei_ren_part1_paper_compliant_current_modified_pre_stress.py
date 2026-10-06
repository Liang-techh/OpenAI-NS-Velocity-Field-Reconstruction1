"""Callable signed stress columns and full local remainder for modified O2/O3."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_O3_repaired_histories import (
    CurrentO3RepairedHistories,HERE,PREFIX,NAME as SOURCE_NAME,RECEIPT as SOURCE_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes,GATES as SOURCE_GATES,OPEN as SOURCE_OPEN)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_modified_pre_stress_operator import (
    modified_pre_stress_rows,modified_pre_remainder_sectors,exact_modified_pre_stress_theorem)

NAME=PREFIX+'current_modified_pre_stress.json'
RECEIPT=PREFIX+'current_modified_pre_stress_check.json'
VIEWS_NAME=PREFIX+'current_modified_pre_stress_views.json.gz'
GATES=('current_modified_full_paper_theta_axial_stress3_signed_sectors_available',
       'current_modified_full_source_remainder2_signed_sectors_available',
       'current_modified_stress_remainder_actual_Pstar_source_modes_certified')
OPEN=SOURCE_OPEN+('current_modified_own_complete_energy_history_installed',
    'current_modified_diagonal_radial_tensor_and_divergence_completion_installed')

class CurrentModifiedPreStress:
    @source_precision
    def __init__(self,source=None,require_checked=True):
        self.source=source if source is not None else CurrentO3RepairedHistories()
        if type(self.source) is not CurrentO3RepairedHistories or not self.source.acceptance_loaded:
            raise ValueError('Checked own repaired field source required')
        self.source.histories.candidate.direction.assert_graph();self.ctx=self.source.ctx
        self.delta=self.source.histories.delta;self.N=self.source.N;self.mu=self.source.mu
        self.theorem=exact_modified_pre_stress_theorem()
        self.hashes={**self.source.hashes,**self.theorem['input_hashes'],SOURCE_NAME:sha(SOURCE_NAME),SOURCE_RECEIPT:sha(SOURCE_RECEIPT)}
        for stem in ('current_modified_pre_stress_operator','current_modified_pre_stress'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_source_definition=self.source.modified_source_definition_sha256,
          operator_sha256=sha(PREFIX+'current_modified_pre_stress_operator.py'),source_sha256=sha(PREFIX+'current_modified_pre_stress.py'))
        self.modified_stress_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Modified local column/remainder receipt scope changed')
            if receipt['modified_stress_definition_sha256']!=self.modified_stress_definition_sha256:
                raise ValueError('Foreign repaired source or signed stress operator')
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=json.loads((HERE/NAME).read_bytes()):raise ValueError('Modified stress source manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def stress(self,region,Z,coordinate):
        view=self.source.history(region,Z,coordinate);c=self.ctx;z=IntervalTaylor.variable(c,c.mpf(Z),5)
        if not all(view[k] for k in SOURCE_GATES) or any(view[k] for k in SOURCE_OPEN):
            raise ValueError('Signed columns require the accepted scoped modified source')
        columns=modified_pre_stress_rows(c,self.delta,z,view)
        remainder=modified_pre_remainder_sectors(c,self.delta,z,view)
        return dict(region=region,coordinate=c.mpf(coordinate),Z=c.mpf(Z),
          modified_source_definition_sha256=self.source.modified_source_definition_sha256,
          same_unique_implicit_repair_definition_sha256=self.source.repair.repair_definition_sha256,
          modified_stress_definition_sha256=self.modified_stress_definition_sha256,
          actual_source=view,full_signed_paper_theta_axial_stress3_sectors=columns,
          full_signed_source_remainder2_sectors=remainder,
          stress_common_positive_radius_half_normalization='R^mode[0]*Pstar^mode[1]/sqrt(2)',
          remainder_original_physical_beta_radius_and_half_modes_preserved=True,
          no_actual_Pstar_or_radius_cap_materialized=True,
          actual_own_absolute_pressure_and_five_histories_retained=True,
          raw_physical_Uz_Pstar_sector_retained=True,
          completed_tensor_energy_physical_NS_and_cone_are_separate_open_work=True,
          **{k:self.acceptance_loaded for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.source.histories.family,
          implicit_source_sha256=self.source.histories.source,datum_enclosure_sha256=self.source.histories.datum_sha,
          parent_modified_source_definition_sha256=self.source.modified_source_definition_sha256,
          modified_stress_definition_sha256=self.modified_stress_definition_sha256,
          finite_integer_N=self.N,exact_signed_modified_stress_remainder_source_theorem=self.theorem,
          actual_signed_modified_source_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Actual modified O2 buffer/O3 variable slope/quiet repaired power signed original paper theta/axial columns through ordinary logR3 and full source remainder through logR2, with physical Uz=Pstar*Vhat, m/k Pstar0+1 and radial Pstar0+1, all nonlinear cross/square terms and physical beta/radius modes. Own complete energy, diagonal/radial tensor, completed divergence/physical/heat joins, common cone N, modified cones/global tensor, recursion/waves/full NS remain open.',
          input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

QUERIES=(('O2_full_buffer','O2_buffer',(0,11)),('O3_full_transition','O3_slope_mu',(0,1)),
  ('quiet_full_repair','quiet_O3_power',(1,2)),('first_overlap','quiet_O3_power','1.201337'),
  ('middle_swirl','quiet_O3_power','1.501337'),('last_overlap','quiet_O3_power','1.801337'),
  ('exact_exit','quiet_O3_power',2))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedPreStress(require_checked=False)
    result=field.manifest();views={key:field.stress(region,(-1,1),q) for key,region,q in QUERIES}
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual modified signed paper stress15/remainder11 sectors constructed; own completed energy/tensor/cone open',flush=True)
    return result

if __name__=='__main__':run()
