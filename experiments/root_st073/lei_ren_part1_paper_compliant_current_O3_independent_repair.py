"""Source-bound independent five-moment implicit repair in quiet O3 power.

Certifies new fixed bumps, the actual scaled nonlinear map and its unique
constant controls. Installing partial repaired histories/tensors is next.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_O3_modulated_histories import (
    CurrentO3ModulatedHistories,HERE,PREFIX,NAME as HISTORY_NAME,RECEIPT as HISTORY_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import (
    CENTERS,RADIUS,ROWS,CONTROLS,exact_repair_theorem,bump_weights,inverse_and_bounds,
    correlated_scaled_defects,contraction_certificate,coefficient_enclosure,scaled_quadratic)

NAME=PREFIX+'current_O3_independent_repair.json'
RECEIPT=PREFIX+'current_O3_independent_repair_check.json'
GATES=('current_O3_quiet_power_independent_fixed_bump_map_available',
       'current_O3_actual_correlated_axial_repair_inverse_certified',
       'current_O3_modulated_unique_implicit_constant_repair_controls_certified')
OPEN=('current_O3_transition_modified_five_moment_repair_certified',
      'current_modified_repaired_partial_histories_pressure_radial_source_installed',
      'current_O3_transition_finite_N_modified_profiles_installed',
      'current_whole_O3_transition_signed_two_vector_cone_certified',
      'current_O2_buffer_modified_whole_cone_certified',
      'current_modified_complete_tensor_and_physical_source_joins_certified',
      'current_modified_analytic_preheat_exterior_compatible_after_repair',
      'current_modified_global_admissible_tensor_and_full_NS_certified')

class CurrentO3IndependentRepair:
    @source_precision
    def __init__(self,histories=None,require_checked=True,bump_cells=256):
        self.histories=histories if histories is not None else CurrentO3ModulatedHistories()
        if type(self.histories) is not CurrentO3ModulatedHistories or not self.histories.acceptance_loaded:
            raise ValueError('Checked actual modified histories required')
        self.histories.candidate.direction.assert_graph();self.ctx=c=self.histories.ctx
        if isinstance(bump_cells,bool) or not isinstance(bump_cells,int) or bump_cells<64 or bump_cells>2048 or bump_cells&(bump_cells-1):
            raise ValueError('Fixed power-of-two bump range subdivisions in[64,2048] required')
        self.bump_cells=bump_cells;self.mu=self.histories.mu;self.N=self.histories.N
        if not 0<endpoints(self.mu)[0]<=endpoints(self.mu)[1]<mp.mpf('1/6'):
            raise ValueError('Actual strictly positive small mu required')
        Tw=c.mpf(self.histories.pre.params.Tw)
        if endpoints(Tw)[0]<=2:raise ValueError('Independent quiet power band must fit actual source')
        original=self.histories.pre.slope_mu(0,1)
        if endpoints(original['original_transition_kernels']['J'])!=(mp.mpf('.5'),mp.mpf('.5')):
            raise ValueError('Actual source reflected J(1)=1/2 required')
        self.normalization=c.mpf(endpoints(self.histories.pre.initial.repair.normalization))
        if endpoints(self.normalization)[0]<=0:raise ValueError('Same exact raw bump normalization required')
        self.theorem=exact_repair_theorem()
        self.weights=bump_weights(c,self.mu,self.normalization,cells=bump_cells)
        self.matrix=inverse_and_bounds(c,self.mu,self.weights)
        self.certificate=contraction_certificate(c,self.mu,self.N,self.matrix)
        self.defects=correlated_scaled_defects(self.histories)
        bounds=self.certificate['scaled_defect_abs_bounds_all_integer_N_at_least1']
        for name,value,bound in zip(ROWS,self.defects,bounds):
            if max(abs(x) for x in endpoints(value))>endpoints(bound)[1]:
                raise ArithmeticError('Actual source defect exceeds scaled uniform bound: '+name)
        self.controls=coefficient_enclosure(c,self.mu,self.N,self.matrix,self.defects,self.certificate)
        self.hashes={**self.histories.hashes,**self.theorem['input_hashes'],HISTORY_NAME:sha(HISTORY_NAME),
          HISTORY_RECEIPT:sha(HISTORY_RECEIPT)}
        for stem in ('current_O3_independent_repair_operator','current_O3_independent_repair'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_modified_source_definition=self.histories.modified_source_definition_sha256,
          source_family=self.histories.family,source=self.histories.source,datum=self.histories.datum_sha,
          finite_integer_N=self.N,centers=CENTERS,radius=RADIUS,bump_cells=self.bump_cells,
          operator_sha256=sha(PREFIX+'current_O3_independent_repair_operator.py'),
          implementation_sha256=sha(PREFIX+'current_O3_independent_repair.py'))
        self.repair_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Independent implicit repair receipt scope differs')
            if receipt['repair_definition_sha256']!=self.repair_definition_sha256:
                raise ValueError('Foreign repair map or modified source')
            raw=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.manifest()))!=raw:raise ValueError('Actual independent implicit repair source differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    @source_precision
    def coefficients(self):
        c=self.ctx
        return dict(repair_definition_sha256=self.repair_definition_sha256,
          parent_modified_source_definition_sha256=self.histories.modified_source_definition_sha256,
          finite_integer_N=self.N,scaled_control_order=CONTROLS,scaled_constant_control_enclosures=self.controls,
          unique_exact_implicit_definition='B*h+Q_N(h,h)=-same_exact_scaled_source_defects; h in strict certified ball',
          exact_controls_are_Z_independent=True,exact_defects_not_replaced_by_interval_midpoints=True,
          physical_axial_coefficient_factor_over_same_A=c.sqrt(self.mu)/self.N,
          physical_swirl_coefficient_factor_over_same_A=self.mu/self.N,
          original_incoming_moments_preserved=True,
          all_five_terminal_increments_vanish_by_defining_map=True,
          terminal_identity_not_yet_connected_to_repaired_partial_history_source=True,
          **{k:self.acceptance_loaded for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.histories.family,
          implicit_source_sha256=self.histories.source,datum_enclosure_sha256=self.histories.datum_sha,
          parent_modified_source_definition_sha256=self.histories.modified_source_definition_sha256,
          repair_definition_sha256=self.repair_definition_sha256,finite_integer_N=self.N,
          fixed_bump_subdivisions=self.bump_cells,same_raw_bump_normalization=self.normalization,
          exact_source_and_new_five_bump_map_theorem=self.theorem,
          new_independent_bump_weights=self.weights,new_scaled_linear_inverse_quadratic_map=self.matrix,
          actual_signed_scaled_defect_enclosures=self.defects,
          exact_unique_scaled_constant_control_enclosures=self.controls,
          uniform_repair_only_contraction_certificate=self.certificate,
          scope='New fixed quiet O3 power log-translated 2+3 bumps, actual correlated (J-M)/mu source defects, common A(Z) normalization, scaled exact linear+quadratic map and unique constant implicit controls for the actual N. Source integrals are not midpoint substitutes. Partial repaired histories/pressure/radial/tensor joins and full physical repair installation remain open; repair-only N does not admit common modified cones/NS/recursion.',
          input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3IndependentRepair(require_checked=False)
    result=field.manifest()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('New independent five-bump map and actual source implicit constant controls constructed; field installation open',flush=True)
    return result

if __name__=='__main__':run()
