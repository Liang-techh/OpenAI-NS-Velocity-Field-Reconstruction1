"""Local source profile and density checks; no modified global field admission."""
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import (
    HERE,NAME,RECEIPT,OPEN,sha,pack,encode,endpoints,source_precision,_verify_hashes,cutoff_rows)

@source_precision
def run(field):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest()
    expected['local_candidate_examples']=[field.profile(region,'.7',v,10**12) for region,v in (
       ('O2_buffer',9),('O2_buffer',10),('O2_buffer',11),
       ('O3_slope_mu',0),('O3_slope_mu','.1337000000001337'),('O3_slope_mu','.4'),('O3_slope_mu','.5'),('O3_slope_mu',1))]
    if encode(pack(expected))!=raw or any(raw[k] for k in OPEN):raise ValueError('Actual finite profile source/scope differs')
    if not field.theorem['passed'] or not all(field.theorem['identities'].values()):raise ValueError('Exact finite-N theorem required')
    if endpoints(field.supported_loop_vs_excess_lower)[0]<=0:raise ValueError('Supported periodic loop floor required')
    c=field.ctx
    for example in expected['local_candidate_examples']:
        if any(example[k] for k in OPEN) or not example['moments_pressure_radial_recovery_and_completed_tensor_not_installed']:
            raise ValueError('Local profiles cannot admit a complete modified field')
        if not example['phase_is_actual_translated_N_logR'] or not example['original_axial_is_exact_zero']:
            raise ValueError('Actual phase and zero axial source required')
        if set(example['actual_five_moment_increment_densities_per_dX'])!={'M','I','J','S','Cp'}:
            raise ValueError('All five increment densities required')
        e=example['modified_theta_over_Pstar_ordinary_logR_axial5']
        v=example['modified_axial_over_Pstar_ordinary_logR_axial5']
        if endpoints(e[0][0])[0]<=0:raise ValueError('Actual finite-N swirl positivity unresolved')
        shear=example['finite_N_shear']
        for name,residual in (
          ('a',1-2*e[1][0]/e[0][0]-2-shear['a_minus2']),
          ('b',2*v[1][0]/e[0][0]-shear['b'])):
            lo,hi=endpoints(residual)
            if not lo<=0<=hi:raise ValueError('Finite-N numerical chain-rule consistency failed: '+name)
        if example['original_profile_jets_unchanged_outside_support']:
            if any(endpoints(row)!=(0,0) for row in example['spatial_cutoff_ordinary_logR_rows']):
                raise ValueError('All cutoff jets must be exactly flat at support edges')
            if encode(pack(e))!=encode(pack(example['original_E0_over_Pstar_rows'])):
                raise ValueError('Original swirl jets changed off support')
            if any(any(endpoints(a)!=(0,0) for a in row.coefficients) for row in v):
                raise ValueError('Axial modulation not zero off support')
    # The source-function join follows from the inherited original O2/O3
    # join and the exact common-coordinate/phase theorem. Byte equality
    # below additionally checks this callable implementation at a fresh Z.
    left=field.profile('O2_buffer','.537',11,11);right=field.profile('O3_slope_mu','.537',0,11)
    for key in ('spatial_cutoff_ordinary_logR_rows','modified_theta_over_Pstar_ordinary_logR_axial5',
                'modified_axial_over_Pstar_ordinary_logR_axial5','finite_N_shear'):
        if encode(pack(left[key]))!=encode(pack(right[key])):raise ValueError('Seam-crossing profile implementation differs: '+key)
    for bad in (0,-1,True,1.5):
        try:field.profile('O3_slope_mu',0,0,bad)
        except ValueError:pass
        else:raise ValueError('Invalid frequency admitted')
    for region,z,v in (('foreign',0,0),('O3_slope_mu',0,-.1),('O2_buffer',0,12),('O3_slope_mu',1.1,0)):
        try:field.profile(region,z,v,11)
        except ValueError:pass
        else:raise ValueError('Invalid original source chart admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
      datum_enclosure_sha256=field.datum_sha,exact_finite_frequency_chain_rule_moment_and_coordinate_identities=len(field.theorem['identities']),
      local_profile_example_count=8,
      actual_seam_crossing_periodic_phase_and_original_flat_profile_edges_checked=True,
      same_current_source_Pstar_factored_before_profile_and_density_evaluation=True,
      original_five_moments_and_pressure_not_copied_for_modified_profiles=True,
      numerical_chain_rule_consistency_checked_after_exact_function_identity=True,
      source_samples_not_a_uniform_finite_N_cone_proof=True,
      integrated_five_moments_pressure_radial_recovery_and_uniform_N_remain_open=True,
      current_strict_nonzero_whole_regions_including_inherited=15,remaining_registry_regions_without_current_whole_cone=17,
      input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual finite-N seam-crossing profile/density construction PASS; no modified field/cone admission',flush=True)
    return result
