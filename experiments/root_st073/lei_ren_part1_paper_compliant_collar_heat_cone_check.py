"""Focused source/proof check of the whole actual heat-part collar cone."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_heat_cone import heat_cone_identities,source_bindings
from lei_ren_part1_paper_compliant_collar_stress_C3 import CompliantCollarStressC3
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def run():
    name=PREFIX+'collar_heat_cone.json'; record=json.loads((HERE/name).read_bytes())
    hashes=dict(record['input_hashes'])
    for source,digest in hashes.items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Heat cone dependency changed: '+source)
    stress=CompliantCollarStressC3()
    if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
        raise ValueError('Heat cone source family changed')
    if heat_cone_identities()!=record['exact_identities']:raise ValueError('Heat cone exact proof changed')
    if source_bindings(stress)!=record['source_bridge']:raise ValueError('Heat cone source bridge changed')
    c=MPIntervalContext(); c.dps=280; bounds=record['bounds']; margins=0
    for label in ('H_lower','angular_future_integrand_lower','angular_shear_perturbation_lower','L_lower','epsilon','delta'):
        if endpoints(read_interval(c,bounds[label]))[0]<=0:raise ArithmeticError('Nonpositive cone prerequisite: '+label)
        margins+=1
    if endpoints(read_interval(c,bounds['epsilon']))[1]>=1:raise ArithmeticError('Nonpositive cutoff factor')
    upper=endpoints(read_interval(c,bounds['log_directional_cone_term_upper']))[1]
    threshold=endpoints(read_interval(c,bounds['log_directional_threshold']))[0]
    if not mp.isfinite(upper) or upper>=threshold:raise ArithmeticError('Whole-domain direction margin failed')
    if record['domain']!=dict(offset=[1,3],Z=[-1,1],strict_upper_offset_excluded=True):raise ValueError('Cone scope changed')
    for flag in ('actual_heat_part_collar_cone_certified','actual_nonzero_stress_theta_positive',
                 'actual_directional_cone_uniform_margin_certified','actual_zero_stress_join_uniform_direction_e_theta_verified',
                 'source_defined_positive_eps_phi_S_not_replaced_by_caps','endpoint3_exact_zero_stress_not_subject_to_strict_cone'):
        if not record[flag]:raise ValueError('Heat-part proof missing: '+flag)
    for flag in ('sigmoid_transition_0_to_1_cone_certified','collar_cone_certified','completed_full_tensor_cone_certified',
                 'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder','temporal_recursion'):
        if record[flag]:raise ValueError('Heat-part cone scope overclaimed: '+flag)
    hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                actual_heat_part_collar_cone_certified=True,positive_source_margins_checked=margins,
                exact_identities=record['exact_identities'],source_bridge=record['source_bridge'],
                whole_domain_log_direction_bound_checked=True,uniform_zero_stress_join_direction_verified=True,
                sigmoid_transition_0_to_1_cone_certified=False,collar_cone_certified=False,
                completed_full_tensor_cone_certified=False,global_admissible_stress_lift_constructed=False,
                independently_bounded_global_flat_remainder=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual heat-part collar cone and uniform flat join direction PASS; whole collar pending',flush=True)
    return result


if __name__=='__main__':run()
