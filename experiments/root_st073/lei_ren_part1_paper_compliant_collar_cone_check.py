"""Focused current-source check of the whole original collar cone."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_cone import sigmoid_cone_identities,sigmoid_source_bridge
from lei_ren_part1_paper_compliant_collar_heat_cone import heat_cone_identities
from lei_ren_part1_paper_compliant_collar_stress_C3 import CompliantCollarStressC3
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def run():
    name=PREFIX+'collar_cone.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
    for source,digest in hashes.items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Whole collar cone source changed: '+source)
    stress=CompliantCollarStressC3()
    if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
        raise ValueError('Whole collar cone family mismatch')
    if sigmoid_cone_identities()!=record['sigmoid_exact_identities']:raise ValueError('Sigmoid exact identities changed')
    if heat_cone_identities()!=record['heat_exact_identities']:raise ValueError('Heat exact identities changed')
    if sigmoid_source_bridge(stress)!=record['source_bridge']:raise ValueError('Actual sigmoid source bridge changed')
    c=MPIntervalContext(); c.dps=280; bounds=record['sigmoid_bounds']; count=0
    for label,value in bounds['positive_margins'].items():
        lo,hi=endpoints(read_interval(c,value))
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Nonpositive actual sigmoid margin: '+label)
        count+=1
    if endpoints(read_interval(c,bounds['positive_margins']['epsilon']))[1]>=1:raise ArithmeticError('Original positive cutoff factor lost')
    upper=endpoints(read_interval(c,bounds['log_directional_term_upper']))[1]
    threshold=endpoints(read_interval(c,bounds['log_threshold']))[0]
    if not mp.isfinite(upper) or upper>=threshold:raise ArithmeticError('Whole sigmoid directional inequality failed')
    if record['domain']!=dict(offset=[0,3],Z=[-1,1],strict_upper_offset_excluded=True):raise ValueError('Whole collar domain changed')
    for flag in ('sigmoid_transition_0_to_1_cone_certified','actual_heat_part_collar_cone_certified',
                 'same_source_offset1_cone_proofs_composed','collar_cone_certified',
                 'whole_collar_stress_theta_positive_before_zero_join','zero_stress_endpoint3_excluded_from_strict_inequalities',
                 'whole_collar_uniform_terminal_direction_verified','fixed_positive_viscosity_and_lambda_factors_preserve_two_vector_cone'):
        if not record[flag]:raise ValueError('Whole collar cone fact missing: '+flag)
    for flag in ('completed_full_tensor_cone_certified','global_admissible_stress_lift_constructed',
                 'independently_bounded_global_flat_remainder','temporal_recursion'):
        if record[flag]:raise ValueError('Whole collar source cone scope overclaimed: '+flag)
    hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                collar_cone_certified=True,sigmoid_transition_0_to_1_cone_certified=True,
                actual_heat_part_collar_cone_certified=True,positive_source_margins_checked=count,
                whole_original_collar_source_cone_domain_verified=record['domain'],
                sigmoid_exact_identities=record['sigmoid_exact_identities'],source_bridge=record['source_bridge'],
                source_offset1_join_and_endpoint3_uniform_direction_verified=True,
                fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                completed_full_tensor_cone_certified=False,global_admissible_stress_lift_constructed=False,
                independently_bounded_global_flat_remainder=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual whole original collar two-vector cone PASS; global stress/flat remainder pending',flush=True)
    return result


if __name__=='__main__':run()
