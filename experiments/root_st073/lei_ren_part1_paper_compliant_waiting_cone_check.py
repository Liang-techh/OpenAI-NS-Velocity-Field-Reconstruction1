"""Whole waiting cone: current common sources, mode theorem and fresh bounds."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_waiting_cone import (
    CompliantWaitingStressC3,CompliantGlobalPhysicalAssembly,DOMAIN,
    waiting_cone_identities,waiting_cone_source_bridge,whole_waiting_bounds)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def run():
    with mp.workdps(300):
        name=PREFIX+'waiting_cone.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting cone source changed: '+source)
        stress=CompliantWaitingStressC3(); assembly=CompliantGlobalPhysicalAssembly(); c=stress.ctx
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Waiting cone actual source family differs')
        if waiting_cone_identities()!=record['exact_identities']:raise ValueError('Waiting monotone mode/cone identities changed')
        if waiting_cone_source_bridge(stress)!=record['source_bridge']:raise ValueError('Actual full-history or log source bridge changed')
        # Re-evaluate the whole source Z box and certified bound arithmetic;
        # the functional theorem covers every q, without phase sampling.
        bounds=whole_waiting_bounds(stress,assembly)
        if encode(pack(bounds))!=record['bounds']:raise ValueError('Current whole-source waiting margins changed')
        count=0
        for label,value in bounds['positive_margins'].items():
            if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Waiting positive margin failed: '+label)
            count+=1
        if endpoints(bounds['theta_uniform_positive_lower'])[0]<=0:raise ArithmeticError('Waiting stress theta sign not strict')
        if endpoints(bounds['log_directional_term_upper'])[1]>=endpoints(c.ln(2))[0]:raise ArithmeticError('Waiting directional inequality failed')
        if record['domain']!=DOMAIN:raise ValueError('Whole original waiting cone domain changed')
        for flag in ('waiting_cone_certified','actual_waiting_nonzero_stress_theta_positive',
                     'actual_waiting_source_shear_strictly_negative','actual_waiting_directional_cone_uniform_margin',
                     'waiting_collar_cone_proofs_same_source_composed','waiting_and_collar_outer_tail_two_vector_cone_certified',
                     'fixed_positive_viscosity_two_vector_cone_transfer_verified'):
            if not record[flag]:raise ValueError('Waiting cone fact missing: '+flag)
        for flag in ('completed_full_tensor_cone_certified','whole_outer_cone_certified',
                     'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
                     'physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Waiting cone scope overclaimed: '+flag)
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    waiting_cone_certified=True,domain=DOMAIN,positive_source_margins_checked=count,
                    whole_Z_source_bounds_recomputed=True,all_waiting_phases_covered_by_exact_mode_theorem=True,
                    actual_waiting_source_shear_strictly_negative=True,
                    actual_Bmax_is_original_source_inlet_verified=True,
                    waiting_and_collar_outer_tail_two_vector_cone_certified=True,
                    same_source_waiting_collar_stress_mixed3_join_consumed=True,
                    fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                    completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS whole original waiting cone and same-source collar composition; global cone pending',flush=True)
        return result


if __name__=='__main__':run()
