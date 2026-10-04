"""Focused original-flatten physical source/right join; consume general oracle."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_flatten_physical_C2 import (
    CompliantFlattenStressC3,flatten_physical_identities,flatten_adapter_binding)
from lei_ren_part1_paper_compliant_outer_power_physical_C2_check import accepted_general_Cartesian_oracle
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def run():
    with mp.workdps(300):
        name=PREFIX+'flatten_physical_C2.json';record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for path,digest in hashes.items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten physical source changed: '+path)
        provider=CompliantFlattenStressC3()
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(provider.family,provider.source):
            raise ValueError('Flatten physical family changed')
        if record['flatten_physical_identities']!=flatten_physical_identities():raise ValueError('Flatten physical source identities changed')
        if record['flatten_adapter_binding']!=flatten_adapter_binding(provider):raise ValueError('Flatten physical adapter/right join changed')
        c=MPIntervalContext();c.dps=300;finite=zeros=0
        def rowcheck(row):
            nonlocal finite,zeros
            lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
            if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Nonfinite flatten physical row')
            if row['exact_zero']:
                if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid flatten physical exact zero row')
                zeros+=1
            else:
                if any(not mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):
                    raise ArithmeticError('Nonfinite flatten physical log bound')
                finite+=1
        for point in record['samples']+[record['whole_flatten']]:
            for flag in ('actual_K_Z_and_K_ZZ_and_full_moment_axial_dependence_retained',
                'actual_regional_physical_flatten_stress_remainder_identity_verified',
                'full_axial_moment_history_and_pressure_baselines_cancelled_before_enclosure',
                'flatten_power_physical_stress_mixed3_and_remainder_mixed2_join_verified','ordinary_t_equals_q_derivatives'):
                if not point[flag]:raise ValueError('Flatten physical gate missing: '+flag)
            for group,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6)):
                for grid in point[group].values():
                    if len(grid)!=count:raise ValueError('Flatten physical mixed coverage incomplete')
                    for row in grid.values():rowcheck(row)
            for group in ('completed_theta_theta_stress_mixed2','physical_angular_axial_viscosity_remainder_mixed2'):
                if len(point[group])!=6:raise ValueError('Flatten physical mixed2 operator coverage incomplete')
                for row in point[group].values():rowcheck(row)
            for label in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
                if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Flatten physical structural zero lost')
                zeros+=1
            for flag in ('flatten_cone_certified','flatten_regional_remainder_exact_zero','whole_outer_cone_certified',
                'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
                'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion','left_pulse_physical_stress_join_verified'):
                if point[flag]:raise ValueError('Flatten physical scope overclaimed: '+flag)
        oracle=accepted_general_Cartesian_oracle(hashes)
        oracle.pop('actual_outer_power_source_admitted_by_new_adapter_and_stress_receipt',None)
        oracle['actual_flatten_source_admitted_by_new_adapter_and_stress_receipt']=True
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
            implicit_source_sha256=record['implicit_source_sha256'],
            actual_regional_physical_flatten_stress_remainder_identity_verified=True,
            flatten_power_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            actual_K_Z_and_K_ZZ_and_full_moment_axial_dependence_retained=True,ordinary_t_equals_q_derivatives=True,
            actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
            accepted_general_Cartesian_oracle=oracle,
            actual_flatten_source_specialization_identities=len(record['flatten_physical_identities']['source_specialization_identities']),
            actual_native_velocity_and_remainder_source_join_identities=len(record['flatten_adapter_binding']['actual_native_velocity_pressure_and_remainder_source_join']['identities']),
            actual_flatten_power_K4_full_moment_stress3_pressure4_source_join_consumed=True,
            flatten_regional_remainder_exact_zero=False,flatten_cone_certified=False,
            left_pulse_physical_stress_join_verified=False,global_admissible_stress_lift_constructed=False,
            independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS original100-unit flatten completed physical stress, variable K, retained remainder and power right join; cone/left pulse pending',flush=True)
        return result


if __name__=='__main__':run()
