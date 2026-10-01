"""Refine the dominant accepted pressure mass without changing its source.

Keeps the old accepted receipts intact. A tighter enclosure is intersected
with the previous bound only after independent directed integration.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures, endpoints
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run(panels=128,order=24):
    base=Path(__file__).parent
    name='lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json'
    old=json.loads((base/name).read_text())
    start=time.monotonic()
    with mp.workdps(340):
        profile,alignment=accepted_profile()
        if old['accepted_schedule_sha256']!=alignment['accepted_schedule']['sha256']:
            raise AssertionError('Accepted source mismatch')
        e=ScheduleEndpointEnclosures(profile.schedule,precision=300)
        ctx=e.iv
        print('Refining slope_transition_ref pressure mass:',panels,order,flush=True)
        refined=HighOrderPreheatIntegrals(e).integrate_stage(
            'slope_transition_ref',panels=panels,order=order)
        rows={}
        old_sum=ctx.mpf(0);new_sum=ctx.mpf(0)
        for stage,row in old['stages'].items():
            if row['beta']!=2:continue
            original=read_interval(ctx,row['true_mass_interval'])
            new=original
            if stage=='slope_transition_ref':
                candidate=refined['normalized_mass_at_Z0_interval']
                left=max(endpoints(original)[0],endpoints(candidate)[0])
                right=min(endpoints(original)[1],endpoints(candidate)[1])
                if left>right:raise AssertionError('Independent mass enclosures disagree')
                new=ctx.mpf([left,right])
            old_sum+=original;new_sum+=new
            rows[stage]=dict(old_mass=original,refined_mass=new,
                width=endpoints(new)[1]-endpoints(new)[0])
        z=ctx.mpf('.3');dt=ctx.mpf('1e-200');q=1+z*z;L=1-dt*z*z
        # Exact pressure contribution to U1=Psi1 in the physical gauge.
        sensitivity=-ctx.exp(28)*((1-z*z)*(-4*z/q**3)-2*(1+dt)*z/q**2)/(2*L)
        radius=ctx.mpf('4.1');target=ctx.mpf('1e-12')
        def budget(mass):
            lo,hi=endpoints(mass)
            return sensitivity*ctx.mpf((hi-lo)/2)*radius
        previous_budget=budget(old_sum);new_budget=budget(new_sum)
        report=dict(accepted_schedule_sha256=alignment['accepted_schedule']['sha256'],
            input_hashes={name:hashlib.sha256((base/name).read_bytes()).hexdigest()},
            panels=panels,Taylor_order=order,precision=e.precision,
            stages=rows,old_beta2_mass=old_sum,refined_beta2_mass=new_sum,
            Psi_first_radial_coefficient_mass_sensitivity=sensitivity,
            scaled_radius=radius,Z=z,normalized_value_error_target=target,
            old_correlated_beta2_midpoint_value_error_upper=previous_budget,
            refined_correlated_beta2_midpoint_value_error_upper=new_budget,
            first_coefficient_source_budget_pass=endpoints(new_budget)[1]<=endpoints(target)[0],
            first_coefficient_budget_is_not_total_core_error=True,
            estimate_is_uncertainty_upper_bound_not_actual_error_lower_bound=True,
            pressure_source_unchanged=True,accepted_pressure_jets_updated=False,
            original_parameter_errors_enclosed=False,core_to_collar_matching_certified=False,
            temporal_recursion=False,elapsed_seconds=time.monotonic()-start)
        path=Path(__file__).with_suffix('.json')
        path.write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('Old/new first-coefficient source value budgets:',
              mp.nstr(endpoints(previous_budget)[1],18),
              mp.nstr(endpoints(new_budget)[1],18),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--panels',type=int,default=128)
    parser.add_argument('--order',type=int,default=24)
    args=parser.parse_args()
    run(args.panels,args.order)
