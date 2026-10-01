"""Retained finite-core shear trace versus the prescribed exit shear.

An exactly stress-free theorem core cannot be inferred from a small finite
coefficient residual. Missing pressure/radial orders remain separate gates.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    inlet,source=read_inlet();ctx=inlet['ctx']
    with mp.workdps(ctx.dps+40):
        def value(name):return inlet[name].value.evaluate(width=0).value
        R=value('R');F=value('F');Itheta=value('I_theta');Iz=value('I_z')
        if endpoints(F)[0]<=0:raise AssertionError('Finite inlet amplitude must be positive')
        angular_shear=2*R*value('F_R');axial_shear=ctx.sqrt(2*R)*value('Uz_R')
        normalized_gap=(Itheta+angular_shear)/F
        lo,hi=endpoints(normalized_gap)
        if lo<=0:raise AssertionError('Retained angular gap no longer separated from zero; revisit diagnosis')
        # At chi=1 the bridge prescribes S=-I. Its logF_y trace therefore
        # differs from the finite core by -T_theta/(2F).
        log_slope_gap=-normalized_gap/2
        # Resolve the retained stress by pressure-parameter degree without
        # changing the datum or recomputing the source. These are coefficients
        # divided by the pressure-free inlet F, not coefficients of T/F.
        Rring=inlet['R'].value
        stress_ring=inlet['I_theta'].value+2*Rring*inlet['F_R'].value
        Fbase=inlet['F'].value.component(0,0).value
        pressure_degree_stress={str(p):stress_ring.component(p,0).value/Fbase
                               for p in range(stress_ring.pressure_order+1)}
        logh=ctx.mpf(-100)-100*ctx.mpf('1e152')
        report=dict(source=source,Z='.3',R=R,radial_degree=18,retained_pressure_order=9,
            core_angular_shear=angular_shear,prescribed_bridge_angular_shear=-Itheta,
            normalized_finite_core_angular_stress=normalized_gap,
            bridge_minus_core_logF_y=log_slope_gap,
            finite_core_axial_stress_interval=Iz+axial_shear,
            normalized_gap_excludes_zero=True,log_gap_over_collar_width=ctx.log(normalized_gap)-logh,
            retained_stress_by_pressure_degree_over_pressure_free_F=pressure_degree_stress,
            pressure_degree_breakdown_is_not_a_missing_order_bound=True,
            conditional_pressure_integral_error_included=True,
            core_to_bridge_C1_matching_certified=False,exact_core_stress_free_identity_certified=False,
            diagnosis_scope='Retained degree18 / pressure9 inlet; not a theorem counterexample',
            required_next_dependency='Enclose infinite core remainder and restore compatible stress-free trace before zero-stress edge cone proof',
            omitted_pressure_orders_enclosed=False,infinite_radial_remainder_enclosed=False,
            full_K_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('finite normalized angular stress',mp.nstr(lo,22),mp.nstr(hi,22),flush=True)


if __name__=='__main__':run()
