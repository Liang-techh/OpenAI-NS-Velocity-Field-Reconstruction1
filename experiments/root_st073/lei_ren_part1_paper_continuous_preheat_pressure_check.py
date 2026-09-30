"""Replay the complete analytic preheat target on the current source schedule.

No nonlinear core or completed pressure is changed by this diagnostic.
"""
import json
from pathlib import Path
import mpmath as mp

from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
from lei_ren_part1_paper_continuous_angular_schedule import install_continuous_angular_schedule


def source_profile():
    with mp.workdps(260):
        log_r = mp.nstr(mp.log(110) + 10 * (mp.mpf('5e151') + 14), 260)
        schedule = PaperOuterSchedule(logPstar='14', logRref=log_r,
            delta='1e-200', Md='.5', c_mu='.001', c_delta='.001', c_epsilon='.01')
        profile = CorrectedSourceProfile(schedule=schedule, match_waiting=True, precision=260)
        install_continuous_angular_schedule(profile.schedule,
            precision=260, primitive_precision=100)
        return profile


def signed(value):
    if value == 0:
        return dict(sign=0, log_abs=None)
    return dict(sign=int(mp.sign(value)), log_abs=mp.nstr(mp.log(abs(value)), 70))


def serialize(value):
    if isinstance(value, mp.mpf):
        return dict(signed(value), arbitrary_exponent_value=mp.nstr(value, 80))
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value


def independent_centered_coefficients(beta, center, degree):
    """Direct mpmath differentiation, independent of adapter recurrence."""
    return mp.taylor(lambda z: (1 + z*z)**(-beta), center, degree)


def run():
    from lei_ren_part1_paper_continuous_preheat_pressure import (
        ContinuousPreheatPressure, _q_power_taylor)
    print('constructing source schedule and complete preheat datum', flush=True)
    profile = source_profile()
    adapter = ContinuousPreheatPressure(profile, quadrature_order=192)
    with mp.workdps(260):
        errors=[]
        for beta in ('0', '.001', '.7', '2'):
            for center in ('0', '.3', '1'):
                got=_q_power_taylor(mp.mpf(beta),mp.mpf(center),12)
                direct=independent_centered_coefficients(mp.mpf(beta),mp.mpf(center),12)
                err=max(abs(a-b) for a,b in zip(got,direct))
                if err>mp.mpf('1e-230'):
                    raise ArithmeticError('Centered Taylor recurrence disagrees with direct differentiation')
                errors.append(mp.nstr(err,25))
        fine=adapter.pressure_jet('.3',degree=8)
        coarse=adapter.pressure_jet('.3',degree=8,quadrature_order=128)
        tail=fine['post_Rv_pressure_coefficients']
        h=mp.mpf('1e-30')
        plus=adapter.post_Rv_components_jet(0,center=mp.mpf('.3')+h)['pressure_value']
        minus=adapter.post_Rv_components_jet(0,center=mp.mpf('.3')-h)['pressure_value']
        derivative=(plus-minus)/(2*h)
        rel=abs((derivative-tail[1])/tail[1])
        if rel>mp.mpf('1e-50'):
            raise ArithmeticError('Actual-source post-Rv derivative check failed')
        stage_checks=[]
        for name,row in fine['post_Rv_components'].items():
            old=coarse['post_Rv_components'][name]
            change=abs((row['pressure_value']-old['pressure_value'])/row['pressure_value'])
            if row['Z_independent'] and any(row['pressure_coefficients'][1:]):
                raise ArithmeticError('Post-flatten preheat stage has a spurious Z derivative')
            stage_checks.append(dict(stage=name,pressure_Taylor_coefficients=[signed(v) for v in row['pressure_coefficients']],pressure=signed(row['pressure_value']),
                pressure_Z=signed(row['derivative']),order128_to192_relative_change=mp.nstr(change,40)))
        prefix=fine['dominant_pressure_coefficients'][0]
        loss=(fine['pressure_value']-prefix)==0 and tail[0]!=0
        report=dict(source_version='2609.35406v2',Z='.3',schedule_inputs=profile.schedule.metadata()['inputs'],
            normalized_prefix=signed(prefix),post_Rv_pressure=signed(tail[0]),
            post_Rv_pressure_Z=signed(tail[1]),post_Rv_to_prefix_log_abs=mp.nstr(mp.log(abs(tail[0]/prefix)),70),
            tail_lost_by_subtracting_nominal_total=bool(loss),stages=stage_checks,
            centered_Taylor_direct_differentiation_errors=errors,
            actual_post_Rv_derivative_relative_error=mp.nstr(rel,40),
            all_post_flatten_Z_derivatives_zero=True,raw_heat_pressure_used=False,
            core_rebuilt_with_complete_target=False,actual_angular_bumps_restored=False,
            quadrature_error_enclosed=False,pressure_terminal_compatibility_certified=False,
            finite_energy_certified=False,stress_cone_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('post_Rv_pressure','post_Rv_pressure_Z',
        'tail_lost_by_subtracting_nominal_total','actual_post_Rv_derivative_relative_error')},indent=2),flush=True)
    return report


if __name__ == '__main__':
    run()
