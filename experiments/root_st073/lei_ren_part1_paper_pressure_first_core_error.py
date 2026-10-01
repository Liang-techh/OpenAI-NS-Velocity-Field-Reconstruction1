"""Directed pressure-error propagation into the first radial core coefficient.

For fixed F0 and U0, Section 8.2 gives
Delta U1 = [(1-Z^2) Delta P0_Z - 2(1+delta) Z Delta P0]
           / [2(1-delta Z^2)].
Only value and first axial derivative of U1 can be bounded from C2 pressure.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from mpmath.ctx_iv import MPIntervalContext


def run():
    iv=MPIntervalContext();iv.dps=300
    r=json.loads(Path(__file__).with_name('lei_ren_part1_paper_uniform_preheat_error_budget.json').read_text())
    assert r['uniform_Z_error_enclosed'] and r['all_pressure_stages_included']
    with mp.workdps(330):
        bounds=[iv.mpf(mp.make_mpf(tuple(r['normalized_uniform_derivative_error_upper_bounds'][name]['exact_mpf_tuple']))) for name in ('value','first_Z','second_Z')]
        e0,e1,e2=bounds;a=iv.mpf(r['radius']);delta=iv.mpf('1e-200')
        c=2*(1+delta);Lmin=1-delta*a*a
        assert endpoints(Lmin)[0]>0
        B=e1+c*a*e0
        value=B/(2*Lmin)
        first=(e2+(2+c)*a*e1+c*e0)/(2*Lmin)+delta*a*B/(Lmin*Lmin)
        report=dict(radius=r['radius'],delta='1e-200',pressure_receipt='lei_ren_part1_paper_uniform_preheat_error_budget.json',
            normalized_U1_error_upper=endpoints(value)[1],normalized_U1_Z_error_upper=endpoints(first)[1],
            normalization='U1 and its derivative divided by Pstar squared',
            F1_pressure_error=0,P1_pressure_error=0,
            U1_second_axial_derivative_enclosed=False,higher_radial_coefficients_enclosed=False,
            axis_F0_U0_held_fixed=True,core_RK_error_enclosed=False,
            original_parameter_errors_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('normalized U1 error',mp.nstr(endpoints(value)[1],16),'U1_Z error',mp.nstr(endpoints(first)[1],16))


if __name__=='__main__':run()
