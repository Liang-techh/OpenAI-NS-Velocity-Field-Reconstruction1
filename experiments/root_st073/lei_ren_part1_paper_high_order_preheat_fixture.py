"""Independent scalar switch derivatives versus interval Taylor coefficients."""
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_high_order_preheat_integrals import switch_taylor,HighOrderPreheatIntegrals


def run():
    e=ScheduleEndpointEnclosures(SimpleNamespace(decimal_precision=100));iv=e.iv
    with mp.workdps(120):
        errors=[]
        for x in (mp.mpf('.3'),mp.mpf('.5'),mp.mpf('.7')):
            jet=switch_taylor(iv,iv.mpf(x),8)
            cell=switch_taylor(iv,iv.mpf([x-mp.mpf('.001'),x+mp.mpf('.001')]),8)
            for n in range(9):
                direct=mp.diff(lambda z:ContinuousAxialPulse.sigma_pair(z)[0],x,n)/mp.factorial(n)
                lo,hi=endpoints(jet.coefficients[n]);middle=(lo+hi)/2
                error=abs(middle-direct)/max(1,abs(direct));assert error<mp.mpf('1e-100')
                cl,ch=endpoints(cell.coefficients[n]);assert cl<=direct<=ch
                errors.append(error)
        integrator=HighOrderPreheatIntegrals(e)
        grid=integrator.primitive_grid(64,12)
        point_errors=[]
        for k in (16,32,48):
            x=mp.mpf(k)/64
            direct=mp.quad(lambda z:ContinuousAxialPulse.sigma_pair(z)[0],[0,x])
            lo,hi=endpoints(grid[k]);assert lo<=direct<=hi
            point_errors.append(mp.nstr(hi-lo,40))
        report=dict(independent_switch_derivative_checks=len(errors),
            maximum_scaled_discrepancy=mp.nstr(max(errors),40),
            primitive_grid_independent_integral_checks=3,primitive_interval_widths=point_errors,
            source_certificate=False,directed_interval_arithmetic=True)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('high-order switch derivative and primitive checks passed',flush=True)


if __name__=='__main__':run()
