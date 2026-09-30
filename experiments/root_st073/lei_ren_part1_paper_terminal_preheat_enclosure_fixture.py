"""Independent resolved H=1 collar integration and stored normalization audit."""
import json
from pathlib import Path
from types import SimpleNamespace
from decimal import Decimal,localcontext
import mpmath as mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure


def run():
    # Resolved schedule data with the same normalization formula, held fixed
    # as decimal strings. A deliberate later shift must not be silently reset.
    with mp.workdps(120):
        logtail=mp.mpf(20);dt=mp.mpf('.001');eps=mp.mpf('.1');logp=mp.mpf(2)
        mu=mp.mpf('.01');yd=mp.mpf(12);yr=mp.mpf(100);Ts=mp.mpf(8);yt=mp.mpf(110)
        J=lambda x:mp.mpf(0) if x<=0 else x-mp.mpf('.5') if x>=1 else (_ for _ in ()).throw(ValueError('partial fixture'))
        ell=mp.mpf('.1')*yt-mp.mpf('.6')*J(yt)-mu*J(yt-yd)-(1-mu)*J(yt-yr)+(1-dt/2)*J(yt-yr-1-Ts)
        logc=logp+ell+(1+dt)*logtail/2-mp.log(2*(1-eps))
        dec=lambda x:Decimal(mp.nstr(x,110))
        s=SimpleNamespace(decimal_precision=110,delta=dec(dt),epsilon=dec(eps),logPstar=dec(logp),
            _log_c_inf=dec(logc),logR_tail=dec(logtail),mu=dec(mu),y_d=dec(yd),y_rel=dec(yr),Ts=dec(Ts),y_tail=dec(yt),y_b=dec(yt+3))
        e=ScheduleEndpointEnclosures(s);r=e.terminal_preheat_bounds()
        scale=mp.exp(2*(mp.mpf(str(s._log_c_inf))-logp)-mp.log(2)-(1+dt)*logtail)
        collar=scale*mp.quad(lambda t:mp.exp(-(1+dt)*t)*ContinuousPreheatPressure._heat_collar_factor(t,eps)**2,[0,1,2,3])
        lo,hi=endpoints(r['heat_collar_mass_interval']);assert lo<=collar<=hi
        tail=scale*mp.exp(-3*(1+dt))/(1+dt)
        lo,hi=endpoints(r['exterior_mass_interval'])
        # Independent scalar rounding is less precise than the interval;
        # compare to its midpoint with an explicit numerical tolerance.
        assert abs((lo+hi)/2-tail)/tail<mp.mpf('1e-110')
        old=endpoints(r['stored_heat_scale'])[0]
        with localcontext() as context:
            context.prec=110
            s._log_c_inf+=Decimal('.01')
        shifted=e.terminal_preheat_bounds()
        ratio=endpoints(shifted['stored_heat_scale'])[0]/old
        assert abs(ratio-mp.exp(mp.mpf('.02')))<mp.mpf('1e-110')
        report=dict(independent_collar_inside_interval=True,exterior_formula_agrees=True,
            stored_normalization_shift_preserved=True,heat_scale_shift_ratio=mp.nstr(ratio,60),
            quadrature_error_used_as_certificate=False,full_heat_velocity_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('terminal preheat integral enclosure checks passed',flush=True)


if __name__=='__main__':run()
