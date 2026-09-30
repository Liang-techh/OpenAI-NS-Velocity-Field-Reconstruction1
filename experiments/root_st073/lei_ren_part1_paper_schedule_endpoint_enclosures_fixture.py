"""Resolved independent primitive/positive-stage integration checks."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints


def run():
    s=PaperOuterSchedule(logPstar='2',logRref='10',delta='.001',Md='.5',c_mu='.001',c_delta='.001',c_epsilon='.01')
    e=ScheduleEndpointEnclosures(s)
    with mp.workdps(e.precision+40):
        for x in ('-.1','0','.3','1','1.2'):
            xx=mp.mpf(x);v=e.primitive(e.scalar(x));lo,hi=endpoints(v)
            direct=mp.mpf(0) if xx<=0 else xx-mp.mpf('.5') if xx>=1 else mp.quad(lambda z:ContinuousAxialPulse.sigma_pair(z)[0],[0,xx])
            assert lo<=direct<=hi
        # Direct integration over the first variable-slope stage; J itself
        # uses independent adaptive quadrature, not the interval primitive.
        def ell(y):
            J=mp.quad(lambda z:ContinuousAxialPulse.sigma_pair(z)[0],[0,y])
            return mp.mpf('.1')*y-mp.mpf('.6')*J
        # The envelope is loose; nested quadrature needs only resolved
        # precision here, not hundreds of digits of endpoint arithmetic.
        with mp.workdps(25):
            actual=mp.quad(lambda y:mp.exp(2*ell(y))/2,[0,mp.mpf('.5'),1])
        upper=e.stage_pressure_upper('slope_transition_ref')['mass_upper']
        assert actual<=upper
        # Boundary-crossing interval at J(1) must remain narrow, rather than
        # revert to [0,1] when the argument encloses both sides of 1.
        narrow=e.primitive(e.iv.mpf(['.999999999999','1.000000000001']))
        lo,hi=endpoints(narrow);assert hi-lo<mp.mpf('3e-12')
        report=dict(primitive_checks_passed=True,boundary_crossing_enclosure_passed=True,
            independent_variable_stage_mass=mp.nstr(actual,25),independent_integration_dps=25,stage_mass_upper=mp.nstr(upper,60),
            parameter_scope='stored Decimal schedule values',directed_interval_arithmetic=True,
            full_source_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('resolved endpoint/mass enclosure checks passed',flush=True)


if __name__=='__main__':run()
