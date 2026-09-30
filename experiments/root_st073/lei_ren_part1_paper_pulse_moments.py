"""Five actual source moments from the reference through the end of pulse.

This combines independent incoming, angular and axial cumulative integrals.
Terminal moments are replayed and may remain nonzero; no cone is presumed.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_reference_moments import PaperReferenceMoments
from lei_ren_part1_paper_angular_cumulative import AngularCumulative
from lei_ren_part1_paper_pulse_cumulative import PaperPulseCumulative
from lei_ren_part1_paper_pulse_energy_cumulative import PulseEnergyCumulative
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


class SourcePulseMoments:
    def __init__(self,profile,*,order=128):
        self.profile=profile; self.precision=profile.precision
        self.reference=PaperReferenceMoments(profile,quadrature_order=order)
        self.angular=AngularCumulative(profile,order=order)
        self.axial=PaperPulseCumulative(profile)
        self.energy=PulseEnergyCumulative(profile)

    @lru_cache(maxsize=64)
    def axis_pressure(self,Z):
        p=self.profile
        with mp.workdps(self.precision):
            atref=from_signed_log(p.pressure_at_log_radius(p.schedule.logRref,Z)['corrected_pressure_nominal'])
            swirl=p.values(p.schedule.logRref,Z)['Utheta']
            return atref-mp.mpf('2.5')*swirl**2

    def __call__(self,logR,Z):
        p=self.profile; s=p.schedule
        if p.offset(logR,s.logRref)<=s.y_d:return self.reference(logR,Z)
        if p.offset(logR,s.logR_v)>0:
            raise ValueError('Five source moments currently stop at Rv')
        with mp.workdps(self.precision):
            angular=self.angular.moments(logR,Z)
            if p.offset(logR,s.logR_p)<0:
                incoming=p.coefficients(float(Z))['incoming']['dimensionless_integrals']
                Rref=mp.exp(mp.mpf(str(s.logRref)))
                mz=Rref*mp.mpf(incoming['I_z'])
                mixed=mp.sqrt(2)*Rref**mp.mpf('1.5')*mp.mpf(incoming['I_theta_z'])
                axialenergy=Rref*mp.mpf(incoming['I_uz2'])
            else:
                pulse=self.axial.evaluate(logR,Z)['moments']
                mz=from_signed_log(pulse['Mz']); mixed=from_signed_log(pulse['Mtheta_z'])
                axialenergy=self.energy.moment(logR,Z)
            pressure=from_signed_log(p.pressure_at_log_radius(logR,Z)['corrected_pressure_nominal'])
            return {'theta':angular['theta'],'z':mz,'theta_z':mixed,
                    'z_theta':axialenergy-angular['swirl_energy']/2,
                    'p':pressure-self.axis_pressure(float(Z))}


def run():
    from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
    from lei_ren_part1_paper_source_stress import SourceStress
    p=CorrectedSourceProfile(); provider=SourcePulseMoments(p)
    with mp.workdps(p.precision):
        radius=p.log_at(p.schedule.logR_p,5/mp.mpf(str(p.schedule.mu)))
        actual=provider(radius,.3); h=mp.mpf('.001')
        plus=provider(p.log_at(radius,h),.3); minus=provider(p.log_at(radius,-h),.3)
        R=mp.exp(mp.mpf(str(radius))); values=p.values(radius,.3)
        rhs=R*(values['Uz']**2-values['Utheta']**2/2)
        defect=abs((plus['z_theta']-minus['z_theta'])/(2*h)-rhs)/abs(rhs)
        stress=SourceStress(p,provider).receipt(radius,.3)
        terminal=provider(p.schedule.logR_v,.3)
        assert defect<mp.mpf('1e-7')
        return {'bulk_logR':str(radius),'Z':'.3',
            'bulk_moments':{k:signed_log(v,p.precision) for k,v in actual.items()},
            'bulk_quadratic_radial_derivative_relative_error':mp.nstr(defect,30),
            'bulk_actual_stress':stress,
            'terminal_moments':{k:signed_log(v,p.precision) for k,v in terminal.items()},
            'scope':'Actual cumulative moments through Rv; numerical quadrature and heat uncertainties remain.',
            'terminal_moments_forced_zero':False,'stress_cone_certified':False,
            'regular_axis_core':False,'global_finite_energy_certified':False,
            'scale_recursion_established':False}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'bulk_quadratic_radial_derivative_relative_error':result['bulk_quadratic_radial_derivative_relative_error'],
                      'five_actual_moments':list(result['bulk_moments'])}))
