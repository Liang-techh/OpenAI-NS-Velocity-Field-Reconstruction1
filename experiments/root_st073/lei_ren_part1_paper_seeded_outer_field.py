"""Install the actual Rh-seeded axial solve in velocity AND mass transport.

The shared angular schedule/pressure are retained. A terminal mean is replayed
from the solved pulse and end bumps, never assigned zero. Finite quadrature,
float-backed Z input and inherited source uncertainties remain explicit.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
from lei_ren_part1_paper_joined_outer import (
    JoinedOuterField, build_joined_field, _finite_difference, _signed_log, _mp)
from lei_ren_part1_paper_seeded_axial_inputs import solve_joined_seed


class SeededAxialProfile(CorrectedSourceProfile):
    def __init__(self, source):
        self.seed_source = source
        original = source.outer
        super().__init__(precision=source.precision,order=original.order,
            schedule=original.schedule,angular=original.angular,tail=original.tail)
        # Axial changes do not change swirl's radial pressure integral.
        self.pressure = original.pressure
        self.pulse = original.pulse
        self.original = original

    @lru_cache(maxsize=64)
    def seed_solve(self,Z):
        return solve_joined_seed(self.seed_source,str(Z))

    @lru_cache(maxsize=64)
    def coefficients(self,Z):
        solved = self.seed_solve(float(Z))
        prior = self.original.coefficients(float(Z))
        return dict(Z=float(Z),incoming=solved['incoming'],
            axial=solved['axial'],angular=prior['angular'])

    def axial_average(self,logR,Z):
        with mp.workdps(self.precision):
            start = self.offset(logR,self.schedule.logR_p)
            if start < 0:
                # The earlier Uz is unchanged; only its incoming primitive
                # includes the actual corrected inner mass.
                offsets = self.seed_source.terminal_offsets(
                    self.seed_source._z_key(_mp(Z)))
                return self.original.axial_average(logR,Z)+offsets['mass_offset']*mp.exp(-_mp(logR))
            # Updated m1 already contains the inner mass here. Replaying the
            # pulse and end bumps must not add that mass a second time.
            return super().axial_average(logR,Z)


class SeededOuterField(JoinedOuterField):
    def __init__(self,source):
        self.original_join = source
        super().__init__(source.inner,derivative_step=source.derivative_step)
        self.outer = SeededAxialProfile(source)
        self.axial_mean_already_seeded = True

    def _outer_average(self,log_radius,z):
        with mp.workdps(self.precision):
            logR = mp.nstr(log_radius,self.precision)
            # Supply the actual seeded mean directly. Subtracting the inner
            # offset and adding it back can erase a smaller exterior mean.
            offset = self.terminal_offsets(self._z_key(z))
            inverse_R = mp.exp(-log_radius)
            if self.outer.offset(logR,self.schedule.logR_p)<0:
                average,derivative=self.original_join._outer_average(log_radius,z)
                return (average+offset['mass_offset']*inverse_R,
                    derivative+offset['mass_offset_Z']*inverse_R)
            average = self.outer.axial_average(logR,z)
            derivative = _finite_difference(
                lambda zz:self.outer.axial_average(logR,float(zz)),z,
                step=self.derivative_step)
            return average,derivative

    def _outer_state(self,log_radius,z,*,carry=True):
        if not carry:
            raise ValueError('Seeded mass already includes the inner offset; uncarried mode is undefined')
        result = super()._outer_state(log_radius,z,carry=True)
        result['region']='seeded_axial_outer_profile'
        result['seeded_axial_coefficients_installed']=True
        result['continuous_mean_closure_certified']=False
        return result


def run():
    print('building corrected shared candidate',flush=True)
    original = build_joined_field()
    field = SeededOuterField(original)
    with mp.workdps(field.precision):
        z=mp.mpf('.3'); schedule=field.schedule
        # Keep stage coordinates in full MP precision (ordinary float
        # offsets cannot resolve these source radii).
        points=[('Rp',_mp(str(schedule.logR_p))),
            ('pulse',_mp(str(schedule.logR_p))+5/_mp(str(schedule.mu))),
            ('first_end_bump',_mp(str(schedule.logR_v))-3),
            ('Rv',_mp(str(schedule.logR_v))),
            ('heat',_mp(str(schedule.logR_b))+1)]
        rows=[]
        for name,logR in points:
            print('replaying '+name,flush=True)
            state=field._outer_state(logR,z)
            R=state['R']
            rows.append(dict(name=name,
                mass_over_Rh=_signed_log(state['mass']/field.Rh),
                mass_Z_over_Rh=_signed_log(state['mass_Z']/field.Rh),
                Uz=_signed_log(state['Uz']),Ur=_signed_log(state['Ur']),
                seeded_axial_coefficients_installed=True))
        solved=field.outer.seed_solve(.3)
        Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').write_text(
            json.dumps(solved,indent=2)+'\n',encoding='utf-8')
        logRp=_mp(str(schedule.logR_p))
        at_Rp=field._outer_state(logRp,z)
        before_Rp=original._outer_state(logRp,z)
        rp_mean_relative_jump=abs((at_Rp['mass']-before_Rp['mass'])/before_Rp['mass'])
        logRv=_mp(str(schedule.logR_v))
        # Beyond Rv, Uz vanishes and the cumulative mean must stay constant.
        rv=field._outer_state(logRv,z)
        after=field._outer_state(logRv+1,z)
        rv_mean_jump_over_Rh=(rv['mass']-after['mass'])/field.Rh
        # Independent radial differences of Ur and Uz plus axial differences
        # of the ACTUAL seeded velocity, rather than uncorrected schedule jets.
        logPulse=points[1][1]; center=field._outer_state(logPulse,z)
        divergence=[]
        for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
            states={i:field._outer_state(logPulse+i*h,z) for i in (-2,-1,1,2)}
            radial=(states[-2]['Ur']-8*states[-1]['Ur']+8*states[1]['Ur']-states[2]['Ur'])/(12*h)
            axial_y=(states[-2]['Uz']-8*states[-1]['Uz']+8*states[1]['Uz']-states[2]['Uz'])/(12*h)
            axial_Z=_finite_difference(lambda zz:field.outer.values(
                mp.nstr(logPulse,field.precision),float(zz))['Uz'],z,step=field.derivative_step)
            R=center['R']; root=mp.sqrt(2*R); dt=field.delta
            terms=[root*radial/R,center['Ur']/root,
                ((1-z*z)*axial_Z-2*z*axial_y-(1+dt)*z*center['Uz'])/(1-dt*z*z)]
            divergence.append(dict(log_radius_step=mp.nstr(h,20),
                mapped_q_divergence=_signed_log(sum(terms)),
                relative_cancellation=mp.nstr(abs(sum(terms))/sum(abs(v) for v in terms),40)))
        old=original._tail_state(z); new=field._tail_state(z)
        old_coefficient=old['tail_radial_transport_coefficient']
        new_coefficient=new['tail_radial_transport_coefficient']
        report=dict(Z='.3',precision=field.precision,rows=rows,
            canonical_linear_replay=solved['canonical_quadrature_replay']['row_relative_differences'],
            energy_equation_replay=solved['axial']['energy_relative_replay'],
            Rp_mean_relative_jump=mp.nstr(rp_mean_relative_jump,40),
            Rv_mean_jump_over_Rh=_signed_log(rv_mean_jump_over_Rh),
            pulse_independent_divergence=divergence,
            old_tail_coefficient_over_Rh=_signed_log(old_coefficient/field.Rh),
            new_tail_coefficient_over_Rh=_signed_log(new_coefficient/field.Rh),
            tail_coefficient_ratio=_signed_log(new_coefficient/old_coefficient),
            terminal_mean_was_assigned_zero=False,
            source_schedule_identity=field.outer.schedule is original.outer.schedule,
            source_angular_identity=field.outer.angular is original.outer.angular,
            source_pressure_identity=field.outer.pressure is original.outer.pressure,
            seeded_axial_coefficients_installed=True,
            continuous_mean_closure_certified=False,finite_energy_certified=False,
            outer_velocity_jets_complete=False,
            scope='Actual numerical velocity/primitive installation; sampled float-Z derivative and inherited quadrature uncertainty remain.')
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({k:report[k] for k in ('canonical_linear_replay','energy_equation_replay','tail_coefficient_ratio')}),flush=True)
        return report


if __name__=='__main__':run()
