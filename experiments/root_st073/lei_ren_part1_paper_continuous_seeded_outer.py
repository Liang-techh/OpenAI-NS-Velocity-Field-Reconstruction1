"""Install continuous axial corrections into the retained joined exterior.

Both Uz and its cumulative mass use one component. Incoming primitives,
float-backed Z derivatives and integral enclosures remain uncertified.
"""
from decimal import Decimal, localcontext
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_solve import (
    solve_continuous, ContinuousAxialCorrection)
from lei_ren_part1_paper_seeded_outer_field import SeededAxialProfile, SeededOuterField
from lei_ren_part1_paper_joined_outer import build_joined_field, _mp, _finite_difference


class ContinuousSeededAxialProfile(SeededAxialProfile):
    def __init__(self,source,*,prepared=None):
        super().__init__(source)
        self.prepared=prepared or {}

    @lru_cache(maxsize=64)
    def seed_solve(self,Z):
        z=float(Z)
        if z in self.prepared:
            seeded,continuous=self.prepared[z]
        else:
            seeded=super().seed_solve(z)
            continuous=solve_continuous(seeded,precision=self.precision)
        if _mp(continuous['input_mu'])!=_mp(str(self.schedule.mu)):
            # Serialized receipts may round only beyond their declared digits.
            with mp.workdps(continuous['algebra_precision']):
                difference=abs(_mp(continuous['input_mu'])/_mp(str(self.schedule.mu))-1)
                if difference>mp.mpf(10)**(-continuous['algebra_precision']+5):
                    raise ValueError('Prepared continuous solve uses a different schedule mu')
        return dict(incoming=seeded['incoming'],axial=continuous,
            energy_target=seeded.get('energy_target',seeded['axial'].get('energy_target')),
            seed_origin=seeded.get('seed_origin'),
            continuous_mean_closure_certified=False,finite_energy_certified=False)

    @lru_cache(maxsize=64)
    def coefficients(self,Z):
        result=self.seed_solve(float(Z))
        return dict(Z=float(Z),incoming=result['incoming'],axial=result['axial'],
                    angular=self.angular.coefficients(float(Z)))

    @lru_cache(maxsize=64)
    def component(self,Z):
        data=self.coefficients(float(Z))
        base=[data['incoming']['row_normalization'][key]
              for key in ('scaled_base_m1','scaled_base_m2')]
        return ContinuousAxialCorrection(data['axial'],base)

    def values(self,logR,Z):
        base=self.schedule.at_log_radius(logR,Z)
        with mp.workdps(self.precision):
            amplitude=mp.exp(_mp(str(base['log_angular_amplitude'])))
            angular_offset=self.offset(logR,self.schedule.logR_rel)
            correction=mp.mpf(0)
            if Decimal('-3.15')<=angular_offset<=Decimal('-.85'):
                correction=from_signed_log(self.angular.relative_bump(
                    float(angular_offset),Z,coefficients=self.coefficients(float(Z))['angular']))
            Uz=_mp(str(base['Uz']))
            start=self.offset(logR,self.schedule.logR_p)
            end=self.offset(logR,self.schedule.logR_v)
            if start>=0 and end<=0:
                component=self.component(float(Z))
                with localcontext() as ctx:
                    ctx.prec=max(self.schedule.decimal_precision,self.precision)
                    xi=self.schedule.mu*start
                if Decimal('-3.15')<=end<=Decimal('-.85'):
                    factor=component.value_jet(str(end))['beta']
                else:
                    factor=component.pulse_value_jet(str(xi))['value']
                Uz=amplitude*factor
            return dict(Utheta=amplitude*(1+correction),Uz=Uz,
                angular_relative_correction=correction,
                base_log_Utheta=str(base['log_angular_amplitude']),
                angular_log1p_correction=mp.nstr(mp.log1p(correction),self.precision),
                continuous_axial_values_installed=True)

    def axial_average_receipt(self,logR,Z):
        with mp.workdps(self.precision):
            start=self.offset(logR,self.schedule.logR_p)
            if start<0:
                value=super().axial_average(logR,Z)
                return dict(nominal=signed_log(value,self.precision),region='inherited_incoming',
                            incoming_uncertainty_enclosed=False)
            end=self.offset(logR,self.schedule.logR_v)
            component=self.component(float(Z))
            if end<=0:
                with localcontext() as ctx:
                    ctx.prec=max(self.schedule.decimal_precision,self.precision)
                    xi=self.schedule.mu*start
                atom=component.cumulative_row(1,xi=str(xi)) if xi<11 else component.cumulative_row(1,end_offset=str(end))
                logE=_mp(str(self.schedule.at_log_radius(logR,Z)['log_angular_amplitude']))
                lam=mp.mpf('.5')-component.mu
                logscale=logE-lam*_mp(str(end))
            else:
                atom=component.cumulative_row(1,end_offset=0)
                logEv=_mp(str(self.schedule.at_log_radius(self.schedule.logR_v,Z)['log_angular_amplitude']))
                # Past Rv cumulative mass is constant: M/R decays as exp(-end),
                # independently of any later heat amplitude evolution.
                logscale=logEv-_mp(str(end))
            nominal=from_signed_log(atom['nominal'])*mp.exp(logscale)
            logbound=atom.get('log_pulse_omitted_absolute_bound')
            return dict(nominal=signed_log(nominal,self.precision),
                log_pulse_omitted_absolute_bound=(mp.nstr(_mp(logbound)+logscale,self.precision)
                                                if logbound is not None else None),
                pulse_omitted_correction_sign=atom.get('pulse_omitted_correction_sign'),
                incoming_uncertainty_enclosed=False,quadrature_enclosure_certified=False,
                terminal_mean_forced_zero=False,region='continuous_seeded_exterior')

    def axial_average(self,logR,Z):
        with mp.workdps(self.precision):
            return from_signed_log(self.axial_average_receipt(logR,Z)['nominal'])


class ContinuousSeededOuterField(SeededOuterField):
    def __init__(self,source,*,prepared=None):
        super().__init__(source)
        self.outer=ContinuousSeededAxialProfile(source,prepared=prepared)

    def _outer_state(self,*args,**kwargs):
        result=super()._outer_state(*args,**kwargs)
        result.update(region='continuous_seeded_axial_outer',
                      continuous_axial_values_and_means_installed=True,
                      continuous_mean_closure_certified=False)
        return result


def run():
    print('building retained joined candidate',flush=True)
    source=build_joined_field()
    folder=Path(__file__).parent
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    continuous=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousSeededOuterField(source,prepared={.3:(seeded,continuous)})
    profile=field.outer
    with mp.workdps(field.precision):
        z=mp.mpf('.3');schedule=field.schedule;mu=_mp(str(schedule.mu))
        points=[('Rp',_mp(str(schedule.logR_p))),
            ('startup',_mp(str(schedule.logR_p))+mp.mpf('.015')/mu),
            ('pulse',_mp(str(schedule.logR_p))+5/mu),
            ('cutoff',_mp(str(schedule.logR_p))+mp.mpf('10.75')/mu),
            ('first_end_bump',_mp(str(schedule.logR_v))-3),
            ('Rv',_mp(str(schedule.logR_v))),
            ('after_Rv',_mp(str(schedule.logR_v))+1)]
        rows=[]
        for name,logR in points:
            print('continuous profile '+name,flush=True)
            values=profile.values(mp.nstr(logR,field.precision),z)
            average=profile.axial_average_receipt(mp.nstr(logR,field.precision),z)
            rows.append(dict(name=name,Uz=signed_log(values['Uz'],field.precision),
                             axial_average=average))
        Rv=_mp(str(schedule.logR_v))
        a=profile.axial_average(mp.nstr(Rv,field.precision),z)
        b=profile.axial_average(mp.nstr(Rv+1,field.precision),z)
        relative=abs(b*mp.e/a-1) if a else None
        print('recovering radial velocity from continuous seeded mean',flush=True)
        logPulse=next(radius for name,radius in points if name=='pulse')
        center=field._outer_state(logPulse,z)
        divergence=[]
        for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
            states={i:field._outer_state(logPulse+i*h,z) for i in (-2,-1,1,2)}
            radial=(states[-2]['Ur']-8*states[-1]['Ur']+8*states[1]['Ur']-states[2]['Ur'])/(12*h)
            axial_y=(states[-2]['Uz']-8*states[-1]['Uz']+8*states[1]['Uz']-states[2]['Uz'])/(12*h)
            axial_Z=_finite_difference(lambda zz:profile.values(
                mp.nstr(logPulse,field.precision),float(zz))['Uz'],z,step=field.derivative_step)
            R=center['R'];root=mp.sqrt(2*R);delta=field.delta
            terms=[root*radial/R,center['Ur']/root,
                ((1-z*z)*axial_Z-2*z*axial_y-(1+delta)*z*center['Uz'])/(1-delta*z*z)]
            divergence.append(dict(log_radius_step=mp.nstr(h,20),
                relative_cancellation=mp.nstr(abs(sum(terms))/sum(abs(v) for v in terms),40),
                mapped_q_divergence=signed_log(sum(terms),field.precision)))
        report=dict(Z='.3',rows=rows,
            source_schedule_identity_preserved=profile.schedule is source.outer.schedule,
            source_angular_identity_preserved=profile.angular is source.outer.angular,
            source_tail_identity_preserved=profile.tail is source.outer.tail,
            source_pressure_identity_preserved=profile.pressure is source.outer.pressure,
            post_Rv_constant_mass_relative_replay=mp.nstr(relative,40) if relative is not None else None,
            continuous_axial_values_and_means_installed=True,
            physical_radial_velocity_sampled=True,
            recovered_Ur=signed_log(center['Ur'],field.precision),
            pulse_divergence_diagnostics=divergence,
            scope='Installed joined adapter with point/mean and independent pulse radial divergence diagnostics; float-backed Z derivatives, incoming integration and global energy remain uncertified.',
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('post_Rv_constant_mass_relative_replay',
        'source_schedule_identity_preserved','continuous_axial_values_and_means_installed')}),flush=True)
    return report


if __name__=='__main__':run()
