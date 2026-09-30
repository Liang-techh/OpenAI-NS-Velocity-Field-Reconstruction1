"""Actual seeded mixed axial-angular moment and its Z jet.

The incoming row and both pulse/end atoms use the installed continuous
providers. Terminal materialized residuals are preserved. The two linear
moments do not certify quadratic energy, pressure, stress or finite energy.
"""
from functools import lru_cache
from decimal import localcontext
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import _mp
from lei_ren_part1_paper_axial_correction import from_signed_log,signed_log


class ContinuousMixedMoments:
    def __init__(self,profile):
        self.profile=profile;self.schedule=profile.schedule
        self.precision=profile.precision
        self.angular=profile.angular_schedule_provider.primitive

    @lru_cache(maxsize=128)
    def _incoming_offset(self,z_key):
        with mp.workdps(self.precision):
            z=_mp(z_key);source=self.profile.seed_source
            value=source.inner.evaluate_x(mp.e,z)
            base=self.schedule.at_log_radius(mp.nstr(source.logRh,self.precision),z)
            theta=mp.mpf(5)/8*mp.sqrt(2)*mp.exp(mp.mpf('1.5')*source.logRh+_mp(base['log_angular_amplitude']))
            LZ=-2*z/(1+z*z)
            return (value['moments']['theta_z']-4*z*theta,
                    value['momentsZ']['theta_z']-4*theta-4*z*LZ*theta)

    @lru_cache(maxsize=128)
    def _partial_mixed_factor(self,y_key):
        with mp.workdps(self.angular.precision):
            y=_mp(y_key);p=self.angular
            if y>=p.cutoff_y:
                return p.mixed_factor(order=self.profile.incoming_order)
            if y<=0:return 4*mp.exp(p.logPstar+mp.mpf('1.6')*y)/mp.mpf('1.6')
            unit=mp.quad(lambda t:mp.exp(mp.mpf('1.6')*t-mp.mpf('.6')*p.J(t)),[0,min(y,mp.mpf(1))])
            transition=mp.mpf(0)
            if y>1:
                transition=mp.exp(mp.mpf('.3'))*mp.quad(
                    lambda t:mp.exp(t)*self.profile.incoming_provider.c(t),[1,y])
            return 4*mp.exp(p.logPstar)*(1/mp.mpf('1.6')+unit+transition)

    def _incoming(self,logR,z):
        y=_mp(str(self.profile.offset(logR,self.schedule.logRref)))
        # The full incoming atom is the same row installed in the live solve.
        if y>=self.angular.cutoff_y:
            incoming=self.profile.seed_solve(float(z))['incoming']
            factor=_mp(incoming['continuous_incoming']['reference_linear_factors']['I_theta_z_times_1plusZ2_over_Z'])
        else:
            factor=self._partial_mixed_factor(mp.nstr(y,self.precision))
        scale=mp.sqrt(2)*mp.exp(mp.mpf('1.5')*_mp(str(self.schedule.logRref)))
        offsets=self._incoming_offset(mp.nstr(z,self.precision))
        return (scale*factor*z/(1+z*z)+offsets[0],
                scale*factor*(1-z*z)/(1+z*z)**2+offsets[1])

    def normalized_row_jet(self,logR,Z):
        """N2 in the frozen Rv/Ev normalization, with its actual Z jet."""
        p=self.profile;s=self.schedule
        with mp.workdps(max(self.precision,p.jet_precision)):
            z=_mp(Z)
            start=p.offset(logR,s.logR_p)
            if start<0:raise ValueError('Normalized pulse row starts at Rp')
            end=min(_mp(str(p.offset(logR,s.logR_v))),mp.mpf(0))
            runtime=p.runtime(float(z));component=runtime.component
            tangent=p.coefficient_tangent(float(z))
            baseZ=from_signed_log(tangent['input']['base_Z'][1])
            with localcontext() as context:
                context.prec=self.precision
                xi=_mp(str(s.mu*start))
            if end>=mp.mpf('-3.15'):
                balance=component.terminal_balance(2,self.precision)
                integral=balance['pulse_integral']
                N=balance['value']-component.weighted_tail(2,end)
                NZ=(baseZ+tangent['a_Z']*integral+runtime.full_end(tangent['c_Z'],2)
                    -component.weighted_tail(2,end,coefficients=tangent['c_Z']))
                method='retained_terminal_balance_minus_owned_end_tail'
            else:
                if xi>=11:
                    integral=runtime.p[1]
                elif xi<=0:
                    integral=mp.mpf(0)
                else:
                    atom=runtime.pulse.partial_row(runtime.mu,2,xi)
                    integral=mp.exp(_mp(atom['log_normalized_pulse_integral']))
                N=runtime.base[1]+runtime.a*integral
                NZ=baseZ+tangent['a_Z']*integral
                method='owned_forward_pulse_atom'
            return dict(value=N,derivative=NZ,method=method,
                terminal_residual_retained=True,terminal_mean_forced_zero=False,
                inherited_input_uncertainty_enclosed=False,
                runtime_atoms_shared=True,input_Z_model_uses_float_backed_future_difference=True)

    def moments_jet(self,logR,Z):
        with mp.workdps(self.precision):
            z=_mp(Z);p=self.profile;s=self.schedule
            if abs(z)>=1:raise ValueError('Actual inner seed requires |Z|<1')
            boundary=_mp(mp.nstr(p.seed_source.logRh,self.precision))
            if _mp(logR)<boundary:raise ValueError('Actual mixed primitive starts at Rh')
            if p.offset(logR,s.logR_p)<0:
                mixed,mixedZ=self._incoming(logR,z)
                method='continuous_incoming_with_actual_inner_mixed_offset'
                normalized=None
            else:
                normalized=self.normalized_row_jet(logR,z)
                Ev=_mp(s.at_log_radius(s.logR_v,z)['log_angular_amplitude'])
                scale=mp.sqrt(2)*mp.exp(mp.mpf('1.5')*_mp(str(s.logR_v))+2*Ev)
                LZ=-2*z/(1+z*z)
                mixed=scale*normalized['value']
                mixedZ=scale*(normalized['derivative']+2*LZ*normalized['value'])
                method=normalized['method']
            mass=p.axial_average_jet(logR,z)
            R=mp.exp(_mp(logR))
            return dict(z=R*mass['value'],z_Z=R*mass['derivative'],
                theta_z=mixed,theta_z_Z=mixedZ,normalized_row2=normalized,
                method=method,actual_inner_seed_installed=True,
                terminal_mean_forced_zero=False,complete_five_moments=False,
                finite_energy_certified=False,scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field
    print('building shared field for mixed moments',flush=True)
    source=build_joined_field();folder=Path(__file__).parent
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    engine=field.outer.mixed_moment_provider
    with mp.workdps(field.precision):
        z=mp.mpf('.3');p=field.outer;s=source.schedule
        at=engine.moments_jet(str(s.logR_p),z)
        before=engine.moments_jet(p.log_at(s.logR_p,'-1'),z)
        jumps=[mp.nstr(abs(at[key]/before[key]-1),40) for key in ('theta_z','theta_z_Z')]
        if max(map(mp.mpf,jumps))>mp.mpf('1e-70'):raise ArithmeticError('Rp mixed matching failed')
        samples=[];errors=[]
        for label,radius in (('bulk',p.log_at(s.logR_p,mp.nstr(5/_mp(str(s.mu)),field.precision))),
                             ('end_bump',p.log_at(s.logR_v,'-2.96')),
                             ('after_Rv',p.log_at(s.logR_v,'2'))):
            row=engine.moments_jet(radius,z)
            samples.append(dict(label=label,theta_z=signed_log(row['theta_z'],field.precision),
                theta_z_Z=signed_log(row['theta_z_Z'],field.precision),method=row['method']))
            if label=='after_Rv':continue
            h=mp.mpf('1e-5');center=engine.normalized_row_jet(radius,z)
            neighbors={i:engine.normalized_row_jet(p.log_at(radius,i*h),z) for i in (-2,-1,1,2)}
            end=_mp(str(p.offset(radius,s.logR_v)))
            lam=mp.mpf('.5')-2*_mp(str(s.mu))
            base=s.at_log_radius(radius,z)
            velocity=p.values_with_jets(radius,z)
            amplitude=mp.exp(_mp(base['log_angular_amplitude']))
            LZ=-2*z/(1+z*z)
            targets=[mp.exp(lam*end)*velocity['Uz']/amplitude,
                     mp.exp(lam*end)*(velocity['Uz_Z']-LZ*velocity['Uz'])/amplitude]
            local=[]
            for key,target in zip(('value','derivative'),targets):
                difference=sum(w*neighbors[i][key] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*h)
                local.append(mp.nstr(abs(difference/target-1),40))
            if max(map(mp.mpf,local))>mp.mpf('1e-8'):raise ArithmeticError('Mixed moment integrand failed')
            errors.append(dict(label=label,relative_errors=local))
        terminal=engine.normalized_row_jet(str(s.logR_v),z)
        report=dict(Rp_mixed_relative_matching=jumps,integrand_checks=errors,samples=samples,
            terminal_row2=signed_log(terminal['value'],field.precision),
            terminal_row2_Z=signed_log(terminal['derivative'],field.precision),
            runtime_atoms_shared=True,terminal_mean_forced_zero=False,
            complete_five_moments=False,finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'Rp_mixed_relative_matching':jumps,'integrand_checks':errors}),flush=True)
    return report


if __name__=='__main__':run()
