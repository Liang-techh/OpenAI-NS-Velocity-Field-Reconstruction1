"""Finite heat cumulative moments with separate tiny corrections.

Uses the same quadratic small-x heat polynomial installed in point jets.
Collar atoms use declared Gauss nodes; exterior polynomial atoms integrate
analytically. Exact heat/input uncertainty and infinite exterior closure are
not certified by these finite-radius cumulative moments.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import _mp
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_axial_correction import signed_log


class ContinuousHeatMoments:
    def __init__(self,angular):
        self.angular=angular;self.profile=angular.profile
        self.schedule=angular.schedule;self.precision=angular.precision

    @staticmethod
    def atom(rate,left,right):
        if right<=left:return mp.mpf(0)
        return mp.exp(rate*left)*mp.expm1(rate*(right-left))/rate if rate else right-left

    @lru_cache(maxsize=128)
    def _collar(self,t_key,z_key):
        with mp.workdps(self.precision):
            t=_mp(t_key);z=_mp(z_key);s=self.schedule
            heat=self.profile.angular_schedule_provider.heat_jet(s.logR_tail,z)
            x=heat['xi'];xz=heat['xi_Z'];h=_mp(s.delta)/2
            if x>mp.mpf('1e-8'):
                raise ValueError('Heat moments require installed small-x polynomial branch throughout collar/exterior')
            c1=h*(1+h);c2=h*(h+1)**2*(h+2)
            rates=[1-h,-2*h];result=[mp.mpf(0)]*6
            eps=_mp(s.epsilon)
            for unit,weight in self.angular.nodes:
                v=t*unit;sw=ContinuousAxialPulse.sigma_pair(v)[0]
                edge=(3-v)/2
                phi=mp.exp(-1/edge**2) if edge>0 else mp.mpf(0)
                factor=1-eps*phi;K0=(1-sw)*(1-eps)+sw*factor
                xx=x*mp.exp(-v);xxz=xz*mp.exp(-v)
                dk=sw*factor*(-c1*xx+c2*xx*xx/2)
                dkz=sw*factor*(-c1+c2*xx)*xxz
                vals=[K0,K0*K0,dk,2*K0*dk+dk*dk,dkz,2*(K0+dk)*dkz]
                for i,val in enumerate(vals):result[i]+=t*weight*mp.exp(rates[i%2]*v)*val
            return tuple(result)

    def increments(self,t,Z):
        with mp.workdps(self.precision):
            t=_mp(t);z=_mp(Z);s=self.schedule
            if t<0:raise ValueError('Heat starts at Rtail')
            row=list(self._collar(mp.nstr(min(t,mp.mpf(3)),self.precision),mp.nstr(z,self.precision)))
            if t>3:
                heat=self.profile.angular_schedule_provider.heat_jet(s.logR_tail,z)
                x=heat['xi'];xz=heat['xi_Z'];h=_mp(s.delta)/2
                c1=h*(1+h);c2=h*(h+1)**2*(h+2)
                b1=-c1*x;b2=c2*x*x/2;b1z=-c1*xz;b2z=c2*x*xz
                rates=[1-h,-2*h]
                for i,rate in enumerate(rates):row[i]+=self.atom(rate,mp.mpf(3),t)
                for j,b,bz in ((1,b1,b1z),(2,b2,b2z)):
                    atom=self.atom(rates[0]-j,mp.mpf(3),t)
                    row[2]+=b*atom;row[4]+=bz*atom
                coeffs=[2*b1,2*b2+b1*b1,2*b1*b2,b2*b2]
                coeffz=[2*b1z,2*b2z+2*b1*b1z,2*(b1z*b2+b1*b2z),2*b2*b2z]
                for j,(b,bz) in enumerate(zip(coeffs,coeffz),1):
                    atom=self.atom(rates[1]-j,mp.mpf(3),t)
                    row[3]+=b*atom;row[5]+=bz*atom
            return row

    def moments_jet(self,logR,Z,*,include_inner=True):
        with mp.workdps(self.precision):
            p=self.profile;s=self.schedule;z=_mp(Z)
            t=_mp(str(p.offset(logR,s.logR_tail)))
            if abs(z)>=1:raise ValueError('Actual inner seed requires |Z|<1')
            anchor=self.angular.moments_jet(str(s.logR_tail),z,include_inner=include_inner)
            inc=self.increments(t,z)
            rlog=_mp(str(s.logR_tail));a=(1+_mp(s.delta))/2
            E=_mp(s._log_c_inf)-a*rlog
            scales=[mp.sqrt(2)*mp.exp(mp.mpf('1.5')*rlog+E),mp.exp(rlog+2*E)]
            base=[scales[i]*inc[i] for i in range(2)]
            correction=[scales[i%2]*inc[i+2] for i in range(4)]
            return dict(theta=anchor['theta']+base[0]+correction[0],
                swirl_energy=anchor['swirl_energy']+base[1]+correction[1],
                theta_Z=anchor['theta_Z']+correction[2],
                swirl_energy_Z=anchor['swirl_energy_Z']+correction[3],
                heat_reference_increments=base,heat_correction_increments=correction,
                heat_normalized_increments=inc,heat_anchor=anchor,
                inner_offsets_reapplied_once=include_inner,heat_supported=True,
                heat_kernel_shared=True,heat_polynomial_truncation_enclosed=False,
                quadrature_order=self.angular.order,quadrature_error_enclosed=False,
                complete_five_moments=False,infinite_heat_tail_integrated=False,
                finite_energy_certified=False,scale_recursion_certified=False)

    def complete_swirl_heat_integral(self,Z):
        """Integral Rtail..infinity Utheta^2 dR, split before summation.

        The declared delta>0 makes the reference exterior integrable. This
        does not cover radial velocity energy or the full physical Z weight.
        """
        with mp.workdps(self.precision):
            s=self.schedule;z=_mp(Z);delta=_mp(s.delta);h=delta/2
            if delta<=0 or h>mp.mpf('.5'):
                raise ValueError('Tail formula requires 0<delta<=1')
            collar=self._collar('3',mp.nstr(z,self.precision))
            heat=self.profile.angular_schedule_provider.heat_jet(s.logR_tail,z)
            x=heat['xi'];xz=heat['xi_Z']
            c1=h*(1+h);c2=h*(h+1)**2*(h+2);c3=h*(h+1)**2*(h+2)**2*(h+3)
            b1=-c1*x;b2=c2*x*x/2;b1z=-c1*xz;b2z=c2*x*xz
            coeff=[2*b1,2*b2+b1*b1,2*b1*b2,b2*b2]
            coeffZ=[2*b1z,2*b2z+2*b1*b1z,2*(b1z*b2+b1*b2z),2*b2*b2z]
            reference=collar[1]+mp.exp(-3*delta)/delta
            correction=collar[3];correctionZ=collar[5]
            for j,(b,bz) in enumerate(zip(coeff,coeffZ),1):
                atom=mp.exp(-3*(delta+j))/(delta+j)
                correction+=b*atom;correctionZ+=bz*atom
            rlog=_mp(s.logR_tail);E=_mp(s._log_c_inf)-(1+delta)*rlog/2
            scale=mp.exp(rlog+2*E)
            # H and the small-x polynomial lie in (0,1], so the square
            # discrepancy is <= 2*|H-P| <= c3*x^3/3 on the whole heat region.
            truncation_bound=scale*c3*x**3/(3*(delta+3))
            return dict(reference=scale*reference,heat_correction=scale*correction,
                heat_correction_Z=scale*correctionZ,
                swirl_heat_integral=scale*(reference+correction),
                analytic_heat_truncation_absolute_bound=truncation_bound,
                infinite_swirl_heat_tail_integrated=True,
                swirl_heat_integrable_for_declared_delta=True,
                physical_Z_energy_weight_integrated=False,radial_energy_included=False,
                quadrature_error_enclosed=False,arithmetic_error_enclosed=False,
                finite_energy_certified=False,scale_recursion_certified=False)


def run(*,tail_only=False):
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field
    folder=Path(__file__).parent
    print('building shared field for cumulative heat moments',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;engine=p.angular_moment_provider.heat_provider;reports=[]
    if tail_only:
        reports=json.loads(Path(__file__).with_suffix('.json').read_text())['samples']
    with mp.workdps(field.precision):
        z=mp.mpf('.3');step=mp.mpf('1e-4');s=p.schedule
        for t in (() if tail_only else (mp.mpf('.6'),mp.mpf(4),mp.mpf(30))):
            inc=engine.increments(t,z)
            neighbors={i:engine.increments(t+i*step,z) for i in (-2,-1,1,2)}
            actual=s.at_log_radius(p.log_at(s.logR_tail,t),z)
            sw=ContinuousAxialPulse.sigma_pair(t)[0];eps=_mp(s.epsilon)
            edge=(3-t)/2;phi=mp.exp(-1/edge**2) if edge>0 else mp.mpf(0)
            factor=1-eps*phi;K0=(1-sw)*(1-eps)+sw*factor
            heat=p.angular_schedule_provider.heat_jet(p.log_at(s.logR_tail,t),z)
            dk=-sw*factor*heat['deficit'];dkz=sw*factor*heat['H_prime']*heat['xi_Z']
            vals=[K0,K0*K0,dk,2*K0*dk+dk*dk,dkz,2*(K0+dk)*dkz]
            rates=[1-_mp(s.delta)/2,-_mp(s.delta)]
            errors=[]
            for k,val in enumerate(vals):
                derivative=sum(w*neighbors[i][k] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
                target=mp.exp(rates[k%2]*t)*val
                errors.append(abs(derivative/target-1))
            if max(errors)>mp.mpf('1e-8'):raise ArithmeticError('Heat cumulative integrand mismatch')
            full=p.angular_moments_jet(p.log_at(s.logR_tail,t),z)
            quadratic=p.quadratic_moments_jet(p.log_at(s.logR_tail,t),z)
            reports.append(dict(t=mp.nstr(t,20),integrand_relative_errors=[mp.nstr(v,40) for v in errors],
                separate_heat_correction_retained=all(v!=0 for v in full['heat_correction_increments']),
                quadratic_heat_dispatch_finite=mp.isfinite(quadratic['z_theta'])))
        zero=engine.increments(0,z)
        if any(zero):raise ArithmeticError('Nonzero heat increment at anchor')
        tail=engine.complete_swirl_heat_integral(z)
        if not(tail['reference']>0 and tail['heat_correction']<0 and tail['analytic_heat_truncation_absolute_bound']>0):
            raise ArithmeticError('Infinite swirl heat tail sign/retention failed')
        tail_neighbors={i:engine.complete_swirl_heat_integral(z+i*step)['heat_correction'] for i in (-2,-1,1,2)}
        tail_dZ=sum(w*tail_neighbors[i] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
        tail_Z_error=abs(tail_dZ/tail['heat_correction_Z']-1)
        if tail_Z_error>mp.mpf('1e-8'):raise ArithmeticError('Infinite swirl heat Z jet mismatch')
        report=dict(samples=reports,Rtail_increment_zero=True,heat_kernel_shared=True,
            actual_inner_anchor_once=True,finite_radius_heat_supported=True,
            quadrature_error_enclosed=False,heat_polynomial_truncation_enclosed=False,
            infinite_swirl_heat_tail_integrated=True,
            infinite_swirl_heat_Z_relative_error=mp.nstr(tail_Z_error,40),
            swirl_heat_tail={key:signed_log(tail[key],50) for key in ('reference','heat_correction','heat_correction_Z','analytic_heat_truncation_absolute_bound')},
            physical_Z_energy_weight_integrated=False,radial_energy_included=False,
            complete_five_moments=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    import sys
    run(tail_only='--tail-only' in sys.argv)
