"""One source candidate for corrected swirl, axial and radial components.

Uses actual coefficient inputs; radial velocity comes from the SAME axial
primitive. The temporary reference is not a regular axis core. Arbitrary
exponents and Decimal stage coordinates are mandatory on the source route.
"""
from decimal import Decimal, localcontext
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from lei_ren_part1_paper_axial_energy_tail import build_default_tail
from lei_ren_part1_paper_axial_incoming import incoming_axial_moments
from lei_ren_part1_paper_axial_pulse_moments import normalized_pulse_integral
from lei_ren_part1_paper_axial_correction import (
    solve_actual_axial, from_signed_log, signed_log, axial_factor)
from lei_ren_part1_paper_outer_closure import _paper_raw_bump
from lei_ren_part1_paper_corrected_pressure import CorrectedPressureAdapter


class CorrectedSourceProfile:
    def __init__(self,*,precision=160,order=128,schedule=None,angular=None,tail=None,
                 match_waiting=False):
        self.precision=precision; self.order=order
        if tail is not None:
            if schedule is not None and tail.schedule is not schedule:
                raise ValueError('explicit schedule and tail must be the same object')
            if angular is not None and tail.correction is not angular:
                raise ValueError('explicit angular correction and tail must be shared')
            if match_waiting:
                raise ValueError('match_waiting cannot be used with an explicit tail')
            self.tail=tail
            construction_mode='explicit_tail'
        else:
            self.tail=build_default_tail(precision=precision,quadrature_order=order,
                    schedule=schedule,correction=angular,match_waiting=match_waiting)
            construction_mode='explicit_schedule' if schedule is not None or angular is not None else 'default'
        self.schedule=self.tail.schedule; self.angular=self.tail.correction
        if schedule is not None and self.schedule is not schedule and tail is None and not match_waiting:
            raise ValueError('tail factory returned a different schedule without matching requested')
        if angular is not None and self.angular is not angular:
            raise ValueError('tail factory did not preserve the supplied angular correction')
        self.pressure=CorrectedPressureAdapter(self.schedule,self.angular,
                    precision=precision,quadrature_order=order)
        self.pulse=[normalized_pulse_integral(self.schedule.mu,i,
                    precision=precision,order=order) for i in (1,2)]
        waiting_was_matched = (
            tail is None
            and schedule is None
            and angular is None
        ) or bool(match_waiting)
        self._construction_metadata={
            'mode':construction_mode,
            'explicit_schedule':schedule is not None,
            'explicit_angular_correction':angular is not None,
            'explicit_tail':tail is not None,
            'match_waiting':waiting_was_matched,
            'schedule_identity_preserved':tail is not None or not match_waiting,
            'tail_schedule_identity':self.schedule is self.tail.schedule,
            'tail_angular_identity':self.angular is self.tail.correction,
            'source_scaled_factory_supported':True,
        }

    def offset(self,logR,checkpoint):
        with localcontext() as ctx:
            ctx.prec=max(self.schedule.decimal_precision,self.precision)
            return Decimal(str(logR))-checkpoint

    def log_at(self,checkpoint,offset):
        with localcontext() as ctx:
            ctx.prec=max(self.schedule.decimal_precision,self.precision)
            return checkpoint+Decimal(str(offset))

    @lru_cache(maxsize=64)
    def coefficients(self,Z):
        incoming=incoming_axial_moments(self.schedule,Z,order=64,precision=self.precision)
        future=self.tail.evaluate(Z,quadrature_order=self.order)
        base=[incoming['row_normalization'][key] for key in ('scaled_base_m1','scaled_base_m2')]
        with mp.workdps(self.precision):
            prior=from_signed_log(incoming['E_prior_mu_Mztheta_over_RpEp2'])
            target=(1-mp.exp(-26))/4-prior+mp.mpf(future['energy_target_contribution_nominal'])
            axial=solve_actual_axial(self.schedule.mu,base,self.pulse,
                    mp.nstr(target,self.precision),precision=self.precision,order=self.order)
        return {'Z':float(Z),'incoming':incoming,'axial':axial,
                'angular':self.angular.coefficients(Z)}

    def values(self,logR,Z):
        """Return mp Utheta, Uz and a separate angular correction term."""
        base=self.schedule.at_log_radius(logR,Z)
        with mp.workdps(self.precision):
            amplitude=mp.exp(mp.mpf(str(base['log_angular_amplitude'])))
            angular_offset=self.offset(logR,self.schedule.logR_rel)
            correction=mp.mpf(0)
            if Decimal('-3.15')<=angular_offset<=Decimal('-.85'):
                row=self.angular.relative_bump(float(angular_offset),Z,
                            coefficients=self.coefficients(float(Z))['angular'])
                correction=from_signed_log(row)
            Uz=mp.mpf(str(base['Uz']))
            start=self.offset(logR,self.schedule.logR_p)
            end=self.offset(logR,self.schedule.logR_v)
            if start>=0 and end<=0:
                with localcontext() as ctx:
                    ctx.prec=max(self.schedule.decimal_precision,self.precision)
                    xi=self.schedule.mu*start
                receipt=self.coefficients(float(Z))['axial']
                if Decimal('-3.15')<=end<=Decimal('-.85'):
                    factor=axial_factor(receipt,end_offset=mp.mpf(str(end)),precision=self.precision)
                    Uz=amplitude*from_signed_log(factor)
                else:
                    # Use the complementary flat step directly. Subtracting
                    # sigma from one would erase its nonzero endpoint tail.
                    from lei_ren_part1_paper_axial_primitive import source_pulse_value_mp
                    x=mp.mpf(str(xi))
                    gp=source_pulse_value_mp(str(xi),precision=self.precision)
                    if 0<x<11 and gp==0:
                        raise ArithmeticError('Nonzero startup pulse is below the primitive tolerance; a bounded log evaluation is required')
                    Uz=amplitude*mp.mpf(receipt['a_p'])*gp
            return {'Utheta':amplitude*(1+correction),'Uz':Uz,
                    'angular_relative_correction':correction,
                    'base_log_Utheta':str(base['log_angular_amplitude']),
                    'angular_log1p_correction':mp.nstr(mp.log1p(correction),self.precision)}

    def _end_bump_integrals(self,end,*,order=192):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        raw=[mp.mpf(str(float(_paper_raw_bump(float(n))))) for n in nodes]
        norm=sum(mp.mpf(str(float(w)))*r for w,r in zip(weights,raw))
        lam=mp.mpf('.5')-mp.mpf(str(self.schedule.mu)); ell=mp.mpf('.15')
        integrals=[]
        for center in (-3,-1):
            left=mp.mpf(center)-ell; right=min(mp.mpf(str(end)),mp.mpf(center)+ell)
            if right<=left:
                integrals.append(mp.mpf(0)); continue
            midpoint=(left+right)/2; half=(right-left)/2
            total=mp.mpf(0)
            for index,(n,w) in enumerate(zip(nodes,weights)):
                t=midpoint+half*mp.mpf(str(float(n)))
                # Full supports use exactly the solve's canonical nodes.
                # Reconstructing a node through (t-center)/ell introduces
                # another floating-point perturbation of the same bump.
                bump=raw[index] if right==mp.mpf(center)+ell else mp.mpf(
                    str(float(_paper_raw_bump(float((t-center)/ell)))))
                beta=bump/(ell*norm)
                total+=half*mp.mpf(str(float(w)))*mp.exp(lam*t)*beta
            integrals.append(total)
        return integrals

    def axial_average(self,logR,Z):
        """A=M^z/R from actual reference, pulse and end bumps.

        Pulse primitive approximation/bounds are supplied by its dedicated
        module. Exterior mass is reintegrated; it is never forced to zero.
        """
        from lei_ren_part1_paper_axial_primitive import normalized_pulse_primitive
        with mp.workdps(self.precision):
            y=self.offset(logR,self.schedule.logRref)
            if y<=0:return mp.mpf(4)*Z
            start=self.offset(logR,self.schedule.logR_p)
            data=self.coefficients(float(Z)); incoming=data['incoming']; axial=data['axial']
            if start<0:
                # Beyond axial turnoff the incoming primitive is constant.
                cutoff_end=mp.exp(mp.mpf(str(self.schedule.Md)))
                if mp.mpf(str(y))>=cutoff_end:
                    Iz=mp.mpf(incoming['dimensionless_integrals']['I_z'])
                    return Iz*mp.exp(-mp.mpf(str(y)))
                nodes,weights=np.polynomial.legendre.leggauss(self.order)
                total=mp.mpf(4)*Z
                # Split at y=1 where the slow axial turnoff starts.
                intervals=[(mp.mpf(0),min(mp.mpf(str(y)),mp.mpf(1)))]
                if y>1:intervals.append((mp.mpf(1),mp.mpf(str(y))))
                for left,right in intervals:
                    half=(right-left)/2; midpoint=(right+left)/2
                    for node,weight in zip(nodes,weights):
                        v=midpoint+half*mp.mpf(str(float(node)))
                        row=self.schedule.at_log_radius(self.log_at(self.schedule.logRref,mp.nstr(v,self.precision)),Z)
                        total+=half*mp.mpf(str(float(weight)))*mp.exp(v)*mp.mpf(str(row['Uz']))
                return total*mp.exp(-mp.mpf(str(y)))
            end=self.offset(logR,self.schedule.logR_v)
            lam=mp.mpf('.5')-mp.mpf(str(self.schedule.mu))
            ap=mp.mpf(axial['a_p'])
            if end<=0:
                xi=mp.mpf(str(self.schedule.mu))*mp.mpf(str(start))
                if xi>=11:
                    # After pulse support, use exactly the full-row input
                    # used by the coefficient solve. Recomputing it at a
                    # different precision can create a jump at Rv even if
                    # both evaluations are individually very accurate.
                    J=mp.exp(mp.mpf(self.pulse[0]['log_normalized_pulse_integral'])
                        -lam*mp.mpf(str(end)))
                else:
                    primitive=normalized_pulse_primitive(self.schedule.mu,1,mp.nstr(xi,self.precision),
                                precision=self.precision,order=self.order)
                    J=mp.exp(mp.mpf(primitive['log_J'])) if primitive['log_J'] is not None else mp.mpf(0)
                m1=from_signed_log(incoming['m1_Mz_over_RpEp'])
                Y=m1*mp.exp(-lam*mp.mpf(str(start)))+ap*J
                if end>=Decimal('-3.15'):
                    c=[from_signed_log(row) for row in axial['c']]
                    Y+=mp.exp(-lam*mp.mpf(str(end)))*sum(
                        coefficient*integral for coefficient,integral in zip(c,
                            self._end_bump_integrals(end,order=axial['quadrature_order'])))
                E=mp.exp(mp.mpf(str(self.schedule.at_log_radius(logR,Z)['log_angular_amplitude'])))
                return E*Y
            c=[from_signed_log(row) for row in axial['c']]
            norm=from_signed_log(incoming['row_normalization']['scaled_base_m1'])\
                 +ap*mp.exp(mp.mpf(self.pulse[0]['log_normalized_pulse_integral']))\
                 +sum(coefficient*integral for coefficient,integral in zip(c,
                    self._end_bump_integrals(0,order=axial['quadrature_order'])))
            Ev=mp.exp(mp.mpf(str(self.schedule.at_log_radius(self.schedule.logR_v,Z)['log_angular_amplitude'])))
            return norm*Ev*mp.exp(-mp.mpf(str(end)))

    def radial_flux(self,logR,Z,*,step=1e-4):
        """V/R from source (3.9), with Z derivative of the same primitive."""
        if abs(float(Z))+2*step>=1:
            raise ValueError('Derivative stencil must remain inside |Z|<1')
        with mp.workdps(self.precision):
            z=mp.mpf(str(Z)); h=mp.mpf(str(step))
            average=self.axial_average(logR,float(Z))
            derivative=(self.axial_average(logR,float(z-2*h))-8*self.axial_average(logR,float(z-h))
                        +8*self.axial_average(logR,float(z+h))-self.axial_average(logR,float(z+2*h)))/(12*h)
            Uz=self.values(logR,float(Z))['Uz']; delta=mp.mpf(str(self.schedule.delta))
            return (2*z*Uz-(1-delta)*z*average-(1-z*z)*derivative)/(1-delta*z*z)

    def cylindrical_chart(self,logR,Z,logq,*,nu='.01',step=1e-4):
        with mp.workdps(self.precision):
            root_nu=mp.sqrt(mp.mpf(str(nu))); qlog=mp.mpf(str(logq)); rlog=mp.mpf(str(logR))
            if root_nu<=0:raise ValueError('nu must be positive')
            h=mp.mpf(str(self.schedule.delta))/2
            values=self.values(logR,Z)
            radial=root_nu*mp.exp((rlog-qlog)/2)/mp.sqrt(2)*self.radial_flux(logR,Z,step=step)
            scale=root_nu*mp.exp((-mp.mpf('.5')-h)*qlog)
            return {'u_r':signed_log(radial,self.precision),
                    'u_theta':signed_log(scale*values['Utheta'],self.precision),
                    'u_z':signed_log(scale*values['Uz'],self.precision),
                    'scope':'Unlocalized source chart candidate; temporary axis reference and unresolved closure uncertainties remain.'}

    def cartesian_chart(self,logR,Z,logq,phi=0,*,nu='.01',step=1e-4):
        """Cartesian velocity and physical coordinates, as signed logs."""
        with mp.workdps(self.precision):
            cylindrical=self.cylindrical_chart(logR,Z,logq,nu=nu,step=step)
            radial=from_signed_log(cylindrical['u_r']); swirl=from_signed_log(cylindrical['u_theta'])
            axial=from_signed_log(cylindrical['u_z']); angle=mp.mpf(str(phi))
            root_nu=mp.sqrt(mp.mpf(str(nu))); qlog=mp.mpf(str(logq))
            h=mp.mpf(str(self.schedule.delta))/2; z=mp.mpf(str(Z))
            radius=root_nu*mp.sqrt(2)*mp.exp((qlog+mp.mpf(str(logR)))/2)
            coordinates=[radius*mp.cos(angle),radius*mp.sin(angle),
                         root_nu*mp.exp((mp.mpf('.5')-h)*qlog)*z]
            velocity=[radial*mp.cos(angle)-swirl*mp.sin(angle),
                      radial*mp.sin(angle)+swirl*mp.cos(angle),axial]
            return {'xyz':[signed_log(v,self.precision) for v in coordinates],
                    'tau':signed_log(mp.exp(qlog)*(1-z*z),self.precision),
                    'uvw':[signed_log(v,self.precision) for v in velocity],
                    'scope':cylindrical['scope']}

    def _inverse_axial_chart(self,z,tau,viscosity):
        delta=mp.mpf(str(self.schedule.delta)); a=(1-delta)/2
        if z==0:return mp.mpf(0),mp.log(tau)
        # Solve in logit(Z^2), keeping the thin d=1-Z^2 separately.
        log_s=mp.log(abs(z)/mp.sqrt(viscosity))-a*mp.log(tau)
        target=2*log_s
        w=target/(1-delta) if target>0 else target
        for _ in range(8):
            softplus=w+mp.log1p(mp.exp(-w)) if w>0 else mp.log1p(mp.exp(w))
            logistic=1/(1+mp.exp(-w))
            w-=(w-delta*softplus-target)/(1-delta*logistic)
        softplus=w+mp.log1p(mp.exp(-w)) if w>0 else mp.log1p(mp.exp(w))
        d=mp.exp(-softplus)
        return mp.sign(z)*mp.sqrt(1-d),mp.log(tau)+softplus

    def velocity_from_tau(self,x,y,z,tau,*,nu='.01',z_inner=.25,z_outer=.5,step=1e-4):
        """Callable physical (u,v,w), with streamfunction axial localization.

        Inputs may be decimal strings with arbitrary exponents. The source
        temporary axis reference remains unaccepted until core replacement.
        """
        from openai_ns_reconstruction.paper_compact_field import cutoff
        with mp.workdps(self.precision):
            x,y,z,tau=(mp.mpf(str(v)) for v in (x,y,z,tau))
            viscosity=mp.mpf(str(nu))
            if tau<=0 or viscosity<=0 or not 0<z_inner<z_outer:
                raise ValueError('Require tau,nu>0 and 0<z_inner<z_outer')
            if abs(z)>=z_outer:
                return [signed_log(mp.mpf(0),self.precision) for _ in range(3)]
            B,derivative=cutoff(float(z*z),z_inner*z_inner,z_outer*z_outer)
            B=mp.mpf(str(B)); Bz=2*z*mp.mpf(str(derivative))
            delta=mp.mpf(str(self.schedule.delta))
            Z,logq=self._inverse_axial_chart(z,tau,viscosity)
            radius=mp.sqrt(x*x+y*y)
            if radius==0:
                axial=mp.sqrt(viscosity)*mp.exp((-mp.mpf('.5')-delta/2)*logq)*4*Z*B
                return [signed_log(mp.mpf(0),self.precision),signed_log(mp.mpf(0),self.precision),
                        signed_log(axial,self.precision)]
            if abs(float(Z))+2*step>=1:
                raise ValueError('Physical chart is too near |Z|=1 for the current coefficient derivative stencil')
            logR=2*mp.log(radius)-mp.log(2*viscosity)-logq
            logR_string=mp.nstr(logR,self.precision); logq_string=mp.nstr(logq,self.precision)
            cylindrical=self.cylindrical_chart(logR_string,float(Z),logq_string,nu=nu,step=step)
            radial=B*from_signed_log(cylindrical['u_r'])
            average=self.axial_average(logR_string,float(Z))
            psi_over_r=viscosity*mp.exp(-delta*logq/2+logR/2)/mp.sqrt(2)*average
            radial-=Bz*psi_over_r
            swirl=B*from_signed_log(cylindrical['u_theta']); axial=B*from_signed_log(cylindrical['u_z'])
            return [signed_log((x*radial-y*swirl)/radius,self.precision),
                    signed_log((y*radial+x*swirl)/radius,self.precision),
                    signed_log(axial,self.precision)]

    def pressure_at_reference(self,Z):
        return self.pressure.corrected_pressure_at_Rref(Z)

    def pressure_at_log_radius(self,logR,Z):
        return self.pressure.corrected_pressure_at_log_radius(logR,Z)

    def pressure_from_tau(self,x,y,z,tau,*,nu='.01',z_inner=.25,z_outer=.5):
        """Physical nominal pressure of the SAME localized source candidate."""
        from openai_ns_reconstruction.paper_compact_field import cutoff
        with mp.workdps(self.precision):
            x,y,z,tau=(mp.mpf(str(v)) for v in (x,y,z,tau)); viscosity=mp.mpf(str(nu))
            if tau<=0 or viscosity<=0 or not 0<z_inner<z_outer:
                raise ValueError('Require tau,nu>0 and 0<z_inner<z_outer')
            if abs(z)>=z_outer:
                return {'pressure':signed_log(mp.mpf(0),self.precision),'outside_axial_support':True}
            B,_=cutoff(float(z*z),z_inner*z_inner,z_outer*z_outer)
            Z,logq=self._inverse_axial_chart(z,tau,viscosity)
            radius=mp.sqrt(x*x+y*y)
            if radius==0:
                data=self.pressure_at_reference(float(Z))
                nominal=mp.mpf(data['corrected_axis_P0_over_Pstar_squared'])\
                    *mp.exp(2*mp.mpf(str(self.schedule.logPstar)))
            else:
                logR=2*mp.log(radius)-mp.log(2*viscosity)-logq
                data=self.pressure_at_log_radius(mp.nstr(logR,self.precision),float(Z))
                nominal=from_signed_log(data['corrected_pressure_nominal'])
            scale=viscosity*mp.mpf(str(B))**2*mp.exp((-1-mp.mpf(str(self.schedule.delta)))*logq)
            return {'pressure':signed_log(scale*nominal,self.precision),
                    'profile_pressure':data,'physical_scale':signed_log(scale,self.precision),
                    'scope':'Nominal pressure with inherited profile quadrature/heat uncertainty; not a full momentum certificate.'}

    def velocity(self,x,y,z,t,*,T=1,**options):
        with mp.workdps(self.precision):
            time=mp.mpf(str(t)); terminal=mp.mpf(str(T))
            if not 0<=time<terminal:
                raise ValueError('Require 0<=t<T; use the tau interface near T')
            return self.velocity_from_tau(x,y,z,mp.nstr(terminal-time,self.precision),**options)

    def velocity_values_from_tau(self,x,y,z,tau,**options):
        """Numeric (u,v,w) as mpmath scalars instead of receipt dictionaries."""
        with mp.workdps(self.precision):
            return tuple(from_signed_log(row) for row in self.velocity_from_tau(x,y,z,tau,**options))

    def pressure_correction(self,logR,Z):
        offset=self.offset(logR,self.schedule.logR_rel)
        return self.pressure.pressure_correction_at_stage(Z,float(offset))

    def metadata(self):
        return {'source':'https://arxiv.org/html/2609.35406v1',
                'schedule':self.schedule.metadata(),
                'construction':self._construction_metadata,
                'components':'Shared candidate angular/axial coefficients; radial recovered from same axial primitive',
                'pressure':'Same angular correction object, separate baseline and signed-log correction',
                'axis_regular':False,'finite_global_energy_certified':False,
                'source_core_reference_scale_compatible':False if self._construction_metadata['mode']=='default' else None,
                'core_scale_limitation':'The default logPstar=14/logRref=10 fail the core-reference relation. Explicit schedules require separate Cstar/Lambda/domain and contraction checks; compatibility is not inferred from injection.',
                'stress_cone_certified':False,'scale_recursion_established':False,
                'physical_localization_applied':'velocity_from_tau only; chart methods remain unlocalized',
                'full_outer_closed':False}


def run():
    profile=CorrectedSourceProfile()
    schedule=profile.schedule; z=.3
    with localcontext() as ctx:
        ctx.prec=profile.precision
        bulk=schedule.logR_p+Decimal(5)/schedule.mu
    samples=[]
    points=[('reference',schedule.logRref),('bulk_pulse',bulk),
            ('end_bump',profile.log_at(schedule.logR_v,-3)),
            ('angular_bump',profile.log_at(schedule.logR_rel,-3)),
            ('heat_exterior',profile.log_at(schedule.logR_b,2))]
    with mp.workdps(profile.precision):
        for name,logR in points:
            values=profile.values(logR,z)
            sample={'stage':name,'logR':str(logR),
                    'Utheta':signed_log(values['Utheta'],profile.precision),
                    'Uz':signed_log(values['Uz'],profile.precision),
                    'angular_log1p_correction':values['angular_log1p_correction']}
            if name in ('reference','bulk_pulse','angular_bump','heat_exterior'):
                sample['same_candidate_pressure']=profile.pressure_at_log_radius(logR,z)
            if name in ('reference','bulk_pulse','end_bump','heat_exterior'):
                sample['axial_average']=signed_log(profile.axial_average(logR,z),profile.precision)
            samples.append(sample)
        # Independent radial derivative of the SAME primitive on the bulk
        # pulse. M_R=Uz is necessary for the streamfunction construction.
        h=mp.mpf('.001')
        averages=[profile.axial_average(profile.log_at(bulk,mp.nstr(i*h,20)),z)
                  for i in (-2,-1,1,2)]
        derivative=(averages[0]-8*averages[1]+8*averages[2]-averages[3])/(12*h)
        center=profile.axial_average(bulk,z); Uz=profile.values(bulk,z)['Uz']
        primitive_relative=float(abs((derivative+center-Uz)/Uz))
        if primitive_relative>1e-8:
            raise ArithmeticError('Same axial primitive failed independent radial identity')
        radial=profile.radial_flux(bulk,z)
        radial_fine=profile.radial_flux(bulk,z,step=5e-5)
        radial_difference=float(abs((radial_fine-radial)/radial_fine)) if radial_fine else 0.
        physical=profile.cylindrical_chart(bulk,z,'-4')
        cartesian=profile.cartesian_chart(bulk,z,'-4',phi='.4')
        xyz=[from_signed_log(row) for row in cartesian['xyz']]
        tau=from_signed_log(cartesian['tau'])
        recovered=profile.velocity_from_tau(*(mp.nstr(v,profile.precision) for v in xyz),
                                            mp.nstr(tau,profile.precision))
        expected=[from_signed_log(row) for row in cartesian['uvw']]
        actual=[from_signed_log(row) for row in recovered]
        roundtrip=[float(abs((v-u)/u)) if u else float(abs(v)) for u,v in zip(expected,actual)]
        if max(roundtrip)>1e-8:
            raise ArithmeticError('Physical-coordinate Cartesian velocity roundtrip failed')
        physical_pressure=profile.pressure_from_tau(*(mp.nstr(v,profile.precision) for v in xyz),
                                                    mp.nstr(tau,profile.precision))
        nominal=from_signed_log(samples[1]['same_candidate_pressure']['corrected_pressure_nominal'])
        pressure_expected=mp.mpf('.01')*mp.exp(4*(1+mp.mpf(str(schedule.delta))))*nominal
        pressure_actual=from_signed_log(physical_pressure['pressure'])
        pressure_roundtrip=float(abs((pressure_actual-pressure_expected)/pressure_expected))
        if pressure_roundtrip>1e-8:
            raise ArithmeticError('Same-candidate physical pressure roundtrip failed')
    report={'candidate':profile.metadata(),'samples':samples,
            'bulk_primitive_radial_identity_relative_error':primitive_relative,
            'bulk_radial_Z_stencil_relative_difference':radial_difference,
            'bulk_physical_chart_velocity':physical,
            'bulk_cartesian_chart':cartesian,
            'physical_velocity_roundtrip_relative_errors':roundtrip,
            'bulk_physical_pressure':physical_pressure,
            'physical_pressure_roundtrip_relative_error':pressure_roundtrip,
            'same_candidate_reference_pressure':profile.pressure_at_reference(z),
            'scope':'Unified candidate with independent bulk primitive check, not a global divergence/energy/core/PDE certificate.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({'primitive_identity_relative':primitive_relative,
                      'radial_Z_refinement_relative':radial_difference}))
    return report


if __name__=='__main__':run()
