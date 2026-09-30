"""Install a continuous MP pre-heat angular schedule in one shared object.

Heat point values use one continuous kernel and retain tiny deficits and
log-amplitude corrections separately. Heat-dependent integral targets still
use their declared Taylor approximations; full heat moments remain open.
"""
from decimal import Decimal, InvalidOperation
import types
import mpmath as mp
from lei_ren_part1_paper_continuous_incoming_angular import ContinuousIncomingAngular
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_continuous_heat_kernel import ContinuousHeatKernel


def _mp(value):
    return value if isinstance(value,mp.mpf) else mp.mpf(str(value))


class ContinuousAngularSchedule:
    def __init__(self,schedule,*,precision=443,primitive_precision=100):
        self.schedule=schedule
        self.precision=max(int(precision),schedule.decimal_precision)
        self.primitive_precision=int(primitive_precision)
        self.primitive=ContinuousIncomingAngular(str(schedule.logPstar),str(schedule.Md),
                                                 precision=self.primitive_precision)
        self.original_evaluate=schedule._evaluate_log_radius
        with mp.workdps(self.precision+20):
            self.heat_kernel=ContinuousHeatKernel(_mp(schedule.delta)/2,precision=self.precision+20)

    def heat_jet(self,log_radius,Z):
        with mp.workdps(self.precision+20):
            z=_mp(Z);invR=mp.exp(-_mp(log_radius))
            xi=2*(1-z*z)*invR
            row=self.heat_kernel.jet(xi)
            row.update(xi=xi,inv_R=invR,xi_Z=-4*z*invR)
            return row

    def J_mp(self,argument):
        with mp.workdps(self.precision+20):
            x=_mp(argument)
            if x<=0:return mp.mpf(0)
            if x>=1:return x-mp.mpf('.5')
            # This is numerical continuous quadrature, not certified accuracy.
            return self.primitive.J(x)

    def logA_mp(self,y):
        s=self.schedule
        with mp.workdps(self.precision+20):
            y=_mp(y);mu=_mp(s.mu);delta=_mp(s.delta)
            J=self.J_mp
            value=_mp(s.logPstar)+mp.mpf('.1')*y-mp.mpf('.6')*J(y)
            value-=mu*J(y-_mp(s.y_d))
            value-=(1-mu)*J(y-_mp(s.y_rel))
            value+=(1-delta/2)*J(y-_mp(s.y_rel)-1-_mp(s.Ts))
            return value

    def logA_decimal(self,y):
        # Combine MP primitives before Decimal conversion: an individual
        # flat primitive can have an exponent outside Decimal's range.
        with mp.workdps(self.precision+20):
            return Decimal(mp.nstr(self.logA_mp(y),self.precision))

    def slope_decimal(self,y):
        s=self.schedule
        with mp.workdps(self.precision+20):
            y=_mp(y);mu=_mp(s.mu);delta=_mp(s.delta)
            switch=lambda x:ContinuousAxialPulse.sigma_pair(x)[0]
            value=mp.mpf('.1')-mp.mpf('.6')*switch(y)-mu*switch(y-_mp(s.y_d))
            value-=(1-mu)*switch(y-_mp(s.y_rel))
            value+=(1-delta/2)*switch(y-_mp(s.y_rel)-1-_mp(s.Ts))
            return Decimal(mp.nstr(value,self.precision))

    def evaluate(self,log_radius,Z):
        s=self.schedule
        with mp.workdps(self.precision+20):
            log_decimal=s._coerce_log_radius(log_radius)
            # Keep exact/logarithmic stage offsets at schedule precision.
            from decimal import localcontext
            with localcontext() as context:
                context.prec=max(self.precision,s.decimal_precision)
                y_decimal=log_decimal-s.logRref
            y=_mp(y_decimal);z=_mp(Z)
            if abs(z)>=1:raise ValueError('Require |Z| < 1')
            # Legacy call supplies stage metadata and axial/heat source data.
            row=self.original_evaluate(log_decimal,float(z))
            heat_branch=y>=_mp(s.y_tail)
            heat_extra={}
            if heat_branch:
                heat=self.heat_jet(log_decimal,z)
                H=_mp(heat['H']);Hp=_mp(heat['H_prime']);xi=_mp(heat['xi'])
                t=y-_mp(s.y_tail);eps=_mp(s.epsilon)
                switch,derivative=ContinuousAxialPulse.sigma_pair(t)
                edge=(3-t)/2
                phi=mp.exp(-1/edge**2) if edge>0 else mp.mpf(0)
                phi_prime=2*phi/edge**3 if edge>0 else mp.mpf(0)
                factor=1-eps*phi
                factor_t=eps*phi_prime/2
                K0=(1-switch)*(1-eps)+switch*factor
                Kt0=-derivative*(1-eps)+derivative*factor+switch*factor_t
                dK= switch*factor*heat['H_minus_one']
                dKt=derivative*factor*heat['H_minus_one']
                dKt+=switch*(-xi*Hp*factor+heat['H_minus_one']*factor_t)
                relative=dK/K0
                log_correction=mp.log1p(relative)
                slope_correction=(dKt*K0-Kt0*dK)/(K0*K0*(1+relative))
                K=K0*(1+relative)
                a=(1+_mp(s.delta))/2
                logu=_mp(s._log_c_inf)-a*_mp(log_decimal)+mp.log(K0)+log_correction
                slope=-a+Kt0/K0+slope_correction
                LZ=switch*factor*Hp*(-4*z*mp.exp(-_mp(log_decimal)))/K
                heat_extra=dict(heat_deficit=heat['deficit'],heat_log_amplitude_correction=log_correction,
                    heat_logarithmic_slope_correction=slope_correction,
                    heat_kernel_method=heat['method'],heat_truncation_absolute_bounds=heat['truncation_absolute_bounds'],
                    heat_arithmetic_error_enclosed=False,heat_integral_targets_regenerated=False)
            else:
                logA=_mp(s._log_A(y_decimal))
                flat,flat_prime=ContinuousAxialPulse.sigma_pair((y-_mp(s.y_v))/_mp(s.Tf))
                zlog=mp.log(1+z*z)
                logu=logA-zlog+flat*(zlog-mp.log(2))
                LZ=-2*z/(1+z*z)*(1-flat)
                slope=_mp(self.slope_decimal(y_decimal))+flat_prime/_mp(s.Tf)*(zlog-mp.log(2))
            logF=logu-(mp.log(2)+_mp(log_decimal))/2
            u=mp.exp(logu);F=mp.exp(logF)
            def dec(value):
                try:return Decimal(mp.nstr(value,self.precision))
                except InvalidOperation:
                    # Decimal's exponent storage is bounded even when its
                    # context is wide. Preserve arbitrary-exponent MP jets.
                    return value
            row.update(log_angular_amplitude=dec(logu),logF=dec(logF),
                Utheta=u,F=F,Utheta_Z=u*LZ,F_Z=F*LZ,
                dlogU_dZ=dec(LZ),dlogF_dZ=dec(LZ),
                logarithmic_slope=dec(slope),logF_slope=dec(slope-mp.mpf('.5')),
                continuous_angular_preheat_installed=True,
                angular_heat_kernel_inherited=False,continuous_heat_kernel_installed=True,
                angular_primitive_working_precision=self.primitive_precision,
                angular_primitive_quadrature_enclosed=False)
            row.update(heat_extra)
            return row

    def install(self):
        s=self.schedule
        s._J=types.MethodType(lambda owner,x:self.J_mp(x),s)
        s._log_A=types.MethodType(lambda owner,y:self.logA_decimal(y),s)
        s._slope_s=types.MethodType(lambda owner,y:self.slope_decimal(y),s)
        s._sigma_J1=mp.mpf('.5');s._sigma_J1_decimal=Decimal('.5')
        # Recompute the same analytic heat prefactor with the new primitive.
        # Point heat jets are installed; complete heat integrals are open.
        with mp.workdps(self.precision+20):
            logc=_mp(s._log_A(s.y_tail))+(1+_mp(s.delta))*_mp(s.logR_tail)/2
            logc-=mp.log(2*(1-_mp(s.epsilon)))
            s._log_c_inf=Decimal(mp.nstr(logc,self.precision))
            if s.c_inf_decimal is not None:
                s.c_inf_decimal=Decimal(mp.nstr(mp.exp(logc),self.precision))
        s._evaluate_log_radius=types.MethodType(lambda owner,r,z:self.evaluate(r,z),s)
        # Bypass the legacy float-Z coercion at both public entry points.
        s.at_log_radius=types.MethodType(lambda owner,r,z:self.evaluate(r,z),s)
        s.at_radius=types.MethodType(lambda owner,r,z:self.evaluate(owner._coerce_radius_log(r),z),s)
        s._continuous_angular_provider=self
        return self


def install_continuous_angular_schedule(schedule,*,precision=443,primitive_precision=100):
    existing=getattr(schedule,'_continuous_angular_provider',None)
    if existing is not None:return existing
    return ContinuousAngularSchedule(schedule,precision=precision,
        primitive_precision=primitive_precision).install()


def run():
    import json
    from pathlib import Path
    from lei_ren_part1_paper_joined_outer import build_joined_field
    print('building angular schedule diagnostic',flush=True)
    field=build_joined_field();s=field.schedule
    provider=install_continuous_angular_schedule(s,precision=field.precision)
    with mp.workdps(provider.precision+20):
        z=mp.mpf('.3')
        points=[('reference_transition',mp.mpf('.3')),
            ('mu_transition',_mp(s.y_d)+mp.mpf('.3')),
            ('Z_flatten',_mp(s.y_v)+50),
            ('steep_transition',_mp(s.y_rel)+mp.mpf('.3')),
            ('delta_transition',_mp(s.y_rel)+1+_mp(s.Ts)+mp.mpf('.3'))]
        rows=[]
        for name,y in points:
            radius=Decimal(mp.nstr(_mp(s.logRref)+y,provider.precision))
            center=s.at_log_radius(radius,z)
            analytic=_mp(center['dlogU_dZ']);errors=[]
            for h in (mp.mpf('1e-4'),mp.mpf('1e-30')):
                values={i:_mp(s.at_log_radius(radius,z+i*h)['log_angular_amplitude']) for i in (-2,-1,1,2)}
                derivative=(values[-2]-8*values[-1]+8*values[1]-values[2])/(12*h)
                error=abs(derivative-analytic)
                if error>mp.mpf('1e-12'):raise AssertionError('MP angular Z jet failed independent difference')
                errors.append(dict(step=mp.nstr(h,20),absolute_error=mp.nstr(error,40)))
            rows.append(dict(stage=name,dlogU_dZ=mp.nstr(analytic,40),
                logarithmic_slope=str(center['logarithmic_slope']),Z_difference=errors))
        start=s.at_log_radius(s.logR_tail,z)
        expected=provider.logA_mp(_mp(s.y_tail))-mp.log(2)
        log_jump=abs(_mp(start['log_angular_amplitude'])-expected)
        if log_jump>mp.mpf('1e-240'):raise AssertionError('Heat normalization interface mismatch')
        result=dict(rows=rows,heat_interface_log_amplitude_jump=mp.nstr(log_jump,40),
            source_schedule_identity_preserved=provider.schedule is s,
            primitive_working_precision=provider.primitive_precision,
            arithmetic_working_precision=provider.precision,
            preheat_angular_Z_jets_installed=True,
            heat_kernel_inherited=False,heat_point_kernel_installed=True,
            heat_integral_targets_regenerated=False,angular_relative_correction_Z_complete=False,
            angular_pressure_and_moments_complete=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(preheat_Z_checks=True,heat_interface_log_jump=result['heat_interface_log_amplitude_jump'])),flush=True)
    return result


if __name__=='__main__':run()
