"""Second-Z propagation through the common R100--R110 switches."""
import mpmath as mp
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches


class SecondAxialPressureWidthExitSwitches(PressureWidthExitSwitches):
    def __init__(self,provider,**kwargs):
        if not getattr(provider.comparison,'automatic_Z_second',False):
            raise TypeError('Second switches require a second-Z source')
        super().__init__(provider,**kwargs)

    def _as_dual(self,value,tangent=0):
        if isinstance(value,AxialSecondJet):return value
        raise TypeError('Second switch source must retain all three axial slots')

    def _start(self,Z):
        with mp.workdps(self.precision):
            z=mp.mpf(str(Z));key=self._z_key(z)
            if key in self._start_cache:return self._start_cache[key]
            value=self.provider.evaluate_R(self.R0,z)
            fields=value['second_jet_fields'];moments=value['second_jet_moments']
            raw=value['second_jet_raw_quadratic_integrals']
            F0=fields['F'];U0=fields['Uz'];R0=self.R0_ring;F0sq=F0*F0
            P0=self.provider.bridge.initial(z)['P0']
            state=[self._dual(0),U0,moments['theta']/(F0*R0*R0),moments['z']/R0,
                   moments['theta_z']/(F0*R0*R0),raw['axial']/R0,
                   raw['swirl']/(F0sq*R0*R0),moments['p']/(F0sq*R0)]
            result=dict(z=z,provider=value,F0=F0,U0=U0,P0=P0,
                        P0_provider=fields['P'],moments=moments,
                        raw_axial=raw['axial'],raw_swirl=raw['swirl'],state=state,
                        Fbar=self.provider.auxiliary.endpoint(z)['F'],
                        raw_integrals_source='common second-Z continuation raw quadratic integrals')
            self._start_cache[key]=result
            return result

    def _output(self,radius,Z,state,region,x):
        result=super()._output(radius,Z,state,region,x)
        with mp.workdps(self.precision):
            start=self._start(Z);F0=start['F0'];R0=self.R0_ring
            g,u,theta,mz,mixed,axial,swirl,pressure=state
            F=F0*g.exp()
            moments=dict(theta=F0*R0*R0*theta,z=R0*mz,
                         theta_z=F0*R0*R0*mixed,
                         z_theta=R0*axial-F0*F0*R0*R0*swirl,
                         p=F0*F0*R0*pressure)
            fields=dict(F=F,Uz=u,P=start['P0']+moments['p'])
            for name,v in fields.items():result[name+'_ZZ']=v.second
            result['moments_ZZ']={k:v.second for k,v in moments.items()}
            raw=dict(axial=R0*axial,swirl=F0*F0*R0*R0*swirl)
            result['raw_quadratic_integrals_ZZ']={k:v.second for k,v in raw.items()}
            result.update(second_jet_fields=fields,second_jet_moments=moments,
                          second_jet_raw_quadratic_integrals=raw,P0_ZZ=start['P0'].second)
            z=mp.mpf(Z);dt=self.comparison.delta;L=1-dt*z*z;root=mp.sqrt(2*radius)
            m=result['moments']['z'];mz1=result['moments_Z']['z'];mz2=result['moments_ZZ']['z']
            N=2*z*radius*result['Uz']-(1-dt)*z*m-(1-z*z)*mz1
            Nz=2*radius*(result['Uz']+z*result['Uz_Z'])-(1-dt)*(m+z*mz1)+2*z*mz1-(1-z*z)*mz2
            result['Ur_Z']=Nz/(L*root)+2*dt*z*N/(L*L*root)
            result['Ur_ZZ']=None
            result.update(second_Z_moments_available=True,radial_first_Z_available=True,
                          second_Z_remainder_enclosed=False)
            return result

    def metadata(self):
        result=super().metadata()
        result.update(second_Z_installed=True,second_Z_remainder_enclosed=False)
        return result
