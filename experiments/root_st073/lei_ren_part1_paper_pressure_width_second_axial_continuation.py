"""Second-Z slots of the analytic post-collar exponential continuation."""
import mpmath as mp
from lei_ren_part1_paper_pressure_width_continuation import PressureWidthExitContinuation


class SecondAxialPressureWidthExitContinuation(PressureWidthExitContinuation):
    def __init__(self,bridge,**kwargs):
        if not getattr(bridge.comparison,'automatic_Z_second',False):
            raise TypeError('Second continuation requires second-Z comparison')
        super().__init__(bridge,**kwargs)

    def evaluate_R(self,R,Z):
        with mp.workdps(self.precision):
            result=super().evaluate_R(R,Z)
            radius=mp.mpf(str(R));z=mp.mpf(str(Z));fun=self.functions(z)
            y=(self.jet(radius)/fun['Rb']).log()
            raw=dict(F=fun['F'].evaluate(y),Uz=fun['Uz'].evaluate(y),R=self.jet(radius))
            moments={k:v.evaluate(y) for k,v in fun['moments'].items()}
            raw['P']=fun['start']['P']+moments['p']-fun['start']['moments']['p']
            for key,value in raw.items():result[key+'_ZZ']=value.second
            result['moments_ZZ']={k:v.second for k,v in moments.items()}
            quadratics={k:v.evaluate(y) for k,v in fun['raw_quadratic'].items()}
            result['raw_quadratic_integrals_ZZ']={k:v.second for k,v in quadratics.items()}
            datum=self.bridge.initial(z)['P0']
            result.update(P0=datum.value,P0_Z=datum.tangent,P0_ZZ=datum.second,
                          second_jet_fields=raw,second_jet_moments=moments,
                          second_jet_raw_quadratic_integrals=quadratics)
            # Differentiate the radial moment recovery at fixed R.
            dt=self.comparison.delta;L=1-dt*z*z;root=mp.sqrt(2*radius)
            m=result['moments']['z'];mz=result['moments_Z']['z'];mzz=result['moments_ZZ']['z']
            N=2*z*radius*result['Uz']-(1-dt)*z*m-(1-z*z)*mz
            Nz=2*radius*(result['Uz']+z*result['Uz_Z'])-(1-dt)*(m+z*mz)+2*z*mz-(1-z*z)*mzz
            result['Ur_Z']=Nz/(L*root)+2*dt*z*N/(L*L*root)
            result['Ur_ZZ']=None
            result.update(second_Z_moments_available=True,radial_first_Z_available=True,
                          second_Z_remainder_enclosed=False)
            return result

    def metadata(self):
        result=super().metadata()
        result.update(second_Z_installed=True,radial_first_Z_available=True,
                      second_Z_remainder_enclosed=False)
        return result
