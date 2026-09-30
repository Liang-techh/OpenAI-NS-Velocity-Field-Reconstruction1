"""Second-Z moments and first-Z radial velocity through finite exit RK."""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_axial_bridge import AxialPressureWidthExitBridge


class SecondAxialPressureWidthExitBridge(AxialPressureWidthExitBridge):
    @lru_cache(maxsize=64)
    def evaluate(self,s,Z):
        with mp.workdps(self.precision):
            result=dict(super().evaluate(s,Z))
            raw=PressureWidthExitBridge.evaluate(self,s,Z)
            for name in ('y','R','F','Uz','P','log_F_over_Fa','F_R_over_F','Uz_R','chi'):
                if not hasattr(raw[name],'second'):
                    raise TypeError('Second-Z bridge requires second-jet comparison source')
                result[name+'_ZZ']=raw[name].second
            result['moments_ZZ']={k:v.second for k,v in raw['moments'].items()}
            result['second_jet_fields']={k:raw[k] for k in ('R','F','Uz','P')}
            result['second_jet_moments']=raw['moments']
            z=mp.mpf(str(Z));dt=self.comparison.delta;R=result['R'];root=(2*R).sqrt()
            if raw['R'].tangent.atoms or raw['R'].second.atoms:
                raise ValueError('Radial first-Z recovery here requires fixed R')
            m=result['moments']['z'];mz=result['moments_Z']['z'];mzz=result['moments_ZZ']['z']
            N=2*z*R*result['Uz']-(1-dt)*z*m-(1-z*z)*mz
            Nz=2*R*(result['Uz']+z*result['Uz_Z'])-(1-dt)*(m+z*mz)+2*z*mz-(1-z*z)*mzz
            L=1-dt*z*z
            result['Ur_Z']=Nz/(L*root)+2*dt*z*N/(L*L*root)
            result['Ur_ZZ']=None
            result.update(second_Z_moments_available=True,radial_first_Z_available=True,
                          radial_second_Z_available=False,second_Z_remainder_enclosed=False)
            return result

    def metadata(self):
        result=super().metadata()
        result.update(second_Z_installed=True,radial_first_Z_available=True,
                      derivative_of_finite_RK_construction=True,
                      second_Z_remainder_enclosed=False)
        return result
