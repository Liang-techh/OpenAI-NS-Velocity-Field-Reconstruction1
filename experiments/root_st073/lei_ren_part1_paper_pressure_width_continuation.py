"""Analytic finite-ring prescribed continuation after the tiny exit collar.

Frozen auxiliary drivers have exact exponential powers. Their prescribed
positive-epsilon changes are integrated, not discarded below MP precision.
The exponential series terminates in the retained width ring; it is not an
enclosure of the untruncated field, core or collar RK construction.
"""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_exponential_polynomial import ExponentialPolynomial
from lei_ren_part1_paper_pressure_width_frozen_comparison import PressureWidthFrozenComparison
from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


class PressureWidthExitContinuation:
    def __init__(self, bridge, *, max_R=100):
        if not getattr(bridge.comparison,'automatic_Z_tangent',False):
            raise TypeError('Derivative-aware continuation requires an axial comparison')
        if not bridge.epsilon_width_tied:
            raise ValueError('Analytic finite-width continuation requires epsilon=W')
        self.bridge=bridge; self.comparison=bridge.comparison
        self.auxiliary=PressureWidthFrozenComparison(self.comparison)
        self.precision=bridge.precision; self.width_order=bridge.width_order
        self.pressure_order=bridge.pressure_order
        with mp.workdps(self.precision):
            self.max_R=mp.mpf(str(max_R))
            if self.max_R<=self.comparison.Ra:
                raise ValueError('max_R must exceed the core exit radius')
        self.jet=self.comparison.jet

    @lru_cache(maxsize=16)
    def functions(self,Z):
        with mp.workdps(self.precision):
            # Raw bridge values retain their automatically differentiated slots.
            start=PressureWidthExitBridge.evaluate(self.bridge,2,Z)
            rb=start['R']; epsilon=self.comparison.W
            samples=[self.auxiliary.evaluate(rb*k,Z) for k in (1,2,4)]
            def recover(key,powers):
                matrix=mp.matrix([[mp.mpf(k)**p for p in powers] for k in (1,2,4)])
                inverse=matrix**-1
                values=[v[key] for v in samples]
                if key=='I_z':
                    values=[(v['R']/2).sqrt()*v['I_z'] for v in samples]
                return [sum((inverse[i,j]*values[j] for j in range(3)),self.jet(0)) for i in range(3)]
            dc=recover('D',(1,0,-1)); jc=recover('I_z',(2,1,0))
            D=ExponentialPolynomial({(p,0):v for p,v in zip((1,0,-1),dc)})
            J=ExponentialPolynomial({(p,0):v for p,v in zip((2,1,0),jc)})
            ginc=D.integral_polynomial()*(-epsilon/2)
            exponential=ExponentialPolynomial({(0,0):self.jet(1)})
            term=exponential
            for n in range(1,self.width_order+1):
                term=term*ginc/n; exponential=exponential+term
            F=exponential*start['F']
            U=ExponentialPolynomial({(0,0):start['Uz']}) + (
                exponential*J*(-epsilon*start['F']/samples[0]['F'])).integral_polynomial()
            R=ExponentialPolynomial({(1,0):rb})
            integrands=dict(theta=2*R*R*F,z=R*U,theta_z=2*R*R*F*U,
                z_theta=R*U*U-R*R*F*F,p=R*F*F)
            moments={k:ExponentialPolynomial({(0,0):start['moments'][k]})+
                     v.integral_polynomial() for k,v in integrands.items()}
            collar_state=start['normalized_state']; Ra=self.comparison.Ra
            Fa=start['Fa']
            raw_quadratic=dict(
                axial=ExponentialPolynomial({(0,0):Ra*collar_state[5]})+(R*U*U).integral_polynomial(),
                swirl=ExponentialPolynomial({(0,0):Fa*Fa*Ra*Ra*collar_state[6]})+(R*R*F*F).integral_polynomial())
            return dict(start=start,Rb=rb,D=D,J=J,ginc=ginc,F=F,Uz=U,moments=moments,
                raw_quadratic=raw_quadratic)

    def evaluate_R(self,R,Z):
        with mp.workdps(self.precision):
            radius=mp.mpf(str(R)); z=mp.mpf(str(Z))
            if not self.comparison.Ra<radius<=self.max_R or abs(z)>=1:
                raise ValueError('Require post-collar R<=max_R and |Z|<1')
            fun=self.functions(z); r=self.jet(radius)
            y=(r/fun['Rb']).log()
            if y.value.constant<0:
                raise ValueError('Radius precedes collar endpoint')
            raw=dict(F=fun['F'].evaluate(y),Uz=fun['Uz'].evaluate(y),R=r)
            raw['moments']={k:v.evaluate(y) for k,v in fun['moments'].items()}
            raw['P']=fun['start']['P']+raw['moments']['p']-fun['start']['moments']['p']
            gy=-self.comparison.W*fun['D'].evaluate(y)/2
            uy=-self.comparison.W*(raw['F']/fun['start']['F'])*(
                fun['start']['F']/self.auxiliary.endpoint(z)['F'])*fun['J'].evaluate(y)
            result={k:v.value for k,v in raw.items() if k!='moments'}
            result.update({k+'_Z':v.tangent for k,v in raw.items() if k!='moments'})
            result['moments']={k:v.value for k,v in raw['moments'].items()}
            result['moments_Z']={k:v.tangent for k,v in raw['moments'].items()}
            quadratics={k:v.evaluate(y) for k,v in fun['raw_quadratic'].items()}
            result['raw_quadratic_integrals']={k:v.value for k,v in quadratics.items()}
            result['raw_quadratic_integrals_Z']={k:v.tangent for k,v in quadratics.items()}
            result['g_y']=gy.value; result['Uz_y']=uy.value
            root=(2*result['R']).sqrt(); F=result['F']
            convert=lambda value:PressureWidthJet(value,pressure_order=self.pressure_order,width_order=self.width_order)
            stress=evaluate_mp_stress(mp.log(radius),z,self.comparison.delta,
                Utheta=root*F,Uz=result['Uz'],Utheta_y=root*F*(mp.mpf('.5')+gy.value),
                Utheta_Z=root*result['F_Z'],Uz_y=uy.value,Uz_Z=result['Uz_Z'],
                moments=result['moments'],moments_Z=result['moments_Z'],P=result['P'],P_Z=result['P_Z'],
                precision=self.precision,radius_override=result['R'],scalar_converter=convert,
                shear_theta=2*F*gy.value,shear_z=root*uy.value/result['R'],include_components=True)
            result.update(Ur=stress['U_r'],Utheta=root*F,stress=stress,
                region='prescribed_post_collar_exit',metadata=self.metadata())
            # Differentiate the actual analytic moment primitive, independently
            # of substituting the moment ODE into the divergence identity.
            mz_y=fun['moments']['z'].derivative().evaluate(y)
            uz_y=fun['Uz'].derivative().evaluate(y)
            dt=self.comparison.delta; L=1-dt*z*z; d=1-z*z
            flux_R=(2*z*result['Uz']+2*z*uz_y.value
                -(1-dt)*z*mz_y.value/result['R']-d*mz_y.tangent/result['R'])/L
            result['Ur_R']=flux_R/root-result['Ur']/(2*result['R'])
            result['analytic_divergence_numerator']=(L*flux_R+d*result['Uz_Z']
                -(1+dt)*z*result['Uz']-2*z*uz_y.value)
            result['moment_radial_identity_defect']=mz_y.value/result['R']-result['Uz']
            result['Uz_y_driver_defect']=uz_y.value-uy.value
            return result

    def metadata(self):
        return dict(frozen_auxiliary_driver=True,prescribed_values_frozen=False,
            positive_epsilon_changes_retained=True,continuation_integrals='analytic exponential polynomials',
            pressure_order=self.pressure_order,width_order=self.width_order,
            post_collar_quadrature_error=False,finite_ring_remainder_enclosed=False,
            collar_RK_error_enclosed=False,functional_terminal_moments_closed=False,
            finite_energy_certified=False,cone_certified=False,global_matching_complete=False)
