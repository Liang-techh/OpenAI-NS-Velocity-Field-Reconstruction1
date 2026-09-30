"""Post-collar source exit, retaining a quantified tiny-positive-shear defect.

Frozen actual values and analytic moment continuation are allowed only when
the prescribed positive-epsilon ODE changes are below the working precision.
This does not implement the subsequent R=100..110 shear changes.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


class ExitContinuation:
    def __init__(self,tangent,*,max_R=100):
        self.tangent=tangent; self.core=tangent.core
        self.comparison=tangent.comparison; self.precision=tangent.precision
        self.r=tangent.r; self.hb=tangent.hb
        with mp.workdps(self.precision):self.max_R=mp.mpf(max_R)

    @lru_cache(maxsize=16)
    def endpoint(self,Z):return self.tangent.evaluate(2*self.hb,Z)

    def evaluate(self,y,Z):
        with mp.workdps(self.precision):
            yy=mp.mpf(y); z=mp.mpf(Z); R=self.r*mp.exp(yy)
            if yy<=2*self.hb:return self.tangent.evaluate(yy,z)
            if R>self.max_R*(1+mp.power(10,-self.precision+10)):
                raise ValueError('Step 2 exit currently ends at R=100')
            a=self.endpoint(z); Re=a['R']; deltaR=R-Re
            F=a['F']; FZ=a['FZ']; u=a['Uz']; uZ=a['UZ']
            bar=self.comparison.evaluate(yy,z)
            coefficients=self.comparison.frozen_driver_coefficients(z)
            dc=coefficients['D']; ic=coefficients['I_z']
            dy=mp.log(R/Re)
            d_terms=[dc[1]*deltaR,dc[0]*dy,dc[-1]*(1/Re-1/R)]
            j_terms=[ic[3]*(R*R-Re*Re)/(2*mp.sqrt(2)),
                     ic[1]*deltaR/mp.sqrt(2),ic[-1]*dy/mp.sqrt(2)]
            eps=self.tangent.epsilon
            g_bound=eps*sum(abs(v) for v in d_terms)/2
            u_bound=eps*abs(F/bar['F'])*sum(abs(v) for v in j_terms)*mp.exp(g_bound)
            # Sufficient bounds for freezing values at this finite precision.
            # expm1 retains a positive error even at astronomical exponents.
            f_relative_bound=mp.expm1(g_bound)
            floor=mp.power(10,-self.precision+20)
            if f_relative_bound>=floor or u_bound>=floor*abs(u):
                raise ArithmeticError('Frozen approximation not below working precision; integrate full ODE')
            moments=dict(a['moments']); mz=dict(a['momentsZ'])
            deltaR2=R*R-Re*Re
            f_bound=abs(F)*f_relative_bound
            quadratic_f_bound=2*abs(F)*f_bound+f_bound*f_bound
            moment_bounds={'theta':f_bound*deltaR2,'z':u_bound*deltaR,
                'theta_z':(abs(F)*u_bound+abs(u)*f_bound+f_bound*u_bound)*deltaR2,
                'z_theta':(2*abs(u)*u_bound+u_bound*u_bound)*deltaR+quadratic_f_bound*deltaR2/2,
                'p':quadratic_f_bound*deltaR}
            moments['theta']+=F*deltaR2
            moments['z']+=u*deltaR
            moments['theta_z']+=F*u*deltaR2
            moments['z_theta']+=u*u*deltaR-F*F*deltaR2/2
            moments['p']+=F*F*deltaR
            mz['theta']+=FZ*deltaR2
            mz['z']+=uZ*deltaR
            mz['theta_z']+=(FZ*u+F*uZ)*deltaR2
            mz['z_theta']+=2*u*uZ*deltaR-F*FZ*deltaR2
            mz['p']+=2*F*FZ*deltaR
            P=a['P']+F*F*deltaR; PZ=a['PZ']+2*F*FZ*deltaR
            # Retain the nonzero prescribed shear even though value changes
            # cannot be added to their larger endpoint values at this precision.
            g_y=-eps*bar['D']/2
            u_y=-eps*mp.sqrt(R/2)*(F/bar['F'])*bar['I_z']
            root=mp.sqrt(2*R)
            stress=evaluate_mp_stress(mp.log(R),z,self.core.delta,
                Utheta=root*F,Uz=u,Utheta_y=root*(F/2+F*g_y),Utheta_Z=root*FZ,
                Uz_y=u_y,Uz_Z=uZ,moments=moments,moments_Z=mz,P=P,P_Z=PZ,
                shear_theta=2*F*g_y,shear_z=root*u_y/R,precision=self.precision)
            return {'R':R,'y':yy,'Z':z,'F':F,'FZ':FZ,'Uz':u,'UZ':uZ,
                'Ur':stress['U_r'],'moments':moments,'momentsZ':mz,'P':P,'PZ':PZ,
                'g_y':g_y,'u_y':u_y,'stress':stress,
                'region':'post_collar_exit','log_F_increment_linear':-eps*sum(d_terms)/2,
                'Uz_increment_linear':-eps*(F/bar['F'])*sum(j_terms),
                'F_relative_value_error_bound':f_relative_bound,'Uz_value_error_bound':u_bound,
                'moment_value_error_bounds':moment_bounds,
                'bound_scope':'Conditional on supplied frozen comparison and endpoint; excludes endpoint RK and driver-Z errors',
                'approximation':'Frozen value/moment continuation; prescribed positive shear retained',
                'Z_error_bound_certified':False,'cone_certified':False,
                'outer_matching_complete':False}


def run():
    from lei_ren_part1_paper_exit_comparison import build_comparison
    from lei_ren_part1_paper_exit_tangents import ExitTangents
    from lei_ren_part1_paper_exit_field import cone_receipt
    comparison=build_comparison()
    continuation=ExitContinuation(ExitTangents(comparison,steps=32))
    rows=[]
    with mp.workdps(comparison.precision):
        def signed_log(x):return {'sign':int(mp.sign(x)),'log_abs':mp.nstr(mp.log(abs(x)),35) if x else None}
        for Z in (comparison.axis.Z0,mp.mpf('.3')):
            for radius in ('1','10','100'):
                R=mp.mpf(radius); a=continuation.evaluate(mp.log(R/continuation.r),Z)
                rows.append({'Z':mp.nstr(Z,35),'R':radius,'Ur':signed_log(a['Ur']),
                    'Uz':mp.nstr(a['Uz'],35),'cone':cone_receipt(a),
                    'F_relative_value_error_bound':signed_log(a['F_relative_value_error_bound']),
                    'Uz_value_error_bound':signed_log(a['Uz_value_error_bound']),
                    'log_F_increment_linear':signed_log(a['log_F_increment_linear']),
                    'Uz_increment_linear':signed_log(a['Uz_increment_linear']),
                    'moment_value_error_bounds':{k:signed_log(v) for k,v in a['moment_value_error_bounds'].items()},
                    'moments':{k:signed_log(v) for k,v in a['moments'].items()}})
    receipt={'rows':rows,'domain':'2hb<log(R/Ra), R<=100',
        'approximation':'Analytic frozen moments, positive ODE shear, explicit subprecision value bounds',
        'Z_error_bound_certified':False,'full_source_constants_certified':False,
        'outer_matching_complete':False,'scale_recursion_established':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'rows':len(rows),'relaxed_passes':[r['cone']['relaxed_passed'] for r in rows]},indent=2))
    return receipt


if __name__=='__main__':run()
