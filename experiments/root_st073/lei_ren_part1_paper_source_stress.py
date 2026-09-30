"""Same-candidate source stress, with explicit cumulative-moment provider.

The provider defines the supported radial domain; queries outside it fail.
Numerical derivatives and pressure heat bounds remain separate uncertainty.
"""
from decimal import Decimal
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


class SourceStress:
    def __init__(self, profile, moments):
        self.profile=profile
        self.moments=moments
        self.precision=profile.precision

    def evaluate(self, logR, Z, *, step='0.0001'):
        p=self.profile
        logR=Decimal(str(logR))
        with mp.workdps(self.precision):
            h=mp.mpf(str(step)); z=mp.mpf(str(Z))
            if h<=0 or abs(z)+2*h>=1:
                raise ValueError('Stencil must remain inside |Z|<1')
            values=p.values(logR,Z)
            def derivative(fn):
                return (fn(-2)-8*fn(-1)+8*fn(1)-fn(2))/(12*h)
            radial={k:derivative(lambda i:p.values(p.log_at(logR,i*h),Z)[k])
                    for k in ('Utheta','Uz')}
            axial={k:derivative(lambda i:p.values(logR,float(z+i*h))[k])
                   for k in ('Utheta','Uz')}
            m=self.moments(logR,Z)
            mz={k:derivative(lambda i:self.moments(logR,float(z+i*h))[k])
                for k in m}
            def pressure_at(zq):
                receipt=p.pressure_at_log_radius(logR,zq)
                return from_signed_log(receipt['corrected_pressure_nominal'])
            pressure_receipt=p.pressure_at_log_radius(logR,Z)
            pressure=from_signed_log(pressure_receipt['corrected_pressure_nominal'])
            pressure_Z=derivative(lambda i:pressure_at(float(z+i*h)))
            stress=evaluate_mp_stress(logR,Z,p.schedule.delta,
                Utheta=values['Utheta'],Uz=values['Uz'],
                Utheta_y=radial['Utheta'],Utheta_Z=axial['Utheta'],
                Uz_y=radial['Uz'],Uz_Z=axial['Uz'],moments=m,
                moments_Z=mz,P=pressure,P_Z=pressure_Z,precision=self.precision)
            return {'stress':stress,'moments':m,'moments_Z':mz,
                'pressure_receipt':pressure_receipt,'derivative_step':str(step),
                'stress_cone_certified':False,'global_moments_closed':False,
                'PDE_recursion_established':False}

    def receipt(self, logR, Z, *, step='0.0001'):
        data=self.evaluate(logR,Z,step=step)
        with mp.workdps(self.precision):
            return {'logR':str(logR),'Z':str(Z),
                'stress':{k:signed_log(v,self.precision) for k,v in data['stress'].items()},
                'moments':{k:signed_log(v,self.precision) for k,v in data['moments'].items()},
                'derivative_step':data['derivative_step'],
                'heat_pressure_bound':data['pressure_receipt']['heat_deficit_bound_in_P'],
                'stress_cone_certified':False,'global_moments_closed':False,
                'PDE_recursion_established':False}


def run():
    from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
    from lei_ren_part1_paper_reference_moments import PaperReferenceMoments
    profile=CorrectedSourceProfile()
    stress=SourceStress(profile,PaperReferenceMoments(profile,quadrature_order=128))
    points=[]
    with mp.workdps(profile.precision):
        for offset in ('-1','.5','2'):
            logR=profile.log_at(profile.schedule.logRref,offset)
            a=stress.evaluate(logR,.3)
            b=stress.evaluate(logR,.3,step='.00005')
            refinement={k:mp.nstr(abs(a['stress'][k]-b['stress'][k])/
                max(abs(b['stress'][k]),mp.mpf('1e-100')),30)
                for k in ('I_theta','I_z','T_theta','T_z')}
            h=mp.mpf('.001')
            neighbors={i:stress.evaluate(profile.log_at(logR,i*h),.3)['stress']
                       for i in (-2,-1,1,2)}
            identities={}
            for ik,nk,factor in (('I_theta','N_theta',1),('I_z','N_z',mp.mpf('.5'))):
                dy=(neighbors[-2][ik]-8*neighbors[-1][ik]
                    +8*neighbors[1][ik]-neighbors[2][ik])/(12*h)
                defect=dy+factor*a['stress'][ik]-a['stress'][nk]
                identities[ik]=mp.nstr(abs(defect)/max(abs(a['stress'][nk]),
                                                         mp.mpf('1e-100')),30)
            points.append({'offset':offset,'logR':str(logR),'Z':'.3',
                'stress':{k:signed_log(v,profile.precision) for k,v in a['stress'].items()},
                'Z_stencil_refinement':refinement,
                'radial_inertial_identity_relative_defect':identities,
                'heat_pressure_bound':a['pressure_receipt']['heat_deficit_bound_in_P']})
        checks=all(mp.mpf(v)<mp.mpf('1e-6') for row in points
            for v in row['radial_inertial_identity_relative_defect'].values())
    return {'points':points,'radial_identity_checks_passed':checks,
        'scope':'Actual same-candidate moments and pressure through initial axial turnoff only.',
        'stress_cone_certified':False,'global_moments_closed':False,
        'PDE_recursion_established':False}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'radial_identity_checks_passed':result['radial_identity_checks_passed'],
        'identities':[p['radial_inertial_identity_relative_defect'] for p in result['points']]}))
