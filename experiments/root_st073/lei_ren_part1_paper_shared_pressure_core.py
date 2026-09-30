"""Actual source-pressure nonlinear local regular core prefix.

Necessary-scale parameters are shared with the exterior. This does not
prove the full core contraction, extend to Ra, join to Rh, or close moments.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets
from lei_ren_part1_paper_core_recursion import core_coefficients,pressure_taylor,core_equation_defects,evaluate_core_jets
from lei_ren_part1_paper_axial_correction import from_signed_log,signed_log


def run(precision=120):
    with mp.workdps(precision):
        lam=mp.mpf(2500); logC=2*mp.log(lam); logP=mp.mpf(14)
        logRref=mp.log(110)+10*(logC+logP)
        schedule=PaperOuterSchedule(logPstar=str(logP),logRref=mp.nstr(logRref,precision),
            delta='1e-32',Md='.5',c_mu='.001',c_delta='.001',c_epsilon='.01')
        profile=CorrectedSourceProfile(schedule=schedule,match_waiting=True,precision=precision)
        @lru_cache(maxsize=128)
        def pressure_float(z):
            data=profile.pressure_at_log_radius(profile.schedule.logRref,z)
            atref=from_signed_log(data['corrected_pressure_nominal'])
            swirl=profile.values(profile.schedule.logRref,z)['Utheta']
            return atref-mp.mpf('2.5')*swirl**2
        def pressure(z):return pressure_float(float(z))
        axis=RegularCoreAxisJets(j='.02',Lambda=lam,logC=logC,
            delta=profile.schedule.delta,pressure=pressure,precision=precision)
        z=mp.mpf('.3'); angular=axis.axis_data(z); K=6
        grad=mp.taylor(lambda w:-lam*(1-axis.delta*w*w)*
            ((1-axis.delta)*w/2+(1-w*w)*(4*w+axis.j))/
            (((1-axis.delta)*w/2+(1-w*w)*(4*w+axis.j))**2+axis.sigma0**2),z,K-2)
        f=[angular['F0']]
        for k in range(K-1):f.append(sum(grad[i]*f[k-i] for i in range(k+1))/(k+1))
        u=[4*z+axis.j,mp.mpf(4)]+[mp.mpf(0)]*(K-2)
        pcoarse=pressure_taylor(pressure,z,spacing='.02',precision=precision)
        pfine=pressure_taylor(pressure,z,spacing='.01',precision=precision)
        sensitivity=[mp.nstr(abs(a-b)/max(abs(b),mp.mpf('1e-100')),30)
                     for a,b in zip(pcoarse[:K],pfine[:K])]
        coefficients=core_coefficients(z,axis.delta,F0_Z_taylor=f,U0_Z_taylor=u,
            P0_Z_taylor=pfine[:K],radial_degree=3,precision=precision)
        samples=[]
        for radius in ('1e-15','5e-16','2.5e-16'):
            defects=core_equation_defects(coefficients,radius)
            values=evaluate_core_jets(coefficients,radius)
            samples.append({'R':radius,
                'defects':{k:signed_log(v,precision) for k,v in defects.items()},
                'angular_defect_over_F0':mp.nstr(abs(defects['angular']/f[0]),40),
                'axial_defect_over_axis_slope':mp.nstr(abs(defects['axial']/coefficients['Uz'][1][0]),40),
                'F_positive':bool(values['F']>0),'F_R_negative':bool(values['F_R']<0),
                'Uz':mp.nstr(values['Uz'],precision)})
        exponents={}
        for key in ('angular_defect_over_F0','axial_defect_over_axis_slope'):
            exponents[key]=mp.nstr(mp.log(mp.mpf(samples[0][key])/mp.mpf(samples[1][key]))/mp.log(2),30)
        first=axis.axis_jets(z)
        slope_agreement={'F_R':mp.nstr(abs(coefficients['F'][1][0]/first.F_R-1),30),
            'Uz_R':mp.nstr(abs(coefficients['Uz'][1][0]/first.Uz_R-1),30)}
        assert all(row['F_positive'] and row['F_R_negative'] for row in samples)
        assert all(mp.mpf(v)>mp.mpf('2.9') for v in exponents.values())
        return {'parameters':{'j':'.02','Lambda':'2500','logCstar':mp.nstr(logC,precision),
                  'logPstar':'14','logRref':mp.nstr(logRref,precision)},
            'actual_axis_pressure':signed_log(pfine[0],precision),
            'axis_data':{k:signed_log(v,precision) for k,v in angular.items()},
            'radial_coefficients':{key:[signed_log(row[0],precision) for row in coefficients[key]]
                                   for key in ('F','Uz','P')},
            'pressure_Taylor_spacing_sensitivity':sensitivity,
            'first_jet_agreement':slope_agreement,'samples':samples,'local_error_exponents':exponents,
            'same_outer_pressure':True,'necessary_reference_scale_relation_satisfied':True,
            'source_absolute_constants_certified':False,'extended_to_Ra':False,
            'matched_to_outer':False,'scale_recursion_established':False,
            'scope':'Third-degree local nonlinear regular-core prefix with actual nominal source pressure; numerical pressure derivatives and heat bounds remain inherited.'}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'local_error_exponents':result['local_error_exponents'],
                      'first_jet_agreement':result['first_jet_agreement']}))
