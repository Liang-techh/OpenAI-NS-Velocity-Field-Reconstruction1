"""Finite sampled residual decomposition, not a PDE or flatness certificate."""
import json
import math
from pathlib import Path
import numpy as np
from lei_ren_part1_physical_remainder import PhysicalRemainder


def run():
    diagnostic=PhysicalRemainder(); field=diagnostic.field
    rows=[]
    for tau in (1/8,1/32,1/128):
        for R,Z in ((.002,.15),(1.,-.3)):
            r,z=field.from_similarity(R,Z,tau)
            point=np.array([r,0.,z])
            measurements=[]
            for fraction in (2e-4,1e-4):
                measurements.append(diagnostic.evaluate(point,tau,
                    spatial_step=min(math.sqrt(field.nu*tau)*fraction,r*fraction),
                    tau_step=tau*fraction))
            rows.append({'tau':tau,'R':R,'Z':Z,'physical_point':point.tolist(),
                         'leading_momentum_scale':math.sqrt(field.nu)*(tau/(1-Z*Z))**(-1.5-field.h),
                         'measurements':measurements,
                         'refinement_max_abs_difference':{key:float(np.max(abs(
                             np.array(measurements[0][key])-np.array(measurements[1][key]))))
                             for key in ('R_B','D_T_B','E_B','axial_viscosity_residual','E_B_without_axial_viscosity')}})
            print(json.dumps({'tau':tau,'R':R,'E_B':measurements[-1]['E_B']}),flush=True)
    fits=[]
    for R,Z in ((.002,.15),(1.,-.3)):
        selected=[row for row in rows if row['R']==R and row['Z']==Z]
        for key in ('R_B','D_T_B','E_B','axial_viscosity_residual','E_B_without_axial_viscosity'):
            norms=[float(np.linalg.norm(row['measurements'][-1][key])) for row in selected]
            normalized=[norm/row['leading_momentum_scale'] for norm,row in zip(norms,selected)]
            fits.append({'R':R,'Z':Z,'quantity':key,
                         'raw_norm_tau_exponent':float(np.polyfit(np.log([row['tau'] for row in selected]),np.log(norms),1)[0]),
                         'norm_divided_by_leading_momentum_scale':normalized,
                         'normalized_tau_exponent':float(np.polyfit(np.log([row['tau'] for row in selected]),np.log(normalized),1)[0])})
    report={'candidate':field.metadata(),'rows':rows,'finite_scale_fits':fits,
            'sampled_component_maxima':{key:max(max(abs(x) for x in row['measurements'][-1][key])
                                              for row in rows) for key in ('R_B','D_T_B','E_B','axial_viscosity_residual','E_B_without_axial_viscosity')},
            'radial_identity_passed':all(row['measurements'][-1]['radial_remainder_equals_radial_residual']
                                         for row in rows),
            'scope':'Six fixed-sector points and two finite-difference steps. Maxima are sampled only; no volume L2, global norms, flatness, cone or corrected-field claim.',
            'volume_L2_measured':False,'flatness_established':False,
            'axial_viscosity_breakdown_is_not_a_field_correction':True,
            'whole_cone_validated':False,'full_PDE_gate_passed':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report['sampled_component_maxima']),flush=True)
    return report


if __name__=='__main__': run()
