"""Same-source uniform normalized C1 preheat pressure constant.

This Kp is the Section 9 pressure constant, not the unrelated pulse
integral called K_p. It does not certify corrected heat pressure.
"""
# Recomputed for the distinct compliant pressure source.
# Formula origin: lei_ren_part1_paper_shared_pressure_Kp.py; legacy source/receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum as LogarithmicPressureDatum
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent

def run():
    datum=LogarithmicPressureDatum('40');c=datum.ctx
    with mp.workdps(c.dps+40):
        # All Md>=2 in the declared family have yd>=exp(2)+11,
        # outer epsilon<=.001 exp(-4(exp(2)+11)-30), even without the
        # min(1e-200,...) clamp. Use this larger conservative bound.
        epsilon_max=c.mpf('.001')*c.exp(-4*(c.exp(2)+11)-30)
        T=c.exp(c.mpf('.6')-c.exp(2)-11)/(2*(1-epsilon_max)**2)
        slope=datum.stages['slope_transition_ref']['mass']
        m2=c.mpf('2.5')+slope+c.exp(c.mpf('-.4'))/2+3*T
        m0=7*T
        flatten=T/(1-c.mpf('.25')**2)**2
        # On real [-1,1], q^-2<=1, |d(q^-2)/dZ|<=4.
        # Flatten derivative uses its certified complex Cauchy radius1/4.
        value=m2+m0+T;derivative=4*m2+4*flatten
        norm=value+derivative
        Kp=max(1,int(mp.ceil(endpoints(norm)[1])))
        epsilon0=1/(c.mpf(10)**6*(1+1000*(1+Kp)))
        result=dict(implicit_source_sha256=datum.source_sha,datum_enclosure_sha256=datum.datum_sha,
            pressure_Kp=Kp,uniform_preheat_C1_normalized_bound=norm,
            normalized_value_bound=value,normalized_first_Z_derivative_bound=derivative,
            uniform_Md_lower=2,real_axial_domain=['-1','1'],
            pressure_units='P0/Pstar^2',C1_convention='sup|f|+sup|f_Z| (also bounds max convention)',
            universal_post_yd_atom_upper=T,uniform_outer_epsilon_upper=epsilon_max,reference_slope_atom=slope,
            KN=1000*(1+Kp),epsilon0=epsilon0,eta_tol_pressure_cap=epsilon0/4,
            same_source_preheat_pressure_Kp_certified=True,
            uniform_over_declared_logarithmic_outer_family=True,
            pulse_integral_Kp_not_used=True,
            corrected_heat_pressure_Kp_certified=False,
            source_j_eta_tol_relation_verified=False,full_Section9_parameter_admission=False,
            input_hashes={**datum.input_hashes,Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('Same-source normalized preheat pressure C1 constant Kp =',Kp)
    return result
if __name__=='__main__':run()
