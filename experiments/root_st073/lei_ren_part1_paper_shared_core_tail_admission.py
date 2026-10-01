"""Tail admission for the new analytic core in actual scaled field variables.

Phi is the normalized swirl shape. Uz=4Z+j+epsilon*Psi; its tail includes
the actual epsilon prefactor. Raw Psi tails are also reported and are not
asserted small. Radial derivatives use r=Lambda*R, not physical R.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    name='lei_ren_part1_paper_shared_core_majorant.json'
    major=json.loads((HERE/name).read_bytes())
    for path,digest in major['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('majorant dependency changed: '+path)
    if not major['contraction_proved'] or len(major['terms'])!=20:
        raise ValueError('complete new-source contraction required')
    tube_name='lei_ren_part1_paper_shared_analytic_tube.json'
    linear_name='lei_ren_part1_paper_shared_linear_resolvent.json'
    tube=json.loads((HERE/tube_name).read_bytes());linear=json.loads((HERE/linear_name).read_bytes())
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        correction=get(major,'scaled_map_size_upper');eps=get(major,'epsilon')
        phi=get(linear,'Phi_model_Xh_norm_upper')+correction
        psi=get(major,'fresh_Psi_model_Xh_norm_upper')+correction
        uz=eps*psi;h=get(tube,'Xh_parameter')
        target=c.mpf('1e-12');attempts=[]
        for N in range(64,257,2):
            rows=[];maximum=mp.mpf(0);raw_max=mp.mpf(0)
            for total in range(4):
                for i in range(total+1):
                    k=total-i
                    factor=tail_factor(c,degree=N,radial_order=i,axial_order=k,radius='4.1',h=h)['tail_per_Xh_norm']
                    pb=phi*factor;ub=uz*factor;rb=psi*factor
                    maximum=max(maximum,endpoints(pb)[1],endpoints(ub)[1])
                    raw_max=max(raw_max,endpoints(rb)[1])
                    rows.append(dict(scaled_radial_order=i,axial_order=k,
                        Phi_tail=pb,Uz_tail=ub,raw_Psi_tail=rb))
            attempts.append(dict(degree=N,maximum_scaled_field_tail=c.mpf(maximum)))
            if maximum<endpoints(target)[0]:break
        else:raise ValueError('field tail degree not admitted')
        report=dict(implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            analytic_fixed_point_exists_for_new_source=True,
            logLambda=get(major,'logLambda'),logC_definition=major['logC_definition'],
            j=get(major,'required_j'),source_j_eta_tol_relation_verified=True,
            analytic_core_family_sha256=major['analytic_core_family_sha256'],
            full_Section9_parameter_admission=False,actual_delta_definition='min(1e-200,exp(-4logPstar-30))',
            radial_degree=N,ordinary_initial_axial_order=N+3,
            analytic_Phi_Xh_norm_upper=phi,analytic_raw_Psi_Xh_norm_upper=psi,
            analytic_Uz_correction_Xh_norm_upper=uz,Xh_parameter=h,
            epsilon_prefactor_retained=True,
            tail_variables=['Phi','Uz=4Z+j+epsilon*Psi'],
            derivative_coordinates=['scaled_radial=Lambda*R','axial=Z'],
            maximum_mixed_C3_scaled_field_tail=c.mpf(maximum),tail_rows=rows,
            raw_Psi_tail_target_met=raw_max<endpoints(target)[0],
            raw_Psi_tail_is_not_the_axial_velocity_error=True,
            degree_selection=attempts,target=target,target_met=True,
            real_scaled_radial_domain=['0','4.1'],analytic_axis_domain=['-1','1'],
            finite_core_generated=False,old_finite_coefficients_reused=False,
            physical_R_derivative_bounds_certified=False,full_NS_residual_bound=False,
            temporal_recursion=False,
            input_hashes={**major['input_hashes'],name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('New analytic field tail degree',N,'mixed scaled C3 upper',mp.nstr(maximum,14),
            'raw Psi target',report['raw_Psi_tail_target_met'],flush=True)
        return report


if __name__=='__main__':run()
