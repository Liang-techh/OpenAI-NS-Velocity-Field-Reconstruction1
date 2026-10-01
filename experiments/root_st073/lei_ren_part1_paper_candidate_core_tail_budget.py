"""Actual conditional infinite-core norm/tail bounds for Lambda1e48.

Uses the certified fixed-datum contraction. Norm tails compare the analytic
solution with its exact Taylor truncation, not any unverified finite adapter.
"""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run(candidate_name='lei_ren_part1_paper_nonlinear_candidate_Lambda48.json', output_name=None):
    base=Path(__file__).parent
    names=[candidate_name,
           'lei_ren_part1_paper_linear_resolvent_bound.json',
           'lei_ren_part1_paper_complex_pressure_bound.json',
           'lei_ren_part1_paper_analytic_radial_tail.json']
    candidate,linear,pressure,tube=[json.loads((base/n).read_text()) for n in names]
    if not candidate['contraction_proved'] or not candidate['analytic_fixed_point_exists_for_fixed_datum']:
        raise AssertionError('Candidate analytic fixed point not established')
    ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def read(record,name):
            row=record[name]
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        for name in names[1:]:
            if candidate['input_hashes'][name]!=hashlib.sha256((base/name).read_bytes()).hexdigest():
                raise AssertionError('Candidate norm dependency changed')
        correction=read(candidate,'scaled_map_size_upper')
        phi_norm=read(linear,'Phi_model_Xh_norm_upper')+correction
        psi_norm=read(pressure,'Psi_model_Xh_norm_upper')+correction
        h=read(tube,'Xh_parameter');radius=ctx.mpf('4.1')
        # Phi-Phi0 has zero axis value. Its value-series starts at n=1.
        # Square radial weights bound that whole series by this geometric sum.
        radial_ratio=radius/20
        zero_axis_value_factor=radial_ratio/(4*(1-radial_ratio))
        correction_value=correction*zero_axis_value_factor
        qmax=radius/2
        reference_lower=1-qmax/2+qmax**2/12-qmax**3/144
        phi_lower=reference_lower-correction_value
        phi_upper=1+correction_value
        if endpoints(phi_lower)[0]<=0:raise AssertionError('Candidate global Phi positivity not proved')
        indices=[(i,total-i) for total in range(4) for i in range(total+1)]
        target=ctx.mpf('1e-12');required=None
        for N in range(20,513):
            worst=ctx.mpf(0)
            for i,k in indices:
                factor=tail_factor(ctx,degree=N,radial_order=i,axial_order=k,radius=radius,h=h)['tail_per_Xh_norm']
                # Psi has the larger bound; enclose both regardless of ordering.
                norm_upper=max(endpoints(phi_norm)[1],endpoints(psi_norm)[1])
                upper=endpoints(factor*ctx.mpf([0,norm_upper]))[1]
                worst=ctx.mpf([0,max(endpoints(worst)[1],upper)])
            if endpoints(worst)[1]<=endpoints(target)[0]:required=N;break
        if required is None:raise AssertionError('Mixed C3 degree search exceeded bound')
        rows=[]
        for N in (20,required):
            for i,k in indices:
                factor=tail_factor(ctx,degree=N,radial_order=i,axial_order=k,radius=radius,h=h)
                rows.append(dict(degree=N,scaled_radial_order=i,axial_order=k,
                    Phi_tail_upper=phi_norm*factor['tail_per_Xh_norm'],
                    Psi_tail_upper=psi_norm*factor['tail_per_Xh_norm'],**factor))
        report=dict(input_hashes={n:hashlib.sha256((base/n).read_bytes()).hexdigest() for n in names},
             accepted_schedule_sha256=candidate['accepted_schedule_sha256'],precision=160,Lambda=candidate['Lambda']['lower'],
             Xh_parameter=h,analytic_correction_norm_upper=correction,
             analytic_Phi_Xh_norm_upper=phi_norm,analytic_Psi_Xh_norm_upper=psi_norm,
             real_domain=dict(scaled_R=['0','4.1'],Z=['-1','1']),
             Bessel_reference_lower=reference_lower,zero_axis_value_factor=zero_axis_value_factor,
             Phi_correction_value_upper=correction_value,Phi_real_lower=phi_lower,Phi_real_upper=phi_upper,
             global_candidate_Phi_positive_relative_to_accepted_data=True,
             infinite_radial_tail_bounds=rows,
             normalized_scaled_mixed_C3_tail_target=target,required_radial_degree=required,
             required_axis_Taylor_length_for_C3=required+4,
             tail_target_is_not_full_NS_residual_target=True,
             exact_Taylor_truncation_tails_enclosed_relative_to_accepted_data=True,
             candidate_finite_coefficients_generated=False,finite_adapter_errors_enclosed=False,
             physical_R_derivative_conversion='multiply radial order i by Lambda^i; Uz correction additionally has epsilon',
             original_parameter_errors_enclosed=False,core_to_collar_matching_certified=False,
             global_stress_cone_certified=False,temporal_recursion=False,
             next_dependency='Generate the coupled candidate coefficients and compare them with these exact Taylor tails; restore moments/pressure and coherent collar data')
        output=base/output_name if output_name else Path(__file__).with_suffix('.json')
        output.write_text(json.dumps(encode(report),indent=2)+'\n')
        print('candidate Phi global lower',mp.nstr(endpoints(phi_lower)[0],22),
              'required radial degree',required,'axis length',required+4,flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--candidate-file',default='lei_ren_part1_paper_nonlinear_candidate_Lambda48.json')
    parser.add_argument('--output-name')
    args=parser.parse_args()
    run(args.candidate_file,args.output_name)
