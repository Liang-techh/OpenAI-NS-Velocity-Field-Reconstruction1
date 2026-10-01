"""Whole-real-axis analytic normalized core bounds, not sampled fields.

The real leading model has chi in [0,1]. Its alternating series supplies
uniform positivity; the fresh fixed-point norm supplies the correction.
Exit input tests follow from coefficient norm embeddings and a signed
square-completion bound. Physical C3 K norms remain separate.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent


def embedding(c, h, r, i, k):
    # Drop (n+1)^-2 from the X_h coefficient bound. Differentiate
    # sum_n binomial(n+k,k)*(r/20)^n=(1-r/20)^(-k-1).
    return (c.mpf(math.factorial(i+k)) / ((k+1)**2 * h**k * 20**i)
            / (1-r/20)**(i+k+1))


def run():
    names = ['lei_ren_part1_paper_shared_core_majorant.json',
             'lei_ren_part1_paper_shared_analytic_tube.json',
             'lei_ren_part1_paper_shared_linear_resolvent.json',
             'lei_ren_part1_paper_shared_bump_constants.json']
    major, tube, linear, fixed = [json.loads((HERE/n).read_bytes()) for n in names]
    for r in (major,tube,fixed):
        for n,d in r['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:
                raise ValueError('Analytic bound dependency changed: '+n)
    if (not major['contraction_proved'] or not tube['common_complex_axis_poles_excluded']
            or major['analytic_core_family_sha256']!=tube['analytic_core_family_sha256']):
        raise ValueError('Fresh analytic admission missing')
    c = MPIntervalContext()
    c.dps = 160
    with mp.workdps(c.dps+40):
        get = lambda r,n: restore_value(c,r[n])
        radius = c.mpf('4.1')
        qmax = radius/2
        # Phi0(q)=sum (-q)^n/[n!(n+1)!], q=chi*r/2 in [0,2.05].
        # Odd cubic is a lower bound: tail terms starting at n=4
        # alternate and decrease. p3'(q)=-1/2+q/6-q^2/48<0 for all q.
        # Even quadratic is an upper bound and <=1 on q<=6.
        model_lower = 1-qmax/2+qmax*qmax/12-qmax**3/144
        tail_start_ratio = qmax/(5*6)
        if not (endpoints(qmax)[1]<3 and endpoints(tail_start_ratio)[1]<1):
            raise ArithmeticError('Alternating model bound not valid')
        h = get(tube,'Xh_parameter')
        eps = get(major,'epsilon')
        correction_norm = get(major,'scaled_map_size_upper')
        phi_norm = get(linear,'Phi_model_Xh_norm_upper')+correction_norm
        psi_norm = get(major,'fresh_Psi_model_Xh_norm_upper')+correction_norm
        correction_value = correction_norm*embedding(c,h,radius,0,0)
        phi_lower = model_lower-correction_value
        phi_upper = 1+correction_value
        if endpoints(phi_lower)[0]<=0:
            raise ArithmeticError('Uniform actual normalized swirl positivity unresolved')
        derivative_rows=[]
        for total in range(4):
            for i in range(total+1):
                k=total-i
                factor=embedding(c,h,radius,i,k)
                derivative_rows.append(dict(scaled_radial_order=i,axial_order=k,
                    Phi_absolute_bound=phi_norm*factor,
                    Phi_minus_model_absolute_bound=correction_norm*factor,
                    Uz_minus_U0_absolute_bound=eps*psi_norm*factor))
        # For Mz/R, coefficient averaging multiplies each radial
        # coefficient by 1/(n+1)<=1. The same positive embedding applies.
        j=get(major,'required_j')
        correction_C2=sum((eps*psi_norm*embedding(c,h,radius,0,k)
                           for k in range(3)),c.mpf(0))
        uz_C2=j+correction_C2
        average_C2=j+correction_C2
        combined=uz_C2+average_C2
        eps0=get(fixed,'epsilon0')
        closeness=endpoints(combined)[1]<endpoints(eps0)[0]
        if not closeness:
            raise ArithmeticError('Two full-axis C2 closeness tests unresolved')
        # Ha=H0+eps*d*Psi, g_axis=L*H0/(H0^2+sigma^2),
        # partial_Z log F=-Lambda*g_axis+Phi_Z/Phi.
        # On the real axis |B'(q)|<=1/2 from its decreasing alternating
        # derivative series; |H0*chi_Z|<=2*20*chi. The model term is
        # therefore <=M*chi. The potentially large axial cross term has
        # size <=B*sqrt(chi), and is controlled by the negative Lambda term.
        sigma=get(tube,'sigma');delta_upper=get(tube,'delta')
        lam=get(major,'Lambda')
        psi_value=psi_norm*embedding(c,h,radius,0,0)
        phi_correction_Z=correction_norm*embedding(c,h,radius,0,1)
        M=20*radius/(2*phi_lower)
        A=lam*(1-delta_upper)-M
        if endpoints(A)[0]<=0:
            raise ArithmeticError('Signed log-gradient negative quadratic not protected')
        B=psi_value/sigma
        square_bound=B*B/(4*A)
        remaining=(get(tube,'H0_modulus_upper')*phi_correction_Z/phi_lower
                   +eps*psi_value*(20*radius/sigma+phi_correction_Z)/phi_lower)
        signed_bound=square_bound+remaining
        signed_pass=endpoints(signed_bound)[1]<mp.mpf(1)/20
        if not signed_pass:
            raise ArithmeticError('Signed Ha partial_Z log F gate unresolved')
        result=dict(
            analytic_core_family_sha256=major['analytic_core_family_sha256'],
            implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            real_axial_domain=['-1','1'],scaled_radial_domain=['0','4.1'],
            leading_model='sum_n (-chi*r/2)^n/[n!(n+1)!], 0<=chi<=1',
            leading_model_proof=dict(odd_cubic='lower bound; p3 derivative negative',
                even_quadratic='upper bound 1-q/2+q^2/12<=1 on q<=6',
                derivative='alternating series: -1/2<=B_prime<=-1/2+q/6<0 on q<3',
                real_chi='H0^2/(H0^2+sigma^2) in [0,1]'),
            leading_model_lower=model_lower,actual_Phi_lower=phi_lower,actual_Phi_upper=phi_upper,
            fixed_point_correction_Xh_norm=correction_norm,
            model_correction_value_bound=correction_value,derivative_bounds=derivative_rows,
            Xh_embedding='(i+k)!/[20^i h^k (k+1)^2 (1-r/20)^(i+k+1)]',
            Uz_minus_4Z_C2_sum_bound=uz_C2,Mz_over_R_minus_4Z_C2_sum_bound=average_C2,
            combined_exit_closeness_C2_bound=combined,epsilon0=eps0,
            C2_convention='sum of separate suprema of derivatives order0,1,2 (no factorial discount)',
            full_real_axis_normalized_swirl_positive=True,
            full_real_axis_Uz_and_average_C2_input_closeness_certified=closeness,
            unique_analytic_field_bounded_not_sampled=True,
            positive_physical_swirl_conditional_on_implicit_positive_F0=True,
            finite_coefficient_state_used=False,whole_axis_finite_velocity_evaluator_built=False,
            signed_Ha_logF_Z_upper=signed_bound,
            signed_square_completion=dict(negative_quadratic_A=A,
                model_chi_coefficient=M,linear_sqrt_chi_B=B,
                square_completion_upper=square_bound,other_correction_upper=remaining),
            signed_Ha_logF_Z_input_test_certified=signed_pass,
            supplied_core_exit_inputs_4_35_and_9_9_certified=True,
            all_supplied_core_exit_inputs_certified=True,
            all_supplied_exit_inputs_certified=False,
            full_physical_C3_K_norms_certified=False,
            frozen_Df_Ef_tests_certified=False,full_Section9_parameter_admission=False,
            physical_R_derivative_bounds_certified=False,temporal_recursion=False,
            proof='odd cubic model bound, Xh embedding/moment averaging, and -A*chi+B*sqrt(chi)<=B^2/(4A)',
            input_hashes={**major['input_hashes'],**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                for n in names+[Path(__file__).name]}})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Whole-axis analytic Phi lower:',mp.nstr(endpoints(phi_lower)[0],14),
              'combined C2 inlet closeness:',mp.nstr(endpoints(combined)[1],14),
              'signed Ha logF_Z upper:',mp.nstr(endpoints(signed_bound)[1],14))
        return result


if __name__=='__main__':
    run()
