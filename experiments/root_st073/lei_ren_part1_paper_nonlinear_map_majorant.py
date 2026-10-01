"""Termwise size/Lipschitz majorants for the full coupled paper (8.50).

All pressure and swirl terms are retained. Bounds are conservative and
conditional on accepted stored data; failure of an upper-bound gate is not
a proof of nonexistence or divergence.
"""
import hashlib
import json
import argparse
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run(Lambda='1e36',output_name=None):
    base=Path(__file__).parent
    names=['lei_ren_part1_paper_analytic_radial_tail.json',
           'lei_ren_part1_paper_linear_resolvent_bound.json',
           'lei_ren_part1_paper_complex_pressure_bound.json',
           'lei_ren_part1_paper_commuting_resolvent_bound.json']
    tube,linear,pressure,commuting=[json.loads((base/n).read_text()) for n in names]
    ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def get(record,name):
            row=record[name]
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        if not linear['analytic_linear_inverse_bound_certified'] or not pressure['all_true_pressure_stages_included']:
            raise AssertionError('Analytic input gates missing')
        eta=get(tube,'complex_tube_radius');h=get(tube,'Xh_parameter');a=1+eta
        dt=ctx.mpf('1e-200');lam=ctx.mpf(Lambda);eps=1/lam;j=ctx.mpf('1e-14')
        if endpoints(lam)[0]<500:raise ValueError('Lambda must satisfy the model gate')
        L=get(tube,'L_modulus_lower');pole=get(tube,'denominator_factor_modulus_lower')
        H=get(tube,'H0_modulus_upper');U=4*a+j;d=1+a*a
        weight=get(linear,'Cauchy_weight_supremum_upper')
        Bphi=get(linear,'Phi_model_Xh_norm_upper')+1
        Bpsi=get(pressure,'Psi_model_Xh_norm_upper')+1
        if commuting['input_sha256']!=hashlib.sha256((base/names[0]).read_bytes()).hexdigest():
            raise AssertionError('Commuting inverse uses different analytic tube')
        Rnorm=get(commuting,'resolvent_norm_upper')
        fixed_multiplier=get(commuting,'axial_convolution_factor_upper')
        # Fixed analytic coefficient norms follow by Cauchy on eta/2 disks.
        coefficients=dict(beta=get(tube,'beta_modulus_upper')*weight,
            W0_over_L=(1+(1+dt)*a*U+4*d)/L*weight,
            H0_over_L=H/L*weight,az=(1+dt)*a/L*weight,
            d_over_L=d/L*weight,z_over_L=a/L*weight,
            axial_linear=((1+dt)/2*(1+4*a*U)+4*d)/L*weight,
            cross_swirl=d/pole*weight)
        Gupper=get(tube,'G_modulus_upper');logC=ctx.mpf('5e151')
        Cstar_margin=logC-2*ctx.log(lam)-lam*Gupper
        if endpoints(Cstar_margin)[0]<=0:raise AssertionError('Candidate fails the complex Cstar guard')
        Fsquare=ctx.exp(2*(-logC+lam*Gupper))*weight
        product=ctx.mpf(256);J1=ctx.mpf(80);J2=ctx.mpf(40)
        axial=ctx.mpf(20480)/h;radial=ctx.mpf(20480)
        Pcal=80*product**2*Fsquare
        rows=[]
        def term(component,name,c,p,q):
            # All displayed formulas have one outer fixed axial multiplier.
            # Cauchy convolution bounds it by S(r), rather than generic256.
            # The internal products of the unknown fields retain256.
            c=c*fixed_multiplier/product
            # c phi^p psi^q is a positive norm majorant after radial inversion.
            value=c*Bphi**p*Bpsi**q
            lip=ctx.mpf(0)
            if p:lip+=c*p*Bphi**(p-1)*Bpsi**q
            if q:lip+=c*q*Bphi**p*Bpsi**(q-1)
            rows.append(dict(component=component,term=name,phi_power=p,psi_power=q,
                             coefficient_upper=c,size_upper=value,Lipschitz_upper=lip))
        # J2 Etheta; combine the fixed non-derivative terms as -beta Phi.
        term('theta','-beta Phi',product*coefficients['beta']*J2,1,0)
        term('theta','W0/L scaledR Phi_R',product*coefficients['W0_over_L']*radial,1,0)
        term('theta','H0/L Phi_Z',product*coefficients['H0_over_L']*axial,1,0)
        term('theta','-eps az M(Psi) Phi',eps*product*coefficients['az']*J2*product,1,1)
        term('theta','-eps az M(Psi) scaledR Phi_R',eps*product*coefficients['az']*radial,1,1)
        term('theta','-eps d/L d_Z M(Psi) Phi',eps*product*coefficients['d_over_L']*axial,1,1)
        term('theta','-eps d/L d_Z M(Psi) scaledR Phi_R',eps*product*coefficients['d_over_L']*axial,1,1)
        term('theta','-eps delta z/L Psi Phi',eps*dt*product*coefficients['z_over_L']*J2*product,1,1)
        term('theta','eps d/L Psi Phi_Z',eps*product*coefficients['d_over_L']*axial,1,1)
        term('theta','-d H0/(H0^2+sigma^2) Psi Phi',product*coefficients['cross_swirl']*J2*product,1,1)
        # J1 Ez; Pcal = F0^2 V(Phi^2), never replaced by an independent datum.
        term('z','W0/L scaledR Psi_R',product*coefficients['W0_over_L']*radial,0,1)
        term('z','axial linear coefficient Psi',product*coefficients['axial_linear']*J1,0,1)
        term('z','H0/L Psi_Z',product*coefficients['H0_over_L']*axial,0,1)
        term('z','-eps az M(Psi) scaledR Psi_R',eps*product*coefficients['az']*radial,0,2)
        term('z','-eps d/L d_Z M(Psi) scaledR Psi_R',eps*product*coefficients['d_over_L']*axial,0,2)
        term('z','-eps (1+delta) z/L Psi^2',eps*(1+dt)*product*coefficients['z_over_L']*J1*product,0,2)
        term('z','eps d/L Psi Psi_Z',eps*product*coefficients['d_over_L']*axial,0,2)
        term('z','d/L d_Z Pcal',product*coefficients['d_over_L']*axial*Pcal,2,0)
        term('z','-2(1+delta)z/L Pcal',2*(1+dt)*product*coefficients['z_over_L']*J1*Pcal,2,0)
        term('z','-2z/L scaledR F0^2 Phi^2',2*product*coefficients['z_over_L']*J1*80*product**2*Fsquare,2,0)
        theta_size=sum((r['size_upper'] for r in rows if r['component']=='theta'),ctx.mpf(0))
        theta_lip=sum((r['Lipschitz_upper'] for r in rows if r['component']=='theta'),ctx.mpf(0))
        z_size=sum((r['size_upper'] for r in rows if r['component']=='z'),ctx.mpf(0))
        z_lip=sum((r['Lipschitz_upper'] for r in rows if r['component']=='z'),ctx.mpf(0))
        size=(Rnorm*theta_size+z_size)/2;lip=(Rnorm*theta_lip+z_lip)/2
        map_size=eps*size;map_lip=eps*lip
        size_gate=endpoints(map_size)[1]<=mp.mpf('.5')
        lip_gate=endpoints(map_lip)[1]<=mp.mpf('.5')
        report=dict(input_hashes={n:hashlib.sha256((base/n).read_bytes()).hexdigest() for n in names},
            accepted_schedule_sha256=pressure['accepted_schedule_sha256'],precision=160,
            Lambda=lam,epsilon=eps,Cstar_log_margin=Cstar_margin,
            pressure_datum_held_fixed=True,
            candidate_requires_shared_core_regeneration=Lambda!='1e36',
            ball_center='paper X0=(Phi0,Psi0)',ball_radius=1,
            Phi_ball_norm_upper=Bphi,Psi_ball_norm_upper=Bpsi,
            product_constant=product,J1_norm_upper=J1,J2_norm_upper=J2,
            fixed_analytic_multiplier_factor_upper=fixed_multiplier,
            angular_resolvent_norm_upper=Rnorm,
            integrated_axial_derivative_constant=axial,integrated_radial_derivative_constant=radial,
            fixed_coefficient_norm_upper=coefficients,F0_squared_Xh_norm_upper=Fsquare,
            restored_pressure_coefficient_norm_upper=Pcal,terms=rows,
            nonlinear_map_size_upper=size,nonlinear_map_Lipschitz_upper=lip,
            scaled_map_size_upper=map_size,scaled_map_Lipschitz_upper=map_lip,
            scaled_map_size_log_upper=ctx.log(map_size),scaled_map_Lipschitz_log_upper=ctx.log(map_lip),
            size_gate_proved=size_gate,Lipschitz_gate_proved=lip_gate,
            contraction_proved=size_gate and lip_gate,
            analytic_fixed_point_exists_for_fixed_datum=size_gate and lip_gate,
            all_paper_8_50_terms_included=True,pressure_and_swirl_couplings_retained=True,
            bounds_relative_to_stored_accepted_parameters=True,
            failed_upper_bound_gate_is_not_nonexistence=True,
            infinite_core_remainder_enclosed=False,original_parameter_errors_enclosed=False,
            temporal_recursion=False,
            next_dependency=('Regenerate shared axis/gauge core at the candidate, verify Ra pressure conventions and apply analytic tails before rebuilding matching layers'
                             if size_gate and lip_gate else
                             'Sharpen preconditioned angular inverse and fixed-data multiplier bounds; do not infer divergence from these upper bounds'))
        output=Path(__file__).with_suffix('.json') if output_name is None else base/output_name
        output.write_text(json.dumps(encode(report),indent=2)+'\n')
        print('all20 terms included; scaled size log upper',mp.nstr(endpoints(ctx.log(map_size))[1],20),
              'scaled Lipschitz log upper',mp.nstr(endpoints(ctx.log(map_lip))[1],20),
              'contraction proved',size_gate and lip_gate,flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--Lambda',default='1e36')
    parser.add_argument('--output-name')
    args=parser.parse_args()
    if args.Lambda!='1e36' and args.output_name is None:
        parser.error('A candidate Lambda needs a separate output-name to preserve the current receipt')
    run(args.Lambda,args.output_name)
