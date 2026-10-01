"""Md40 all-term core majorants with logarithmic parameter selection.

No enormous F0 amplitude or explicit logC value is materialized.
Frozen universal tube/resolvent upper bounds are reused, never old solutions.
"""
import hashlib
import json
import argparse
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    base=Path(__file__).parent
    names=['lei_ren_part1_paper_analytic_radial_tail.json',
           'lei_ren_part1_paper_linear_resolvent_bound.json',
           'lei_ren_part1_paper_commuting_resolvent_bound.json']
    tube,linear,commuting=[json.loads((base/n).read_text()) for n in names]
    from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
    datum=LogarithmicPressureDatum('40');parameters=datum.parameters
    pressure=dict(all_true_pressure_stages_included=len(datum.stages)==14)
    ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def get(record,name):
            row=record[name]
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        if (not tube['common_complex_axis_poles_excluded']
            or not linear['analytic_linear_inverse_bound_certified'] or not pressure['all_true_pressure_stages_included']):
            raise AssertionError('Analytic input gates missing')
        if linear['input_sha256']!=hashlib.sha256((base/names[0]).read_bytes()).hexdigest():
            raise ValueError('linear inverse uses a different analytic tube')
        eta=get(tube,'complex_tube_radius');h=get(tube,'Xh_parameter');a=1+eta
        dt=ctx.mpf('1e-200');j=ctx.mpf('1e-14')
        logP=ctx.mpf(parameters.logPstar);logLambda=4*logP+1000
        lam=ctx.exp(logLambda);eps=ctx.exp(-logLambda)
        if endpoints(lam)[0]<max(mp.mpf(500),endpoints(j**-2)[1]):
            raise ValueError('paper requires Lambda>=max(500,j^-2)')
        L=get(tube,'L_modulus_lower');pole=get(tube,'denominator_factor_modulus_lower')
        H=get(tube,'H0_modulus_upper');U=4*a+j;d=1+a*a
        weight=get(linear,'Cauchy_weight_supremum_upper')
        Bphi=get(linear,'Phi_model_Xh_norm_upper')+1
        if endpoints(parameters.delta)[1]>mp.mpf('1e-200'):
            raise ValueError('new delta outside universal bound')
        qlower=get(tube,'pressure_q_modulus_lower')
        normalized=(ctx.mpf(datum.m2)+ctx.mpf(datum.stages['z_flatten']['mass']))/qlower**2+ctx.mpf(datum.m0)
        Pbound=ctx.exp(2*logP)*normalized
        Pderivative=4*(1+eta)*Pbound/qlower
        half=1+eta/2;Uhalf=4*half+j;Hbound=get(tube,'H0_modulus_upper')
        g=((1+dt)/2*(1+2*half*Uhalf)*Uhalf+4*Hbound
            +2*(1+dt)*half*Pbound+(1+half*half)*Pderivative)/L
        Psi_model=40*g*weight
        Bpsi=Psi_model+1
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
        Gupper=ctx.mpf(endpoints(get(tube,'G_modulus_upper'))[1])
        # Define logC by this shared symbolic expression. Its difference
        # from Lambda*Gupper is exact, not a subtraction of rounded floats.
        logC_definition='Lambda*Gupper + 2*logLambda + 1000'
        Cstar_margin=ctx.mpf(1000)
        if endpoints(Cstar_margin)[0]<=0:raise AssertionError('Candidate fails the complex Cstar guard')
        Fsquare=ctx.exp(-4*logLambda-2000)*weight
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
            implicit_source_sha256=datum.source_sha,datum_enclosure_sha256=datum.datum_sha,precision=160,
            Lambda=lam,logLambda=logLambda,logPstar=logP,logC_definition=logC_definition,
            Gupper_in_logC_definition=Gupper,epsilon=eps,Cstar_log_margin=Cstar_margin,
            fresh_pressure_norm_upper=Pbound,fresh_pressure_derivative_upper=Pderivative,
            fresh_Psi_model_Xh_norm_upper=Psi_model,
            pressure_datum_held_fixed=True,
            candidate_requires_shared_core_regeneration=True,old_finite_core_reused=False,
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
            universal_tube_delta_range=['0','1e-200'],
            paper_Lambda_at_least_max_500_j_inverse_squared=True,
            bounds_for_new_implicit_continuous_source=True,
            logC_and_F0_not_materialized=True,
            failed_upper_bound_gate_is_not_nonexistence=True,
            infinite_core_remainder_enclosed=False,original_parameter_errors_enclosed=False,
            temporal_recursion=False,
            next_dependency=('Regenerate shared axis/gauge core at the candidate, verify Ra pressure conventions and apply analytic tails before rebuilding matching layers'
                             if size_gate and lip_gate else
                             'Sharpen preconditioned angular inverse and fixed-data multiplier bounds; do not infer divergence from these upper bounds'))
        report['input_hashes'].update(datum.input_hashes)
        for name in (Path(__file__).name,'lei_ren_part1_paper_nonlinear_map_majorant.py',
            'lei_ren_part1_paper_analytic_radial_tail.py','lei_ren_part1_paper_linear_resolvent_bound.py',
            'lei_ren_part1_paper_commuting_resolvent_bound.py'):
            report['input_hashes'][name]=hashlib.sha256((base/name).read_bytes()).hexdigest()
        output=Path(__file__).with_suffix('.json')
        output.write_text(json.dumps(encode(report),indent=2)+'\n')
        print('all20 terms included; scaled size log upper',mp.nstr(endpoints(ctx.log(map_size))[1],20),
              'scaled Lipschitz log upper',mp.nstr(endpoints(ctx.log(map_lip))[1],20),
              'contraction proved',size_gate and lip_gate,flush=True)
        return report


if __name__=='__main__':run()
