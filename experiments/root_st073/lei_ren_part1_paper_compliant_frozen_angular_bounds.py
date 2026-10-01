"""Whole-axis frozen angular direction bounds from the fresh analytic core.

Proves Df>0 through Ra..110 and Df>=4 through100..110. Does not prove
the separate Hf=(Df^2+Ef^2)/Df gate or the admissible stress lift.
"""
# Recomputed for the distinct compliant pressure source.
# Formula origin: lei_ren_part1_paper_shared_frozen_angular_bounds.py; legacy source/receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=['lei_ren_part1_paper_compliant_core_uniform_bounds.json',
           'lei_ren_part1_paper_compliant_core_transfer.json',
           'lei_ren_part1_paper_shared_analytic_tube.json']
    uniform,major,tube=[json.loads((HERE/n).read_bytes()) for n in names]
    for r in (uniform,major,tube):
        for n,d in r['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:
                raise ValueError('Analytic direction dependency changed: '+n)
    if (not uniform['supplied_core_exit_inputs_4_35_and_9_9_certified'] or
        any(r['analytic_core_family_sha256']!=major['analytic_core_family_sha256'] for r in (uniform,tube))):
        raise ValueError('Fresh core input identity mismatch')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(c.dps+40):
        get=lambda r,n:restore_value(c,r[n])
        lo=lambda x:c.mpf(endpoints(x)[0])
        hi=lambda x:c.mpf(endpoints(x)[1])
        radius=c.mpf('4.1');eps=get(major,'epsilon');dt=hi(get(tube,'delta'))
        Lmin=lo(1-dt);sigma=lo(get(tube,'sigma'));j=hi(get(major,'required_j'))
        phiMin=lo(get(uniform,'actual_Phi_lower'));phiMax=hi(get(uniform,'actual_Phi_upper'))
        h=get(tube,'Xh_parameter');corr=get(major,'scaled_map_size_upper')
        psiNorm=get(major,'fresh_Psi_model_Xh_norm_upper')+corr
        psi=hi(psiNorm*embedding(c,h,radius,0,0))
        psiZ=hi(psiNorm*embedding(c,h,radius,0,1))
        correctionZ=hi(corr*embedding(c,h,radius,0,1))
        Hmax=hi(get(tube,'H0_modulus_upper'))
        betaMin=lo(3-dt/2-j)
        # Move the W*r*Phi_r term to the left. Everything below bounds
        # the forcing without using the sign of Phi_r.
        A=phiMin-eps*10*radius/Lmin
        if endpoints(A)[0]<=0:raise ArithmeticError('Angular negative quadratic lost')
        B=eps*psi*phiMax/sigma
        square=B*B/(4*A)
        average_error=eps**2*phiMax*((1+dt)*psi+psiZ)/Lmin
        derivative_error=eps*Hmax*correctionZ/Lmin
        axial_gradient_error=eps**2*psi*(20*radius/sigma+correctionZ)/Lmin
        error=square+average_error+derivative_error+axial_gradient_error
        margin=lo(betaMin*phiMin-error/eps)
        Wmax=hi(3+4*dt+j+eps*(psi+psiZ))
        odeCoefficient=4+eps*Wmax*radius/Lmin
        odeCoefficientLower=4-eps*Wmax*radius/Lmin
        if not (endpoints(margin)[0]>0 and endpoints(Wmax)[1]<4 and
                endpoints(odeCoefficientLower)[0]>0 and
                endpoints(odeCoefficient)[1]<5 and
                endpoints(betaMin/4-margin/5)[0]>0):
            raise ArithmeticError('Uniform angular derivative barrier failed')
        derivativeFloor=eps*margin/5
        Qentry=lo(2*Lmin*margin/(5*phiMax))
        signed=hi(get(uniform,'signed_Ha_logF_Z_upper'))
        S0=lo(3-dt/2-j-eps*(psi+psiZ)-signed)
        # The constant j cancels in mu/Ra=eps*(average Psi-Psi).
        muBound=hi(2*eps*(psi+psiZ))
        Qreference=lo((S0-muBound)/2)
        Qglobal=c.mpf(min(endpoints(Qentry)[0],endpoints(Qreference)[0]))
        thetaUpper=hi(4*eps/100)
        Qterminal=lo(Qreference*(1-thetaUpper**2))
        Dterminal=lo(100*Qterminal)
        Dglobal=4*eps*Qglobal
        if not (endpoints(Qglobal)[0]>0 and endpoints(Dglobal)[0]>0 and
                endpoints(Dterminal)[0]>=4):
            raise ArithmeticError('Frozen angular gates failed')
        result=dict(
            analytic_core_family_sha256=major['analytic_core_family_sha256'],
            implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            real_axial_domain=['-1','1'],core_scaled_radial_domain=['0','4.1'],
            frozen_physical_radial_domain=['Ra=4/Lambda','110'],
            force_error=dict(square_completion=square,average=average_error,
                Phi_Z_correction=derivative_error,axial_gradient=axial_gradient_error),
            negative_force_margin=margin,W_absolute_bound=Wmax,
            derivative_ode_coefficient_upper=odeCoefficient,
            derivative_ode_coefficient_lower=odeCoefficientLower,
            negative_Phi_scaled_radial_derivative_lower=derivativeFloor,
            frozen_Q_entry_lower=Qentry,frozen_S0_lower=S0,
            normalized_entry_moment_error_bound=muBound,frozen_Q_reference_lower=Qreference,
            normalized_entry_moment_error_units='pointwise |A(mu_z)/Ra|; mu_z=Mz(Ra)-Ra*v',
            frozen_Q_whole_domain_lower=Qglobal,frozen_D_whole_domain_lower=Dglobal,
            terminal_theta_upper=thetaUpper,frozen_Q_R100_R110_lower=Qterminal,
            frozen_D_R100_R110_lower=Dterminal,
            frozen_D_positive_full_axis_Ra_R110_certified=True,
            frozen_D_at_least4_full_axis_R100_R110_certified=True,
            proof=dict(
                angular_ODE='2r Phi_rr+[4-eps W r/L]Phi_r=forcing',
                force='forcing<=-eps*margin from -A*chi+B*sqrt(chi)<=B^2/(4A)',
                initial='-Phi_r(0)=(chi+eps beta)/4 >=eps betaMin/4',
                barrier='y=-Phi_r; 2r y_prime+a(r)y>=eps margin, a(r)<=5; y>=eps margin/5',
                stress_free_entry='Qentry=-2 L Phi_r/(eps Phi), Dentry=-8 Phi_r/Phi at rscaled=4',
                frozen_identity='Qf=theta^2 Qentry+(S0/2)(1-theta^2)+(A mu/Ra)theta(1-theta)',
                frozen_lower='theta(1-theta)<=.5(1-theta^2); Qf>=theta^2 Qentry+(1-theta^2)Qreference',
                units='theta=Ra/R, Qf=L Df/R, Df>=R Qf since 0<L<=1'),
            actual_Phi_radial_monotonicity_certified=True,
            finite_coefficient_state_used=False,frozen_E_or_H_test_certified=False,
            full_frozen_profile_test_9_14_certified=False,
            admissible_stress_lift_constructed=False,full_Section9_parameter_admission=False,
            full_physical_C3_K_norms_certified=False,temporal_recursion=False,
            input_hashes={**uniform['input_hashes'],**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                for n in names+[Path(__file__).name]}})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Whole-axis frozen angular gates: Q lower',mp.nstr(endpoints(Qglobal)[0],14),
              'D100..110 lower',mp.nstr(endpoints(Dterminal)[0],14),
              'full H gate remains open')
        return result


if __name__=='__main__':run()
