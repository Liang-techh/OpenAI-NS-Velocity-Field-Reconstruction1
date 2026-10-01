"""Whole-axis analytic cone bounds for the exact Section9 exit prescription.

The field is specified by exact integrations of the admitted analytic core;
this certificate is not a sampled or finite-atlas velocity evaluator. Its
cone is relaxed outside the inner admissible collar. Stress lifting, outer
matching, terminal repairs, and temporal recursion remain separate.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def run():
    names=['shared_K1_ledger','shared_physical_norm_family',
           'shared_frozen_H_bounds','shared_frozen_angular_bounds',
           'shared_core_uniform_bounds','shared_fixed_step_bound']
    data={n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
    hashes={}
    for record in data.values():
        for n,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
                raise ValueError('Global exit dependency changed: '+n)
            hashes[n]=digest
    ledger=data['shared_K1_ledger'];norms=data['shared_physical_norm_family']
    frozen=data['shared_frozen_H_bounds'];angular=data['shared_frozen_angular_bounds']
    core=data['shared_core_uniform_bounds']
    if not (ledger['K1_numeric_bound_certified'] and all(ledger['smallness_checks'].values())
            and ledger['source_j_equals_adopted_eta_tol_over_8']
            and ledger['eta_tol_min_branch_verified']
            and ledger['fixed_K_K1_cstar_and_radius_inner_gates_certified']):
        raise ValueError('Fixed-coefficient parameter gates missing')
    if ledger['uniform_Cstar_family_sha256']!=norms['uniform_Cstar_family_sha256']:
        raise ValueError('Selected analytical family mismatch')
    if any(r['analytic_core_family_sha256']!=ledger['base_analytic_core_family_sha256']
           for r in (frozen,angular,core)):
        raise ValueError('Whole-axis input family mismatch')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        hi=lambda x:c.mpf(endpoints(x)[1]);lo=lambda x:c.mpf(endpoints(x)[0])
        gamma=get(ledger,'gamma');eps0=get(ledger,'epsilon0')
        eta=get(ledger,'eta_tol');K1=c.mpf(ledger['K1']);ch=hi(get(ledger,'cstar'))
        Kmin=c.mpf(1000000)
        small=ledger['decreasing_K_smallness_bounds']
        h_error=get(small,'comparison_H_error')
        point_error=get(small,'actual_error_over_omega')
        comparison_error_times_K=get(small,'comparison_error_times_K')
        # The physical K sum charges both Df and Ef, hence |q_f|<=K.
        # K*|qbar-q_f|<.5 implies |qbar|<=K+.5/K<=2K for K>=1.
        # This stronger pointwise bound is needed for the weak cone; the
        # printed C2 bound 2K^3 would not imply the same coefficient.
        if endpoints(comparison_error_times_K)[1]>=mp.mpf('.5'):
            raise ArithmeticError('Comparison q point-bound implication failed')
        weak_error=get(small,'weak_cone_error')
        strong_ratio=get(small,'strong_cone_error_ratio')
        frozen_D_terminal=lo(get(angular,'frozen_D_R100_R110_lower'))
        if endpoints(frozen_D_terminal)[0]<4:
            raise ValueError('Terminal frozen D input below4')
        # Uniform core budgets hold for the whole increased-Cstar family.
        rho=get(ledger,'rho_core_C2_bound');j=get(ledger,'required_j')
        rho_average=get(ledger,'rho_core_average_C2_bound')
        rho_bridge=3*K1*ch*Kmin**-80
        axial_error=hi(j+rho+rho_bridge)
        inherited_mean_error=hi(j+rho_average)
        # Mz(R)/R=(Ra/R)*Mz(Ra)/Ra+(1/R)*int_Ra^R V(s)ds.
        # The weights are positive, sum to1, and are Z independent.
        # rho_bridge includes both short-switch velocity increments, and
        # their exact integral remains in the mean; moments are never reset.
        mean_error=c.mpf(max(endpoints(inherited_mean_error)[1],endpoints(axial_error)[1]))
        if endpoints(axial_error)[1]>=endpoints(eta*c.mpf('5')/8)[1]:
            # The upper bound j+rho+rho_bridge is <=5eta/8; use an
            # endpoint-safe exact inequality if roundoff straddles equality.
            if endpoints(rho)[1]>endpoints(eta/4)[0] or endpoints(rho_bridge)[1]>endpoints(eta/4)[0]:
                raise ArithmeticError('Axial connection budget failed')
        if (endpoints(axial_error)[1]>=endpoints(eps0)[0]
                or endpoints(mean_error)[1]>=endpoints(eps0)[0]
                or endpoints(axial_error+mean_error)[1]>=endpoints(2*eps0)[0]):
            raise ArithmeticError('Velocity and averaged moment not within9.28 budget')

        # Early core condition plus the actual signed-gradient perturbation.
        signed_gradient_upper=hi(get(core,'signed_Ha_logF_Z_upper')+get(small,'signed_log_gradient_error'))
        if endpoints(signed_gradient_upper)[1]>=mp.mpf('.1'):
            raise ArithmeticError('Signed gradient condition failed')
        # Terminal q is within .5 of the frozen input and actual q within .25.
        comparison_D_floor=c.mpf('3.5')
        actual_D_floor=c.mpf('3.25')
        first_switch_e2_squared_upper=hi(point_error**2)
        first_switch_e2_squared_margin=lo(comparison_D_floor-first_switch_e2_squared_upper)
        if (endpoints(point_error)[1]>=mp.mpf('.25')
                or endpoints(first_switch_e2_squared_margin)[0]<=0):
            raise ArithmeticError('First-switch e2 squared comparison-D gate failed')
        weak_margin=lo(3*gamma-h_error-weak_error)
        strong_margin=lo(1-strong_ratio)
        switch_margin=lo(actual_D_floor-c.mpf('.25')-2)
        if min(endpoints(v)[0] for v in (weak_margin,strong_margin,switch_margin))<=0:
            raise ArithmeticError('Actual cone margin failed')
        # After b becomes zero, V is radial-constant and its averaged axial
        # moment is a convex average. The source term in9.32 stays positive.
        delta_max=c.mpf(1)/200
        SQ=lo(c.mpf('1.8')-c.mpf('7.9')*delta_max-6*eps0-signed_gradient_upper)
        if endpoints(SQ)[0]<mp.mpf('1.4'):
            raise ArithmeticError('Constant-power angular source bound failed')
        barrier=lo(100*SQ-3)  # 1-a/2<=1, R>=100, L<=1
        if endpoints(barrier)[0]<=0:
            raise ArithmeticError('D=3 crossing barrier failed')
        prescription=dict(
            comparison='alpha=1-sigma((y-h)/h); d_y logFbar=alpha*d_y logFc, d_y Vbar=alpha*d_y Vc until2h, then freeze their values; continue all five moments',
            actual='chi=1-(1-epsilon)*sigma(y/h); logF=logf-.5 integralchi*Dbar dy; V=v-integralchi sqrt(R/2)F Ebar dy, Ra<R<=100',
            first_switch='100<=R<=100exp(h): a=epsilon Dbar,b=-epsilon Ebar*(1-sigma(log(R/100)/h))',
            second_switch='100exp(h)<=R<=100exp(2h): b=0,a interpolates epsilon Dbar to4/5 using the same flat step',
            last_interval='a=4/5,b=0 until110; integrate d_y logF=-a/2 andd_y V=b sqrt(2R)F/2',
            recovery='all moments start with the same core data atRa and integrate the actual field; P=P0+Mp; inertial stress is9.13, not reset or fitted',
            parameters='same positive h=epsilon=cstar*K^-100; K is the exact physical9.16 norm sum, with admitted logarithmic bounds')
        proof=dict(
            domain='Z in[-1,1], Ra<R<=110, exact analytic integrations',
            comparison='9.14 inputs, K1 ledger and shared width giveDbar>=1/(2K),|q|<=2K,Hbar>=2+2gamma,Hbar<=K^10; Dbar>=3.5 on100..110',
            early_weight='error e=I/F-q satisfies|e|<=K1K^20(h+epsilon)omega,omega=1-chi>0 for y>0; omega is retained even as y tends to0',
            weak_cone='if kappa=chi Hbar<=2, Hbar+(e dotq)/Dbar >=Hbar-4K^2|e|>2',
            strong_cone='if kappa>2, |q|>=1/(2K),sqrt(1+Hbar)<=2K^5 give sufficient error thresholdomega/(40K^6); the stored strong ratio is below1',
            admissible_collar='t_an=sigma^-1(gamma/(10K^10)) lies in(0,1); Ran=Ra exp(h*t_an); chi Hbar>2 forRa<R<=Ran; nonzero admissible stress there',
            smooth_axis_exit='the cutoff is flat at0 and comparison equals the stress-free core throughh, so actual field matches every core derivative atRa',
            switch_cone='first kappa<=epsilon Hbar<1, D>=3.25, Ebar*E>=-e2^2/4 ande2^2<=Dbar giveD+tEbarE/Dbar>2; second b=0,0<a<=.8 andD>3',
            q_point_bound='the physical K sum contains both Df andEf norms, so|q_f|<=K; K*|qbar-q_f|<.5 implies|qbar|<=K+.5/K<=2K forK>=1',
            first_switch_squared_gate='|e2|<=point_error<1/4 gives e2^2<1/16<3.5<=Dbar; thus Ebar E>=-e2^2/4>=-Dbar/4 throughout the first switch',
            axial_budget='rho_bridge<=3K1*cstar*K^-80<=eta_tol/4 includes both short-switch C2 velocity increments; withj=eta/8 andrho_core<=eta/4, the actual velocity stays within5eta/8<epsilon0',
            inherited_axial_mean='atRa the positive Xh coefficient averaging gives ||Mz/Ra-(4Z+j)||C2<=rho_core_average; its n+1 denominator cannot enlarge the core embedding',
            exact_cumulative_mean='Mz(R)/R=(Ra/R)*Mz(Ra)/Ra+int_Ra^R V(s)ds/R is a positive Z-independent average of the inherited mean and actual velocity, including both short switches; its summed C2 error is at most the maximum of the two stored bounds, no moment reset',
            W_budget='separate velocity and exact cumulative-mean C2 errors are each<epsilon0 and their sum<2epsilon0; formula9.31 then gives W=-3+4delta Z^2+error with C1 error<=10epsilon0',
            gradient='core signedHv partial_Z logf plus the admitted signed perturbation remains<.1; aftertheconstantpowersegment zeta is unchanged inR',
            angular_barrier='withb=0, SQ=-W(1-a/2)-delta(1-2ZV)/2-Hv*zeta>=7/5; d_yD+(1-a/2)D=R SQ/L haspositive derivative atD=3, soD>3 is preserved to110',
            amplitude_budget='log(Cstar*u1*(1+Z^2)) C2<=2A: entrylog(Cstar*f) C2<=A-10, relative exitlog corrections<1, andsqrt220 plus(1+Z^2) logarithms fit the remaining fixed slack',
            limitation='these are analytic bounds for the actual prescribed field, not finite-atlas evaluations or a stress lift; the cone outsideRan is relaxed')
        result=dict(
            admitted_inner_parameter_family_sha256=ledger['admitted_inner_parameter_family_sha256'],
            uniform_Cstar_family_sha256=norms['uniform_Cstar_family_sha256'],
            base_analytic_core_family_sha256=ledger['base_analytic_core_family_sha256'],
            real_axial_domain=['-1','1'],physical_radial_domain=['Ra','110'],
            field_prescription=prescription,proof=proof,
            weak_cone_normalized_lower_margin=weak_margin,strong_cone_sufficient_ratio_margin=strong_margin,
            terminal_comparison_D_floor=comparison_D_floor,short_switch_actual_D_floor=actual_D_floor,
            comparison_q_point_bound_2K_derived=True,
            first_switch_e2_squared_upper=first_switch_e2_squared_upper,
            first_switch_e2_squared_comparison_D_lower_margin=first_switch_e2_squared_margin,
            first_switch_e2_squared_le_comparison_D_certified=True,
            first_switch_relaxed_lower_margin=switch_margin,
            final_constant_power_SQ_lower=SQ,angular_D3_crossing_derivative_lower=barrier,
            rho_core_C2_bound=rho,rho_core_average_C2_bound=rho_average,
            inherited_axial_mean_C2_error_upper=inherited_mean_error,
            rho_bridge_C2_upper=rho_bridge,
            axial_velocity_C2_error_upper=axial_error,axial_cumulative_mean_C2_error_upper=mean_error,
            combined_axial_velocity_and_mean_C2_error_upper=axial_error+mean_error,
            short_switch_moments_included_by_exact_velocity_averaging=True,
            separate_axial_velocity_and_mean_9_28_gates_certified=True,
            axial_velocity_and_average_individual_C2_error_upper=axial_error,
            signed_Hv_logF_Z_upper=signed_gradient_upper,
            actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified=True,
            nonempty_inner_admissible_collar_analytically_certified=True,
            exact_implicit_exit_field_specified=True,whole_axis_finite_velocity_evaluator_built=False,
            finite_atlas_or_sampling_used_as_global_proof=False,
            admissible_stress_lift_constructed=False,whole_transition_admissible_cone=False,
            five_terminal_moment_identities_repaired=False,corrected_outer_at_selected_radius_built=False,
            full_Section9_parameter_admission=False,temporal_recursion=False,
            input_hashes={**hashes,**{PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in names},
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Whole-axis exact exit: analytic relaxed cone throughR110 and inner admissible collar PASS',flush=True)
        print('Constant-power SQ lower=',mp.nstr(endpoints(SQ)[0],14),'D3 barrier=',mp.nstr(endpoints(barrier)[0],14),flush=True)
        return result


if __name__=='__main__':
    run()
