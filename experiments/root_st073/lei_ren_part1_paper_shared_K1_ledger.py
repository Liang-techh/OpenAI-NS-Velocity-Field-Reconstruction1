"""Constructive fixed coefficients for paper Section 9.4, not fitted norms.

K1 is generated from the fixed radius110, derivative product rules, and
explicit five-moment stress formulas. No j, Lambda, pressure, or input
profile is used to choose K1. Analytic source receipts are bound only after
the integer ledger has been constructed.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def fixed_ledger():
    R = 110
    # Sum C3(F) uses derivatives K,K^2,2K^3,5K^4. Componentwise
    # convolution keeps F^2 in K^5, rather than the naive K^8 bound.
    products = dict(F_C3=9,F_times_V_C3=23,F_squared_C3=31,V_squared_C3=4)
    moments = dict(theta=1+R*R*products['F_C3'],z=1+R,
        theta_z=1+R*R*products['F_times_V_C3'],
        ztheta=1+R*products['V_squared_C3']+(R*R//2)*products['F_squared_C3'],
        p=1+R*products['F_squared_C3'])
    pressure = 1+moments['p']
    B = R+14*moments['z']
    angular = 3*moments['theta']+14*moments['theta_z']
    axial = 2*B+3*moments['z']+11*moments['ztheta']+20*R*pressure
    # Sum C2 product factor2. Each reciprocal F and Linv is <=4K,4.
    # Angular term has scalar1/(2R)<=K/2; axial1/sqrt(2R)<=K.
    Cq = 8*B+32*angular+64*axial
    # Difference fractions use the entire R<=110 domain, including switches.
    # X is the summed C1 difference of F,V; C1 is a Banach algebra.
    moment_differences = dict(theta=R*R,z=R,theta_z=6*R*R,
        ztheta=4*R+4*R*R,p=8*R,P=8*R)
    CM = ((sum(moment_differences.values())+9999)//10000)*10000
    Cabs = 250000
    delta_axial = 334*CM+331
    stress = (6*CM+16*Cabs)+8*delta_axial+96*Cabs
    CX = 256*(1+Cq)
    switch_fields = 128*(1+Cq)
    switch_moments = 2000000
    signed_log = 128*(1+Cq)
    coefficients = dict(comparison_velocity=28,comparison_stress=28*stress,
        comparison_H=700*stress,q_C2=Cq,actual_weighted_velocity=CX,
        actual_weighted_stress=stress*CX,
        stress_after_short_switches=stress*(CX+switch_fields),
        switch_fields=switch_fields,switch_moments=switch_moments,
        signed_log_gradient=signed_log,axial_bridge=CX+switch_fields)
    K1 = 1+sum(coefficients.values())
    return dict(K1=K1,derivative_products=products,
        comparison_moment_C3_coefficients=moments,comparison_pressure_C3_coefficient=pressure,
        stress_C2_numerator_coefficients=dict(B=B,angular=angular,axial=axial),
        moment_difference_C1_radius_upper=R,
        moment_difference_C1_individual_coefficients=moment_differences,
        moment_difference_C1_coefficient=CM,absolute_moment_pressure_C1_coefficient=Cabs,
        normalized_stress_difference_coefficient=stress,
        coefficients=coefficients)


def run():
    ledger = fixed_ledger()  # constructed before reading any input family
    names = ['shared_physical_norm_family','shared_Cstar_envelope_transfer',
             'shared_bump_constants','shared_core_uniform_bounds',
             'shared_frozen_H_bounds','shared_tolerance_exit_parameters',
             'shared_fixed_step_bound']
    records = {n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
    hashes = {}
    for record in records.values():
        for n,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
                raise ValueError('K1 dependency changed: '+n)
            hashes[n]=digest
    norms=records['shared_physical_norm_family']; fixed=records['shared_bump_constants']
    uniform=records['shared_core_uniform_bounds']; frozen=records['shared_frozen_H_bounds']
    transfer=records['shared_Cstar_envelope_transfer']
    if (not records['shared_fixed_step_bound']['global_derivative_bound_certified']
            or records['shared_fixed_step_bound']['global_derivative_upper']!=8):
        raise ValueError('Fixed smooth step hypothesis not certified')
    if not (norms['full_physical_C3_K_norms_certified_for_uniform_analytic_family']
            and norms['radius_9_17_compatibility_certified_for_selected_analytic_family']
            and transfer['larger_Cstar_core_and_local_exit_envelopes_bound']
            and frozen['full_frozen_profile_test_9_14_certified']):
        raise ValueError('Uniform physical norm/radius/input gates missing')
    if transfer['selected_uniform_Cstar_family_sha256']!=norms['uniform_Cstar_family_sha256']:
        raise ValueError('Selected radius/envelope family mismatch')
    base=norms['base_analytic_core_family_sha256']
    if any(r['analytic_core_family_sha256']!=base for r in (uniform,frozen)):
        raise ValueError('Frozen/core input family mismatch')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        hi=lambda x:c.mpf(endpoints(x)[1])
        lo=lambda x:c.mpf(endpoints(x)[0])
        gamma=c.mpf('.01'); eps0=get(fixed,'epsilon0')
        j=get(fixed,'required_j'); eta=get(fixed,'eta_tol')
        if endpoints(8*j)!=endpoints(eta):
            raise ValueError('Source j is not the adopted eta_tol/8')
        e_star=get(fixed,'e_star')
        if not (endpoints(e_star/100)[1]<endpoints(eps0/4)[0]
                and endpoints(eta)==endpoints(e_star/100)
                and records['shared_tolerance_exit_parameters']['definition']['eta_tol']
                    =='min(epsilon0/4,e_star/100)'):
            raise ValueError('Adopted eta_tol min branch not bound')
        a=eps0*gamma**2/1000000
        b=eta/12
        if endpoints(b)[1]>=endpoints(a)[0]:
            raise ValueError('Resolve cstar min branch with tighter arithmetic')
        cstar=b/(2*ledger['K1'])
        if not 0<endpoints(cstar)[0]<=endpoints(cstar)[1]<mp.mpf('.5'):
            raise ValueError('Positive fixed cstar missing')
        Kmin=c.mpf(1000000); ch=hi(cstar)
        coef=ledger['coefficients']; Cq=coef['q_C2']; st=ledger['normalized_stress_difference_coefficient']
        # Every expression decreases with K. Evaluate at the fixed minimum.
        small=dict(comparison_log_C2=4*ch*Kmin**-99,
            actual_log_C2=2*(Cq+1)*ch*Kmin**-92,
            axial_velocity_change=128*Cq*ch*Kmin**-91,
            comparison_error_times_K=coef['comparison_stress']*ch*Kmin**-91,
            comparison_H_error=coef['comparison_H']*ch*Kmin**-88,
            signed_log_gradient_error=2*coef['signed_log_gradient']*ch*Kmin**-90,
            actual_error_over_omega=2*ledger['K1']*ch*Kmin**-80,
            weak_cone_error=8*ledger['K1']*ch*Kmin**-78,
            strong_cone_error_ratio=80*ledger['K1']*ch*Kmin**-74,
            axial_bridge_error_over_eta=(3*ledger['K1']*ch*Kmin**-80)/lo(eta),
            initial_switch_a=2*ch*Kmin**-99,
            small_shear_kappa=ch*Kmin**-90,
            double_h=2*ch*Kmin**-100)
        thresholds=dict(comparison_log_C2='.25',actual_log_C2='.25',
            axial_velocity_change='1',comparison_error_times_K='.5',
            comparison_H_error='.01',signed_log_gradient_error='.05',
            actual_error_over_omega='.25',weak_cone_error='.01',strong_cone_error_ratio='1',
            axial_bridge_error_over_eta='.25',initial_switch_a='.8',small_shear_kappa='1',double_h='.01')
        checks={n:endpoints(v)[1]<endpoints(c.mpf(thresholds[n]))[0] for n,v in small.items()}
        if not all(checks.values()):
            raise ArithmeticError('Fixed coefficient/width smallness failed')
        rho_core=sum((get(row,'Uz_minus_U0_absolute_bound') for row in uniform['derivative_bounds']
                      if row['scaled_radial_order']==0 and row['axial_order']<=2),c.mpf(0))
        if endpoints(rho_core)[1]>endpoints(eta/4)[0]:
            raise ValueError('Full-axis source core axial budget not admitted')
        # Mz/R divides the nth radial coefficient by n+1. The absolute
        # Xh embedding is positive, so it also bounds this exact average.
        rho_average=rho_core
        if endpoints(get(uniform,'Mz_over_R_minus_4Z_C2_sum_bound'))!=endpoints(j+rho_average):
            raise ValueError('Inherited core moment-average embedding changed')
        logC=get(norms,'selected_logCstar'); logKbar=get(norms,'log_Kbar_upper')
        log_h_bounds=c.mpf([endpoints(c.ln(cstar)-100*(logC+logKbar))[0],
                            endpoints(c.ln(cstar)-100*logC)[1]])
        family_definition=dict(uniform_Cstar_family_sha256=norms['uniform_Cstar_family_sha256'],
            K1=ledger['K1'],cstar='eta_tol/(24*K1), admitted min branch',
            K='actual physical norm sum in9.16, bounded by Cstar*Kbar',
            h_b='cstar*K^-100',epsilon_b='h_b, the identical scalar',
            logCstar=norms['selected_logCstar'],logRref='log110+10(logCstar+logPstar)')
        parameter_sha=hashlib.sha256(json.dumps(family_definition,sort_keys=True).encode()).hexdigest()
        proof=dict(
            norms='unweighted sum norms; C1 product factor1, C2 factor2, C3 factor4; a Z derivative consumes one more order',
            comparison='integration by parts gives positive convex weights for log(Cstar F),V in C_Z^3; s<=2 gives relative log and V differences<=4Kh; finite exp derivatives give norm<=2 and difference<=3T for T<=1/4',
            comparison_reciprocal='f and1/f C2<=K; multiplication by exp(+-relative log) gives each comparison F and1/F C2<=4K',
            cubic_products='F derivatives K,K^2,2K^3,5K^4 give FV C3<=23K^5 andF^2 C3<=31K^5; squaring a full C3 bound would incorrectly lose theK^7 target',
            q='Aop:C3->C2 coefficient14 andPop:C3->C2 coefficient20; use inherited moments and current integrals; R^-1<=K,LinvC2<=4,FinvC2<=4K give displayed CqK^7',
            differences='on R<=110 the displayed six difference coefficients sum to135410, rounded up to140000: five moment and pressure C1 differences<=140000K X; absolute C1 moment/pressure sum<=250000K^2; exact9.13 fractions and reciprocal difference<=16K^2X give Cstress K^6 X',
            early_weight='omega nondecreasing and Z independent; ell<=2Kh omega,chi expell-1 C2<=2omega; Fdiff<=48K^2h omega,Vdiff<=8Kh omega',
            later_weight='integralchi<=K(h+epsilon),log110/r<=K; log(F/f)C2<=CqK^8(h+epsilon)/2,V-vC2<=64CqK^9(h+epsilon); omega>=1/2 after h yields CXK^10(h+epsilon)omega',
            actual_stress='the same difference fractions give error<=Cstress CXK^16(h+epsilon)omega<=K1K^20(h+epsilon)omega; the flat omega factor is retained',
            H_transfer='gradient ofD+E^2/D bounded by25K^4 onD>=1/(2K),|D|+|E|<=2K; comparisonH error<=700Cstress K^12h',
            switches='C2 velocity increments<=128(1+Cq)K^8h; physical length of2h switches<=440h gives summed C1 moments andpressure increments<=2000000K^2h',
            signed_gradient='actualHv<=6 after the axial budget; signedHv*zeta changes byHv_actual*Delta_zeta+d*DeltaV*zeta_entry, with|zeta_entry|<=K; coefficient128(1+Cq)K^10(h+epsilon) suffices',
            constant_independence='fixed_ledger reads no profile data; its integers depend only on radius110, derivative product rules, andfixed cutoffs/auxiliary hypotheses',
            inherited_average='Mz/R core coefficient averaging divides by n+1<=1, so the same positive Xh embedding bounds both Uz-U0 andMz/R-U0 in summed C2',
            min_branch='eta_tol is the adopted min(epsilon0/4,e_star/100), equals8j; eta_tol/12 < epsilon0 gamma^2/1e6; cstar=eta_tol/(24K1) satisfies both9.4 restrictions')
        result=dict(K1=ledger['K1'],fixed_coefficient_ledger=ledger,proof=proof,
            K1_numeric_bound_certified=True,K1_chosen_independently_of_input_profiles=True,
            base_analytic_core_family_sha256=base,
            uniform_Cstar_family_sha256=norms['uniform_Cstar_family_sha256'],
            admitted_inner_parameter_family_sha256=parameter_sha,parameter_definition=family_definition,
            cstar=cstar,gamma=gamma,epsilon0=eps0,eta_tol=eta,required_j=j,
            rho_core_C2_bound=rho_core,rho_core_average_C2_bound=rho_average,
            source_j_equals_adopted_eta_tol_over_8=True,eta_tol_min_branch_verified=True,
            shared_positive_width_log_enclosure=log_h_bounds,
            positive_width_not_materialized=True,h_b_equals_epsilon_b_by_definition=True,
            decreasing_K_smallness_bounds=small,smallness_thresholds=thresholds,smallness_checks=checks,
            fixed_K_K1_cstar_and_radius_inner_gates_certified=True,
            global_actual_exit_cone_certified=False,full_Section9_parameter_admission=False,
            corrected_outer_at_selected_radius_built=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,
            input_hashes={**hashes,**{PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in names},
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Fixed input-independent K1 =',ledger['K1'],'Cq =',Cq,flush=True)
        print('PASS:',len(checks),'shared-width smallness gates; full actual exit certificate is separate',flush=True)
        return result


if __name__=='__main__':
    run()
