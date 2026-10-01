"""Same-family Section 9.5 log field and whole-axis analytic join bounds.

The selected bound Abar replaces the exact A in the conservative choice
T=400Abar; F27 imposed the radius restrictions using this same Abar.
Fields remain implicit analytic profiles. No exp(Abar), physical radius,
physical swirl or inherited moment is rounded to zero or reset.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_functional_defects import logjet
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def sigma8(c,s):
    """Same flat step, with the separately proved global derivative bound8."""
    if endpoints(s)[1]<=0:
        return c.mpf(0),c.mpf(0)
    if endpoints(s)[0]>=1:
        return c.mpf(1),c.mpf(0)
    sig=1-alpha_box(c,s+1,c.mpf(1))
    if endpoints(s)[0]<=0 or endpoints(s)[1]>=1:
        return sig,c.mpf([0,8])
    ds=sig*(1-sig)*(2/s**3+2/(1-s)**3)
    return sig,c.mpf([max(mp.mpf(0),endpoints(ds)[0]),min(mp.mpf(8),endpoints(ds)[1])])


def coordinate_square_plus_one(c,z):
    # z is the coordinate jet; square its real interval as a power to
    # preserve 1+Z^2>=1, rather than interval-multiplying two copies.
    if len(z.coefficients)<2 or endpoints(z[1])!=(mp.mpf(1),mp.mpf(1)):
        raise ValueError('The identity Z coordinate jet is required')
    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in z.coefficients[2:]):
        raise ValueError('Nonlinear coordinate jet supplied')
    return IntervalTaylor(c,[1+z[0]**2,2*z[0]]
                          +([c.mpf(1)] if len(z.coefficients)>2 else [])
                          +[c.mpf(0)]*max(0,len(z.coefficients)-3))


def shape_log_jet(c,phase,Abar,logC,logP,z,B):
    """Enclose log(u/Pstar), a and zeta from supplied ordinary inlet jets.

    B is log(Cstar*u1*(1+Z^2)), not an independently fitted profile.
    This function does not supply the inlet jets or recover moments.
    """
    phase=c.mpf(phase)
    if endpoints(phase)[0]<0 or endpoints(phase)[1]>1 or endpoints(Abar)[0]<10:
        raise ValueError('shape phase in[0,1] and Abar>=10 required')
    sig,ds=sigma8(c,phase)
    logu=B*(1-sig)-logjet(coordinate_square_plus_one(c,z))+(40*Abar*phase-logC-logP)
    a=B*(ds/(200*Abar))+c.mpf('.8')
    return dict(log_Utheta_over_Pstar_coefficients=list(logu.coefficients),
                angular_shear_a=a[0],angular_shear_a_coefficients=list(a.coefficients),
                shape_cutoff=sig,shape_cutoff_derivative=ds,
                physical_radius_or_amplitude_materialized=False)


def reference_log_jet(c,offset,z):
    """Reference log(u/Pstar) at the exact offset log(R/Rref)."""
    return -logjet(coordinate_square_plus_one(c,z))+c.mpf(offset)/10


def restore_axial_jet(c,t,z,v1):
    """(9.38); t is the exact relative offset log(R/Rz), not logR-logRz."""
    t=c.mpf(t)
    if endpoints(t)[0]<0 or endpoints(t)[1]>1:
        raise ValueError('restore offset in[0,1] required')
    sig,ds=sigma8(c,t)
    v=v1+(z*4-v1)*sig
    vy=(z*4-v1)*ds
    return dict(Uz_coefficients=list(v.coefficients),Uz_y_coefficients=list(vy.coefficients),
                restoration_cutoff=sig,restoration_cutoff_derivative=ds)


def run():
    names=['shared_global_exit_certificate','shared_K1_ledger',
           'shared_physical_norm_family','shared_pressure_Kp','shared_fixed_step_bound',
           'shared_core_majorant','shared_analytic_tube','shared_core_uniform_bounds',
           'shared_core_seed_check']
    data={n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
    hashes={}
    for record in data.values():
        for n,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
                raise ValueError('Reference join dependency changed: '+n)
            hashes[n]=digest
    exit=data['shared_global_exit_certificate'];ledger=data['shared_K1_ledger']
    norms=data['shared_physical_norm_family'];pressure=data['shared_pressure_Kp']
    if not (exit['actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified']
            and exit['separate_axial_velocity_and_mean_9_28_gates_certified']
            and norms['radius_9_17_compatibility_certified_for_selected_analytic_family']
            and pressure['same_source_preheat_pressure_Kp_certified']
            and data['shared_fixed_step_bound']['global_derivative_upper']==8):
        raise ValueError('Same-family exit/radius/pressure/cutoff inputs required')
    if (exit['admitted_inner_parameter_family_sha256']!=ledger['admitted_inner_parameter_family_sha256']
            or exit['uniform_Cstar_family_sha256']!=norms['uniform_Cstar_family_sha256']
            or norms['implicit_source_sha256']!=pressure['implicit_source_sha256']
            or norms['datum_enclosure_sha256']!=pressure['datum_enclosure_sha256']):
        raise ValueError('Connection family or analytic pressure changed')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        hi=lambda x:c.mpf(endpoints(x)[1]);lo=lambda x:c.mpf(endpoints(x)[0])
        A=get(norms,'A_upper');logC=get(norms,'selected_logCstar')
        logP=get(data['shared_core_majorant'],'logPstar')
        tube=data['shared_analytic_tube'];j=get(tube,'j')
        seed=data['shared_core_seed_check'];dt=get(tube,'delta')
        if (tube['core_family_definition']['G_anchor']!='unique H0 root in [-j,0]'
                or endpoints(j)[1]>=mp.mpf('.001')
                or endpoints(get(data['shared_core_uniform_bounds'],'actual_Phi_upper'))[1]>=2
                or not seed['unique_axis_anchor_root_proved']
                or not seed['seed_acceptance_passed']
                or seed['primitive_branch']!='G(Z0)=0 at unique H root'
                or seed['analytic_core_family_sha256']!=exit['base_analytic_core_family_sha256']):
            raise ValueError('Normalized swirl/root normalization inputs changed')
        half_coeff=lo((1-dt)/2)
        G_sign_gates=dict(L_floor=lo(1-dt),H_negative_region_upper=hi(-half_coeff*j),
            H_positive_region_lower=c.mpf(min(endpoints(c.mpf('.75')*j)[0],endpoints(half_coeff/2)[0])),
            H_middle_derivative_floor=lo((9-dt)/2-14*j*j),
            H_root_left_upper=hi(-(1-dt)*j/2-3*j*(1-j*j)),H_root_right_lower=lo(j))
        if (min(endpoints(G_sign_gates[n])[0] for n in
                ('L_floor','H_positive_region_lower','H_root_right_lower'))<=0
                or endpoints(G_sign_gates['H_middle_derivative_floor'])[0]<=4
                or max(endpoints(G_sign_gates[n])[1] for n in
                    ('H_negative_region_upper','H_root_left_upper'))>=0):
            raise ArithmeticError('Real-axis G sign/branch gates failed')
        eps0=get(ledger,'epsilon0');ch=hi(get(ledger,'cstar'));Kmin=c.mpf(1000000)
        Cq=ledger['fixed_coefficient_ledger']['coefficients']['q_C2']
        actual_log=get(ledger['decreasing_K_smallness_bounds'],'actual_log_C2')
        # Sum of two switches, each with length h and the factor1/2 in
        # d_ylogF=-a/2. Z-independent interpolation has C2 norm<=.8+epsilon*CqK^7.
        switch_log=hi(ch*Kmin**-100*(c.mpf('.8')+Cq*ch*Kmin**-93))
        relative_log=hi(actual_log+switch_log)
        extra=hi(relative_log+c.mpf('.4')*c.ln(c.mpf('1.1'))
                 +c.ln(220)/2+c.ln(2)+1+2)
        if endpoints(extra)[1]>=10 or endpoints(relative_log)[1]>=mp.mpf('.5'):
            raise ArithmeticError('Recomputed normalized R110 B budget failed')
        # Abar includes10+the unweighted core log/velocity C3 upper. Thus
        # ||B||C2<=Abar-10+extra<Abar<=2Abar. The shape uses T=400Abar.
        T=400*A
        delta=c.mpf(1)/200
        signed_entry=get(exit,'signed_Hv_logF_Z_upper')
        signed_shape=hi(signed_entry+4*eps0)
        a_min=c.mpf('.7');a_max=c.mpf('.9')
        W_error=10*eps0
        SQ=lo((3-4*delta-W_error)*(1-a_max/2)-c.mpf('5.5')*delta-signed_shape)
        q_one_unit=lo((c.mpf('1.4')/c.mpf('1.65'))*(1-c.exp(c.mpf('-1.65'))))
        barrier=lo(110*SQ-3)
        if (endpoints(SQ)[0]<=mp.mpf('1.4') or endpoints(q_one_unit)[0]<=mp.mpf('.5')
                or endpoints(barrier)[0]<=0):
            raise ArithmeticError('Long-shape angular source/barrier failed')

        # All large scales remain logarithmic. Derive these expressions
        # from Rref=110(Cstar Pstar)^10, not from a rounded logRref value.
        shape_to_z_gap=lo(10*(logC+logP)-T-8)
        # log_Kbar can itself be too large to exponentiate. Use the
        # established log-sum majorant instead (Kbar>=1 => ln(1+Kbar)<=lnKbar+ln2).
        inherited_radius_margin=lo(8*logC+12*logP-8
                                   -2*(get(norms,'log_Kbar_upper')+c.ln(2)))
        if endpoints(shape_to_z_gap)[0]<=1 or endpoints(inherited_radius_margin)[0]<=0:
            raise ArithmeticError('Shape/reference/restore radius ordering failed')

        early_mp_coefficient=1760*(3+2*A)
        early_mp_exp_margin=lo(4*A-c.ln(early_mp_coefficient))
        early_mp_log_upper=hi(4*A-2*logC)
        theta_increment_log_upper=hi(c.ln(48400)+c.ln(1+A)-logC)
        mixed_increment_log_upper=hi(c.ln(435600)+c.ln(1+A)-logC)
        if (endpoints(early_mp_exp_margin)[0]<=0 or endpoints(early_mp_log_upper)[1]>=0
                or endpoints(theta_increment_log_upper)[1]>=0
                or endpoints(mixed_increment_log_upper)[1]>=0):
            raise ArithmeticError('Relative-amplitude early moment bounds failed')
        Kp=pressure['pressure_Kp'];KN=pressure['KN']
        # The shape pressure integral <=(5+30Abar)u_sh^2, and9.17 gives
        # u_sh^2/Pstar^2<=exp(-2)/(1+Abar)^2. No enormous amplitude is formed.
        shape_pressure=hi(35*c.exp(-2)/(1+A))
        reference_pressure=c.mpf('7.5')
        pressure_C1=hi(Kp+1+shape_pressure+reference_pressure)
        if endpoints(pressure_C1)[1]>=Kp+100:
            raise ArithmeticError('Normalized whole-join C1 pressure bound failed')
        # Each inherited core moment C1<=K; additional Mz,V^2 integrals
        # <=990,8910 and angular/pressure increments<1. Thus moment and
        # pressure C1 norms<=2K. Formula9.13 atR110 then gives N110<=20K.
        inherited_N_over_K=hi(5/Kmin+c.mpf(25)/110+c.mpf('4.02'))
        if endpoints(inherited_N_over_K)[1]>=20:
            raise ArithmeticError('Inherited axial stress bound failed')
        N_source=hi(1+c.mpf('2.01')*(Kp+100)+64*eps0+c.mpf('27.6375')+c.mpf('27.5'))
        if endpoints(N_source)[1]>=500+4*Kp or 501+4*Kp>=KN:
            raise ArithmeticError('Axial stress normalization source bound failed')
        u_restore_min=c.exp(c.mpf('-.8'))/2  # u/Pstar at Rz and throughout restoration
        J_over_KN=hi(1/u_restore_min**2)
        if endpoints(J_over_KN)[1]>=20:
            raise ArithmeticError('Normalized axial J bound failed')
        vy=hi(16*eps0);b_times_P=hi(144*eps0)
        bw=hi(1280*KN*eps0)
        kappa=hi(c.mpf('.8')+c.mpf('1.25')*b_times_P**2)
        restore_cone_margin=lo(3*(c.mpf('.8')-bw)-c.mpf('1.6'))
        if (endpoints(bw)[1]>=mp.mpf('.01') or endpoints(kappa)[1]>=1
                or endpoints(restore_cone_margin)[0]<=0):
            raise ArithmeticError('Axial-restoration relaxed cone failed')

        # Callable enclosure examples use the full admitted B jet box.
        # They are log fields, not selected point solutions or moment evaluations.
        z=IntervalTaylor(c,[c.mpf([-1,1]),c.mpf(1),c.mpf(0)])
        B=IntervalTaylor(c,[c.mpf([-endpoints(2*A)[1],endpoints(2*A)[1]]),
                            c.mpf([-endpoints(2*A)[1],endpoints(2*A)[1]]),
                            c.mpf([-endpoints(A)[1],endpoints(A)[1]])])
        shapes=[dict(phase=p,**shape_log_jet(c,p,A,logC,logP,z,B)) for p in ('0','.5','1')]
        reference={p:list(reference_log_jet(c,p,z).coefficients) for p in ('-8','-7','-6','-5')}
        velocity_error=hi(get(exit,'axial_velocity_C2_error_upper'))
        v1=z*4+IntervalTaylor(c,[c.mpf([-endpoints(velocity_error)[1],endpoints(velocity_error)[1]])
                               for _ in range(3)])
        restores=[dict(relative_offset=t,**restore_axial_jet(c,t,z,v1)) for t in ('0','.5','1')]
        definition=dict(admitted_inner_parameter_family_sha256=ledger['admitted_inner_parameter_family_sha256'],
            Abar='proved upper for exact9.16 A, fixed once for selected family',T='400*Abar',
            Rref='110*(Cstar*Pstar)^10, exact expression',Rsh='110*exp(T)',
            offsets=dict(Rz=-8,restore_end=-7,Rm=-6,Rh=-5),
            swirl='log(u/Pstar)=y/10-logCstar-logPstar-log(1+Z^2)+(1-sigma(y/T))*B',
            B='log(Cstar*u1*(1+Z^2)) of the actual R110 exit',
            axial='v1 toRz; v1+(4Z-v1)*sigma(log(R/Rz)) to eRz;4Z afterwards',
            moments='continue all five exact axis primitives, P=P0+Mp; no resetting')
        sha=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        proof=dict(
            G_nonnegative='H0=(4Z+j)(1-Z^2)+(1-delta)Z/2 is negative forZ<=-j andpositive forZ>=0; on[-j,0] H0_prime>4 forj<.001,delta<=.005, so the unique root is a global minimum for G anchored atzero there',
            G_primitive_branch='G(Z)=integral_Z0^Z L(s)H0(s)/(H0(s)^2+sigma^2)ds is the implicit primitive definition; seed acceptance binds ell=-G_prime jets andG(Z0)=0. The positive denominator andL sign giveG>=0. No pole-log primitive is used for this family.',
            amplitude='Abar-10 bounds the core log(Cstar F) C3; integrated actual andswitch relativelogs plus constant-power, sqrt220 andlog(1+Z^2) C2 add<10; hence B C2<Abar<=2Abar',
            T_choice='the conservative fixed Abar replaces exactA only in T; F27 used the same Abar in its stronger radius restriction, so no accepted core/outer data is changed',
            slope='sigma_prime<=8 and|B|<=2Abar give|a-.8|<=.08, safely[.7,.9]; logu_y=(1-a)/2>=.05',
            signed_gradient='zeta convexly interpolates entryzeta and-2Z/(1+Z^2); forV=4Z+error, Hv*(-2Z/(1+Z^2)) has nonpositive leadingterm anderror<=epsilon0; same bound holds throughrestoration',
            mean='Mz/R exact positive Z-independent average of inherited mean andsubsequent velocities; both stay<epsilon0 from4Z, including restoration, so Werror<=10epsilon0',
            angular='9.32 SQ has the stored lower bound>7/5; damping[31/20,33/20], positiveD3crossing andone-unitQ>.5 prove relaxedcone before restoration',
            early_pressure='uniform normalized Phi<=2 andG>=0, relativelog<.5 implyF<=4/Cstar, |logF_Z|<=Abar+1; integrate110*16*(3+2Abar)/Cstar^2; the displayed exponentmargin bounds this byexp4Abar/Cstar^2<=1',
            pressure='shapeC1<= (5+30Abar)u_sh^2, compatibilityu_sh^2/Pstar^2<=exp(-2)/(1+Abar)^2, referenceC1<=7.5Pstar^2, sameP0C1<=KpPstar^2; no pressure-tail addition',
            inherited_N='angular/mixed/pressure increments<1; MzC1increments<=990,MzthetaC1increments<=8911; inheritedcoremoments/P0 chargedtoK, K>=1e6 => each actualmoment/P C1<=2K; formula9.13 givesN110/K<20',
            N_source='9.35 source bound: Zu^2<=Pstar^2, PopP<=2.01(Kp+100)Pstar^2, |W Vy|<=64epsilon0, otheraxialterms<=27.6375+27.5; thus source< (500+4Kp)Pstar^2',
            inherited_stress='integratingN_y+N retains110N110/R; the admitted Rz radius makes it<=Pstar^2, and501+4Kp<KN',
            restoration='u/Pstar>=exp(-.8)/2, J=N/u^2<=20KN, Q>=.5; bw=2J Vy/Q<=1280KNepsilon0, b<=144epsilon0/Pstar, kappa<1, D(a-bw)>2a',
            topology='offsets -8,-7,-6,-5 are exact from the formal logRref expression; finite differences are never computed by subtracting huge rounded absolute logs',
            limitations='exact implicit prescription and analyticbounds; logjet examples are enclosures of input families, not completed physical evaluation, five-moment repair, outerheatfield orstresslift')
        result=dict(reference_join_family_sha256=sha,definition=definition,proof=proof,
            admitted_inner_parameter_family_sha256=ledger['admitted_inner_parameter_family_sha256'],
            uniform_Cstar_family_sha256=norms['uniform_Cstar_family_sha256'],
            real_axial_domain=['-1','1'],radial_domain=['110','Rh=exp(-5)*Rref'],
            Abar=A,T=T,B_C2_upper=2*A,normalized_R110_log_extra_C2_upper=extra,
            G_real_axis_sign_gates=G_sign_gates,
            G_branch_definition='real integral ofL H0/(H0^2+sigma^2) anchored atuniqueH0root',
            pole_log_primitive_used=False,
            short_switch_relative_log_C2_upper=switch_log,
            shape_to_Rz_log_gap_lower=shape_to_z_gap,inherited_N_radius_log_margin=inherited_radius_margin,
            shape_shear_interval=[a_min,a_max],signed_Hv_zeta_upper=signed_shape,
            angular_SQ_lower=SQ,D3_crossing_derivative_lower=barrier,Q_after_one_log_unit_lower=q_one_unit,
            early_Mp_log_upper=early_mp_log_upper,early_Mp_exp_coefficient_margin=early_mp_exp_margin,
            early_theta_increment_log_upper=theta_increment_log_upper,
            early_mixed_increment_log_upper=mixed_increment_log_upper,
            normalized_shape_pressure_C1_upper=shape_pressure,normalized_join_pressure_C1_upper=pressure_C1,
            inherited_N_R110_over_K_upper=inherited_N_over_K,normalized_N_source_upper=N_source,
            pressure_Kp=Kp,KN=KN,N_after_Rz_over_Pstar_squared_upper=KN,
            J_after_Rz_upper=20*KN,axial_restoration_Vy_upper=vy,
            axial_restoration_b_times_Pstar_upper=b_times_P,axial_restoration_bw_upper=bw,
            axial_restoration_kappa_upper=kappa,axial_restoration_relaxed_lower_margin=restore_cone_margin,
            log_shape_enclosure_examples=shapes,exact_offset_reference_log_jet_examples=reference,
            axial_restore_enclosure_examples=restores,
            normalized_R110_amplitude_budget_recomputed=True,
            real_axis_G_nonnegative_normalization_verified=True,
            callable_log_shape_and_axial_restoration_installed=True,
            same_family_R110_Rh_relaxed_cone_analytically_certified=True,
            same_preheat_axis_pressure_retained=True,actual_moments_continuously_inherited=True,
            physical_radius_or_amplitude_materialized=False,whole_axis_finite_velocity_evaluator_built=False,
            five_terminal_moment_identities_repaired=False,corrected_outer_at_selected_radius_built=False,
            heat_exterior_matched=False,admissible_stress_lift_constructed=False,
            full_Section9_parameter_admission=False,temporal_recursion=False,
            input_hashes={**hashes,**{PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in names},
                **{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                   (Path(__file__).name,PREFIX+'interval_taylor.py',PREFIX+'interval_functional_defects.py')}})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Same-family log reshape/reference/axial restore: analytic relaxed cone throughRh PASS',flush=True)
        print('SQ lower=',mp.nstr(endpoints(SQ)[0],14),'restoration bw<=',mp.nstr(endpoints(bw)[1],14),flush=True)
        return result


if __name__=='__main__':
    run()
