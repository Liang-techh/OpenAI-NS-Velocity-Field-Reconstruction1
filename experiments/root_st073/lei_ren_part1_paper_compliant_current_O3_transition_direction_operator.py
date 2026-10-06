"""Actual variable O3 transition direction and strict-shear endpoint obstruction."""
import importlib
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_power_cone_operator import (
    current_O3_power_source_theorem,SourceAST,pack,encode,endpoints,source_precision,upper,absolute)

DOMAIN=(0,1)


def actual_variable_transition_theorem():
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('Actual O3 transition identity failed: '+name)
        checks[name]=True
    z,mu,delta,v,J,Q=s.symbols('Z mu delta offset same_J same_Q',real=True)
    Ua,Ma,da,EZa,EQa,KE=s.symbols('actual_Ua actual_Ma actual_da actual_EZa actual_EQa actual_KE',real=True)
    C=1/(1+z*z);L=1-delta*z*z;k=1-delta/2
    A=k*C+(1-delta)*z*z*C*C
    N=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp)
    KT=s.exp(v-mu*J)-1+mu*Q
    native=dict(Utheta_over_Pstar=Ua*C,Mz_over_R=Ma*z,
        Mtheta_over_sqrt2_R_3half_Pstar=Ua*(1-da)*C,
        Mtheta_z_over_sqrt2_R_3half_Pstar=Ua*(Ma-4*da)*z*C,
        Mztheta_over_R_Pstar_squared=EZa*z*z+EQa*C*C)
    env=dict(c=c,t=v,f=s.exp(-v/2-mu*J),decay=s.exp(-v),d3=s.exp(-3*v/2),u1=Ua*C,
        get=lambda key:native[key],kernels=dict(theta=KT,energy=KE))
    vals={key:asts.evaluate(asts.expression('outer_buffer','slope_mu',key),env) for key in ('u','m','h','k','e')}
    U=vals['u']/C;D=da*s.exp(-v+mu*J);V=Q*s.exp(-v+mu*J);M=Ma*s.exp(-v)
    zero('actual_variable_X_deficit_and_positive_Q',vals['h']/vals['u'],1-D+mu*V)
    zero('actual_variable_shared_M_K_X_correlation',vals['k']/(z*C*U),M*s.exp(mu*J)-4*D)
    AZ=EZa/Ua**2*s.exp(2*mu*J)
    AQ=(EQa/Ua**2-KE/2)*s.exp(2*mu*J)
    CE=z*z*AZ+C*C*AQ
    zero('actual_variable_full_energy_canonical_shape',vals['e']/U**2,CE)
    sig=s.Function('same_sigma')(v);jf=s.Function('same_J')(v);uf=Ua*s.exp(-v/2-mu*jf)
    # The jet helper and the positive-integral cutoff are two evaluations
    # of the same defining function, including the reflected right half.
    asts.expression('flat_pulse_derivatives','_sigma_left','odds',wanted='1/(1-x)**2-1/x**2')
    jet_sigma=asts.expression('flat_pulse_derivatives','_sigma_left','value',wanted='e/(1+e)')
    integral_sigma=asts.expression('outer_initial','stable_sigma','value',wanted='a/(a+b)')
    aa,bb=s.symbols('same_positive_left_exp same_positive_right_exp',positive=True)
    cutoff=asts.evaluate(integral_sigma,dict(a=aa,b=bb))
    zero('actual_jet_and_integral_sigma_same_positive_fraction',
        asts.evaluate(jet_sigma,dict(e=aa/bb)),cutoff)
    zero('actual_jet_sigma_right_reflection_same_fraction',
        1-cutoff.xreplace({aa:bb,bb:aa}),cutoff)
    zero('actual_variable_logU_derivative',s.diff(uf,v).subs(s.diff(jf,v),sig),-(s.Rational(1,2)+mu*sig)*uf)
    zero('actual_variable_shear_a_minus2',1-2*(-s.Rational(1,2)-mu*sig)-2,2*mu*sig)
    # Replay the full actual raw variable-amplitude stress. Ordinary
    # derivative rows are symbolic; no fixed-power derivative is inserted.
    shifted=importlib.import_module('lei_ren_part1_paper_compliant_collar_stress_C3').shifted_rows
    product=importlib.import_module('lei_ren_part1_paper_compliant_collar_Gamma_C4').product_rows
    op=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(
        axial_derivative=lambda value:s.diff(value,z),product_rows=product,shifted_rows=shifted))
    H=U*(1-D+mu*V);K=U*(M*s.exp(mu*J)-4*D)
    # Pabs is the same analytic absolute pressure in raw Pstar^2 units.
    Hfuture=s.symbols('same_positive_future_swirl_integral',real=True)
    Pv=s.Function('same_actual_signed_Rv')(z);p=1+2*mu
    Tw=s.symbols('same_Tw',positive=True)
    U1=Ua*s.exp(-s.Rational(1,2)-mu/2)
    Pmemory=U1**2*(Pv+C*C/(2*p))*s.exp(-p*(13/mu+Tw))
    Pabs=-U**2*C*C*Hfuture/2+Pmemory
    rows=lambda expr:[expr]+[s.Function('ordinary_'+str(j))(z) for j in range(1,5)]
    hist=dict(m=[M*z*(-1)**j for j in range(5)],h=rows(H*C),k=rows(K*z*C),e=rows(U**2*CE))
    urows=[U*C,-(s.Rational(1,2)+mu*sig)*U*C]+[s.Function('variable_u'+str(j))(z) for j in range(2,5)]
    native_rows=op(c,delta,z,urows,[s.Integer(0)]*5,hist,rows(Pabs))
    inertial=sum(part['shape'][0] for name,part in native_rows['theta'].items() if name!='variable_radial_shear')/U
    theta=(A*(1-D+mu*V)-C)/L+C*M+N*(M*s.exp(mu*J)-4*D)
    zero('actual_full_variable_theta_operator',inertial,theta)
    zero('actual_variable_theta_correlated_positive_drift_split',theta,
        (A-C)/L+(-A/L-4*N)*D+mu*A*V/L+(C+N)*M+N*M*(s.exp(mu*J)-1))
    zero('actual_variable_radial_shear_operator',native_rows['theta']['variable_radial_shear']['shape'][0]/U,
        -(2+2*mu*sig)*C)
    IE=(2*delta*z*CE-(1-z*z)*s.diff(CE,z))/L
    Pbase=-C*C*Hfuture/2
    IP=(2*(1+delta)*z*Pbase-(1-z*z)*s.diff(Pbase,z))/L
    zero('actual_full_energy_future_pressure_odd_Z_factor',IE+IP,
        z*((2*delta*z*z-2*(1-z*z))*AZ+(2*delta*C*C+4*(1-z*z)*C**3)*AQ
            -Hfuture*((1+delta)*C*C+2*(1-z*z)*C**3))/L)
    full_axial=(native_rows['axial']['retained_full_energy']['shape'][0]
        +native_rows['axial']['actual_absolute_pressure']['shape'][0])/U**2
    memop=(2*(1+delta)*z*Pmemory-(1-z*z)*s.diff(Pmemory,z))/L
    zero('actual_full_energy_absolute_functional_pressure_operator',full_axial,IE+IP+memop/U**2)
    for name in ('local_axial_transport','nonlinear_meridional_transport','retained_linear_axial_moment','axial_radial_shear'):
        for j,expr in enumerate(native_rows['axial'][name]['shape']):zero('actual_'+name+'_ordinary_row'+str(j)+'_zero',expr,s.Integer(0))
    zero('actual_variable_time_AMGM_factor',U**2/D,Ua**2*s.exp(-3*mu*J)/da)
    # Original absolute pressure FTC from the same phase1 datum.
    # f(s)^2/f(v)^2 lies between exp(-p*(s-v)) and exp(-(s-v)).
    # The future integral H is therefore 1/p <= H <= 1.
    jj=s.Function('same_J')(v);integral=s.Function('same_backward_swirl_integral')(v)
    Pm1=s.symbols('same_absolute_Rw_pressure',real=True)
    f2=Ua**2*s.exp(-v-2*mu*jj)
    backward=Pm1-C*C*Ua**2*integral/2
    zero('actual_absolute_transition_pressure_FTC',
        s.diff(backward,v).subs(s.diff(integral,v),-s.exp(-v-2*mu*jj)),C*C*f2/2)
    phase=s.Symbol('periodic_phase',real=True);E0=s.Symbol('same_current_E0',positive=True)
    a0=2+2*mu*sig;aL=a0+mu*s.cos(4*s.pi*phase);bL=2*s.sqrt(mu)*s.sin(2*s.pi*phase)
    Ap=-mu*s.sin(4*s.pi*phase)/(8*s.pi)
    Bp=-E0*s.sqrt(mu)*s.cos(2*s.pi*phase)/(2*s.pi)
    zero('actual_shear_loop_a_mean',s.integrate(aL,(phase,0,1)),a0)
    zero('actual_shear_loop_b_weighted_mean',s.integrate(E0*bL,(phase,0,1)),s.Integer(0))
    zero('actual_shear_loop_A_primitive',s.diff(Ap,phase),-(aL-a0)/2)
    zero('actual_shear_loop_B_primitive',s.diff(Bp,phase),E0*bL/2)
    zero('actual_shear_loop_A_zero_mean',s.integrate(Ap,(phase,0,1)),s.Integer(0))
    zero('actual_shear_loop_B_zero_mean',s.integrate(Bp,(phase,0,1)),s.Integer(0))
    y=s.Symbol('same_sin_squared',real=True)
    zero('actual_shear_loop_vs_positive_lower_expression',
        mu*(1-2*y)+4*mu*y/(2+3*mu),mu-6*mu**2*y/(2+3*mu))
    oldT,oldZ,daL=s.symbols('actual_Ttheta_over_F actual_Tz_over_F shear_a_increment',real=True)
    zero('actual_shear_loop_uses_same_inviscid_vector_not_unchanged_stress',
        (oldZ+bL)/(oldT-daL),((oldZ/oldT)+bL/oldT)/(1-daL/oldT))
    aa,bb,ww=s.symbols('actual_loop_a actual_loop_b actual_new_stress_ratio',real=True)
    dd=1-bb*ww/aa;cross=ww+bb/aa;km=aa-2+bb*bb/aa
    zero('actual_loop_direction_bracket_equals_original_two_vector_quadratic',
        2*dd*dd-km*cross*cross,
        (1+bb*bb/(aa*aa))*(2-((aa-2)*ww*ww+2*bb*ww+bb*bb/aa)))
    return dict(identities=checks,
        formulas=dict(D='da*exp(-offset+mu*J)',X='1-D+mu*Q*exp(-offset+mu*J)',
            K_over_U='M*exp(mu*J)-4D',pressure_future_H='[f(1)^2/p+integral_offset^1 f(s)^2 ds]/f(offset)^2',
            pressure_memory='U1^2*(Pv+C^2/(2p))*exp(-p*(13/mu+Tw))',
            shear_a_minus2='2*mu*sigma(offset)'),
        continuous_inequalities=dict(nonnegative_original_integrals='0<=J<=1/2, 0<=KE<=offset<=1, Q>=0',
            theta='Theta>=cX*Z^2+cD*D(offset), D>=da*exp(-1)>0',
            pressure='1/p<=H<=1 by 0<=sigma<=1 on the same original future integral',
            drift='|N*M*(exp(mu*J)-1)|<=Ncap*Ma*mu*exp(mu/2)/2',
            energy_direction='2*mu*w_energy^2<=mu*Pstar^2*Ua^2*Lz^2/(2*cD*cX*da)'),
        strict_shear_positive_only_on='offset in (0,1]',
        explicit_mean_preserving_shear_loop=dict(a='a0+mu*cos(4*pi*phase)',b='2*sqrt(mu)*sin(2*pi*phase)',
            A='-mu*sin(4*pi*phase)/(8*pi)',B='-E0*sqrt(mu)*cos(2*pi*phase)/(2*pi)',
            same_inviscid_vector_target='Ttheta_loop/F=Ttheta0/F-(aL-a0); Tz_loop/F=Tz0/F+bL',
            spatial_support_and_finite_N_and_five_moment_repair_not_installed=True),
        exact_offset0_shear_excess_zero=True,offset0_nonzero_stress_requires_original_shear_repair=True,
        full_closed_transition_admissible_cone=False,input_hashes=asts.hashes,passed=True)


def current_variable_transition_bindings(powercone):
    powercone.assert_graph();owner=powercone.registry.owners['o3']
    if not powercone.acceptance_loaded or not owner.acceptance_loaded:raise ValueError('Checked current O3 power/transition required')
    source=current_O3_power_source_theorem(powercone.entrancecone);asts=SourceAST()
    rawmod=importlib.import_module('lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator')
    if owner.chart.__func__.__wrapped__.__globals__['raw_pre_stress_rows'] is not rawmod.raw_pre_stress_rows:
        raise ValueError('Original variable raw stress operator required')
    sigmod=importlib.import_module('lei_ren_part1_paper_compliant_flat_pulse_derivatives')
    if owner.physical.pre.slope_mu.__func__.__globals__['sigma_jets'] is not sigmod.sigma_jets:
        raise ValueError('Original variable cutoff jet function required')
    asts.expression('pre_pulse_mixed_C4','slope_mu','K',wanted='transition_kernels(c,t,mu,self.cells)')
    asts.expression('pre_pulse_mixed_C4','slope_mu','logU',
        wanted="[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]")
    asts.expression('current_O3_transition_background_tensor','chart','P',
        wanted="[copy_jet(c,raw['p'][0])+datum]+[copy_jet(c,row) for row in raw['p'][1:]]")
    return dict(current_actual_power_composed_function_and_pressure_proof=source,
        current_variable_transition_complete_operator_and_join=owner.proof,
        same_current_original_variable_kernels_and_all_logU_derivatives=True,
        same_exact_positive_cutoff_for_jets_and_history_integrals=True,
        same_signed_absolute_Rv_pressure_and_power_phase0_function=True,
        no_fixed_power_rate_substituted_for_variable_transition=True,
        input_hashes={**source['input_hashes'],**asts.hashes},passed=True)


def validate_whole_transition(powercone,view):
    if tuple(view[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) \
            !=(powercone.family,powercone.source,powercone.datum_sha):raise ValueError('Foreign transition source')
    raw=view['original_complete_view']
    if raw['chart']!='O3_slope_mu' or endpoints(raw['Z'])!=(-1,1) or endpoints(raw['coverage_coordinate'])!=(0,1):
        raise ValueError('Whole variable transition source required')
    for flag in ('actual_full_stress_not_local_difference','all_variable_log_amplitude_ordinary_derivatives_retained',
        'full_nonzero_histories_and_radial_remainder_retained_when_local_V_zero'):
        if not raw[flag]:raise ValueError('Complete original variable source required: '+flag)
    for key in ('actual_upstream_original_pre_transition_source','current_raw_five_history_rows',
        'current_absolute_pressure_ordinary_y_rows','current_source_three_component_velocity_rows'):
        if key not in raw:raise ValueError('Current full transition source omitted: '+key)
    expected={'theta':{'local_transport':(.5,1,0,0),'retained_angular_moment':(.5,1,0,0),
        'retained_mixed_moment':(.5,1,0,0),'meridional_transport':(.5,1,0,0),'variable_radial_shear':(-.5,1,0,0)},
        'axial':{'local_axial_transport':(.5,0,0,0),'nonlinear_meridional_transport':(.5,0,0,0),
            'retained_linear_axial_moment':(.5,0,0,0),'retained_full_energy':(.5,2,0,0),
            'actual_absolute_pressure':(.5,2,0,0),'axial_radial_shear':(-.5,0,0,0)}}
    sectors=raw['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if set(sectors)!=set(expected):raise ValueError('Original transition stress layout required')
    for label,parts in expected.items():
        if set(sectors[label])!=set(parts):raise ValueError('Full variable signed sectors required')
        for name,mode in parts.items():
            part=sectors[label][name];row=raw['physical_cylindrical_stress_mixed3'][label][name]['r0_z0']
            if tuple(part['mode'])!=mode or encode(pack(row['signed_coefficient']))!=encode(pack(part['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])) \
                    or row['actual_source_log_parts']!=part['exact_source_log_parts']:
                raise ValueError('Actual transition signed coefficient/mode/log changed')
    return raw


@source_precision
def whole_current_transition_bounds(powercone,view,inlet):
    raw=validate_whole_transition(powercone,view);c=powercone.ctx
    owner=powercone.registry.owners['o3'];mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta);p=1+2*mu
    positive={}
    def require(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Actual transition bound unresolved: '+name)
        positive[name]=value
    for name,value in dict(mu=mu,delta=delta,one_minus_delta=1-delta,delta_below_one_eighth=c.mpf('.125')-delta).items():require(name,value)
    if inlet['chart']!='O2_11_unit_buffer' or endpoints(inlet['Z'])!=(0,0) or endpoints(inlet['coverage_coordinate'])!=(11,11):
        raise ValueError('Actual canonical O2 buffer endpoint source required')
    iv=lambda value:c.mpf(endpoints(value));history=inlet['actual_normalized_primitive_y_derivative_axial5']
    Ua=iv(inlet['Utheta_over_Pstar_axial5_coefficients'][0]);Ma=iv(history['m'][0][1]);require('actual_Ua',Ua);require('actual_Ma',Ma)
    source=powercone.bounds;ds=source['actual_O2_slope_deficit'];tau=c.mpf(owner.physical.pre.params.logPstar)-1
    logda=c.ln(ds)-tau;da=c.exp(logda);dmin=da*c.exp(-1)
    cX=(1-delta)/4;cD=(c.mpf(15)/8-14*delta)/8
    delta_cap=delta/(2*(1-delta));drift=2*Ma*mu*c.exp(mu/2)/(2*(1-delta))
    require('positive_deficit_coefficient',cD);require('deficit_reserve_dominates_delta_drift',cD*dmin-delta_cap-drift)
    logr=iv(raw['current_actual_source_stress_packet']['exact_logR'])
    radial_relative_log=c.ln(2+2*mu)-logr-c.ln(cD*dmin)
    require('radial_relative_log_gap',-1000-radial_relative_log)
    radial_cap=cD*dmin*c.exp(-1000)
    reserve=cD*dmin-delta_cap-drift-radial_cap;require('full_theta_positive_reserve',reserve)
    theta_floor=cD*dmin+reserve;require('whole_theta_lower',theta_floor)
    EQa=iv(history['e'][0][0]);EZa=iv(history['e'][0][2])+2*EQa;growth=c.exp(mu)
    AQ=upper(c,(absolute(c,EQa/Ua**2)+c.mpf('.5'))*growth);AZ=upper(c,absolute(c,EZa/Ua**2)*growth)
    Lz=upper(c,((2+2*delta)*AZ+(4+2*delta)*AQ+3+delta)/(1-delta));require('full_energy_future_pressure_Z_factor',Lz)
    lp=c.mpf(owner.physical.logP);lu=c.mpf(inlet['log_Utheta_over_Pstar_base_source'])
    energy_log=c.ln(mu)+2*lp+2*lu+2*c.ln(Lz)-c.ln(2*cD*cX)-logda
    require('whole_energy_direction_log_gap',-1000-energy_log)
    Pv=powercone.entrancecone.whole_view['original_complete_view']['current_signed_absolute_Rv_pressure']
    pm=upper(c,(2*(1+delta)*(absolute(c,Pv[0])+1/(2*p))+absolute(c,Pv[1])+2/p)/(1-delta))
    Tw=c.mpf(owner.physical.pre.params.Tw);require('actual_Tw',Tw)
    pressure_log=lp+lu+c.ln(pm)-p*(13/mu+Tw)-c.ln(theta_floor)
    require('whole_functional_pressure_memory_log_gap',-1000-pressure_log)
    ecap=c.exp(-1000);wE=c.sqrt(ecap/(2*mu));wP=c.exp(-1000);direction=2*mu*(wE+wP)**2
    require('direction_below_one_point8',c.mpf('1.8')-direction);require('full_directional_bracket',2-direction)
    # A constructive periodic shear loop for the SAME original inviscid
    # vector ps. Stress changes with shear, so w is not held fixed.
    require('loop_mu_below_one_sixth',c.mpf(1)/6-mu)
    minimum_stress_over_F_log=logr+c.ln(theta_floor)
    require('original_stress_over_F_inverse_log_gap',minimum_stress_over_F_log-1000)
    inverse_t0=c.exp(-1000);w0=c.sqrt(direction/(2*mu))
    require('loop_stress_denominator_positive',c.mpf('.5')-mu*inverse_t0)
    loop_w=2*(w0+2*c.sqrt(mu)*inverse_t0)
    loop_vs_lower=mu-3*mu**2;require('loop_vs_exceeds_half_mu',loop_vs_lower-mu/2)
    loop_align=2-mu-2*c.sqrt(mu)*loop_w;require('loop_signed_alignment',loop_align)
    loop_direction=3*mu*loop_w**2+4*c.sqrt(mu)*loop_w+4*mu/(2-mu)
    require('loop_direction_below_one_point8',c.mpf('1.8')-loop_direction)
    return dict(positive_margins=positive,mu=mu,delta=delta,actual_Ua=Ua,actual_Ma=Ma,actual_da_log=logda,
        actual_da=da,minimum_actual_D=dmin,cX=cX,cD=cD,full_theta_uniform_lower=theta_floor,
        full_theta_lower_formula='cX*Z^2+cD*D(offset)+positive_reserve',positive_theta_reserve=reserve,
        full_energy_EQa=EQa,full_energy_EZa=EZa,full_energy_AQ_upper=AQ,full_energy_AZ_upper=AZ,
        original_future_pressure_H_bounds=dict(lower=1/p,upper=c.mpf(1)),
        full_energy_and_future_pressure_Z_factor_upper=Lz,
        full_energy_direction_log_upper=energy_log,full_pressure_memory_log_upper=pressure_log,
        original_Pv_and_Z_derivatives=Pv,functional_pressure_memory_coefficient_upper=pm,
        variable_radial_relative_log=radial_relative_log,variable_radial_error_upper=radial_cap,
        full_directional_expression_upper=direction,whole_closed_direction_and_positive_stress_certified=True,
        source_correlated_periodic_shear_loop=dict(original_stress_over_F_log_lower=minimum_stress_over_F_log,
            inverse_original_stress_over_F_upper=inverse_t0,new_stress_ratio_absolute_upper=loop_w,
            shear_vs_minus2_lower=loop_vs_lower,alignment_lower=loop_align,directional_expression_upper=loop_direction,
            same_inviscid_vector_with_changed_loop_target=True,strict_loop_for_all_phase_and_current_Z_offset=True,
            spatial_support_and_finite_frequency_and_moment_repair_installed=False),
        shear_excess_formula='vs-2=2*mu*sigma(offset)',exact_offset0_shear_excess=0,
        offset0_stress_nonzero_by_same_positive_theta_lower=True,
        strict_shear_on_open_offset_domain_by_original_cutoff_positivity=True,
        closed_transition_strict_cone_certified=False,original_shear_repair_required=True,
        complete_original_signed_stress_pressure_velocity_and_remainder_preserved=True,
        defining_source_not_caps_or_phase_samples=True)
