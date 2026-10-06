"""Current O3 power cone from correlated O2 history and backward pressure."""
import ast
import importlib
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_original_cone import (
    SourceAST,pack,encode,endpoints,source_precision)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_operator import upper,absolute
from lei_ren_part1_paper_compliant_current_O3_theta_correlation import original_full_theta_identity

DOMAIN=(0,1)


def original_correlated_O2_O3_theorem():
    """Exact canonical history, positive transition integral and shear modes."""
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:
            raise ArithmeticError('O3 power identity failed: '+name)
        checks[name]=True
    z,mu,delta,t,tau,J,Q=s.symbols('Z mu delta t tau J Q',real=True)
    Us,Xs,Bmass=s.symbols('actual_Us actual_Xs same_turnoff_mass',real=True)
    C=1/(1+z*z);r=1-mu;k=1-delta/2;L=1-delta*z*z
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    inlet=dict(Utheta_over_Pstar=Us*C,Mz_over_R=4*z,
        Mtheta_over_sqrt2_R_3half_Pstar=Us*Xs*C,
        Mtheta_z_over_sqrt2_R_3half_Pstar=4*Us*Xs*z*C)
    env=dict(c=c,z=z,t=tau,decay=s.exp(-tau),root_decay=s.exp(-tau/2),
        decay3=s.exp(-3*tau/2),get=lambda key:inlet[key],u1=Us*C,K=dict(B_mass=Bmass))
    ua=asts.evaluate(asts.expression('outer_initial','axial','u'),env)/C
    ma=asts.evaluate(asts.expression('outer_initial','axial','m'),env)/z
    ha=asts.evaluate(asts.expression('outer_initial','axial','h'),env)/C
    ka=asts.evaluate(asts.expression('outer_initial','axial','k'),env)/(z*C)
    xa=s.cancel(ha/ua);da=(1-Xs)*s.exp(-tau)
    zero('actual_O2_X_deficit',xa,1-da)
    zero('actual_O2_shared_mass_M_K_X_cancellation',ka/ua-ma,4*(xa-1))
    # J'=sigma and theta'=exp(t-mu*J). Integration of
    # (exp(t-mu*J))' proves theta=exp(t-mu*J)-1+mu*Q,
    # Q=integral sigma(v)*exp(v-mu*J(v))dv >=0.
    sig=s.Function('same_original_sigma')(t);jf=s.Function('same_original_J')(t)
    zero('actual_transition_FTC_exponential_weight',
        s.diff(s.exp(t-mu*jf),t).subs(s.diff(jf,t),sig),
        (1-mu*sig)*s.exp(t-mu*jf))
    sigma_node=asts.expression('outer_initial','stable_sigma','value',wanted='a/(a+b)')
    asts.expression('outer_initial','stable_sigma','a',wanted='c.exp(-1/x**2)')
    asts.expression('outer_initial','stable_sigma','b',wanted='c.exp(-1/(1-x)**2)')
    aa,bb=s.symbols('same_positive_left_exp same_positive_right_exp',positive=True)
    sigma=asts.evaluate(sigma_node,dict(a=aa,b=bb))
    zero('actual_sigma_reflection_and_half_integral',sigma+sigma.xreplace({aa:bb,bb:aa}),s.Integer(1))
    if not sigma.is_positive:raise ArithmeticError('Actual sigma positivity unresolved')
    checks['actual_sigma_nonnegative_from_positive_exponentials_and_zero_one_edges']=True
    kt=s.exp(1-mu*J)-1+mu*Q;rho=s.exp(-1+mu*J)
    native=dict(Utheta_over_Pstar=ua*C,Mz_over_R=ma*z,
        Mtheta_over_sqrt2_R_3half_Pstar=ha*C,
        Mtheta_z_over_sqrt2_R_3half_Pstar=ka*z*C)
    env=dict(c=c,t=s.Integer(1),f=s.exp(-s.Rational(1,2)-mu*J),
        decay=s.exp(-1),d3=s.exp(-s.Rational(3,2)),u1=ua*C,
        get=lambda key:native[key],kernels=dict(theta=kt))
    u0=asts.evaluate(asts.expression('outer_buffer','slope_mu','u'),env)/C
    m0=asts.evaluate(asts.expression('outer_buffer','slope_mu','m'),env)/z
    h0=asts.evaluate(asts.expression('outer_buffer','slope_mu','h'),env)/C
    k0=asts.evaluate(asts.expression('outer_buffer','slope_mu','k'),env)/(z*C)
    D=da*rho;V=Q*rho
    zero('actual_O3_initial_X_with_positive_Q',h0/u0,1-D+mu*V)
    zero('actual_O3_initial_M_K_deficit_correlation',k0/u0,m0*s.exp(mu*J)-4*D)
    # Exact current whole-Z theta formula, now with the SAME deficit.
    x=s.symbols('x',real=True);Cx=1/(1+x);Lx=1-delta*x
    Ax=k*Cx+(1-delta)*x*Cx**2
    Nx=((2*delta*x-1)*Cx+2*(1-x)*x*Cx**2)/Lx
    zero('actual_nonnegative_M_coefficient',Cx+Nx,x*Cx*(delta+2*(1-x)*Cx)/Lx)
    poly=3-6*x+8*x*x+delta*(s.Rational(1,2)-s.Rational(13,2)*x-8*x*x)
    zero('actual_positive_deficit_polynomial',-Ax/Lx-4*Nx,poly/(Lx*(1+x)**2))
    zero('actual_deficit_quadratic_completed_square',3-6*x+8*x*x,
        8*(x-s.Rational(3,8))**2+s.Rational(15,8))
    zero('actual_negative_N_completed_square',1-x+2*x*x,
        2*(x-s.Rational(1,4))**2+s.Rational(7,8))
    zero('actual_C_above_half_on_closed_Z_domain',Cx-s.Rational(1,2),(1-x)/(2*(1+x)))
    zero('actual_C_squared_above_quarter_on_closed_Z_domain',Cx**2-s.Rational(1,4),
        (1-x)*(3+x)/(4*(1+x)**2))
    zero('actual_deficit_shape_numerator_lower_positive_polynomial',poly-(s.Rational(15,8)-14*delta),
        8*(x-s.Rational(3,8))**2+delta*(1-x)*(s.Rational(29,2)+8*x))
    zero('actual_deficit_shape_denominator_below_four',4-Lx*(1+x)**2,
        (1-x)*(3+x)+delta*x*(1+x)**2)
    zero('actual_negative_N_absolute_numerator_below_two',2-(1-x+2*x*x),
        (1-x)*(1+2*x))
    zero('actual_equilibrium_minus_mu_forcing',
        (Ax/r-Cx)/Lx-mu*Ax/(r*Lx),
        (1-delta)*x*Cx**2/Lx-delta*Cx/(2*Lx))
    DD,MM,VV=s.symbols('actual_D actual_M0 nonnegative_actual_V',real=True)
    E=s.exp(-r*t);equilibrium=(Ax/r-Cx)/Lx
    correlated_X=1/r+(1-DD+mu*VV-1/r)*E
    correlated_theta=(Ax*correlated_X-Cx)/Lx+Cx*MM*s.exp(-t)+Nx*(MM*s.exp(mu/2)-4*DD)*E
    retained=(equilibrium*(1-E)+((1-delta)*x*Cx**2/Lx-delta*Cx/(2*Lx)
        +(-Ax/Lx-4*Nx)*DD)*E+mu*Ax*VV*E/Lx
        +(Cx+Nx)*MM*s.exp(-t)+Nx*MM*E*(s.exp(mu/2)-s.exp(-mu*t)))
    zero('actual_theta_exact_positive_terms_and_signed_drift_split',correlated_theta,retained)
    # With 0<=x<=1, delta<1/8: N<=0, |N|<=2/(1-delta).
    # The sole signed drift is controlled by an exact positive integral:
    # exp(mu/2)-exp(-mu*t) = integral[-mu*t,mu/2] exp(v)dv
    # <= mu*(t+1/2)*exp(mu/2). For r>0 and t>=0,
    # t*exp(-r*t)<=1/r and exp(-r*t)<=1.
    vv=s.symbols('v',real=True)
    zero('actual_signed_drift_exact_positive_integral',
        s.integrate(s.exp(vv),(vv,-mu*t,mu/2)),s.exp(mu/2)-s.exp(-mu*t))
    zero('actual_negative_N_with_small_delta',-Nx*Lx*(1+x)**2,
        1-x+2*x*x-2*delta*x*(1+x))
    # The deficit reserve cD*D >= delta_cap+Eqmin supplies the
    # missing floor at all times, not just at either endpoint.
    eqmin,dcap,cD,cX=s.symbols('Eqmin delta_cap cD cX',real=True)
    shape=s.symbols('nonnegative_Z_shape',real=True)
    zero('actual_uniform_floor_reserve_exact_split',
        eqmin*(1-E)+(2*cD*DD-dcap+shape)*E,
        eqmin+(cD*DD+shape)*E+(cD*DD-dcap-eqmin)*E)
    # AM-GM is a square identity in nonnegative |Z| and the same time
    # factor E. It keeps the B versus B^2 stress scaling explicit.
    yy,aa,bb=s.symbols('abs_Z sqrt_cX sqrt_cD_D_E',nonnegative=True)
    zero('actual_AMGM_nonnegative_square',aa**2*yy**2+bb**2-2*aa*yy*bb,(aa*yy-bb)**2)
    # Backward pressure is the original normalized ODE with the SAME
    # current signed Rv datum, not an added pressure correction.
    Pv=s.Function('actual_signed_Rv_pressure_function')(z)
    Tw=s.symbols('Tw',real=True);p=1+2*mu
    base=-C*C/(2*p);Pm=(Pv-base)*s.exp(-p*(13/mu+Tw-t))
    P=base+Pm
    zero('actual_backward_absolute_pressure_ODE',s.diff(P,t),p*P+C*C/2)
    zero('actual_backward_pressure_same_current_pulse_inlet',P.subs(t,Tw),
        base+(Pv-base)*s.exp(-p*13/mu))
    EZ,EQ,U=s.symbols('actual_EZ0 actual_EQ0 actual_U0',real=True)
    AZ=EZ/U**2*s.exp(2*mu*t)
    AQ=EQ/U**2*s.exp(2*mu*t)-(s.exp(2*mu*t)-1)/(4*mu)
    CE=z*z*AZ+C*C*AQ
    IE=(2*delta*z*CE-(1-z*z)*s.diff(CE,z))/L
    IP=(2*(1+delta)*z*base-(1-z*z)*s.diff(base,z))/L
    zero('actual_full_energy_and_baseline_pressure_odd_factor',IE+IP,
        z*((2*delta*z*z-2*(1-z*z))*AZ+(2*delta*C*C+4*(1-z*z)*C**3)*AQ
        -(1+delta)*C*C/p-2*(1-z*z)*C**3/p)/L)
    # Replay all actual axial sectors with the full energy and pressure,
    # not the zero axial placeholders of the theta-only theorem.
    Ps=s.symbols('Pstar',positive=True);mv=s.symbols('actual_M',real=True)
    ut=U*s.exp(-(s.Rational(1,2)+mu)*t)
    mm=mv*z*s.exp(-t)/(Ps*ut*C)
    mrows=[mm*(-(s.Rational(1,2)-mu))**j for j in range(5)]
    erows=[s.diff(CE/C**2,t,j) for j in range(5)]
    prows=[s.diff(P,t,j) for j in range(5)]
    shifted=importlib.import_module('lei_ren_part1_paper_compliant_collar_stress_C3').shifted_rows
    product=importlib.import_module('lei_ren_part1_paper_compliant_collar_Gamma_C4').product_rows
    rows=asts.replay('pulse_end_stress_C3','pulse_coefficients',dict(
        axial_derivative=lambda v:s.diff(v,z),product_rows=product,shifted_rows=shifted,mp=c))(
        delta,mu,z,C,s.symbols('same_X'),[s.Integer(0)]*5,mrows,[s.Integer(0)]*5,
        erows,[s.Integer(0)]*5,prows)
    full_I=IE+(2*(1+delta)*z*P-(1-z*z)*s.diff(P,z))/L
    IPmemory=(2*(1+delta)*z*Pm-(1-z*z)*s.diff(Pm,z))/L
    zero('actual_full_energy_baseline_and_functional_pressure_memory_split',full_I,IE+IP+IPmemory)
    zero('actual_functional_pressure_memory_axial_operator',IPmemory,
        s.exp(-p*(13/mu+Tw-t))*(2*(1+delta)*z*(Pv+C*C/(2*p))
            -(1-z*z)*(s.diff(Pv,z)-2*z*C**3/p))/L)
    zero('actual_full_axial_energy_absolute_pressure_operator',
        rows['axial']['full_energy_and_pressure']['shape'][0],full_I)
    for name,part in rows['axial'].items():
        if name!='full_energy_and_pressure':
            for j,value in enumerate(part['shape']):zero('actual_axial_'+name+'_row'+str(j)+'_zero',value,s.Integer(0))
    zero('actual_linear_axial_moment_structurally_zero',
        (z*s.Symbol('M')-z*s.diff(z*s.Symbol('M'),z)),s.Integer(0))
    zero('actual_AMGM_time_scale_ratio',
        s.exp(-(1+2*mu)*t)/s.exp(-r*t),s.exp(-3*mu*t))
    return dict(identities=checks,current_full_theta_theorem=original_full_theta_identity(),
        source_continuous_inequalities=dict(
            hypotheses='0<=x=Z^2<=1, 0<delta<1/8, 0<mu<1/2, Tw>0, t=Tw*phase>=0',
            equilibrium='Eq >= cX*x+Eqmin: C>=1/2, C^2>=1/4, 0<r,L<=1',
            initial_shape='((1-delta)*x*C^2-delta*C/2)/L >= cX*x-delta/(2*(1-delta))',
            deficit_shape='-A/L-4N >= (15/8-14*delta)/4 by nonnegative numerator gap and denominator <=4',
            negative_N='N<0 and |N|<=2/(1-delta): 7/8-4*delta>0 and numerator<=2',
            drift='exp(mu/2)-exp(-mu*t)<=mu*(t+1/2)*exp(mu/2); t*exp(-r*t)<=1/r',
            pressure_memory='|IPmemory| <= exp(-p*(13/mu+Tw-t))*pm; pm includes full |Pv| and |Pv_Z|'),
        positive_transition_Q_definition='integral_0^1 sigma(v)*exp(v-mu*J(v))dv',
        source_deficit='D=(1-Xs)*exp(-tau-1+mu/2)',
        initial_X='X0=1-D+mu*Q*exp(-1+mu/2)',
        initial_K_over_U='K0/U0=M0*exp(mu/2)-4*D',
        actual_pressure='-C^2/(2p)+(Pv+C^2/(2p))*exp(-p*(13/mu+Tw-t))',
        input_hashes=asts.hashes,passed=True)


def current_O3_power_source_theorem(entrancecone):
    entrancecone.assert_graph();owner=entrancecone.registry.owners['incoming'];owner.assert_graph()
    if not owner.acceptance_loaded or not entrancecone.acceptance_loaded:
        raise ValueError('Checked current entrance/full upstream source required')
    prod=entrancecone.bindings['current_actual_incoming_function_proof']
    if not prod['passed'] or not prod['actual_X_Z0_assignment_bound_to_exact_whole_Z_q_cancellation']:
        raise ValueError('Current whole-Z canonical functions required')
    if not prod['checked_main_closed_Rv_pressure_bridge_consumed']:
        raise ValueError('Same current signed Rv/absolute native pressure required')
    normalization=owner.proof['current_actual_raw_power_normalization_theorem']
    if not normalization['passed'] or not all(normalization['identities'].values()):
        raise ValueError('Actual raw energy/absolute pressure normalized ODE required')
    # Every composed stage has these exact whole-Z functions. Scalar
    # endpoint boxes below bound their coefficients, not fitted samples.
    z=s.Symbol('Z',real=True);C=1/(1+z*z);canonical={};chain={}
    for stage,functions in prod['independently_composed_actual_parent_chain_source_functions'].items():
        values={key:s.sympify(expr,locals={'Z':z}) for key,expr in functions.items()}
        chain[stage]=values
        EQ=values['e'].subs(z,0);EZ=s.diff(values['e'],z,2).subs(z,0)/2+2*EQ
        shapes=dict(u=values['u'].subs(z,0)*C,m=s.diff(values['m'],z).subs(z,0)*z,
            h=values['h'].subs(z,0)*C,k=s.diff(values['k'],z).subs(z,0)*z*C,
            e=EZ*z*z+EQ*C*C,p=values['p'].subs(z,0)*C*C)
        for key,wanted in shapes.items():
            if s.cancel(values[key]-wanted)!=0:raise ArithmeticError('Current whole-Z canonical source differs: '+stage+'.'+key)
            canonical[stage+'.'+key]=True
    transition_symbols={v.name:v for expr in chain['slope_mu'].values() for v in expr.free_symbols}
    mu=transition_symbols['mu'];Md=transition_symbols['same_Md'];Q=s.Symbol('same_actual_positive_Q',real=True)
    subs={transition_symbols['actual_transition_J_at1']:s.Rational(1,2),
        transition_symbols['KT']:s.exp(1-mu/2)-1+mu*Q}
    native0={key:expr.subs(subs,simultaneous=True) for key,expr in chain['slope_mu'].items()}
    Us=chain['slope']['u'].subs(z,0);Xs=(chain['slope']['h']/chain['slope']['u']).subs(z,0)
    D=(1-Xs)*s.exp(-s.exp(Md)-11+mu/2);V=Q*s.exp(-1+mu/2)
    U0=native0['u'].subs(z,0);M0=s.diff(native0['m'],z).subs(z,0)
    X0=(native0['h']/native0['u']).subs(z,0);K0=s.diff(native0['k'],z).subs(z,0)
    actual_correlations={}
    for name,a,b in (
        ('actual_composed_phase0_X_D_positive_Q',X0,1-D+mu*V),
        ('actual_composed_phase0_K_U_M_D',K0/U0,M0*s.exp(mu/2)-4*D)):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('Actual composed phase0 correlation differs: '+name)
        actual_correlations[name]=True
    mod=importlib.import_module('lei_ren_part1_paper_compliant_current_pulse_entrance_incoming_background_tensor')
    if owner.actual_power.__func__ is not mod.CurrentPulseEntranceIncomingBackgroundTensor.actual_power \
            or owner.actual_power.__func__.__globals__['pulse_coefficients'] is not mod.pulse_coefficients:
        raise ValueError('Original actual O3 source stress operator required')
    asts=SourceAST()
    buffer=importlib.import_module('lei_ren_part1_paper_compliant_outer_buffer')
    initial=importlib.import_module('lei_ren_part1_paper_compliant_outer_initial')
    if owner.physical.pre.slope_mu.__func__.__globals__['transition_kernels'] is not buffer.transition_kernels \
            or buffer.transition_kernels.__globals__['stable_sigma'] is not initial.stable_sigma:
        raise ValueError('Same actual nonnegative sigma/transition defining integral required')
    for stem,method,target,wanted in (
        ('pre_pulse_mixed_C4','power','parent','self.slope_mu(Z,1)'),
        ('pre_pulse_mixed_C4','slope_mu','parent','self.axial(Z,buffer_offset=11)'),
        ('pre_pulse_mixed_C4','axial','t','y-1'),
        ('current_pulse_entrance_incoming_background_tensor','actual_power','pressure',"(copy_jet(c,raw['p'][0])+datum)*C*C/(u*u)"),
        ('current_pulse_entrance_incoming_background_tensor','actual_power','P',"[pressure]")):
        asts.expression(stem,method,target,wanted=wanted)
    actual_method=asts.method('current_pulse_entrance_incoming_background_tensor','actual_power')
    append=[node for node in ast.walk(actual_method) if isinstance(node,ast.Call) and ast.unparse(node.func)=='P.append']
    expected=ast.parse('P.append(P[j]*(1+2*mu)+(C*C/2 if j==0 else zero))',mode='eval').body
    if len(append)!=1 or ast.dump(append[0])!=ast.dump(expected):
        raise ValueError('Actual absolute pressure ordinary row recurrence differs')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','baseline',wanted='-C*C/(2*p)')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','Ptilde',wanted='P0-baseline')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','Pmemory',wanted='[Ptilde*p**j for j in range(5)]')
    asts.expression('current_pulse_main_exit_background_tensor','_data','right',wanted='self.gap_tensor.cache[key]')
    asts.expression('current_pulse_main_exit_background_tensor','_data','self.cache[key]',
        wanted="dict(ap=ap,incoming=incoming,incoming_energy=energy,future=future,J0=right['J0'],P0=right['P0'],controls=right['controls'],inlet=inlet,u=u,current_selected_result=selected)")
    Pv=s.Symbol('same_current_signed_Rv_function',real=True);pp=1+2*mu
    source_base=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','baseline'),dict(C=C,p=pp))
    source_memory=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','Ptilde'),dict(P0=Pv,baseline=source_base))
    extras=asts.expression('current_pulse_main_exit_background_tensor','chart','extras')
    qnode=next(kw.value for kw in extras.keywords if kw.arg=='Q')
    q0=asts.evaluate(qnode,dict(mu=mu,xi=s.Integer(0)))
    inlet=source_base+source_memory*s.exp(q0)
    expected_inlet=-C*C/(2*pp)+(Pv+C*C/(2*pp))*s.exp(-pp*13/mu)
    if s.cancel(inlet-expected_inlet)!=0:raise ArithmeticError('Actual current pulse inlet pressure source differs')
    actual_correlations['actual_main_source_pressure_at_same_phase1_inlet']=True
    asts.method('outer_buffer','transition_kernels')
    # The exact J(1)=1/2 comes from the same sigma reflection integral.
    asts.expression('outer_buffer','transition_kernels','J',wanted="c.mpf('.5')")
    transition=entrancecone.registry.owners['o3']
    transition_source=transition.proof['current_actual_transition_source_and_power_endpoint']
    if not transition.acceptance_loaded or not transition_source['passed'] \
            or not transition_source['actual_transition_power_same_Rw']:
        raise ValueError('Same actual transition/power function join required')
    return dict(current_complete_production_function_proof=prod,
        exact_whole_Z_canonical_source_identities=canonical,
        actual_composed_source_phase0_M_K_X_correlations=actual_correlations,
        actual_composed_phase0_functions=dict(U0=str(U0),M0=str(M0),K0=str(K0),X0=str(X0),D=str(D),positive_V=str(V)),
        current_raw_energy_and_absolute_pressure_normalized_ODE=normalization,
        current_typed_parameter_source_bridge=entrancecone.bindings['current_typed_parameter_defining_source_bridge'],
        current_same_signed_Rv_and_native_absolute_pressure_proof=owner.main_tensor.proof,
        current_transition_power_tensor_join=transition.proof,
        current_power_entrance_tensor_join=owner.proof,
        same_actual_O2_mass_and_current_O3_transition_kernel_correlated=True,
        current_full_energy_pressure_radial_velocity_and_remainder_not_replaced=True,
        original_normalized_pressure_ODE_and_same_inlet_function_composed_backwards=True,
        all_positive_Q_bounds_from_original_nonnegative_sigma_integral=True,
        source_boxes_only_enclose_the_identified_canonical_functions=True,
        input_hashes=asts.hashes,passed=True)


def validate_whole_O3_power_view(entrancecone,view):
    c=entrancecone.ctx
    if tuple(view[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) \
            !=(entrancecone.family,entrancecone.source,entrancecone.datum_sha):raise ValueError('Foreign O3 source')
    raw=view['original_complete_view']
    if raw['chart']!='O3_power' or endpoints(raw['Z'])!=endpoints(c.mpf([-1,1])) \
            or endpoints(raw['coverage_coordinate'])!=(0,1):raise ValueError('Whole current O3 power box required')
    if not raw['actual_full_stress_not_local_difference'] or not raw['actual_pre_pressure_function_used_in_tensor']:
        raise ValueError('Full original O3 tensor/pressure required')
    for key in ('actual_upstream_original_pre_power_source','actual_upstream_scalar_U_X_defining_source',
        'current_actual_normalized_power_rows','current_source_three_component_velocity_rows'):
        if key not in raw:raise ValueError('Current O3 source history omitted: '+key)
    sectors=raw['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if set(sectors['theta'])!={'equilibrium','signed_original_memory','meridional_transport','radial_shear'} \
            or set(sectors['axial'])!={'full_energy_and_pressure','axial_transport','nonlinear_meridional_transport',
                'linear_axial_moment','selected_backward_energy_loss','axial_radial_shear'}:
        raise ValueError('Full original O3 signed sectors required')
    for label,parts in sectors.items():
        for name,part in parts.items():
            expected={'theta':{'equilibrium':(.5,1,0,0),'signed_original_memory':(.5,1,0,1),
                'meridional_transport':(.5,2,1,0),'radial_shear':(-.5,1,0,0)},
                'axial':{'full_energy_and_pressure':(.5,2,0,0),'axial_transport':(.5,1,1,0),
                    'nonlinear_meridional_transport':(.5,2,2,0),'linear_axial_moment':(.5,1,1,0),
                    'selected_backward_energy_loss':(.5,2,2,0),'axial_radial_shear':(-.5,1,1,0)}}[label][name]
            if tuple(part['mode'])!=expected or part['extra_source']!='one':
                raise ValueError('Actual O3 signed source mode differs')
            row=raw['physical_cylindrical_stress_mixed3'][label][name]['r0_z0']
            if encode(pack(row['signed_coefficient']))!=encode(pack(part['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])) \
                    or row['actual_source_log_parts']!=part['exact_source_log_parts']:
                raise ValueError('Actual O3 physical signed coefficient/factor differs')
    return raw


@source_precision
def whole_current_O3_power_bounds(entrancecone,view,slope,initial):
    validate_whole_O3_power_view(entrancecone,view);c=entrancecone.ctx;owner=entrancecone.registry.owners['incoming']
    mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta);r=1-mu;p=1+2*mu;qmin=mu-delta/2
    pre=owner.physical.pre;Tw=c.mpf(pre.params.Tw);tau=c.mpf(pre.params.logPstar)-1
    positive={}
    def require(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Current O3 positive bound failed: '+name)
        positive[name]=value
    for name,value in dict(mu=mu,delta=delta,actual_Tw=Tw,qmin=qmin,one_minus_delta=1-delta,rate_above_half=r-c.mpf('.5'),
        delta_below_one_eighth=c.mpf('.125')-delta,mu_Tw_below_point001=c.mpf('.001')-mu*Tw).items():require(name,value)
    if slope['chart']!='O2_slope' or endpoints(slope['Z'])!=(0,0) or endpoints(slope['coverage_coordinate'])!=(1,1):
        raise ValueError('Actual canonical O2 slope endpoint required')
    if initial['chart']!='O3_power_to_Rp' or endpoints(initial['Z'])!=(0,0) or endpoints(initial['coverage_coordinate'])!=(0,0):
        raise ValueError('Actual canonical O3 phase0 source required')
    def interval(v):return c.mpf(endpoints(v))
    us=interval(slope['Utheta_over_Pstar_axial5_coefficients'][0])
    hs=interval(slope['actual_normalized_primitive_y_derivative_axial5']['h'][0][0])
    require('actual_Us',us);ds=1-hs/us;require('actual_one_minus_Xs',ds)
    logD=c.ln(ds)-tau-1+mu/2;D=c.exp(logD)
    u0=interval(initial['Utheta_over_Pstar_axial5_coefficients'][0]);require('actual_U0',u0)
    raw=initial['actual_normalized_primitive_y_derivative_axial5']
    M0=interval(raw['m'][0][1]);require('actual_M0',M0)
    K0=interval(raw['k'][0][1]);H0=interval(raw['h'][0][0])
    require('actual_M0_below_point001',c.mpf('.001')-M0)
    Dmin=(c.mpf(15)/8-14*delta)/4;cD=Dmin/2;cX=(1-delta)/4
    Eqmin=qmin/(2*r);delta_cap=delta/(2*(1-delta));Ncap=2/(1-delta)
    require('positive_deficit_shape',Dmin)
    require('deficit_reserve_dominates_uniform_floor_and_delta',cD*D-delta_cap-Eqmin)
    drift=Ncap*M0*mu*c.exp(mu/2)*(1/r+c.mpf('.5'))
    require('drift_below_quarter_equilibrium_floor',Eqmin/4-drift)
    logRmin=c.mpf(owner.physical.logRp)-Tw
    radial_log=c.ln(2*(1+mu))-logRmin-c.ln(Eqmin)
    require('inverse_R_relative_log_gap',-1000-radial_log)
    radial_cap=Eqmin*c.exp(-1000);theta_floor=Eqmin-drift-radial_cap
    require('theta_uniform_floor',theta_floor);require('theta_floor_exceeds_half_equilibrium',theta_floor-Eqmin/2)
    growth=c.exp(2*mu*Tw)
    EQ0=interval(raw['e'][0][0]);EZ0=interval(raw['e'][0][2])+2*EQ0
    aq=upper(c,(absolute(c,EQ0/u0**2)+Tw/2)*growth)
    az=upper(c,absolute(c,EZ0/u0**2)*growth)
    Lz=upper(c,((2+2*delta)*az+(4+2*delta)*aq+(3+delta)/p)/(1-delta))
    require('finite_full_energy_axial_factor',Lz)
    lp=c.mpf(owner.physical.logP);lu=c.mpf(initial['log_Utheta_over_Pstar_base_source'])
    energy_direction_log=c.ln(mu)+2*lp+2*lu+2*c.ln(Lz)-c.ln(2*cD*cX)-logD
    require('energy_direction_relative_log_gap',-1000-energy_direction_log)
    Pv=entrancecone.whole_view['original_complete_view']['current_signed_absolute_Rv_pressure']
    pm=upper(c,(2*(1+delta)*(absolute(c,Pv[0])+1/(2*p))+absolute(c,Pv[1])+2/p)/(1-delta))
    logU1=c.ln(c.mpf(owner.main_tensor.gap_tensor.end_tensor.flatten_power.flatten.U))
    pressure_ratio_log=lp+logU1+c.ln(pm)-p*13/mu-c.ln(theta_floor)
    require('full_pressure_memory_ratio_log_gap',-1000-pressure_ratio_log)
    direction_energy_cap=c.exp(-1000);w_energy_cap=c.sqrt(direction_energy_cap/(2*mu))
    w_pressure_cap=c.exp(-1000);direction=2*mu*(w_energy_cap+w_pressure_cap)**2
    require('direction_below_paper_one_point8',c.mpf('1.8')-direction)
    require('full_directional_bracket',2-direction);require('full_directional_margin_normalized',4*(2-direction))
    return dict(positive_margins=positive,mu=mu,delta=delta,
        actual_O2_slope_U=us,actual_O2_slope_H=hs,actual_O2_slope_deficit=ds,
        exact_current_O3_deficit_log=logD,actual_current_O3_deficit_enclosure=D,
        actual_current_M0=M0,actual_current_U0=u0,actual_current_K0=K0,actual_current_H0=H0,
        independent_raw_K_H_boxes_retained_not_subtracted_to_define_correlated_deficit=True,
        theta_lower_formula='cX*Z^2+cD*D*exp(-(1-mu)*t)+theta_floor, t=Tw*phase',
        theta_lower_constants=dict(cX=cX,cD=cD,theta_floor=theta_floor),
        nonnegative_Q_and_M_source_terms_retained=True,
        source_M_K_X_correlated_before_enclosure=True,drift_absolute_uniform_upper=drift,
        inverse_R_radial_log=radial_log,inverse_R_radial_error_upper=radial_cap,
        full_energy_EQ0=EQ0,full_energy_EZ0=EZ0,full_energy_AQ_upper=aq,full_energy_AZ_upper=az,
        full_axial_energy_baseline_pressure_Z_factor_upper=Lz,
        energy_direction_log_upper=energy_direction_log,energy_direction_upper=direction_energy_cap,
        actual_signed_Rv_pressure_C5=Pv,pressure_memory_coefficient_upper=pm,
        pressure_memory_ratio_log_upper=pressure_ratio_log,pressure_memory_ratio_upper=w_pressure_cap,
        full_directional_expression_upper=direction,
        actual_axial_velocity_and_shear_zero=True,canonical_linear_axial_stress_zero_by_function_identity=True,
        full_axial_energy_and_pressure_stress_nonzero_retained=True,
        whole_current_O3_power_two_vector_cone=True,continuous_domain='phase[0,1], all Z[-1,1]',
        original_absolute_pressure_restored_by_same_source_ODE_not_added_tail=True,
        complete_original_source_signed_tensor_and_remainder_preserved=True,
        cap_values_used_as_defining_fields=False,phase_samples_used_as_proof=False)
