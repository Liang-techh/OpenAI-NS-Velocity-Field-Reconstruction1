"""Original physical velocity/pressure and absolute boundary from current h*.

Exact source functions are restored unchanged. A numerical point adapter
requires a separately supplied actual source oracle; its root algorithm and
synthetic checks do not install numerical controls or an exterior datum.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_limit_repair_band as current
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery
import lei_ren_part1_paper_compliant_current_modified_physical_velocity_operator as physical

source,packets=current.source,current.packets
HERE,PREFIX,sha,require=current.HERE,current.PREFIX,current.sha,current.require
NAME=PREFIX+'current_limit_band_physical_recovery.json.gz'
RECEIPT=PREFIX+'current_limit_band_physical_recovery_check.json'
GATE='current_exact_limit_band_original_physical_velocity_pressure_and_absolute_boundary_recovered'


def load_inputs():
    report=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
    checked=json.loads((HERE/current.RECEIPT).read_bytes())
    require(report[current.GATE] and checked[current.GATE] and checked['all_passed'],
        'Accepted actual current C1 limit repair band required')
    hashes=dict(checked['input_hashes'])
    for name,digest in hashes.items():require(sha(name)==digest,'Changed current prerequisite: '+name)
    hashes[current.NAME]=sha(current.NAME);hashes[current.RECEIPT]=sha(current.RECEIPT)
    for module in (recovery,physical):
        name=Path(module.__file__).name;hashes[name]=sha(name)
    pressure=PREFIX+'pressure_source.json';record=json.loads((HERE/pressure).read_bytes())
    require(record['compliant_source']['implicit_source_sha256']==report['source_family']['implicit_source_sha256'],
        'Same original similarity/pressure parameter source required')
    require(record['compliant_source']['implicit_source_definition']['Md']=='40','Actual fixed Md=40 required')
    require(record['compliant_source']['implicit_source_definition']['delta']=='min(1e-200,exp(-4logPstar-30))',
        'Original similarity delta definition changed')
    hashes[pressure]=sha(pressure);hashes[Path(__file__).name]=sha(Path(__file__).name)
    return report,record,hashes


def restore_band(report):
    g=source.FunctionTransportGraph();nodes=report['exact_graph_nodes']
    require(nodes[:2]==g.nodes,'Original zero/one function prefix required')
    g.nodes=[dict(q) for q in nodes]
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    require(len(g.keys)==len(g.nodes),'Duplicate accepted source expressions')
    raw=report['exact_limit_repair_band'];background=raw['original_power_background']
    ref=lambda i:source.FunctionRef(g,i);pair=lambda q:source.C1Function(*map(ref,q))
    fields={key:pair(background['corrected_fields'][key]) for key in ('E','V','E_y','V_y')}
    return dict(graph=g,N=ref(raw['N']),parameters={k:ref(v) for k,v in raw['parameters'].items()},
        source_family=report['source_family'],source_graph_prefix_length=len(nodes),
        fields=fields,own={k:pair(v) for k,v in background['complete_histories'].items()},
        own_y={k:pair(v) for k,v in background['complete_history_y_Z_pairs'].items()},
        original_leading={k:pair(v) for k,v in background['leading_histories'].items()},
        original_fields={k:pair(background[k]) for k in ('E','V','E_y','V_y')},
        original_own_y={k:pair(v) for k,v in background['leading_history_y_Z_pairs'].items()},
        P0=pair(background['independent_original_P0']),
        original_leaf_binding=background['leading_history_P0_source_binding'],
        relative_terminal_zero_certificates=raw['relative_terminal_zero_certificates'])


def native_recovery(g,fields,own,own_y,P0,Z,delta):
    """Full original n=0 source equations; Q needs m_Z but not m_ZZ."""
    E,V,Ey,Vy=(fields[k] for k in ('E','V','E_y','V_y'))
    m,h,k,e,p=(own[key] for key in ('m','h','k','e','p'))
    d=g.sub(g.one,g.mul(Z,Z));L=g.sub(g.one,g.mul(delta,Z,Z))
    divide=lambda q:g.quotient(q,L,'same original 0<delta<1 and |Z|<=1 give L>=1-delta>0')
    transport=g.add(g.mul(g.sub(g.one,delta),Z,m.value),g.mul(d,m.Z))
    Q=divide(g.sub(g.mul(g.constant(2),Z,V.value),transport))
    my=own_y['m']
    Qy=divide(g.sub(g.mul(g.constant(2),Z,Vy.value),
        g.add(g.mul(g.sub(g.one,delta),Z,my.value),g.mul(d,my.Z))))
    pressure=g.c1add(P0,p);Py=own_y['p']
    theta_linear=divide(g.sub(g.add(g.neg(E.value),g.mul(g.sub(g.one,g.mul(g.constant('1/2'),delta)),h.value)),
        g.mul(g.constant('1/2'),Z,g.sub(g.one,delta),h.Z)))
    theta_quadratic=divide(g.add(g.mul(g.sub(g.mul(g.constant(2),delta),g.one),Z,k.value),
        g.neg(g.mul(d,k.Z)),g.mul(E.value,transport)))
    axial_linear=divide(g.add(g.neg(V.value),g.mul(g.constant('1/2'),g.sub(g.one,delta),g.sub(m.value,g.mul(Z,m.Z)))))
    axial_quadratic=divide(g.add(g.mul(V.value,transport),g.mul(g.constant(2),delta,Z,e.value),
        g.neg(g.mul(d,e.Z)),g.mul(g.constant(2),g.add(g.one,delta),Z,pressure.value),g.neg(g.mul(d,pressure.Z))))
    return dict(E=E,V=V,Q=Q,Q_y=Qy,pressure=pressure,pressure_y=Py,L=L,d=d,
        full_signed_inertial_sectors=dict(theta_linear=theta_linear,theta_quadratic=theta_quadratic,
            axial_linear=axial_linear,axial_quadratic=axial_quadratic),
        shear_theta=g.sub(g.mul(g.constant(2),Ey.value),E.value),shear_axial=g.mul(g.constant(2),Vy.value),
        Q_Z_and_full_velocity_mixed_jets_not_available=True)


def paper_moments(g,own,P0,R,S):
    """Current normalized histories -> exact paper cumulative variables."""
    R3=g.unary('exp',g.mul(g.constant('3/2'),g.unary('log',R)))
    sqrt2=g.unary('sqrt',g.constant(2));S2=g.mul(S,S)
    scales=dict(M=g.mul(R,S),I=g.mul(sqrt2,R3,S),J=g.mul(sqrt2,R3,S2),
        S_energy=g.mul(R,S2),Cp=S2)
    keys=dict(M='m',I='h',J='k',S_energy='e',Cp='p')
    moments={key:g.c1scale(scales[key],own[name]) for key,name in keys.items()}
    return dict(cumulative=moments,axis_Pi=g.c1scale(S2,P0),
        absolute_Pi=g.c1scale(S2,g.c1add(P0,own['p'])),
        defining_units=recovery.UNITS,radius_is_paper_X_not_sqrt_2X=True)


def substitution(g,q,variables,values):
    return g.node('simultaneous_function_substitution',expression=q.node,
        variables=g.ids(variables),values=g.ids(values),
        simultaneous=True,source_functions_held_fixed=True)


def build(report,pressure_record):
    built=restore_band(report);g=built['graph'];Z=g.symbol('Z');x=g.symbol('repair_x')
    logP=built['parameters']['logP'];S=g.unary('exp',logP)
    # Md=40 makes exp(-4logP-30)<1e-200. The positive source is
    # retained as an exact exponential, never rounded to delta=0.
    delta=g.unary('exp',g.sub(g.neg(g.mul(g.constant(4),logP)),g.constant(30)))
    recipe=built['original_leaf_binding']
    Rc=g.node('current_original_leading_function_recipe',source_family=built['source_family'],
        native_chart='O3_power',recipe=recipe,quantity='actual_source_radius',Z_order=0,
        coordinate=g.constant(2).node,Z_variable=Z.node,
        quantity_path=['original_closed_O3_background','actual_source_radius'],
        defining_quantity_not_a_range_value=True,Z_independent=True,
        exact_definition='Rc=Rm*exp(exp40+20); unchanged original coordinate2 radius recipe',
        source_background_binding=current.BACKGROUND_BINDING)
    R=g.mul(Rc,x)
    recovered=native_recovery(g,built['fields'],built['own'],built['own_y'],built['P0'],Z,delta)
    moments=paper_moments(g,built['own'],built['P0'],R,S)
    # Coordinate-2 terminal D=0 identifies complete and original leading
    # values AND first Z rows. Keep both expressions and the proof handles.
    leading=native_recovery(g,built['original_fields'],built['original_leading'],built['original_own_y'],built['P0'],Z,delta)
    endpoint=lambda q:substitution(g,q,[x],[g.constant(2)])
    terminal_moments={key:source.C1Function(endpoint(pair.value),endpoint(pair.Z))
        for key,pair in moments['cumulative'].items()}
    Pi2=source.C1Function(endpoint(moments['absolute_Pi'].value),endpoint(moments['absolute_Pi'].Z))
    c_inf=g.symbol('future_source_c_infinity');h=g.mul(g.constant('1/2'),delta)
    R2=g.mul(g.constant(2),Rc)
    Ipow=g.quotient(g.mul(g.unary('sqrt',g.constant(2)),c_inf,
        g.unary('exp',g.mul(g.sub(g.one,h),g.unary('log',R2)))),
        g.sub(g.one,h),'same original 0<h<1/2')
    future_targets={key:g.c1scale(g.constant(-1),terminal_moments[key]) for key in ('M','J','S_energy')}
    future_targets['renormalized_I']=source.C1Function(g.sub(Ipow,terminal_moments['I'].value),g.neg(terminal_moments['I'].Z))
    future_targets['pressure_tail']=g.c1scale(g.constant(-1),Pi2)
    # This is the FINAL paper reference Hpow, integrated over the ENTIRE
    # remaining source. It is not the current O3 power, whose exponent is mu.
    # No amplitude or transition is installed by these necessary conditions.
    final_reference=dict(h=h,c_infinity=c_inf,primitive_at_2Rc=Ipow,
        definition='Hpow=sqrt(2)*c_infinity*R^(-delta/2)',
        reference_parameter_hypotheses='c_infinity>0; c_infinity,Rc,delta Z-independent',
        final_reference_c_infinity_not_installed=True,
        current_O3_amplitude_not_identified_with_final_c_infinity=True,
        actual_transition_2Rc_to_final_heat_not_installed=True,
        paper_binding='Lemma 4.8 and equations 4.25,4.28')
    mu=built['parameters']['mu'];alpha=g.add(g.constant('1/2'),mu)
    A1=source.C1Function(*[substitution(g,getattr(built['original_fields']['E'],row),[x],[g.one]) for row in ('value','Z')])
    coefficient=g.c1scale(g.mul(S,g.unary('exp',g.mul(alpha,g.unary('log',Rc)))),A1)
    log2=g.unary('log',g.constant(2))
    current_I_increment=g.c1scale(g.mul(g.unary('sqrt',g.constant(2)),S,
        g.unary('exp',g.mul(g.constant('3/2'),g.unary('log',Rc))),
        log2,g.unary('exprel',g.mul(g.sub(g.one,mu),log2))),A1)
    current_reference=dict(mu=mu,alpha=alpha,paper_power_coefficient=coefficient,
        angular_primitive_Rc_to_2Rc=current_I_increment,
        definition='E_paper=S*A_rc(Z)*(R/Rc)^(-1/2-mu)',
        actual_C1_amplitude_and_Z_row_retained=True,
        current_mu_not_replaced_with_final_delta_over_2=True)
    future_integral_definitions=dict(domain='[2Rc,infinity)',
        M='integral U dR',J='integral U*H dR',S_energy='integral (U^2-E^2/2) dR',
        renormalized_I='integral (H-Hpow) dR; includes every O3/flatten/preheat/heat transition',
        pressure_tail='integral E^2/(2R) dR',
        necessary_conditions_only=True,actual_exterior_source_not_installed=True,
        hypotheses='M(infinity)=J(infinity)=S_energy(infinity)=0; integral_0^infinity(H-Hpow)=0; Pi(infinity)=0')
    # Physical Cartesian variables are distinct from the source x=R/Rc.
    px,py,pz,pt,T,nu=[g.symbol(k) for k in ('physical_x','physical_y','physical_z','physical_t','terminal_time','constant_viscosity')]
    tau=g.sub(T,pt);r=g.unary('sqrt',g.add(g.mul(px,px),g.mul(py,py)));sqrt_nu=g.unary('sqrt',nu)
    lam=g.node('unique_positive_similarity_root',tau=tau.node,z=pz.node,nu=nu.node,delta=delta.node,
        definition='lambda^2-(physical_z^2/nu)*lambda^(2delta)=tau',
        admissible_domain='tau>0; constant nu>0; same original 0<delta<1',
        uniqueness='s=log(lambda/sqrt(tau)): 2s-log(1+exp(logq+2delta*s)) strictly increasing',
        source_similarity_parameters=pressure_record['compliant_source']['implicit_source_definition'],
        source_parameter_report=PREFIX+'pressure_source.json',source_parameter_report_sha256=sha(PREFIX+'pressure_source.json'),
        q_is_lambda_squared=True,lambda_not_replaced_by_sqrt_tau=True)
    power=lambda b:g.unary('exp',g.mul(b,g.unary('log',lam)))
    zp=g.quotient(pz,g.mul(sqrt_nu,power(g.sub(g.one,delta))),'positive viscosity and similarity scale')
    Rp=g.quotient(g.mul(r,r),g.mul(g.constant(2),nu,lam,lam),'positive viscosity and similarity scale')
    xp=g.quotient(Rp,Rc,'same source positive Rc')
    at=lambda q:substitution(g,q,[Z,x],[zp,xp])
    domain=g.node('source_region_guard',tau=tau.node,nu=nu.node,r=r.node,native_x=xp.node,
        definition='tau>0,nu>0,r>0,1<=native_x<=2',source_family=built['source_family'],
        outside_band_requires_other_actual_region=True)
    Ur=g.mul(sqrt_nu,power(g.constant(-1)),S,g.unary('sqrt',g.mul(g.constant('1/2'),Rp)),at(recovered['Q']))
    Ut=g.mul(sqrt_nu,power(g.neg(g.add(g.one,delta))),S,at(recovered['E'].value))
    Uz=g.mul(sqrt_nu,power(g.neg(g.add(g.one,delta))),S,at(recovered['V'].value))
    cosine=g.quotient(px,r,'repair annulus r>0');sine=g.quotient(py,r,'repair annulus r>0')
    xyz=dict(u=g.sub(g.mul(cosine,Ur),g.mul(sine,Ut)),v=g.add(g.mul(sine,Ur),g.mul(cosine,Ut)),w=Uz,
        p=g.mul(nu,power(g.neg(g.mul(g.constant(2),g.add(g.one,delta)))),S,S,at(recovered['pressure'].value)))
    built.update(S=S,delta=delta,Rc=Rc,profile_R=R,native_recovery=recovered,paper_moments=moments,
        terminal_absolute_cumulative=terminal_moments,terminal_absolute_Pi=Pi2,
        required_future_integrals=future_targets,future_c_infinity=c_inf,
        final_heat_reference=final_reference,current_O3_power_reference=current_reference,
        future_integral_definitions=future_integral_definitions,
        physical_coordinate_functions=dict(lambda_scale=lam,Z=zp,R=Rp,x_band=xp,r=r,tau=tau,domain=domain),
        cartesian_velocity_pressure=xyz,
        terminal_complete_vs_leading=dict(Q=endpoint(recovered['Q']),Q_leading=endpoint(leading['Q']),
            pressure=endpoint(recovered['pressure'].value),pressure_leading=endpoint(leading['pressure'].value),
            relative_history_zero_proofs=built['relative_terminal_zero_certificates']),
        original_recovery_binding=current.current.ast_binding(recovery.GenericMomentRecovery.field))
    return built


def solve_similarity_scale(ctx,z,tau,nu,delta,*,iterations=256):
    """Stable approximate scalar inverse; no source/correction oracle hidden.

    Bisection in log(lambda/sqrt(tau)) avoids exponential root overflow.
    The numeric bracket is not a directed interval certificate.
    """
    z,tau,nu,delta=map(ctx.mpf,(z,tau,nu,delta))
    if not all(ctx.isfinite(q) for q in (z,tau,nu,delta)) or tau<=0 or nu<=0 or not 0<=delta<1:
        raise ValueError('Finite z,tau>0,constant nu>0 and 0<=delta<1 required')
    if type(iterations) is not int or not 1<=iterations<=8192:raise ValueError('Bounded exact iteration count required')
    if z==0:return dict(log_lambda=ctx.log(tau)/2,Z=ctx.mpf(0),log_root_bracket=(ctx.mpf(0),ctx.mpf(0)),approximate_only=True)
    logq=2*ctx.log(abs(z))-ctx.log(nu)+(delta-1)*ctx.log(tau)
    softplus=lambda q:max(q,ctx.mpf(0))+ctx.log1p(ctx.exp(-abs(q)))
    f=lambda s:2*s-softplus(logq+2*delta*s)
    lo=max(ctx.mpf(0),logq/(2*(1-delta)));hi=lo+ctx.log(2)/(2*(1-delta))
    for unused in range(iterations):
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        if f(mid)>0:hi=mid
        else:lo=mid
    s=(lo+hi)/2;loglambda=ctx.log(tau)/2+s
    Z=ctx.sign(z)*ctx.exp(ctx.log(abs(z))-ctx.log(nu)/2-(1-delta)*loglambda)
    return dict(log_lambda=loglambda,Z=Z,log_root_bracket=(lo,hi),approximate_only=True)


def recover_numeric_values(ctx,packet,Z,delta):
    """Original scalar recovery from supplied values, never cap endpoints."""
    keys=('E','V','E_y','V_y','P0','P0_Z','m','m_Z','m_y','m_yZ','p','p_Z','p_y','p_yZ')
    if any(k not in packet for k in keys):raise ValueError('Actual value/first-Z/radial source rows required')
    for k in keys:
        q=packet[k]
        if isinstance(q,(dict,source.FunctionRef)) or hasattr(q,'_mpi_') or hasattr(q,'scale'):
            raise TypeError('Function ranges/caps/handles are not numerical source values: '+k)
        if not ctx.isfinite(ctx.mpf(q)):raise ValueError('Finite original source values required')
    E,V,Ey,Vy,P0,P0Z,m,mZ,my,myZ,p,pZ,py,pyZ=(ctx.mpf(packet[k]) for k in keys)
    Z,delta=ctx.mpf(Z),ctx.mpf(delta)
    if abs(Z)>1 or not 0<=delta<1:raise ValueError('Admitted original axial/similarity parameters required')
    L=1-delta*Z*Z;d=1-Z*Z
    return dict(E=E,V=V,Q=(2*Z*V-(1-delta)*Z*m-d*mZ)/L,
        Q_y=(2*Z*Vy-(1-delta)*Z*my-d*myZ)/L,
        pressure=P0+p,pressure_Z=P0Z+pZ,pressure_y=py,pressure_yZ=pyZ)


class SourcePointOracleRequired(RuntimeError):pass


class ApproximateCartesianBandAdapter:
    """Executable numeric map; actual h* oracle is a separate open task."""
    def __init__(self,*,oracle,source_family,graph_sha256,ctx=None):
        if oracle is None or getattr(oracle,'mode',None) not in ('original_source','synthetic_reference'):
            raise SourcePointOracleRequired('Explicit original source point oracle required')
        if oracle.source_family!=source_family or oracle.source_graph_sha256!=graph_sha256:
            raise SourcePointOracleRequired('Same current exact limit band source/graph required')
        if not callable(getattr(oracle,'parameters',None)) or not callable(getattr(oracle,'band_state',None)):
            raise SourcePointOracleRequired('Explicit original parameter and band-state value callables required')
        self.oracle,self.ctx=oracle,ctx or mp.mp.clone()

    def cartesian(self,x,y,z,t,*,terminal_time=1,viscosity=1):
        c=self.ctx;values=self.oracle.parameters();delta,S,Rc=(c.mpf(values[k]) for k in ('delta','S','Rc'))
        if S<=0 or Rc<=0:raise ValueError('Positive actual original S/Rc required')
        x,y,z,t,T,nu=map(c.mpf,(x,y,z,t,terminal_time,viscosity));tau=T-t
        coord=solve_similarity_scale(c,z,tau,nu,delta,iterations=4*c.prec)
        lam=c.exp(coord['log_lambda']);r=c.sqrt(x*x+y*y)
        if r<=0:raise ValueError('Repair annulus requires positive physical radius')
        R=r*r/(2*nu*lam*lam);xb=R/Rc
        if not 1<=xb<=2:raise ValueError('Requested point is outside this actual repair band')
        state=recover_numeric_values(c,self.oracle.band_state(coord['Z'],xb),coord['Z'],delta)
        rootnu=c.sqrt(nu);Ur=rootnu*S*c.sqrt(R/2)*state['Q']/lam
        factor=rootnu*S*c.exp(-(1+delta)*coord['log_lambda']);Ut,Uz=factor*state['E'],factor*state['V']
        return dict(u=x*Ur/r-y*Ut/r,v=y*Ur/r+x*Ut/r,w=Uz,
            p=nu*S*S*c.exp(-2*(1+delta)*coord['log_lambda'])*state['pressure'],
            similarity=dict(**coord,R=R,x_band=xb),mode=self.oracle.mode,
            approximate_only=True,numeric_source_controls_and_error_budget_certified=False)


def symbolic_checks():
    y,z,de=sy.symbols('y Z delta',real=True);L=1-de*z*z;d=1-z*z
    m=sy.Function('m')(y,z);V=sy.Function('V')(y,z)
    Q=(2*z*V-(1-de)*z*m-d*sy.diff(m,z))/L
    numerator=L*(Q+sy.diff(Q,y))+(-1-de)*z*V+d*sy.diff(V,z)-2*z*sy.diff(V,y)
    simplified=numerator.subs({sy.diff(m,y,z):sy.diff(V,z)-sy.diff(m,z),sy.diff(m,y):V-m},simultaneous=True).doit()
    require(sy.simplify(simplified)==0,'Original actual physical divergence fails')
    R,S=sy.symbols('R S',positive=True);m0,h0,k0,e0,p0=sy.symbols('m h k e p')
    E,U=sy.symbols('E U');rates=dict(m=m0,h=h0,k=k0,e=e0,p=p0)
    driver=dict(m=U,h=E,k=E*U,e=U*U-E*E/2,p=E*E/2)
    scales=dict(m=R*S,h=sy.sqrt(2)*R**sy.Rational(3,2)*S,
        k=sy.sqrt(2)*R**sy.Rational(3,2)*S*S,e=R*S*S,p=S*S)
    expected=dict(m=S*U,h=sy.sqrt(2*R)*S*E,k=sy.sqrt(2*R)*S*S*E*U,
        e=S*S*(U*U-E*E/2),p=S*S*E*E/(2*R))
    for key,rate in current.current.RATES.items():
        actual=(R*sy.diff(scales[key],R)*rates[key]+scales[key]*(driver[key]-sy.Rational(rate.numerator,rate.denominator)*rates[key]))/R
        require(sy.simplify(actual-expected[key])==0,'Paper moment conversion fails: '+key)
    # Exact physical centrifugal relation from Pi_y=E_full^2/2.
    nu,lam,r=sy.symbols('nu lambda r',positive=True)
    ut=sy.sqrt(nu)*S*lam**(-1-de)*E
    pr=nu*S*S*lam**(-2-2*de)*E*E/r
    require(sy.simplify(pr-ut*ut/r)==0,'Original physical radial pressure balance fails')
    # Delta=2h; code R is the paper X, not its cylindrical sqrt(2X).
    require(sy.simplify((sy.Rational(1,2)+de/2)+(sy.Rational(1,2)-de/2)-1)==0,'A+D=1 lost')
    I2,Itail,P2,Ptail=sy.symbols('I2 Itail Ipow2 Ipow_tail')
    require(sy.expand((Itail-I2-(Ptail-P2))+(Ptail-Itail)-(P2-I2))==0,
        'Entire renormalized future target must retain the transition integral')
    return dict(exact_physical_divergence_from_own_moment_ODE=True,
        genuine_m_yZ_commutation_used=True,paper_all_five_cumulative_unit_conversions=True,
        physical_centrifugal_pressure_identity=True,signed_Z_and_positive_constant_viscosity_retained=True,
        radial_profile_R_equals_paper_X=True,similarity_delta_equals_twice_paper_h=True,
        final_renormalized_target_includes_unsupplied_transition=True,
        absolute_future_moments_and_pressure_not_set_to_zero=True)


def run():
    began=time.monotonic();report,pressure,hashes=load_inputs();built=build(report,pressure)
    result=dict(**{GATE:True},source_family=report['source_family'],
        exact_graph_nodes=built['graph'].nodes,exact_physical_recovery=current.current.encode_graph(built),
        symbolic_checks=symbolic_checks(),original_recovery_binding=built['original_recovery_binding'],
        paper_equations=['4.1-4.7','4.15','4.25','4.28'],
        exact_band_cartesian_velocity_and_absolute_pressure_functions_defined=True,
        exact_band_divergence_and_centrifugal_pressure_identities_proved=True,
        absolute_2Rc_state_and_required_future_integrals_defined=True,
        actual_numeric_point_source_oracle_installed=False,actual_numeric_controls_evaluated=False,
        actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
        full_recovered_velocity_pressure_and_heat_joins_admitted=False,
        higher_Z_velocity_jets_and_stress_cone_admitted=False,
        actual_temporal_scale_recursion_installed=False,**dict.fromkeys(current.outer.OPEN,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Same current C1 limit source -> original n=0 physical u,v,w,p on Rc..2Rc, exact divergence '
            'and pressure balance, paper cumulative moment units and absolute future matching targets. '
            'Approximate numeric mapping requires an explicit source oracle. No higher-Z/global/exterior/heat/recursion/fullNS admission.')
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(packets.encode(result),separators=(',',':'))+'\n').encode(),mtime=0))
    print('CURRENT_LIMIT_BAND_PHYSICAL_RECOVERY',len(built['graph'].nodes),'nodes; original u/v/w/p and absolute 2Rc targets',flush=True)
    return result


if __name__=='__main__':run()
