"""Independent full kernels, common-unit source recovery and local Duhamel."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_long_reshape_finite_N as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as scalar
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def flow(c):return MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def contains(row,value,allowance='1e-60',max_log=1000):
    lo,hi=ep(row.finite_interval(max_log=max_log)+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def symbolic():
    y=sy.symbols('y');l=sy.Function('l')(y);E=sy.Function('E')(y);V,S=sy.symbols('V S',nonzero=True)
    names=('theta','theta_z','pressure','swirl','mean','axial');shapes={name:sy.Function(name)(y) for name in names}
    ode=dict(theta=1-(sy.Rational(3,2)+l)*shapes['theta'],theta_z=V-(sy.Rational(3,2)+l)*shapes['theta_z'],
        pressure=1-2*l*shapes['pressure'],swirl=1-(1+2*l)*shapes['swirl'],
        mean=V-shapes['mean'],axial=V*V-shapes['axial'])
    replace={sy.diff(E,y):l*E,**{sy.diff(shapes[name],y):value for name,value in ode.items()}}
    raw=dict(m=shapes['mean'],h=E*shapes['theta'],k=E*shapes['theta_z'],
        e=shapes['axial']/S**2-E**2*shapes['swirl']/2,p=E**2*shapes['pressure']/2)
    expected=dict(m=V-raw['m'],h=E-sy.Rational(3,2)*raw['h'],k=E*V-sy.Rational(3,2)*raw['k'],
        e=V**2/S**2-raw['e']-E**2/2,p=E**2/2)
    for name in current.RATES:assert sy.simplify(sy.diff(raw[name],y).xreplace(replace)-expected[name])==0
    B,T,ds=sy.symbols('B T ds')
    assert sy.simplify(1-2*(sy.Rational(1,10)-B*ds/T)-(sy.Rational(4,5)+2*B*ds/T))==0
    lam,width,suffix=sy.symbols('lam width suffix',positive=True)
    mass=(1-sy.exp(-lam*width))/lam
    assert sy.simplify(sy.integrate(sy.exp(-lam*(width-y)),(y,0,width))-mass)==0
    assert sy.limit(mass,lam,0)==width
    z,d,m=sy.symbols('z d m');bs=sy.symbols('b0:6')
    exact=sy.exp(m*d*sum(bs[j]*z**j for j in range(6)))
    polynomials=[sy.simplify(sy.diff(exact,z,n).subs(z,0)/sy.factorial(n)/sy.exp(m*d*bs[0])) for n in range(6)]
    return dict(passed=True,all_five_original_raw_y_ODEs_and_Pstar_units_verified=True,
        original_correlated_C_identity_verified=True,full_finite_own_rate_mass_and_zero_rate_limit_verified=True,
        independent_full_integrand_Z_coefficients_through5_verified=True),polynomials


def kernel_references(polynomials):
    c=MPIntervalContext();c.dps=150;f=flow(c);p=mp.mp.clone();p.dps=85
    z,d,m=sy.symbols('z d m');bs=sy.symbols('b0:6');compared=0
    for sign in (-1,1):
        coeff=[p.mpf(sign)*p.mpf('.02'),p.mpf('.003'),p.mpf('-.002'),p.mpf('.001'),p.mpf('-.0007'),p.mpf('.0002')]
        B=current.IntervalTaylor(c,[c.mpf(v) for v in coeff]);T=20;rho=p.mpf('.5');y=T*rho
        for kind,(k,m0) in dict(theta=('1.6',1),pressure=('.2',2),swirl=('1.2',2)).items():
            rows,proof=current.finite_kernel(f,B,c.mpf(T),c.mpf(y),kind)
            for n in range(6):
                formula=polynomials[n].subs({m:m0,**{b:sy.Rational(str(v)) for b,v in zip(bs,coeff)}})
                polynomial=sy.lambdify(d,formula,'mpmath')
                def integrand(t):
                    delta=scalar.flat_step(p,rho)-scalar.flat_step(p,(y-t)/T)
                    return p.exp(-p.mpf(k)*t+m0*coeff[0]*delta)*polynomial(delta)
                actual=p.quad(integrand,[0,y/4,y/2,3*y/4,y]);contains(rows[n],actual);compared+=1
            assert not proof['full_finite_lower_tail'].zero and not proof['full_finite_upper_tail'].zero
            zero,_=current.finite_kernel(f,B,c.mpf(T),c.mpf(0),kind)
            assert all(row.zero for row in zero)
    rejects=0
    bad=current.IntervalTaylor(c,[c.mpf(2)]+[c.mpf(0)]*5)
    try:current.finite_kernel(f,bad,c.mpf(1),c.mpf('.5'),'theta')
    except ValueError:rejects+=1
    else:raise AssertionError('Unadmitted finite kernel rate accepted')
    return dict(passed=True,independent_full_finite_original_kernel_Z_coefficients=compared,
        signed_nonzero_B_with_axial_rows_through5=True,finite_integral_not_terminal_limit=True,
        exact_zero_inlet_and_nonzero_finite_tails_verified=True,rate_guard_rejections=rejects)


def manufactured_source():
    c=MPIntervalContext();c.dps=170;f=flow(c);p=mp.mp.clone();p.dps=90
    z=sy.symbols('z',real=True);R=sy.Rational;T=20;y=10;rho=p.mpf('.5')
    Bexpr=R(1,50)+z*R(3,1000)-z*z/500+z**3/1000
    Vexpr=R(2,5)+z*R(3,100)+z*z/100
    initial=dict(theta=R(7,10)+z/50,theta_z=R(1,4)+z/100,pressure=R(5)+z/20,
        swirl=R(5,6)-z/30,mean=R(1,3)+z/40,axial=R(1,5)+z/50)
    P0expr=R(3,5)+z/10+z*z/30
    def jet(expr):return current.IntervalTaylor(c,[c.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),155))) for n in range(6)])
    B=jet(Bexpr);V=f.jet(jet(Vexpr));P0=f.jet(jet(P0expr));zjet=current.IntervalTaylor.variable(c,0,5)
    long=SimpleNamespace(flow=f,c=c,B=B,T=c.mpf(T),logC=c.ln(4),p0=jet(P0expr),V=V,
        inlets={name:f.jet(jet(expr)) for name,expr in initial.items()})
    reference=SimpleNamespace(flow=f,T=long.T,logC=long.logC,P0=P0,z=zjet,zrows=f.jet(zjet),delta=c.mpf('.01'))
    radius=f.factor((0,0,0,0,0),c.ln(110)+y)
    op,packet,recovered,a,record=current.source_packet(long,reference,radius,c.mpf('.5'),c.mpf(y))
    kernels={};sig=scalar.flat_step(p,rho);ds=p.diff(lambda value:scalar.flat_step(p,value),rho)
    for kind,(k,m) in dict(theta=('1.6',1),pressure=('.2',2),swirl=('1.2',2)).items():
        rows=[]
        for n in range(3):
            def integrand(t):
                delta=sig-scalar.flat_step(p,(y-t)/T);b0=p.mpf('0.02');b1=p.mpf('.003');b2=-p.mpf('.002')
                factor=1 if n==0 else m*delta*b1 if n==1 else m*delta*b2+(m*delta*b1)**2/2
                return p.exp(-p.mpf(k)*t+m*b0*delta)*factor
            rows.append(p.quad(integrand,[0,y/4,y/2,3*y/4,y]))
        kernels[kind]=sum(sy.Float(str(value),100)*z**n for n,value in enumerate(rows))
    sigsym=sy.Float(str(sig),100);dssym=sy.Float(str(ds),100)
    E=sy.exp(Bexpr*(1-sigsym)-sy.log(1+z*z)+R(y,10)-sy.log(4)-sy.log(3))
    shapes={name:initial[name]*sy.exp(-sy.Rational(k)*y+m*Bexpr*sigsym)+
        (Vexpr if name=='theta_z' else 1)*kernels[kind]
        for name,kind,m,k in (('theta','theta',1,'1.6'),('theta_z','theta',1,'1.6'),
                            ('pressure','pressure',2,'.2'),('swirl','swirl',2,'1.2'))}
    shapes['mean']=Vexpr+(initial['mean']-Vexpr)*sy.exp(-y)
    shapes['axial']=Vexpr**2+(initial['axial']-Vexpr**2)*sy.exp(-y)
    vv=Vexpr/3;hist=dict(m=shapes['mean']/3,h=E*shapes['theta'],k=E*shapes['theta_z']/3,
        e=shapes['axial']/9-E**2*shapes['swirl']/2,p=E**2*shapes['pressure']/2)
    de=R(1,100);L=1-de*z*z;d=1-z*z;pressure=P0expr+hist['p'];m,h,k,e,pp=(hist[name] for name in current.RATES)
    transport=m*z*(1-de)+sy.diff(m,z)*d
    sectors=dict(theta_linear=(-E+h*(1-de/2)-sy.diff(h,z)*z*(1-de)/2)/L,
        theta_quadratic=(k*z*(2*de-1)-sy.diff(k,z)*d+E*transport)/L,
        axial_linear=(-vv+(m-sy.diff(m,z)*z)*(1-de)/2)/L,
        axial_quadratic=(vv*transport+e*z*2*de-sy.diff(e,z)*d+pressure*z*2*(1+de)-sy.diff(pressure,z)*d)/L)
    expected_a=1-2*(R(1,10)-Bexpr*dssym/T)
    p2=110*sy.exp(y)*(sectors['axial_linear']+3*sectors['axial_quadratic'])/E
    scales=scalar.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1e12',dps=100)
    eta=sy.Float(scales.ctx.nstr(scales.eta,105),100)
    source,qr,qrecord=current.source_quotients(op,recovered,a,c.ln(c.mpf(scales.ctx.nstr(scales.eta,105))))
    qexpr=sy.sqrt((2*eta-(expected_a-2))/(2*expected_a));compared=0
    def reference_value(expr,n):return mp.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),95)))
    for name,expr in dict(E=E,V=vv).items():
        for n in range(3):contains(recovered['common_velocity_'+name+'_axial5'][n],reference_value(expr,n));compared+=1
    for name,expr in hist.items():
        for n in range(3):contains(recovered['common_own_five_histories_axial5'][name][n],reference_value(expr,n));compared+=1
    for name,expr in sectors.items():
        for n in range(2):contains(recovered['full_signed_inertial_sectors_axial4'][name][n],reference_value(expr,n));compared+=1
    for name,expr in dict(a=expected_a,E=E,p2=p2,q=qexpr).items():
        rows=qr if name=='q' else source['roots'][name]
        for n in (0,1):contains(rows[current.C0 if n==0 else current.Z],reference_value(expr,n));compared+=1
    # Independent original scalar inverse diagnostic at fixed phi, using the
    # genuine source a,p2,E derivatives. Auxiliary p1 is unused by A/B.
    dstar=c.ln(c.mpf(scales.ctx.nstr(scales.d_star,105)))
    got=current.primitives.all_u_primitive_bounds(f,source,qr,dstar,c.mpf([0,1]))
    density=current.phase.densities.density_Z_kernels(source['roots']['E'][current.C0],source['roots']['E'][current.Z],
        recovered['common_velocity_V_axial5'][0],recovered['common_velocity_V_axial5'][1],got['values'],257)
    ps=scales.ctx;functions=[sy.lambdify(z,expr,'mpmath') for expr in (expected_a,p2,E,vv)];hstep=ps.mpf('1e-7');density_checks=0
    for phi in ('.137','.554'):
        def scalar_density(zz):
            aa,pp2,ee,v0=(fun(zz) for fun in functions)
            loop=scalar.GenericShearLoop(scales,a=aa,b=0,p1=200,p2=pp2,Utheta=ee)
            row=loop.at_angle(loop.angle_at_phase(ps.mpf(phi)));inc=ee*ps.expm1(row['A']/257);dv=row['B']/257
            return dict(m=dv,h=inc,k=v0*inc+ee*dv+inc*dv,
                e=2*v0*dv+dv*dv-ee*inc-inc*inc/2,p=ee*inc+inc*inc/2)
        center=scalar_density(0);minus=scalar_density(-hstep);plus=scalar_density(hstep)
        for name in current.RATES:
            contains(density['kernels'][name],center[name],'1e-30')
            contains(density['Z_derivatives'][name],(plus[name]-minus[name])/(2*hstep),'1e-14');density_checks+=2
    return dict(passed=True,independent_full_kernel_normalized_source_and_inertial_coefficients=compared,
        independent_original_finite_N_density_C0_Z_comparisons=density_checks,
        exact_current_P0_preserved_and_Pstar_divided_once=True,nonzero_R110_leading_memory_and_pressure_terms_used=True,
        original_scalar_inverse_Z_difference_is_finite_diagnostic_only=True,
        scalar_inverse_auxiliary_p1_not_a_source_inertia_or_cone_claim=True)


def weight_references():
    c=MPIntervalContext();c.dps=150;f=flow(c);p=mp.mp.clone();p.dps=90;comparisons=0
    width=p.mpf(7)/13;suffix=p.mpf(2)/5;N=257;alpha=p.mpf('.2')
    for name,rate in current.RATES.items():
        mass,decay,downstream=current.kernel_weight(f,c.mpf(width),c.mpf(suffix),rate)
        lam=p.mpf(rate.numerator)/rate.denominator
        actual_mass=width if not rate else -p.expm1(-lam*width)/lam
        contains(mass,actual_mass);contains(decay,p.exp(-lam*width));contains(downstream,p.exp(-lam*suffix));comparisons+=3
        # Independent analytic integral of signed finite-frequency forcing.
        complex_mass=(p.exp(1j*N*width)-p.exp(-lam*width))/(lam+1j*N)
        exact=p.exp(-lam*suffix)*(complex_mass.imag/N+alpha*actual_mass/(N*N))
        cap=mass*downstream*c.mpf(1/N+alpha/(N*N))
        contains(current.bounds.symmetric(f,cap),exact);comparisons+=1
    tiny=p.mpf('1e-1000');tiny_checks=0
    for rate in current.RATES.values():
        mass,_,_=current.kernel_weight(f,c.mpf(tiny),c.mpf(0),rate)
        lam=p.mpf(rate.numerator)/rate.denominator
        actual=tiny if not rate else -p.expm1(-lam*tiny)/lam
        contains(mass,actual,'1e-1080',3000);assert ep(mass.finite_interval(max_log=3000))[0]>0;tiny_checks+=1
    B=current.IntervalTaylor.constant(c,0,5)
    for kind,k in dict(theta='1.6',pressure='.2',swirl='1.2').items():
        rows,_=current.finite_kernel(f,B,c.mpf(20),c.mpf(tiny),kind)
        contains(rows[0],-p.expm1(-p.mpf(k)*tiny)/p.mpf(k),'1e-1080',3000)
        assert ep(rows[0].finite_interval(max_log=3000))[0]>0;tiny_checks+=1
    rejects=0
    for width0,suffix0 in ((0,0),(-1,0),(1,-1)):
        try:current.kernel_weight(f,c.mpf(width0),c.mpf(suffix0),current.RATES['m'])
        except ValueError:rejects+=1
        else:raise AssertionError('Invalid Duhamel cell accepted')
    return dict(passed=True,independent_own_rate_mass_decay_and_signed_integral_comparisons=comparisons,
        positive_microscopic_width_mass_and_kernel_comparisons=tiny_checks,
        both_N_levels_and_nonzero_mean_forcing_retained=True,weight_guard_rejections=rejects)


def native(owner,report):
    cells=densities=source_rows=memory=guard_rejects=0
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(current.serialized(packet))==report['frames'][label]
        long=owner.long_wrapper.owner(label);reference=owner.reference_wrapper.owner(label);f,c=long.flow,long.c
        assert reference.flow is f and packet['exact_common_P0_axial5'] is reference.P0
        assert packet['actual_R110_correction_is_still_unsupplied'] and not packet['actual_Rm_incoming_correction_supplied']
        for cell in packet['actual_full_window_source_cells']:
            source=cell['source'];generic=source['original_generic_source'];record=source['original_background_function_record']
            assert source['exact_common_P0_axial5'] is reference.P0 and generic['common_original_P0_axial5'] is reference.P0
            assert record['exact_common_P0_axial5'] is reference.P0 and record['exact_V_y_and_V_yZ_zero']
            assert record['one_Pstar_conversion_for_raw_m_k_V'] and source['actual_radius_phase_Z_exact_zero']
            assert not source['genuine_finite_N_R110_boundary_correction_supplied']
            assert source['current_R110_source_binding']['exact_fixed_inlet_radius']==[110,1]
            assert source['current_R110_source_binding']['same_actual_anchored_B_axial5']
            assert generic['source_frame_conditional_on_same_current_R110_background']
            assert 'source_frame_conditional_on_same_actual_Rm_inlet' not in generic
            for proof in record['actual_full_finite_kernel_proofs'].values():
                assert proof['full_source_T_not_shortened'] and proof['source_kernel_enclosures_not_selected_values']
                assert not proof['full_finite_lower_tail'].zero and not proof['full_finite_upper_tail'].zero
            for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                for row in rows:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;source_rows+=1
            for part in ('kernels','Z_derivatives'):
                for row in source['original_signed_five_density_C0_Z'][part].values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;densities+=1
            for name,proof in cell['original_own_rate_weights'].items():
                assert proof['physical_log_radius_measure_once'] and proof['nonlinear_source_driver_not_zeroed']
                if name=='p':assert ep(proof['cell_incoming_decay'].coefficient)==(1,1)
            assert all(source[key] is False for key in fields.previous.OPEN);cells+=1
        for row in packet['actual_nonzero_incoming_memory_decays'].values():assert not row.zero;memory+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        original=owner.long_wrapper.saved['post_power_packets'][label][-1]
        for field,bad in (('source_family','wrong-source'),('source_frame','wrong-frame'),('radius',[111,1])):
            altered=dict(original)
            if field=='radius':
                altered['post_power_function_evaluation']=dict(original['post_power_function_evaluation'])
                altered['post_power_function_evaluation']['geometry']=dict(original['post_power_function_evaluation']['geometry'],exact_fixed_radius=bad)
            else:altered[field]=bad
            mock=SimpleNamespace(c=c,family=owner.family,
                long_wrapper=SimpleNamespace(saved=dict(post_power_packets={label:[altered]})))
            try:current.OriginalLongReshapeFiniteN.inlet_binding(mock,label,long)
            except ValueError:guard_rejects+=1
            else:raise AssertionError('Wrong current R110 source binding accepted')
        print('Current original full long-reshape finite-N source checked',label,flush=True)
    return dict(passed=True,actual_current_source_cells=cells,genuine_current_source_profile_and_history_coefficients=source_rows,
        original_finite_N_signed_density_C0_Z_rows=densities,nonzero_incoming_memory_rows=memory,
        current_R110_family_frame_radius_guard_rejections=guard_rejects,
        full_finite_T_and_current_P0_basis_ledger_verified=True,
        real_R110_Rm_inlet_whole_axis_cone_global_N_recursion_not_admitted=True)


def run():
    began=time.monotonic();theorem,polys=symbolic()
    with mp.workdps(240):
        kernels=kernel_references(polys);source=manufactured_source();weights=weight_references()
    print('Independent original long-reshape source, finite-N and Duhamel references PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalLongReshapeFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_source_and_weights=theorem,independent_full_finite_kernels=kernels,
        independent_original_common_unit_source_and_density=source,independent_signed_Duhamel=weights,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Current original long-reshape finite-N driver and local affine Duhamel PASS',flush=True);return result


if __name__=='__main__':run()
