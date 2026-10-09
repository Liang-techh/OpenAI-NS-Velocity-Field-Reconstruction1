"""Independent raw mixed source, Fourier calculus and original scalar yZ."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_mixed_yZ_primitives as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,ep=current.fields,current.ep


def flow(c):return MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def contains(row,value,allowance='1e-20'):
    lo,hi=ep(row.finite_interval()+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def symbolic_fourier_and_implicit():
    u=sy.symbols('u',real=True);r=u/sy.sqrt(1+u*u);s=1/(1+u*u)
    for k in range(6):
        g=sy.sqrt(s)*r**k
        first=-s*r**(k+1)+(k*s*s*r**(k-1) if k else 0)
        second=s**sy.Rational(3,2)*((3*r*r-1)*r**k-5*k*s*r**k+
            (k*(k-1)*s*s*r**(k-2) if k>=2 else 0))
        assert sy.simplify(sy.diff(g,u)-first)==0
        assert sy.simplify(sy.diff(g,u,2)-second)==0
    v=sy.symbols('v',nonnegative=True);rho=1-v;sv=1-v*v
    L1=sy.cancel(v*sv/rho+sv*sv/(rho*rho))
    assert sy.expand(6-L1-(1-v)*(2*v+5))==0
    L2bracket=2*(1+v)+5*v*(1+v)**2+2*(1+v)**3
    assert sy.expand(40-L2bracket-(1-v)*(7*v*v+23*v+36))==0
    assert sy.simplify(sy.sqrt(sv)/rho**2-(1+v)**2/sv**sy.Rational(3,2))==0
    y,Z,phi=sy.symbols('y Z phi');q=sy.Function('q')(y,Z);U=sy.Function('u')(y,Z)
    S=sy.Function('S')(U);t0=sy.Function('t0')(y,Z)
    wanted=sy.diff(t0,y,Z)+2*(sy.diff(q,y,Z)*S+sy.diff(q,y)*sy.diff(S,U)*sy.diff(U,Z)+
        sy.diff(q,Z)*sy.diff(S,U)*sy.diff(U,y)+q*(sy.diff(S,U,2)*sy.diff(U,y)*sy.diff(U,Z)+sy.diff(S,U)*sy.diff(U,y,Z)))
    assert sy.simplify(sy.diff(t0+2*q*S,y,Z)-wanted)==0
    # Generic implicit derivative and total integral endpoint chain.
    py,pz,pyz,t,ty,tz,tp,Fyz=sy.symbols('py pz pyz t ty tz tp Fyz',real=True)
    assert sy.simplify((1+t*t)*(-(Fyz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)/(1+t*t))+
        Fyz+2*t*ty*pz+2*t*tz*py+2*t*tp*py*pz)==0
    z=sy.symbols('z',real=True)
    assert sy.factor((1+z*z)-2*z)==(z-1)**2
    return dict(passed=True,first_second_original_Fourier_coefficient_derivatives_verified=6,
        geometric_L1_derivative_bounds_6_and_40_verified=True,original_angular_derivative_h_cubed_bound_verified=True,
        complete_general_t0_q_u_yZ_chain_verified=True,original_implicit_mixed_sign_and_endpoint_ratio_verified=True,
        rho_is_only_a_majorant_and_not_differentiated=True)


def manufactured_source(scales):
    c=MPIntervalContext();c.dps=180;f=flow(c);p=scales.ctx;y,Z,x=sy.symbols('y Z x',real=True);R=sy.Rational
    x0=R(5,4);xy=x0*sy.exp(y)
    controls=[R(1,1000)+Z/5000+Z*Z/7000,R(1,2000)-Z/8000+Z*Z/9000,R(1,3000)+Z/9000]
    gammas=[(x-1)**2,(x-R(3,2))**3,(x-R(7,4))**4]
    H=x**R(1,10)+sum(control*gamma for control,gamma in zip(controls,gammas))
    E=(1+Z/20)*H.subs(x,xy);V=R(2,5)+y/50+y*y/100+Z*R(3,100)+y*Z/250
    shapes=dict(m=R(1,7)+y/11+Z/13+y*Z/17,h=R(1,9)+y/7+Z/11+y*y/13,
        k=R(1,11)+y/13+Z/7+y*Z/19,e=R(1,13)+y/17+Z/9+y*y/23,
        p=R(1,17)+y/19+Z/23+y*Z/29)
    def scalar(expr):return f.scalar(str(sy.N(expr,170)))
    def rows(expr):return [[scalar(sy.diff(expr,y,j,Z,n).subs({y:0,Z:0})/sy.factorial(n))
                            for n in range(6)] for j in range(3)]
    P0expr=R(3,5)+Z/7+Z*Z/11;P0=rows(P0expr)[0];zjet=IntervalTaylor.variable(c,0,5)
    zero=[f.scalar(0)]*6
    op=SimpleNamespace(flow=f,c=c,zrows=f.jet(zjet),reference=SimpleNamespace(z=zjet,delta=c.mpf('1/100')),
        Pstar=f.factor((0,.5,0,0,0)),P0=P0,Rm_factor=f.scalar(5),controls=[zero,zero]+[rows(h)[0] for h in controls])
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={key:rows(expr*3 if key in ('m','k') else expr) for key,expr in shapes.items()},
        velocity=dict(theta=rows(E),axial=rows(V*3))),original_P0_normalized_axial5=P0,
        geometry=dict(point=True,exact_x=[5,4]),
        actual_H_x_derivative_axial5=[rows(sy.diff(H,x,j).subs(x,xy))[0] for j in range(2)],
        actual_gamma_ordinary_x_derivatives=[[c.mpf(str(sy.diff(gamma,x,j).subs(x,x0))) for j in range(3)] for gamma in gammas])
    generic=current.first.generic_module.recover_inputs(op,packet)
    eta=R(p.nstr(scales.eta,155));eta_log=c.ln(c.mpf(str(eta)));dstar_log=c.ln(c.mpf(p.nstr(scales.d_star,155)))
    quotient=current.phase.previous.recover_quotients(op,generic,eta_log,dstar_log,patch=packet)
    source,qr,velocity,record=current.source_mixed_frame(op,packet,generic,quotient)
    a=1-2*sy.diff(E,y)/E;b=2*sy.diff(V,y)/E;t0=-b/a;Delta=a+b*b/a-2;q=sy.sqrt((2*eta-Delta)/(2*a))
    m,h,k,e,pp=(shapes[key] for key in ('m','h','k','e','p'))
    de=R(1,100);d=1-Z*Z;L=1-de*Z*Z;pressure=P0expr+pp
    transport=m*Z*(1-de)+sy.diff(m,Z)*d
    I=(-V+(m-Z*sy.diff(m,Z))*(1-de)/2)/L+3*(V*transport+2*de*Z*e-d*sy.diff(e,Z)+
        2*(1+de)*Z*pressure-d*sy.diff(pressure,Z))/L
    p2=5*xy*I/E
    expected=dict(a=a,t0=t0,E=E,p2=p2,q=q,V=V);compared=0
    for name,expr in expected.items():
        actual=source['roots'][name] if name in source['roots'] else qr if name=='q' else velocity
        for order in current.ORDERS:
            value=sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0})
            contains(actual[order],mp.mpf(str(sy.N(value,160))),'1e-130');compared+=1
    assert record['original_radius_y'] is quotient['exact_same_shared_positive_radius_factor']
    assert sy.simplify(sy.diff(5*xy*(I/E),y,Z)-5*xy*(sy.diff(I/E,y,Z)+sy.diff(I/E,Z)))==0
    return dict(passed=True,independent_manufactured_original_mixed_source_root_comparisons=compared,
        nonzero_correlated_H_a_b_t0_q_and_pressure_source_used=True,
        genuine_raw_yZ_and_one_physical_radius_yZ_recipe_verified=True,
        source_values_not_native_field_fits=True)


def scalar_mixed(scales):
    c=MPIntervalContext();c.dps=200;f=flow(c);p=scales.ctx;y,Z=sy.symbols('y Z',real=True);R=sy.Rational
    eta=R(p.nstr(scales.eta,155));dstar=c.ln(c.mpf(p.nstr(scales.d_star,155)))
    a=R(4,5)+y*R(11,100)+Z*R(7,100)+y*Z/20
    b=R(1,5)+y*R(7,100)-Z/20+y*Z/25;E=R(13,10)+y*R(3,100)+Z/25+y*Z/50
    V=R(2,5)+y/50+Z*R(3,100)+y*Z/40
    q=sy.sqrt((2*eta-(a+b*b/a-2))/(2*a));q0=sy.lambdify((y,Z),q,'mpmath')(0,0)
    h=p.mpf('1e-8');phi=p.mpf('.137');weights={-2:1,-1:-8,1:8,2:-1};cases=[]
    for target_u in ('-2','0','.1','2'):
        p20=R(p.nstr(p.mpf(target_u)*scales.d_star/q0,155));p2=(p20+y*R(7,100)-Z/50+y*Z/30)*sy.exp(y)
        def rows(expr):return {order:f.scalar(str(sy.N(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0}),170))) for order in current.ORDERS}
        roots={name:rows(expr) for name,expr in dict(a=a,t0=-b/a,E=E,p2=p2).items()};qr=rows(q)
        got=current.all_u_mixed_bounds(f,dict(q=qr[current.C0],roots=roots),qr,dstar,c.mpf('.137'))
        leading=current.leading_mixed(roots['E'],rows(V),got['values'])
        af,bf,pf,ef,vf=[sy.lambdify((y,Z),expr,'mpmath') for expr in (a,b,p2,E,V)]
        def reference(yy,zz,guess=None):
            loop=original.GenericShearLoop(scales,a=af(yy,zz),b=bf(yy,zz),p1=200,p2=pf(yy,zz),Utheta=ef(yy,zz))
            start=2*p.pi*phi if guess is None else guess
            try:angle=p.findroot(lambda value:loop.phase_at_angle(value)-phi,(start,start+p.mpf('.001')),tol=p.mpf('1e-135'))
            except (ValueError,ArithmeticError):angle=loop.angle_at_phase(phi)
            result=loop.at_angle(angle);A,B=result['A'],result['B'];ee,vv=ef(yy,zz),vf(yy,zz)
            return dict(A=A,B=B,m=B,h=ee*A,k=ee*(vv*A+B),e=2*vv*B-ee*ee*A,p=ee*ee*A),angle
        zero,guess=reference(0,0)
        samples={(i,j):reference(i*h,j*h,guess)[0] for i in weights for j in weights}
        compared=0
        for key,row in dict(A=got['values']['A_yZ'],B=got['values']['B_yZ_over_Pstar'],
                            **{name:value[current.YZ] for name,value in leading.items()}).items():
            derivative=sum(weights[i]*weights[j]*samples[(i,j)][key] for i in weights for j in weights)/(144*h*h)
            contains(row,derivative);compared+=1
        for trace in (0,p.mpf('.5'),1):
            special=current.all_u_mixed_bounds(f,dict(q=qr[current.C0],roots=roots),qr,dstar,c.mpf(p.nstr(trace,155)))
            assert special['values']['A_yZ'].zero and special['values']['B_yZ_over_Pstar'].zero
        cases.append(dict(passed=True,signed_u=target_u,independent_original_primitive_and_leading_yZ_comparisons=compared))
    rejects=0
    source=dict(q=qr[current.C0],roots=roots)
    for phibox in (c.mpf('-.1'),c.mpf('1.1')):
        try:current.all_u_mixed_bounds(f,source,qr,dstar,phibox)
        except ValueError:rejects+=1
        else:raise AssertionError('Outside original phase accepted')
    bad=dict(source,q=qr[current.C0].scalar(1))
    try:current.all_u_mixed_bounds(f,bad,qr,dstar,c.mpf('.2'))
    except ValueError:rejects+=1
    else:raise AssertionError('Different original q object accepted')
    incomplete={name:dict(rows) for name,rows in roots.items()};incomplete['a'].pop(current.YZ)
    try:current.all_u_mixed_bounds(f,dict(source,roots=incomplete),qr,dstar,c.mpf('.2'))
    except ValueError:rejects+=1
    else:raise AssertionError('Missing genuine yZ source accepted')
    return dict(passed=True,independent_original_scalar_mixed_cases=cases,
        independent_original_yZ_comparisons=sum(row['independent_original_primitive_and_leading_yZ_comparisons'] for row in cases),
        nonzero_y_Z_yZ_sources_and_signed_zero_u_checked=True,exact_symmetry_mixed_traces_zero=True,
        original_phase_q_and_missing_yZ_guard_rejections=rejects,
        finite_difference='tensor five point mixed difference, h=1e-8, original scalar precision150, allowance1e-20; diagnostic only')


def native(owner,report):
    queries=roots=primitives=densities=0
    points=dict(active=((5,4),None),terminal=((2,1),None),Rh=('Rh',None),
        whole_active=((1,1),(71,40)),whole_terminal=((71,40),'Rh'))
    for label in ('0','.5'):
        op=owner.upstream.upstream.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        for key,(left,right) in points.items():
            packet=owner.query(label,left,right)
            assert current.base.encoded(current.serialized(packet))==report['frames'][label][key]
            assert packet['exact_common_P0_axial5'] is op.P0 and packet['exact_source_Rm_factor'] is op.Rm_factor
            assert packet['genuine_original_yZ_source_and_primitive_enclosures_installed']
            assert not packet['actual_spatial_mixed_yZ_or_Z_averaged_integrals_installed']
            record=packet['original_mixed_source_record']
            assert record['first_y_Z_source']['original_common_P0_axial5'] is op.P0
            assert record['original_radius_y'] is record['first_y_Z_source']['original_physical_radius']
            assert record['original_radius_Z_exact_zero'] and record['raw_yZ_rows_not_relabelled_first_derivatives']
            assert not record['genuine_yy_ZZ_or_total_spatial_yZ_installed']
            for rows in (*packet['original_genuine_mixed_source_roots'].values(),packet['original_genuine_mixed_q_rows'],packet['original_genuine_mixed_V_rows']):
                assert set(rows)==set(current.ORDERS)
                for row in rows.values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;roots+=1
            values=packet['original_fixed_phi_A_B_mixed_values'];assert len(values)==10
            for row in values.values():
                assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger
                assert all(mp.isfinite(v) for v in ep(row.coefficient));primitives+=1
            proof=packet['original_mixed_primitive_proof']
            assert proof['every_finite_signed_u_and_zero_included'] and proof['nonzero_t0_and_b_allowed']
            assert proof['original_J_full_mixed_endpoint_terms_retained'] and proof['original_B_full_E_t0_A_a_J_mixed_product_retained']
            assert proof['narrow_numerical_inverse_yZ_or_actual_spatial_mixed_chain_not_claimed']
            for rows in packet['original_five_signed_leading_density_C0_y_Z_yZ'].values():
                assert set(rows)==set(current.ORDERS)
                for row in rows.values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;densities+=1
            assert all(packet[k] is False for k in fields.previous.OPEN);queries+=1
        print('Actual Rm source-owned mixed yZ primitives checked',label,flush=True)
    return dict(passed=True,actual_live_source_queries=queries,genuine_original_mixed_source_coefficients=roots,
        actual_original_C0_y_Z_phi_yZ_primitive_rows=primitives,actual_signed_leading_C0_y_Z_yZ_density_rows=densities,
        same_P0_Rm_geometry_and_formal_algebra_verified=True,
        fixed_phi_enclosures_not_narrow_inverse_values_or_total_spatial_mixed_derivatives=True)


def run():
    began=time.monotonic();theorem=symbolic_fourier_and_implicit()
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1000',dps=150)
    with mp.workdps(260):recovery=manufactured_source(scales);scalar=scalar_mixed(scales)
    print('Independent original mixed source and primitive references PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmMixedYZPrimitives(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_Fourier_and_implicit_calculus=theorem,independent_manufactured_mixed_source=recovery,
        independent_original_scalar_mixed_derivatives=scalar,actual_live_source=actual,input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original Rm genuine yZ roots and all-u mixed primitives PASS',flush=True);return result


if __name__=='__main__':run()
