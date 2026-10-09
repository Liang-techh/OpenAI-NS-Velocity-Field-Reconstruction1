"""Independent symbolic recovery, scalar spatial derivatives and actual jets."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_first_spatial_jets as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,ep=current.fields,current.ep


def flow(c):
    return MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def contains(row,value,allowance='1e-20'):
    lo,hi=ep(row.finite_interval()+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def independent_recovery():
    c=MPIntervalContext();c.dps=180;f=flow(c);y,Z=sy.symbols('y Z',real=True);R=sy.Rational
    shapes=dict(m=R(1,7)+y/11+Z/13+y*Z/17,
        h=R(1,9)+y/7+Z/11+y*y/13,
        k=R(1,11)+y/13+Z/7+y*Z/19,
        e=R(1,13)+y/17+Z/9+y*y/23,
        p=R(1,17)+y/19+Z/23+y*Z/29,
        E=R(13,10)+y/7+y*y/11+Z/13+y*Z/17,
        V=R(2,5)+y/11+y*y/13+Z/17+y*Z/19)
    def rows(expr):
        return [[f.scalar(str(sy.diff(expr,y,j,Z,n).subs({y:0,Z:0})/sy.factorial(n)))
                 for n in range(6)] for j in range(3)]
    P0expr=R(3,5)+Z/7+Z*Z/11
    P0=rows(P0expr)[0];zjet=IntervalTaylor.variable(c,0,5)
    op=SimpleNamespace(flow=f,c=c,zrows=f.jet(zjet),reference=SimpleNamespace(z=zjet,delta=c.mpf('1/100')),
        Pstar=f.factor((0,.5,0,0,0)),P0=P0,Rm_factor=f.scalar(5))
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={key:rows(expr*3 if key in ('m','k') else expr)
                   for key,expr in shapes.items() if key not in ('E','V')},
        velocity=dict(theta=rows(shapes['E']),axial=rows(shapes['V']*3))),
        original_P0_normalized_axial5=P0,geometry=dict(point=True,exact_x=[5,4]))
    dual=current.lift_generic_y(op,packet)
    m,h,k,e,p,E,V=(shapes[key] for key in ('m','h','k','e','p','E','V'))
    de=R(1,100);d=1-Z*Z;L=1-de*Z*Z;pressure=P0expr+p
    transport=m*Z*(1-de)+sy.diff(m,Z)*d
    expressions=dict(E=E,C=E-2*sy.diff(E,y),B=2*sy.diff(V,y),
        theta_linear=(-E+h*(1-de/2)-sy.diff(h,Z)*Z*(1-de)/2)/L,
        theta_quadratic=(k*Z*(2*de-1)-sy.diff(k,Z)*d+E*transport)/L,
        axial_linear=(-V+(m-Z*sy.diff(m,Z))*(1-de)/2)/L,
        axial_quadratic=(V*transport+2*de*Z*e-d*sy.diff(e,Z)+2*(1+de)*Z*pressure-d*sy.diff(pressure,Z))/L)
    comparisons=0
    for key,expr in expressions.items():
        exported=dual['actual_generic_source_numerators'][key] if key in ('E','C','B') else dual['full_signed_inertial_sectors_axial4'][key]
        for n,row in enumerate(exported[:2]):
            for order,value in ((0,row.v),(1,row.dy)):
                expected=sy.diff(expr,y,order,Z,n).subs({y:0,Z:0})/sy.factorial(n)
                contains(value,mp.mpf(str(sy.N(expected,160))),'1e-145');comparisons+=1
    I=expressions['axial_linear']+3*expressions['axial_quadratic'];p2bar=I/E
    sector=dual['full_signed_inertial_sectors_axial4']
    I0=sector['axial_linear'][0].v+sector['axial_quadratic'][0].v*op.Pstar
    Iy=sector['axial_linear'][0].dy+sector['axial_quadratic'][0].dy*op.Pstar
    E0=dual['common_velocity_E_axial5'][0].v;Ey=dual['common_velocity_E_axial5'][0].dy
    proof=current.phase.previous.positive_source(f,E0,'fixture E')
    bar=I0.positive_divide(E0,proof['source_log_lower'])
    bary=(Iy-bar*Ey).positive_divide(E0,proof['source_log_lower'])
    radius=op.Rm_factor*c.mpf('5/4');full_y=radius*(bary+bar)
    expected=sy.diff(R(25,4)*sy.exp(y)*p2bar,y).subs({y:0,Z:0})
    contains(full_y,mp.mpf(str(sy.N(expected,160))),'1e-145')
    assert sy.simplify(sy.diff(sy.exp(y)*p2bar,y)-sy.exp(y)*(sy.diff(p2bar,y)+p2bar))==0
    # Independent exact quotient and correlated shear differentiation.
    H,F,aa,bb,Delta=[sy.Function(name)(y) for name in ('H','F','a','b','Delta')]
    assert sy.simplify(sy.diff(R(4,5)-2*F/H,y)+2*(sy.diff(F,y)-F*sy.diff(H,y)/H)/H)==0
    eta=sy.symbols('eta',positive=True);q2=(2*eta-Delta)/(2*aa)
    assert sy.simplify(sy.diff(q2,y)-(-sy.diff(Delta,y)-2*q2*sy.diff(aa,y))/(2*aa))==0
    assert sy.simplify(sy.diff(aa+bb*bb/aa-2,y)-(sy.diff(aa,y)+(2*bb*sy.diff(bb,y)-bb*bb/aa*sy.diff(aa,y))/aa))==0
    x=sy.symbols('x',positive=True);g=sy.Function('g')(x)
    assert sy.simplify(x*sy.diff(x*sy.diff(g,x)-g/10,x)-x*(R(9,10)*sy.diff(g,x)+x*sy.diff(g,x,2)))==0
    return dict(passed=True,independent_symbolic_full_signed_recovery_value_y_Z_comparisons=comparisons,
        velocity_second_y_rows_and_history_y_rows_used=True,one_and_only_one_R_y_term_verified=True,
        exact_correlated_a_and_Delta_q_squared_derivatives_verified=True,
        fixtures_are_manufactured_sources_not_native_field_values=True)


def scalar_spatial():
    c=MPIntervalContext();c.dps=200;f=flow(c);y,Z=sy.symbols('y Z',real=True);R=sy.Rational
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1000',dps=150)
    p=scales.ctx;eta=R(p.nstr(scales.eta,155));dstar=c.ln(c.mpf(p.nstr(scales.d_star,155)))
    a=R(4,5)+y*R(11,100)+Z*R(7,100)
    b=R(1,5)+y*R(7,100)-Z/20;E=R(13,10)+y*R(3,100)+Z/25;V=R(2,5)+y/50+Z*R(3,100)
    Delta=a+b*b/a-2;q=sy.sqrt((2*eta-Delta)/(2*a))
    q0=sy.lambdify((y,Z),q,'mpmath')(0,0);N=257;phi=p.mpf('.137');h=p.mpf('1e-12');cases=[]
    for target_u in ('-2','0','.1','2'):
        p20=R(p.nstr(p.mpf(target_u)*scales.d_star/q0,155))
        p2=(p20+y*R(7,100)-Z/50)*sy.exp(y)
        def roots_of(expression):
            return {order:f.scalar(c.mpf(str(sy.N(sy.diff(expression,y,order[0],Z,order[1]).subs({y:0,Z:0}),170))))
                    for order in ((0,0),(1,0),(0,1))}
        roots={name:roots_of(expr) for name,expr in dict(a=a,t0=-b/a,E=E,p2=p2).items()};qr=roots_of(q)
        source=dict(q=qr[(0,0)],roots=roots)
        af,bf,pf,ef,vf=[sy.lambdify((y,Z),expr,'mpmath') for expr in (a,b,p2,E,V)]
        def reference(yy,zz,target,guess=None):
            loop=original.GenericShearLoop(scales,a=af(yy,zz),b=bf(yy,zz),p1=200,p2=pf(yy,zz),Utheta=ef(yy,zz))
            start=2*p.pi*target if guess is None else guess
            try:angle=p.findroot(lambda v:loop.phase_at_angle(v)-target,(start,start+p.mpf('.001')),tol=p.mpf('1e-135'))
            except (ValueError,ArithmeticError):angle=loop.angle_at_phase(target)
            return loop.at_angle(angle),angle
        zero,guess=reference(0,0,phi);comparisons=0
        got=current.phase.first.conditioned_first_jets(source,qr,dstar,c.mpf('.137'))
        assert got['values'] is not None
        bounded=current.all_u_first_bounds(f,source,qr,dstar,c.mpf('.137'))
        for result in (got,bounded):
            values=result['values'];assert set(values)==set(current.phase.first.OUTPUTS)
            for key,rkey in (('A','A'),('B_over_Pstar','B')):
                contains(values[key],zero[rkey]);comparisons+=1
            for direction in ('y','Z','phi'):
                samples=[reference(shift*h if direction=='y' else 0,shift*h if direction=='Z' else 0,
                    phi+(shift*h if direction=='phi' else 0),guess)[0] for shift in (-2,-1,1,2)]
                for primitive in ('A','B'):
                    derivative=(samples[0][primitive]-8*samples[1][primitive]+8*samples[2][primitive]-samples[3][primitive])/(12*h)
                    key=primitive+'_'+direction+('_over_Pstar' if primitive=='B' else '')
                    contains(values[key],derivative);comparisons+=1
            spatial=dict(values,A_y=values['A_y']+N*values['A_phi'],
                         B_y_over_Pstar=values['B_y_over_Pstar']+N*values['B_phi_over_Pstar'])
            def density_reference(yy,zz,total):
                primitive,_=reference(yy,zz,phi+(N*yy if total else 0),guess)
                ee,vv=ef(yy,zz),vf(yy,zz);dE=ee*p.expm1(primitive['A']/N);dV=primitive['B']/N
                pressure=ee*dE+dE*dE/2
                return dict(m=dV,h=dE,k=vv*dE+ee*dV+dE*dV,e=2*vv*dV+dV*dV-pressure,p=pressure)
            for direction,total,primitives in (('y',False,values),('y',True,spatial),('Z',True,values)):
                Ed=roots['E'][(1,0) if direction=='y' else (0,1)]
                Vd=roots_of(V)[(1,0) if direction=='y' else (0,1)]
                density=current.density_derivative(roots['E'][(0,0)],Ed,roots_of(V)[(0,0)],Vd,primitives,N,direction)
                samples=[density_reference(shift*h if direction=='y' else 0,shift*h if direction=='Z' else 0,total)
                         for shift in (-2,-1,1,2)]
                for key in current.previous.RATES:
                    derivative=(samples[0][key]-8*samples[1][key]+8*samples[2][key]-samples[3][key])/(12*h)
                    contains(density['derivatives'][key],derivative);comparisons+=1
            if target_u=='0':assert not values['A_y'].zero and not values['B_y_over_Pstar'].zero
        cases.append(dict(signed_u=target_u,comparisons=comparisons,passed=True))
    return dict(passed=True,independent_original_scalar_spatial_cases=cases,
        independent_comparisons=sum(row['comparisons'] for row in cases),
        both_original_conditioned_first_graph_and_all_u_genuine_y_Z_bounds_checked=True,
        fixed_phi_y_total_spatial_y_and_Z_five_density_derivatives_checked=True,
        actual_fast_phase_chain_independently_finite_differenced=True,
        finite_difference='five point h=1e-12, original scalar precision150, allowance1e-20; diagnostic only')


def native(owner,report):
    queries=cells=primitive_rows=density_rows=0;routes={}
    points=dict(active=((5,4),None),terminal=((2,1),None),Rh=('Rh',None),
                whole_active=((1,1),(71,40)),whole_terminal=((71,40),'Rh'))
    for label in ('0','.5'):
        op=owner.phase.upstream.upstream.upstream.owner(label).op;f=op.flow
        for key,(left,right) in points.items():
            packet=owner.query(label,left,right)
            assert current.base.encoded(fields.serialized(packet))==report['frames'][label][key]
            assert packet['genuine_first_y_Z_source_and_density_chain_installed']
            assert not packet['genuine_mixed_second_derivatives_or_sharp_integrals_installed']
            assert all(packet[k] is False for k in fields.previous.OPEN)
            radius=packet['actual_original_Rm_radius_phase'];assert radius['exact_source_Rm_factor'] is op.Rm_factor
            assert len(packet['actual_first_spatial_source_cells'])==len(radius['phase_boxes'])
            for cell in packet['actual_first_spatial_source_cells']:
                assert cell['source_family']==owner.family and cell['source_geometry']==radius['source_geometry']
                assert cell['exact_common_P0_axial5'] is op.P0
                record=cell['source_first_y_Z']
                assert record['original_common_P0_axial5'] is op.P0 and record['original_P0_y_exactly_zero']
                assert record['original_physical_radius_y'] is record['original_physical_radius']
                assert record['original_physical_radius'] is cell['exact_same_shared_radius_factor']
                assert record['genuine_y_from_mixed4_raw_rows_not_Z_surrogate'] and record['q_y_from_original_active_branch_with_fixed_eta']
                assert cell['actual_phase_y_exactly_N'] and cell['actual_phase_Z_exactly_zero']
                assert cell['genuine_first_y_Z_only_not_mixed_second_jets']
                values=cell['original_fixed_phi_first_primitive_values'];spatial=cell['actual_spatial_first_primitive_values']
                assert set(values)==set(current.phase.first.OUTPUTS)
                for name in ('A','B'):
                    ykey=name+'_y'+('_over_Pstar' if name=='B' else '')
                    phikey=name+'_phi'+('_over_Pstar' if name=='B' else '')
                    expected=values[ykey]+values[phikey]*packet['candidate_N']
                    assert current.base.encoded(fields.serialized(expected))==current.base.encoded(fields.serialized(spatial[ykey]))
                for row in values.values():
                    assert row.ctx is op.c and row.scale.bases is f.logs and row.ledger is f.ledger
                    assert all(mp.isfinite(v) for v in ep(row.coefficient));primitive_rows+=1
                for name in ('original_five_signed_density_C0','original_five_signed_density_fixed_phi_y',
                             'actual_five_signed_density_y','actual_five_signed_density_Z'):
                    assert set(cell[name])==set(current.previous.RATES)
                    for row in cell[name].values():
                        assert row.ctx is op.c and row.scale.bases is f.logs and row.ledger is f.ledger
                        assert all(mp.isfinite(v) for v in ep(row.coefficient));density_rows+=1
                assert all(cell[k] is False for k in fields.previous.OPEN)
                route=cell['original_primitive_route'];routes[route]=routes.get(route,0)+1;cells+=1
            queries+=1
        print('Actual Rm genuine first source/spatial density jets checked',label,flush=True)
    return dict(passed=True,actual_live_source_queries=queries,actual_phase_source_cells=cells,
        actual_original_C0_y_Z_phi_primitive_rows=primitive_rows,actual_five_C0_fixed_y_total_y_Z_density_rows=density_rows,
        original_function_routes=routes,same_P0_Rm_geometry_and_source_algebra_verified=True,
        whole_source_cells_have_first_derivative_enclosures_not_mixed_second_or_sharp_integral_closure=True)


def run():
    began=time.monotonic()
    with mp.workdps(260):recovery=independent_recovery();scalar=scalar_spatial()
    print('Independent genuine recovery and original spatial density derivatives PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmFirstSpatialJets(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_recovery_and_radius_derivatives=recovery,independent_original_scalar_spatial_derivatives=scalar,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original Rm genuine first y/Z and five density spatial chains PASS',flush=True);return result


if __name__=='__main__':run()
