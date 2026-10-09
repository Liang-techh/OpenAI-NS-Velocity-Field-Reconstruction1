"""Independent original scalar references and actual all-u edge replay."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_all_u_density_integrals as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def symbolic_integral_theorem():
    r,x=sy.symbols('r x',real=True);k=sy.symbols('k',integer=True,nonnegative=True)
    n=sy.symbols('n',integer=True,positive=True);angle=sy.symbols('angle',real=True)
    geometric=sy.summation(x**k,(k,0,sy.oo))
    assert isinstance(geometric,sy.Piecewise)
    assert sy.simplify(geometric.args[0][0]-1/(1-x))==0
    assert geometric.args[0][1]==(sy.Abs(x)<1)
    S0=1/(1-r*r)
    S1=sy.diff(S0,r)/2
    S2=sy.diff(r*S1,r)/(2*r)
    s=1-r*r
    assert sy.simplify(s*S0-1)==0
    norm_r=r*r/s*S0-2*r*S1+s*S2
    assert sy.simplify(norm_r-1/s**2)==0
    u=sy.symbols('u',real=True);ru=u/sy.sqrt(1+u*u)
    assert sy.simplify(1-ru*ru-1/(1+u*u))==0
    assert sy.simplify(sy.diff(ru,u)**2-(1+u*u)**-3)==0
    assert sy.simplify(norm_r*s**3-s)==0
    # The geometric Fourier real part is the original Poisson direction.
    numerator=sy.re((sy.cos(angle)+sy.I*sy.sin(angle))*(1-r*sy.cos(angle)+sy.I*r*sy.sin(angle)))
    assert sy.trigsimp(numerator-(sy.cos(angle)-r))==0
    assert sy.integrate(sy.cos(n*angle),(angle,0,2*sy.pi))==0
    assert sy.integrate(sy.cos(n*angle)**2,(angle,0,2*sy.pi))==sy.pi
    m=sy.symbols('m',integer=True,nonzero=True)
    assert sy.integrate(sy.cos(m*angle),(angle,0,2*sy.pi))==0
    aa,bb=sy.symbols('aa bb',real=True)
    assert sy.trigsimp(sy.cos(aa)*sy.cos(bb)-(sy.cos(aa-bb)+sy.cos(aa+bb))/2)==0
    t0,q=sy.symbols('t0 q',real=True)
    assert sy.simplify(2*sy.pi*t0*t0+4*sy.pi*q*q*s*S0-2*sy.pi*(t0*t0+2*q*q))==0
    assert sy.simplify(4*sy.pi*q*q*norm_r*s**3-4*sy.pi*q*q*s)==0
    M=sy.symbols('M',nonnegative=True);nu=1+M
    assert sy.simplify((nu-sy.sqrt(M))*(nu+sy.sqrt(M))-(M*M+M+1))==0
    v=sy.symbols('v',nonnegative=True)
    assert sy.simplify(sy.Rational(1,2)-v/(1+v*v)-(v-1)**2/(2*(1+v*v)))==0
    z=sy.symbols('Z',real=True)
    E,t,A,a,J=[sy.Function(name)(z) for name in ('E','t0','A','a','J')]
    product=sy.diff(E,z)*(t*A-a*J)+E*(sy.diff(t,z)*A+t*sy.diff(A,z)-sy.diff(a,z)*J-a*sy.diff(J,z))
    assert sy.simplify(sy.diff(E*(t*A-a*J),z)-product)==0
    return dict(passed=True,geometric_series_sum_verified_on_abs_r_less_than_one=True,
        original_Poisson_direction_real_part_verified=True,
        full_period_cosine_mean_norm_and_cross_orthogonality_verified=True,
        normalized_Fourier_coefficient_norm_and_u_derivative_norm_verified=True,
        all_finite_u_has_s_positive_and_abs_r_less_than_one=True,
        original_t_tq_tu_L2_constants_verified=True,
        Cauchy_implicit_Z_majorant_uses_sqrt_M_leq_nu=True,
        endpoint_J_Z_factor_from_abs_t_over_one_plus_t_squared_leq_half=True,
        complete_B_Z_product_rule_verified=True,
        identities_are_source_theorems_not_finite_native_value_fits=True)


def scalar_fixtures():
    c=MPIntervalContext();c.dps=200;z=sy.Symbol('Z',real=True);R=sy.Rational
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1000',dps=150)
    p=scales.ctx;eta=R(p.nstr(scales.eta,155));dstar=c.ln(c.mpf(p.nstr(scales.d_star,155)))
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    a=R(4,5)+R(7,100)*z;b=R(1,5)-z/20;E=R(13,10)+z/25
    Delta=a+b*b/a-2;q=sy.sqrt((2*eta-Delta)/(2*a))
    q0=sy.lambdify(z,q,'mpmath')(mp.mpf(0))
    tolerance=mp.mpf('1e-140');allowance=c.mpf(['-1e-24','1e-24']);h=p.mpf('1e-7');evidence=[]
    for target_u in ('-100','-2','-.1','0','.1','2','100'):
        p20=R(p.nstr(p.mpf(target_u)*scales.d_star/q0,155));p2=p20+R(7,100)*z
        def pair(expression):
            return {(0,n):f.scalar(c.mpf([value-tolerance,value+tolerance]))
                    for n in (0,1) for value in [sy.lambdify(z,sy.diff(expression,z,n),'mpmath')(mp.mpf(0))]}
        roots={name:pair(expression) for name,expression in dict(a=a,t0=-b/a,E=E,p2=p2).items()};qr=pair(q)
        source=dict(roots=roots,q=qr[(0,0)])
        af,bf,pf,ef=[sy.lambdify(z,expression,'mpmath') for expression in (a,b,p2,E)]
        def reference(zz,target,guess=None):
            loop=original.GenericShearLoop(scales,a=af(zz),b=bf(zz),p1=200,p2=pf(zz),Utheta=ef(zz))
            start=2*p.pi*target if guess is None else guess
            try:angle=p.findroot(lambda value:loop.phase_at_angle(value)-target,(start,start+p.mpf('.001')),tol=p.mpf('1e-135'))
            except (ValueError,ArithmeticError):angle=loop.angle_at_phase(target)
            return loop.at_angle(angle),angle
        comparisons=0
        for phi in (p.mpf('.137'),p.mpf('.499'),p.mpf('.863')):
            got=current.all_u_primitive_bounds(f,source,qr,dstar,c.mpf(p.nstr(phi,155)))
            assert set(got['values'])==set(current.phase.OUTPUTS)
            assert got['record']['original_genuine_Z_rows_retained']
            ref,guess=reference(p.mpf(0),phi)
            for key,rkey in (('A','A'),('B_over_Pstar','B')):
                lo,hi=ep(got['values'][key].finite_interval()+allowance)
                assert lo<=ref[rkey]<=hi,(target_u,phi,key);comparisons+=1
            for direction in ('Z','phi'):
                samples=[reference(shift*h if direction=='Z' else p.mpf(0),
                                   phi+(shift*h if direction=='phi' else 0),guess)[0] for shift in (-2,-1,1,2)]
                for primitive in ('A','B'):
                    value=(samples[0][primitive]-8*samples[1][primitive]+8*samples[2][primitive]-samples[3][primitive])/(12*h)
                    key=primitive+'_'+direction+('_over_Pstar' if primitive=='B' else '')
                    lo,hi=ep(got['values'][key].finite_interval()+allowance)
                    assert lo<=value<=hi,(target_u,phi,key);comparisons+=1
            N=257;density=current.phase.densities.density_Z_kernels(roots['E'][(0,0)],roots['E'][(0,1)],
                f.scalar('.4'),f.scalar('.03'),got['values'],N)
            def density_reference(zz):
                primitive,_=reference(zz,phi,guess);ee=ef(zz);vv=p.mpf('.4')+p.mpf('.03')*zz
                dE=ee*p.expm1(primitive['A']/N);dV=primitive['B']/N;cross=ee*dE+dE*dE/2
                return dict(m=dV,h=dE,k=vv*dE+ee*dV+dE*dV,e=2*vv*dV+dV*dV-cross,p=cross)
            zeroth=density_reference(p.mpf(0));samples=[density_reference(shift*h) for shift in (-2,-1,1,2)]
            for key in current.RATES:
                derivative=(samples[0][key]-8*samples[1][key]+8*samples[2][key]-samples[3][key])/(12*h)
                for row,value in ((density['kernels'][key],zeroth[key]),(density['Z_derivatives'][key],derivative)):
                    lo,hi=ep(row.finite_interval()+allowance)
                    assert lo<=value<=hi,(target_u,phi,key);comparisons+=1
            if target_u=='0':
                assert got['record']['original_u']['sign']=='undetermined'
                assert got['record']['original_u_Z']['sign']!='zero'
                assert not got['values']['A_Z'].zero and not got['values']['B_Z_over_Pstar'].zero
        evidence.append(dict(passed=True,target_signed_u=target_u,
            independent_original_scalar_primitive_Z_phi_and_signed_density_comparisons=comparisons,
            three_nonvertex_phase_values_checked=True,
            finite_difference='five point, h=1e-7, precision150, allowance1e-24; diagnostic only'))
    rejects=0
    for phi in (c.mpf('-.1'),c.mpf('1.1')):
        try:current.all_u_primitive_bounds(f,source,qr,dstar,phi)
        except ValueError:rejects+=1
        else:raise AssertionError('Outside closed original phase interval accepted')
    bad=dict(q=qr[(0,0)].scalar(1),roots=roots)
    try:current.all_u_primitive_bounds(f,bad,qr,dstar,c.mpf('.2'))
    except ValueError:rejects+=1
    else:raise AssertionError('Different q source object accepted')
    return dict(passed=True,independent_original_cases=evidence,
        total_independent_original_scalar_comparisons=sum(row['independent_original_scalar_primitive_Z_phi_and_signed_density_comparisons'] for row in evidence),
        phase_and_exact_q_source_guard_rejections=rejects,
        signed_small_large_and_zero_u_with_nonzero_u_Z_checked=True,
        finite_fixture_values_not_native_function_definitions=True)


def native(owner,report):
    closed=phase_boxes=primitives=densities=active_cells=cell_rows=whole_rows=0;signs={}
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(fields.serialized(packet))==report['frames'][label]
        op=owner.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        assert packet['exact_common_P0_axial5'] is op.P0
        partition=packet['exact_actual_active_partition'];current.previous.validate_active_partition(c,partition)
        cells=packet['actual_complete_active_cells_to_join']
        assert len(partition)==len(cells)+1
        for left,right,cell in zip(partition,partition[1:],cells):
            assert cell['source_geometry']==current.previous.previous.exact_geometry(left,right)
            assert cell['source_family']==owner.family and cell['target_coordinate']==current.previous.JOIN
            assert cell['actual_phase_source']['exact_source_Rm_factor'] is op.Rm_factor
            assert cell['actual_phase_source']['actual_Rm_radius_phase_Z_independent']
            assert ep(cell['source_radial_logarithmic_cell_width'])[0]>0
            assert ep(cell['source_logarithmic_suffix_to_target'])[0]>=0
            for key,rate in current.RATES.items():
                assert ep(cell['original_positive_own_rate_masses'][key])[0]>0
                if not rate:assert ep(cell['original_own_rate_suffix_decays'][key])==(1,1)
                for row in cell['actual_cell_to_target_integral_C0_Z'][key]:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger
                    assert all(mp.isfinite(v) for v in ep(row.coefficient));cell_rows+=1
            active_cells+=1
        assert len(packet['actual_all_u_closed_edge_sources'])==(4 if label=='0' else 0)
        for replacement in packet['actual_all_u_closed_edge_sources']:
            source=replacement['actual_all_u_closed_density_source'];radius=source['actual_original_Rm_radius_phase']
            assert source['actual_original_spatial_Z_density_interface_installed'] and source['actual_closed_radial_source_cell']
            assert radius['exact_source_Rm_factor'] is op.Rm_factor and radius['candidate_N']==257
            assert radius['source_geometry']==replacement['original_unresolved_source_geometry']
            assert len(source['actual_source_bound_phase_density_cells'])==len(radius['phase_boxes'])
            for cell in source['actual_source_bound_phase_density_cells']:
                record=cell['original_phase_Z_only_result'];values=cell['actual_original_primitive_Z_values']
                assert record['status']=='enclosed_by_original_all_u_integral_theorem'
                assert record['original_inverse_identity_and_Poisson_L2_theorem']==current.THEOREM
                assert record['original_u']['sign']=='undetermined' and record['original_u_Z']['sign']!='zero'
                assert record['original_genuine_Z_rows_retained'] and record['numerical_inverse_bracket_or_selected_field_value_not_claimed']
                assert cell['original_common_P0_axial5'] is op.P0
                assert set(values)==set(current.phase.OUTPUTS) and not any('_y' in key for key in values)
                for row in values.values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;primitives+=1
                density=cell['actual_original_five_signed_density_C0_Z'];assert density is not None
                for name in ('kernels','Z_derivatives'):
                    assert set(density[name])==set(current.RATES)
                    for row in density[name].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;densities+=1
                phase_boxes+=1
            closed+=1
        active=packet['actual_active_local_defect_integral_C0_Z']
        whole=packet['actual_whole_patch_local_defect_integral_C0_Z']
        expected=current.previous.previous.affine_transport(f,active,packet['actual_terminal_local_defect_integral_C0_Z'],
                                                           packet['actual_terminal_logarithmic_interval_width'])
        assert current.base.encoded(fields.serialized(expected))==current.base.encoded(fields.serialized(whole))
        for key in current.RATES:
            for n in (0,1):
                total=sum((cell['actual_cell_to_target_integral_C0_Z'][key][n] for cell in cells),f.scalar(0))
                assert current.base.encoded(fields.serialized(total))==current.base.encoded(fields.serialized(active[key][n]))
                whole_rows+=1
        assert packet['unknown_local_source_integral_count']==0 and packet['complete_local_whole_patch_C0_Z_integral_bounds_installed']
        assert packet['finite_N_Rm_incoming_defect_is_unsupplied_affine_argument']
        assert packet['integral_bound_completeness_not_sharp_defect_or_five_moment_closure']
        assert all(packet[key] is False for key in fields.previous.OPEN)
        try:owner.transport_supplied_incoming(label,None)
        except ValueError:pass
        else:raise AssertionError('Actual Rm incoming correction silently reset')
        signs[label]={key:[row.record()['sign'] for row in rows] for key,rows in whole.items()}
    return dict(passed=True,actual_previous_unknown_edge_cells_closed=closed,
        actual_all_u_source_phase_boxes=phase_boxes,original_all_u_C0_Z_phi_primitive_rows=primitives,
        original_all_u_five_density_C0_Z_rows=densities,complete_actual_active_source_cells=active_cells,
        complete_actual_active_cell_to_join_C0_Z_rows=cell_rows,complete_whole_patch_C0_Z_integral_rows=whole_rows,
        conditional_complete_local_integral_signs=signs,
        missing_actual_prefix_sharp_moment_closure_whole_Z_repair_cone_global_N_and_recursion_not_admitted=True)


def run():
    began=time.monotonic()
    theorem=symbolic_integral_theorem()
    with mp.workdps(260):fixture=scalar_fixtures()
    print('Independent original all-u C0/Z/phi and five signed density references PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmAllUDensityIntegrals(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family and report['original_all_u_theorem']==current.base.encoded(current.THEOREM)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        original_all_u_symbolic_source_theorem=current.THEOREM,independent_original_scalar_fixtures=fixture,
        independent_symbolic_Poisson_integral_and_derivative_proof=theorem,
        actual_live_source=evidence,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original all-u bounds and complete actual Rm local C0/Z integrals PASS',flush=True);return result


if __name__=='__main__':run()
