"""Independent exact split / frozen original quadrature and actual Rm IBP."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_weighted_averaging as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def symbolic():
    E,V,A,B,RE,F,N=sy.symbols('E V A B RE F N',nonzero=True)
    dE=E*A/N+RE/N**2;dV=B/N
    exact=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    leading=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)
    remainder=dict(m=0,h=RE,k=V*RE+F*B,e=-E*RE+B*B-F*F/2,p=E*RE+F*F/2)
    for key in current.RATES:
        assert sy.expand(exact[key]-leading[key]/N-sy.sympify(remainder[key]).subs(F,E*A+RE/N)/N**2)==0
    y,phi=sy.symbols('y phi');lam=sy.symbols('lam',nonnegative=True)
    G=sy.Function('G')(y,phi);K=sy.Function('K')(y)
    total=lam*K*G+K*(sy.diff(G,y)+N*sy.diff(G,phi))
    assert sy.expand(K*sy.diff(G,phi)/N-(total-K*(sy.diff(G,y)+lam*G))/N**2)==0
    psi,a,b,t0,T1=sy.symbols('psi a b t0 T1',nonzero=True)
    aa=a/2*(phi-psi/(2*sy.pi));bb=E/2*(-a*T1/(2*sy.pi)-b*phi)
    assert sy.simplify(aa.xreplace({phi:1-phi,psi:2*sy.pi-psi})+aa)==0
    assert sy.simplify((bb.xreplace({phi:1-phi,T1:2*sy.pi*t0-T1})+bb).subs(t0,-b/a))==0
    x,t=sy.symbols('x t',real=True)
    assert sy.simplify(sy.diff(sy.exp(x),x,2)-sy.exp(x))==0
    assert sy.integrate(1-t,(t,0,1))==sy.Rational(1,2)
    return dict(passed=True,all_five_exact_signed_leading_and_second_remainder_splits_verified=True,
        original_general_nonzero_t0_reflection_verified=True,weighted_IBP_sign_and_Nminus2_verified=True,
        second_exponential_remainder_integral_mass_half_verified=True,
        full_nonlinear_density_mean_not_assumed_zero=True)


def frozen_original_quadrature():
    c=MPIntervalContext();c.dps=160
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1000',dps=90)
    p=scales.ctx;N=257;E=p.mpf('1.3');V=p.mpf('.4');cases=[];zero=f.scalar(0)
    for sign in (-1,1):
        loop=original.GenericShearLoop(scales,a='.8',b='.2',p1=200,p2=sign*p.mpf('.2'),Utheta=E)
        constant=lambda value:{(0,0):f.scalar(p.nstr(value,95)),(1,0):zero,(0,1):zero}
        roots={name:constant(value) for name,value in dict(a=loop.a,t0=loop.t0,E=E,p2=loop.p2).items()}
        qr=constant(loop.q);source=dict(q=qr[(0,0)],roots=roots)
        values=current.previous.all_u_first_bounds(f,source,qr,c.ln(c.mpf(p.nstr(scales.d_star,95))),c.mpf([0,1]))['values']
        A,B=values['A'],values['B_over_Pstar'];e=f.scalar('1.3');v=f.scalar('.4')
        pairs=current.leading_pairs(e,zero,v,zero,A,values['A_y'],B,values['B_y_over_Pstar'])
        capA=current.magnitude(f,A);capB=current.magnitude(f,B)
        exponent=c.exp(current.previous.phase.first.phase.bounded_value(capA*(c.mpf(1)/N)))
        Q=e*current.previous.phase.first.current.square(capA)*exponent*c.mpf('.5');F=e*capA*exponent
        rem=dict(m=zero,h=Q,k=v*Q+F*capB,
            e=e*Q+current.previous.phase.first.current.square(capB)+current.previous.phase.first.current.square(F)*c.mpf('.5'),
            p=e*Q+current.previous.phase.first.current.square(F)*c.mpf('.5'))
        primitive_source=dict(original_leading_zero_phase_mean=True,phase_primitive_endpoints_exactly_zero=True,
            original_leading_signed_density_pairs=pairs,source_phase_primitive_definition='original integral',
            source_phase_primitive_y_definition='original slow derivative integral')
        def data(angle):
            row=loop.at_angle(angle);aa,bb=row['A'],row['B'];de=E*p.expm1(aa/N);dv=bb/N
            lead=dict(m=bb,h=E*aa,k=E*(V*aa+bb),e=2*V*bb-E*E*aa,p=E*E*aa)
            density=dict(m=dv,h=de,k=V*de+E*dv+de*dv,e=2*V*dv+dv*dv-E*de-de*de/2,p=E*de+de*de/2)
            return row,lead,density
        phase0=p.mpf('.137');phase1=p.mpf('.554');width=(2+phase1-phase0)/N
        left_angle=loop.angle_at_phase(phase0);right_angle=loop.angle_at_phase(phase1)
        comparisons=0
        for key,rate in current.RATES.items():
            lam=p.mpf(rate.numerator)/rate.denominator
            def integral_piece(lower,upper,turn):
                def integrand(angle):
                    row,_,density=data(angle);yy=(turn+row['phase']-phase0)/N
                    return p.exp(-lam*(width-yy))*density[key]*row['phase_angle_derivative']/N
                return p.quad(integrand,[lower,(lower+upper)/2,upper])
            exact=integral_piece(left_angle,2*p.pi,0)+integral_piece(0,2*p.pi,1)+integral_piece(0,right_angle,2)
            cap,proof=current.weighted_cap(f,*pairs[key],rem[key],c.mpf(p.nstr(width,95)),c.mpf(0),rate,N,
                c.mpf('.137'),c.mpf('.446'))
            lo,hi=ep(current.symmetric(f,cap).finite_interval()+c.mpf(['-1e-70','1e-70']))
            assert lo<=exact<=hi,(sign,key,exact,lo,hi);comparisons+=1
            assert ep(proof['endpoint_primitive_cap'].coefficient)[1]>0
            # A separately integrated same original phase primitive, rather than its cap.
            G=p.quad(lambda angle:data(angle)[1][key]*data(angle)[0]['phase_angle_derivative'],[0,right_angle/2,right_angle])
            enclosed=current.phase_primitive_enclosure(f,primitive_source,c.mpf('.554'))
            lo,hi=ep(enclosed['original_phase_primitive_C0'][key].finite_interval()+c.mpf(['-1e-70','1e-70']))
            assert lo<=G<=hi,(sign,key,G,lo,hi);comparisons+=1
            mean=p.quad(lambda angle:data(angle)[1][key]*data(angle)[0]['phase_angle_derivative'],[0,p.pi,2*p.pi])
            assert abs(mean)<p.mpf('1e-75'),(sign,key,mean);comparisons+=1
            if not rate:
                assert ep(proof['incoming_cell_decay'])==(1,1) and ep(proof['downstream_suffix_decay'])==(1,1)
                bias=p.quad(lambda angle:data(angle)[2]['p']*data(angle)[0]['phase_angle_derivative'],[0,p.pi,2*p.pi])
                assert bias>0
        for phi in (0,1):
            enclosed=current.phase_primitive_enclosure(f,primitive_source,c.mpf(phi))
            assert all(row.zero for name in ('original_phase_primitive_C0','original_phase_primitive_slow_y') for row in enclosed[name].values())
        cases.append(dict(passed=True,signed_p2=sign,independent_weighted_integral_primitive_and_mean_comparisons=comparisons,
            nonlinear_pressure_mean_is_strictly_nonzero=True))
    rejects=0
    for boxes in ([],[c.mpf('-.1')],[c.mpf('1.1')]):
        try:current.phase_distance(c,boxes)
        except ValueError:rejects+=1
        else:raise AssertionError('Invalid phase primitive support accepted')
    assert ep(current.phase_distance(c,[c.mpf(['.4','.6'])]))==(mp.mpf('.5'),mp.mpf('.5'))
    return dict(passed=True,independent_original_signed_frozen_cases=cases,
        independent_original_comparisons=sum(row['independent_weighted_integral_primitive_and_mean_comparisons'] for row in cases),
        original_angle_to_phase_Jacobian_used_for_independent_quadrature=True,nonvertex_endpoint_terms_retained=True,
        exact_endpoint_primitive_zero_and_phase_guards_verified=True,phase_guard_rejections=rejects,
        numeric_diagnostic_precision90_allowance1e_minus70=True)


def native(owner,report):
    cells=endpoint_rows=leading_rows=IBP_uses=reductions=0;comparison={}
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(fields.serialized(packet))==report['frames'][label]
        op=owner.upstream.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['endpoint_retaining_C0_Nminus2_averaging_installed'] and not packet['ordinary_Z_Nminus2_averaging_installed']
        assert packet['real_finite_N_Rm_incoming_correction_is_still_unsupplied']
        assert len(packet['actual_source_weighted_cells'])==len(current.PARTITION)-1
        for cell in packet['actual_source_weighted_cells']:
            source=cell['source'];assert source['source_family']==owner.family
            assert source['exact_common_P0_axial5'] is op.P0 and source['exact_source_Rm_factor'] is op.Rm_factor
            assert source['original_slow_y_not_total_spatial_N_chain'] and source['original_full_density_mean_not_zeroed']
            assert source['phase_primitive_endpoints_exactly_zero'] and source['primitive_caps_are_bounds_not_function_values']
            assert ep(cell['source_log_width'])[0]>0 and ep(cell['source_suffix_to_Rh'])[0]>=0
            for pair in source['original_leading_signed_density_pairs'].values():
                for row in pair:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;leading_rows+=1
            for boxes in source['actual_endpoint_phase_primitive_enclosures'].values():
                for box in boxes:
                    for name in ('original_phase_primitive_C0','original_phase_primitive_slow_y'):
                        for row in box[name].values():
                            assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;endpoint_rows+=1
            for key,detail in cell['original_own_rate_averaging'].items():
                proof=detail['original_IBP']
                assert proof['all_cell_endpoints_retained'] and proof['nonlinear_mean_not_zeroed']
                assert proof['original_log_radius_measure_once'] and proof['internal_endpoint_cancellation_not_claimed']
                assert ep(proof['positive_own_rate_mass'])[0]>0
                if key=='p':assert ep(proof['incoming_cell_decay'])==(1,1) and ep(proof['downstream_suffix_decay'])==(1,1)
                cap=detail['used_valid_cap'];assert cap.zero or ep(cap.coefficient)[0]>=0
                if detail['used_route']=='left_valid_bound':IBP_uses+=1
            assert all(source[k] is False for k in fields.previous.OPEN);cells+=1
        baseline=owner.baseline['frames'][label]['actual_whole_patch_local_defect_integral_C0_Z']
        for key,rows in packet['actual_averaged_whole_patch_local_integral_C0_Z'].items():
            assert current.base.encoded(fields.serialized(rows[1]))==baseline[key][1]
            comp=packet['actual_C0_comparisons'][key]
            assert comp['ordinary_Z_accepted_direct_bound_unchanged']
            if comp['strict_absolute_upper_reduction']:
                assert ep(comp['log_absolute_upper_reduction'])[0]>0;reductions+=1
        assert all(packet[k] is False for k in fields.previous.OPEN)
        comparison[label]={key:dict(strict_reduction=v['strict_absolute_upper_reduction'],
            logarithmic_reduction=v['log_absolute_upper_reduction']) for key,v in packet['actual_C0_comparisons'].items()}
        print('Actual endpoint-retaining weighted Rm source integrals checked',label,flush=True)
    return dict(passed=True,actual_source_cells=cells,genuine_original_leading_C0_slow_y_rows=leading_rows,
        actual_endpoint_phase_primitive_rows=endpoint_rows,cell_rows_using_IBP_cap=IBP_uses,
        strict_whole_patch_C0_reductions=reductions,actual_conditional_C0_comparisons=comparison,
        accepted_ordinary_Z_bounds_bitwise_unchanged=True,actual_prefix_closure_global_N_recursion_not_admitted=True)


def run():
    began=time.monotonic()
    theorem=symbolic()
    with mp.workdps(240):independent=frozen_original_quadrature()
    print('Independent original signed phase primitive and weighted integral references PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmWeightedAveraging(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_leading_remainder_and_IBP=theorem,independent_original_frozen_quadrature=independent,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm source phase primitives and weighted C0 averaging PASS',flush=True);return result


if __name__=='__main__':run()
