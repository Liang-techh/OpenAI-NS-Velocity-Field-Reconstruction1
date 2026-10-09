"""Independent signed axial split, original weighted quadrature and live Z IBP."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_weighted_Z_averaging as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def contains(row,value,allowance='1e-65'):
    lo,hi=ep(row.finite_interval()+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def flow(c):return MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def symbolic():
    y,Z,phi=sy.symbols('y Z phi');N=sy.symbols('N',positive=True)
    E,V,A,B,RE,F=[sy.Function(name)(Z) for name in ('E','V','A','B','RE','F')]
    dE=E*A/N+RE/N**2;dV=B/N
    exact=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    leading=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)
    remainder=dict(m=sy.Integer(0),h=RE,k=V*RE+F*B,e=-E*RE+B*B-F*F/2,p=E*RE+F*F/2)
    for key in current.RATES:
        identity=exact[key]-leading[key]/N-remainder[key].subs(F,E*A+RE/N)/N**2
        assert sy.expand(sy.diff(identity,Z))==0
    expected=dict(m=sy.Integer(0),h=sy.diff(RE,Z),
        k=sy.diff(V,Z)*RE+V*sy.diff(RE,Z)+sy.diff(F,Z)*B+F*sy.diff(B,Z),
        e=-sy.diff(E,Z)*RE-E*sy.diff(RE,Z)+2*B*sy.diff(B,Z)-F*sy.diff(F,Z),
        p=sy.diff(E,Z)*RE+E*sy.diff(RE,Z)+F*sy.diff(F,Z))
    for key in current.RATES:assert sy.expand(sy.diff(remainder[key],Z)-expected[key])==0
    exactF=N*E*(sy.exp(A/N)-1)
    assert sy.simplify(sy.diff(exactF,Z)-(N*sy.diff(E,Z)*(sy.exp(A/N)-1)+E*sy.exp(A/N)*sy.diff(A,Z)))==0
    x,t=sy.symbols('x t',real=True)
    exprel=(sy.exp(x)-1)/x
    assert sy.simplify(exprel+x*sy.diff(exprel,x)-sy.exp(x))==0
    assert sy.integrate(1-t,(t,0,1))==sy.Rational(1,2)
    assert sy.integrate(t*(1-t),(t,0,1))==sy.Rational(1,6)
    G=sy.Function('G')(y,Z,phi);K=sy.Function('K')(y);lam=sy.symbols('lam',nonnegative=True)
    GZ=sy.diff(G,Z);total=lam*K*GZ+K*(sy.diff(GZ,y)+N*sy.diff(GZ,phi))
    assert sy.expand(K*sy.diff(GZ,phi)/N-(total-K*(sy.diff(GZ,y)+lam*GZ))/N**2)==0
    return dict(passed=True,all_five_exact_signed_density_Z_splits_verified=True,
        all_original_second_remainder_Z_product_terms_verified=True,
        original_F_N_Z_including_exprel_prime_verified=True,R2_and_R2_prime_integral_masses_verified=True,
        genuine_G_Z_G_yZ_endpoint_retaining_IBP_sign_and_Nminus2_verified=True,
        phase_Z_exact_zero_required=True)


def independent_remainders():
    c=MPIntervalContext();c.dps=160;f=flow(c);p=mp.mp.clone();p.dps=100;N=257;compared=0
    for av in ('-.3','0','.2'):
        E={current.C0:f.scalar('1.3'),current.Z:f.scalar('-.04')}
        V={current.C0:f.scalar('.4'),current.Z:f.scalar('.03')}
        values=dict(A=f.scalar(av),A_Z=f.scalar('.07'),B_over_Pstar=f.scalar('-.2'),B_Z_over_Pstar=f.scalar('.05'))
        caps,proof=current.remainder_Z_caps(f,E,V,values,N)
        def exact(zz):
            ee=p.mpf('1.3')-p.mpf('.04')*zz;vv=p.mpf('.4')+p.mpf('.03')*zz
            aa=p.mpf(av)+p.mpf('.07')*zz;bb=-p.mpf('.2')+p.mpf('.05')*zz
            RE=N*N*ee*(p.expm1(aa/N)-aa/N);F=N*ee*p.expm1(aa/N)
            return dict(Q=RE,F=F,m=p.mpf(0),h=RE,k=vv*RE+F*bb,
                e=-ee*RE+bb*bb-F*F/2,p=ee*RE+F*F/2)
        for name,row in proof.items():
            if name in ('original_R_E_cap','original_F_N_cap'):
                contains(current.averaging.symmetric(f,row),exact(0)['Q' if name=='original_R_E_cap' else 'F']);compared+=1
            elif name in ('original_R_E_Z_cap','original_F_N_Z_cap'):
                key='Q' if name=='original_R_E_Z_cap' else 'F'
                contains(current.averaging.symmetric(f,row),p.diff(lambda zz:exact(zz)[key],0));compared+=1
        for key,row in caps.items():
            contains(current.averaging.symmetric(f,row),p.diff(lambda zz:exact(zz)[key],0));compared+=1
    return dict(passed=True,independent_exact_exponential_remainder_and_Z_comparisons=compared,
        signed_A_and_nonzero_E_Z_A_Z_B_Z_V_Z_used=True,nonlinear_Z_bias_not_zeroed=True)


def original_weighted_quadrature():
    c=MPIntervalContext();c.dps=170;f=flow(c)
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='1000',p2_abs_max='1000',dps=90)
    p=scales.ctx;N=257;cases=[];zero=f.scalar(0)
    phase0=p.mpf('.137');phase1=p.mpf('.554');width=(2+phase1-phase0)/N
    alpha,beta,gamma,delta,epsilon=[p.mpf(v) for v in ('.07','.03','.02','.04','.05')]
    ybox=c.mpf([0,p.nstr(width,95)])
    e=f.scalar('1.3')*c.exp(c.mpf('.07')*ybox)
    E={current.C0:e,current.Y:e*c.mpf('.07'),current.Z:e*c.mpf('.03'),current.YZ:e*c.mpf('.0021')}
    V={current.C0:f.scalar('.4')+f.scalar('.02')*ybox,current.Y:f.scalar('.02'),
       current.Z:f.scalar('.04')+f.scalar('.05')*ybox,current.YZ:f.scalar('.05')}
    constants=lambda value:{key:f.scalar(p.nstr(value,95)) if key==current.C0 else zero for key in current.previous.ORDERS}
    for sign in (-1,1):
        loop=original.GenericShearLoop(scales,a='.8',b='.2',p1=200,p2=sign*p.mpf('.2'),Utheta=1)
        roots=dict(a=constants(loop.a),t0=constants(loop.t0),E=E,p2=constants(loop.p2));qr=constants(loop.q)
        got=current.previous.all_u_mixed_bounds(f,dict(q=qr[current.C0],roots=roots),qr,c.ln(c.mpf(p.nstr(scales.d_star,95))),c.mpf([0,1]))
        leading=current.previous.leading_mixed(E,V,got['values'])
        remainder,_=current.remainder_Z_caps(f,E,V,got['values'],N)
        source=dict(original_five_signed_leading_density_C0_y_Z_yZ=leading,original_radius_phase_Z_exact_zero=True,
            original_leading_zero_phase_mean=True,phase_primitive_endpoints_exactly_zero=True)
        left_angle=loop.angle_at_phase(phase0);right_angle=loop.angle_at_phase(phase1)
        def data(angle,yy):
            row=loop.at_angle(angle);A,Bbar=row['A'],row['B']
            ee=p.mpf('1.3')*p.exp(alpha*yy);ez=beta*ee;ey=alpha*ee;eyz=alpha*ez
            vv=p.mpf('.4')+gamma*yy;vz=delta+epsilon*yy;vy=gamma;vyz=epsilon
            B=ee*Bbar;By=ey*Bbar;Bz=ez*Bbar;Byz=eyz*Bbar
            leadZ=dict(m=Bz,h=ez*A,k=vz*ee*A+vv*ez*A+ez*B+ee*Bz,
                e=2*vz*B+2*vv*Bz-2*ee*ez*A,p=2*ee*ez*A)
            leadYZ=dict(m=Byz,h=eyz*A,
                k=(vyz*ee+vz*ey+vy*ez+vv*eyz)*A+eyz*B+ez*By+ey*Bz+ee*Byz,
                e=2*(vyz*B+vz*By+vy*Bz+vv*Byz)-2*(ey*ez+ee*eyz)*A,
                p=2*(ey*ez+ee*eyz)*A)
            de=ee*p.expm1(A/N);dez=ez*p.expm1(A/N);dv=B/N;dvz=Bz/N
            densityZ=dict(m=dvz,h=dez,k=vz*de+vv*dez+ez*dv+ee*dvz+dez*dv+de*dvz,
                e=2*vz*dv+2*vv*dvz+2*dv*dvz-ez*de-ee*dez-de*dez,
                p=ez*de+ee*dez+de*dez)
            return row,leadZ,leadYZ,densityZ
        comparisons=0
        for key,rate in current.RATES.items():
            lam=p.mpf(rate.numerator)/rate.denominator
            def piece(lower,upper,turn):
                def integrand(angle):
                    row=loop.at_angle(angle);yy=(turn+row['phase']-phase0)/N
                    return p.exp(-lam*(width-yy))*data(angle,yy)[3][key]*row['phase_angle_derivative']/N
                return p.quad(integrand,[lower,(lower+upper)/2,upper])
            exact=piece(left_angle,2*p.pi,0)+piece(0,2*p.pi,1)+piece(0,right_angle,2)
            cap,proof=current.averaging.weighted_cap(f,leading[key][current.Z],leading[key][current.YZ],remainder[key],
                c.mpf(p.nstr(width,95)),c.mpf(0),rate,N,c.mpf('.137'),c.mpf('.446'))
            contains(current.averaging.symmetric(f,cap),exact);comparisons+=1
            enclosed=current.phase_primitive_Z_enclosure(f,source,c.mpf('.554'))
            for index,name in ((1,'original_phase_primitive_Z'),(2,'original_phase_primitive_slow_yZ')):
                G=p.quad(lambda angle:data(angle,width/2)[index][key]*loop.at_angle(angle)['phase_angle_derivative'],
                    [0,right_angle/2,right_angle])
                contains(enclosed[name][key],G);comparisons+=1
            assert ep(proof['endpoint_primitive_cap'].coefficient)[1]>0
            if key=='p':assert ep(proof['incoming_cell_decay'])==(1,1) and ep(proof['downstream_suffix_decay'])==(1,1)
        for phi in (0,1):
            enclosed=current.phase_primitive_Z_enclosure(f,source,c.mpf(phi))
            assert all(row.zero for name in ('original_phase_primitive_Z','original_phase_primitive_slow_yZ') for row in enclosed[name].values())
        cases.append(dict(passed=True,signed_p2=sign,original_weighted_Z_G_Z_G_yZ_comparisons=comparisons))
    rejects=0
    for bad in (dict(source,original_radius_phase_Z_exact_zero=False),dict(source,original_leading_zero_phase_mean=False),
                dict(source,phase_primitive_endpoints_exactly_zero=False),dict(source,original_five_signed_leading_density_C0_y_Z_yZ={})):
        try:current.phase_primitive_Z_enclosure(f,bad,c.mpf('.2'))
        except ValueError:rejects+=1
        else:raise AssertionError('Invalid axial primitive source accepted')
    return dict(passed=True,independent_original_signed_cases=cases,
        independent_original_weighted_Z_and_primitive_comparisons=sum(row['original_weighted_Z_G_Z_G_yZ_comparisons'] for row in cases),
        nonzero_genuine_E_yZ_V_yZ_and_slow_y_used=True,original_angle_to_phase_Jacobian_used=True,
        all_nonvertex_endpoints_retained=True,exact_endpoint_zeros_verified=True,source_guard_rejections=rejects,
        quadrature_precision90_allowance1e_minus65_diagnostic_only=True)


def native(owner,report):
    cells=leading_rows=endpoint_rows=IBP_uses=reductions=0;comparison={}
    try:owner.contribution('0',owner.baseline['candidate_N']+1)
    except ValueError:pass
    else:raise AssertionError('Cross-N comparison with a fixed baseline accepted')
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(current.previous.serialized(packet))==report['frames'][label]
        op=owner.upstream.upstream.upstream.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['endpoint_retaining_genuine_Z_Nminus2_averaging_installed']
        assert not packet['actual_spatial_mixed_yZ_chain_installed']
        assert packet['real_finite_N_Rm_incoming_correction_is_still_unsupplied']
        assert len(packet['actual_source_weighted_Z_cells'])==len(current.PARTITION)-1
        for cell in packet['actual_source_weighted_Z_cells']:
            source=cell['source'];current.require_phase_Z_zero(source)
            assert source['source_family']==owner.family and source['exact_common_P0_axial5'] is op.P0
            assert source['exact_source_Rm_factor'] is op.Rm_factor and source['original_slow_yZ_not_total_spatial_yZ']
            mixed=source['genuine_original_mixed_source']
            assert mixed['source_family']==owner.family and mixed['exact_common_P0_axial5'] is op.P0
            assert mixed['exact_source_Rm_factor'] is op.Rm_factor
            assert mixed['original_mixed_source_record']['original_radius_Z_exact_zero']
            assert source['original_remainder_Z_proof']['full_signed_remainder_derivatives_retained']
            for rows in source['original_five_signed_leading_density_C0_y_Z_yZ'].values():
                for row in rows.values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;leading_rows+=1
            for boxes in source['actual_endpoint_phase_primitive_Z_enclosures'].values():
                for box in boxes:
                    for name in ('original_phase_primitive_Z','original_phase_primitive_slow_yZ'):
                        for row in box[name].values():
                            assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;endpoint_rows+=1
            for detail in cell['original_own_rate_Z_averaging'].values():
                proof=detail['original_Z_IBP']
                assert proof['all_cell_endpoints_retained'] and proof['nonlinear_mean_not_zeroed']
                assert proof['original_log_radius_measure_once'] and proof['internal_endpoint_cancellation_not_claimed']
                if detail['used_route']=='left_valid_bound':IBP_uses+=1
            assert all(source[key] is False for key in fields.previous.OPEN);cells+=1
        old=owner.baseline['frames'][label]['actual_averaged_whole_patch_local_integral_C0_Z']
        for key,rows in packet['actual_averaged_whole_patch_local_integral_C0_Z'].items():
            assert current.base.encoded(fields.serialized(rows[0]))==old[key][0]
            comp=packet['actual_Z_comparisons'][key];assert comp['accepted_C0_bound_unchanged']
            used,current_old=comp['used_Z_cap'],comp['accepted_direct_whole_Z_cap']
            ratio=(used.scale-current_old.scale).evaluate()+c.ln(used.coefficient)-c.ln(current_old.coefficient)
            assert ep(ratio)[1]<mp.mpf('1e-450')
            if comp['strict_absolute_upper_reduction']:
                assert ep(comp['log_absolute_upper_reduction'])[0]>0;reductions+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        comparison[label]={key:dict(strict_reduction=row['strict_absolute_upper_reduction'],
            logarithmic_reduction=row['log_absolute_upper_reduction']) for key,row in packet['actual_Z_comparisons'].items()}
        print('Actual genuine weighted Rm Z integrals checked',label,flush=True)
    return dict(passed=True,actual_source_cells=cells,genuine_original_leading_C0_y_Z_yZ_rows=leading_rows,
        actual_endpoint_G_Z_G_yZ_rows=endpoint_rows,cell_rows_using_Z_IBP_cap=IBP_uses,
        strict_whole_patch_Z_reductions=reductions,actual_conditional_Z_comparisons=comparison,
        accepted_C0_bounds_bitwise_unchanged=True,all_Z_bounds_no_worse_than_accepted=True,
        cross_N_baseline_comparison_rejected=True,
        actual_prefix_closure_global_N_recursion_not_admitted=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(240):
        remainders=independent_remainders();independent=original_weighted_quadrature()
    print('Independent original weighted Z and nonlinear references PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmWeightedZAveraging(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_Z_split_and_IBP=theorem,independent_nonlinear_remainder_Z=remainders,
        independent_original_weighted_Z_quadrature=independent,actual_live_source=actual,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original Rm G_Z G_yZ and weighted Z averaging PASS',flush=True);return result


if __name__=='__main__':run()
