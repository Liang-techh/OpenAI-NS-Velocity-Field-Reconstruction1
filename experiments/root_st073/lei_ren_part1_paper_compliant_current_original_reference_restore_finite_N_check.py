"""Independent original source, general q union and full downstream driver checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as scalar
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def contains(row,value,allowance='1e-60'):
    lo,hi=ep(row.finite_interval(max_log=2000)+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def test_flow(c):
    return MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def symbolic():
    a,b,az,bz,eta,D,Dz=sy.symbols('a b az bz eta D Dz',real=True)
    z=sy.symbols('z');aa=a+az*z;bb=b+bz*z
    assert sy.simplify(sy.diff(aa+bb*bb/aa-2,z).subs(z,0)-(az*(1-b*b/a**2)+2*b*bz/a))==0
    root=sy.sqrt((2*eta-D-Dz*z)/(2*aa))
    assert sy.simplify(sy.diff(root,z).subs(z,0)/root.subs(z,0)+Dz/(2*(2*eta-D))+az/(2*a))==0
    t0,p,q=sy.symbols('t0 p q',real=True);nu=1+t0*t0+2*q*q
    assert sy.simplify(((1+t0*t0)*2*sy.pi*p/nu/(2*sy.pi)).subs(q,0)-p)==0
    x,l1,l2,l3=sy.symbols('x l1 l2 l3',positive=True)
    assert sy.simplify(sy.exp(-x*l1)*sy.exp(-x*l2)*sy.exp(-x*l3)-sy.exp(-x*(l1+l2+l3)))==0
    t,T,logC,logP,logq,log110=sy.symbols('t T logC logP logq log110',real=True)
    logref=10*(logC+logP);gap=logref-T-8
    originalR=log110+T*(1-t)+(logref-8)*t
    currentR=log110+T+gap*t
    factoredR=5*(2*logP)+10*logC+log110-8-gap*(1-t)
    assert sy.simplify(originalR-currentR)==0 and sy.simplify(originalR-factoredR)==0
    logE=-logq+(T/10-logC-logP)*(1-t)-sy.Rational(4,5)*t
    assert sy.simplify(sy.diff(logE,t)/sy.diff(currentR,t)-sy.Rational(1,10))==0
    for offset in (-8+t,-7+t):
        assert sy.simplify(sy.diff(-logq+offset/10,t)/sy.diff(log110+logref+offset,t)-sy.Rational(1,10))==0
    z,V0,alpha,alphay=sy.symbols('z V0 alpha alphay',real=True)
    vel=4*z+(V0-4*z)*alpha
    assert sy.simplify(vel.subs(alpha,1)-V0)==0 and sy.simplify(vel.subs(alpha,0)-4*z)==0
    assert sy.simplify(sy.diff(vel,alpha)*alphay-(V0-4*z)*alphay)==0
    return dict(passed=True,full_b_squared_kappa_Z_and_q_root_derivative_verified=True,
        exact_flat_phase_inverse_and_primitive_identity_verified=True,
        full_window_memory_composition_verified=True,
        original_and_new_reference_log_radius_and_formal_factor_identity_verified=True,
        all_three_exact_physical_log_E_y_slopes_and_correlated_a_verified=True,
        original_restoration_velocity_forcing_and_nonzero_y_derivative_verified=True)


def q_union():
    c=MPIntervalContext();c.dps=170;f=test_flow(c);p=mp.mp.clone();p.dps=100
    eta=p.mpf('.2');az=p.mpf('.04');dz=p.mpf('.06');comparisons=unions=flats=0
    cases=[('-1.2','-1.2'),('-.05','-.05'),('.1','.1'),('.2','.2'),('.4','.4'),('.1','.3'),('-.1','.3')]
    for lo,hi in cases:
        a=[f.scalar(c.mpf('.8')),f.scalar(c.mpf(az))]
        Delta=[f.scalar(c.mpf([lo,hi])),f.scalar(c.mpf(dz))]
        rows,proof=current.general_q_C0_Z(f,a,Delta,c.ln(c.mpf(eta)))
        for D in sorted(set((p.mpf(lo),p.mpf(hi),(p.mpf(lo)+p.mpf(hi))/2))):
            def original(zz):
                d=D+dz*zz;aa=p.mpf('.8')+az*zz
                if d>=eta:return p.mpf(0)
                return scalar.flat_step(p,1-d/eta)*p.sqrt((2*eta-d)/(2*aa))
            contains(rows[current.C0],original(0));contains(rows[current.Z],p.diff(original,0));comparisons+=2
        if proof.get('mixed_branch_union'):
            assert proof['active_gamma_positive_intersection_scope_only']
            assert proof['exact_flat_q_and_q_Z_zero_only_on_flat_subdomain']
            assert not proof['q_globally_positive'] and proof['derivative_of_selected_cap_not_used'];unions+=1
        if rows[current.C0].zero:assert rows[current.Z].zero;flats+=1
    return dict(passed=True,independent_original_q_C0_Z_comparisons=comparisons,
        explicit_mixed_source_domain_unions=unions,exact_flat_cells=flats,
        active_domain_positive_theorem_not_globalized=True)


def manufactured_source():
    c=MPIntervalContext();c.dps=180;f=test_flow(c);p=mp.mp.clone();p.dps=100
    z=sy.symbols('z',real=True);R=sy.Rational
    initial=dict(theta=R(7,10)+z/50,theta_z=R(1,4)+z/100,pressure=R(5)+z/20,
        swirl=R(5,6)-z/30,mean=R(1,3)+z/40,axial=R(1,5)+z/50)
    V0=R(2,5)+z*R(3,100)+z*z/100;P0=R(3,5)+z/10+z*z/30
    def jet(expr):return current.long.IntervalTaylor(c,[c.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),165))) for n in range(6)])
    def full_kernels(t):
        return {name:p.quad(lambda s:p.exp(-p.mpf(rate)*(t-s))*(1-scalar.flat_step(p,s))**power,[0,t/4,t/2,3*t/4,t])
            for name,rate,power in (('mean',1,1),('mixed','1.6',1),('square',1,2))}
    end=full_kernels(p.mpf(1))
    endpoint_integrals={part:c.mpf(value)*c.exp(c.mpf(rate)) for part,value,rate in
        (('linear_1',end['mean'],1),('linear_8_over_5',end['mixed'],'1.6'),('square_1',end['square'],1))}
    provider=current.background.RestorationKernelProvider(c,endpoint_integrals)
    op=current.background.ActualReferenceRestoreFunctions(f,c.mpf(0),c.mpf('.01'),
        {name:f.jet(jet(expr)) for name,expr in initial.items()},f.jet(jet(V0)),f.jet(jet(P0)),c.mpf(2),c.mpf('1.5'),provider)
    forcing=V0-4*z;gap=10*(p.mpf('1.5')+p.log(3))-2-8
    cent0=dict(mean_error=initial['mean']-4*z,angular_error=initial['theta']-R(5,8),
        mixed_error=initial['theta_z']-4*z*initial['theta'],axial_square=initial['axial']-8*z*initial['mean']+16*z*z,
        swirl_error=initial['swirl']-R(5,6),pressure_error=initial['pressure']-5)
    sf=lambda value:sy.Float(str(value),110)
    def reference_cent(g):
        decay={name:sf(p.exp(-p.mpf(rate)*g)) for name,rate in current.background.RATES.items()}
        return dict(mean_error=forcing+(cent0['mean_error']-forcing)*decay['mean_error'],
            angular_error=cent0['angular_error']*decay['angular_error'],
            mixed_error=cent0['mixed_error']*decay['mixed_error']+forcing*(1-decay['mixed_error'])*R(5,8),
            axial_square=forcing*forcing+(cent0['axial_square']-forcing*forcing)*decay['axial_square'],
            swirl_error=cent0['swirl_error']*decay['swirl_error'],pressure_error=cent0['pressure_error']*decay['pressure_error'])
    centRz=reference_cent(gap)
    def restore_cent(t):
        decay={name:sf(p.exp(-p.mpf(rate)*t)) for name,rate in current.background.RATES.items()}
        result={name:centRz[name]*decay[name] for name in centRz};K=full_kernels(t)
        result['mean_error']+=forcing*sf(K['mean']);result['mixed_error']+=forcing*sf(K['mixed'])
        result['axial_square']+=forcing*forcing*sf(K['square']);return result
    centExit=restore_cent(p.mpf(1));comparisons=kernel_comparisons=nonzero_shear=0
    for chart in current.CHARTS:
        for tstr,left,right in (('.137',0,1),('.554',2,3)):
            t=p.mpf(tstr);source=current.background_cell(op,chart,c.mpf([c.mpf(left)/4,c.mpf(right)/4]))
            proxy,packet,recovered=current.recover_cell(op,source)
            roots,qrows,proof=current.quotient_cell(proxy,recovered,c.ln(c.mpf('.2')))
            if chart=='reference':
                cent=reference_cent(gap*t);V=V0;Vy=sy.Integer(0)
                offset=sf((p.mpf('.2')-p.mpf('1.5')-p.log(3))*(1-t)-p.mpf('.8')*t)
                logradius=sf(p.log(110)+2+gap*t)
            elif chart=='restoration':
                cent=restore_cent(t);alpha=1-scalar.flat_step(p,t);dy=-p.diff(lambda s:scalar.flat_step(p,s),t)
                V=4*z+forcing*sf(alpha);Vy=forcing*sf(dy);offset=sf((-8+t)/10)
                logradius=sf(p.log(110)+10*(p.mpf('1.5')+p.log(3))-8+t);nonzero_shear+=1
                for name,value in full_kernels(t).items():
                    lo,hi=ep(source['full_original_restoration_kernels'][name]);assert lo<=value<=hi;kernel_comparisons+=1
            else:
                cent={name:centExit[name]*sf(p.exp(-p.mpf(current.background.RATES[name])*t)) for name in centExit}
                V=4*z;Vy=sy.Integer(0);offset=sf((-7+t)/10)
                logradius=sf(p.log(110)+10*(p.mpf('1.5')+p.log(3))-7+t)
            H=cent['angular_error']+R(5,8)
            shapes=dict(theta=H,theta_z=4*z*H+cent['mixed_error'],mean=4*z+cent['mean_error'],
                axial=16*z*z+8*z*cent['mean_error']+cent['axial_square'],
                swirl=cent['swirl_error']+R(5,6),pressure=cent['pressure_error']+5)
            E=sy.exp(offset)/(1+z*z);vv=V/3
            hist=dict(m=shapes['mean']/3,h=E*shapes['theta'],k=E*shapes['theta_z']/3,
                e=shapes['axial']/9-E*E*shapes['swirl']/2,p=E*E*shapes['pressure']/2)
            de=R(1,100);L=1-de*z*z;d=1-z*z;m,h,k,e,pp=(hist[name] for name in current.RATES)
            transport=m*z*(1-de)+sy.diff(m,z)*d;pressure=P0+pp
            sectors=dict(theta_linear=(-E+h*(1-de/2)-sy.diff(h,z)*z*(1-de)/2)/L,
                theta_quadratic=(k*z*(2*de-1)-sy.diff(k,z)*d+E*transport)/L,
                axial_linear=(-vv+(m-sy.diff(m,z)*z)*(1-de)/2)/L,
                axial_quadratic=(vv*transport+e*z*2*de-sy.diff(e,z)*d+pressure*z*2*(1+de)-sy.diff(pressure,z)*d)/L)
            def actual(expr,n):return mp.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),105)))
            for name,expr in dict(E=E,V=vv).items():
                for n in range(3):contains(recovered['common_velocity_'+name+'_axial5'][n],actual(expr,n));comparisons+=1
            for name,expr in hist.items():
                for n in range(3):contains(recovered['common_own_five_histories_axial5'][name][n],actual(expr,n));comparisons+=1
            for name,expr in sectors.items():
                for n in (0,1):contains(recovered['full_signed_inertial_sectors_axial4'][name][n],actual(expr,n));comparisons+=1
            b=2*Vy/(3*E);t0=-b*R(5,4);Delta=R(4,5)+R(5,4)*b*b-2
            p2=sy.exp(logradius)*(sectors['axial_linear']+3*sectors['axial_quadratic'])/E
            for name,expr in dict(b=b,t0=t0,Delta=Delta,p2=p2).items():
                rows=(proof['actual_nonzero_b_axial5'] if name=='b' else proof['actual_Delta_axial5'] if name=='Delta'
                      else proof['actual_t0_axial5'] if name=='t0' else proof['full_signed_p2_axial4'])
                for n in (0,1):contains(rows[n],actual(expr,n));comparisons+=1
    return dict(passed=True,independent_full_original_source_and_inertia_coefficients=comparisons,
        independent_full_finite_restoration_kernel_comparisons=kernel_comparisons,
        nonzero_restoration_V_y_source_cases=nonzero_shear,
        whole_cell_source_covers_tested_at_nonendpoint_points=True,actual_P0_separate_and_one_Pstar_conversion=True)


def native(owner,report):
    cells=rows=densities=memory=nonzero_shear=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.reference.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['genuine_current_finite_N_R110_correction_still_unsupplied']
        assert not packet['actual_finite_N_Rm_incoming_correction_supplied']
        assert packet['accepted_long_driver_rehydrated_in_same_live_source_algebra']
        for chart,window in packet['actual_source_windows'].items():
            for cell in window['actual_full_source_cells']:
                source=cell['source'];generic=source['original_generic_source'];back=source['original_background_interval_source']
                assert source['exact_common_P0_axial5'] is op.P0 and generic['common_original_P0_axial5'] is op.P0
                assert generic['exact_E_y_equals_E_over10_and_a_four_fifths']
                assert back['radial_rows_use_physical_log_radius'] and back['original_full_source_window_not_shortened']
                for values in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                    for row in values:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;rows+=1
                for part in ('kernels','Z_derivatives'):
                    for row in source['original_signed_five_density_C0_Z'][part].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;densities+=1
                if chart=='restoration':
                    assert back['full_original_restoration_kernels'] is not None
                    assert any(not row.zero for row in back['actual_velocity_V_y']);nonzero_shear+=1
                else:assert all(row.zero for row in back['actual_velocity_V_y'])
                assert all(source[key] is False for key in fields.previous.OPEN);cells+=1
            for row in window['actual_nonzero_incoming_memory'].values():assert not row.zero;memory+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        assert all(not row.zero for row in packet['actual_R110_Rm_incoming_memory'].values())
        print('Current original downstream source windows checked',label,flush=True)
    rejected=0
    op=owner.reference.owner('0')
    for chart,t in (('wrong',0),('reference',-1),('restoration',2)):
        try:current.background_cell(op,chart,op.c.mpf(t))
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid actual source chart/cell accepted')
    try:owner.contribution('0',N=258)
    except ValueError:rejected+=1
    else:raise AssertionError('Cross-N long-driver composition accepted')
    return dict(passed=True,actual_full_downstream_source_cells=cells,actual_common_source_coefficients=rows,
        original_finite_N_signed_density_C0_Z_rows=densities,retained_window_memory_rows=memory,
        genuine_nonzero_restoration_V_y_cells=nonzero_shear,source_chart_range_and_cross_N_rejections=rejected,
        complete_downstream_drivers_not_real_R110_or_Rm_inlet=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(240):q=q_union();source=manufactured_source()
    print('Independent original downstream source and general-q union PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalReferenceRestoreFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_interval_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_source_and_q=theorem,independent_original_general_q=q,
        independent_full_original_downstream_source=source,actual_live_source=actual,
        original_interval_source_bindings=owner.bindings,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Current original full downstream finite-N source and Duhamel PASS',flush=True);return result


if __name__=='__main__':run()
