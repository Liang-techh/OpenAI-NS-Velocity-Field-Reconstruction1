"""Independent switch source, physical conversion and finite-N driver checks."""
import copy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_switch_finite_N as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as scalar
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def contains(row,value,allowance='1e-35'):
    lo,hi=ep(row.finite_interval(max_log=2000)+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def test_flow(c):
    return MacroFlow(c,c.ln(c.mpf('.001')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))


def symbolic():
    y=sy.symbols('y');R,F,S=sy.symbols('R F S',positive=True)
    phi,V,H,M,K,A,B,C=sy.symbols('phi V H M K A B C')
    raw=dict(m=M/S,h=sy.sqrt(R/2)*F*H/S,k=sy.sqrt(R/2)*F*K/S**2,
        e=A/S**2-R*F**2*B/S**2,p=R*F**2*C/S**2)
    ode={R:R,H:2*phi-2*H,M:V-M,K:2*phi*V-2*K,A:V**2-A,B:phi**2-2*B,C:phi**2-C}
    E=sy.sqrt(2*R)*F*phi/S;v=V/S
    wanted=dict(m=v-raw['m'],h=E-sy.Rational(3,2)*raw['h'],
        k=E*v-sy.Rational(3,2)*raw['k'],e=v*v-raw['e']-E*E/2,p=E*E/2)
    for name,expr in raw.items():
        derivative=sum(sy.diff(expr,var)*rhs for var,rhs in ode.items())
        assert sy.simplify(derivative-wanted[name])==0
    phiy,a=sy.symbols('phiy a');Ey=sy.sqrt(2*R)*F*(phiy+phi/2)/S
    assert sy.simplify((E-2*Ey).subs(phiy,-a*phi/2)-a*E)==0
    h,Y,s,l,r,rate=sy.symbols('h Y s l r rate',positive=True)
    assert sy.simplify(h+h+(Y-2*h)-Y)==0
    assert sy.simplify(sy.exp(h*(1+s))-(sy.exp(h)*sy.exp(h*s)))==0
    assert sy.simplify(sy.exp(-rate*h)**2*sy.exp(-rate*(Y-2*h))-sy.exp(-rate*Y))==0
    sigma=sy.Function('sigma');Dbar=sy.Function('Dbar');u=sy.symbols('u')
    first_log=-h*h*sy.Integral(Dbar(u),(u,0,s))/2
    second_log=-h*h*sy.Integral((1-sigma(u))*Dbar(u),(u,0,s))/2-sy.Rational(2,5)*h*sy.Integral(sigma(u),(u,0,s))
    assert sy.simplify(-2*sy.diff(first_log,s)/h-h*Dbar(s))==0
    assert sy.simplify(-2*sy.diff(second_log,s)/h-(h*Dbar(s)*(1-sigma(s))+sy.Rational(4,5)*sigma(s)))==0
    az,D,Dz,eta=sy.symbols('az D Dz eta',real=True);z=sy.symbols('z')
    root=sy.sqrt((2*eta-D-Dz*z)/(2*(a+az*z)))
    assert sy.simplify(sy.diff(root,z).subs(z,0)/root.subs(z,0)+Dz/(2*(2*eta-D))+az/(2*a))==0
    return dict(passed=True,all_five_cumulative_macro_to_common_ODEs_verified=True,
        physical_sqrt_radius_and_one_Pstar_conversion_verified=True,exact_correlated_C_equals_a_E_verified=True,
        true_micro_jacobian_and_total_R100_R110_width_verified=True,
        full_three_window_memory_identity_verified=True,true_variable_a_q_root_Z_identity_verified=True)


def manufactured_source():
    c=MPIntervalContext();c.dps=160;f=test_flow(c);p=mp.mp.clone();p.dps=55
    z=sy.symbols('z',real=True);Q=sy.Rational
    phi0=Q(7,10)+z/50;V0=Q(2,5)+3*z/100+z*z/100
    initial={name:Q(n,10)+z/100 for n,name in enumerate(current.first.moments.RATES,1)}
    P0=Q(3,5)+z/10+z*z/30;amplitude=sy.exp(z*Q(11,100)+z*z*Q(3,100))
    def jet(expr):
        return current.long.IntervalTaylor(c,[c.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),145))) for n in range(6)])
    zero=jet(sy.Integer(0));one=jet(sy.Integer(1))
    f.set_sources([zero,jet(Q(7,10)),zero],{part:[zero,zero,zero] for part in fields.PARTS},
        one,[f.scalar(0)]*6,one,jet(V0),[f.scalar(0)]*6)
    first=current.first.FirstSwitchFunctions(f,dict(phi=f.jet(jet(phi0)),V=f.jet(jet(V0))),
        {name:f.jet(jet(expr)) for name,expr in initial.items()})
    exit1=first.evaluate((1,1));second=current.second.SecondSwitchFunctions(first,exit1['actual_fields'],exit1['actual_six_histories'])
    ref=SimpleNamespace(flow=f,c=c,P0=f.jet(jet(P0)),z=jet(z),zrows=f.jet(jet(z)),delta=c.mpf('.01'))
    inlet=dict(original_F0_derivative_ratios_ordinary=f.jet(jet(amplitude)),
        original_F0_squared_derivative_ratios_ordinary=f.jet(jet(amplitude**2)))
    h=p.mpf('.001');D=p.mpf('3.5')*p.exp(2*h);k=h*h*D/2;Y=p.log(p.mpf('1.1'));L=Y-2*h
    sf=lambda x:sy.Float(str(x),65)
    def conv(rate,t,decay):
        return h*(p.exp(-decay*t)-p.exp(-rate*h*t))/(rate*h-decay)
    def first_shapes(t):
        ph=phi0*sf(p.exp(-k*t));src=dict(H=2*phi0,M=V0,K=2*phi0*V0,A=V0*V0,B=phi0*phi0,C=phi0*phi0)
        result={}
        for name,rate in current.first.moments.RATES.items():
            decay=k if name in ('H','K') else 2*k if name in ('B','C') else p.mpf(0)
            result[name]=initial[name]*sf(p.exp(-rate*h*t))+src[name]*sf(conv(rate,t,decay))
        return ph,result
    ph1,I1=first_shapes(p.mpf(1));cache={}
    def complement(t):
        key=str(t)
        if key not in cache:
            cache[key]=p.quad(lambda u:1-scalar.flat_step(p,u),[0,t/4,t/2,3*t/4,t])
        return cache[key]
    def factor2(t):
        comp=complement(t)
        return p.exp(-k*comp-p.mpf('.4')*h*(t-comp))
    def second_shapes(t):
        ph=ph1*sf(factor2(t));src=dict(H=2*ph1,M=V0,K=2*ph1*V0,A=V0*V0,B=ph1*ph1,C=ph1*ph1)
        result={}
        for name,rate in current.first.moments.RATES.items():
            power=1 if name in ('H','K') else 2 if name in ('B','C') else 0
            mass=p.quad(lambda u:p.exp(-rate*h*(t-u))*factor2(u)**power*h,[0,t/4,t/2,3*t/4,t]) if power else (1-p.exp(-rate*h*t))/rate
            result[name]=I1[name]*sf(p.exp(-rate*h*t))+src[name]*sf(mass)
        return ph,result
    ph2,I2=second_shapes(p.mpf(1));comparisons=kernel_comparisons=0
    for chart,tstr,left,right in (('first_switch','.137',(0,1),(1,4)),('second_switch','.137',(0,1),(1,4)),
                                ('post_power','.554',(1,2),(3,4))):
        t=p.mpf(tstr)
        source=(current.post_cell(first,second,current.downstream.scalar_hull(c,c.mpf('.5'),c.mpf('.75')))
                if chart=='post_power' else current.prefix_cell(first,second,chart,left,right))
        proxy,packet,recovered,a=current.recover_source(ref,inlet,source)
        if chart=='first_switch':
            ph,hist=first_shapes(t);logR=p.log(100)+h*t;av=h*D
        elif chart=='second_switch':
            ph,hist=second_shapes(t);logR=p.log(100)+h*(1+t)
            sig=scalar.flat_step(p,t);av=h*D*(1-sig)+p.mpf('.8')*sig
        else:
            theta=p.exp(-L*t);ph=ph2*sf(theta**p.mpf('.4'));logR=p.log(100)+2*h+L*t;av=p.mpf('.8')
            hist=dict(H=I2['H']*sf(theta**2)+ph2*sf((theta**p.mpf('.4')-theta**2)*p.mpf('1.25')),
                M=I2['M']*sf(theta)+V0*sf(1-theta),K=I2['K']*sf(theta**2)+ph2*V0*sf((theta**p.mpf('.4')-theta**2)*p.mpf('1.25')),
                A=I2['A']*sf(theta)+V0*V0*sf(1-theta),B=I2['B']*sf(theta**2)+ph2*ph2*sf((theta**p.mpf('.8')-theta**2)*p.mpf(5)/6),
                C=I2['C']*sf(theta)+ph2*ph2*sf((theta**p.mpf('.8')-theta)*5))
        E=sf(p.sqrt(2*p.exp(logR))*p.mpf(2)/3)*amplitude*ph;V=V0/3
        own=dict(m=hist['M']/3,h=sf(p.sqrt(p.exp(logR)/2)*p.mpf(2)/3)*amplitude*hist['H'],
            k=sf(p.sqrt(p.exp(logR)/2)*p.mpf(2)/9)*amplitude*hist['K'],
            e=hist['A']/9-sf(p.exp(logR)*p.mpf(4)/9)*amplitude**2*hist['B'],
            p=sf(p.exp(logR)*p.mpf(4)/9)*amplitude**2*hist['C'])
        def actual(expr,n):return p.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),50)))
        for name,expr in dict(phi=ph,V=V0).items():
            for n in range(3):contains(source['fields'][name][n],actual(expr,n));comparisons+=1
        for name,expr in hist.items():
            for n in range(3):contains(source['histories'][name][n],actual(expr,n));comparisons+=1
        for name,expr in dict(E=E,V=V).items():
            for n in range(3):contains(recovered['common_velocity_'+name+'_axial5'][n],actual(expr,n));comparisons+=1
        for name,expr in own.items():
            for n in range(3):contains(recovered['common_own_five_histories_axial5'][name][n],actual(expr,n));comparisons+=1
        contains(a[0],av);contains(packet['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][1][0],actual(E*(1-sf(av))/2,0));comparisons+=2
        assert recovered['common_original_P0_axial5'] is ref.P0
        for name,rate in current.RATES.items():
            mass,decay,tail=current.own_weights(first,chart,left,right,rate)
            lo=p.mpf(left[0])/left[1];hi=p.mpf(right[0])/right[1];length=h if chart!='post_power' else L
            lam=p.mpf(rate.numerator)/rate.denominator;w=length*(hi-lo)
            contains(mass,w if not rate else -p.expm1(-lam*w)/lam)
            contains(decay,p.exp(-lam*w));contains(tail,p.exp(-lam*length*(1-hi)));kernel_comparisons+=3
    # Separate original nonzero hydro drive; its variable V is not replaced
    # by the constant-V fixture used for the independent six-history test.
    f.set_sources([zero,jet(Q(7,10)),zero],{part:[zero,jet(Q(1,25)) if part=='hydro' else zero,zero]
        for part in fields.PARTS},one,[f.scalar(0)]*6,one,jet(V0),[f.scalar(0)]*6)
    driven=current.first.FirstSwitchFunctions(f,dict(phi=f.jet(jet(phi0)),V=f.jet(jet(V0))),
        {name:f.jet(jet(expr)) for name,expr in initial.items()})
    source=current.prefix_cell(driven,None,'first_switch',(0,1),(1,4));t=p.mpf('.137');G=p.mpf('.2')*p.exp(2*h)
    mass=p.quad(lambda u:(1-scalar.flat_step(p,u))*p.exp(-k*u),[0,t/4,t/2,3*t/4,t])
    expectedV=V0-phi0*sf(h*h*G*mass)
    expectedVy=-phi0*sf(h*G*p.exp(-k*t)*(1-scalar.flat_step(p,t)))
    _,_,recovered,_=current.recover_source(ref,inlet,source);nonzero=0
    for rows,expr in ((source['fields']['V'],expectedV),(source['V_y'],expectedVy),
                      (recovered['actual_generic_source_numerators']['B'],2*expectedVy/3)):
        for n in range(3):contains(rows[n],actual(expr,n));nonzero+=1
    assert any(not row.zero for row in source['V_y'])
    return dict(passed=True,independent_original_switch_source_and_conversion_coefficients=comparisons,
        independent_true_measure_Duhamel_weight_comparisons=kernel_comparisons,
        independent_nonzero_first_switch_drive_and_V_y_coefficients=nonzero,
        actual_nonendpoint_whole_cell_values_tested=True,cumulative_histories_not_instantaneous_reset=True,
        separate_analytic_P0_and_nonzero_F0_axial_derivatives=True,
        independently_quadratured_original_second_switch_cutoff_prefix=True)


def variable_q_and_bounded_exponent():
    c=MPIntervalContext();c.dps=150;f=test_flow(c);p=mp.mp.clone();p.dps=80;comparisons=unions=flats=0
    for av,blo,bhi in (('.003','.01','.01'),('.8','.8','.8'),('.8','1.1','1.1'),('.8','1','1.1')):
        az=p.mpf('.00002') if av=='.003' else p.mpf('.04');bz=p.mpf('.03');eta=p.mpf('.2')
        a=[f.h*3 if av=='.003' else f.scalar(c.mpf(av)),f.scalar(c.mpf(az))]+[f.scalar(0)]*4
        E=[f.scalar(1)]+[f.scalar(0)]*5;b=[f.scalar(c.mpf([blo,bhi])),f.scalar(c.mpf(bz))]+[f.scalar(0)]*4
        rec=dict(actual_generic_source_numerators=dict(E=E,B=b),full_signed_inertial_sectors_axial4={name:[f.scalar(0)]*5 for name in ('theta_linear','theta_quadratic','axial_linear','axial_quadratic')})
        proxy=SimpleNamespace(flow=f,c=c,Pstar=f.scalar(3),source_radius=f.scalar(100))
        _,qr,proof=current.general_quotients(proxy,rec,a,c.ln(c.mpf(eta)))
        for bb in sorted(set((p.mpf(blo),p.mpf(bhi),(p.mpf(blo)+p.mpf(bhi))/2))):
            def original(zz):
                aa=p.mpf(av)+az*zz;b0=bb+bz*zz;D=aa+b0*b0/aa-2
                return p.mpf(0) if D>=eta else scalar.flat_step(p,1-D/eta)*p.sqrt((2*eta-D)/(2*aa))
            contains(qr[current.C0],original(0));contains(qr[current.Z],p.diff(original,0));comparisons+=2
        if proof['original_q_Z_scope'].get('mixed_branch_union'):unions+=1
        if qr[current.C0].zero:assert qr[current.Z].zero;flats+=1
    original=f.factor((0,0,0,0,0),c.mpf(['-1e100',str(mp.log(mp.mpf('.4')))]))*c.mpf([-1,1])
    az=f.scalar('.03');got=dict(values=dict(A=original,A_Z=az),record={})
    bounded=current.bounded_exponent_cover(f,got,257)
    contains(bounded['values']['A'],p.mpf('.4'),allowance='1e-14');contains(bounded['values']['A'],-p.mpf('.4'),allowance='1e-14')
    assert bounded['values']['A_Z'] is az
    return dict(passed=True,independent_variable_a_nonzero_b_q_C0_Z_comparisons=comparisons,
        explicit_mixed_active_flat_cells=unions,exact_flat_cells=flats,
        nonzero_formal_a_scale_accepted=True,bounded_A_cover_does_not_replace_source_a_or_A_Z=True)


def native(owner,report):
    cells=rows=densities=memory=nonzero_shear=bounded=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.downstream.reference.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.downstream.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['genuine_finite_N_R100_correction_still_unsupplied']
        assert not packet['actual_finite_N_Rm_incoming_correction_supplied']
        assert packet['accepted_downstream_driver_rehydrated_in_same_live_source_algebra']
        for chart,window in packet['actual_source_windows'].items():
            for cell in window['actual_full_source_cells']:
                source=cell['source'];generic=source['original_generic_source'];back=source['actual_background_source']
                assert source['exact_common_P0_axial5'] is op.P0 and generic['common_original_P0_axial5'] is op.P0
                assert back['source_interval_functions_not_endpoint_hulls']
                for values in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                    for row in values:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;rows+=1
                for part in ('kernels','Z_derivatives'):
                    for row in source['original_signed_five_density_C0_Z'][part].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;densities+=1
                if chart=='first_switch':assert any(not row.zero for row in back['V_y']);nonzero_shear+=1
                else:assert all(row.zero for row in back['V_y'])
                if source['original_primitive_proof'].get('bounded_A_C0_cover_only_not_source_or_A_Z_replacement'):bounded+=1
                assert all(source[key] is False for key in fields.previous.OPEN);cells+=1
            for row in window['retained_actual_incoming_memory'].values():assert not row.zero;memory+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        print('Current original switch source/driver checked',label,flush=True)
    rejected=0;op=owner.first.owner('0');two=owner.second.owner('0')
    for chart,left,right in (('wrong',(0,1),(1,1)),('first_switch',(-1,4),(1,4)),('second_switch',(3,4),(1,4))):
        try:current.prefix_cell(op,two,chart,left,right)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original source phase accepted')
    try:owner.contribution('0',N=258)
    except ValueError:rejected+=1
    else:raise AssertionError('Cross-N driver composition accepted')
    accepted=owner.saved_downstream['frames']['0'];P0=owner.downstream.reference.owner('0').P0
    for key,value in (('source_family','wrong'),('source_frame','.5'),('candidate_N',258),('exact_common_P0_axial5',[])):
        bad=copy.deepcopy(accepted);bad[key]=value
        try:current.guard_saved_source(owner.family,'0',257,P0,bad)
        except ValueError:rejected+=1
        else:raise AssertionError('Mismatched downstream source accepted')
    return dict(passed=True,actual_full_switch_source_cells=cells,actual_common_source_coefficients=rows,
        original_finite_N_signed_density_C0_Z_rows=densities,retained_window_memory_rows=memory,
        genuine_nonzero_first_switch_V_y_cells=nonzero_shear,bounded_A_C0_exponent_cells=bounded,
        source_range_family_frame_P0_and_cross_N_rejections=rejected,
        complete_R100_Rm_local_operator_not_real_R100_or_Rm_inlet=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(210):q=variable_q_and_bounded_exponent();source=manufactured_source()
    print('Independent switch source/measure/conversion and variable-a q PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalSwitchFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_source_conversion_and_measure=theorem,independent_variable_a_q=q,
        independent_original_switch_source=source,actual_live_source=actual,
        original_source_bindings=owner.bindings,input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Current original R100/R110 finite-N source and complete local R100/Rm drivers PASS',flush=True);return result


if __name__=='__main__':run()
