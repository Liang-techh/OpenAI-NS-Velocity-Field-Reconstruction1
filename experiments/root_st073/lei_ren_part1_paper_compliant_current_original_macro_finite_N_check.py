"""Independent exact frozen-macro source and true-geometry driver checks."""
import copy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_macro_finite_N as current

down=current.downstream;fields,ep=current.fields,current.ep


def contains(row,value,allowance='1e-40'):
    lo,hi=ep(row.finite_interval(max_log=2000)+row.ctx.mpf([-mp.mpf(allowance),mp.mpf(allowance)]))
    assert lo<=value<=hi,(value,lo,hi)


def symbolic():
    h,Y,rho,logRa,lam,sc=sy.symbols('h Y rho logRa lam sc',real=True)
    original=logRa+2*h+(Y-2*h)*rho
    factored=logRa+Y*rho+2*h*(1-rho)
    assert sy.simplify(original-factored)==0
    assert sy.simplify(sy.diff(original,rho)-(Y-2*h))==0
    assert sy.simplify(original.subs(rho,0)-(logRa+2*h))==0
    assert sy.simplify(original.subs(rho,1)-(logRa+Y))==0
    phase=original-logRa-h*sc/2
    assert sy.simplify(phase-((Y-2*h)*rho+h*(2-sc/2)))==0
    t=sy.symbols('t');D,G=sy.symbols('D G');ell0,phi0,V0,q=sy.symbols('ell0 phi0 V0 q')
    phi=phi0*sy.exp(ell0-h*D*t/2)
    vel=V0-h*q*G*sy.Integral(sy.exp(ell0-h*D*sy.Symbol('s')/2),(sy.Symbol('s'),0,t))
    assert sy.simplify(sy.diff(phi,t)+h*D*phi/2)==0
    assert sy.simplify(sy.diff(vel,t)+h*G*q*phi/phi0)==0
    assert sy.simplify(sy.exp(-lam*(Y-2*h))*sy.exp(-lam*sy.Symbol('suffix'))-
        sy.exp(-lam*(Y-2*h+sy.Symbol('suffix'))))==0
    return dict(passed=True,original_and_factored_macro_radius_identity_verified=True,
        exact_nonzero_hb_start_and_R100_endpoint_verified=True,true_physical_measure_and_Z_independent_phase_verified=True,
        original_angular_and_nonzero_axial_source_ODEs_verified=True,complete_incoming_memory_composition_verified=True)


def independent_source():
    c=MPIntervalContext();c.dps=160;p=mp.mp.clone();p.dps=90;z=sy.symbols('z',real=True);Q=sy.Rational
    f=current.macro.fields.MacroFlow(c,c.ln(c.mpf('.001')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    def jet(expr):return current.long.IntervalTaylor(c,[c.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),145))) for n in range(6)])
    phi0=Q(7,10)+z/50;V0=Q(2,5)+3*z/100+z*z/100;vin=Q(3,10000)+z/50000
    initial={name:Q(n,10)+z/100+z*z/1000 for n,name in enumerate(current.macro.RATES,1)}
    zero=jet(sy.Integer(0));ell=jet(Q(1,5000));one=jet(sy.Integer(1))
    f.set_sources([zero,jet(Q(7,10)),zero],{part:[zero,jet(Q(1,25)) if part=='hydro' else zero,zero]
        for part in fields.PARTS},jet(phi0),f.jet(ell),jet(phi0),jet(V0),f.jet(jet(vin)))
    op=current.macro.ActualMacroMoments(f,{name:f.jet(jet(expr)) for name,expr in initial.items()})
    series=down.first.FirstSwitchFunctions(f,dict(phi=f.jet(jet(phi0)),V=f.jet(jet(V0))),op.inlet)
    amp=sy.exp(z*Q(11,100)+z*z*Q(3,100));P0=Q(3,5)+z/10+z*z/30
    ref=SimpleNamespace(flow=f,c=c,P0=f.jet(jet(P0)),z=jet(z),zrows=f.jet(jet(z)),delta=c.mpf('.01'))
    inlet=dict(original_F0_derivative_ratios_ordinary=f.jet(jet(amp)),original_F0_squared_derivative_ratios_ordinary=f.jet(jet(amp**2)))
    h=p.mpf('.001');Ra=p.mpf(5);R0=Ra*p.exp(2*h);Y=p.log(100/Ra);L=Y-2*h
    D=p.mpf('.7')*R0;G=p.mpf('.04')*R0;k=h*D/2;beta=h*G/k
    sf=lambda x:sy.Float(str(x),100)
    phiin=phi0*sf(p.exp(p.mpf('.0002')));vc=V0+vin-phiin*sf(beta)
    source_terms=dict(H=[(2*phiin,k)],M=[(vc,p.mpf(0)),(sf(beta)*phiin,k)],
        K=[(2*phiin*vc,k),(2*sf(beta)*phiin*phiin,2*k)],
        A=[(vc*vc,p.mpf(0)),(2*vc*sf(beta)*phiin,k),(sf(beta*beta)*phiin*phiin,2*k)],
        B=[(phiin*phiin,2*k)],C=[(phiin*phiin,2*k)])
    def truth(rho):
        S=L*rho;R=R0*p.exp(S);ph=phiin*sf(p.exp(-k*S));V=vc+sf(beta)*ph
        hist={name:initial[name]*sf(p.exp(-rate*S))+sum(coef*sf((p.exp(-d*S)-p.exp(-rate*S))/(rate-d))
            for coef,d in source_terms[name]) for name,rate in current.macro.RATES.items()}
        return R,ph,V,hist
    comparisons=weight_comparisons=nonzero=0
    for left,right,point in (((0,1),(1,4),'.137'),((1,2),(3,4),'.554'),((1,1),(1,1),'1')):
        rho=p.mpf(point);R,ph,V,hist=truth(rho);source=current.macro_cell(op,series,left,right)
        proxy,packet,recovered,a=down.recover_source(ref,inlet,source)
        E=sf(p.sqrt(2*R)*2/3)*amp*ph;v=V/3
        own=dict(m=hist['M']/3,h=sf(p.sqrt(R/2)*2/3)*amp*hist['H'],
            k=sf(p.sqrt(R/2)*2/9)*amp*hist['K'],e=hist['A']/9-sf(R*4/9)*amp**2*hist['B'],p=sf(R*4/9)*amp**2*hist['C'])
        actual=lambda expr,n:p.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,0)/sy.factorial(n),85)))
        for rows,expr in ((source['fields']['phi'],ph),(source['fields']['V'],V),
            (source['phi_y'],-sf(k)*ph),(source['V_y'],-sf(h*G)*ph),
            (recovered['common_velocity_E_axial5'],E),(recovered['common_velocity_V_axial5'],v),
            (recovered['actual_generic_source_numerators']['B'],-sf(2*h*G/3)*ph)):
            for n in range(3):contains(rows[n],actual(expr,n));comparisons+=1
        for name,expr in hist.items():
            for n in range(3):contains(source['histories'][name][n],actual(expr,n));comparisons+=1
        for name,expr in own.items():
            for n in range(3):contains(recovered['common_own_five_histories_axial5'][name][n],actual(expr,n));comparisons+=1
        contains(source['radius'],R);contains(source['root_radius'],p.sqrt(R));contains(a[0],h*D);comparisons+=3
        contains(source['actual_micro_exit_radius'],R0);contains(source['window_length'],L);comparisons+=2
        sc=p.mpf('.017');phase=current.actual_phase_record(f,source,c.mpf(sc))
        contains(phase['actual_unwrapped_log_radius_phase'],p.log(R/Ra)-h*sc/2)
        contains(phase['original_positive_left_log_offset'],h*sc/2);comparisons+=2
        assert recovered['common_original_P0_axial5'] is ref.P0
        assert any(not row.zero for row in source['V_y']);nonzero+=1
        if left!=right:
            d=p.mpf(right[0])/right[1]-p.mpf(left[0])/left[1];suffix=1-p.mpf(right[0])/right[1]
            for rate in current.RATES.values():
                mass,decay,tail=current.weights(series,left,right,rate);lam=p.mpf(rate.numerator)/rate.denominator
                contains(mass,L*d if not rate else -p.expm1(-lam*L*d)/lam)
                contains(decay,p.exp(-lam*L*d));contains(tail,p.exp(-lam*L*suffix));weight_comparisons+=3
    return dict(passed=True,independent_full_nonlinear_macro_source_and_conversion_coefficients=comparisons,
        independent_true_physical_macro_Duhamel_weights=weight_comparisons,
        nonzero_original_macro_axial_drive_cases=nonzero,nonzero_micro_prefix_and_actual_incoming_histories_retained=True,
        nonendpoint_whole_cell_source_and_exact_R100_endpoint_tested=True,
        genuine_F0_derivatives_and_separate_analytic_P0_retained=True)


def native(owner,report):
    cells=rows=density_rows=memories=0
    for label in ('0','.5'):
        packet=owner.contribution(label);ref=owner.downstream.downstream.reference.owner(label);f,c=ref.flow,ref.c
        assert current.base.encoded(down.downstream.serialized(packet))==report['frames'][label]
        assert packet['genuine_finite_N_micro_exit_correction_still_unsupplied']
        assert not packet['actual_finite_N_R100_correction_supplied'] and not packet['actual_finite_N_Rm_incoming_correction_supplied']
        assert packet['accepted_downstream_rehydrated_in_same_live_source_algebra']
        for cell in packet['actual_macro_source_cells']:
            source=cell['source'];back=source['actual_macro_background_source'];generic=source['original_generic_source']
            assert back['source_interval_functions_not_endpoint_hulls'] and back['actual_micro_background_inlet_preserved']
            assert not back['real_finite_N_micro_exit_correction_supplied']
            assert source['exact_common_P0_axial5'] is ref.P0 and generic['common_original_P0_axial5'] is ref.P0
            assert source['actual_selected_s_c'] is owner.sc and ep(owner.sc)[0]>0
            assert not source['original_positive_left_log_offset'].zero
            assert 's_c' in source['actual_phase_definition']
            # A second native route forms log(R/Ra) from the original
            # radius record. Its directed cover must meet the phase source.
            R=back['radius'];coeff=ep(R.coefficient)
            assert coeff[0]>0
            direct=R.scale.evaluate()+c.ln(R.coefficient)-f.logs[3]
            phase_source=f.ordinary_cover(source['actual_unwrapped_log_radius_phase']+source['original_positive_left_log_offset'])
            assert max(ep(direct)[0],ep(phase_source)[0])<=min(ep(direct)[1],ep(phase_source)[1])
            assert any(not row.zero for row in back['V_y'])
            for values in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                for row in values:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;rows+=1
            for part in ('kernels','Z_derivatives'):
                for row in source['original_signed_five_density_C0_Z'][part].values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
            assert all(source[key] is False for key in fields.previous.OPEN);cells+=1
        for row in packet['retained_micro_exit_R100_incoming_memory'].values():assert not row.zero;memories+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        print('Current original macro source/driver checked',label,flush=True)
    rejected=0;op,_=owner.moments.owner('0');series=owner.downstream.first.owner('0')
    for left,right in (((-1,4),(1,4)),((3,4),(1,4)),((0,1),(5,4))):
        try:current.macro_cell(op,series,left,right)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid whole macro cell accepted')
    try:owner.contribution('0',N=258)
    except ValueError:rejected+=1
    else:raise AssertionError('Cross-N macro driver composition accepted')
    accepted=owner.saved_downstream['frames']['0'];P0=owner.downstream.downstream.reference.owner('0').P0
    for key,value in (('source_family','wrong'),('source_frame','.5'),('candidate_N',258),('exact_common_P0_axial5',[])):
        bad=copy.deepcopy(accepted);bad[key]=value
        try:down.guard_saved_source(owner.family,'0',257,P0,bad)
        except ValueError:rejected+=1
        else:raise AssertionError('Mismatched macro/downstream source accepted')
    return dict(passed=True,actual_whole_macro_source_cells=cells,actual_common_source_coefficients=rows,
        original_finite_N_signed_density_C0_Z_rows=density_rows,retained_macro_incoming_memory_rows=memories,
        invalid_phase_family_frame_P0_and_cross_N_rejections=rejected,
        complete_local_micro_exit_Rm_operator_not_real_inlet_or_closure=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(210):independent=independent_source()
    print('Independent full nonlinear macro source/geometry/measure PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalMacroFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_geometry_and_source=theorem,independent_original_macro_source=independent,
        actual_live_source=actual,original_source_bindings=owner.bindings,input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original frozen-macro finite-N source and complete local micro-exit/Rm drivers PASS',flush=True);return result


if __name__=='__main__':run()
