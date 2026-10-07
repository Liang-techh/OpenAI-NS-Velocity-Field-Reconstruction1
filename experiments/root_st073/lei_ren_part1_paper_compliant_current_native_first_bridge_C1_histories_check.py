"""Independent original loop/cutoff references and whole first-bridge replay."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_first_bridge_C1_histories as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks
import lei_ren_part1_paper_compliant_current_native_q_slow_jets_check as qchecks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=checks.require;contains=checks.contains;original=checks.original;density=current.density


def independent_original_loop_checks(c):
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='3',p1_abs_max='1e9',p2_abs_max='1e8',dps=120)
    p=scales.ctx;Z=sy.symbols('Z');eta_text=p.nstr(scales.eta,125);eta_s=sy.Rational(eta_text)
    shape=sy.Rational(1,10)*Z;aa=sy.Rational(4,5)+shape;bb=sy.Rational(1,5)-shape/2
    cases=[('body_signed_positive',aa,bb,sy.Rational(1,5)+shape),
        ('body_signed_negative',aa,-bb,-sy.Rational(1,5)+shape),
        ('signed_r_zero_crossing',aa,bb,shape),
        ('large_positive_u',aa,bb,sy.Integer(10**7)+shape),
        ('large_negative_u',aa,-bb,-sy.Integer(10**7)+shape),
        ('all_three_cutoff_branches',2+eta_s*(sy.Rational(1,2)+Z),sy.Integer(0),shape)]
    Eexpr=sy.Rational(13,10)+shape;count=0;qcount=0;regions={};h=p.mpf('1e-9');allow=c.mpf(('-1e-28','1e-28'))
    for label,aexpr,bexpr,p2expr in cases:
        bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger()
        scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
        zbox=c.mpf(('-1','1'))
        expressions=dict(a=aexpr,b=bexpr,p2=p2expr,E=Eexpr,kappa_minus2=aexpr+bexpr*bexpr/aexpr-2)
        def interval(expr):
            if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
            if expr==Z:return zbox
            if expr.is_Add:return sum((interval(term) for term in expr.args),c.mpf(0))
            if expr.is_Mul:
                value=c.mpf(1)
                for term in expr.args:value*=interval(term)
                return value
            if expr.is_Pow and expr.exp.is_Integer:return interval(expr.base)**int(expr.exp)
            raise ValueError('Unsupported exact interval fixture expression: '+str(expr))
        roots={name:{order:scalar(interval(sy.diff(expr,Z,order[1]))) for order in (current.ZERO,current.DZ)} for name,expr in expressions.items()}
        eta_log=c.ln(c.mpf(eta_text));dstar_log=c.ln(c.mpf(p.nstr(scales.d_star,125)))
        got=current.whole_period_C1(roots,eta_log,c.ln(c.mpf('.7')),dstar_log)
        qgot=current.whole_cutoff_C1(roots,eta_log,c.ln(c.mpf('.7')))
        # Explicit rational evaluation in the scalar reference context; no
        # fixture value is injected into the native source graph.
        def eval_expr(expr,z):
            terms=sy.Poly(expr,Z).terms()
            return sum((p.mpf(int(coefficient.p))/int(coefficient.q)*z**power[0] for power,coefficient in terms),p.mpf(0))
        def ref(z,phi,guess=None):
            av,bv,pv,ev=[eval_expr(expr,z) for expr in (aexpr,bexpr,p2expr,Eexpr)]
            loop=original.GenericShearLoop(scales,a=av,b=bv,p1=10**8,p2=pv,Utheta=ev)
            if loop.q==0 or phi in (0,p.mpf('.5'),1):psi=2*p.pi*phi
            else:
                start=2*p.pi*phi if guess is None else guess
                try:psi=p.findroot(lambda x:loop.phase_at_angle(x)-phi,(start,start+p.mpf('.001')),solver='secant',tol=p.mpf('1e-100'))
                except (ValueError,ArithmeticError):psi=loop.angle_at_phase(phi)
            return loop.at_angle(psi),psi,loop.q
        local=0;seen=set()
        for point in ('-.75','0','.75'):
            z=p.mpf(point)
            av,bv=[eval_expr(expr,z) for expr in (aexpr,bexpr)];Delta=av+bv*bv/av-2
            seen.add('flat' if Delta>=scales.eta else 'body' if Delta<=0 else 'transition')
            for phase in ('.137','.537'):
                phi=p.mpf(phase);value,guess,q=ref(z,phi)
                for key,name in (('A','A'),('B_over_Pstar','B')):
                    require(contains(got['values'][key].finite_interval()+allow,c.mpf(p.nstr(value[name],125))),
                        'Original whole-period value not enclosed: '+label+' '+key);count+=1;local+=1
                rows=[ref(z+shift*h,phi,guess)[0] for shift in (-2,-1,1,2)]
                for key,name in (('A_Z','A'),('B_Z_over_Pstar','B')):
                    derivative=(rows[0][name]-8*rows[1][name]+8*rows[2][name]-rows[3][name])/(12*h)
                    require(contains(got['values'][key].finite_interval()+allow,c.mpf(p.nstr(derivative,125))),
                        'Original implicit fixed-phase Z derivative not enclosed: '+label+' '+key);count+=1;local+=1
            qreference=[ref(z+shift*h,p.mpf('.5'))[2] for shift in (-2,-1,1,2)]
            derivative=(qreference[0]-8*qreference[1]+8*qreference[2]-qreference[3])/(12*h)
            require(contains(qgot['q'].finite_interval()+allow,c.mpf(p.nstr(q,125))),'Original body/transition/flat q union failed: '+label)
            require(contains(qgot['q_Z'].finite_interval()+allow,c.mpf(p.nstr(derivative,125))),'Original body/transition/flat q_Z union failed: '+label)
            qcount+=2
        if label=='all_three_cutoff_branches':require(seen=={'flat','body','transition'},'Independent fixture must genuinely cross all three cutoff branches')
        regions[label]=dict(original_value_and_implicit_Z_comparisons=local,actual_scalar_cutoff_branches=sorted(seen))
    # The cap is deliberately a cover. Its fixed-phase Z rows need not be
    # zero when q, u or signed r crosses zero inside the original source box.
    require(not got['values']['A_Z'].zero,'Cutoff crossing cannot erase whole-box first derivatives')
    # Very small eta remains a positive formal source in body/transition.
    bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger();scalar=lambda v:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),v,ledger)
    tiny_eta=c.mpf('-1e40');roots=dict(a={current.ZERO:scalar(2),current.DZ:scalar(0)},b={current.ZERO:scalar(0),current.DZ:scalar(0)},
        kappa_minus2={current.ZERO:scalar(0),current.DZ:scalar(0)})
    tiny=current.whole_cutoff_C1(roots,tiny_eta,c.ln(2));require(not tiny['q'].zero,'Tiny positive original cutoff q became zero')
    return dict(passed=True,independent_original_loop_whole_period_value_and_implicit_Z_comparisons=count,
        independent_original_cutoff_q_C0_Z_union_comparisons=qcount,regions=regions,
        genuine_body_transition_flat_and_signed_r_zero_crossings_checked=True,
        large_positive_and_negative_Poisson_parameters_checked=True,unmaterializable_positive_eta_q_retained=True,
        scalar_fixtures_do_not_define_native_field_values=True,finite_difference_reference_allowance='1e-28')


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2,'Actual whole first bridge C1 records required')
    require(not any(saved.get(k) for k in packets.OPEN),'First bridge cannot complete global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed whole first-bridge prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(current.first.NativePhaseFirstJets(current.slow.NativeQSlowJets(current.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))))
        owner=current.NativeFirstBridgeC1Histories(current.transfer.NativeTrueChartC1Transfer(current.history.NativeC1HistoryTransfer(c1)))
        independent=independent_original_loop_checks(owner.ctx);rows=0;regions={}
        for name,old in saved['actual_whole_first_bridge_C1_history_records'].items():
            got=owner.first_bridge(packets.interval(owner.ctx,old['Z_box']),saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual complete first-bridge C1 history changed: '+name)
            require(old['no_interval_skipped_between_true_inlet_and_first_bridge_exit'] and old['actual_active_loop_retained_across_cutoff_crossings'],'Actual first-bridge route cannot skip its active field')
            require(not old['global_inlet_to_Rc_histories_admitted'] and not any(old.get(k) for k in packets.OPEN),'First-bridge bounds cannot admit terminal/global gates')
            require(not got['geometry']['width'].zero and ep(got['geometry']['width'].coefficient)[0]>0,'True microscopic bridge width lost')
            require(old['original_phase1_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history'],'Original P0/P0_Z must remain separate')
            for key in current.RATES:
                require(packets.encode(got['correction']['values'][key].record())==packets.encode(got['contributions'][key].record()),'Known zero initial correction must give its full integrated contribution')
                require(packets.encode(got['correction']['Z_derivatives'][key].record())==packets.encode(got['Z_derivatives'][key].record()),'Actual first-Z incoming correction mismatch')
            for value in [*got['own'].values(),*got['own_Z'].values(),*got['correction']['values'].values(),*got['correction']['Z_derivatives'].values()]:
                require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,'Whole bridge histories must share actual arithmetic basis/ledger')
            row=got['primitives']['record'];require(row['actual_periodic_phase_covered_without_sample_or_selected_inverse'] and row['q_and_phase_crossings_have_genuine_C1_covers'],'Whole function-domain phase/C1 cover required')
            require(row['original_three_branch_cutoff_C1_cover']['flat_Delta_ge_eta_q_and_q_Z_exact_zero'],'Original flat branch must be included exactly')
            require(old['first_bridge_covers_are_conservative_not_small_error_or_terminal_closure'],'Wide covers cannot imply tight closure')
            rows+=20;regions[name]=dict(actual_correction_C0_Z_rows=10,actual_own_history_C0_Z_rows=10,
                whole_remaining_phase_scope='3sc/4 ->1; known prior initial collar completes sc/2 ->1')
            print('Actual complete first-bridge C1 exit checked:',name,flush=True)
        try:owner.first_bridge(N=159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe candidate frequency accepted for analytic primitive cap')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        actual_whole_first_bridge_correction_and_own_C0_Z_rows_checked=rows,native_Z_queries_checked=2,
        original_periodic_primitive_and_implicit_derivative_theorem=current.original_periodic_theorem(),
        independent_original_loop_and_cutoff_C1_checks=independent,regions=regions,
        actual_original_inlet_to_phase1_no_gap_C1_history_covers_checked=True,
        covers_are_conservative_and_do_not_prove_tight_error_or_terminal_closure=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Whole first-bridge C1 PASS:',rows,'actual history rows',flush=True)
    return result


if __name__=='__main__':run()
