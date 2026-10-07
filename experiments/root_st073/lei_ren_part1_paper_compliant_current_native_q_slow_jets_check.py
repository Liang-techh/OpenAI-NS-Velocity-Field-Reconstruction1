"""Original scalar derivative references and native q mixed-jet replay."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_q_slow_jets as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

prior=current.prior;packets=current.packets;HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
ep=current.ep;require=checks.require;contains=checks.contains;original=checks.original


def new_ledger():
    return dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)


def scalar_jet_checks(c):
    p=mp.mp.clone();p.dps=150;y,Z=sy.symbols('y Z');eta_s=sy.Rational(1,100)
    shape=sy.Rational(11,100)*y+sy.Rational(7,100)*Z+sy.Rational(31,1000)*y*y+sy.Rational(1,20)*y*Z+sy.Rational(29,1000)*y*y*Z
    aa=sy.Rational(4,5)+shape
    bb=sy.Rational(1,5)+sy.Rational(7,100)*y-sy.Rational(1,20)*Z+sy.Rational(1,50)*y*y*Z
    cases=[('plateau_positive_b',aa,bb,None),('plateau_negative_b',aa,-bb,None),
        ('transition',2+eta_s*(sy.Rational(37,100)+shape),sy.Integer(0),None),
        ('flat',2+eta_s*(sy.Rational(13,10)+shape),sy.Integer(0),None),
        ('Delta_zero',2+eta_s*shape,sy.Integer(0),'zero')]
    count=0;branches={}
    for label,aexpr,bexpr,boundary in cases:
        bases=tuple(c.mpf(0) for _ in range(5));ledger=new_ledger()
        scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
        eta=prior.ScaledEnclosure(prior.FormalScale(bases,offset=c.ln(c.mpf('.01'))),1,ledger)
        Delta=aexpr+bexpr*bexpr/aexpr-2
        jets={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0}))) for order in current.ORDERS}
              for name,expr in (('a',aexpr),('kappa_minus2',Delta))}
        if boundary=='zero':jets['kappa_minus2'][current.ZERO]=scalar(0)
        loop=current.current.q_enclosure(jets['a'][current.ZERO],jets['kappa_minus2'][current.ZERO],c.ln(c.mpf('.01')),c.ln(c.mpf('.7')))
        got=current.original_q_jet(jets,c.ln(c.mpf('.01')),c.ln(c.mpf('.7')),loop)
        require(got['status']=='enclosed','Bounded original q derivative fixture unresolved: '+label)
        def polynomial(expr):
            terms=sy.Poly(expr,y,Z).terms()
            return lambda yv,zv:sum((p.mpf(int(coefficient.p))/int(coefficient.q)*yv**i*zv**k
                for (i,k),coefficient in terms),p.mpf(0))
        af=polynomial(aexpr);bf=polynomial(bexpr)
        def reference(yv,zv):
            a,b=af(yv,zv),bf(yv,zv);excess=a+b*b/a-2;eta=p.mpf('.01')
            if excess>=eta:return p.mpf(0)
            return original.flat_step(p,1-excess/eta)*p.sqrt((2*eta-excess)/(2*a))
        for order,value in got['rows'].items():
            ref=p.diff(reference,(p.mpf(0),p.mpf(0)),order)
            covered=value.finite_interval()+c.mpf(('-1e-110','1e-110'))
            require(contains(covered,c.mpf(p.nstr(ref,145))),'Independent original q derivative not enclosed: '+label+str(order))
            count+=1
        if label=='flat':
            require(not got['active_body_enclosure_evaluated'] and all(v.zero for v in got['rows'].values()),'Flat cutoff derivatives must stay lazy/exact zero')
        branches[label]=dict(branch=got['branch'],cutoff=got['cutoff_scope'])
    # Exact shared-source Delta=eta at a point. A singleton log keeps this
    # identity visible to the existing formal-scale cancellation, unlike an
    # independent enclosing log(.01) interval. eta=exp(-1)<1/2 is admissible.
    bases=tuple(c.mpf(0) for _ in range(5));ledger=new_ledger()
    scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
    eta_log=c.mpf(-1);eta=prior.ScaledEnclosure(prior.FormalScale(bases,offset=eta_log),1,ledger)
    Delta={order:eta*(1 if order==current.ZERO else c.mpf('.1')) for order in current.ORDERS}
    a={order:Delta[order]+(2 if order==current.ZERO else 0) for order in current.ORDERS}
    loop=current.current.q_enclosure(a[current.ZERO],Delta[current.ZERO],eta_log,c.ln(2))
    got=current.original_q_jet(dict(a=a,kappa_minus2=Delta),eta_log,c.ln(2),loop)
    require(got['status']=='enclosed' and not got['active_body_enclosure_evaluated'] and all(v.zero for v in got['rows'].values()),'Exact shared-source cutoff boundary needs all six flat derivatives zero')
    count+=6;branches['Delta_eta']=dict(branch='flat',cutoff=got['cutoff_scope'],exact_shared_eta_source=True)
    # Positive unmaterializable q at Delta=0 has genuine zero slow derivatives
    # for constant a/Delta, rather than q itself being replaced by zero.
    bases=tuple(c.mpf(0) for _ in range(5));ledger=new_ledger();scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
    roots={'a':{o:scalar(2 if o==current.ZERO else 0) for o in current.ORDERS},
           'kappa_minus2':{o:scalar(0) for o in current.ORDERS}}
    eta_log=c.mpf('-1e40');loop=current.current.q_enclosure(roots['a'][current.ZERO],roots['kappa_minus2'][current.ZERO],eta_log,c.ln(2))
    got=current.original_q_jet(roots,eta_log,c.ln(2),loop)
    require(not got['rows'][current.ZERO].zero and all(v.zero for o,v in got['rows'].items() if o!=current.ZERO),'Tiny constant positive q slow jets lost')
    eta=prior.ScaledEnclosure(prior.FormalScale(bases,offset=eta_log),1,ledger)
    ratio=current.current.square(got['rows'][current.ZERO]).positive_divide(eta,eta_log).finite_interval()
    require(contains(ratio,c.mpf('.5')),'Tiny Delta=0 q must retain q^2/eta=1/a')
    # Regression for original O2 a=2: tiny positive b^2/a must not disappear
    # in (2+tiny)-2. This directly controls the true cutoff branch when eta
    # is even smaller, rather than only improving an interval's display.
    tiny_b=prior.ScaledEnclosure(prior.FormalScale(bases,offset=c.mpf('-1e40')),1,ledger)
    b={order:tiny_b if order==current.ZERO else scalar(0) for order in current.ORDERS}
    excess=current.correlated_excess_rows(dict(a=roots['a'],b=b),c.ln(2))
    require(ep(excess[current.ZERO].coefficient)[0]>0,'Same-source a=2 positive b^2/a lower bound lost')
    refined_loop=current.current.q_enclosure(roots['a'][current.ZERO],excess[current.ZERO],c.mpf('-3e40'),c.ln(2))
    require(refined_loop['branch']=='flat' and not refined_loop['active_body_enclosure_evaluated'],'Tiny positive excess exceeding smaller eta must select exact flat source branch')
    roots['kappa_minus2'][current.ZERO]=scalar((0,'.03'))
    mixed=current.current.q_enclosure(roots['a'][current.ZERO],roots['kappa_minus2'][current.ZERO],c.ln(c.mpf('.01')),c.ln(2))
    result=current.original_q_jet(roots,c.ln(c.mpf('.01')),c.ln(2),mixed)
    require(result['rows'] is None and result['status']=='requires_source_branch_subdivision','Mixed source box cannot invent q derivatives')
    return dict(passed=True,independent_original_scalar_mixed_q_derivative_comparisons=count,
        plateau_signed_b_transition_flat_and_both_cutoff_boundaries_checked=branches,
        reference_roundoff_allowed='1e-110',original_scalar_cutoff_used=True,
        unmaterializable_Delta_zero_positive_q_and_exact_constant_derivatives_retained=True,
        exact_O2_constant_cancelled_before_positive_tiny_quotient=True,
        tiny_positive_excess_controls_actual_smaller_eta_flat_branch=True,
        unresolved_branch_box_returns_no_invented_derivatives=True,fixtures_are_not_native_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_query_chart_count']==17,'17 native original source q-jet queries required')
    require(not any(saved.get(k) for k in packets.OPEN),'Slow q derivatives cannot complete global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed q slow-jet prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeQSlowJets(current.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))
        independent=scalar_jet_checks(owner.ctx);regions={};enclosed=0;flat=0;unresolved=0;row_count=0
        records=dict(saved['actual_native_q_slow_jet_records']);records['actual_integral_cell']=saved['actual_integral_cell_q_slow_jet_record']
        records.update({'Z_subdivision_'+key:value for key,value in saved['actual_source_Z_subdivision_q_slow_jet_records'].items()})
        for name,old in records.items():
            provenance=old['source_provenance'];Z=packets.interval(owner.ctx,provenance['Z_box']);coordinate=packets.interval(owner.ctx,provenance['coordinate_box'])
            got=owner.query(old['chart'],Z,coordinate)
            require(packets.encode(got['record'])==old,'Native q slow derivative source changed: '+name)
            require(got['source']['record']['original_q_enclosure']['branch']==got['record']['branch'],'Returned source record has stale pre-refinement branch')
            require(packets.encode(got['source']['record']['original_q_enclosure']['q'])==packets.encode(got['source']['q'].record()),'Returned source record has stale pre-refinement q')
            require(not any(got['record'].get(k) for k in packets.OPEN),'Native q slow-jet query cannot admit global stages')
            if got['rows'] is None:
                require(got['record']['status'].startswith('requires_source_'),'Missing rows require explicit source subdivision')
                unresolved+=1
            else:
                require(set(got['rows'])==set(current.ORDERS),'All six ordinary q rows required')
                require(got['record']['native_width_and_Pstar_conversion_not_reapplied'],'Original ordinary source conversions must occur once')
                enclosed+=1;row_count+=len(got['rows'])
                if got['record']['branch']=='flat':
                    require(not got['record']['active_body_enclosure_evaluated'] and all(v.zero for v in got['rows'].values()),'Flat source q jets must all be lazy zero')
                    flat+=1
                if old['chart']=='O2_buffer':
                    require(not got['rows'][current.ZERO].zero and all(v.zero for o,v in got['rows'].items() if o!=current.ZERO),'Actual constant buffer q is tiny positive with zero slow derivatives')
            regions[name]=dict(chart=old['chart'],status=got['record']['status'],branch=got['record']['branch'],ordinary_rows=0 if got['rows'] is None else len(got['rows']))
            print('Native original q slow jets checked:',name,flush=True)
    require(enclosed==18 and unresolved==2 and flat==4 and row_count==108,'Expected15 original enclosed charts plus integral cell and2 actual Z subdivisions;4 flat and2 unresolved')
    require(regions['actual_integral_cell']['status']=='enclosed','Actual whole integral radial/Z cell q jets required')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_source_queries_checked=20,
        native_original_chart_queries_checked=17,enclosed_original_chart_queries=15,actual_integral_cell_q_jet_query_checked=True,
        original_q_ordinary_derivative_enclosures_checked=row_count,exact_lazy_flat_source_queries=flat,
        actual_source_Z_subdivision_queries_checked=2,
        unresolved_source_queries_retained=unresolved,regions=regions,independent_original_q_derivative_checks=independent,
        phase_A_B_slow_jets_or_global_C1_histories_or_common_N_admission=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original q slow jets PASS:',row_count,'ordinary derivative enclosures',flush=True)
    return result


if __name__=='__main__':run()
