"""Focused checks of original shear correlations and the lazy q kernel."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_correlated_shear_q as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as scalar_original

prior=current.prior;packets=current.packets;HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def contains(outer,inner):
    a,b=current.ep(outer);c,d=current.ep(inner)
    return a<=c and d<=b


def q_kernel_check(c):
    bases=tuple(c.mpf(0) for _ in range(5))
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
    # The original scalar constructor implements the paper independently of
    # this new interval kernel; all inputs here are labeled bounded fixtures.
    original=scalar_original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=90)
    cases=0
    for aa,bb,pp1,pp2 in (('.8','.2','5','.2'),('.8','-.2','5','-.2'),('2.019','0','5','0')):
        loop=scalar_original.GenericShearLoop(original,a=aa,b=bb,p1=pp1,p2=pp2,Utheta='1.3')
        eta=c.mpf(str(original.eta));a=c.mpf(aa);b=c.mpf(bb)
        Delta=a+b*b/a-2
        q=current.q_enclosure(scalar(a),scalar(Delta),c.ln(eta),c.ln(c.mpf('.7')))
        actual=q['q'].finite_interval()
        # Scalar finite precision is an independent numerical comparison, not
        # an interval proof. Allow its90-digit arithmetic error explicitly.
        reference=c.mpf(str(loop.q));error=c.mpf('1e-85')
        require(contains(actual+ c.mpf((-current.ep(error)[1],current.ep(error)[1])),reference),
            'q interval disagrees with original scalar source formula')
        if aa=='2.019':require(q['q'].zero and not q['active_body_enclosure_evaluated'],'Flat q must avoid its active body')
        cases+=1
    # Delta=0 is active with q=sqrt(eta/a); eta here cannot be materialized.
    tiny_eta=c.mpf('-1e40')
    tiny=current.q_enclosure(scalar(2),scalar(0),tiny_eta,c.ln(2))
    require(tiny['branch']=='active' and not tiny['q'].zero and tiny['positive_q_lower_retained'],
        'Unmaterializable positive eta must produce positive q, not zero')
    q2=current.square(tiny['q'])
    ratio=q2.positive_divide(prior.ScaledEnclosure(prior.FormalScale(bases,offset=tiny_eta),1,ledger),tiny_eta)
    require(contains(ratio.finite_interval(),c.mpf('.5')),'Delta=0 must satisfy q^2/eta=1/a')
    # On an unresolved branch box the full analytic q range includes both a
    # positive active value and the exact flat value; gamma must stay lazy.
    eta=c.mpf('.01');mixed=current.q_enclosure(scalar(2),scalar((0,'.03')),c.ln(eta),c.ln(2))
    require(mixed['branch']=='requires_source_box_refinement' and current.ep(mixed['q'].coefficient)[0]==0,
        'Mixed cutoff box must include the exact flat branch')
    require(contains(mixed['q'].finite_interval(),c.sqrt(eta/2)),
        'Mixed cutoff box must retain its Delta=0 active value')
    crossed=current.square(scalar((-2,3)))
    require(current.ep(crossed.coefficient)[0]>=0 and contains(crossed.finite_interval(),c.mpf((0,9))),
        'Correlated signed square must retain nonnegative values')
    return dict(passed=True,original_scalar_q_fixture_comparisons=cases,
        positive_tiny_eta_q_not_zeroed=True,exact_flat_active_body_skipped=True,
        mixed_cutoff_hull_checked=True,correlated_signed_square_checked=True,
        fixtures_are_not_current_source_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['live_native_chart_count']==17,'All17 live correlated source queries required')
    require(not any(saved.get(k) for k in packets.OPEN),'Correlated q stage cannot complete global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed correlated q prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge)))
        kernel=q_kernel_check(owner.ctx);regions={};direct=0;changed=0;delta0=0
        for chart,record in saved['actual_correlated_shear_and_q_records'].items():
            p=record['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);v=packets.interval(owner.ctx,p['coordinate_box'])
            live=owner.query(chart,Z,v)
            require(packets.encode(live['record'])==record,'Stored correlated shear/q differs from live original source: '+chart)
            roots=live['roots'];loop=live['loop'];direct+=live['evidence']['correlated_log_source_available']
            require(current.ep(live['q'].coefficient)[0]>=0,'Original q must have a nonnegative enclosure')
            if loop['branch']=='flat':require(live['q'].zero and not loop['active_body_enclosure_evaluated'],'Flat query cannot evaluate active q')
            expected=None
            if chart in ('switch_power','inner_reference','axial_restore','restore_buffer','Rh_reference'):
                expected={(j,k):owner.ctx.mpf('.8') if (j,k)==(0,0) else owner.ctx.mpf(0) for j,k in current.ORDERS}
            elif chart in ('O2_axial','O2_buffer'):
                expected={(j,k):owner.ctx.mpf(2) if (j,k)==(0,0) else owner.ctx.mpf(0) for j,k in current.ORDERS}
            elif chart=='O2_slope':
                sigma=prior.sigma_jets(owner.ctx,v)
                expected={(j,k):owner.ctx.mpf(0) if k else (owner.ctx.mpf('.8')+sigma[0]*owner.ctx.mpf('1.2')
                    if j==0 else sigma[j]*(owner.ctx.mpf('1.2')*math.factorial(j))) for j,k in current.ORDERS}
            if expected:
                for order,value in expected.items():
                    require(contains(roots['a'][order].finite_interval(),value),'Original ordinary log-shear identity not enclosed: '+chart+str(order))
            if chart=='O2_buffer':
                require(roots['kappa_minus2'][(0,0)].zero,'Original zero-excess buffer must retain Delta=0 exactly')
                require(loop['branch']=='active' and loop['positive_q_lower_retained'] and not live['q'].zero,
                    'Original zero-excess buffer has tiny positive q')
                delta0+=1
            changed+=loop['branch']!=record['previous_independent_root_branch']
            regions[chart]=dict(branch=loop['branch'],previous_branch=record['previous_independent_root_branch'],
                correlated_source=live['evidence']['correlated_log_source_available'],positive_q_lower=loop['positive_q_lower_retained'])
            print('Live original correlated q checked:',chart,flush=True)
    require(direct==16 and delta0==1,'Expected16 log-shear chart sources and exact-zero-excess buffer')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        live_q_enclosure_queries_checked=len(regions),original_correlated_log_shear_sources_checked=direct,
        branch_classifications_sharpened=changed,zero_excess_tiny_positive_q_source_checked=delta0,
        focused_original_q_kernel_check=kernel,regions=regions,
        source_midpoints_or_caps_used_as_field_values=False,actual_phase_inverse_or_integral_admission=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Live correlated source shear/q PASS:',len(regions),'queries;',direct,'original log-shear sources',flush=True)
    return result


if __name__=='__main__':run()
