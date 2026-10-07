"""Independent original scalar first derivatives and same-source native replay."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_phase_first_jets as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks
import lei_ren_part1_paper_compliant_current_native_q_slow_jets_check as qchecks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
prior=current.prior;slow=current.slow;packets=current.packets;ep=current.ep
require=checks.require;contains=checks.contains;original=checks.original


def scalar_first_checks(c):
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    p=scales.ctx;y,Z=sy.symbols('y Z');shape=sy.Rational(11,100)*y+sy.Rational(7,100)*Z+sy.Rational(1,20)*y*Z
    eta=p.nstr(scales.eta,105);eta_s=sy.Rational(eta)
    aa=sy.Rational(4,5)+shape;bb=sy.Rational(1,5)+sy.Rational(7,100)*y-sy.Rational(1,20)*Z
    ee=sy.Rational(13,10)+sy.Rational(3,100)*y+sy.Rational(1,25)*Z
    cases=[('positive_Mobius',aa,bb,sy.Rational(1,5)+shape),
        ('negative_Mobius',aa,-bb,-sy.Rational(1,5)+shape),
        ('small_nonzero',aa,bb,sy.Rational(1,50)+shape),
        ('p2_zero_crossing',aa,bb,shape),
        ('transition',2+eta_s*(sy.Rational(37,100)+shape),sy.Integer(0),sy.Rational(1,5)+shape),
        ('flat',sy.Rational(2019,1000)+shape,sy.Integer(0),shape)]
    comparisons=0;regions={};h=p.mpf('1e-8');allowance=c.mpf(('-1e-27','1e-27'))
    for label,aexpr,bexpr,p2expr in cases:
        bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger()
        scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
        expressions=dict(a=aexpr,b=bexpr,p2=p2expr,E=ee,t0=-bexpr/aexpr,kappa_minus2=aexpr+bexpr*bexpr/aexpr-2)
        roots={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0}))) for order in slow.ORDERS}
            for name,expr in expressions.items()}
        eta_log=c.ln(c.mpf(eta));dstar=c.ln(c.mpf(p.nstr(scales.d_star,105)))
        loop=current.current.q_enclosure(roots['a'][(0,0)],roots['kappa_minus2'][(0,0)],eta_log,c.ln(c.mpf('.7')))
        qjet=slow.original_q_jet(roots,eta_log,c.ln(c.mpf('.7')),loop)
        require(qjet['rows'] is not None,'Independent first-jet fixture q unresolved: '+label)
        source=dict(q=qjet['rows'][(0,0)],roots=roots)
        def polynomial(expr):
            terms=sy.Poly(expr,y,Z).terms()
            return lambda yy,zz:sum((p.mpf(int(coefficient.p))/int(coefficient.q)*yy**i*zz**k
                for (i,k),coefficient in terms),p.mpf(0))
        af,bf,pf,ef=[polynomial(expr) for expr in (aexpr,bexpr,p2expr,ee)]
        def reference(yy,zz,phi,guess=None):
            ref=original.GenericShearLoop(scales,a=af(yy,zz),b=bf(yy,zz),p1=5,p2=pf(yy,zz),Utheta=ef(yy,zz))
            if ref.q==0 or phi in (0,p.mpf('.5'),1):psi=2*p.pi*phi
            else:
                start=2*p.pi*phi if guess is None else guess
                try:psi=p.findroot(lambda x:ref.phase_at_angle(x)-phi,(start,start+p.mpf('.001')),solver='secant',tol=p.mpf('1e-92'))
                except (ValueError,ArithmeticError):psi=ref.angle_at_phase(phi)
            result=ref.at_angle(psi);return result,psi
        phi=p.mpf('.137');ref,guess=reference(p.mpf(0),p.mpf(0),phi)
        got=current.conditioned_first_jets(source,qjet['rows'],dstar,'.137')
        require(got['record']['status']=='enclosed','Independent first-jet fixture inverse unresolved: '+label)
        for name,rname in (('A','A'),('B_over_Pstar','B')):
            require(contains(current.phase.bounded_value(got['values'][name])+allowance,c.mpf(p.nstr(ref[rname],105))),
                'Independent original C0 primitive not enclosed: '+label+' '+name)
            comparisons+=1
        for direction in ('y','Z','phi'):
            rows=[]
            for shift in (-2,-1,1,2):
                yy=shift*h if direction=='y' else p.mpf(0);zz=shift*h if direction=='Z' else p.mpf(0)
                target=phi+shift*h if direction=='phi' else phi
                rows.append(reference(yy,zz,target,guess)[0])
            for name,rname in (('A','A'),('B','B')):
                numeric=(rows[0][rname]-8*rows[1][rname]+8*rows[2][rname]-rows[3][rname])/(12*h)
                key=name+'_'+direction+('_over_Pstar' if name=='B' else '')
                require(contains(current.phase.bounded_value(got['values'][key])+allowance,c.mpf(p.nstr(numeric,105))),
                    'Independent original first derivative not enclosed: '+label+' '+key)
                comparisons+=1
        if label=='p2_zero_crossing':
            require(got['record']['geometry']=='small_r_series' and not got['values']['A_Z'].zero,
                'r=0 point cannot erase a genuine source derivative across p2=0')
        if label=='flat':require(all(v.zero for v in got['values'].values()),'Verified flat q jets imply exact flat primitive jets')
        regions[label]=dict(geometry=got['record']['geometry'],selected_chart=got['record'].get('selected_chart'),comparisons=8)
    # Special phases have zero value/slow derivatives, but the fast derivative
    # is computed from the original loop, rather than zeroed by symmetry.
    aexpr,bexpr,p2expr=cases[0][1:]
    bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger();scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(bases),v,ledger)
    roots={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0}))) for order in slow.ORDERS}
        for name,expr in dict(a=aexpr,b=bexpr,p2=p2expr,E=ee,t0=-bexpr/aexpr,kappa_minus2=aexpr+bexpr*bexpr/aexpr-2).items()}
    eta_log=c.ln(c.mpf(eta));dstar=c.ln(c.mpf(p.nstr(scales.d_star,105)))
    loop=current.current.q_enclosure(roots['a'][(0,0)],roots['kappa_minus2'][(0,0)],eta_log,c.ln(c.mpf('.7')))
    qjet=slow.original_q_jet(roots,eta_log,c.ln(c.mpf('.7')),loop);source=dict(q=qjet['rows'][(0,0)],roots=roots)
    scalar_loop=original.GenericShearLoop(scales,a='.8',b='.2',p1=5,p2='.2',Utheta='1.3')
    for phi in ('0','.5','1'):
        got=current.conditioned_first_jets(source,qjet['rows'],dstar,phi);psi=2*p.pi*p.mpf(phi)
        t=scalar_loop.direction(psi);ratio=(1+scalar_loop.t0**2+2*scalar_loop.q**2)/(1+t*t)
        references=dict(A_phi=scalar_loop.a/2*(1-ratio),B_phi_over_Pstar=scalar_loop.Utheta*scalar_loop.a/2*(scalar_loop.t0-t*ratio))
        require(all(got['values'][name].zero for name in current.OUTPUTS if '_phi' not in name),'Symmetry primitive/slow derivative not zero')
        for name,value in references.items():
            require(contains(current.phase.bounded_value(got['values'][name])+allowance,c.mpf(p.nstr(value,105))),
                'Original fast derivative lost at symmetry phase '+phi+' '+name)
            comparisons+=1
    unresolved=current.conditioned_first_jets(source,None,dstar,'.137')
    require(unresolved['values'] is None,'Missing source q rows cannot invent first primitive derivatives')
    return dict(passed=True,independent_original_scalar_C0_and_first_derivative_comparisons=comparisons,
        regions=regions,finite_difference_method='five-point centered, h=1e-8, numerical comparison allowance1e-27; scalar root precision100 dps',
        independent_original_phase_and_primitive_functions_used=True,
        p2_zero_crossing_derivatives_not_erased=True,symmetry_fast_derivatives_retained=True,
        unresolved_q_rows_do_not_invent_derivatives=True,fixtures_are_not_native_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and len(saved['native_fixed_phase_first_jet_records'])==5,'Five original native first-jet boxes required')
    require(not any(saved.get(k) for k in packets.OPEN),'First jets cannot complete global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed first-jet prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativePhaseFirstJets(slow.NativeQSlowJets(current.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge)))))
        independent=scalar_first_checks(owner.ctx);queries=0;rows=0;chain_rows=0;regions={}
        def inspect(got,old,name):
            nonlocal queries,rows
            require(packets.encode(got['record'])==old,'Original native first-jet replay changed: '+name)
            require(got['record']['status']=='enclosed' and set(got['values'])==set(current.OUTPUTS),'Eight original primitive first-jet rows required')
            require(not any(got['record'].get(k) for k in packets.OPEN),'Native first jets cannot admit global gates')
            if got['loop'].flat:require(all(v.zero for v in got['values'].values()),'Native flat derivatives must be exact zero')
            queries+=1;rows+=8
        for chart,old in saved['native_fixed_phase_first_jet_records'].items():
            prov=old['source_provenance'];got=owner.query(chart,packets.interval(owner.ctx,prov['Z_box']),packets.interval(owner.ctx,prov['coordinate_box']),'.137')
            inspect(got,old,chart);regions[chart]=dict(geometry=got['record']['geometry'],first_rows=8)
            if chart=='O2_buffer':require(not got['values']['B_over_Pstar'].zero,'Tiny positive native buffer cannot be erased')
            print('Original native A/B first jets checked:',chart,flush=True)
        for phi,old in saved['native_symmetry_phase_first_jet_records'].items():
            got=owner.query('O2_slope',('.5','.5'),'.1337',phi);inspect(got,old,'symmetry '+phi)
            require(all(got['values'][key].zero for key in current.OUTPUTS if '_phi' not in key),'Original native symmetry slow rows must be zero')
        spatial_records=dict(saved['native_actual_spatial_first_jet_records']);spatial_records['actual_whole_integral_cell']=saved['actual_whole_integral_cell_spatial_first_jet_record']
        for name,old in spatial_records.items():
            geometry=old['actual_original_radius_phase'];prov=old['source_provenance']
            coordinate={'original_power_offset':'.537'} if old['chart']=='O3_power' else packets.interval(owner.ctx,prov['coordinate_box'])
            exact=geometry['declared_exact_coordinate']
            if exact is not None:coordinate=str(exact['numerator'])+'/'+str(exact['denominator'])
            got=owner.spatial_query(old['chart'],packets.interval(owner.ctx,prov['Z_box']),coordinate,saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Original native spatial first-jet replay changed: '+name)
            require(geometry['phase_independent_of_Z'] and not any(old.get(k) for k in packets.OPEN),'Actual radius chain requires Z-independent original phase and open global gates')
            for cell in got['cells']:
                require(cell['values'] is not None and set(cell['values'])==set(current.OUTPUTS),'Original spatial first jets unresolved')
                require(set(cell['chain'])=={'A_y_total','B_y_total_over_Pstar','A_Z_total','B_Z_total_over_Pstar'},'Four total spatial first derivative rows required')
                require(cell['chain']['A_Z_total'] is cell['values']['A_Z'] and cell['chain']['B_Z_total_over_Pstar'] is cell['values']['B_Z_over_Pstar'],'Z-independent actual phase must preserve original slow Z rows')
                queries+=1;rows+=8;chain_rows+=4
            print('Original native actual spatial first jets checked:',name,flush=True)
    require(queries==14 and rows==112 and chain_rows==24,'Expected8 fixed/symmetry plus6 actual spatial query cells')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_first_jet_query_cells_checked=queries,
        original_A_B_value_and_first_derivative_rows_checked=rows,actual_total_spatial_first_derivative_rows_checked=chain_rows,
        actual_whole_integral_radial_Z_cell_first_jets_checked=True,independent_original_first_derivative_checks=independent,
        exact_symmetry_phase_queries_checked=3,regions=regions,density_C1_integrals_or_global_C1_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original native A/B first jets PASS:',rows,'primitive rows;',chain_rows,'actual spatial chain rows',flush=True)
    return result


if __name__=='__main__':run()
