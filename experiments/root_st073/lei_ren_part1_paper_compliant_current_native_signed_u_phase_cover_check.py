"""Original scalar reference comparisons and actual broad signed-u replay."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_signed_u_phase_cover as current
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as accepted

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
ep=current.ep;packets=current.packets
require=accepted.require;contains=accepted.contains;original=accepted.original


def independent_branch_checks(c):
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    p=scales.ctx;y,Z=sy.symbols('y Z');h=p.mpf('1e-8');allowance=c.mpf(('-1e-27','1e-27'))
    comparisons=0;branch_count=0;regions={}
    for label,base in (('negative_overlap',sy.Rational(-1,25)),('p2_zero',sy.Integer(0)),('positive_overlap',sy.Rational(1,25))):
        aexpr=sy.Rational(4,5)+sy.Rational(11,100)*y+sy.Rational(7,100)*Z
        bexpr=sy.Rational(1,5)+sy.Rational(7,100)*y-sy.Rational(1,20)*Z
        pexpr=base+sy.Rational(3,100)*y+sy.Rational(1,25)*Z
        eexpr=sy.Rational(13,10)+sy.Rational(3,100)*y+sy.Rational(1,25)*Z
        bases=tuple(c.mpf(0) for _ in range(5));ledger=accepted.qchecks.new_ledger()
        scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
        expressions=dict(a=aexpr,b=bexpr,p2=pexpr,E=eexpr,t0=-bexpr/aexpr,kappa_minus2=aexpr+bexpr*bexpr/aexpr-2)
        roots={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0})))
            for order in current.first.slow.ORDERS} for name,expr in expressions.items()}
        eta_log=c.ln(c.mpf(p.nstr(scales.eta,105)));dstar=c.ln(c.mpf(p.nstr(scales.d_star,105)))
        q0=current.current.q_enclosure(roots['a'][(0,0)],roots['kappa_minus2'][(0,0)],eta_log,c.ln(c.mpf('.7')))
        qjet=current.first.slow.original_q_jet(roots,eta_log,c.ln(c.mpf('.7')),q0)
        require(qjet['rows'] is not None,'Independent original q fixture unresolved')
        source=dict(q=qjet['rows'][(0,0)],roots=roots)
        _,branches,_=current.signed_u_branches(source,dstar)
        names={b['name'] for b in branches}
        require(names==({'central'} if label=='p2_zero' else {'central','negative_tail' if base<0 else 'positive_tail'}),
            'Fixture must exercise the actual signed-u overlap')
        def polynomial(expr):
            terms=sy.Poly(expr,y,Z).terms()
            return lambda yy,zz:sum((p.mpf(int(k.p))/int(k.q)*yy**i*zz**j for (i,j),k in terms),p.mpf(0))
        af,bf,pf,ef=[polynomial(expr) for expr in (aexpr,bexpr,pexpr,eexpr)]
        def reference(yy,zz,phi,guess=None):
            loop=original.GenericShearLoop(scales,a=af(yy,zz),b=bf(yy,zz),p1=5,p2=pf(yy,zz),Utheta=ef(yy,zz))
            start=2*p.pi*phi if guess is None else guess
            psi=p.findroot(lambda x:loop.phase_at_angle(x)-phi,(start,start+p.mpf('.001')),solver='secant',tol=p.mpf('1e-92'))
            return loop.at_angle(psi),psi
        phi=p.mpf('.137');value,guess=reference(p.mpf(0),p.mpf(0),phi)
        references=dict(A=value['A'],B_over_Pstar=value['B'])
        for direction in ('y','Z','phi'):
            rows=[]
            for shift in (-2,-1,1,2):
                yy=shift*h if direction=='y' else p.mpf(0);zz=shift*h if direction=='Z' else p.mpf(0)
                target=phi+shift*h if direction=='phi' else phi
                rows.append(reference(yy,zz,target,guess)[0])
            for name in ('A','B'):
                references[name+'_'+direction+('_over_Pstar' if name=='B' else '')]=(rows[0][name]-8*rows[1][name]+8*rows[2][name]-rows[3][name])/(12*h)
        for branch in branches:
            got=current.branch_first_jets(source,qjet['rows'],dstar,'.137',branch)
            require(got['values'] is not None,'Independent branch first jets unresolved')
            for name,reference_value in references.items():
                require(contains(current.phase.bounded_value(got['values'][name])+allowance,c.mpf(p.nstr(reference_value,105))),
                    'Original scalar reference outside branch '+label+' '+branch['name']+' '+name)
                comparisons+=1
            if label=='p2_zero':
                require(not got['loop'].flat and not got['values']['A_Z'].zero,'p2=0 must retain genuine source derivatives')
            for target in ('0','.5','1'):
                boundary=current.branch_first_jets(source,qjet['rows'],dstar,target,branch)
                require(all(boundary['values'][name].zero for name in current.first.OUTPUTS if '_phi' not in name),
                    'Symmetry primitive/slow rows not exactly zero')
                require(any(not boundary['values'][name].zero for name in current.first.OUTPUTS if '_phi' in name),
                    'Actual fast derivatives must survive symmetry')
            missing=current.branch_first_jets(source,None,dstar,'.137',branch)
            require(missing['values'] is None,'Missing original q rows cannot invent derivatives')
            branch_count+=1
        regions[label]=dict(branches=sorted(names),original_scalar_comparisons=len(references)*len(branches))
    return dict(passed=True,independent_original_scalar_comparisons=comparisons,conditional_branches_checked=branch_count,
        regions=regions,symmetry_branch_queries=3*branch_count,p2_zero_does_not_erase_original_derivatives=True,
        missing_original_q_rows_fail_closed=True,each_overlap_branch_checked_individually=True,
        finite_difference='original GenericShearLoop,100 dps; five-point h=1e-8; allowance1e-27',
        synthetic_reference_fixtures_not_original_field_values=True)


@current.native.inlet.source_precision
def run(original_owner,report=None):
    began=time.monotonic();owner=current.NativeSignedUDensityCover(original_owner);c=owner.ctx
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    require(report[current.GATE] and report['source_family']==owner.family,'Same original signed-u producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed original source prerequisite: '+name)
    require(not any(report.get(key) for key in packets.OPEN),'Local cover cannot complete global construction gates')
    record=report['actual_broad_original_O2_query']
    require(record['status']=='enclosed' and {b['name'] for b in record['conditional_branches']}=={'negative_tail','central','positive_tail'},
        'Full broad original O2 signed cover required')
    for cell in record['spatial_signed_density_branch_cells']:
        require(cell['status']=='enclosed' and all(b['status']=='enclosed' for b in cell['original_conditional_branch_first_jets']),
            'Every nonempty branch must resolve')
        require(cell['original_nonlinear_density_computed_before_branch_union'] and cell['overlapping_branch_ranges_hulled_not_added'],
            'Original nonlinear density must precede overlapping branch hull')
    independent=independent_branch_checks(c)
    # A new original point query still covers the full Z interval and all
    # three u branches; this is not a selected midpoint defining a field.
    query=owner.spatial_query('O2_slope',(-1,1),'.1337',1024)
    require(query['record']['status']=='enclosed','Additional broad original source query unresolved')
    require(all(v.ctx is c and v.ledger is query['source']['q'].ledger
        for cell in query['cells'] for group in ('kernels','Z_derivatives') for v in cell['values'][group].values()),
        'Branch density changed source context or ledger')
    quiet=owner.spatial_query('O3_power',(-1,1),{'original_power_offset':'.537'},2048)
    require(quiet['record']['status']=='enclosed' and all(v.zero for cell in quiet['cells']
        for group in ('kernels','Z_derivatives') for v in cell['values'][group].values()),'Original exact q-flat modulation must vanish')
    require(not quiet['source']['roots']['E'][(0,0)].zero,'Flat modulation must retain original nonzero velocity')
    try:owner.spatial_query('O2_slope',(-1,1),'.1337',159)
    except ValueError:pass
    else:raise AssertionError('N<160 incorrectly accepted')
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        independent_original_overlap_first_jet_checks=independent,
        inherited_actual_broad_original_source=record,additional_original_full_Z_point_query=query['record'],
        original_exact_flat_q_query=quiet['record'],actual_local_C0_Z_integral_rows=10,
        original_source_context_and_ledger_preserved=True,original_integer_N_guard_retained=True,
        working_dependency_hashes_checked=len(report['input_hashes']),
        all_17_chart_original_range_oracle_admitted=False,full_factored_function_graph_evaluator_installed=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Independent original scalar/first-jet comparisons per signed-u overlap branch, exact p2=0 and phase symmetries, original full-Z O2 and exact-flat O3 queries. Broad local integral cover only; global transport/control/closure/recursion remain open.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original signed-u cover focused check PASS:',independent['independent_original_scalar_comparisons'],'independent comparisons',flush=True)
    return result
