"""Independent mixed q squared derivatives and original phase comparisons."""
import ast
import inspect
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_collected_q2_transport as current
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as reference

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets=current.packets;ep=current.ep;require=reference.require;contains=reference.contains


@current.native.inlet.source_precision
def independent_q2_checks(c):
    p=reference.mp.mp.clone();p.dps=110;y,Z=sy.symbols('y Z');eta=p.mpf('.1');eta_log=c.ln(c.mpf('.1'))
    allowance=c.mpf(('-1e-70','1e-70'));comparisons=0;regions={}
    for label,base in (('negative','-.9'),('zero','0'),('transition','.037'),('eta_endpoint','.1'),('flat','.3')):
        bases=tuple(c.mpf(0) for _ in range(5));ledger=reference.qchecks.new_ledger()
        scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
        aa=sy.Rational(4,5)+sy.Rational(11,100)*y+sy.Rational(7,100)*Z+sy.Rational(1,50)*y*y*Z
        dd=sy.Rational(base)+sy.Rational(3,100)*y+sy.Rational(1,25)*Z+sy.Rational(1,20)*y*Z+sy.Rational(1,70)*y*y*Z
        roots={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0})))
            for order in current.ORDERS} for name,expr in (('a',aa),('kappa_minus2',dd))}
        def original_q_squared(yy,zz):
            a=p.mpf('.8')+p.mpf('.11')*yy+p.mpf('.07')*zz+yy*yy*zz/50
            D=p.mpf(base)+p.mpf('.03')*yy+p.mpf('.04')*zz+p.mpf('.05')*yy*zz+yy*yy*zz/70
            if D>=eta:return p.mpf(0)
            q=reference.original.flat_step(p,1-D/eta)*p.sqrt((2*eta-D)/(2*a))
            return q*q
        expected={order:p.diff(original_q_squared,(p.mpf(0),p.mpf(0)),order) for order in current.ORDERS}
        branches,_=current.cutoff.conditional_cutoff_branches(roots['kappa_minus2'][current.ZERO],eta_log)
        for branch in branches:
            got=current.conditional_q2_jet(roots,eta_log,c.ln(c.mpf('.7')),branch)
            for order,value in got.items():
                require(contains(current.cover.phase.bounded_value(value)+allowance,c.mpf(p.nstr(expected[order],115))),
                    'Original q squared derivative outside '+label+' '+branch['name']+' '+str(order));comparisons+=1
            if label=='zero':require(not got[(0,1)].zero,'Delta=0 cannot erase original square/source derivative')
            if branch['name']=='flat':require(all(value.zero for value in got.values()),'Original flat q squared jet changed')
        regions[label]=dict(branches=[branch['name'] for branch in branches],comparisons=len(branches)*len(current.ORDERS))
    # A deliberately huge eta logarithm exercises the shared-factor algebra,
    # separately from original field data. q squared first derivatives must
    # not acquire an artificial exp(-eta_log) after eta/eta cancellation.
    bases=tuple(c.mpf(0) for _ in range(5));ledger=reference.qchecks.new_ledger()
    scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
    huge=c.mpf('-1e100');eta_source=current.prior.ScaledEnclosure(current.prior.FormalScale(bases,offset=huge),1,ledger)
    roots={name:{order:scalar(1 if name=='a' and order==current.ZERO else 0) for order in current.ORDERS}
        for name in ('a','kappa_minus2')}
    roots['kappa_minus2'][(0,1)]=scalar(1)
    branch=dict(name='transition',condition='0<=Delta<=eta',Delta=eta_source*c.mpf('.37'),theta=c.mpf('.37'))
    rows=current.conditional_q2_jet(roots,huge,c.mpf(0),branch)
    require(ep(rows[(0,1)].record()['log_absolute_upper'])[1]<100,
        'First-order eta/eta cancellation still evaluated as independent huge logs')
    require(not rows[current.ZERO].zero,'Tiny positive original q squared cannot become exact zero')
    return dict(passed=True,independent_original_q_squared_ordinary_derivative_comparisons=comparisons,regions=regions,
        independent_original_scalar_q_squared_used=True,ordinary_orders=[list(order) for order in current.ORDERS],
        huge_eta_log_first_derivative_cancellation_and_nonzero_C0_passed=True,
        derivative_method='mpmath multivariate diff110 dps; allowance1e-70',synthetic_fixtures_not_original_field_values=True)


@current.native.inlet.source_precision
def independent_collected_phase_checks(c):
    # Reuse the independently defined scalar roots, original GenericShearLoop,
    # numerical inverse and centered differences. Only the tested consumer is
    # adapted; the reference functions and all comparison rules are unchanged.
    scales=reference.original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    eta_log=c.ln(c.mpf(scales.ctx.nstr(scales.eta,105)))
    def collected_fixture(source,qrows,dstar,phi):
        if qrows is None:return reference.current.conditioned_first_jets(source,None,dstar,phi)
        cutoffs,_=current.cutoff.conditional_cutoff_branches(source['roots']['kappa_minus2'][current.ZERO],eta_log)
        u,branches,empty=current.cover.signed_u_branches(source,dstar);outputs=[]
        for cutoff_branch in cutoffs:
            rows=current.conditional_q2_jet(source['roots'],eta_log,c.ln(c.mpf('.7')),cutoff_branch)
            for branch in branches:
                got=current.compile_collected_first_jets(source,qrows,rows,dstar,phi,branch)
                require(got['values'] is not None,'Collected scalar phase fixture unresolved');outputs.append(got)
        union=current.cover.density.local.same_source_union;got=dict(outputs[0])
        got['values']={key:union([row['values'][key] for row in outputs]) for key in reference.current.OUTPUTS}
        return got
    tree=ast.parse(inspect.getsource(reference.scalar_first_checks))
    changes=[('current.conditioned_first_jets(source,qjet[\'rows\'],dstar,\'.137\')',"collected_fixture(source,qjet['rows'],dstar,'.137')"),
        ('current.conditioned_first_jets(source,qjet[\'rows\'],dstar,phi)',"collected_fixture(source,qjet['rows'],dstar,phi)"),
        ("current.conditioned_first_jets(source,None,dstar,'.137')","collected_fixture(source,None,dstar,'.137')")]
    tree,counts=current.density.replace_expressions(tree,changes)
    scope=dict(vars(reference));scope['collected_fixture']=collected_fixture
    exec(compile(tree,'<independent-original-scalar-checks-of-collected-q-squared-consumer>','exec'),scope)
    got=scope['scalar_first_checks'](c);got['test_consumer_AST_replacements']=counts
    got['independent_original_references_and_allowances_unchanged']=True
    return got


@current.native.inlet.source_precision
def run(report,owner,live):
    began=time.monotonic()
    require(report[current.GATE] and report['source_family']==owner.family,'Same original collected q squared producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed q squared prerequisite: '+name)
    require(report['enclosed_original_whole_cells']==24 and report['full24_original_C1_integral_range_transport_enclosed'],
        'Original24 continuous cells and all first-Z history rows required')
    require(len(live['cells'])==24 and live['history'] is not None,'Actual original route objects required')
    first_jet_rows=0;branch_jet_rows=0;integral_rows=0
    for cell in live['cells']:
        for group in ('values','Z_derivatives'):
            for value in cell[group].values():
                require(value.ctx is owner.ctx and value.ledger is owner.coordinates.ledger,'Original context or ledger changed');integral_rows+=1
        if cell['frame'] is None:continue
        source=cell['frame'].record['actual_original_spatial_source']
        qsource=source['original_q_slow_jet_source']
        require(qsource['original_q_and_six_ordinary_q_rows_unchanged'],'Original linear q rows changed')
        for row in qsource['original_branch_local_q_squared_jets']:
            require(len(row['ordinary_q_squared_rows'])==6 and row['original_a_positive_denominator_only'],
                'Six original q squared mixed rows and original positivity required');branch_jet_rows+=6
        for branch in source['original_cutoff_branch_density_queries']:
            spatial=branch['original_conditional_spatial_density']
            for phase in spatial['spatial_signed_density_branch_cells']:
                for item in phase['original_conditional_branch_first_jets']:
                    adapter=item['original_direct_q_squared_first_jet_adapter']
                    require(adapter['AST_replacement_counts']==[1,1,1,5] and
                        adapter['original_linear_q_q_Z_u_h_and_implicit_phase_rules_unchanged'],
                        'Original five-site quadratic adapter or linear rows changed');first_jet_rows+=1
    for row in report['original_serial_cells']:
        if row['chart']=='O3_power':
            require(row['incoming_C0']['p']==row['outgoing_C0']['p'] and row['incoming_Z']['p']==row['outgoing_Z']['p'],
                'Original quiet pressure memory changed')
    require(not any(report.get(key) for key in packets.OPEN) and not report['useful_repair_contraction_or_actual_controls_established'],
        'Source derivative coordinates cannot complete controls/global N/closure/recursion')
    mixed=independent_q2_checks(owner.ctx);phase=independent_collected_phase_checks(owner.ctx)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=2048,
        original_continuous_route_cells=24,original_live_integral_C0_Z_rows=integral_rows,
        original_branch_local_q_squared_ordinary_rows_checked=branch_jet_rows,
        original_first_jet_five_site_adapters_checked=first_jet_rows,
        independent_original_q_squared_checks=mixed,independent_original_collected_phase_checks=phase,
        strict_target_absolute_upper_reductions=report['strict_target_absolute_upper_reductions'],
        original_linear_q_and_q_Z_and_quiet_pressure_memory_retained=True,
        genuine_higher_order_inverse_eta_factors_retained=True,completed_actual_controls_global_N_closure_or_recursion=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Independent mixed original q squared derivatives, unchanged scalar phase/primitive references and all24 original continuous transport cells. Original linear q derivatives remain; no actual controls/global N/closure/recursion admission.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original collected q squared focused check PASS;',first_jet_rows,'original five-site adapters',flush=True)
    return result
