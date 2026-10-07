"""Original parameter-bound references and actual complete O2 C1 route."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_O2_C1_histories as current
import lei_ren_part1_paper_compliant_current_native_serial_C1_cell as serial
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks
import lei_ren_part1_paper_compliant_current_native_q_slow_jets_check as qchecks
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories_check as previous_checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets;ep=current.ep
require=checks.require;contains=checks.contains;original=checks.original
same=previous_checks.same;overlaps=previous_checks.overlaps

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
        roots={name:{order:scalar(interval(sy.diff(expr,Z,order[1]))) for order in (serial.ZERO,serial.DZ)} for name,expr in expressions.items()}
        eta_log=c.ln(c.mpf(eta_text));dstar_log=c.ln(c.mpf(p.nstr(scales.d_star,125)))
        got=serial.whole_period_C1(roots,eta_log,c.ln(c.mpf('.7')),dstar_log)
        qgot=serial.p.whole_cutoff_C1(roots,eta_log,c.ln(c.mpf('.7')))
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
    tiny_eta=c.mpf('-1e40');roots=dict(a={serial.ZERO:scalar(2),serial.DZ:scalar(0)},b={serial.ZERO:scalar(0),serial.DZ:scalar(0)},
        kappa_minus2={serial.ZERO:scalar(0),serial.DZ:scalar(0)})
    tiny=serial.p.whole_cutoff_C1(roots,tiny_eta,c.ln(2));require(not tiny['q'].zero,'Tiny positive original cutoff q became zero')
    return dict(passed=True,independent_original_loop_whole_period_value_and_implicit_Z_comparisons=count,
        independent_original_cutoff_q_C0_Z_union_comparisons=qcount,regions=regions,
        genuine_body_transition_flat_and_signed_r_zero_crossings_checked=True,
        large_positive_and_negative_Poisson_parameters_checked=True,unmaterializable_positive_eta_q_retained=True,
        scalar_fixtures_do_not_define_native_field_values=True,finite_difference_reference_allowance='1e-28')



def exact_flat_primitive_and_density_check(c):
    bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger()
    scalar=lambda v:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),v,ledger)
    roots={name:{serial.ZERO:scalar(value),serial.DZ:scalar(jet)} for name,value,jet in
        (('a',3,0),('b',0,0),('p2',10,2),('E',2,1),('kappa_minus2',1,0))}
    got=serial.whole_period_C1(roots,c.ln(c.mpf('.01')),c.ln(3),c.ln(c.mpf('.1')))
    require(all(v.zero for v in got['values'].values()),'Original flat support must give exact zero A/B and first Z')
    density=current.density.density_Z_kernels(roots['E'][serial.ZERO],roots['E'][serial.DZ],
        scalar('.75'),scalar('.1'),got['values'],1024)
    require(all(v.zero for v in [*density['kernels'].values(),*density['Z_derivatives'].values()]),
        'Exact flat original loop must give zero increments, including pressure')
    return dict(passed=True,original_flat_primitive_C0_Z_zero=True,all_five_density_C0_Z_exact_zero=True,
        zero_increments_do_not_assert_incoming_correction_zero=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2 and saved['actual_original_O2_cell_count']==8,
        'Actual original eight-cell O2 history route required')
    require(not saved['global_inlet_to_Rc_histories_admitted'] and not any(saved.get(k) for k in packets.OPEN),
        'Conservative O2 histories cannot admit terminal/global gates')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed original O2 source: '+name)
    previous=json.loads((HERE/current.preceding.NAME).read_bytes())
    inherited=json.loads((HERE/current.preceding.RECEIPT).read_bytes())
    require(inherited['all_passed'] and inherited[current.preceding.GATE],'Accepted actual O2-inlet source required')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeO2C1Histories(current.make_middle_owner(bridge))
        independent=independent_original_loop_checks(owner.ctx)
        flat=exact_flat_primitive_and_density_check(owner.ctx)
        require(saved['original_periodic_parameter_C1_theorem']==serial.periodic_parameter_theorem(),
            'Original signed-r parameter proof changed')
        rows=0;incoming_rows=0;local_rows=0;composite_rows=0;terminal_rows=0;regions={}
        expected=[dict(label=label,chart=chart,left=left,right=right) for label,chart,left,right in current.ROUTE]
        for name,old in saved['actual_original_inlet_to_O2_exit_C1_records'].items():
            got=owner.route(packets.interval(owner.ctx,old['Z_box']),saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual original O2 C1 route changed: '+name)
            require(old['known_actual_original_O2_inlet_C1_history']==previous['actual_inlet_to_O2_inlet_C1_history_records'][name],
                'Actual accepted O2 coordinate0 correction was changed')
            require(old['ordered_original_O2_cells']==expected and old['no_original_interval_skipped_from_true_inlet_to_O2_exit'],
                'Full original O2 route must have no gaps')
            incoming=got['initial']['correction'];previous_chart,previous_endpoint='O2_slope',0
            require(any(not v.zero for v in incoming['values'].values()) and any(not v.zero for v in incoming['Z_derivatives'].values()),
                'Actual inherited history boxes cannot be zeroed')
            for index,(label,chart,left,right) in enumerate(current.ROUTE,1):
                cell=got['cells'][label];row=old['actual_O2_serial_cell_C1_records'][label]
                require(row['original_left_endpoint']==dict(chart=previous_chart,coordinate=previous_endpoint),
                    'Wrong actual inherited endpoint: '+label)
                require(cell['geometry']['coordinate']._mpi_==owner.ctx.mpf((left,right))._mpi_,
                    'Entire original O2 cell domain required: '+label)
                geometry=cell['geometry']['record']
                require(geometry['width_and_endpoints_independent_of_Z'] and ep(cell['geometry']['width'].coefficient)[0]>0,
                    'Original positive fixed-Z width required: '+label)
                if chart==previous_chart:
                    require(row['original_background_history_P0_join_receipt']['same_original_chart_and_source_function'],
                        'Same original source function required inside a chart')
                    require(current.density.spatial.exact_coordinate(left)==current.density.spatial.exact_coordinate(previous_endpoint),
                        'Same-chart adjacent endpoint gap')
                else:
                    seam=previous_chart+' -> '+chart
                    interface='slope_axial' if chart=='O2_axial' else 'axial_buffer'
                    require(row['original_same_radius_phase_or_same_chart_seam']==seam
                        and seam in geometry['exact_original_radius_Jacobian_identities']['original_same_radius_periodic_phase_seam_identities'],
                        'Exact original O2 radius/phase seam absent')
                    require(row['original_background_history_P0_join_receipt']==owner.joins[interface],
                        'Original source/P0 join must be separately bound')
                if chart=='O2_axial':
                    M=owner.ctx.mpf(owner.transfer.geometry.binder.seed.params.Md)
                    require(cell['geometry']['regular']._mpi_==owner.ctx.expm1(M)._mpi_,
                        'Full original axial width must be exp(Md)-1')
                if chart=='O2_buffer':
                    require(row['smaller_relaxed_buffer_admission_scope']=='selector[0,9] = shared_offset[-11,-2]'
                        and row['last_two_units_source_covered_without_extended_cone_admission']==(left==9)
                        and not row['global_strict_or_relaxed_cone_admission_from_this_integral'],
                        'Original buffer source domain and smaller cone admission must remain distinct')
                require(row['original_chart_uniform_a_positive_certificate']['source_function_positivity_not_inferred_from_saved_box'],
                    'Original chart-uniform positive theorem required')
                proof=row['original_periodic_parameter_C1_cover']
                require(proof['original_periodic_parameter_C1_theorem']==serial.periodic_parameter_theorem(),
                    'Parameter-uniform whole-period derivative theorem required')
                if not proof.get('original_full_box_q_and_q_Z_exact_zero_implies_primitive_C0_Z_exact_zero'):
                    require(proof['whole_period_derivatives_avoid_Poisson_denominator_powers']
                        and proof['original_full_source_derivatives_not_clipped_or_erased']
                        and proof['actual_fractional_phase_covered_without_samples_or_selected_inverse'],
                        'Original whole-period source derivatives cannot be sampled or clipped')
                require(row['old_unresolved_q_rows_not_used_for_first_Z_primitive_bounds']
                    and row['original_true_width_and_native_ordinary_y_conversion_applied_once'],
                    'Whole source roots and single true-width conversion required')
                require(row['O2_inlet_to_current_endpoint_C1_operator']['steps']==index,'Cumulative operator lost a cell')
                for key in current.RATES:
                    same(cell['incoming']['values'][key],incoming['values'][key],'Actual inherited C0 memory lost: '+label+' '+key)
                    same(cell['incoming']['Z_derivatives'][key],incoming['Z_derivatives'][key],'Actual inherited Z memory lost: '+label+' '+key)
                    decay=cell['factors'][key]['decay']
                    same(cell['correction']['values'][key],decay*incoming['values'][key]+cell['contributions'][key],
                        'Original signed C0 transfer mismatch: '+label+' '+key)
                    same(cell['correction']['Z_derivatives'][key],decay*incoming['Z_derivatives'][key]+cell['Z_derivatives'][key],
                        'Original signed Z transfer mismatch: '+label+' '+key)
                    same(cell['own'][key],cell['background']['originals'][key]+cell['correction']['values'][key],
                        'Original background reset or double-added: '+label+' '+key)
                    same(cell['own_Z'][key],cell['background']['Z_derivatives'][key]+cell['correction']['Z_derivatives'][key],
                        'Original background derivative reset or double-added: '+label+' '+key)
                require(ep(cell['factors']['p']['decay'].coefficient)==(1,1)
                    and ep(cell['factors']['p']['decay'].scale.evaluate())==(0,0),
                    'Rate-zero pressure memory must retain coefficient1')
                require(row['original_right_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history'],
                    'Original analytic P0 and P0_Z must stay separate')
                for value in [*cell['own'].values(),*cell['own_Z'].values(),*cell['correction']['values'].values(),*cell['correction']['Z_derivatives'].values()]:
                    require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,
                        'Every original history must use common arithmetic bases/ledger')
                require(not row['global_inlet_to_Rc_histories_admitted'] and not any(row.get(k) for k in packets.OPEN),
                    'O2 cell cannot admit quantitative global gates')
                rows+=20;incoming_rows+=10;incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
            local=got['local'].apply(got['local_incoming']['values'],got['local_incoming']['Z_derivatives'],owner.family)
            composite=got['cumulative'].apply(got['initial']['correction']['values'],got['initial']['correction']['Z_derivatives'],owner.family)
            require(got['local'].steps==3 and got['cumulative'].steps==8,'Original local/whole operator lengths changed')
            for key in current.RATES:
                for kind in ('values','Z_derivatives'):
                    require(overlaps(local[kind][key],got['cells']['slope_14_15']['correction'][kind][key]),
                        'Actual .12 incoming three-cell transfer disagrees with serial cover: '+key);local_rows+=1
                    require(overlaps(composite[kind][key],got['correction'][kind][key]),
                        'Whole O2 composite transfer disagrees with serial cover: '+key);composite_rows+=1
            require(old['original_terminal_O2_O3_background_history_P0_join']==owner.joins['buffer_transition'],
                'Original terminal buffer/transition source-function join required')
            provenance=got['background']['record']['source_provenance']
            require(provenance['chart']=='O3_slope_mu' and ep(provenance['coordinate_box'])==(0,0),
                'O3 transition0 background must be queried directly')
            require(got['background']['record']['P0_not_merged_into_pressure_history'],'Terminal P0_Z must stay separate')
            for key in current.RATES:
                same(got['own'][key],got['background']['originals'][key]+got['correction']['values'][key],
                    'O3 transition0 own C0 mismatch: '+key)
                same(got['own_Z'][key],got['background']['Z_derivatives'][key]+got['correction']['Z_derivatives'][key],
                    'O3 transition0 own Z mismatch: '+key)
                require(overlaps(got['own'][key],got['cells']['buffer_9_11']['own'][key])
                    and overlaps(got['own_Z'][key],got['cells']['buffer_9_11']['own_Z'][key]),
                    'Same-function original buffer/transition history covers disagree: '+key);terminal_rows+=2
            require(old['original_O3_to_Rc_route_still_missing'] and not old['global_inlet_to_Rc_histories_admitted'],
                'Missing O3/Rc and quantitative closure must stay explicit')
            regions[name]=dict(actual_O2_cells=8,actual_correction_and_own_C0_Z_rows=160,
                actual_incoming_C0_Z_rows=80,local_12_15_comparison_rows=10,whole_O2_comparison_rows=10,
                directly_queried_O3_transition0_own_rows=10)
            print('Actual complete O2 C1 route checked:',name,flush=True)
        try:owner.route(N=159)
        except ValueError:pass
        else:raise ArithmeticError('Unsafe original whole-period frequency accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        native_Z_queries_checked=2,complete_added_original_charts_checked=3,actual_original_O2_cells_checked=8,
        actual_added_cell_correction_and_own_C0_Z_rows_checked=rows,actual_inherited_C0_Z_rows_checked=incoming_rows,
        actual_12_15_local_serial_composite_overlap_rows_checked=local_rows,
        whole_O2_serial_composite_overlap_rows_checked=composite_rows,actual_O3_transition0_own_C0_Z_rows_checked=terminal_rows,
        independent_original_loop_reference_checks=independent,original_exact_flat_check=flat,regions=regions,
        original_periodic_parameter_C1_theorem=serial.periodic_parameter_theorem(),
        accepted_original_O2_background_history_P0_source_joins_bound=4,
        actual_original_O2_12_incoming_and_three_cell_transfer_checked=True,
        original_full_O2_source_domains_true_widths_and_buffer_admission_distinctions_checked=True,
        rate0_pressure_C0_Z_memory_and_separate_P0_retained=True,
        conservative_covers_not_quantitative_terminal_closure=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Complete O2 C1 PASS:',rows,'actual cell rows;',incoming_rows,'inherited rows',flush=True)
    return result


if __name__=='__main__':run()
