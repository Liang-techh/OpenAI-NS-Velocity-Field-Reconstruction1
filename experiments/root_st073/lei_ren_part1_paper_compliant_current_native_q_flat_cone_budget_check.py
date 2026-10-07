"""Independent source identities, restricted margins and local N budgets.

Only saved original native records are read. No ancestor construction or
point-field oracle is needed to check this directed local budget layer.
"""
import ast
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_q_flat_cone_budget as producer

HERE, PREFIX, sha = producer.HERE, producer.PREFIX, producer.sha
packets, ep = producer.packets, producer.ep


def require(condition, message):
    if not condition:raise ArithmeticError(message)


def symbolic_source_and_stability():
    a,b,p1,p2=s.symbols('a b p1 p2',nonzero=True)
    delta=a+b*b/a-2;D=p1-b*p2/a-a-b*b/a;J=p2+b*p1/a
    G=[a,a*a+b*b-2*a,a*p1-b*p2-a*a-b*b]
    G.append(2*a*G[2]**2-G[1]*(b*p1+a*p2)**2)
    expected=[a,a*delta,a*D,a**3*(2*D*D-delta*J*J)]
    for got,want in zip(G,expected):require(s.cancel(got-want)==0,'Independent weighted q-flat G identity differs')
    grad=[]
    for expr in G:
        polys=[s.Poly(s.diff(expr,v),a,b,p1,p2) for v in (a,b,p1,p2)]
        require(all(p.total_degree()<=5 for p in polys),'State gradient degree exceeds five')
        grad.append(int(sum(abs(v) for p in polys for v in p.coeffs())))
    require(grad==[1,6,8,208] and s.Rational(208,512)<s.Rational(1,2),'Whole state stability factor insufficient')
    db,eb,chi,et,ez=s.symbols('Dbar Ebar chi e_theta e_z',nonzero=True)
    hb=(db*db+eb*eb)/db;omega=1-chi
    actualD=D.subs({a:chi*db,b:-chi*eb,p1:db+et,p2:eb+ez})
    actualJ=J.subs({a:chi*db,b:-chi*eb,p1:db+et,p2:eb+ez})
    require(s.cancel(actualD-(omega*hb+et+eb*ez/db))==0,'Original bridge D projection differs')
    require(s.cancel(actualJ-(ez-eb*et/db))==0,'Original bridge transverse projection differs')
    K=s.symbols('K',positive=True);C=s.symbols('C',positive=True)
    error=C*K**-80;ratio=40*K**6*error
    require(s.simplify(s.diff(ratio,K)+74*ratio/K)==0,'Pointwise strong ratio does not decrease with actual K')
    rho=1/(20*K**5)
    require(s.simplify(K**10*rho*rho-s.Rational(1,400))==0,'K dependence was lost before the cone reduction')
    x=s.symbols('x',positive=True)
    log_odds=-1/x**2+1/(1-x)**2
    require(s.simplify(s.diff(log_odds,x)-2/x**3-2/(1-x)**3)==0,'Original flat cutoff monotonicity differs')
    return dict(passed=True,weighted_G_identities=4,original_bridge_projection_identities=2,
        independent_gradient_coefficient_sums=grad,gradient_degree_at_most_five=True,
        pointwise_K_correlation_cancelled_before_independent_ranges=True,
        monotone_original_sigma_on_open_unit_interval=True)


def original_axial_factored_recipe(owner):
    name=PREFIX+'current_generic_shear_uniform_inputs.py'
    require(owner.hashes[name]==sha(name),'Accepted original axial recipe source changed')
    tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(q for q in tree.body if isinstance(q,ast.ClassDef) and q.name=='CurrentUniformInputs')
    fn=next(q for q in cls.body if isinstance(q,ast.FunctionDef) and q.name=='margins')
    assignments={ast.unparse(q.targets[0]):q.value for q in fn.body if isinstance(q,ast.Assign) and len(q.targets)==1}
    expected={
        'axial':"self.rows['current_O2_axial_relaxed_cone']",
        'base':"axial['whole_original_axial_baseline']",
        'proof':"axial['whole_original_axial_relaxed_cone']",
        'axiallog':"self.service.data['logRref'] + 1 + c.ln(read(base['canonical_theta_positive_reserve'])) + c.ln(read(proof['full_D_over_theta_lower']))",
        "conditions['axial_log_D_positive_lower']":"lower(axiallog)"}
    for key,value in expected.items():
        require(key in assignments and ast.dump(assignments[key])==ast.dump(ast.parse(value,mode='eval').body),
                'Accepted source no longer factors the original canonical theta floor: '+key)
    axial=owner.rows['current_O2_axial_relaxed_cone']
    bound=owner.rows['current_generic_shear_uniform_inputs']['whole_actual_original_generic_input_margin']['chart_groups']['O2_axial']['original_source_proof']
    require(bound['baseline']==axial['whole_original_axial_baseline'] and bound['signed_cone']==axial['whole_original_axial_relaxed_cone'],
            'Factored attachment is not the same original axial baseline and signed cone')
    require(ep(packets.interval(owner.ctx,bound['baseline']['canonical_theta_positive_reserve']))[0]>0,
            'Original canonical theta reserve is not positive')
    return dict(passed=True,original_source_AST_factored_recipe_assignments=len(expected),
        canonical_theta_reserve_source_and_signed_D_factor_bound=True,
        theta_floor_not_inferred_from_two_unrelated_lower_bounds=True,
        directed_subtraction_only_recovers_original_product_recipe=True)


def recorded_source_margins(owner,report):
    c=owner.ctx;read=lambda row:packets.interval(c,row)
    uniform=owner.rows['current_generic_shear_uniform_inputs'];groups=uniform['whole_actual_original_generic_input_margin']['chart_groups']
    ledger=owner.rows['K1_ledger'];exit=owner.rows['global_exit_certificate']
    require(ledger['smallness_checks']['strong_cone_error_ratio'] and ledger['smallness_checks']['small_shear_kappa'],
            'Original decreasing-K bridge/switch gates missing')
    K=c.mpf(1000000);ch=c.mpf(ep(read(ledger['cstar']))[1]);S=80*ledger['K1']*ch*K**-74
    storedS=read(ledger['decreasing_K_smallness_bounds']['strong_cone_error_ratio'])
    require(ep(S)[0]<=ep(storedS)[1] and ep(S)[1]>=ep(storedS)[0] and ep(storedS)[1]<1,
            'Strong source ratio is not the original pointwise K^-74 bound at Kmin')
    require(ep(read(exit['strong_cone_sufficient_ratio_margin']))[0]>0,'Original strong gate missing')
    relative=report['whole_strong_bridge_relative_cone_rebinding'];rhomax=1/(20*K**5)
    require(ep(read(relative['rho_upper']))[1]>=ep(rhomax)[1],'Bridge relative error cap too small')
    Drel=read(relative['D_over_omega_Hbar_lower']);Qrel=read(relative['Q_over_omega_Hbar_squared_lower'])
    require(ep(Drel)[0]<=ep(1-rhomax)[1] and ep(Drel)[0]>ep(c.mpf('.95'))[1],'Bridge D reserve invalid')
    require(ep(Qrel)[0]<=ep(2*(1-rhomax)**2-c.mpf(1)/400)[1] and ep(Qrel)[0]>ep(c.mpf('1.8'))[1],'Bridge Q reserve invalid')
    first=report['actual_first_bridge_modification_floor'];sc=read(owner.rows['current_inner_exit_strict_collar']['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint'])
    require(first['selected_sc']==packets.encode(sc),'First phase is not the accepted original sc')
    require(uniform['current_chosen_domain']['formal_radii']['r_minus']=='Ra*exp(hb*s_c/2)' and first['lower_phase']=='sc/2 at r_minus',
            'Positive cutoff floor does not belong to the actual modification inlet')
    require(ep(read(first['log_omega_positive_lower']))[1]<=ep(-4/(sc*sc)-c.ln(4))[1],
            'First bridge stress logarithm lower floor too large')
    require(first['Ra_or_native_phase_zero_excluded'] and ep(c.ln(c.mpf('.5'))-read(ledger['shared_positive_width_log_enclosure']))[0]>0,
            'Ra phase zero or original width condition mishandled')
    require(set(report['actual_original_q_flat_nonempty_charts'])==set(producer.NONEMPTY) and
            set(report['actual_original_q_flat_empty_charts'])==set(producer.EMPTY),'Original seventeen-chart inventory differs')
    for chart,row in report['actual_original_q_flat_empty_charts'].items():
        require(ep(read(row['kappa_upper']))[1]<=2 and row['q_flat_subset_empty'] and row['active_branch_certificate_required_and_retained'],
                'Empty q-flat case loses its original weak/active route: '+chart)
    for chart,key in (('axial_restore','current_restore_relaxed_inputs'),('restore_buffer','current_restore_relaxed_inputs'),('actual_patch','current_patch_relaxed_inputs')):
        require(report['actual_original_q_flat_empty_charts'][chart]['kappa_upper']==groups[key]['original_source_proof']['original_kappa_upper'],
                'Empty subset kappa cap changed from original source: '+chart)
    require(groups['reshape_reference']['original_source_proof']['whole_path_kappa_equals_a_below1'],'Original reshape theorem absent')
    ref=groups['O2_reference_slope']['original_source_proof']
    require(ref['actual_bs']==0 and ep(read(ref['actual_kappa_range'][1]))[1]==2,'Original O2 reference/slope weak limit changed')
    require('first kappa<=epsilon Hbar<1' in exit['proof']['switch_cone'] and '0<a<=.8' in exit['proof']['switch_cone'],
            'Original switch theorem does not remove the three strong subsets')
    # The checker derives the same lower floors from source identities, not
    # from sampled points or by evaluating defining functions at cap values.
    nonempty=report['actual_original_q_flat_nonempty_charts'];sources=owner.rows['current_generic_shear_source_bounds']['current_original_source_log_bound_charts']
    sources={**sources,**owner.rows['current_generic_shear_O3_sources']['original_O3_quotient_log_norms']}
    scales=uniform['current_actual_logarithmic_loop_scales'];eta=read(scales['selected_positive_eta_log'])
    amin=read(scales['logarithmic_selected_positive_lower_constants']['a_min'])
    axial=owner.rows['current_O2_axial_relaxed_cone']['whole_original_axial_relaxed_cone'];attach=uniform['whole_actual_original_generic_input_margin']['source_margin_attachment_conditions']
    checked=0
    for chart,row in nonempty.items():
        margin=row['margin'];loga=amin if chart.startswith('bridge_') else c.ln(2);logdelta=eta
        if chart.startswith('bridge_'):
            logomega=read(first['log_omega_positive_lower']) if chart=='bridge_first' else -c.ln(2)
            scale=logomega+c.ln(2);logD=scale+c.ln(Drel);logQ=2*scale+c.ln(Qrel)
        elif chart=='O2_axial':
            logD=read(attach['axial_log_D_positive_lower']);theta=logD-c.ln(read(axial['full_D_over_theta_lower']))
            logQ=2*theta+c.ln(read(axial['full_Q_over_theta_squared_lower']))
            require(row['theta_floor_recovered_from_original_factored_log_D_recipe'],'Original factored axial recipe not declared')
        else:
            stem='O3_transition' if chart=='O3_slope_mu' else 'O3_power'
            logD=read(attach[stem+'_full_Ttheta_over_F_log_lower']);factor=read(attach[stem+'_full_signed_Q_over_Ttheta_squared_lower'])
            require(ep(factor)[0]>0 and row['full_energy_and_absolute_pressure_retained'],'Full signed original O3 cone not retained')
            logQ=2*logD+c.ln(factor)
            if chart=='O3_power':logdelta=c.ln(2*read(uniform['current_chosen_domain']['right_collar_mu']))
        floors=dict(G1=loga,G2=loga+logdelta,G3=loga+logD,G4=3*loga+logQ)
        for key,value in floors.items():
            require(ep(read(margin['G_log_positive_lowers'][key]))[0]<=ep(value)[1],'Original weighted margin too large: '+chart+'/'+key)
        logg=read(margin['log_g_positive_lower'])
        require(all(ep(logg)[1]<=ep(read(v))[0] for v in margin['G_log_positive_lowers'].values()),'Common G floor too large: '+chart)
        normlogs=[c.ln(3)]
        for key in ('a','b','p1','p2'):
            norm=sources[chart]['admitted_original_quotient_log_norms'][key]['y0_Z0']
            if not norm['exact_zero']:normlogs.append(read(norm['log_absolute_upper']))
        logH=read(margin['log_H_state_upper'])
        require(all(ep(logH-v)[0]>=0 for v in [c.ln(2),*normlogs]),'Whole a,b,p1,p2 cover incomplete: '+chart)
        require(ep(logH-c.mpf(max(ep(v)[1] for v in normlogs))-c.ln(2))[1]>=0,
                'State stability cover does not include the unit perturbation path: '+chart)
        require(ep(read(margin['log_rho_stability_positive']))[1]<=0 and
                ep(read(margin['log_rho_stability_positive'])-logg+c.ln(512)+5*logH)[0]<=0,
                'State tolerance exceeds original weighted reserve: '+chart)
        checked+=4
    return dict(passed=True,source_proved_empty_q_flat_charts=11,nonempty_chart_weighted_positive_floors=checked,
        all_four_original_whole_chart_state_norms_checked=True,source_pressure_energy_and_radial_sectors_preserved=True,
        original_K_decreasing_ratio_and_actual_first_phase_floor_bound=True,
        positive_floor_at_stress_free_Ra_not_claimed=True)


def recorded_coverage_and_frequencies(owner,report):
    c=owner.ctx;read=lambda v:packets.interval(c,v)
    original=owner.rows['current_native_generic_cone_state_errors'];got=report['original_continuous_modulation_q_flat_records']
    require(len(got)==len(original['original_continuous_cell_state_error_records'])==24 and original['actual_source_branches']==35,
            'Original continuous cell/branch route changed')
    count=0;states=0;combined=read(report['source_repair_active_flat_quiet_band_C0_sufficient_log_N_lower'])
    qflat=read(report['q_flat_sufficient_log_N_lower'])
    for old,new in zip(original['original_continuous_cell_state_error_records'],got):
        for key in ('label','chart','original_geometry'):require(old[key]==new[key],'Original cell source geometry changed: '+key)
        require(old['exact_original_zero_initial_collar']==new['original_initial_zero_correction'],'Initial exact-zero correction route changed')
        require(new['all_original_active_conditional_branches_retained']==len(old['original_branch_state_error_records']),
                'Original active conditional branch removed')
        if old['chart'] in producer.EMPTY:
            require(new['q_flat_subset_empty'] and not new['q_flat_branch_records'],'Empty subset unexpectedly assigned a strict cone')
            continue
        require(len(new['q_flat_branch_records'])==len(old['original_branch_state_error_records']),'Nonempty chart omits cutoff branch')
        rho=read(report['actual_original_q_flat_nonempty_charts'][old['chart']]['margin']['log_rho_stability_positive'])
        for prior,row in zip(old['original_branch_state_error_records'],new['q_flat_branch_records']):
            for key,target in (('original_source_provenance','original_source_provenance'),('original_cutoff_branch','original_cutoff_branch'),
                               ('normalized_state_error_log_polynomials','original_continuous_state_polynomials')):
                require(prior[key]==row[target],'Original signed continuous branch/state bound changed: '+key)
            require(row['conditional_only_on_original_kappa_minus2_ge_eta'],'q-flat budget lost subset condition')
            require(set(row['q_flat_requirements'])=={'a','b','p1','p2'},'One normalized state row omitted')
            largest=read(row['q_flat_sufficient_log_N_lower'])
            for key,poly in row['original_continuous_state_polynomials'].items():
                require(all(term['N_power']<=-1 for term in poly['terms']),'State includes nondecaying frequency terms')
                logs=[read(term['coefficient']['log_absolute_upper']) for term in poly['terms']]
                bound=c.mpf(max(ep(v)[1] for v in logs))+c.ln(len(logs)) if logs else None
                threshold=read(row['q_flat_requirements'][key]['sufficient_log_N_lower'])
                if bound is not None:require(ep(threshold-bound+rho)[1]>=0,'Directed sufficient N threshold below log-sum cap')
                require(ep(threshold)[0]>=0 and ep(largest-threshold)[0]>=0,'Branch frequency does not cover all states')
                states+=1
            require(ep(qflat-largest)[0]>=0,'Combined q-flat frequency omits original branch')
            count+=1
    band=read(owner.rows['current_native_repair_band_cone_budget']['source_repair_active_quiet_band_sufficient_log_N_lower'])
    require(count==report['q_flat_original_conditional_branch_count']==13 and states==52,'Original applicable branch count differs')
    require(ep(combined-qflat)[0]>=0 and ep(combined-band)[0]>=0,'Combined local N budget omits existing repair/active/quiet/band condition')
    for key in ('actual_five_controls_installed','certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed',
                'current_whole_N_selected','original_source_ancestor_constructors_called',*packets.OPEN):
        require(report[key] is False,'Local budget incorrectly closes the full construction: '+key)
    return dict(passed=True,original_continuous_cells=24,original_active_conditional_branches_retained=35,
        restricted_q_flat_conditional_branches=13,negative_power_state_requirements=52,
        source_active_repair_quiet_band_frequency_conditions_all_retained=True,
        actual_controls_tail_terminal_closure_global_N_outer_cone_and_recursion_remain_open=True)


def run():
    began=time.monotonic();owner=producer.SavedOriginalQFlatBudget();report=json.loads((HERE/producer.NAME).read_bytes())
    require(report[producer.GATE] and report['source_family']==owner.family,'Current original same-family q-flat output required')
    for name,value in owner.hashes.items():require(report['input_hashes'].get(name)==value,'Producer native dependency differs: '+name)
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_symbolic_source_and_stability=symbolic_source_and_stability(),
        independent_original_axial_factored_recipe=original_axial_factored_recipe(owner),
        actual_original_source_margins=recorded_source_margins(owner,report),
        actual_original_route_and_frequency_coverage=recorded_coverage_and_frequencies(owner,report),
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        checker_uses_saved_native_output_without_reconstructing_original_owners=True,
        local_modulation_and_repair_C0_cone_budget_bound=True,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes={**report['input_hashes'],producer.NAME:sha(producer.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Independent original restricted q-flat margins, axial product recipe and full local C0 route budget PASS',flush=True)
    return result


if __name__=='__main__':run()
