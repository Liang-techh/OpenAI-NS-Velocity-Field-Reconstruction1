"""Quantitative q-flat input margins on the actual modulation domain.

The first bridge excludes Ra: its lower native phase is sc/2 at r_minus.
Original relative bridge errors are rebound on the whole strong subset;
collar-only numerical cone constants are not promoted to a whole theorem.
Saved original whole-source norms and continuous state bounds are reused.
"""
import json
from pathlib import Path
import time
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget as band

state, allN, packets = band.state, band.allN, band.packets
HERE, PREFIX, sha, ep = band.HERE, band.PREFIX, band.sha, band.ep
LogUpper, Poly = band.LogUpper, band.Poly
NAME=PREFIX+'current_native_q_flat_cone_budget.json'
RECEIPT=PREFIX+'current_native_q_flat_cone_budget_check.json'
GATE='current_original_restricted_q_flat_margins_and_modulation_C0_frequency_budget_bound'
NONEMPTY=('bridge_first','bridge_second','bridge_macro','O2_axial','O3_slope_mu','O3_power')
EMPTY=('switch_first','switch_second','switch_power','reshape','inner_reference','axial_restore',
       'restore_buffer','actual_patch','Rh_reference','O2_slope','O2_buffer')


def decode_poly(c,row):
    return Poly(c,{term['N_power']:LogUpper(c,packets.interval(c,term['coefficient']['log_absolute_upper']))
                   for term in row['terms']})


def bridge_relative_constants(c):
    """From actual K>=10^6, e<=omega/(40K^6), |qbar|>=1/(2K).

    rho<=1/(20K^5); kappa-2<=Hbar<=K^10 gives
    (kappa-2)*rho^2<=1/400 BEFORE taking independent ranges.
    """
    Kmin=c.mpf(1000000);rho=1/(20*Kmin**5)
    D=1-rho;Q=2*(1-rho)**2-c.mpf(1)/400
    if ep(D-c.mpf('.95'))[0]<=0 or ep(Q-c.mpf('1.8'))[0]<=0:
        raise ArithmeticError('Whole strong bridge relative reserve failed')
    return dict(actual_K_lower=Kmin,rho_upper=rho,D_over_omega_Hbar_lower=D,
                Q_over_omega_Hbar_squared_lower=Q,correlated_kappa_rho_squared_upper=c.mpf(1)/400,
                domain='all original bridge source points with kappa>2, full Z[-1,1]',
                source_error_gate='e(K)/omega <= 1/(40K^6) for actual K>=K1>=10^6',
                actual_K_correlation_retained=True,collar_numeric_D_Q_record_not_used=True)


def first_bridge_omega_log_floor(c,sc):
    """s>=sc/2>0, hb<1/2, original sigma denominator<=2.

    log omega >= -4/sc^2-log4. This bounds the log of stress;
    it neither evaluates tiny sigma nor changes the original width hb.
    """
    if ep(sc)[0]<=0 or ep(sc)[1]>=ep(c.mpf('.25'))[0]:
        raise ValueError('Accepted original selected 0<sc<1/4 required')
    lower=c.mpf(ep(-4/(sc*sc)-c.ln(4))[0])
    return dict(log_omega_positive_lower=lower,lower_phase='sc/2 at r_minus',
        selected_sc=sc,source_sigma='exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))',
        sigma_monotonic_on_open_unit_interval=True,original_hb_less_half=True,
        Ra_or_native_phase_zero_excluded=True,stress_scale_or_profile_not_materialized=True)


def weighted_margin(c,loga,logdelta,logD,logQ,norms):
    """G=(a,a*Delta,a*D,a^3*Q), Delta=kappa-2>=eta.

    H covers ALL four original state components, including a,b; the
    active-branch bound |a|,|b|<=3 is not imposed on q-flat input.
    """
    logs=dict(G1=loga,G2=loga+logdelta,G3=loga+logD,G4=3*loga+logQ)
    logg=c.mpf(min(ep(v)[0] for v in logs.values()))
    rows=[LogUpper.constant(c,3)]
    for key in ('a','b','p1','p2'):
        record=norms[key]['y0_Z0']
        if not record['exact_zero']:rows.append(LogUpper(c,packets.interval(c,record['log_absolute_upper'])))
    maxlog=c.mpf(max(ep(v.log)[1] for v in rows))
    logH=LogUpper.add(c,[LogUpper.constant(c,2),LogUpper(c,maxlog)]).log
    logrho=c.mpf(min(0,ep(logg-c.ln(512)-5*logH)[0]))
    return dict(G_log_positive_lowers=logs,log_g_positive_lower=logg,
        log_H_state_upper=logH,log_rho_stability_positive=logrho,
        log_D_positive_lower=logD,log_Q_positive_lower=logQ,
        original_a_b_p1_p2_whole_chart_norms_used=True,
        tolerance='max state error <= min(1,g/(512H^5)) gives G1..G4 >= g/2',
        domain='actual modulation-domain source subset kappa-2>=eta',
        stress_free_Ra_not_in_modification_domain=True)


class SavedOriginalQFlatBudget:
    """Use checked native artifacts, not source-ancestor reconstruction."""
    def __init__(self):
        self.ctx=MPIntervalContext();self.ctx.dps=300;self.hashes={};self.rows={}
        checked=json.loads((HERE/band.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[band.GATE]:raise ValueError('Accepted native band budget required')
        self.family=checked['source_family'];self.hashes.update(checked['input_hashes'])
        self.hashes[band.RECEIPT]=sha(band.RECEIPT)
        for stem in ('current_native_repair_band_cone_budget','current_native_generic_cone_state_errors',
            'current_generic_shear_uniform_inputs','current_generic_shear_source_bounds','current_generic_shear_O3_sources',
            'current_inner_relaxed_inputs','current_inner_exit_strict_collar','global_exit_certificate','K1_ledger',
            'current_O2_axial_relaxed_cone','current_O3_transition_direction','current_O3_power_cone'):
            name=PREFIX+stem+'.json';expected=self.hashes.get(name)
            if expected is not None and expected!=sha(name):raise ValueError('Accepted original input changed: '+name)
            self.rows[stem]=json.loads((HERE/name).read_bytes());self.hashes[name]=sha(name)
        for stem in ('current_generic_shear_uniform_inputs','current_generic_shear_source_bounds','current_generic_shear_O3_sources','current_inner_relaxed_inputs'):
            if self.rows[stem]['source_family']!=self.family:raise ValueError('Original source family differs: '+stem)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def margins(self):
        c=self.ctx;read=lambda row:packets.interval(c,row)
        uniform=self.rows['current_generic_shear_uniform_inputs'];scales=uniform['current_actual_logarithmic_loop_scales']
        eta=read(scales['selected_positive_eta_log'])
        amin=read(scales['logarithmic_selected_positive_lower_constants']['a_min'])
        prior=self.rows['current_generic_shear_source_bounds']['current_original_source_log_bound_charts']
        outer=self.rows['current_generic_shear_O3_sources']['original_O3_quotient_log_norms']
        norms={chart:row['admitted_original_quotient_log_norms'] for chart,row in {**prior,**outer}.items() if chart!='core'}
        if set(norms)!=set(NONEMPTY+EMPTY):raise ArithmeticError('All17 original whole-chart norm inventory required')
        ledger=self.rows['K1_ledger'];exit=self.rows['global_exit_certificate']
        if ep(read(exit['strong_cone_sufficient_ratio_margin']))[0]<=0:
            raise ArithmeticError('Original pointwise strong bridge error gate required')
        if ep(1-read(ledger['decreasing_K_smallness_bounds']['strong_cone_error_ratio']))[0]<=0:
            raise ArithmeticError('Original decreasing-K ratio does not imply the uniform strong error gate')
        if ep(c.ln(c.mpf('.5'))-read(ledger['shared_positive_width_log_enclosure']))[0]<=0:
            raise ArithmeticError('Original hb<1/2 required')
        inner=self.rows['current_inner_relaxed_inputs']['current_whole_inner_relaxed_input_and_bounds']
        if ep(read(inner['positive_directed_margins']['original_K_exceeds_one_million_log_gap']))[0]<=0:
            raise ArithmeticError('Actual K>=10^6 source proof required')
        relative=bridge_relative_constants(c)
        collar=self.rows['current_inner_exit_strict_collar']['explicit_current_inner_exit_strict_collar']
        sc=read(collar['selected_first_phase_endpoint']);first=first_bridge_omega_log_floor(c,sc)
        if uniform['current_chosen_domain']['formal_radii']['r_minus']!='Ra*exp(hb*s_c/2)':
            raise ArithmeticError('Positive lower first native phase not bound to actual r_minus')
        Drel=relative['D_over_omega_Hbar_lower'];Qrel=relative['Q_over_omega_Hbar_squared_lower']
        records={}
        for chart in ('bridge_first','bridge_second','bridge_macro'):
            logomega=first['log_omega_positive_lower'] if chart=='bridge_first' else -c.ln(2)
            logscale=logomega+c.ln(2) # Hbar>=2+2gamma>=2.
            records[chart]=dict(margin=weighted_margin(c,amin,eta,logscale+c.ln(Drel),
                2*logscale+c.ln(Qrel),norms[chart]),original_relative_source=relative,
                actual_first_phase_floor=first if chart=='bridge_first' else None,
                original_omega_lower='log bound at sc/2' if chart=='bridge_first' else '1-hb>1/2',
                full_Z_domain=[-1,1],original_source_subset='kappa-2>=eta, actual r_minus<=R<=r_plus')
        attach=uniform['whole_actual_original_generic_input_margin']['source_margin_attachment_conditions']
        axial=self.rows['current_O2_axial_relaxed_cone']['whole_original_axial_relaxed_cone']
        logD=read(attach['axial_log_D_positive_lower']);Dfactor=read(axial['full_D_over_theta_lower'])
        # The accepted uniform producer factors logD as
        # lower(logRref+1+log(canonical_theta_reserve)+log(Dfactor)).
        # Removing that same factor recovers a conservative original theta
        # floor. This is NOT an inference of theta from D>=Dfactor*theta.
        logtheta=logD-c.ln(Dfactor)
        logQ=2*logtheta+c.ln(read(axial['full_Q_over_theta_squared_lower']))
        records['O2_axial']=dict(margin=weighted_margin(c,c.ln(2),eta,logD,logQ,norms['O2_axial']),
            source_subset='b^2/2>=eta; a=2, excludes zero-shear edges and midplane',
            log_normalized_theta_positive_lower=logtheta,
            original_factored_theta_floor_recipe='current_generic_shear_uniform_inputs.margins.axiallog',
            theta_floor_recovered_from_original_factored_log_D_recipe=True,
            full_signed_D_over_theta_lower=Dfactor,
            full_signed_Q_over_theta_squared_lower=read(axial['full_Q_over_theta_squared_lower']),
            full_pressure_energy_and_radial_sectors_retained=True)
        for chart,stem in (('O3_slope_mu','O3_transition'),('O3_power','O3_power')):
            logD=read(attach[stem+'_full_Ttheta_over_F_log_lower'])
            Qfactor=read(attach[stem+'_full_signed_Q_over_Ttheta_squared_lower'])
            if ep(Qfactor)[0]<=0:raise ArithmeticError('Original signed full O3 direction reserve required')
            delta=eta if chart=='O3_slope_mu' else c.ln(2*read(uniform['current_chosen_domain']['right_collar_mu']))
            records[chart]=dict(margin=weighted_margin(c,c.ln(2),delta,logD,2*logD+c.ln(Qfactor),norms[chart]),
                source_subset='2mu*sigma>=eta' if chart=='O3_slope_mu' else '2mu>=eta',
                original_a='2+2mu*sigma' if chart=='O3_slope_mu' else '2+2mu',original_b_exact_zero=True,
                full_signed_Q_over_D_squared_lower=Qfactor,full_energy_and_absolute_pressure_retained=True)
        groups=uniform['whole_actual_original_generic_input_margin']['chart_groups']
        proofs={chart:dict(kappa_upper=c.mpf(1),source_proof='global_exit_certificate.proof.switch_cone',
            source_recipe=exit['field_prescription'][recipe]) for chart,recipe in
                (('switch_first','first_switch'),('switch_second','second_switch'),('switch_power','last_interval'))}
        proofs['switch_first']['kappa_upper']=read(ledger['decreasing_K_smallness_bounds']['small_shear_kappa'])
        for chart in ('switch_second','switch_power'):proofs[chart]['kappa_upper']=c.mpf('.8')
        shape=groups['reshape_reference']['original_source_proof']
        if not shape['whole_path_kappa_equals_a_below1']:raise ArithmeticError('Original reshape weak source gate required')
        for chart in ('reshape','inner_reference'):
            proofs[chart]=dict(kappa_upper=c.mpf(1),source_proof='uniform.chart_groups.reshape_reference',
                exact_b_zero_and_kappa_equals_a=True,original_below_one_theorem=True)
        for chart,key in (('axial_restore','current_restore_relaxed_inputs'),('restore_buffer','current_restore_relaxed_inputs'),
                          ('actual_patch','current_patch_relaxed_inputs')):
            proofs[chart]=dict(kappa_upper=read(groups[key]['original_source_proof']['original_kappa_upper']),
                source_proof='uniform.chart_groups.'+key+'.original_source_proof.original_kappa_upper')
        ref=groups['O2_reference_slope']['original_source_proof']
        if ref['actual_bs']!=0:raise ArithmeticError('Original O2 reference/slope exact zero b required')
        for chart in ('Rh_reference','O2_slope'):
            proofs[chart]=dict(kappa_upper=read(ref['actual_kappa_range'][1]),
                source_proof='uniform.chart_groups.O2_reference_slope',original_b_exact_zero=True)
        proofs['O2_buffer']=dict(kappa_upper=c.mpf(2),source_proof='uniform.chart_groups.O2_buffer_first/O2_buffer_last',
            original_a_exact_two_b_exact_zero=True)
        empty={}
        for chart in EMPTY:
            proof=proofs[chart]
            if ep(proof['kappa_upper'])[1]>2:raise ArithmeticError('Cannot remove a nonempty original q-flat subset: '+chart)
            empty[chart]=dict(q_flat_subset_empty=True,eta_strictly_positive=True,**proof,
                active_branch_certificate_required_and_retained=True)
        return dict(nonempty=records,empty=empty,relative=relative,first=first,norms=norms,log_eta=eta)

    def compute(self):
        c=self.ctx;original=self.rows['current_native_generic_cone_state_errors'];bands=self.rows['current_native_repair_band_cone_budget']
        margins=self.margins();cells=[];largest=c.mpf(0);count=0
        for cell in original['original_continuous_cell_state_error_records']:
            chart=cell['chart'];branches=[]
            for row in cell['original_branch_state_error_records']:
                if chart not in margins['nonempty']:continue
                margin=margins['nonempty'][chart]['margin']
                polynomials={key:decode_poly(c,v) for key,v in row['normalized_state_error_log_polynomials'].items()}
                requirements=state.negative_power_threshold(c,polynomials,margin['log_rho_stability_positive'])
                logN=c.mpf(max(ep(v['sufficient_log_N_lower'])[1] for v in requirements.values()))
                largest=c.mpf(max(ep(largest)[1],ep(logN)[1]));count+=1
                branches.append(dict(original_source_provenance=row['original_source_provenance'],
                    original_cutoff_branch=row['original_cutoff_branch'],
                    original_continuous_state_polynomials=row['normalized_state_error_log_polynomials'],
                    q_flat_requirements=requirements,q_flat_sufficient_log_N_lower=logN,
                    conditional_only_on_original_kappa_minus2_ge_eta=True))
            cells.append(dict(label=cell['label'],chart=chart,original_geometry=cell['original_geometry'],
                q_flat_subset_empty=chart in margins['empty'],q_flat_branch_records=branches,
                original_initial_zero_correction=cell['exact_original_zero_initial_collar'],
                all_original_active_conditional_branches_retained=len(cell['original_branch_state_error_records'])))
        if len(cells)!=24 or original['actual_source_branches']!=35:raise ArithmeticError('Whole24-cell/35-branch modulation route required')
        combined=c.mpf(max(ep(largest)[1],ep(packets.interval(c,bands['source_repair_active_quiet_band_sufficient_log_N_lower']))[1]))
        return dict(margins=margins,cells=cells,q_flat_branch_count=count,q_flat_logN=largest,combined_logN=combined)


def run():
    began=time.monotonic();owner=SavedOriginalQFlatBudget();got=owner.compute()
    report=dict(source_family=owner.family,**{GATE:True},
        actual_original_q_flat_nonempty_charts=got['margins']['nonempty'],
        actual_original_q_flat_empty_charts=got['margins']['empty'],
        whole_strong_bridge_relative_cone_rebinding=got['margins']['relative'],
        actual_first_bridge_modification_floor=got['margins']['first'],
        original_selected_eta_log=got['margins']['log_eta'],
        original_continuous_modulation_q_flat_records=got['cells'],
        q_flat_original_conditional_branch_count=got['q_flat_branch_count'],
        q_flat_sufficient_log_N_lower=got['q_flat_logN'],
        source_repair_active_flat_quiet_band_C0_sufficient_log_N_lower=got['combined_logN'],
        all17_original_q_flat_chart_cases_accounted=True,
        local_modulation_and_repair_C0_cone_budget_bound=True,
        original_whole_source_norms_not_native_point_query_caps_used=True,
        phase_zero_Ra_not_assigned_a_positive_stress_floor=True,
        existing_original_native_graph_and_continuous_budget_records_reused=True,
        original_source_ancestor_constructors_called=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Quantitative original restricted q-flat input margins on actual r_minus..r_plus, first-phase cutoff before Ra, source-proved empty switches/middle/weak charts, and combined all24-cell C0 modulation plus conditional quiet/repair-band frequency budgets. Not whole N, evaluated controls/point field, terminal closure, outer/global cone or recursion.')
    (HERE/NAME).write_text(json.dumps(packets.encode(report),indent=2)+'\n',encoding='utf8')
    print('Original restricted q-flat margins and composed local C0 cone budget produced',flush=True)
    return report


if __name__=='__main__':run()
