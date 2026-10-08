"""Whole-Z original source derivative refinement through the actual24 route.

All active cutoff alternatives are retained. Original C0 operators,
backgrounds, pressure datum, geometry and phase remain unchanged. Directed
derivative bounds refine source functions, never define field coefficients.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_full_predicate_24_cell_controls as base
import lei_ren_part1_paper_compliant_current_native_paired_C1_transport as paired
import lei_ren_part1_paper_compliant_current_native_active_kappa_mixed_conditioning as active

cutoff=paired.cutoff;q2=paired.accepted.accepted;prior=q2.prior
first=base.driver.middle.p;serial=paired.accepted.serial;density=first.density
HERE,PREFIX,sha=base.HERE,base.PREFIX,base.sha;ep=base.ep;KEYS=base.KEYS
ZERO=(0,0);DZ=(0,1);N=1024
NAME=PREFIX+'current_original_full_Z_paired_C1_controls.json.gz'
RECEIPT=PREFIX+'current_original_full_Z_paired_C1_controls_check.json'
GATE='actual_original_full_Z_correlated_cutoff_paired_C1_local_ranges_and_24_cell_control_diagnostics_installed'
ACTIVE_INDICES=(*range(1,13),*range(18,22))


def encode(value):
    """Retain exact scalar endpoint bits as well as the accepted interval codec."""
    if hasattr(value,'_mpf_') and not hasattr(value,'_mpi_'):
        return dict(decimal=mp.nstr(value,85),exact_mpf_tuple=list(value._mpf_))
    if isinstance(value,dict):return {key:encode(row) for key,row in value.items()}
    if isinstance(value,(list,tuple)):return [encode(row) for row in value]
    return base.encode(value)


def magnitude_log(value):
    return None if value.zero else ep(first.absolute_upper(value).scale.evaluate())[1]


def smaller_cover(old,new):
    """Select a complete same-function cover by certified magnitude only."""
    old.coerce(new);a,b=magnitude_log(old),magnitude_log(new)
    changed=b is None and a is not None or a is not None and b is not None and b<a
    return (new if changed else old),dict(old_absolute_upper_log=a,alternative_absolute_upper_log=b,
        selected_absolute_upper_log=b if changed else a,strict_absolute_upper_reduction=changed,
        selected_complete_source_function_cover_not_endpoint_value=True)


def correlated_active_roots(roots,branch,logamin):
    """Original active Delta_Z identity, with conditional C0 support only."""
    c=roots['a'][ZERO].ctx;scalar=roots['a'][ZERO].scalar
    restricted=dict(roots,kappa_minus2=dict(roots['kappa_minus2']))
    restricted['kappa_minus2'][ZERO]=branch['Delta']
    if branch['name']=='flat':return restricted,dict(exact_flat_branch_all_q_jets_zero=True)
    limit=2 if branch['name']=='negative' else c.mpf('2.5')
    a=first.positive_restriction(roots['a'][ZERO],logamin,c.ln(limit))
    if a is None:return None,dict(active_branch_empty_by_original_kappa_ge_a=True,active_a_upper=limit)
    b,bproof=q2.accepted.signed_C0_support_range(roots['b'][ZERO],scalar(c.mpf(('-1.25','1.25'))))
    derived_t0=(-b).positive_divide(a,logamin)
    support_t0=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=(c.ln(3)-logamin)/2),1,a.ledger)
    bound=serial.minimum_upper(derived_t0,roots['t0'][ZERO],support_t0)
    t0=first.symmetric_bound(bound)
    DeltaZ=(scalar(1)-first.current.square(t0))*roots['a'][DZ]-(t0*roots['b'][DZ])*2
    selected,decision=smaller_cover(roots['kappa_minus2'][DZ],first.symmetric_bound(DeltaZ))
    restricted['a']=dict(roots['a']);restricted['a'][ZERO]=a
    restricted['b']=dict(roots['b']);restricted['b'][ZERO]=b
    restricted['t0']=dict(roots['t0']);restricted['t0'][ZERO]=t0
    restricted['kappa_minus2'][DZ]=selected
    return restricted,dict(active_support='a<=kappa=2+Delta<=2+eta<=2.5; b²=a*(kappa-a)<=kappa²/4',
        original_conditional_a_C0=a.record(),original_conditional_b_C0=b.record(),b_support_intersection=bproof,
        original_active_t0_absolute_support=bound.record(),original_correlated_Delta_Z=DeltaZ.record(),
        original_raw_Delta_Z=roots['kappa_minus2'][DZ].record(),selected_same_function_Delta_Z=selected.record(),
        Delta_Z_comparison=decision,exact_Delta_Z_identity='(1-t0²)*a_Z-2*t0*b_Z; t0=-b/a; eta_Z=0',
        all_original_a_Z_b_Z_E_Z_p2_Z_and_higher_rows_retained=True,
        first_Z_function_enclosure_refined_not_derivative_of_support_caps=True,
        original_unconditioned_roots_not_mutated=True)


def whole_Z_paired_derivatives(roots,eta_log,logamin,dstar_log):
    c=roots['a'][ZERO].ctx
    if ep(eta_log)[1]>ep(-c.ln(2))[0]:raise ValueError('Same original constant eta<=1/2 required')
    candidates,empty=cutoff.conditional_cutoff_branches(roots['kappa_minus2'][ZERO],eta_log)
    values=[];rows=[]
    for branch in candidates:
        restricted,correlation=correlated_active_roots(roots,branch,logamin)
        if restricted is None:
            empty.append(dict(name=branch['name'],proof=correlation));continue
        q=cutoff.conditional_q_jet(restricted,eta_log,logamin,branch)
        square=q2.conditional_q2_jet(restricted,eta_log,logamin,branch)
        result=paired.paired_C1_support(restricted,q['rows'],square,logamin,dstar_log)
        values.append(result['values'])
        rows.append(dict(name=branch['name'],condition=branch['condition'],conditional_Delta_C0=branch['Delta'].record(),
            theta=branch.get('theta'),same_original_active_correlation=correlation,
            q_C0=q['rows'][ZERO].record(),q_Z=q['rows'][DZ].record(),direct_q_squared_C0=square[ZERO].record(),direct_q_squared_Z=square[DZ].record(),
            paired_support=result['record'],derivative_covers=base.records(result['values']),
            original_linear_q_Z_retained_no_q_squared_division_by_q=True))
    if not values:raise ValueError('Every original active/flat branch was excluded without a complete cover')
    union=density.local.same_source_union
    hull={key:union([row[key] for row in values]) for key in ('A_Z','B_Z_over_Pstar')}
    return dict(values=hull,record=dict(original_native_log_bases=roots['a'][ZERO].scale.bases,source_context_dps=c.dps,
        original_source_root_C0_Z={name:{str(order):row[order].record() for order in (ZERO,DZ)}
        for name,row in roots.items()},eta_log=eta_log,log_a_positive_lower=logamin,dstar_log=dstar_log,
        original_three_branch_cover='Delta<=0 union 0<=Delta<=eta union Delta>=eta',
        actual_conditional_branches=rows,branches_proved_empty=empty,whole_source_derivative_hulls=base.records(hull),
        complete_branch_union_precedes_whole_source_cover_selection=True,
        original_paired_parameter_theorem=paired.paired_parameter_theorem(),
        first_Z_only_no_higher_jet_or_source_point_oracle_admission=True))


class OriginalFullZPairedC1Controls:
    def __init__(self,bridge=None):
        self.base=base.OriginalFullPredicate24CellControls(bridge);self.target_owner=self.base.target_owner
        self.ctx=self.base.ctx;self.coordinates=self.base.coordinates;self.family=self.base.family;self.hashes=dict(self.base.hashes)
        self.accepted=base.driver.accepted(self.hashes,base,self.family)
        for module in (paired,q2,cutoff,active):base.driver.accepted(self.hashes,module,self.family)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.source=None

    @base.native.inlet.source_precision
    def refine(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Unchanged accepted actual original candidate N1024 required')
        source=self.base.assemble();c=self.ctx;coords=self.coordinates;root_owner=self.target_owner.q_owner.owner.owner
        signed=self.target_owner.transfer.owner.signed_owner;proofs=[]
        eta=base.packets.interval(c,root_owner.scales['selected_positive_eta_log'])
        dstar=base.packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        with mp.workdps(c.dps+40):
            for index in ACTIVE_INDICES:
                cell=source['cells'][index];spec=base.parameters.ROUTE[index];label,chart,left,right=spec
                geom=self.target_owner.geometry(*spec)
                if encode(geom['record'])!=cell['record']['actual_original_geometry']:raise ValueError('Original whole native source geometry changed')
                query=self.target_owner.q_owner.query(chart,(-1,1),geom['coordinate']);roots=query['source']['roots'];packet=query['source']['packet']
                positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
                positivity='whole_actual_source_positive_not_inferred_from_saved_denominator_box' if chart.startswith('O3_') else 'source_function_positivity_not_inferred_from_saved_box'
                if positive.get(positivity) is not True:raise ValueError('Original whole-chart positive a proof required')
                logamin=positive['log_actual_a_positive_lower']
                # Keep genuine old source C0 primitive objects; only first-Z bounds change.
                old=first.whole_period_C1(roots,eta,logamin,dstar) if index==1 else serial.whole_period_C1(roots,eta,logamin,dstar)
                alternative=whole_Z_paired_derivatives(roots,eta,logamin,dstar);selected=dict(old['values']);decisions={}
                for key in ('A_Z','B_Z_over_Pstar'):selected[key],decisions[key]=smaller_cover(old['values'][key],alternative['values'][key])
                if selected['A'] is not old['values']['A'] or selected['B_over_Pstar'] is not old['values']['B_over_Pstar']:
                    raise ValueError('Original C0 A/B primitive objects must remain unchanged')
                def axial(k):
                    row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
                    return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),roots['E'][ZERO].scale.bases,roots['E'][ZERO].ledger)
                got=density.density_Z_kernels(roots['E'][ZERO],roots['E'][DZ],axial(0),axial(1),selected,N)
                factors={key:base.preceding.downstream.transfer.true_width_kernel(coords,geom,rate) for key,rate in base.history.RATES.items()}
                newjets={key:coords.rebase(got['Z_derivatives'][key],self.family)*factors[key]['mass'] for key in KEYS}
                local_decisions={};old_jets=cell['Z_derivatives'];chosen={}
                for key in KEYS:chosen[key],local_decisions[key]=smaller_cover(old_jets[key],newjets[key])
                record=dict(label=label,chart=chart,original_cell_index=index,source_family=self.family,candidate_N=N,Z_box=c.mpf((-1,1)),
                    actual_original_geometry=geom['record'],fresh_original_source_query=query['record'],original_chart_positive_a_theorem=positive,
                    old_whole_period_source_primitive_C0_Z=old['record'],actual_complete_paired_cutoff_C1=alternative['record'],
                    primitive_derivative_comparisons=decisions,selected_original_primitive_C0_Z=base.records(selected),
                    original_density_velocity_and_increment_source_C0_Z=base.records(got['velocities']),
                    original_signed_density_Z_covers=base.records(got['Z_derivatives']),
                    original_true_width_kernel_factors={key:dict(mass=row['mass'].record(),decay=row['decay'].record()) for key,row in factors.items()},
                    accepted_local_Z_contributions=base.records(old_jets),paired_local_Z_contributions=base.records(newjets),
                    selected_local_Z_contributions=base.records(chosen),local_derivative_comparisons=local_decisions,
                    accepted_C0_local_contributions_unchanged=base.records(cell['values']),
                    actual_whole_Z_and_native_cell_requeried_not_strict_sign_archive_relabelled=True,
                    same_original_C0_A_B_E_EZ_V_VZ_mass_radius_phase_background_and_P0=True)
                cell['old_Z_derivatives']=old_jets;cell['Z_derivatives']=chosen;cell['operator'].Z_increments=chosen
                cell['record']['ordinary_Z_refinement_binding']=dict(module=Path(__file__).name,record_index=len(proofs))
                proofs.append(record)
                print('Actual whole-Z paired source:',label,'local Z reductions',sum(q['strict_absolute_upper_reduction'] for q in local_decisions.values()),flush=True)
            incoming=dict(values={key:coords.scalar(0) for key in KEYS},Z_derivatives={key:coords.scalar(0) for key in KEYS});budgets=[]
            for cell in source['cells']:
                oldrecord=cell['record'];op=cell['operator'];out=op.apply(incoming['values'],incoming['Z_derivatives'],self.family)
                bg=cell['background'];own={key:bg['values'][key]+out['values'][key] for key in KEYS}
                ownZ={key:bg['Z_derivatives'][key]+out['Z_derivatives'][key] for key in KEYS}
                for key in KEYS:
                    if encode(out['values'][key].record())!=encode(oldrecord['actual_right_correction_C0'][key]):raise ValueError('Accepted C0 affine history changed')
                cell['incoming']=incoming;cell['correction']=out
                oldrecord.update(actual_cell_C1_operator=op.record(),actual_local_Z_contributions=base.records(cell['Z_derivatives']),
                    actual_inherited_correction_Z=base.records(incoming['Z_derivatives']),actual_right_correction_Z=base.records(out['Z_derivatives']),
                    actual_right_own_history_Z=base.records(ownZ),absolute_pressure_Z=(bg['P0_Z']+ownZ['p']).record(),
                    accepted_original_C0_operator_background_P0_and_phase_unchanged=True)
                budgets.append(dict(label=oldrecord['label'],outgoing_Z_log_caps={key:magnitude_log(v) for key,v in out['Z_derivatives'].items()}))
                incoming=out
            amplitude=source['amplitude']
            target=dict(**base.fixed.fixed_N_target_rows(incoming['values'],incoming['Z_derivatives'],amplitude['A'],amplitude['AZ'],amplitude['logA'],amplitude['mu_source'],amplitude['logmu']),
                **{key:amplitude[key] for key in ('A','AZ','mu','logA','logmu')})
            source.update(history=incoming['values'],Z_derivatives=incoming['Z_derivatives'],target=target)
            source['record'].update(actual_original24_source_cell_records=[row['record'] for row in source['cells']],
                actual_Rc_correction_C0_Z=[base.records(incoming['values']),base.records(incoming['Z_derivatives'])],
                actual_Rc_joint_target_C0_Z=[base.records(target['values']),base.records(target['Z_derivatives'])],
                actual_Rc_joint_numerator_C0_Z=[target['joint_numerator'].record(),target['joint_numerator_Z'].record()],
                actual_full_Z_paired_C1_refinement_records=proofs,actual_updated_24_cell_Z_error_budget=budgets,
                original_C0_history_and_background_P0_and_N_and_geometry_unchanged=True,
                fresh_complete_body_transition_flat_derivative_unions_not_saved_sign_tiles=True)
        self.source=source;self.hashes.update(self.base.hashes);self.hashes.update(self.target_owner.service.hashes)
        self.base.hashes.update(self.hashes)
        return source

    def finite_controls(self,source,*,iterations=3):
        if source is not self.source or source is None:raise ValueError('Issued refined actual original24 function ranges required')
        result=self.base.finite_controls(source,iterations=iterations);self.hashes.update(self.base.hashes)
        return result


@base.native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=base.native.inlet.native_bridge_owner()
    with base.native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalFullZPairedC1Controls(bridge);source=owner.refine();result=owner.finite_controls(source)
    diagnostics=result['report'];previous=owner.accepted['actual_original24_finite_control_diagnostics']
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_full_Z_paired_C1_original24_source_ranges=source['record'],actual_refined_finite_control_diagnostics=diagnostics,
        previous_accepted_full_predicate24_report=dict(filename=base.NAME,sha256=sha(base.NAME)),
        previous_actual_N_scaled_target_C0_Z_ranges=previous['actual_N_scaled_target_C0_Z_ranges'],
        previous_fixed_N_contraction_diagnostic=previous['fixed_N_sufficient_contraction_diagnostic'],
        original_source_construction=construction,original_source_runtime=runtime.record(),
        active_source_cells_freshly_requeried=16,full24_original_C1_integral_range_transport_enclosed=True,
        actual_finite_picard_and_residual_C0_Z_ranges_installed=True,actual_five_controls_installed=False,
        functional_terminal_identity_solved=False,current_whole_N_selected=False,**dict.fromkeys(base.packets.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Complete actual N1024 whole-Z source branch unions refine16 original local ordinary-Z integrals and propagate through unchanged original24 C0/background/P0/geometry/phase. Same original Rc target and finite controls/residuals; no fixed point, global N, terminal, recursion or corrected NS completion.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual full-Z paired source24 finite control diagnostics generated',flush=True);return report


if __name__=='__main__':run()
