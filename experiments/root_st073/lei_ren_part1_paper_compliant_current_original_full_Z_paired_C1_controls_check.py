"""Focused original-source first-Z, complete branch and actual24 acceptance.

Accepted ancestors are reused. Saved complete typed source covers are
replayed as enclosures, never as numerical source owners or field values.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_full_Z_paired_C1_controls as current
import lei_ren_part1_paper_compliant_current_original_full_predicate_24_cell_controls_check as checked

base=current.base;ep=current.ep;iv=base.packets.interval;KEYS=base.KEYS
ZERO,DZ=current.ZERO,current.DZ;prior=current.prior;first=current.first


def encoded(value):return json.loads(json.dumps(current.encode(value)))


def scalar(c,record):
    return c.mpf(mp.make_mpf(tuple(record['exact_mpf_tuple']))) if isinstance(record,dict) else c.mpf(record)


def native_restore(record,bases,ledger):
    s=record['formal_positive_scale'];c=bases[0].ctx
    value=prior.ScaledEnclosure(prior.FormalScale(bases,tuple(s['source_exponents'])+(s['radius_power'],),
        iv(c,s['additional_log_interval'])),iv(c,record['coefficient_interval']),ledger)
    checked.same(record,value)
    return value


def upper(value):
    if value.zero:return None
    return ep(value.scale.evaluate()+value.ctx.ln(value.ctx.mpf(max(abs(v) for v in ep(value.coefficient)))))[1]


def selection(old,new,record):
    old.coerce(new);a,b=upper(old),upper(new)
    improved=b is None and a is not None or a is not None and b is not None and b<a
    assert record['strict_absolute_upper_reduction']==improved
    for name,value in (('old',a),('alternative',b),('selected',b if improved else a)):
        saved=record[name+'_absolute_upper_log']
        assert saved is None if value is None else ep(scalar(old.ctx,saved))==(value,value)
    assert record['selected_complete_source_function_cover_not_endpoint_value']
    return new if improved else old


def symmetric(value):
    lo,hi=ep(value.coefficient)
    assert lo==-hi, 'Primitive derivative choice requires complete symmetric source covers'


def independent_cutoff_first_Z(roots,branch,eta_log,logamin):
    """Ordinary first derivatives only; no division by the cutoff q."""
    a,az=roots['a'][ZERO],roots['a'][DZ];c=a.ctx;zero=a.scalar(0)
    if branch['name']=='flat':return {ZERO:zero,DZ:zero},{ZERO:zero,DZ:zero}
    eta=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=eta_log),1,a.ledger)
    deltaZ=roots['kappa_minus2'][DZ]
    gamma_lower=eta_log+(c.ln(2) if branch['name']=='negative' else 0)
    gamma=(eta*2-branch['Delta']).positive_intersection(gamma_lower)
    ratio=gamma.positive_divide(a*2,logamin+c.ln(2))
    logroot=(gamma_lower-c.mpf(upper(a))-c.ln(2))/2
    ratio=ratio.positive_intersection(logroot*2)
    root=first.current.nonnegative_sqrt(ratio).positive_intersection(logroot)
    ratioZ=(-deltaZ-ratio*az*2).positive_divide(a*2,logamin+c.ln(2))
    rootZ=ratioZ.positive_divide(root*2,logroot+c.ln(2))
    if branch['name']=='negative':return {ZERO:root,DZ:rootZ},{ZERO:ratio,DZ:ratioZ}
    theta=branch['theta'];s,sp=prior.sigma_jets(c,1-theta)[:2]
    sigmaZ=(-deltaZ).positive_divide(eta,eta_log)*sp
    q={ZERO:root*s,DZ:rootZ*s+root*sigmaZ}
    q2=ratio*(s*s)
    Hprime=-2*s*sp*(2-theta)-s*s
    q2Z=(deltaZ*Hprime-q2*az*2).positive_divide(a*2,logamin+c.ln(2))
    return q,{ZERO:q2,DZ:q2Z}


def derivative_identities(c):
    z=sy.Symbol('Z');a,b=sy.Function('a')(z),sy.Function('b')(z);t=-b/a
    assert sy.simplify(sy.diff(a+b*b/a-2,z)-((1-t*t)*sy.diff(a,z)-2*t*sy.diff(b,z)))==0
    theta=sy.Function('theta')(z);eta=sy.Symbol('eta',positive=True);s=sy.Function('sigma')(1-theta)
    R=eta*(2-theta)/(2*a);Q=s*s*R
    Hprime=-2*s*sy.Subs(sy.Derivative(sy.Function('sigma')(sy.Symbol('x')),sy.Symbol('x')),sy.Symbol('x'),1-theta)*(2-theta)-s*s
    assert sy.simplify(sy.diff(Q,z)-(Hprime*eta*sy.diff(theta,z)/(2*a)-Q*sy.diff(a,z)/a))==0
    # Bind the accepted all-order seam theorem and check the admitted jets
    # at both endpoints. The unmodulated root has gamma=2eta at Delta=0.
    theorem=first.original_periodic_theorem()
    assert theorem['cutoff_smooth_seams']=='sigma and all positive derivatives are flat at0 and1; the active/flat union defines the same original smooth q/A/B functions'
    for endpoint in (0,1):
        rows=prior.sigma_jets(c,c.mpf(endpoint))
        assert ep(rows[0])==(endpoint,endpoint) and all(ep(row)==(0,0) for row in rows[1:])
    return dict(passed=True,independent_Delta_Z_and_eta_collected_q_squared_Z_identities=2,
        original_eta_Z_exact_zero=True,cutoff_q_not_used_as_positive_denominator=True,
        accepted_all_order_sigma_seam_theorem=theorem['cutoff_smooth_seams'],exact_sigma_endpoint_jet_checks=10,
        body_transition_root_positive='Delta=0 gives gamma=2eta>0 and a>0; sigma(1)=1, sigma_Z=0',
        transition_flat_C1_gluing='Delta=eta gives gamma=eta>0; sigma(0)=sigma_prime(0)=0; q/q_Z/q2/q2_Z vanish')


def refine_saved_sources(owner,live,saved):
    proofs=saved['actual_full_Z_paired_C1_refinement_records'];c=owner.ctx;coords=owner.coordinates
    assert [r['original_cell_index'] for r in proofs]==list(current.ACTIVE_INDICES)
    root_owner=owner.target_owner.q_owner.owner.owner
    firstgeom=owner.target_owner.geometry(*base.parameters.ROUTE[1])
    fresh=owner.target_owner.q_owner.query('bridge_first',(-1,1),firstgeom['coordinate'])
    freshroots=fresh['source']['roots'];ledger_template=dict(freshroots['a'][ZERO].ledger)
    primitive_count=branch_count=q_count=density_count=mass_count=local_count=reductions=0
    branch_inventory=[]
    for proof in proofs:
        index=proof['original_cell_index'];cell=live['cells'][index];spec=base.parameters.ROUTE[index]
        label,chart,left,right=spec;alt=proof['actual_complete_paired_cutoff_C1'];geom=owner.target_owner.geometry(*spec)
        assert (proof['label'],proof['chart'])==(label,chart)
        assert proof['source_family']==owner.family and proof['candidate_N']==1024 and ep(iv(c,proof['Z_box']))==(-1,1)
        assert encoded(geom['record'])==proof['actual_original_geometry']==cell['record']['actual_original_geometry']
        assert geom['record']['width_and_endpoints_independent_of_Z']
        positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
        assert encoded(positive)==proof['original_chart_positive_a_theorem']
        positivity='whole_actual_source_positive_not_inferred_from_saved_denominator_box' if chart.startswith('O3_') else 'source_function_positivity_not_inferred_from_saved_box'
        assert positive[positivity]
        bases=tuple(iv(c,x) for x in alt['original_native_log_bases']);ledger=dict(ledger_template)
        assert len(bases)==5 and alt['source_context_dps']==c.dps
        signed=owner.target_owner.transfer.owner.signed_owner
        packet=signed.packet(chart,(-1,1),geom['coordinate'])
        packet_bases=tuple(packet.algebra.logs)+(packet.provenance['logR_cover'],)
        assert all(ep(x)==ep(y) for x,y in zip(bases,packet_bases,strict=True))
        restore=lambda r:native_restore(r,bases,ledger)
        roots={name:{order:restore(rows[str(order)]) for order in (ZERO,DZ)} for name,rows in alt['original_source_root_C0_Z'].items()}
        query=proof['fresh_original_source_query']
        assert query['source_family']==owner.family and query['chart']==chart
        provider=query['original_correlated_shear_and_q']['correlated_shear_and_signed_root_enclosures']
        for name,rows in roots.items():
            for order,value in rows.items():checked.same(provider[name]['y%d_Z%d'%order],value)
        if index==1:
            for name,rows in roots.items():
                for order,value in rows.items():checked.same(alt['original_source_root_C0_Z'][name][str(order)],freshroots[name][order])
            assert all(ep(x)==ep(y) for x,y in zip(bases,freshroots['a'][ZERO].scale.bases,strict=True))
        eta,logamin,dstar=[iv(c,alt[k]) for k in ('eta_log','log_a_positive_lower','dstar_log')]
        assert ep(logamin)==ep(positive['log_actual_a_positive_lower']) and ep(eta)[1]<=ep(-c.ln(2))[0]
        candidates,empty=current.cutoff.conditional_cutoff_branches(roots['kappa_minus2'][ZERO],eta)
        rows=alt['actual_conditional_branches'];byname={r['name']:r for r in rows};covers=[]
        assert len(byname)==len(rows)
        for branch in candidates:
            name=branch['name'];restricted={name:dict(row) for name,row in roots.items()}
            restricted['kappa_minus2'][ZERO]=branch['Delta']
            if name!='flat':
                limit=2 if name=='negative' else c.mpf('2.5')
                a=first.positive_restriction(roots['a'][ZERO],logamin,c.ln(limit))
                if a is None:
                    assert name not in byname
                    empty.append(dict(name=name,proof=dict(active_branch_empty_by_original_kappa_ge_a=True,active_a_upper=limit)))
                    continue
                assert name in byname;corr=byname[name]['same_original_active_correlation']
                b,_=current.q2.accepted.signed_C0_support_range(roots['b'][ZERO],a.scalar(c.mpf(('-1.25','1.25'))))
                derived=(-b).positive_divide(a,logamin)
                supported=prior.ScaledEnclosure(prior.FormalScale(bases,offset=(c.ln(3)-logamin)/2),1,ledger)
                cap=current.serial.minimum_upper(derived,roots['t0'][ZERO],supported);t=first.symmetric_bound(cap)
                DeltaZ=(a.scalar(1)-first.current.square(t))*roots['a'][DZ]-(t*roots['b'][DZ])*2
                checked.same(corr['original_conditional_a_C0'],a);checked.same(corr['original_conditional_b_C0'],b)
                checked.same(corr['original_active_t0_absolute_support'],cap);checked.same(corr['original_correlated_Delta_Z'],DeltaZ)
                new=first.symmetric_bound(DeltaZ);symmetric(new)
                selected=selection(roots['kappa_minus2'][DZ],new,corr['Delta_Z_comparison'])
                checked.same(corr['selected_same_function_Delta_Z'],selected)
                restricted['a'][ZERO]=a;restricted['b'][ZERO]=b;restricted['t0'][ZERO]=t;restricted['kappa_minus2'][DZ]=selected
            assert name in byname;row=byname[name];branch_count+=1
            assert row['condition']==branch['condition'];checked.same(row['conditional_Delta_C0'],branch['Delta'])
            if name=='transition':assert ep(iv(c,row['theta']))==ep(branch['theta'])
            q={ZERO:restore(row['q_C0']),DZ:restore(row['q_Z'])}
            q2={ZERO:restore(row['direct_q_squared_C0']),DZ:restore(row['direct_q_squared_Z'])}
            iq,ir=independent_cutoff_first_Z(restricted,branch,eta,logamin)
            for actual,independent in ((q,iq),(q2,ir)):
                for order in (ZERO,DZ):assert checked.overlaps(actual[order],independent[order]);q_count+=1
            if name=='flat':assert all(v.zero for v in (*q.values(),*q2.values()))
            support=current.paired.paired_C1_support(restricted,q,q2,logamin,dstar)
            for key,value in support['values'].items():checked.same(row['derivative_covers'][key],value);symmetric(value)
            covers.append(support['values'])
        assert set(byname)=={branch['name'] for branch in candidates}-{row['name'] for row in empty}
        assert encoded(empty)==alt['branches_proved_empty']
        assert alt['complete_branch_union_precedes_whole_source_cover_selection']
        hull={key:current.density.local.same_source_union([r[key] for r in covers]) for key in ('A_Z','B_Z_over_Pstar')}
        old=(first.whole_period_C1 if index==1 else current.serial.whole_period_C1)(roots,eta,logamin,dstar)
        chosen=dict(old['values'])
        for key in hull:
            checked.same(alt['whole_source_derivative_hulls'][key],hull[key]);symmetric(hull[key]);symmetric(old['values'][key])
            chosen[key]=selection(old['values'][key],hull[key],proof['primitive_derivative_comparisons'][key]);primitive_count+=1
        assert chosen['A'] is old['values']['A'] and chosen['B_over_Pstar'] is old['values']['B_over_Pstar']
        for key,value in chosen.items():checked.same(proof['selected_original_primitive_C0_Z'][key],value)
        velocity={key:restore(value) for key,value in proof['original_density_velocity_and_increment_source_C0_Z'].items()}
        checked.same(proof['original_density_velocity_and_increment_source_C0_Z']['original_E'],roots['E'][ZERO])
        checked.same(proof['original_density_velocity_and_increment_source_C0_Z']['original_E_Z'],roots['E'][DZ])
        for k,name in ((0,'original_V'),(1,'original_V_Z')):
            row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
            value=signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),bases,ledger)
            checked.same(proof['original_density_velocity_and_increment_source_C0_Z'][name],value)
        got=current.density.density_Z_kernels(velocity['original_E'],velocity['original_E_Z'],velocity['original_V'],velocity['original_V_Z'],chosen,1024)
        selected_jets={}
        for key,rate in base.history.RATES.items():
            checked.same(proof['original_signed_density_Z_covers'][key],got['Z_derivatives'][key]);density_count+=1
            factors=base.preceding.downstream.transfer.true_width_kernel(coords,geom,rate)
            for field in ('mass','decay'):checked.same(proof['original_true_width_kernel_factors'][key][field],factors[field]);mass_count+=1
            contribution=coords.rebase(got['Z_derivatives'][key],owner.family)*factors['mass']
            checked.same(proof['paired_local_Z_contributions'][key],contribution)
            checked.same(proof['accepted_local_Z_contributions'][key],cell['Z_derivatives'][key])
            chosen_local=selection(cell['Z_derivatives'][key],contribution,proof['local_derivative_comparisons'][key])
            checked.same(proof['selected_local_Z_contributions'][key],chosen_local)
            checked.same(proof['accepted_C0_local_contributions_unchanged'][key],cell['values'][key])
            selected_jets[key]=chosen_local;local_count+=1;reductions+=int(proof['local_derivative_comparisons'][key]['strict_absolute_upper_reduction'])
        cell['Z_derivatives']=selected_jets;cell['operator'].Z_increments=selected_jets
        branch_inventory.append(dict(label=label,admitted_branches=list(byname),proved_empty_branches=[r['name'] for r in empty]))
    return dict(passed=True,fresh_full_first_bridge_source_anchor=True,fresh_original_native_packet_basis_and_axial_velocity_anchors=16,source_cells=16,
        original_whole_Z_branch_inventory=branch_inventory,conditional_branch_replays=branch_count,
        independent_linear_q_and_direct_q_squared_C0_Z_overlaps=q_count,complete_primitive_Z_cover_choices=primitive_count,
        original_signed_density_Z_replays=density_count,original_mass_and_decay_comparisons=mass_count,
        complete_local_Z_cover_choices=local_count,strict_local_Z_absolute_upper_reductions=reductions,
        original_C0_primitives_and_local_integrals_retained=True,only_first_Z_covers_refined=True)


def propagate(owner,live,saved,old):
    incoming=dict(values={key:owner.coordinates.scalar(0) for key in KEYS},Z_derivatives={key:owner.coordinates.scalar(0) for key in KEYS})
    count=0
    for cell,row,previous in zip(live['cells'],saved['actual_original24_source_cell_records'],old['actual_original24_source_cell_records'],strict=True):
        op=cell['operator'];out=op.apply(incoming['values'],incoming['Z_derivatives'],owner.family)
        for key in KEYS:
            for suffix,field in (('C0','values'),('Z','Z_derivatives')):
                checked.same(row['actual_inherited_correction_'+suffix][key],incoming[field][key])
                checked.same(row['actual_local_'+suffix+'_contributions'][key],cell[field][key])
                independent=incoming[field][key]*op.coefficients[key]+cell[field][key]
                checked.same(row['actual_right_correction_'+suffix][key],out[field][key]);assert checked.overlaps(independent,out[field][key]);count+=1
                checked.same(row['actual_right_own_history_'+suffix][key],cell['background'][field][key]+out[field][key])
            assert row['actual_right_correction_C0'][key]==previous['actual_right_correction_C0'][key]
            assert row['actual_local_C0_contributions'][key]==previous['actual_local_C0_contributions'][key]
        for key in ('actual_original_geometry','actual_original_source_cell_binding','actual_right_background_C0','actual_right_background_Z',
            'original_separate_P0','original_separate_P0_Z','absolute_pressure_C0'):
            assert row[key]==previous[key],key
        bg=cell['background'];checked.same(row['absolute_pressure_Z'],bg['P0_Z']+(bg['Z_derivatives']['p']+out['Z_derivatives']['p']))
        cell['incoming']=incoming;cell['correction']=out;cell['record']=row;incoming=out
    for j,field in enumerate(('values','Z_derivatives')):
        for key in KEYS:checked.same(saved['actual_Rc_correction_C0_Z'][j][key],incoming[field][key])
    a=live['amplitude']
    target=dict(**base.fixed.fixed_N_target_rows(incoming['values'],incoming['Z_derivatives'],a['A'],a['AZ'],a['logA'],a['mu_source'],a['logmu']),
        **{key:a[key] for key in ('A','AZ','mu','logA','logmu')})
    live.update(history=incoming['values'],Z_derivatives=incoming['Z_derivatives'],target=target,
        record=dict(saved,Z_box=iv(owner.ctx,saved['Z_box'])))
    result=checked.joint_target(owner,live,saved)
    result.update(independent_actual24_C0_Z_affine_checks=count,all_original_C0_operator_background_P0_geometry_bindings_unchanged=True)
    return result


def finite_controls(owner,live,result,saved,previous):
    assert encoded(result['report'])==saved,checked.first_difference(encoded(result['report']),saved)
    for key in ('original_exact_weight_enclosures','original_exact_matrix_enclosures','exact_control_graph_nodes','exact_original_source_roles',
        'exact_function_source_graph_sha256','exact_finite_control_function_roots','exact_repair_residual_function_roots'):
        assert saved[key]==previous[key],key
    built=result['live'];matrix=built['exact_matrix_enclosures'];c=owner.ctx;d=built['N_scaled_target_ranges'];zero=d[0].value.scalar(0)
    updates=residuals=quadratics=0
    for depth,row in enumerate(built['ranges']):
        Q=checked.independent_Q(c,live['target'],matrix,row)
        production=base.controls.range_quadratic(c,live['target']['mu'],live['target']['logmu'],matrix,row)
        for j,label in enumerate(base.controls.ROWS):
            for field,suffix in (('value','C0'),('Z','Z')):
                assert checked.overlaps(getattr(Q[j],field),getattr(production[j],field));quadratics+=1
                residual=getattr(d[j],field)+getattr(Q[j],field)*(c.mpf(1)/1024)
                residual+=sum((getattr(row[k],field)*matrix['linear_enclosure'][j][k] for k in reversed(range(5))),zero)
                assert checked.overlaps(residual,owner.restore(saved['actual_finite_control_and_residual_ranges'][depth]['actual_repair_equation_residual_C0_Z_ranges'][label][suffix]));residuals+=1
                if depth<saved['finite_picard_depth']:
                    rhs=[getattr(d[k],field)+getattr(Q[k],field)*(c.mpf(1)/1024) for k in range(5)]
                    nxt=-sum((rhs[k]*matrix['inverse_enclosure'][j][k] for k in reversed(range(5))),zero)
                    assert checked.overlaps(nxt,getattr(built['ranges'][depth+1][j],field));updates+=1
    diagnostic=base.controls.fixed_N_contraction_diagnostic(built)
    assert encoded(diagnostic)==saved['fixed_N_sufficient_contraction_diagnostic']
    assert not saved['certified_fixed_point_tail_installed'] and not diagnostic['globally_compatible_N_admitted']
    return dict(passed=True,exact_original_source_graph_and_finite_function_definitions_unchanged=True,
        independent_quadratic_C0_Z_overlaps=quadratics,independent_finite_Picard_C0_Z_update_overlaps=updates,
        actual_linear_B_numerical_residual_C0_Z_overlaps=residuals,finite_depth=saved['finite_picard_depth'],
        sufficient_half_contraction_bound_passed=diagnostic['half_contraction_sufficient_bound_passed'],finite_ranges_not_fixed_point=True)


@base.native.inlet.source_precision
def run():
    began=time.monotonic();path=base.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    flags=('actual_five_controls_installed','functional_terminal_identity_solved','current_whole_N_selected',*base.packets.OPEN)
    assert manifest[current.GATE] and manifest['candidate_N']==1024 and all(manifest[k] is False for k in flags)
    for name,digest in manifest['input_hashes'].items():assert base.sha(name)==digest,name
    previous=json.loads(gzip.decompress((base.HERE/base.NAME).read_bytes()))
    assert manifest['previous_accepted_full_predicate24_report']['sha256']==base.sha(base.NAME)
    bridge,_=base.native.inlet.native_bridge_owner();checks={}
    with base.native.inlet.CheckedSourceRuntime():
        service=current.OriginalFullZPairedC1Controls(bridge);owner=service.base;live=owner.assemble()
        assert owner.family==manifest['source_family']
        with mp.workdps(owner.ctx.dps+40):
            checks['source_derivative_identities_and_C1_seams']=derivative_identities(owner.ctx)
            saved=manifest['actual_full_Z_paired_C1_original24_source_ranges']
            checks['complete_whole_Z_paired_source_refinement']=refine_saved_sources(owner,live,saved)
            print('Full-Z original first derivative and complete source branch replay PASS',flush=True)
            checks['actual24_affine_and_new_joint_target']=propagate(owner,live,saved,previous['actual_original24_full_predicate_source_range_bridge'])
            print('Unchanged C0 and actual24 refined Z propagation PASS',flush=True)
        result=owner.finite_controls(live)
        with mp.workdps(owner.ctx.dps+40):
            checks['same_exact_functions_actual_finite_controls_and_residuals']=finite_controls(owner,live,result,
                manifest['actual_refined_finite_control_diagnostics'],previous['actual_original24_finite_control_diagnostics'])
            count=0
            for callback in (lambda:service.refine(N=True),lambda:service.refine(N=2048),lambda:service.finite_controls(None),lambda:service.finite_controls(dict(live))):
                try:callback()
                except (TypeError,ValueError):count+=1
                else:raise AssertionError('Invalid frequency or unissued refined source admitted')
            assert count==4;checks['invalid_N_and_unissued_source_guards']=dict(passed=True,rejections=count)
            print('Same original finite control and actual numerical residual replay PASS',flush=True)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        full24_original_C1_integral_range_transport_enclosed=True,actual_finite_picard_and_residual_C0_Z_ranges_installed=True,
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:base.sha(current.NAME),Path(__file__).name:base.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original full-Z source first-derivative refinement with complete branch unions, unchanged C0/background/P0/geometry, original24 affine propagation, new joint targets and same exact finite control functions. Ancestor suites reused. No global N, fixed point, terminal, heat/cone/recursion/pulse/corrected NS completion.')
    (base.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(receipt),indent=2).encode()+b'\n')
    print('Actual full-Z paired original24 controls focused acceptance PASS',flush=True)
    return receipt


if __name__=='__main__':run()
