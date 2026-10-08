"""Centered branch applications, independent all-N densities and affine targets."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_centered_all_N_C1_controls as new
import lei_ren_part1_paper_compliant_current_original_centered_first_bridge_C1_controls_check as centered_check
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls_check as graph_check

current=new.current;centered=new.centered;ep=new.ep;iv=new.packets.interval
checked=centered_check.checked;restore_native=centered_check.preceding.native_restore
ZERO,DZ=new.ZERO,new.DZ;encode=new.encode


def branch_caps(branch,c,bases,ledger):
    whole=branch['primitives']['alternative']['record'];restore=lambda row:restore_native(row,bases,ledger)
    n=vchecks=overlaps=0
    assert whole['all_original_body_transition_flat_branches_hulled_before_primitive_selection']
    for row in whole['complete_original_cutoff_branches']:
        proof=row['original_centered_first_Z'];n+=1
        if row['name']=='flat':
            assert proof['original_flat_branch_exact_primitive_Z_zero']
            assert all(restore(v).zero for key in ('original_q_C0_Z','original_q_squared_C0_Z') for v in row[key])
            continue
        assert proof['v_positive_lower']=='2' and proof['v_upper']=='3'
        assert proof['fixed_free_angle_derivatives_not_total_inverse_rows']
        assert proof['original_eta_Z_exact_zero'] and proof['no_cutoff_q_floor_or_source_derivative_clipping']
        roots=whole['original_source_root_C0_Z'];DeltaZ=restore(row['original_active_correlation']['selected_same_function_Delta_Z'])
        if row['name']=='negative':
            assert restore(proof['original_v_Z_direct']).zero and restore(proof['original_v_Z_upper']).zero
            assert checked.overlaps(restore(proof['original_aq2_Z_direct']),DeltaZ*(-c.mpf('.5')))
        else:
            assert row['name']=='transition'
            theta=iv(c,row['theta']);s,sp=centered.prior.sigma_jets(c,1-theta)[:2]
            hp=-(s*sp)*2*(2-theta)-s*s
            assert ep(hp)==ep(iv(c,proof['Hcut_prime_wrt_theta']))
            assert checked.overlaps(restore(proof['original_v_Z_direct']),DeltaZ*(1+hp))
            assert checked.overlaps(restore(proof['original_aq2_Z_direct']),DeltaZ*(hp*c.mpf('.5')))
        vchecks+=1
        get=lambda key:restore(proof[key])
        M2Z=get('original_aq2_Z_upper');vZ=get('original_v_Z_upper')
        JZ=get('original_fixed_centered_Jtilde_Z_upper');KZ=get('original_fixed_centered_K_Z_upper');inv=get('original_inverse_contribution_upper')
        assert checked.overlaps(JZ,get('original_fixed_M2H_Z_upper')+get('original_fixed_bqP_Z_upper')+M2Z*c.pi)
        assert checked.overlaps(KZ,JZ*c.mpf('.5')+vZ*(c.pi/4))
        assert checked.overlaps(inv,KZ*(3/c.pi))
        absolute=centered.first.absolute_upper
        AZ=absolute(restore(roots['a'][str(DZ)]));BZ=absolute(restore(roots['b'][str(DZ)]))
        E=absolute(restore(roots['E'][str(ZERO)]));EZ=absolute(restore(roots['E'][str(DZ)]))
        assert checked.overlaps(get('original_total_A_Z_upper'),AZ*c.mpf('.5')+inv)
        assert checked.overlaps(get('original_total_M_Z_upper'),BZ+inv+get('original_fixed_aqP_Z_upper')*(1/c.pi))
        assert checked.overlaps(get('original_total_B_Z_upper'),get('original_total_M_Z_upper')*E*c.mpf('.5')+EZ*c.mpf('1.5'))
        overlaps+=6
    assert n>0
    p=branch['primitives']
    for key in ('A','B_over_Pstar'):assert p['values'][key] is p['old']['values'][key]
    for key in ('A_Z','B_Z_over_Pstar'):
        old,alt=p['old']['values'][key],p['alternative']['values'][key]
        alog=centered.accepted.magnitude_log(alt);olog=centered.accepted.magnitude_log(old)
        change=alog is None and olog is not None or alog is not None and olog is not None and alog<olog
        assert p['decisions'][key]['strict_absolute_upper_reduction']==change
        assert p['values'][key] is (alt if change else old)
    return n,vchecks,overlaps


def route_checks(owner,live,saved):
    c=owner.ctx;coords=owner.target.coordinates;n=v=capchecks=densitychecks=affine=unchanged=0
    incoming={k:{p:coords.scalar(0) for p in new.ORDERS} for k in new.RATES}
    incomingZ={k:dict(row) for k,row in incoming.items()}
    old=owner.previous['actual_all_N_continuous_cell_range_records']
    for index,(cell,row,prior,spec) in enumerate(zip(live['cells'],saved['actual_all_N_continuous_cell_range_records'],old,new.target_module.ROUTE,strict=True)):
        assert (row['label'],row['chart'])==spec[:2]
        assert row['original_geometry']==prior['original_geometry']
        assert row['original_true_mass_applied_once'] and row['quiet_and_pressure_memories_not_reset']
        for b in cell['branches']:
            roots=b['conditional']['source']['roots'];a=roots['a'][ZERO]
            assert b['conditional']['record']['source_family']==owner.family
            nn,vv,cc=branch_caps(b,c,a.scale.bases,a.ledger);n+=nn;v+=vv;capchecks+=cc
            E,EZ,V,VZ=roots['E'][ZERO],roots['E'][DZ],*b['axial']
            A,AZ,B,BZ=[b['primitives']['values'][k] for k in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar')]
            factor=c.exp(c.mpf(0) if A.zero else c.mpf(('-1.25','1.25'))/160)
            F=(A*factor)*E;FZ=(AZ*factor)*E+(A*factor)*EZ;zero=E.scalar(0)
            expected=dict(m={-1:BZ,-2:zero},h={-1:FZ,-2:zero},
                k={-1:(F*VZ+B*EZ)+(FZ*V+BZ*E),-2:B*FZ+F*BZ},
                e={-1:(B*VZ+BZ*V)*2-(F*EZ+FZ*E),-2:BZ*B*2-FZ*F},
                p={-1:FZ*E+F*EZ,-2:FZ*F})
            for key in new.RATES:
                for p in new.ORDERS:
                    assert checked.overlaps(b['density']['Z_derivatives'][key][p],expected[key][p]);densitychecks+=1
        for key in new.RATES:
            for p in new.ORDERS:
                # Actual C0 remains exactly the previous accepted all-N cover.
                checked.same(prior['coefficient_contributions'][key][str(p)],cell['values'][key][p])
                checked.same(prior['coefficient_outgoing'][key][str(p)],cell['outgoing'][key][p]);unchanged+=2
                for field,actual in (('coefficient_contributions',cell['values']),('coefficient_Z_contributions',cell['Z_derivatives']),
                    ('coefficient_outgoing',cell['outgoing']),('coefficient_Z_outgoing',cell['outgoingZ'])):
                    checked.same(row[field][key][str(p)],actual[key][p])
                out=cell['values'][key][p]+incoming[key][p]*cell['factors'][key]['decay']
                outZ=cell['Z_derivatives'][key][p]+incomingZ[key][p]*cell['factors'][key]['decay']
                assert checked.overlaps(out,cell['outgoing'][key][p]) and checked.overlaps(outZ,cell['outgoingZ'][key][p]);affine+=2
        incoming,incomingZ=cell['outgoing'],cell['outgoingZ']
    A,AZ,mu,logA,logmu=live['amplitude'],live['amplitude_Z'],coords.scalar(live['mu']),live['logA'],live['logmu']
    invratio=AZ.positive_divide(A,logA);quotients=0
    for p in new.ORDERS:
        for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
            denom=A if degree==1 else A*A
            val=incoming[key][p].positive_divide(denom,degree*logA)
            jet=(incomingZ[key][p]-incoming[key][p]*invratio*degree).positive_divide(denom,degree*logA)
            assert checked.overlaps(val,live['targets']['values'][row][p]) and checked.overlaps(jet,live['targets']['Z_derivatives'][row][p]);quotients+=2
        numerator=incoming['k'][p]-incoming['m'][p]*A
        numeratorZ=incomingZ['k'][p]-(incoming['m'][p]*AZ+incomingZ['m'][p]*A)
        divisor=A*A*mu;row=new.repair.ROWS[1]
        val=numerator.positive_divide(divisor,2*logA+logmu)
        jet=(numeratorZ-numerator*invratio*2).positive_divide(divisor,2*logA+logmu)
        assert checked.overlaps(val,live['targets']['values'][row][p]) and checked.overlaps(jet,live['targets']['Z_derivatives'][row][p]);quotients+=2
    return dict(passed=True,original24_cells=len(live['cells']),outer_original_branches=live['source_branch_count'],
        complete_centered_inner_branches=n,independent_centered_v_derivative_checks=v,independent_centered_cap_overlaps=capchecks,
        independently_rearranged_density_Z_overlaps=densitychecks,independent_affine_C0_Z_overlaps=affine,
        unchanged_old_all_N_C0_contribution_and_history_rows=unchanged,independent_joint_target_C0_Z_quotients=quotients,
        original_C0_primitives_retained_as_same_objects=True,branch_union_and_true_mass_once=True)


def branch_consistency(saved):
    """Check the saved, fresh-query branch inventory without rerunning a route.

    Reconditioning a formal Delta cover can include endpoint neighbors or
    additional unreachable branches when its microscopic positive lower
    bound is not retained. Such extra covers are conservative; they must
    never be mistaken for an identical predicate partition.
    """
    branch_count=extras=0
    conditions={'negative':'Delta<=0','transition':'0<=Delta<=eta','flat':'Delta>=eta'}
    for cell in saved['actual_all_N_continuous_cell_range_records']:
        if cell['label']=='initial_flat_collar':continue
        query=cell['original_source_query'];branches=cell['original_conditional_branch_results']
        assert query['status']=='enclosed' and query['source_family']==saved['source_family'] and branches
        assert query['cutoff_branch_cover']=='Delta<=0 union 0<=Delta<=eta union Delta>=eta covers R'
        expected=query['conditional_branches']
        assert len(expected)==len(branches)
        for outer,item in zip(expected,branches,strict=True):
            got=item['original_cutoff_branch'];assert got==outer
            name=got['conditional_branch'];assert name in conditions
            assert got['status']=='enclosed' and got['condition']==conditions[name]
            assert got['source_family']==saved['source_family'] and got['branch_boundary_derivatives_not_introduced']
            inner=item['original_full_period_primitive_C1']['original_complete_centered_first_Z']
            assert inner['all_original_body_transition_flat_branches_hulled_before_primitive_selection']
            innernames=[r['name'] for r in inner['complete_original_cutoff_branches']]
            emptynames=[r['name'] for r in inner['branches_proved_empty']]
            assert name in innernames or name in emptynames
            for row in inner['complete_original_cutoff_branches']:
                assert row['condition']==conditions[row['name']]
                if row['name']==name:continue
                assert row['name'] in conditions
                extras+=1
            # The complete same-function hull can be wider than the outer
            # predicate. Original C0 and flat-zero rows are still retained;
            # only one of two valid complete derivative covers is chosen.
            primitive=item['original_full_period_primitive_C1']
            for key in ('A','B_over_Pstar'):
                assert primitive['original_previous_primitive_C0_Z'][key]==primitive['selected_original_primitive_C0_Z'][key]
            if name=='flat':
                for key in ('A_Z','B_Z_over_Pstar'):
                    assert primitive['original_previous_primitive_C0_Z'][key]['exact_zero']
                    assert primitive['selected_original_primitive_C0_Z'][key]['exact_zero']
            branch_count+=1
    assert branch_count==saved['original_conditional_coefficient_branch_count']==35
    return dict(passed=True,enclosed_outer_branch_and_centered_union_bindings=branch_count,
        additional_conservative_inner_branch_covers=extras,inner_cover_not_claimed_identical_to_outer_predicate=True,
        original_flat_derivative_zero_and_C0_preserved=True,complete_body_transition_flat_source_partition_preserved=True)


def accept_branch_consistency_extension(*,previous_checker_sha):
    """Incremental acceptance: only add this new inventory check.

    The prior full replay remains bound by its receipt and every unchanged
    producer/dependency byte. No old check is reclassified as newly run.
    """
    began=time.monotonic();path=new.HERE/new.RECEIPT;receipt=json.loads(path.read_bytes())
    checker=Path(__file__).name
    assert receipt['all_passed'] and receipt[new.GATE] and receipt['input_hashes'][checker]==previous_checker_sha
    for name,digest in receipt['input_hashes'].items():
        if name!=checker:assert new.sha(name)==digest,name
    raw=gzip.decompress((new.HERE/new.NAME).read_bytes())
    assert hashlib.sha256(raw).hexdigest()==receipt['compressed_producer_report']['lossless_original_json_sha256']
    receipt['original_source_branch_consistency']=branch_consistency(json.loads(raw))
    receipt['incremental_branch_inventory_acceptance']=dict(previous_focused_checker_sha256=previous_checker_sha,
        unchanged_full_source_replay_and_algebra_results_reused=True,execution_seconds=time.monotonic()-began)
    receipt['input_hashes'][checker]=new.sha(checker)
    path.write_bytes(json.dumps(receipt,indent=2).encode()+b'\n')
    print('Original all-N source branch consistency extension PASS; unchanged full replay reused',flush=True)
    return receipt


@new.paired.native.inlet.source_precision
def run():
    began=time.monotonic();path=new.HERE/new.NAME;raw=gzip.decompress(path.read_bytes());saved=json.loads(raw)
    flags=('actual_five_controls_installed','certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed',
        'functional_terminal_identity_solved','current_whole_N_selected',*new.packets.OPEN)
    assert saved[new.GATE] and all(saved[k] is False for k in flags)
    for name,digest in saved['input_hashes'].items():assert new.sha(name)==digest,name
    _,owner,built,live=new.run(return_live=True,write=False)
    assert saved['source_family']==owner.family and saved['original_uniform_integer_N_lower']==160
    old=owner.previous
    assert built['graph'].nodes==old['exact_function_graph_nodes']==saved['exact_function_graph_nodes']
    for key,actual in (('exact_coefficient_history_roots',current.encode_orders(built['coefficient_history'])),
        ('exact_coefficient_target_roots',current.encode_orders(built['coefficient_target_orders'])),
        ('exact_reconstructed_N_scaled_target_roots',new.controls.pair_roots(built['N_scaled_targets'].values(),built['N_scaled_targets'].keys()))):
        assert saved[key]==old[key]==actual
    checks=dict(centered_source_algebra=centered_check.identities(),
        original_seam_binding=centered_check.preceding.derivative_identities(owner.ctx),
        original_source_branch_consistency=branch_consistency(saved),
        fresh_original24_branch_density_affine_target=route_checks(owner,live,saved),
        actual_N_phase_and_coefficient_integrals=graph_check.original_phase_and_live_ranges(owner,built,live))
    # Parent exact-function/scalar suites were already accepted; graph bytes,
    # roots and uniform_coefficients are unchanged, so do not repeat them.
    checks['unchanged_accepted_exact_all_N_function_graph']=dict(passed=True,node_count=len(built['graph'].nodes),
        accepted_parent_density_C1_and_scalar_uniform_suites_reused=True,exact_phase_shared_N_and_roots_unchanged=True)
    conditions=owner.frequency_conditions(live)
    assert encode(conditions)==saved['original_source_and_repair_combined_log_conditions']
    assert conditions['remaining_global_cone_higher_jet_join_N_conditions_not_certified'] and not conditions['current_whole_N_selected']
    rejections=0
    for fn in (lambda:owner.frequency_conditions(None),lambda:owner.frequency_conditions(dict(live)),lambda:owner.route(Z=(0,1))):
        try:fn()
        except (ValueError,TypeError):rejections+=1
        else:raise AssertionError('Unissued/all-Z source guard failed')
    assert rejections==3
    checks['issued_source_and_global_N_scope_guards']=dict(passed=True,rejections=rejections,
        frequency_is_conditional_source_repair_lower_bound_only=True)
    receipt=dict(all_passed=True,**{new.GATE:True},source_family=owner.family,**checks,
        strict_all_N_target_C1_reductions=sum(row['strict_log_upper_reduction'] for row in new.comparison(owner,live).values()),
        compressed_producer_report=dict(filename=new.NAME,compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw),
            lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        **dict.fromkeys(flags,False),read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        input_hashes={**saved['input_hashes'],new.NAME:new.sha(new.NAME),Path(__file__).name:new.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Independent centered branch algebra/applications, all-N density/affine/joint C1 ranges, preserved original C0/exact phase/function graph, and conditional frequency guard. No global N/actual control/terminal/recursion admission.')
    (new.HERE/new.RECEIPT).write_bytes(json.dumps(encode(receipt),indent=2).encode()+b'\n')
    print('Original centered all-N C1 source/frequency focused acceptance PASS',flush=True)
    return receipt


if __name__=='__main__':run()
