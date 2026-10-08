"""Use the original active kappa relation before conditional q/q² Z jets."""
from contextlib import contextmanager
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_upstream_paired_branch_C1 as previous

baseline=previous.baseline;common=previous.common;middle=previous.middle
cutoff=previous.cutoff;q2=previous.q2;paired=previous.paired;prior=previous.prior
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_native_upstream_active_kappa_C1.json'
RECEIPT=PREFIX+'current_native_upstream_active_kappa_C1_check.json'
REPLAY=PREFIX+'current_original_O2_active_kappa_upstream_replay.json'
GATE='original_active_kappa_conditional_Delta_Z_and_actual_upstream_C1_histories_executed'
N=previous.N;KEYS=previous.KEYS;LABELS=previous.LABELS
ZERO=previous.ZERO;DZ=previous.DZ;ep=previous.ep;encode=previous.encode
BASE_BRANCH_COVER=previous.branch_paired_cover


def active_kappa_theorem():
    a,b,az,bz=sy.symbols('a b a_Z b_Z',real=True)
    k=a+b*b/a
    if sy.cancel(sy.diff(k,a)*az+sy.diff(k,b)*bz-(az*(1-b*b/(a*a))+2*b*bz/a))!=0:
        raise ArithmeticError('Original kappa ordinary-Z identity changed')
    return dict(passed=True,original_identity='Delta_Z=a_Z*(1-b²/a²)+2*b*b_Z/a',
        original_source_relation='kappa=a+b²/a; Delta=kappa-2; a>0',
        conditional_relation='kappa<=K => a<=K, b²/a<=K, |b/a|<=sqrt(K/a)',
        same_derivative_majorant='|Delta_Z|<=|a_Z|*(1+K/a)+2*|b_Z|*sqrt(K/a)',
        body_K=2,transition_K='2+eta<=5/2',eta_Z_exact_zero_original_global_constant=True,
        bound_valid_only_on_nonflat_conditional_branch=True,
        original_a_Z_b_Z_retained_not_derivatives_of_range_selectors=True,
        sigma_flatness_preserves_Delta_zero_and_eta_seam_jets=True)


def active_root_cover(roots,log_a_lower,branch):
    """Narrow covers of the same a/Delta_Z, never their defining functions."""
    if branch['name']=='flat':return roots,dict(flat_branch_no_active_constraint_or_body_evaluation=True)
    if branch['name'] not in ('negative','transition'):raise ValueError('Original nonflat cutoff branch required')
    a=roots['a'][ZERO];c=a.ctx;K=c.mpf(2 if branch['name']=='negative' else '2.5')
    bounded_a=middle.p.positive_restriction(a,log_a_lower,c.ln(K))
    if bounded_a is None:raise ArithmeticError('Conditional a positivity/support inconsistent; certify empty branch explicitly')
    ratio=a.scalar(K).positive_divide(bounded_a,log_a_lower)
    root=q2.current.nonnegative_sqrt(ratio)
    az=baseline.serial.absolute_upper(roots['a'][DZ]);bz=baseline.serial.absolute_upper(roots['b'][DZ])
    bound=baseline.serial.symmetric_bound(az*(a.scalar(1)+ratio)+bz*root*2)
    before=baseline.magnitude_log(roots['kappa_minus2'][DZ]);after=baseline.magnitude_log(bound)
    tighter=after is None or before is not None and after<before
    selected=bound if tighter else roots['kappa_minus2'][DZ]
    result=dict(roots,a=dict(roots['a']),kappa_minus2=dict(roots['kappa_minus2']))
    result['a'][ZERO]=bounded_a;result['kappa_minus2'][DZ]=selected
    if result['a'][DZ] is not roots['a'][DZ] or result['b'][DZ] is not roots['b'][DZ]:
        raise ValueError('Original source a_Z/b_Z objects must survive')
    return result,dict(original_active_kappa_theorem=active_kappa_theorem(),conditional_K=K,
        original_a_C0=a.record(),conditional_same_source_positive_a_C0=bounded_a.record(),
        original_a_Z=roots['a'][DZ].record(),original_b_Z=roots['b'][DZ].record(),
        K_over_a=ratio.record(),sqrt_K_over_a=root.record(),
        inherited_original_Delta_Z=roots['kappa_minus2'][DZ].record(),
        same_original_derivative_conditional_majorant=bound.record(),selected_Delta_Z_cover=selected.record(),
        inherited_absolute_upper_log=before,conditional_majorant_absolute_upper_log=after,
        strict_Delta_Z_upper_reduction=tighter and before!=after,
        original_a_Z_b_Z_and_remaining_mixed_root_rows_retained=True,
        only_same_source_a_C0_and_Delta_Z_enclosure_coordinates_refined=True,
        root_map_copied_original_roots_not_mutated=True)


def active_branch_cover(roots,eta_log,log_a_lower,dstar_log):
    original=BASE_BRANCH_COVER(roots,eta_log,log_a_lower,dstar_log)
    a=roots['a'][ZERO];c=a.ctx
    if ep(eta_log)[1]>ep(-c.ln(2))[0]:raise ValueError('Original eta<=1/2 required')
    candidates,empty=cutoff.conditional_cutoff_branches(roots['kappa_minus2'][ZERO],eta_log)
    branches=[];values=[]
    for branch in candidates:
        restricted=dict(roots,kappa_minus2=dict(roots['kappa_minus2']))
        restricted['kappa_minus2'][ZERO]=branch['Delta']
        conditional,proof=active_root_cover(restricted,log_a_lower,branch)
        q=cutoff.conditional_q_jet(conditional,eta_log,log_a_lower,branch)
        square=q2.conditional_q2_jet(conditional,eta_log,log_a_lower,branch)
        got=paired.paired_C1_support(conditional,q['rows'],square,log_a_lower,dstar_log)
        values.append(got['values'])
        branches.append(dict(name=branch['name'],condition=branch['condition'],theta=branch.get('theta'),
            original_conditional_Delta_C0=branch['Delta'].record(),same_source_active_root_refinement=proof,
            original_q_jet={('y%d_Z%d'%order):value.record() for order,value in q['rows'].items()},
            original_direct_q_squared_jet=q2.q2_record(square,branch),
            original_q_branch_derivation={k:v for k,v in q.items() if k!='rows'},
            original_paired_whole_period_support=got['record'],
            paired_derivative_covers={k:v.record() for k,v in got['values'].items()},
            original_a_Z_b_Z_and_remaining_mixed_root_source_rows_not_changed=True))
    union=middle.density.local.same_source_union
    hull={k:union([row[k] for row in values]) for k in ('A_Z','B_Z_over_Pstar')}
    selected={};comparisons={}
    for k,value in hull.items():
        before=baseline.magnitude_log(original['values'][k]);after=baseline.magnitude_log(value)
        tighter=after is None or before is not None and after<before
        selected[k]=value if tighter else original['values'][k]
        comparisons[k]=dict(accepted_paired_absolute_upper_log=before,active_kappa_absolute_upper_log=after,
            strict_upper_reduction=tighter and before!=after,
            selected_backend='active_kappa_branch_hull' if tighter else 'accepted_paired_branch_hull')
    record=dict(original['record'],conditional_cutoff_paired_branches=branches,
        whole_source_paired_derivative_hulls={k:v.record() for k,v in selected.items()},
        new_active_kappa_branch_derivative_hulls={k:v.record() for k,v in hull.items()},
        accepted_paired_cover_on_same_original_roots=original['record'],active_kappa_comparisons=comparisons,
        actual_Delta_Z_ranges_refined_from_original_identity_not_selector_derivatives=True)
    return dict(values=selected,record=record)


@contextmanager
def active_firstbridge_primitives():
    original=previous.branch_paired_cover;previous.branch_paired_cover=active_branch_cover
    try:
        with previous.firstbridge_paired_primitives() as traces:yield traces
    finally:
        changed=previous.branch_paired_cover is not active_branch_cover
        previous.branch_paired_cover=original
        if changed:raise RuntimeError('Active kappa primitive runtime binding changed unexpectedly')


@middle.native.inlet.source_precision
def run():
    begin=time.monotonic();bridge,seed=middle.native.inlet.native_bridge_owner()
    with middle.native.inlet.CheckedSourceRuntime():
        owner=baseline.build_owner(bridge);hashes,binding=baseline.accepted_hashes(owner.family)
        common.attach_receipt(hashes,previous,owner.family)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Original source closures disagree: '+name)
            hashes[name]=digest
        archives=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            with active_firstbridge_primitives() as traces:route=owner.route(Z,N=N)
            if len(traces)!=12:raise ValueError('All original twelve active chart calls required')
            cells=baseline.live_active_cells(route)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                complete_actual_pre_O2_tile_route=route['record'],exact_upstream_source_function_incoming_binding=binding,
                twelve_source_primitive_cover_records=traces,
                per_chart_actual_true_width_contribution_attribution={label:dict(C0={k:v.record() for k,v in cell['contributions'].items()},
                    Z={k:v.record() for k,v in cell['Z_derivatives'].items()}) for label,cell in cells},
                full_source_bounds_requeried_on_this_tile_not_relabelled_whole_Z=True,
                original_source_normalization_densities_masses_and_incoming_equations_unchanged=True,
                only_firstbridge_ordinary_Z_primitive_cover_tightened=True)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_native_upstream_active_kappa_C1_'+tag+'.json.gz'
            if len(compressed)>=100*1024*1024:raise ArithmeticError('Split native evidence losslessly')
            (HERE/filename).write_bytes(compressed)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()))
            print('Genuine active kappa source derivative histories:',tag,flush=True)
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            genuine_active_kappa_upstream_archives=archives,original_active_kappa_theorem=active_kappa_theorem(),
            exact_upstream_source_function_incoming_binding=binding,source_seed_evidence=seed,
            original_a_Z_b_Z_and_C0_densities_background_P0_native_owner_retained=True,
            full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False),
            input_hashes=hashes,execution_seconds=time.monotonic()-begin,
            scope='Genuine original13-cell upstream same-N1024 strict-sign histories, firstbridge source-equivalent conditional a/Delta_Z range refinement from original active kappa relation. Original a_Z/b_Z, other mixed roots, C0 A/B/densities, native owner, widths and pressure datum unchanged. No full atlas/Rc matching/repair/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
