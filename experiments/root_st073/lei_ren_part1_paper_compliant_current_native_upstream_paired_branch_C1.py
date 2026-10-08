"""Genuine firstbridge cutoff branches and exact paired Poisson Z supports.

Only complete derivative covers are tightened. Original C0 A/B objects,
the native density owner, true widths and all thirteen upstream cells stay.
"""
from contextlib import contextmanager
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_upstream_tile_uniform_C1 as baseline
import lei_ren_part1_paper_compliant_current_native_paired_C1_transport as paired

common=baseline.common;middle=baseline.middle;cutoff=paired.cutoff
q2=paired.accepted.accepted;prior=q2.prior;ep=baseline.ep
HERE,PREFIX,sha=baseline.HERE,baseline.PREFIX,baseline.sha
NAME=PREFIX+'current_native_upstream_paired_branch_C1.json'
RECEIPT=PREFIX+'current_native_upstream_paired_branch_C1_check.json'
GATE='genuine_firstbridge_conditional_cutoff_paired_C1_source_histories_executed'
REPLAY=PREFIX+'current_original_O2_paired_branch_upstream_replay.json'
N=baseline.N;KEYS=baseline.KEYS;LABELS=baseline.LABELS
ZERO=(0,0);DZ=(0,1);ORDERS=cutoff.ORDERS
encode=lambda v:middle.packets.encode(common.base.encoded(v))


def branch_paired_cover(roots,eta_log,log_a_lower,dstar_log):
    """Evaluate each nonempty original cutoff branch, then hull its Z covers."""
    a=roots['a'][ZERO];c=a.ctx
    if ep(eta_log)[1]>ep(-c.ln(2))[0]:raise ValueError('Original eta<=1/2 support proof required')
    candidates,empty=cutoff.conditional_cutoff_branches(roots['kappa_minus2'][ZERO],eta_log)
    branches=[];values=[]
    for branch in candidates:
        restricted=dict(roots,kappa_minus2=dict(roots['kappa_minus2']))
        restricted['kappa_minus2'][ZERO]=branch['Delta']
        q=cutoff.conditional_q_jet(restricted,eta_log,log_a_lower,branch)
        square=q2.conditional_q2_jet(restricted,eta_log,log_a_lower,branch)
        got=paired.paired_C1_support(restricted,q['rows'],square,log_a_lower,dstar_log)
        values.append(got['values'])
        branches.append(dict(name=branch['name'],condition=branch['condition'],
            original_conditional_Delta_C0=branch['Delta'].record(),theta=branch.get('theta'),
            original_q_jet={('y%d_Z%d'%order):value.record() for order,value in q['rows'].items()},
            original_direct_q_squared_jet=q2.q2_record(square,branch),
            original_q_branch_derivation={k:v for k,v in q.items() if k!='rows'},
            original_paired_whole_period_support=got['record'],
            paired_derivative_covers={k:v.record() for k,v in got['values'].items()},
            only_Delta_C0_restricted_all_original_derivative_objects_retained=True))
    union=middle.density.local.same_source_union
    hull={k:union([row[k] for row in values]) for k in ('A_Z','B_Z_over_Pstar')}
    return dict(values=hull,record=dict(
        native_source_log_bases=a.scale.bases,
        original_root_rows={name:{('y%d_Z%d'%order):value.record() for order,value in rows.items()} for name,rows in roots.items()},
        original_eta_log=eta_log,original_log_a_positive_lower=log_a_lower,original_dstar_log=dstar_log,
        complete_cutoff_branch_cover='Delta<=0 union 0<=Delta<=eta union Delta>=eta',
        conditional_cutoff_paired_branches=branches,branches_proved_empty=empty,
        whole_source_paired_derivative_hulls={k:v.record() for k,v in hull.items()},
        paired_parameter_theorem=paired.paired_parameter_theorem(),
        active_a_b_support_constraints_applied_only_inside_conditional_branches=True,
        original_linear_q_Z_and_nonzero_p2_Z_retained=True,
        conditional_branches_hulled_not_added_or_selected=True,
        no_derivative_of_C0_caps_or_branch_selectors=True))


@contextmanager
def firstbridge_paired_primitives():
    traces=[]
    with baseline.uniform_upstream_primitives() as inherited:
        original=middle.p.whole_period_C1
        def replacement(roots,eta_log,log_a_lower,dstar_log):
            old=original(roots,eta_log,log_a_lower,dstar_log)
            index=len(traces)
            if index>=len(LABELS):raise ValueError('Expected exactly twelve original upstream primitive calls')
            if index!=0:
                traces.append(old['record']);return old
            got=branch_paired_cover(roots,eta_log,log_a_lower,dstar_log)
            selected=dict(old['values']);comparisons={}
            for key in ('A_Z','B_Z_over_Pstar'):
                before=baseline.magnitude_log(old['values'][key]);after=baseline.magnitude_log(got['values'][key])
                tighter=after is None or (before is not None and after<before)
                if tighter:selected[key]=got['values'][key]
                comparisons[key]=dict(accepted_uniform_absolute_upper_log=before,paired_branch_union_absolute_upper_log=after,
                    selected_backend='paired_conditional_branch_union' if tighter else 'accepted_uniform_cover',
                    strict_absolute_upper_reduction=tighter and before!=after)
            if selected['A'] is not old['values']['A'] or selected['B_over_Pstar'] is not old['values']['B_over_Pstar']:
                raise ValueError('Original C0 primitive enclosure objects must remain unchanged')
            trace=dict(label=LABELS[index],accepted_uniform_baseline_cover=old['record'],
                genuine_same_source_conditional_paired_cover=got['record'],paired_derivative_comparisons=comparisons,
                selected_primitive_C0_Z_covers={k:v.record() for k,v in selected.items()},
                original_C0_A_B_objects_and_native_density_owner_unchanged=True,
                full_cutoff_branch_union_encloses_original_source_derivatives=True)
            traces.append(trace)
            return dict(values=selected,record=trace)
        middle.p.whole_period_C1=replacement
        try:yield traces
        finally:
            changed=middle.p.whole_period_C1 is not replacement
            middle.p.whole_period_C1=original
            if changed:raise RuntimeError('Firstbridge paired primitive runtime binding changed unexpectedly')
        if len(inherited)!=len(traces):raise ValueError('Every original cell must retain its baseline source cover')


@middle.native.inlet.source_precision
def run():
    begin=time.monotonic();bridge,seed=middle.native.inlet.native_bridge_owner()
    with middle.native.inlet.CheckedSourceRuntime():
        owner=baseline.build_owner(bridge);hashes,binding=baseline.accepted_hashes(owner.family)
        for module in (baseline,paired,q2,cutoff):common.attach_receipt(hashes,module,owner.family)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Same original source closures disagree: '+name)
            hashes[name]=digest
        archives=[];comparisons={}
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            with firstbridge_paired_primitives() as traces:route=owner.route(Z,N=N)
            if len(traces)!=12:raise ValueError('Every original active upstream cell required')
            reductions=sum(r['strict_absolute_upper_reduction'] for r in traces[0]['paired_derivative_comparisons'].values())
            cells=baseline.live_active_cells(route)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                complete_actual_pre_O2_tile_route=route['record'],exact_upstream_source_function_incoming_binding=binding,
                twelve_source_primitive_cover_records=traces,firstbridge_paired_derivative_strict_reductions=reductions,
                per_chart_actual_true_width_contribution_attribution={label:dict(C0={k:v.record() for k,v in cell['contributions'].items()},
                    Z={k:v.record() for k,v in cell['Z_derivatives'].items()}) for label,cell in cells},
                full_source_bounds_requeried_on_this_tile_not_relabelled_whole_Z=True,
                original_source_normalization_densities_masses_and_incoming_equations_unchanged=True,
                only_firstbridge_ordinary_Z_primitive_cover_tightened=True)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_native_upstream_paired_branch_C1_'+tag+'.json.gz'
            if len(compressed)>=100*1024*1024:raise ArithmeticError('Split genuine native source evidence losslessly')
            (HERE/filename).write_bytes(compressed)
            archive=dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),firstbridge_paired_derivative_strict_reductions=reductions)
            archives.append(archive);comparisons[tag]=traces[0]['paired_derivative_comparisons']
            print('Genuine firstbridge paired branch histories:',tag,'strict primitive derivative reductions',reductions,flush=True)
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            genuine_firstbridge_paired_branch_upstream_archives=archives,firstbridge_derivative_bound_comparisons=comparisons,
            exact_upstream_source_function_incoming_binding=binding,source_seed_evidence=seed,
            same_original_C0_functions_all13_cells_actual_pressure_memory_and_native_owner_retained=True,
            full_Rc_terminal_closure_or_global_N_or_stress_or_recursion_claimed=False,
            **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-begin,
            scope='Genuine same-N1024 strict-sign original13-cell upstream source histories; only firstbridge ordinary-Z primitive covers refined with complete conditional cutoff union/direct q²/paired Poisson parameter identities. Original C0 A/B, source densities/units/widths/owner unchanged. No full axial/Rc matching/controls/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
