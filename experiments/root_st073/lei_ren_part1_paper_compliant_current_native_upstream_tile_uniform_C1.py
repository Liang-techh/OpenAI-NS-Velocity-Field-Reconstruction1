"""Genuine pre-O2 axial tiles with accepted parameter-uniform periodic jets.

Only the primitive enclosure algorithm is selected in a serial context.
Original source queries, cutoff branches, densities, masses and all actual
incoming histories are unchanged. No integral from another N is inserted.
"""
from contextlib import contextmanager
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories as middle
import lei_ren_part1_paper_compliant_current_native_serial_C1_cell as serial
import lei_ren_part1_paper_compliant_current_native_O2_C1_histories as theorem
import lei_ren_part1_paper_compliant_current_original_O2_source_incoming_common_N as common

HERE,PREFIX,sha=middle.HERE,middle.PREFIX,middle.sha
NAME=PREFIX+'current_native_upstream_tile_uniform_C1.json'
RECEIPT=PREFIX+'current_native_upstream_tile_uniform_C1_check.json'
GATE='actual_pre_O2_strict_sign_tiles_with_uniform_periodic_C1_enclosures_executed'
N=1024;KEYS=tuple(middle.RATES);ep=middle.ep
LABELS=('active_first_bridge','second_bridge','bridge_macro','switch_first','switch_second',
    'switch_power','reshape','inner_reference','axial_restore','restore_buffer','actual_patch','Rh_reference')


def magnitude_log(value):
    if value.zero:return None
    return ep(serial.absolute_upper(value).scale.evaluate())[1]


@contextmanager
def uniform_upstream_primitives():
    """Choose complete certified function covers, never source point values."""
    original=middle.p.whole_period_C1;calls=[]
    def replacement(roots,eta_log,log_a_lower,dstar_log):
        old=original(roots,eta_log,log_a_lower,dstar_log)
        new=serial.whole_period_C1(roots,eta_log,log_a_lower,dstar_log)
        selected={};decisions={}
        for key in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar'):
            a,b=old['values'][key],new['values'][key]
            a.coerce(b)
            la,lb=magnitude_log(a),magnitude_log(b)
            tighter=lb is None or (la is not None and lb<la)
            selected[key]=b if tighter else a
            decisions[key]=dict(old_absolute_upper_log=la,uniform_absolute_upper_log=lb,
                selected_backend='parameter_uniform' if tighter else 'inherited_whole_period',
                strict_absolute_upper_reduction=tighter and la!=lb,
                complete_certified_source_function_cover_selected_not_field_value=True)
        trace=dict(original_periodic_parameter_C1_theorem=serial.periodic_parameter_theorem(),
            same_original_source_root_objects=True,
            old_whole_period_primitive_cover=old['record'],
            parameter_uniform_primitive_cover=new['record'],
            selected_primitive_C0_Z_covers={k:v.record() for k,v in selected.items()},
            same_source_bound_comparisons=decisions,
            all_original_cutoff_branches_and_q_Z_retained=True,
            no_derivative_of_interval_selector_or_C0_cap=True)
        calls.append(trace)
        return dict(values=selected,record=trace)
    middle.p.whole_period_C1=replacement
    try:yield calls
    finally:
        changed=middle.p.whole_period_C1 is not replacement
        middle.p.whole_period_C1=original
        if changed:raise RuntimeError('Serial upstream primitive binding changed unexpectedly')


def build_owner(bridge):
    p=middle.p
    c1=middle.density.NativeDensityC1LocalIntegrals(p.first.NativePhaseFirstJets(p.slow.NativeQSlowJets(
        p.current.NativeCorrelatedShearQ(middle.prior.NativeSignedInputEnclosures(middle.native.NativeGenericSourcePackets(bridge))))))
    second=middle.preceding.preceding.NativeSecondBridgeC1Histories(p.NativeFirstBridgeC1Histories(
        middle.transfer.NativeTrueChartC1Transfer(middle.history.NativeC1HistoryTransfer(c1))))
    return middle.NativeMiddleO2InletC1Histories(middle.preceding.NativeBridgeSwitchC1Histories(second))


def live_active_cells(route):
    switch=route['initial'];second=switch['initial'];first=second['incoming']
    result=[('active_first_bridge',first),('second_bridge',second)]
    result.extend(switch['cells'].items());result.extend(route['cells'].items())
    if tuple(label for label,cell in result)!=LABELS:raise ValueError('All twelve original active pre-O2 cells required')
    return result


def accepted_hashes(family):
    hashes={}
    for module in (middle,theorem,common):
        report,receipt=common.attach_receipt(hashes,module,family)
    binding=report['exact_upstream_source_function_incoming_binding']
    if binding['ordered_pre_O2_function_cells']!=['initial_flat_collar',*LABELS]:
        raise ValueError('Exact original thirteen-cell function prefix required')
    if receipt['input_hashes'].get(Path(serial.__file__).name)!=sha(Path(serial.__file__).name):
        raise ValueError('Accepted original parameter-uniform theorem source must be bound')
    return hashes,binding


def final_summary(route,baseline):
    rows={}
    c=route['correction']['values']['m'].ctx
    with mp.workdps(c.dps+40):
        for kind,field in (('values','actual_original_inlet_to_O2_inlet_correction_C0'),
                ('Z_derivatives','actual_original_inlet_to_O2_inlet_correction_Z')):
            for key,value in route['correction'][kind].items():
                old=common.interval(c,baseline[field][key]['log_absolute_upper']) if baseline[field][key]['log_absolute_upper'] is not None else None
                before=None if old is None else ep(old)[1];after=magnitude_log(value)
                rows[kind+'_'+key]=dict(accepted_whole_Z_absolute_upper_log=before,
                    actual_tile_uniform_absolute_upper_log=after,
                    strict_upper_reduction_vs_whole_Z=after is None and before is not None or
                        after is not None and before is not None and after<before,
                    comparison_is_bound_not_observed_physical_amplitude=True)
    return rows


@middle.native.inlet.source_precision
def run():
    begin=time.monotonic();bridge,seed_evidence=middle.native.inlet.native_bridge_owner()
    with middle.native.inlet.CheckedSourceRuntime():
        owner=build_owner(bridge);hashes,binding=accepted_hashes(owner.family)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Original source closures disagree: '+name)
            hashes[name]=digest
        baseline=json.loads((HERE/middle.NAME).read_bytes())['actual_inlet_to_O2_inlet_C1_history_records']['whole_Z']
        archives=[];summaries={}
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            with uniform_upstream_primitives() as calls:
                route=owner.route(Z,N=N)
            if len(calls)!=12:raise ValueError('All twelve active cells must use reviewed uniform primitive selection')
            reductions=sum(d['strict_absolute_upper_reduction'] for r in calls for d in r['same_source_bound_comparisons'].values())
            if not reductions:raise ArithmeticError('No actual same-source primitive bound was tightened')
            cells=live_active_cells(route)
            attribution={label:dict(
                C0={k:v.record() for k,v in cell['contributions'].items()},
                Z={k:v.record() for k,v in cell['Z_derivatives'].items()}) for label,cell in cells}
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                complete_actual_pre_O2_tile_route=route['record'],
                exact_upstream_source_function_incoming_binding=binding,
                twelve_same_source_uniform_primitive_comparisons=calls,
                per_chart_actual_true_width_contribution_attribution=attribution,
                same_source_primitive_strict_upper_reductions=reductions,
                full_source_bounds_requeried_on_this_tile_not_relabelled_whole_Z=True,
                original_source_normalization_densities_masses_and_incoming_equations_unchanged=True)
            comparison=final_summary(route,baseline);summaries[tag]=comparison
            raw=json.dumps(middle.packets.encode(common.base.encoded(payload)),indent=2).encode()+b'\n'
            compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            if len(compressed)>=100*1024*1024:raise ArithmeticError('Split complete native source evidence losslessly')
            filename=PREFIX+'current_native_upstream_tile_uniform_C1_'+tag+'.json.gz'
            (HERE/filename).write_bytes(compressed)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),
                uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
                same_source_primitive_strict_upper_reductions=reductions))
            print('Genuine upstream tile with uniform periodic jets:',tag,'same-source reductions',reductions,flush=True)
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        hashes[Path(serial.__file__).name]=sha(Path(serial.__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            actual_strict_sign_upstream_tile_archives=archives,
            exact_upstream_source_function_incoming_binding=binding,
            final_actual_incoming_bound_comparisons=summaries,
            source_seed_evidence=seed_evidence,parameter_uniform_theorem=serial.periodic_parameter_theorem(),
            all_twelve_active_charts_and_true_initial_flat_collar_retained=True,
            original_P0_and_actual_correction_memory_separate=True,
            no_original_O2_integral_recomputed_or_N7_contribution_reused=True,
            global_terminal_closure_or_global_N_or_stress_or_recursion_claimed=False,
            **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,
            execution_seconds=time.monotonic()-begin,
            scope='Actual original sc/2-to-O2-slope0 correction C0/Z histories on continuous strict-sign Z tiles at N1024, all13 original cells, accepted parameter-uniform whole-period primitive bounds selected against inherited bounds on the same source. No original field/source/cutoff altered, no O2 integral rerun, no global closure/N/stress/recursion/full NS claim.')
        (HERE/NAME).write_text(json.dumps(middle.packets.encode(common.base.encoded(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
