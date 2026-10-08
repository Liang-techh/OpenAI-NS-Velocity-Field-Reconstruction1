"""Continue genuine refined O2-slope exits through the original six Rc cells.

The native source is queried on the actual strict-sign tiles. These are full
source-function covers, not selected numerical values or terminal closure.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_inner_reference_weighted_pressure as inlet
import lei_ren_part1_paper_compliant_current_native_Rc_C1_histories as rc

common=inlet.common;baseline=inlet.upstream.baseline;serial=rc.serial
HERE,PREFIX,sha=inlet.HERE,inlet.PREFIX,inlet.sha;N=inlet.N;KEYS=inlet.KEYS
ep=inlet.ep;iv=inlet.iv;encode=inlet.encode
NAME=PREFIX+'current_weighted_O2_to_Rc_tail.json'
RECEIPT=PREFIX+'current_weighted_O2_to_Rc_tail_check.json'
GATE='actual_weighted_original_O2_exits_continued_through_all_six_tail_cells_to_Rc'
LABELS=('axial','buffer_0_9','buffer_9_11','transition','power_to_r_plus','power_to_Rc')


def inlet_functions(record,coordinates):
    coordinates.require_family(record['source_family'])
    if record['original_P0_datum_sha256']!=coordinates.family['datum_enclosure_sha256']:
        raise ValueError('Same accepted original axis pressure datum required')
    if record['normalized_own_units']!=common.previous.five.UNITS:
        raise ValueError('Same original five normalized units required')
    if record['candidate_N']!=N or record['exact_y_window']!=['0','1']:
        raise ValueError('Same genuine N1024 entire slope-chart exit required')
    return dict(values={k:common.restore_common_source(v,coordinates)
        for k,v in record['propagated_actual_correction_C0'].items()},
        Z_derivatives={k:common.restore_common_source(v,coordinates)
        for k,v in record['propagated_actual_correction_Z'].items()})


def source_route(owner,Z,incoming):
    """No upstream owner.route call: the actual accepted slope exit is injected."""
    o2=owner.owner;c=owner.ctx;cells=[];previous_chart='O2_slope';previous_endpoint=1
    for label,chart,left,right in rc.preceding.ROUTE[5:]:
        geometry=owner.transfer.geometry.cell(chart,left,right)
        if previous_chart==chart:
            seam='same '+chart+' chart at '+str(left)
            join=dict(same_original_chart_and_source_function=True,coordinate=str(left))
        else:
            seam=previous_chart+' -> '+chart
            if seam not in owner.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
                raise ValueError('Original O2 tail radius/phase seam required')
            join=o2.joins['slope_axial' if chart=='O2_axial' else 'axial_buffer']
        cell=serial.serial_cell(o2,chart=chart,geometry=geometry,endpoint=c.mpf(right),Z=Z,N=N,
            incoming=incoming,left_record=dict(chart=previous_chart,coordinate=previous_endpoint),
            right_record=dict(chart=chart,coordinate=right),seam=seam,join=join)
        cell['record']['original_tail_label']=label;cells.append(cell)
        incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
        print('Original refined tail source:',Z,label,flush=True)
    Tw=owner.transfer.geometry.binder.fixed['Tw']
    for label,chart,left,right in rc.ROUTE:
        if chart=='O3_power':
            geometry=owner.transfer.geometry.cell(chart,dict(original_power_offset=left),dict(original_power_offset=right))
            endpoint=c.mpf(right)/Tw
            left_record=dict(chart=previous_chart,coordinate=previous_endpoint) if left==0 else dict(chart=chart,original_power_offset=left)
            right_record=dict(chart=chart,original_power_offset=right,original_phase_expression=str(right)+'/Tw')
        else:
            geometry=owner.transfer.geometry.cell(chart,left,right);endpoint=c.mpf(right)
            left_record=dict(chart=previous_chart,coordinate=previous_endpoint)
            right_record=dict(chart=chart,coordinate=right)
        if previous_chart==chart:
            seam='same original '+chart+' source';join=dict(same_original_chart_and_source_function=True,original_power_offset=left)
        else:
            seam=previous_chart+' -> '+chart
            if seam not in owner.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
                raise ValueError('Original O3 tail radius/phase seam required')
            join=o2.joins['buffer_transition'] if chart=='O3_slope_mu' else owner.join
        cell=rc.original_O3_cell(owner,chart=chart,geometry=geometry,endpoint=endpoint,Z=Z,N=N,
            incoming=incoming,left_record=left_record,right_record=right_record,seam=seam,join=join)
        if chart=='O3_power':
            if not all(v.zero for v in [*cell['contributions'].values(),*cell['Z_derivatives'].values()]):
                raise ValueError('Original canonical power flat source must have zero own increments')
            cell['record'].update(original_whole_power_q_q_Z_and_density_C0_Z_exact_zero=True,
                original_actual_power_offsets=[left,right],quiet_local_density_does_not_reset_inherited_histories=True)
        cell['record']['original_tail_label']=label;cells.append(cell)
        incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
        print('Original refined tail source:',Z,label,flush=True)
    if tuple(cell['record']['original_tail_label'] for cell in cells)!=LABELS:
        raise ValueError('All original six tail cells required')
    return cells


def replay_saved(coordinates,input_record,rows):
    incoming=inlet_functions(input_record,coordinates);path=[]
    for row in rows:
        values={};jets={}
        for k in KEYS:
            factor=common.restore_common_source(row['original_true_width_kernel_factors'][k]['decay'],coordinates)
            addition=common.restore_common_source(row['actual_true_width_C0_contributions'][k],coordinates)
            additionZ=common.restore_common_source(row['actual_true_width_Z_contributions'][k],coordinates)
            values[k]=factor*incoming['values'][k]+addition
            jets[k]=factor*incoming['Z_derivatives'][k]+additionZ
        background=row['original_right_background_and_separate_P0_Z']
        originals={k:common.restore_common_source(v,coordinates) for k,v in background['original_normalized_history_C0_enclosures'].items()}
        originalZ={k:common.restore_common_source(v,coordinates) for k,v in background['original_normalized_history_Z_enclosures'].items()}
        path.append(dict(label=row['original_tail_label'],actual_inherited_C0={k:v.record() for k,v in incoming['values'].items()},
            actual_inherited_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
            actual_right_correction_C0={k:v.record() for k,v in values.items()},actual_right_correction_Z={k:v.record() for k,v in jets.items()},
            actual_right_own_history_C0={k:(originals[k]+values[k]).record() for k in KEYS},
            actual_right_own_history_Z={k:(originalZ[k]+jets[k]).record() for k in KEYS}))
        incoming=dict(values=values,Z_derivatives=jets)
    P0=common.restore_common_source(background['original_separate_P0_over_Pstar_squared'],coordinates)
    P0Z=common.restore_common_source(background['original_separate_P0_Z_over_Pstar_squared'],coordinates)
    own={k:originals[k]+values[k] for k in KEYS};ownZ={k:originalZ[k]+jets[k] for k in KEYS}
    return dict(source_family=coordinates.family,candidate_N=N,exact_Z_range=input_record['exact_Z_range'],
        original_O2_ordered_source_cells=input_record['ordered_source_cells'],
        genuine_slope_exit_binding=dict(manifest=inlet.NAME,manifest_sha256=sha(inlet.NAME),receipt=inlet.RECEIPT,
            receipt_sha256=sha(inlet.RECEIPT),exact_Z_range=input_record['exact_Z_range'],
            ordered_source_cells=input_record['ordered_source_cells'],whole_original_slope_y_window=['0','1'],
            input_is_actual_correction_not_original_plus_correction=True),
        ordered_original_tail_labels=list(LABELS),six_original_cell_actual_histories=path,
        actual_Rc_correction_C0={k:v.record() for k,v in values.items()},actual_Rc_correction_Z={k:v.record() for k,v in jets.items()},
        actual_Rc_own_history_C0={k:v.record() for k,v in own.items()},actual_Rc_own_history_Z={k:v.record() for k,v in ownZ.items()},
        original_Rc_background_and_separate_P0_Z=background,
        original_Rc_endpoint=dict(radius_expression='Rw*exp(2)',original_power_offset=2,original_power_phase_expression='2/Tw'),
        actual_absolute_pressure_over_Pstar_squared=(P0+own['p']).record(),
        actual_absolute_pressure_Z_over_Pstar_squared=(P0Z+ownZ['p']).record(),
        final_correction_absolute_upper_logs={k:common.absolute_log_upper(v) for k,v in values.items()},
        final_correction_Z_absolute_upper_logs={k:common.absolute_log_upper(v) for k,v in jets.items()},
        original_source_covers_not_selected_field_values_or_terminal_defects=True,
        all_six_cells_keep_actual_incoming_and_rate_zero_pressure_memory=True,
        numerical_original_source_oracle_installed=False,actual_original_numerical_tail_integrals_evaluated=False,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))


@rc.native.inlet.source_precision
def run():
    begin=time.monotonic();accepted=json.loads((HERE/inlet.NAME).read_bytes())
    hashes={};common.attach_receipt(hashes,inlet,accepted['source_family'])
    common.attach_receipt(hashes,rc,accepted['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(baseline.build_owner(bridge)))
        if owner.family!=accepted['source_family']:raise ValueError('Same accepted original source family required')
        archives=[];replays=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            inputs=[r for r in accepted['genuine_original_O2_weighted_pressure_replays'] if inlet.same_Z(r['exact_Z_range'],Z)]
            if len(inputs)!=2:raise ValueError('Both accepted genuine O2 refinements required')
            primary=next(r for r in inputs if r['ordered_source_cells']==2048)
            incoming=inlet_functions(primary,owner.coordinates)
            cells=source_route(owner,Z,incoming)
            rows=encode([cell['record'] for cell in cells])
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                genuine_primary_incoming_O2_record=primary,ordered_original_six_cell_source_records=rows,
                exact_native_original_radius_Jacobian_and_seam_binding=owner.transfer.geometry.binder.identity,
                original_O3_transition_power_join=owner.join,original_Rc_reservation=owner.reservation,
                source_requeried_on_actual_axial_tile_not_relabelled_whole_Z=True,
                actual_original_source_density_geometries_and_backgrounds_unchanged=True,
                only_original_six_tail_cells_executed_upstream_route_not_recomputed=True)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_weighted_O2_to_Rc_tail_'+tag+'.json.gz'
            if len(compressed)>=100*1024*1024:raise ValueError('Split large source evidence losslessly')
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archive=dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest());archives.append(archive)
            for input_record in inputs:
                row=replay_saved(owner.coordinates,input_record,rows)
                row['actual_six_cell_source_binding']=dict(archive=filename,sha256=sha(filename),
                    source_is_genuine_original_same_N_tile_function_cover=True)
                replays.append(row)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Original source closures disagree: '+name)
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            original_six_cell_tail_source_archives=archives,actual_refined_O2_to_Rc_replays=replays,
            input_hashes=hashes,source_seed_evidence=seed,execution_seconds=time.monotonic()-begin,
            previous_weighted_pressure_and_O2_integrals_reused_without_recomputation=True,
            full_original_pre_O2_prefix_plus_entire_O2_slope_plus_all_six_tail_cells_connected=True,
            numerical_original_source_oracle_installed=False,actual_original_numerical_tail_integrals_evaluated=False,
            full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False),
            scope='Genuine original same-N1024 weighted pressure/active-kappa O2 exits, both strict-sign tiles, all6 remaining cells to actual Rc offset2/Tw. Full source-domain covers and actual inherited five C0/Z rows, not numerical tail integral values/target closure/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
