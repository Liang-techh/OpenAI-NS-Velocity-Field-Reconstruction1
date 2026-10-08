"""Batch centered C1 supports on the eleven remaining original upstream charts."""
import copy,gzip,hashlib,json,time
from pathlib import Path
import lei_ren_part1_paper_compliant_current_firstbridge_centered_C1 as previous

common,weighted,terminal,slope,centered=previous.common,previous.weighted,previous.terminal,previous.slope,previous.centered
HERE,PREFIX,sha,encode,iv,ep=previous.HERE,previous.PREFIX,previous.sha,previous.encode,previous.iv,previous.ep
N,KEYS=previous.N,previous.KEYS
NAME=PREFIX+'current_upstream_centered_C1.json';RECEIPT=PREFIX+'current_upstream_centered_C1_check.json'
GATE='actual_all12_upstream_centered_C1_transport_and_Rc_targets_executed'
LABELS=weighted.upstream.LABELS


def live_chart(owner,payload,label,index):
    c=owner.ctx;row=weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])[1][index][1]
    sourcekey=weighted.accepted.row_fields(label)[1];prov=row[sourcekey]['source_provenance']
    query=owner.q_owner.query(prov['chart'],iv(c,prov['Z_box']),iv(c,prov['coordinate_box']))
    assert encode(query['record']['source_provenance'])==prov
    roots=query['source']['roots'];E=roots['E'][(0,0)];packet=query['source']['packet'];signed=owner.transfer.owner.signed_owner
    root_owner=owner.q_owner.owner.owner
    positive=root_owner.decode(root_owner.inventory[prov['chart']]['actual_positive_denominator_theorem'])
    velocities=[]
    for k in (0,1):
        leaf=slope.prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:slope.prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)})
        velocities.append(signed.leaf(leaf,E.scale.bases,E.ledger).record())
    return encode(dict(label=label,original_source_provenance=prov,native_source_bases=E.scale.bases,
        original_root_rows={name:{'y%d_Z%d'%order:v.record() for order,v in rows.items()} for name,rows in roots.items()},
        original_positive_a_certificate=positive,original_eta_log=root_owner.scales['selected_positive_eta_log'],
        original_dstar_log=root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'],V=velocities[0],V_Z=velocities[1],
        original_packet_queried_with_exact_saved_chart_Z_coordinate=True))


def chart_cover(c,payload,live,index,sigma):
    label=live['label'];bases=tuple(iv(c,v) for v in live['native_source_bases'])
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    restore=lambda v:slope.native_restore(v,bases,ledger)
    roots={name:{(int(k[1]),int(k[-1])):restore(v) for k,v in rows.items()} for name,rows in live['original_root_rows'].items()}
    trace=payload['twelve_source_primitive_cover_records'][index];old={k:restore(v) for k,v in trace['selected_primitive_C0_Z_covers'].items()}
    caps=centered.primitive_caps(c,roots,live['original_positive_a_certificate'],iv(c,live['original_eta_log']),iv(c,live['original_dstar_log']),sigma,dict(values=old))
    selected=dict(old);decisions={}
    for key,kind in (('A_Z','A'),('B_Z_over_Pstar','B')):
        cap=caps[kind][centered.DZ]
        candidate=slope.prior.ScaledEnclosure(slope.prior.FormalScale(bases,offset=cap.log),c.mpf(('-1','1')),ledger) if cap.log is not None else old[key].scalar(0)
        before=weighted.upstream.baseline.magnitude_log(old[key]);after=weighted.upstream.baseline.magnitude_log(candidate)
        tighter=after is None or before is not None and after<before;selected[key]=candidate if tighter else old[key]
        decisions[key]=dict(old_log_absolute_upper=before,new_log_absolute_upper=after,strict_upper_reduction=tighter and before!=after)
    assert selected['A'] is old['A'] and selected['B_over_Pstar'] is old['B_over_Pstar']
    den=weighted.upstream.middle.density.density_Z_kernels(roots['E'][(0,0)],roots['E'][(0,1)],restore(live['V']),restore(live['V_Z']),selected,N)
    row=weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])[1][index][1];geo=row[weighted.accepted.row_fields(label)[0]];coords=weighted.replay.inputs(c,payload)[0]
    g=dict(width=common.restore_common_source(geo['positive_true_log_radius_width'],coords),regular=iv(c,geo['regular_true_log_radius_width_cover']),scalar_cover=iv(c,geo['scalar_width_cover_used_only_for_directed_kernel_bounds']))
    values={};comparisons={}
    for key,rate in slope.rc.RATES.items():
        mass=slope.rc.transfer.true_width_kernel(coords,g,rate)['mass'];candidate=coords.rebase(den['Z_derivatives'][key],coords.family)*mass
        accepted=common.restore_common_source(row['true_log_radius_signed_Z_contributions'][key],coords)
        before=weighted.upstream.baseline.magnitude_log(accepted);after=weighted.upstream.baseline.magnitude_log(candidate)
        tighter=after is None or before is not None and after<before;values[key]=(candidate if tighter else accepted).record()
        comparisons[key]=dict(old_log_absolute_upper=before,new_log_absolute_upper=after,strict_upper_reduction=tighter and before!=after)
    return dict(label=label,selected_Z_contributions=encode(values),record=dict(original_live_source=live,
        centered_primitive_proof=caps['record'],primitive_comparisons=decisions,contribution_comparisons=comparisons,
        selected_original_primitive_C0_Z={k:v.record() for k,v in selected.items()},all_original_C0_and_native_width_conversion_unchanged=True))


def transport(c,payload,charts,pressure):
    changed=copy.deepcopy(payload);flat,rows=weighted.accepted.route_rows(changed['complete_actual_pre_O2_tile_route']);coords=weighted.replay.inputs(c,payload)[0]
    jets={k:common.restore_common_source(flat['actual_correction_Z_enclosures'][k],coords) for k in KEYS};path=[]
    assert [label for label,row in rows]==list(LABELS)
    for label,row in rows:
        row['true_log_radius_signed_Z_contributions']=charts[label]['selected_Z_contributions']
        geo=row[weighted.accepted.row_fields(label)[0]]
        g=dict(width=common.restore_common_source(geo['positive_true_log_radius_width'],coords),regular=iv(c,geo['regular_true_log_radius_width_cover']),scalar_cover=iv(c,geo['scalar_width_cover_used_only_for_directed_kernel_bounds']))
        jets={k:slope.rc.transfer.true_width_kernel(coords,g,rate)['decay']*jets[k]+common.restore_common_source(row['true_log_radius_signed_Z_contributions'][k],coords) for k,rate in slope.rc.RATES.items()}
        path.append(dict(label=label,actual_right_correction_Z={k:v.record() for k,v in jets.items()}))
    half=slope.native.HalfPstarCoordinates(c,coords.logP_squared,payload['source_family'])
    return dict(jets=encode({k:v.record() for k,v in jets.items()}),trace=slope.pre_slope_attribution(half,changed,pressure),
        record=dict(ordered_actual_Z_transport=path,all_original_C0_P0_true_widths_and_zero_rate_pressure_retained=True,
        old_saved_output_fields_not_relabelled_as_new=True))


@slope.rc.native.inlet.source_precision
def run():
    began=time.monotonic();m=json.loads((HERE/previous.NAME).read_bytes());hashes={};common.attach_receipt(hashes,previous,m['source_family'])
    wm=json.loads((HERE/weighted.NAME).read_bytes());cm=json.loads((HERE/slope.source.NAME).read_bytes());sig=json.loads((HERE/(PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
    bridge,seed=slope.rc.native.inlet.native_bridge_owner();archives=[]
    with slope.rc.native.inlet.CheckedSourceRuntime():
        owner=weighted.upstream.baseline.build_owner(bridge);c=owner.ctx;sigma={int(k):centered.previous.read_cap(c,v) for k,v in sig['sigma_global_derivative_log_caps'].items()}
        for old_archive in m['actual_firstbridge_centered_C1_archives']:
            old=terminal.load_archive(old_archive);payload=terminal.load_archive(old['accepted_original_pre_slope_archive']);sd=terminal.load_archive(old['accepted_original_slope_archive']);Z=old['exact_Z_range']
            charts={'active_first_bridge':dict(label='active_first_bridge',selected_Z_contributions=old['selected_firstbridge_Z_contributions'],record=dict(accepted_firstbridge_centered_archive=old_archive))}
            for index,label in enumerate(LABELS[1:],start=1):
                live=live_chart(owner,payload,label,index);charts[label]=chart_cover(c,payload,live,index,sigma)
                print('Actual centered upstream source:',Z,label,'executed',flush=True)
            pressure=next(r for r in wm['actual_weighted_pressure_tiles'] if weighted.same_Z(r['exact_Z_range'],Z));upstream=transport(c,payload,charts,pressure)
            primary=next(r for r in wm['genuine_original_O2_weighted_pressure_replays'] if r['ordered_source_cells']==2048 and weighted.same_Z(r['exact_Z_range'],Z))
            report=next(r for r in cm['original_complete_same_N_source_incoming_O2_integral_refinements'] if r['ordered_source_cells']==2048 and weighted.same_Z(r['exact_Z_range'],Z));output=previous.updated_terminal(c,sd,upstream,primary,report)
            result=dict(source_family=owner.family,candidate_N=N,exact_Z_range=Z,actual_centered_upstream_chart_sources=charts,
                actual_updated_upstream_Z=upstream['jets'],actual_updated_upstream_attribution=upstream['trace'],actual_twelve_chart_transport=upstream['record'],actual_updated_Rc_terminal_defects=output,
                accepted_firstbridge_centered_archive=old_archive,actual_global_function_control_N_heat_stress_recursion_admitted=False,**dict.fromkeys(common.current.FLAGS,False))
            raw=json.dumps(encode(result),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0);tag='positive' if slope.Fraction(Z[0])>0 else 'negative';name=PREFIX+'current_upstream_centered_C1_'+tag+'.json.gz';(HERE/name).write_bytes(compressed);hashes[name]=sha(name)
            archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),exact_Z_range=Z))
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Original source closure disagreement: '+name)
            hashes[name]=digest
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=m['source_family'],candidate_N=N,actual_all_upstream_centered_C1_archives=archives,
        only22_remaining_upstream_native_packets_queried=True,accepted_firstbridge_slope_and_downstream_producers_and_inverse_solvers_not_rerun=True,
        actual_global_function_controls_or_recursion_admitted=False,**dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
