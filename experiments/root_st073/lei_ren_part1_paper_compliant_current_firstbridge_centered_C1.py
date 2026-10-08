"""Apply accepted centered phase identities to the actual firstbridge C1.

Only two original firstbridge packets are queried for their axial velocity.
All other upstream charts, slope integrals and downstream inverse data are
reused. Original C0/P0 stay unchanged; only tighter derivative covers enter.
"""
import copy,gzip,hashlib,json,time
from pathlib import Path
import lei_ren_part1_paper_compliant_current_slope_joint_source as slope
import lei_ren_part1_paper_compliant_current_native_centered_phase_conditioning as centered

common=slope.common;weighted=slope.weighted;terminal=slope.terminal
HERE,PREFIX,sha,encode,iv,ep=slope.HERE,slope.PREFIX,slope.sha,slope.encode,slope.iv,slope.ep
N=slope.N;KEYS=slope.KEYS
NAME=PREFIX+'current_firstbridge_centered_C1.json';RECEIPT=PREFIX+'current_firstbridge_centered_C1_check.json'
GATE='actual_firstbridge_centered_C1_covers_and_updated_Rc_targets_executed'


def live_velocity(first,payload):
    cell=weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])[1][0][1]
    provenance=cell['original_whole_cell_source_and_cutoff_C1']['source_provenance'];c=first.ctx
    query=first.q_owner.query('bridge_first',iv(c,provenance['Z_box']),iv(c,provenance['coordinate_box']))
    roots=query['source']['roots'];packet=query['source']['packet'];E=roots[( 'E')][(0,0)]
    saved=payload['twelve_source_primitive_cover_records'][0]['genuine_same_source_conditional_paired_cover']
    assert encode({k:{'y%d_Z%d'%order:v.record() for order,v in rows.items()} for k,rows in roots.items()})==saved['original_root_rows']
    signed=first.owner.owner.signed_owner;values=[]
    for k in (0,1):
        row=slope.prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
        values.append(signed.leaf(slope.prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger).record())
    return encode(dict(source_family=first.family,native_source_bases=E.scale.bases,V=values[0],V_Z=values[1],
        original_source_provenance=provenance,all_original_root_rows_exactly_match_accepted_archive=True,
        original_firstbridge_packet_queried_without_route_or_producer=True))


def tightened_source(c,payload,velocity,sigma):
    trace=payload['twelve_source_primitive_cover_records'][0];saved=trace['genuine_same_source_conditional_paired_cover']
    bases=tuple(iv(c,v) for v in saved['native_source_log_bases']);ledger=dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    restore=lambda rec:slope.native_restore(rec,bases,ledger)
    assert velocity['source_family']==payload['source_family'] and tuple(iv(c,v)._mpi_ for v in velocity['native_source_bases'])==tuple(v._mpi_ for v in bases)
    roots={name:{(int(k[1]),int(k[-1])):restore(v) for k,v in rows.items()} for name,rows in saved['original_root_rows'].items()}
    old={k:restore(v) for k,v in trace['selected_primitive_C0_Z_covers'].items()}
    caps=centered.primitive_caps(c,roots,dict(log_actual_a_positive_lower=saved['original_log_a_positive_lower']),
        iv(c,saved['original_eta_log']),iv(c,saved['original_dstar_log']),sigma,dict(values=old))
    values=dict(old);decisions={}
    for key,capkey in (('A_Z','A'),('B_Z_over_Pstar','B')):
        cap=caps[capkey][centered.DZ]
        candidate=slope.prior.ScaledEnclosure(slope.prior.FormalScale(bases,offset=cap.log),c.mpf(('-1','1')),ledger) if cap.log is not None else old[key].scalar(0)
        before=weighted.upstream.baseline.magnitude_log(old[key]);after=weighted.upstream.baseline.magnitude_log(candidate)
        tighter=after is None or before is not None and after<before
        values[key]=candidate if tighter else old[key]
        decisions[key]=dict(old_log_absolute_upper=before,new_log_absolute_upper=after,strict_upper_reduction=tighter and before!=after)
    assert values['A'] is old['A'] and values['B_over_Pstar'] is old['B_over_Pstar']
    density=weighted.upstream.middle.density.density_Z_kernels(roots['E'][(0,0)],roots['E'][(0,1)],restore(velocity['V']),restore(velocity['V_Z']),values,N)
    coords=weighted.replay.inputs(c,payload)[0];cell=weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])[1][0][1]
    geometry=cell[weighted.accepted.row_fields('active_first_bridge')[0]]
    g=dict(width=common.restore_common_source(geometry['positive_true_log_radius_width'],coords),regular=iv(c,geometry['regular_true_log_radius_width_cover']),scalar_cover=iv(c,geometry['scalar_width_cover_used_only_for_directed_kernel_bounds']))
    contributions={};comparison={}
    for key,rate in slope.rc.RATES.items():
        mass=slope.rc.transfer.true_width_kernel(coords,g,rate)['mass']
        candidate=coords.rebase(density['Z_derivatives'][key],coords.family)*mass
        accepted=common.restore_common_source(cell['true_log_radius_signed_Z_contributions'][key],coords)
        before=weighted.upstream.baseline.magnitude_log(accepted);after=weighted.upstream.baseline.magnitude_log(candidate)
        tighter=after is None or before is not None and after<before
        contributions[key]=(candidate if tighter else accepted).record()
        comparison[key]=dict(old_log_absolute_upper=before,new_log_absolute_upper=after,strict_upper_reduction=tighter and before!=after)
    return dict(selected_Z_contributions=encode(contributions),record=dict(original_live_velocity_source=velocity,
        centered_same_original_source_primitive_proof=caps['record'],original_centered_phase_theorem=centered.exact_kernel_theorem(),
        selected_native_primitive_C0_Z=encode({k:v.record() for k,v in values.items()}),primitive_bound_comparisons=decisions,
        original_native_signed_density_Z={k:v.record() for k,v in density['Z_derivatives'].items()},
        firstbridge_contribution_comparisons=comparison,original_C0_source_primitives_and_contributions_unchanged=True,
        whole_original_period_and_cutoff_union_bounded_no_phase_sample=True,actual_new_phase_inverse_or_function_point_oracle_installed=False))


def replay_upstream(c,payload,tightened,pressure):
    changed=copy.deepcopy(payload);route=changed['complete_actual_pre_O2_tile_route']
    flat,rows=weighted.accepted.route_rows(route);coords=weighted.replay.inputs(c,payload)[0]
    jets={k:common.restore_common_source(flat['actual_correction_Z_enclosures'][k],coords) for k in KEYS}
    path=[]
    for label,row in rows:
        if label=='active_first_bridge':row['true_log_radius_signed_Z_contributions']=tightened['selected_Z_contributions']
        geo=row[weighted.accepted.row_fields(label)[0]]
        g=dict(width=common.restore_common_source(geo['positive_true_log_radius_width'],coords),regular=iv(c,geo['regular_true_log_radius_width_cover']),scalar_cover=iv(c,geo['scalar_width_cover_used_only_for_directed_kernel_bounds']))
        jets={k:slope.rc.transfer.true_width_kernel(coords,g,rate)['decay']*jets[k]+common.restore_common_source(row['true_log_radius_signed_Z_contributions'][k],coords) for k,rate in slope.rc.RATES.items()}
        path.append(dict(label=label,actual_right_correction_Z={k:v.record() for k,v in jets.items()}))
    # Keep the payload as contribution evidence only: its historical saved
    # output fields are not relabelled as new outputs.
    half=slope.native.HalfPstarCoordinates(c,coords.logP_squared,payload['source_family'])
    trace=slope.pre_slope_attribution(half,changed,pressure)
    return dict(jets=encode({k:v.record() for k,v in jets.items()}),trace=trace,
        record=dict(actual_new_twelve_chart_Z_transport=path,all_other11_chart_source_contributions_unchanged=True,
        actual_initial_flat_correction_source_zero=True,original_all_C0_P0_and_pressure_rate_zero_retained=True))


def updated_terminal(c,sd,upstream,primary,report):
    old=terminal.load_archive(sd['accepted_original_terminal_archive']);bd=terminal.load_archive(old['accepted_original_buffer_archive']);ad=terminal.load_archive(bd['accepted_original_axial_archive'])
    td=terminal.load_archive(bd['accepted_original_transition_archive'])
    coords=slope.native.HalfPstarCoordinates(c,iv(c,sd['common_directed_coordinate_theorem']['common_log_bases'][1]),sd['source_family'])
    amp=terminal.amplitude(coords,ad);primary=copy.deepcopy(primary);primary['actual_refined_upstream_Z']=upstream['jets']
    C=terminal.complete.restore_half_source(sd['slope_own_joint_C_at_slope_exit'],coords);CZ=terminal.complete.restore_half_source(sd['slope_own_joint_C_Z_at_slope_exit'],coords)
    output=slope.terminal_composition(coords,old,bd,ad,report,primary,C,CZ,amp,upstream['trace'])
    common_bases=(c.mpf(0),coords.logP_squared,c.mpf(0),c.mpf(0),c.mpf(0))
    lift=lambda rec:coords.rebase(slope.native_restore(rec,common_bases,coords.ledger),coords.family)
    all_width=coords.logP_squared/2+3;after=all_width-1
    jets={k:lift(upstream['jets'][k])*coords.decay(all_width,rate)+coords.scalar(iv(c,report['five_genuine_ordinary_Z_integral_contributions'][k]))*coords.decay(after,rate) for k,rate in slope.rc.RATES.items()}
    post=terminal.component_source_parts(coords,bd,ad,td)
    for name in ('original_axial_joint_source','original_buffer_joint_source','original_transition_joint_source'):
        for k in KEYS:jets[k]=jets[k]+post[name]['Z_derivatives'][k]
    values=terminal.restored(bd['actual_updated_Rd_Rc_correction_and_own_history']['actual_updated_Rc_correction_C0'],coords)
    targets=terminal.normalize(values,jets,terminal.complete.restore_half_source(encode(output['actual_joint_terminal_C']),coords),terminal.complete.restore_half_source(encode(output['actual_joint_terminal_C_Z']),coords),amp)
    output.update(actual_signed_five_terminal_defect_C0=terminal.native.records(targets['values']),actual_signed_five_terminal_defect_Z=terminal.native.records(targets['Z_derivatives']),
        actual_N_scaled_targets_C0=terminal.native.records({k:v*N for k,v in targets['values'].items()}),actual_N_scaled_targets_Z=terminal.native.records({k:v*N for k,v in targets['Z_derivatives'].items()}),
        actual_five_target_diagnostics=terminal.summary(targets),original_five_control_linear_response_enclosures=terminal.linear_control_enclosures(coords,targets,amp),
        unchanged_original_component_histories_used_for_other_four_rows=False,actual_other_four_Z_rows_recomposed_from_new_upstream_and_unchanged_sources=True)
    return output


@slope.rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/slope.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,slope,accepted['source_family']);common.attach_receipt(hashes,centered,accepted['source_family'])
    wm=json.loads((HERE/weighted.NAME).read_bytes());cm=json.loads((HERE/slope.source.NAME).read_bytes())
    sigma_manifest=json.loads((HERE/(PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
    bridge,seed=slope.rc.native.inlet.native_bridge_owner();archives=[]
    with slope.rc.native.inlet.CheckedSourceRuntime():
        owner=weighted.upstream.baseline.build_owner(bridge);first=owner.owner.owner.owner;c=owner.ctx
        sigma={int(k):centered.previous.read_cap(c,v) for k,v in sigma_manifest['sigma_global_derivative_log_caps'].items()}
        for archive in accepted['actual_original_slope_joint_source_archives']:
            sd=terminal.load_archive(archive);Z=sd['exact_Z_range']
            pressure=next(row for row in wm['actual_weighted_pressure_tiles'] if weighted.same_Z(row['exact_Z_range'],Z));payload=terminal.load_archive(pressure['original_source_archive'])
            velocity=live_velocity(first,payload);tight=tightened_source(c,payload,velocity,sigma);upstream=replay_upstream(c,payload,tight,pressure)
            primary=next(row for row in wm['genuine_original_O2_weighted_pressure_replays'] if row['ordered_source_cells']==2048 and weighted.same_Z(row['exact_Z_range'],Z))
            report=next(row for row in cm['original_complete_same_N_source_incoming_O2_integral_refinements'] if row['ordered_source_cells']==2048 and weighted.same_Z(row['exact_Z_range'],Z))
            output=updated_terminal(c,sd,upstream,primary,report)
            result=dict(source_family=owner.family,candidate_N=N,exact_Z_range=Z,firstbridge_centered_source=tight['record'],selected_firstbridge_Z_contributions=tight['selected_Z_contributions'],
                actual_updated_pre_slope_Z=upstream['jets'],actual_updated_pre_slope_C1_attribution=upstream['trace'],actual_twelve_chart_transport=upstream['record'],actual_updated_Rc_terminal_defects=output,
                accepted_original_slope_archive=archive,accepted_original_pre_slope_archive=pressure['original_source_archive'],
                actual_pre_slope_phase_correlation_nonlinear_controls_whole_Z_axis_global_N_heat_stress_recursion_admitted=False,**dict.fromkeys(common.current.FLAGS,False))
            raw=json.dumps(encode(result),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            tag='positive' if slope.Fraction(Z[0])>0 else 'negative';name=PREFIX+'current_firstbridge_centered_C1_'+tag+'.json.gz';(HERE/name).write_bytes(compressed)
            hashes[name]=sha(name);archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),exact_Z_range=Z))
            print('Firstbridge centered derivative source and new Rc targets:',Z,'executed',flush=True)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Accepted source closure disagreement: '+name)
            hashes[name]=digest
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],candidate_N=N,actual_firstbridge_centered_C1_archives=archives,
        only_two_firstbridge_native_packets_queried=True,accepted_other11_charts_slope_and_post_slope_producers_or_inverse_solvers_not_rerun=True,
        actual_global_function_control_stress_recursion_admitted=False,**dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
