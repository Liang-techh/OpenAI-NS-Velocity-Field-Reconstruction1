"""Check genuine N1024 phase evidence and real upstream affine C1 attachment."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_source_incoming_common_N as source
import lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport_check as accepted

base=source.base;ep=source.ep;KEYS=source.KEYS;iv=source.interval
TOTALS=('five_original_C0_integral_contributions','five_genuine_ordinary_Z_integral_contributions',
    'five_unmodulated_original_integral_contributions','five_unmodulated_original_Z_integral_contributions')
HULLS=('five_signed_C0_density_hulls','five_genuine_signed_ordinary_Z_density_hulls',
    'five_original_unmodulated_density_hulls','five_original_unmodulated_density_Z_hulls')
ADDS=('five_C0_contributions','five_ordinary_Z_contributions',
    'five_original_unmodulated_contributions','five_original_unmodulated_Z_contributions')


def same_factored(c,record,value):
    """Compare the full cover, including its factor powers and directed offset."""
    scale=record['formal_positive_scale']
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero
    assert not record['point_value_selected'] and record['encloses_original_source_function']


def actual_attachment_checks(owner,record):
    c=owner.c;coordinates=owner.coordinates
    attached=record['actual_source_incoming_and_own_history_transport']
    binding=attached['actual_incoming_source_binding']
    assert attached['source_family']==owner.family and attached['candidate_N']==source.N
    assert attached['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
    assert attached['exact_Z_range']==record['exact_Z_range']
    assert attached['normalized_own_units']==source.previous.five.UNITS
    assert binding['source_family']==owner.family and binding['candidate_N']==source.N
    assert binding['manifest']==source.middle.NAME and binding['manifest_sha256']==source.sha(source.middle.NAME)
    assert binding['receipt']==source.middle.RECEIPT and binding['receipt_sha256']==source.sha(source.middle.RECEIPT)
    assert binding['source_record_key']=='whole_Z' and binding['actual_upstream_Z_cover']==owner.incoming['Z_box']
    assert binding['same_exact_incoming_function_graph']==owner.exact_binding
    assert binding['input_is_correction_not_original_plus_correction']
    assert binding['wider_whole_Z_cover_used_on_subtile_not_a_refined_function_value']
    assert binding['common_coordinates']==base.encoded(coordinates.record())
    # Independent accepted affine operator: append the actual one-unit O2
    # increments, then apply the complete upstream function covers explicitly.
    modulation={k:coordinates.scalar(iv(c,record[TOTALS[0]][k])) for k in KEYS}
    modulationZ={k:coordinates.scalar(iv(c,record[TOTALS[1]][k])) for k in KEYS}
    operator=source.middle.history.C1DuhamelOperator(coordinates)
    operator.append(c.mpf(1),modulation,modulationZ,owner.family)
    propagated=operator.apply(owner.incoming_values,owner.incoming_Z,owner.family)
    originals={k:coordinates.scalar(iv(c,record['source_defined_original_histories_at_y1'][k])) for k in KEYS}
    originalsZ={k:coordinates.scalar(iv(c,record['source_defined_original_history_Z_at_y1'][k])) for k in KEYS}
    for k in KEYS:
        for field,value in (
            ('actual_source_incoming_C0',owner.incoming_values[k]),
            ('actual_source_incoming_Z',owner.incoming_Z[k]),
            ('genuine_O2_window_C0',modulation[k]),('genuine_O2_window_Z',modulationZ[k]),
            ('original_background_history_C0',originals[k]),('original_background_history_Z',originalsZ[k]),
            ('propagated_actual_correction_C0',propagated['values'][k]),
            ('propagated_actual_correction_Z',propagated['Z_derivatives'][k]),
            ('actual_original_plus_propagated_correction_C0',originals[k]+propagated['values'][k]),
            ('actual_original_plus_propagated_correction_Z',originalsZ[k]+propagated['Z_derivatives'][k])):
            same_factored(c,attached[field][k],value)
        raw0=source.restore_common_source(owner.incoming['actual_original_inlet_to_O2_inlet_correction_C0'][k],coordinates)
        rawZ=source.restore_common_source(owner.incoming['actual_original_inlet_to_O2_inlet_correction_Z'][k],coordinates)
        same_factored(c,attached['actual_source_incoming_C0'][k],raw0)
        same_factored(c,attached['actual_source_incoming_Z'][k],rawZ)
    same_factored(c,attached['separate_P0'],owner.P0)
    same_factored(c,attached['separate_P0_Z'],owner.P0_Z)
    same_factored(c,attached['actual_absolute_pressure'],owner.P0+(originals['p']+propagated['values']['p']))
    same_factored(c,attached['actual_absolute_pressure_Z'],owner.P0_Z+(originalsZ['p']+propagated['Z_derivatives']['p']))
    assert operator.coefficients['p'].scale.powers==(0,0,0,0,0)
    assert ep(operator.coefficients['p'].coefficient)==(1,1) and ep(operator.coefficients['p'].scale.offset)==(0,0)
    assert attached['original_background_added_once_after_correction_transport']
    assert attached['rate_zero_pressure_correction_retains_entire_actual_upstream_memory']
    assert attached['common_factored_context_ledger_retained_no_enormous_scalar_cast']
    assert attached['no_incoming_zero_fallback_or_N7_integral_reused']
    assert attached['full_Rc_route_or_functional_terminal_identity_solved'] is False
    assert all(attached[k] is False for k in source.current.FLAGS)
    bad_records=(dict(record,explicit_candidate_N=7),dict(record,source_family={}),
        dict(record,original_P0_datum_sha256='different-datum'),dict(record,exact_y_window=['.5','1']),
        dict(record,normalized_own_units={}),dict(record,own_rates={**source.previous.five.RATES,'p':1}),
        dict(record,exact_Z_range=['1','2']))
    rejected=0
    for bad in bad_records:
        try:owner.compose(bad)
        except ValueError:rejected+=1
        else:raise AssertionError('Changed N/family/datum/window/unit/rate/Z domain was admitted')
    try:owner.phase_level(record['ordered_source_cells'],7)
    except ValueError:rejected+=1
    else:raise AssertionError('N7 upstream mismatch was admitted')
    assert rejected==8
    return dict(passed=True,complete_actual_source_C0_Z_cover_comparisons=10,
        independently_accepted_affine_operator_C0_Z_matches=10,changed_contract_rejections=rejected,
        original_background_added_once=True,P0_separate_and_rate_zero_pressure_memory_preserved=True,
        upstream_whole_Z_cover_is_conservative_not_refined_subtile_point_value=True)


def rectangle_checks(owner,record,old,archive_hashes):
    c=owner.c;count=record['ordered_source_cells'];zl,zh=record['exact_Z_range']
    assert record['source_family']==owner.family and record['explicit_candidate_N']==source.N
    source.current.require_identity(record,owner.family,owner.family['datum_enclosure_sha256'])
    assert record['phase_cover_regenerated_for_actual_upstream_N']
    assert record['accepted_N_independent_y_source_cache_phase_candidate']==7
    assert record['entire_continuous_Z_interval_not_samples']
    assert record['complete_original_y_window_all_source_and_phase_pieces_integrated']
    level,phase,origin=owner.phase_level(count,source.N)
    encoded_origin=base.encoded(origin)
    assert record['fresh_same_N_true_radius_phase_origin']==encoded_origin
    assert record['original_nonzero_Z_inlet_contract']==owner.parent.inlet_contract
    sums=[{k:c.mpf(0) for k in KEYS} for unused in range(4)]
    p=mp.mp.clone();p.dps=c.dps+40
    expected_index=0;pieces=mass_checks=0
    for archive in record['whole_original_density_and_density_Z_source_archives']:
        name=archive['filename'];compressed=(source.HERE/name).read_bytes();raw=gzip.decompress(compressed)
        assert len(compressed)==archive['compressed_bytes']<100*1024*1024
        assert len(raw)==archive['uncompressed_bytes'] and hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
        chunk=json.loads(raw);left,right=archive['source_cell_index_range']
        assert chunk['source_cell_index_range']==[left,right] and left==expected_index
        assert chunk['source_cell_count']==count and right<=count and right-left<=512
        assert chunk['source_family']==owner.family and chunk['candidate_N']==source.N
        assert chunk['exact_Z_range']==[zl,zh] and chunk['fresh_same_N_true_radius_phase_origin']==encoded_origin
        assert len(chunk['records'])==right-left
        for i,row in enumerate(chunk['records'],start=left):
            original=row['source']['original_source_function_coefficient_record']
            assert row['source']['ordered_source_cell_level']==count and row['source']['ordered_source_cell_index']==i
            assert original['exact_y_cell']==[str(i)+'/'+str(count),str(i+1)+'/'+str(count)]
            assert original['exact_Z_range']==[zl,zh] and original['whole_source_ranges_not_field_points']
            assert row['entire_source_y_Z_rectangle_enclosed'] and row['phase_and_masses_Z_independent']
            phases=phase['whole_source_cells'][i]['true_common_N_phase_boxes']
            proofs=row['actual_source_phase_held_Z_records']
            assert len(proofs)==len(row['source']['conditioned_source_piece_geometry'])*len(phases)>0
            split=row['source']['overlapping_source_coordinate_split']
            if split is not None:
                assert split['entire_native_logabsu_source_range_covered']
                assert split['overlapping_original_small_and_large_u_coordinates']
            for j,proof in enumerate(proofs):
                assert proof['source_piece_index']==j//len(phases)
                assert proof['actual_true_phase_cover']==phases[j%len(phases)]
                target=iv(c,proof['actual_true_phase_cover'])
                assert accepted.previous.fixture.contains(iv(c,proof['original_inverse']['phase_image']),target)
                assert accepted.previous.fixture.contains(iv(c,proof['original_derivative_inverse']['phase_image']),target)
                accepted.previous.carrier_contract(c,proof['correlated_original_p2_Z_carrier'],row['source'])
                contract=proof['original_slow_Z_contract']
                assert contract['original_phase_held_Z_not_finite_difference']
                assert contract['source_Z_derivatives_not_interval_selector_derivatives']
                assert contract['q_Z_and_dstar_Z_exact_zero_source_identity']
                assert proof['density_Z_is_genuine_ordinary_source_derivative']
                assert all(v['encloses_original_source_function'] and not v['point_value_selected']
                    for v in proof['native_original_primitive_Z_enclosures'].values())
                pieces+=1
            for key,rate in source.previous.five.RATES.items():
                mass=iv(c,row['positive_own_rate_final_endpoint_masses'][key])
                l,r=p.mpf(i)/count,p.mpf(i+1)/count;rate=p.mpf(rate)
                exact=r-l if rate==0 else (p.exp(-rate*(1-r))-p.exp(-rate*(1-l)))/rate
                assert 0<ep(mass)[0]<=exact<=ep(mass)[1];mass_checks+=1
                for total,hulls,adds in zip(sums,HULLS,ADDS,strict=True):
                    hull=iv(c,row[hulls][key]);add=iv(c,row[adds][key])
                    assert all(mp.isfinite(v) for v in ep(hull)) and ep(hull*mass)==ep(add)
                    total[key]+=add
        archive_hashes[name]=source.sha(name);expected_index=right
        del chunk,raw,compressed
    assert expected_index==count and record['every_complete_source_cell_saved_losslessly']
    for total,field in zip(sums,TOTALS,strict=True):
        for k in KEYS:assert ep(total[k])==ep(iv(c,record[field][k]))
    # Original source integrals do not depend on N; exact old accepted totals
    # provide a stronger cross-frequency test than rechecking old ancestors.
    for field in TOTALS[2:]:
        for k in KEYS:assert ep(iv(c,record[field][k]))==ep(iv(c,old[field][k]))
    zlo,zhi=(base.point.pressure.exact_Z(v) for v in (zl,zh))
    z=c.mpf((ep(c.mpf(int(zlo.p))/int(zlo.q))[0],ep(c.mpf(int(zhi.p))/int(zhi.q))[1]))
    endpoint,endpointZ=accepted.original_endpoint(c,level,z)
    for k in KEYS:
        assert accepted.overlap(endpoint[k],iv(c,record['source_defined_original_histories_at_y1'][k]))
        assert accepted.overlap(endpointZ[k],iv(c,record['source_defined_original_history_Z_at_y1'][k]))
    attachment=actual_attachment_checks(owner,record)
    return dict(passed=True,exact_Z_range=[zl,zh],candidate_N=source.N,complete_source_cells=count,
        fresh_actual_phase_ordinary_Z_source_pieces=pieces,independent_positive_own_mass_checks=mass_checks,
        original_unmodulated_integrals_exactly_match_N7_source_cache=True,
        source_endpoint_ordinary_jet_defining_mass_intersections=10,actual_incoming_attachment=attachment,
        five_C0_widths={k:ep(iv(c,record[TOTALS[0]][k]))[1]-ep(iv(c,record[TOTALS[0]][k]))[0] for k in KEYS},
        five_Z_widths={k:ep(iv(c,record[TOTALS[1]][k]))[1]-ep(iv(c,record[TOTALS[1]][k]))[0] for k in KEYS})


def run():
    begin=time.monotonic();manifest=json.loads((source.HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert manifest['genuine_N1024_O2_phase_and_integrals_not_N7_relabelling']
    assert manifest['accepted_upstream_whole_Z_function_covers_restored_as_factored_covers']
    assert manifest['upstream_cover_widths_not_improved_by_downstream_O2_integration']
    assert manifest['all_chart_terminal_matching_or_global_N_or_true_recursion_admitted'] is False
    assert all(manifest[k] is False for k in source.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert source.sha(name)==digest,name
    owner=source.OriginalO2SourceIncomingCommonN();c=owner.c
    assert manifest['exact_upstream_source_function_incoming_binding']==owner.exact_binding
    old_manifest=json.loads((source.HERE/source.current.NAME).read_bytes())
    old={(tuple(r['exact_Z_range']),r['ordered_source_cells']):r for r in old_manifest['continuous_Z_integral_archives']}
    summaries=[];archive_hashes={};ranges={}
    with mp.workdps(c.dps+40):
        for record in manifest['original_complete_same_N_source_incoming_O2_integral_refinements']:
            before=old[(tuple(record['exact_Z_range']),record['ordered_source_cells'])]
            summaries.append(rectangle_checks(owner,record,before,archive_hashes))
            for field in TOTALS:
                for k in KEYS:ranges.setdefault((tuple(record['exact_Z_range']),field,k),[]).append(iv(c,record[field][k]))
            print('Same-N1024 genuine source/phase integrals and actual upstream C1 attachment PASS',record['exact_Z_range'],record['ordered_source_cells'],flush=True)
        for values in ranges.values():assert max(ep(v)[0] for v in values)<=min(ep(v)[1] for v in values)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        same_N_actual_upstream_original_O2_checks=summaries,
        complete_N1024_true_phase_and_source_cell_evidence_checked=True,
        accepted_original_scalar_implicit_derivative_checks_reused=True,
        all_original_and_modulation_C0_Z_refinements_have_common_intersection=True,
        exact_13_cell_upstream_function_graph_binding=owner.exact_binding,
        full_Rc_terminal_identity_or_global_N_or_stress_or_recursion_claimed=False,
        **dict.fromkeys(source.current.FLAGS,False),
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),**archive_hashes,
            Path(__file__).name:source.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-begin,scope=manifest['scope'])
    (source.HERE/source.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual same-N original O2 incoming C1 history attachment PASS',flush=True);return result


if __name__=='__main__':run()
