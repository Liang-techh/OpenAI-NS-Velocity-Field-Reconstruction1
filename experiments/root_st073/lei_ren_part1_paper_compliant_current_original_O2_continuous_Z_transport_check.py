"""Check complete continuous rectangles, own masses and actual-inlet C1 wiring."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport as current
import lei_ren_part1_paper_compliant_current_original_O2_density_Z_integrals_check as previous

base=current.base;ep=current.ep;prev=current.current;KEYS=current.KEYS


def iv(c,value):return prev.interval(c,value)


def overlap(a,b):return max(ep(a)[0],ep(b)[0])<=min(ep(a)[1],ep(b)[1])


def original_endpoint(c,level,z):
    """Independent endpoint history formulas from accepted defining masses."""
    node=level['ordered_nodes'][-1];C=1/(1+z*z)
    H,D,P=(iv(c,node[k]) for k in ('H','D','P'))
    # Retain the accepted positive native Pstar inverse via the same formal
    # number system; only the endpoint formulas/defining masses are separate.
    bases=(c.exp(40)+11,c.mpf(0),c.mpf(0),c.mpf(0),c.mpf(0))
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    scalar=lambda v:base.prior.ScaledEnclosure(base.prior.FormalScale(bases),v,ledger)
    inv=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,(-1,0,0,0,0)),c.mpf(1),ledger)
    inv2=current.prior.slow.square(inv)
    values=dict(m=inv*(4*z),h=scalar(C*H),k=inv*(4*z*C*H),
        e=inv2*(16*z*z)-scalar(C*C*D),p=scalar(C*C*P))
    jets=dict(m=inv*4,h=scalar(-2*z*C*C*H),k=inv*(4*(1-z*z)*C*C*H),
        e=inv2*(32*z)+scalar(4*z*C**3*D),p=scalar(-4*z*C**3*P))
    return ({k:prev.bounded(v) for k,v in values.items()},
            {k:prev.bounded(v) for k,v in jets.items()})


def incoming_checks(c,record):
    """Nonzero synthetic data test the API; they are not physical inlets."""
    incoming={k:c.mpf((str(i+1),str(i+2))) for i,k in enumerate(KEYS)}
    incomingZ={k:c.mpf((str(-i-2),str(-i-1))) for i,k in enumerate(KEYS)}
    args=dict(incoming=incoming,incoming_Z=incomingZ,source_family=record['source_family'],
        original_P0_datum_sha256=record['original_P0_datum_sha256'])
    result=current.apply_actual_incoming(c,record,**args)
    assert result['actual_incoming_histories_not_set_to_zero'] and result['P0_not_reset_or_added_to_cumulative_pressure']
    assert not result['functional_terminal_identity_solved'] and all(result[k] is False for k in current.FLAGS)
    for key,rate in prev.five.RATES.items():
        decay=c.exp(-c.mpf(rate))
        old=iv(c,record['source_defined_original_histories_at_y1'][key])
        oldZ=iv(c,record['source_defined_original_history_Z_at_y1'][key])
        assert ep(result['original_histories_at_y1'][key])==ep(old)
        assert ep(result['original_history_Z_at_y1'][key])==ep(oldZ)
        defect=decay*incoming[key]+iv(c,record['five_original_C0_integral_contributions'][key])
        defectZ=decay*incomingZ[key]+iv(c,record['five_genuine_ordinary_Z_integral_contributions'][key])
        assert ep(result['transported_incoming_and_modulation_defects'][key])==ep(defect)
        assert ep(result['transported_defect_Z'][key])==ep(defectZ)
        assert ep(result['own_five_histories_at_y1'][key])==ep(old+defect)
        assert ep(result['own_five_ordinary_Z_at_y1'][key])==ep(oldZ+defectZ)
    rejected=0
    for bad in (dict(args,incoming={'m':c.mpf(1)}),dict(args,incoming_Z={'p':c.mpf(1)}),
                dict(args,original_P0_datum_sha256='different-datum')):
        try:current.apply_actual_incoming(c,record,**bad)
        except ValueError:rejected+=1
        else:raise AssertionError('Missing incoming/source identity was admitted')
    for bad_record in (dict(record,exact_y_window=['.5','1']),
            dict(record,own_rates={**prev.five.RATES,'p':1}),
            dict(record,original_P0_datum_sha256='different-datum')):
        try:current.apply_actual_incoming(c,bad_record,**args)
        except ValueError:rejected+=1
        else:raise AssertionError('Changed recorded window/rate/datum was admitted')
    assert rejected==6
    return dict(passed=True,nonzero_affine_value_and_ordinary_Z_memory_cases=10,
        missing_incoming_or_changed_datum_rejections=rejected,
        synthetic_incoming_wiring_test_not_physical_source_result=True)


def rectangle_checks(c,record,owner,fixed=None):
    count=record['ordered_source_cells'];zl,zh=record['exact_Z_range']
    level,phase=owner.parent.parent.parent.levels[count]
    assert record['source_family']==owner.family and record['explicit_candidate_N']==7
    assert record['exact_y_window']==['0','1'] and record['own_rates']==prev.five.RATES
    assert record['normalized_own_units']==prev.five.UNITS
    assert record['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
    assert record['entire_continuous_Z_interval_not_samples'] and record['actual_incoming_histories_not_set_to_zero']
    assert not record['Z_interval_terminal_identities_or_global_controls_installed'] and not record['current_whole_N_selected']
    assert record['original_nonzero_Z_inlet_contract']==owner.inlet_contract
    rows=record['whole_original_density_and_density_Z_source_records'];assert len(rows)==count
    total_fields=('five_original_C0_integral_contributions','five_genuine_ordinary_Z_integral_contributions',
                  'five_unmodulated_original_integral_contributions','five_unmodulated_original_Z_integral_contributions')
    hull_fields=('five_signed_C0_density_hulls','five_genuine_signed_ordinary_Z_density_hulls',
                 'five_original_unmodulated_density_hulls','five_original_unmodulated_density_Z_hulls')
    add_fields=('five_C0_contributions','five_ordinary_Z_contributions',
                'five_original_unmodulated_contributions','five_original_unmodulated_Z_contributions')
    sums=[{k:c.mpf(0) for k in KEYS} for unused in range(4)]
    inverses=mass_checks=fixed_intersections=0;p=mp.mp.clone();p.dps=c.dps+40
    for i,(row,phase_cell) in enumerate(zip(rows,phase['whole_source_cells'],strict=True)):
        source=row['source'];original=source['original_source_function_coefficient_record']
        assert source['ordered_source_cell_level']==count and source['ordered_source_cell_index']==i
        assert original['exact_y_cell']==[str(i)+'/'+str(count),str(i+1)+'/'+str(count)]
        assert original['exact_Z_range']==[zl,zh] and original['whole_source_ranges_not_field_points']
        assert row['entire_source_y_Z_rectangle_enclosed'] and row['phase_and_masses_Z_independent']
        phases=phase_cell['true_common_N_phase_boxes'];proofs=row['actual_source_phase_held_Z_records']
        assert proofs and len(proofs)==len(source['conditioned_source_piece_geometry'])*len(phases)
        split=source['overlapping_source_coordinate_split']
        if split is not None:
            assert split['entire_native_logabsu_source_range_covered']
            assert split['overlapping_original_small_and_large_u_coordinates']
        for j,proof in enumerate(proofs):
            assert proof['source_piece_index']==j//len(phases) and proof['actual_true_phase_cover']==phases[j%len(phases)]
            target=iv(c,proof['actual_true_phase_cover'])
            assert previous.fixture.contains(iv(c,proof['original_inverse']['phase_image']),target)
            assert previous.fixture.contains(iv(c,proof['original_derivative_inverse']['phase_image']),target)
            previous.carrier_contract(c,proof['correlated_original_p2_Z_carrier'],source)
            contract=proof['original_slow_Z_contract']
            assert contract['original_phase_held_Z_not_finite_difference'] and contract['source_Z_derivatives_not_interval_selector_derivatives']
            assert contract['q_Z_and_dstar_Z_exact_zero_source_identity'] and proof['density_Z_is_genuine_ordinary_source_derivative']
            assert all(v['encloses_original_source_function'] and not v['point_value_selected'] for v in proof['native_original_primitive_Z_enclosures'].values())
            inverses+=1
        for key,rate in prev.five.RATES.items():
            mass=iv(c,row['positive_own_rate_final_endpoint_masses'][key])
            left,right=p.mpf(i)/count,p.mpf(i+1)/count;r=p.mpf(rate)
            exact=right-left if r==0 else (p.exp(-r*(1-right))-p.exp(-r*(1-left)))/r
            assert 0<ep(mass)[0]<=exact<=ep(mass)[1];mass_checks+=1
            for total,hulls,adds in zip(sums,hull_fields,add_fields,strict=True):
                hull=iv(c,row[hulls][key]);add=iv(c,row[adds][key])
                assert all(mp.isfinite(v) for v in ep(hull)) and ep(hull*mass)==ep(add)
                total[key]+=add
            if fixed is not None:
                before=fixed['whole_original_density_and_density_Z_source_records'][i]
                for field in hull_fields[:2]:
                    assert overlap(iv(c,row[field][key]),iv(c,before[field][key]))
                    fixed_intersections+=1
    for total,field in zip(sums,total_fields,strict=True):
        for key in KEYS:assert ep(total[key])==ep(iv(c,record[field][key]))
    zlo,zhi=(base.point.pressure.exact_Z(v) for v in (zl,zh))
    z=c.mpf((ep(c.mpf(int(zlo.p))/int(zlo.q))[0],ep(c.mpf(int(zhi.p))/int(zhi.q))[1]))
    endpoint,endpointZ=original_endpoint(c,level,z)
    for key in KEYS:
        assert overlap(endpoint[key],iv(c,record['source_defined_original_histories_at_y1'][key]))
        assert overlap(endpointZ[key],iv(c,record['source_defined_original_history_Z_at_y1'][key]))
    inlet=record['source_defined_original_inlet']
    assert inlet['microscopic_positive_inverse_Pstar_retained'] and inlet['whole_Z_rational_source_functions_not_axial_samples']
    assert not inlet['native_original_inlet_sources']['m']['exact_zero']
    assert not inlet['native_original_inlet_Z_sources']['m']['exact_zero']
    incoming=incoming_checks(c,record)
    return dict(passed=True,exact_Z_range=[zl,zh],ordered_source_cells=count,
        complete_source_cells=count,true_original_phase_derivative_pieces=inverses,
        independent_positive_own_mass_checks=mass_checks,
        accepted_fixed_Z_density_cover_intersections=fixed_intersections,
        original_endpoint_history_and_ordinary_Z_independent_defining_mass_intersections=10,
        incoming_C1_transport=incoming,
        five_C0_widths={k:ep(iv(c,record[total_fields[0]][k]))[1]-ep(iv(c,record[total_fields[0]][k]))[0] for k in KEYS},
        five_Z_widths={k:ep(iv(c,record[total_fields[1]][k]))[1]-ep(iv(c,record[total_fields[1]][k]))[0] for k in KEYS})


def run():
    begin=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and all(manifest[k] is False for k in current.FLAGS)
    assert manifest['actual_original_nonzero_Z_inlet_functions_installed']
    assert manifest['explicit_actual_incoming_value_and_derivative_interface_installed']
    assert not manifest['functional_terminal_identity_solved']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.OriginalO2ContinuousZTransport();c=owner.c
    saved=json.loads(gzip.decompress((current.HERE/current.prior.NAME).read_bytes()))
    fixed={r['ordered_source_cells']:r for r in saved['actual_original_five_density_Z_integral_refinements']}
    summaries=[];archive_hashes={};ranges={}
    with mp.workdps(c.dps+40):
        for archive in manifest['continuous_Z_integral_archives']:
            path=current.HERE/archive['filename'];compressed=path.read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            record=json.loads(raw);assert record['exact_Z_range']==archive['exact_Z_range']
            count=record['ordered_source_cells'];assert count==archive['ordered_source_cells']
            positive=record['exact_Z_range'][0]=='9/25'
            summaries.append(rectangle_checks(c,record,owner,fixed[count] if positive else None))
            for field in ('five_original_C0_integral_contributions','five_genuine_ordinary_Z_integral_contributions',
                          'five_unmodulated_original_integral_contributions','five_unmodulated_original_Z_integral_contributions'):
                assert record[field]==archive[field]
                for key in KEYS:ranges.setdefault((tuple(record['exact_Z_range']),field,key),[]).append(iv(c,record[field][key]))
            archive_hashes[archive['filename']]=current.sha(archive['filename'])
            print('Continuous whole-Z source, inlet, original endpoint and explicit C1 history PASS',record['exact_Z_range'],count,flush=True)
        for values in ranges.values():assert max(ep(v)[0] for v in values)<=min(ep(v)[1] for v in values)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        continuous_original_Z_rectangle_checks=summaries,original_inlet_source_contract=owner.inlet_contract,
        original_coordinates_velocity_inverse_and_zero_mass_source_assignments_guarded=True,
        complete_value_Z_and_original_history_refinements_have_common_intersection=True,
        independent_original_scalar_derivatives_reused_from_accepted_predecessor=True,
        crossing_Z_functional_terminal_repair_global_N_stress_recursion_not_claimed=True,
        functional_terminal_identity_solved=False,**dict.fromkeys(current.FLAGS,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),**archive_hashes,
            Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-begin,
        scope=manifest['scope'])
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Continuous original O2 five-history C1 transport PASS',flush=True);return result


if __name__=='__main__':run()
