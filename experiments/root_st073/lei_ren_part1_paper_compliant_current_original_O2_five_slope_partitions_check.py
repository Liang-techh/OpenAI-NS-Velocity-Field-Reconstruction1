"""Exact source-bin restriction and original five-cell C0/Z memory checks."""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_five_slope_partitions as current

ep=current.ep;iv=current.iv;KEYS=current.KEYS


def same(record,value):
    scale=record['formal_positive_scale'];c=value.ctx
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero and not record['point_value_selected'] and record['encloses_original_source_function']


def overlaps(a,b):
    assert a.ctx is b.ctx and a.scale.bases is b.scale.bases and a.ledger is b.ledger
    anchor=a.ctx.mpf(max(ep(a.scale.evaluate())[1],ep(b.scale.evaluate())[1]))
    left=a.coefficient*a.bounded_exp(a.scale.evaluate()-anchor)
    right=b.coefficient*b.bounded_exp(b.scale.evaluate()-anchor)
    return max(ep(left)[0],ep(right)[0])<=min(ep(left)[1],ep(right)[1])


def finite_overlap(a,b):return max(ep(a)[0],ep(b)[0])<=min(ep(a)[1],ep(b)[1])


def independent_exact_partition_and_weights(owner,saved):
    c=owner.ctx;seen=[];masses=0;densityrows=0;decays=0
    assert saved['exact_original_slope_edges']==list(map(str,current.EDGES))
    assert len(saved['actual_original_five_slope_cell_records'])==5
    for cell,label,a,b in zip(saved['actual_original_five_slope_cell_records'],current.LABELS,current.EDGES[:-1],current.EDGES[1:],strict=True):
        assert cell['label']==label and cell['chart']=='O2_slope' and cell['candidate_N']==1024
        pieces=cell['source_bin_piece_records'];lo=a
        for piece in pieces:
            index=piece['original_source_bin'];left,right=map(Fraction,piece['exact_y_piece'])
            assert left==lo and a<=left<right<=b and Fraction(index,64)<=left<right<=Fraction(index+1,64)
            assert Fraction(piece['exact_piece_width'])==right-left
            assert piece['source_density_cover_domain']==owner.source_cells[index]['exact_y_cell']
            assert piece['same_source_whole_bin_density_hull_valid_on_contained_piece']
            assert piece['actual_signed_whole_function_C0_Z_not_point_or_quadrature']
            width=current.rational(c,right-left)
            for key,r in current.five.RATES.items():
                rate=c.mpf(r)
                # Independent closed exponential integral, not exp_average.
                local=width if r==0 else -c.expm1(-rate*width)/rate
                suffix=c.exp(-rate*current.rational(c,b-right))
                final=local*c.exp(-rate*current.rational(c,1-right))
                assert finite_overlap(local,iv(c,piece['positive_piece_local_endpoint_masses'][key]))
                assert finite_overlap(final,iv(c,piece['positive_piece_global_endpoint_masses'][key]))
                assert finite_overlap(suffix,iv(c,piece['piece_to_partition_endpoint_decay'][key]));masses+=2
            for j in range(4):
                for key in KEYS:
                    actual=iv(c,owner.source_cells[index]['normalized_four_density_unions'][j][key])
                    expected=(owner.Z_unit if j==1 else owner.coordinates.scalar(1))*actual
                    assert overlaps(expected,owner.densities[index][j][key]);densityrows+=1
            phase=current.five.pressure.phase_boxes(c,owner.origin,1024,current.rational(c,left),current.rational(c,right))
            assert current.encode(phase)==piece['actual_phase_boxes_same_global_origin']
            lo=right;seen.append((index,left,right))
        assert lo==b
        for key,r in current.history.RATES.items():
            record=cell['actual_cell_C1_operator']['incoming_C0_Z_decay_coefficients'][key]
            same(record,owner.coordinates.decay(current.rational(c,b-a),r));decays+=1
    assert len(seen)==68 and saved['exact_source_bin_piece_count']==68
    for index in range(64):
        local=sorted((a,b) for i,a,b in seen if i==index)
        assert local[0][0]==Fraction(index,64) and local[-1][1]==Fraction(index+1,64)
        assert all(a[1]==b[0] for a,b in zip(local,local[1:]))
        assert sum((b-a for a,b in local),Fraction(0))==Fraction(1,64)
    assert len([q for q in seen if q[0]==7])==2 and len([q for q in seen if q[0]==8])==3 and len([q for q in seen if q[0]==9])==2
    return dict(passed=True,exact_original_slope_partitions=5,disjoint_source_measure_pieces=68,
        original_source_bins_complete_exactly_once=64,independent_local_global_positive_mass_comparisons=masses,
        retained_actual_source_density_C0_Z_and_background_rows=densityrows,independent_partition_own_rate_decays=decays,
        actual_global_radius_origin_retained_for_every_piece=True)


def actual_five_source_operator_replay(owner,saved):
    assert saved['accepted_actual_N1024_density_and_upstream_binding']==current.encode(owner.binding)
    live=owner.integrate();assert current.encode(live['report'])==saved
    c=owner.ctx;serial=live['correction'];composite=live['cumulative'].apply(owner.incoming['values'],owner.incoming['Z_derivatives'],owner.family)
    for field in ('values','Z_derivatives'):
        for key in KEYS:assert overlaps(serial[field][key],composite[field][key])
    old=owner.source_owner.saved;overlap_count=0
    oldfields=('changed_C0_integral','changed_ordinary_Z_integral','original_unmodulated_C0_integral','original_unmodulated_ordinary_Z_integral')
    convert=lambda row:owner.source_owner.native_to_common(current.preceding.driver.integrals.restore_value(owner.source_owner.atlas,row))
    for j,field in enumerate(oldfields):
        for key in KEYS:
            assert overlaps(live['global_totals'][j][key],convert(owner.direct[field][key]));overlap_count+=1
    for j,field in enumerate(('values','Z_derivatives')):
        for key in KEYS:
            assert overlaps(serial[field][key],convert(old['propagated_actual_direct_C0_Z_corrections'][j][key]));overlap_count+=1
            assert overlaps(live['background'][field][key],convert(owner.direct['source_defined_original_histories_at_y1' if j==0 else 'source_defined_original_history_Z_at_y1'][key]));overlap_count+=1
    pressurechecks=0
    for cell in saved['actual_original_five_slope_cell_records']:
        restore=lambda record:current.preceding.driver.upstream.restore_common_source(record,owner.coordinates)
        for suffix,P0 in (('C0',owner.P0),('Z',owner.P0_Z)):
            bg=restore(cell['original_right_background_'+suffix]['p']);correction=restore(cell['actual_right_correction_'+suffix]['p'])
            same(cell['actual_right_own_history_'+suffix]['p'],bg+correction)
            same(cell['actual_absolute_pressure_'+suffix],P0+(bg+correction));pressurechecks+=1
            assert not correction.zero
    assert saved['accepted_full_spatial_O2_alternative_remains_separate_global_endpoint_only']
    assert saved['whole_O2_spatial_IBP_bound_not_sliced_as_partial_cell_increment']
    assert not saved['full24_original_C1_integral_range_transport_enclosed']
    return dict(passed=True,actual_five_interval_live_operator_replays=5,
        serial_composite_actual_correction_C0_Z_overlaps=10,accepted_whole_O2_direct_and_background_overlaps=overlap_count,
        original_background_once_separate_P0_nonzero_pressure_checks=pressurechecks,
        upstream_and_density_queries_not_replayed=True,spatial_alternative_kept_global_not_sliced=True,
        original24_control_bridge_not_prematurely_admitted=True)


def guards(owner):
    count=0
    for callback in (lambda:owner.integrate(N=True),lambda:owner.integrate(N=7),lambda:owner.integrate(N=2048),
        lambda:current.exact_pieces(1,0),lambda:current.exact_pieces(-1,0),lambda:current.exact_pieces(0,2)):
        try:callback()
        except ValueError:count+=1
        else:raise AssertionError('Wrong source frequency or piece domain admitted')
    return dict(passed=True,candidate_and_exact_measure_domain_guards=count)


@current.native.inlet.source_precision
def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['candidate_N']==1024
    flags=('full24_original_C1_integral_range_transport_enclosed','actual_five_controls_installed',
        'functional_terminal_identity_solved','current_whole_N_selected',*current.packets.OPEN)
    assert all(manifest[flag] is False for flag in flags)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    bridge,_=current.native.inlet.native_bridge_owner();checks={}
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.OriginalO2FiveSlopePartitions(bridge);saved=manifest['actual_original_O2_five_slope_partitions']
        with mp.workdps(owner.ctx.dps+40):
            for label,callback in (('independent_exact_piece_measures_original_density_and_masses',lambda:independent_exact_partition_and_weights(owner,saved)),
                ('actual_five_source_operator_memory_and_whole_O2_overlap',lambda:actual_five_source_operator_replay(owner,saved)),
                ('candidate_and_exact_domain_guards',lambda:guards(owner))):
                checks[label]=callback();print(label+' PASS',flush=True)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        actual_five_original_O2_slope_partition_operators_installed=True,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Exact five original O2 intervals with68 disjoint measure pieces use accepted actual N1024 full-Z predicate density covers, independent true masses/phase origin, actual upstream and background/P0 memory. Whole-O2 source overlaps. Original24 adapter/control/terminal/global N/recursion/full NS remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(receipt),indent=2).encode()+b'\n')
    print('Actual full-predicate original five O2 slope partitions PASS',flush=True);return receipt


if __name__=='__main__':run()
