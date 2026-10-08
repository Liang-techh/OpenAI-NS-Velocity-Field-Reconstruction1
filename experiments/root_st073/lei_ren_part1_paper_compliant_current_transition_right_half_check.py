"""Focused true flat-collar and explicit unresolved source-atlas audit.

Checks saved genuine source, original q/ordinary jets and zero own integrals;
the original nonzero incoming history remains an explicit affine argument.
No native owner or accepted earlier numerical integrals are rerun.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_transition_right_half as source
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals_check as generic

HERE,sha,ep,iv=source.HERE,source.sha,generic.ep,generic.iv
rc,common=source.rc,source.common;same=generic.same


def archive_payload(archive):
    compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
    assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
    assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    return json.loads(raw)


def atlas_diagnostic_check(c,coordinates,Z,payload):
    assert payload['source_family']==coordinates.family and payload['candidate_N']==source.N
    assert ep(c.mpf(payload['exact_Z_range']))==Z
    assert payload['full_source_partition_preserved_even_when_cutoff_unresolved']
    assert not payload['numerical_complete_tail_or_Rc_targets_admitted']
    queries=0;windows=0;reasons={};charts={}
    for chart,got in payload['original_source_charts'].items():
        summary=got['summary'];rows=got['ordered_original_source_partition'];parents=got['original_refused_parent_source_queries']
        assert summary['chart']==chart and summary['exact_original_native_domain']==['0','1']
        assert list(map(Fraction,summary['original_seed_partition_edges']))==source.atlas_source.original_edges(chart)
        cursor=Fraction(0);covered=0;unresolved=0
        for row in rows:
            left,right=map(Fraction,row['exact_native_endpoints']);assert left==cursor<right;cursor=right
            assert row['chart']==chart and row['original_endpoint_focused_partition']
            assert 0<=row['source_partition_depth']<=2
            parsed,_=generic.cell_check(c,coordinates,row,Z);queries+=1;windows+=1
            if parsed is None:
                unresolved+=1
                q=row['actual_original_phase_inverse_and_density_source']['original_q_slow_jet_source']
                reasons[q['status']]=reasons.get(q['status'],0)+1
                assert row['status']=='requires_original_source_phase_subdivision'
            else:covered+=1
        assert cursor==1 and covered==summary['enclosed_source_cells'] and unresolved==summary['unresolved_source_cells']
        assert summary['accepted_partition_cells']==len(rows)
        assert summary['exact_unresolved_native_windows']==[r['exact_native_endpoints'] for r in rows if r['status']!='enclosed']
        assert not summary['complete_chart_numerical_integral_admitted']
        for row in parents:
            assert row['chart']==chart and row['source_partition_depth']<2
            parsed,_=generic.cell_check(c,coordinates,row,Z);assert parsed is None;queries+=1
            left,right=map(Fraction,row['exact_native_endpoints']);mid=(left+right)/2
            assert any(Fraction(r['exact_native_endpoints'][1])==mid for r in rows)
        assert summary['actual_source_queries']==len(rows)+len(parents)
        assert len(parents)==summary['refined_refused_parent_cells']<=12
        charts[chart]=dict(enclosed_source_cells=covered,explicit_unresolved_source_cells=unresolved)
    return dict(exact_Z_range=payload['exact_Z_range'],source_queries=queries,full_partition_windows=windows,
        charts=charts,actual_unresolved_q_status_counts=reasons,
        source_refusals_preserved_as_refusals_not_zero_or_integral=True)


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert not manifest['actual_numerical_entire_transition_or_tail_or_Rc_histories_admitted']
    assert all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    atlas=json.loads((HERE/source.atlas_source.NAME).read_bytes())
    assert atlas[source.atlas_source.GATE] and atlas['source_family']==manifest['source_family']
    assert not atlas['full_original_axial_transition_function_phase_coverage_certified']
    c=MPIntervalContext();c.dps=240;checks=[];diagnostics=[];phase_count=0
    with mp.workdps(300):
        for archive in manifest['actual_original_transition_right_half_source_archives']:
            payload=archive_payload(archive);Z=ep(c.mpf(payload['exact_Z_range']))
            theorem=payload['common_directed_coordinate_theorem'];bases=theorem['common_log_bases']
            coordinates=rc.history.CommonSourceCoordinates(c,iv(c,bases[1]),manifest['source_family'])
            assert all(ep(iv(c,bases[k]))==(0,0) for k in (0,2,3,4))
            assert payload['source_family']==coordinates.family and payload['candidate_N']==source.N
            assert payload['exact_original_native_window']==['1/2','1'] and payload['exact_true_log_radius_length']=='1/2'
            assert payload['unresolved_transition_left_half_not_skipped_or_assumed_zero']
            assert payload['original_incoming_is_required_argument_and_rate0_memory_is_not_reset']
            assert not payload['actual_numerical_entire_transition_or_Rc_histories_admitted']
            row=payload['original_entire_transition_right_half_source_query']
            assert row['chart']=='O3_slope_mu' and row['exact_native_endpoints']==['1/2','1']
            got,n=generic.cell_check(c,coordinates,row,Z);assert got is not None;phase_count+=n
            q=row['actual_original_phase_inverse_and_density_source']['original_q_slow_jet_source']
            assert q['status']=='enclosed' and q['branch']=='flat' and not q['active_body_enclosure_evaluated']
            assert q['current_refined_original_q_C0']['exact_zero']
            assert all(v['exact_zero'] for v in q['original_q_ordinary_slow_derivative_enclosures'].values())
            assert q['original_correlated_shear_and_q']['original_q_enclosure']['branch']=='flat'
            difference=q['original_correlated_shear_and_q']['original_q_enclosure']['branch_difference']
            assert ep(iv(c,difference['coefficient_interval']))[0]>=0
            assert ep(got['geometry']['regular'])==(mp.mpf('.5'),mp.mpf('.5'))
            assert ep(got['geometry']['scalar_cover'])==(mp.mpf('.5'),mp.mpf('.5'))
            assert all(v.zero for v in (*got['values'].values(),*got['Z_derivatives'].values()))
            operator=payload['original_right_half_C1_integral_operator']
            assert operator['steps']==1 and operator['incoming_argument_not_assumed_or_reset']
            assert operator['quiet_cells_preserve_original_rate0_pressure_memory']
            assert operator['original_rates']=={k:str(r) for k,r in rc.RATES.items()}
            for k in rc.RATES:
                same(c,operator['incoming_C0_Z_decay_coefficients'][k],got['factors'][k]['decay'])
                same(c,operator['cumulative_signed_increment_C0_enclosures'][k],got['values'][k])
                same(c,operator['cumulative_signed_increment_Z_enclosures'][k],got['Z_derivatives'][k])
            pressure=got['factors']['p']['decay'];assert ep(pressure.coefficient)==(1,1) and pressure.scale.powers==(0,0,0,0,0)
            diagnostic_archive=next(a for a in atlas['actual_original_endpoint_focused_source_atlas_archives']
                if ep(c.mpf(a['exact_Z_range']))==Z)
            diagnostics.append(atlas_diagnostic_check(c,coordinates,Z,archive_payload(diagnostic_archive)))
            checks.append(dict(exact_Z_range=payload['exact_Z_range'],original_native_window=['1/2','1'],
                all_original_q_and_ordinary_q_jets_exact_zero=True,all_ten_own_integral_C0_Z_increments_exact_zero=True,
                original_incoming_m_h_k_e_decayed_not_reset=True,original_pressure_incoming_C0_Z_exactly_retained=True))
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        actual_original_transition_right_half_source_checks=checks,actual_flat_phase_density_cells_checked=phase_count,
        actual_exact_zero_own_integral_C0_Z_rows_checked=20,executed_endpoint_atlas_diagnostics=diagnostics,
        source_atlas_diagnostics_not_admitted_as_numerical_integrals=True,
        whole_original_transition_or_tail_targets_or_physical_closure_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original transition right-half C1 PASS:20 exact-zero own rows; incoming memory retained; unresolved atlas explicit',flush=True)
    return result


if __name__=='__main__':run()
