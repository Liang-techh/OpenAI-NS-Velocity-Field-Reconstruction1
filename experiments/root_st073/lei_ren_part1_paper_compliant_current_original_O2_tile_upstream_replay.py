"""Reuse accepted genuine N1024 O2 integrals with refined actual tile inlets."""
from fractions import Fraction
import lei_ren_part1_paper_compliant_current_native_upstream_tile_uniform_C1 as upstream

common=upstream.common;ep=common.ep;iv=common.interval;KEYS=upstream.KEYS
HERE,PREFIX,sha=upstream.HERE,upstream.PREFIX,upstream.sha
NAME=PREFIX+'current_original_O2_tile_upstream_replay.json'
GATE='original_O2_N1024_integrals_replayed_with_genuine_refined_upstream_tiles'


def inputs(c,payload):
    family=payload['source_family'];route=payload['complete_actual_pre_O2_tile_route']
    if route['source_family']!=family or route['candidate_N']!=upstream.N:
        raise ValueError('Same genuine upstream family and candidate N required')
    if not route['no_original_interval_skipped_between_true_inlet_and_O2_inlet']:
        raise ValueError('Complete actual original upstream route required')
    background=route['actual_O2_inlet_original_background_and_separate_P0_Z']
    if background['original_units']!=common.previous.five.UNITS or not background['P0_not_merged_into_pressure_history']:
        raise ValueError('Original five units and separate pressure datum required')
    bases=background['common_directed_coordinate_theorem']['common_log_bases']
    if any(ep(iv(c,bases[i]))!=(0,0) for i in (0,2,3,4)):
        raise ValueError('Same original common Pstar-squared basis required')
    coordinates=common.middle.history.CommonSourceCoordinates(c,iv(c,bases[1]),family)
    values={k:common.restore_common_source(v,coordinates) for k,v in route['actual_original_inlet_to_O2_inlet_correction_C0'].items()}
    jets={k:common.restore_common_source(v,coordinates) for k,v in route['actual_original_inlet_to_O2_inlet_correction_Z'].items()}
    if set(values)!=set(KEYS) or set(jets)!=set(KEYS):raise ValueError('All five actual upstream C0/Z functions required')
    P0=common.restore_common_source(background['original_separate_P0_over_Pstar_squared'],coordinates)
    P0_Z=common.restore_common_source(background['original_separate_P0_Z_over_Pstar_squared'],coordinates)
    return coordinates,values,jets,P0,P0_Z


def replay(c,payload,record,archive):
    family=payload['source_family']
    common.current.require_identity(record,family,family['datum_enclosure_sha256'])
    if record['explicit_candidate_N']!=upstream.N or payload['candidate_N']!=upstream.N:
        raise ValueError('N7 cannot enter the actual N1024 upstream route')
    if tuple(Fraction(v) for v in record['exact_Z_range'])!=tuple(Fraction(v) for v in payload['exact_Z_range']):
        raise ValueError('Same entire original axial tile required')
    coordinates,values,jets,P0,P0_Z=inputs(c,payload)
    scalar=lambda v:coordinates.scalar(iv(c,v))
    original={k:scalar(v) for k,v in record['source_defined_original_histories_at_y1'].items()}
    originalZ={k:scalar(v) for k,v in record['source_defined_original_history_Z_at_y1'].items()}
    modulation={k:scalar(v) for k,v in record['five_original_C0_integral_contributions'].items()}
    modulationZ={k:scalar(v) for k,v in record['five_genuine_ordinary_Z_integral_contributions'].items()}
    correction={};correctionZ={};own={};ownZ={}
    for k,rate in common.previous.five.RATES.items():
        decay=coordinates.decay(1,rate)
        correction[k]=decay*values[k]+modulation[k]
        correctionZ[k]=decay*jets[k]+modulationZ[k]
        own[k]=original[k]+correction[k];ownZ[k]=originalZ[k]+correctionZ[k]
    return dict(source_family=family,candidate_N=upstream.N,exact_Z_range=record['exact_Z_range'],
        ordered_source_cells=record['ordered_source_cells'],exact_y_window=['0','1'],
        original_P0_datum_sha256=family['datum_enclosure_sha256'],normalized_own_units=common.previous.five.UNITS,
        actual_upstream_binding=dict(manifest=upstream.NAME,manifest_sha256=sha(upstream.NAME),
            complete_tile_source_archive=archive['filename'],archive_sha256=sha(archive['filename']),
            exact_incoming_source_function_graph=payload['exact_upstream_source_function_incoming_binding'],
            source_is_actual_correction_not_original_plus_correction=True,
            actual_tile_source_requeried_not_whole_Z_relabelled=True),
        accepted_original_O2_integral_binding=dict(manifest=common.NAME,manifest_sha256=sha(common.NAME),
            receipt=common.RECEIPT,receipt_sha256=sha(common.RECEIPT),candidate_N=upstream.N,
            exact_integral_record_Z_range=record['exact_Z_range'],ordered_source_cells=record['ordered_source_cells'],
            original_integral_cells_and_true_phase_evidence_reused_without_recomputation=True),
        common_coordinates=coordinates.record(),
        actual_refined_upstream_C0={k:v.record() for k,v in values.items()},
        actual_refined_upstream_Z={k:v.record() for k,v in jets.items()},
        propagated_actual_correction_C0={k:v.record() for k,v in correction.items()},
        propagated_actual_correction_Z={k:v.record() for k,v in correctionZ.items()},
        actual_original_plus_propagated_correction_C0={k:v.record() for k,v in own.items()},
        actual_original_plus_propagated_correction_Z={k:v.record() for k,v in ownZ.items()},
        separate_P0=P0.record(),separate_P0_Z=P0_Z.record(),
        actual_absolute_pressure=(P0+own['p']).record(),actual_absolute_pressure_Z=(P0_Z+ownZ['p']).record(),
        original_background_added_once_and_rate_zero_pressure_memory_retained=True,
        no_O2_integral_rerun_no_incoming_zero_fallback_no_enormous_scalar_cast=True,
        full_Rc_terminal_identity_or_global_N_or_stress_or_recursion_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False))
