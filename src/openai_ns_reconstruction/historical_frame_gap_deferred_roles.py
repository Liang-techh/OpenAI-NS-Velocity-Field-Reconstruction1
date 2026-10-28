"""Recover saved historical pulse-gap source roles as deferred substitutions.

Only authenticated archived JSON fields are read; no graph nodes are
traversed, no old constructor is called, and no source value is evaluated.
The unexecuted substitution identifies candidate inherited F roles; it does
not prove the adapter's global pressure decomposes into these local pieces.
"""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path

PREFIX='lei_ren_part1_paper_compliant_current_original_'
RP=PREFIX+'Rp_native_C5_pulse_source.json.gz'
C3=PREFIX+'n1_pulse_pressure_C3_functions.json.gz'
OWN=PREFIX+'n1_pulse_own_fields.json.gz'
RP_CHECK=PREFIX+'Rp_native_C5_pulse_source_check.json'
C3_CHECK=PREFIX+'n1_pulse_pressure_C3_functions_check.json'
OWN_CHECK=PREFIX+'n1_pulse_own_fields_check.json'
COMPARE=PREFIX+'native_Rp_m15_source_comparison.json'
COMPARE_HASH='c38540a019963d2023ba9fc7d4e9c292fd7af4d23fe103018d3f5009687a6586'
HASHES={RP:'442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
 C3:'bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d',
 OWN:'2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b'}


def _load_archive(path: Path, expected: str):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('archive raw hash mismatch: '+path.name)
    return json.loads(gzip.decompress(raw))


def _read_json(path: Path):
    return json.loads(path.read_bytes())


def recover(root: Path):
    folder=root/'experiments/root_st073'
    checks={name:_read_json(folder/name) for name in (RP_CHECK,C3_CHECK,OWN_CHECK)}
    if not all(v.get('all_passed') is True for v in checks.values()):
        raise ValueError('all three archived receipts must have passed')
    for check in (checks[RP_CHECK],checks[C3_CHECK]):
        for name in (C3,OWN):
            if check['input_hashes'].get(name)!=HASHES[name]:
                raise ValueError('archived dependency pin mismatch: '+name)
    if checks[RP_CHECK]['input_hashes'].get(RP)!=HASHES[RP]:
        raise ValueError('Rp archived receipt pin mismatch')
    rp=_load_archive(folder/RP,HASHES[RP])
    c3=_load_archive(folder/C3,HASHES[C3])
    own=_load_archive(folder/OWN,HASHES[OWN])
    if not (rp['source_family']==c3['source_family']==own['source_family']):
        raise ValueError('archive source-family metadata differs')
    frame=rp['Rp_native_C5_pulse_source_bindings']
    compare_path=folder/COMPARE
    if hashlib.sha256(compare_path.read_bytes()).hexdigest()!=COMPARE_HASH:
        raise ValueError('frame pressure comparison report hash mismatch')
    compare=_read_json(compare_path)
    global_frame=compare['component_result']['source_interface']['native_complete_pulse_pressure']['m15']
    if frame['pressure']['m15']!=global_frame:
        raise ValueError('adapter global pressure differs from saved frame role')
    nested=frame['inputs'][C3]
    original=c3['C5_pulse_C3_pressure_bindings']
    # Direct saved-field equality only. This deliberately does not follow or
    # compare graph nodes, evaluate ids, or reconstruct an endpoint mapping.
    fields=('native','leading_fields','endpoints','pieces','pressure')
    for key in fields:
        if nested[key]!=original[key]:
            raise ValueError('saved C3 binding disagrees with Rp input: '+key)
    gap=own['own_main_gap_end_fields']['gap']
    endpoints=original['endpoints']['gap']
    if not isinstance(endpoints,list) or len(endpoints)!=2:
        raise ValueError('expected two historical gap endpoints')
    substitution=frame['input_substitution']
    if not isinstance(substitution,dict) or not substitution:
        raise ValueError('missing saved inlet substitution')
    selected={k:endpoints[0][k] for k in ('R_q_y','q_plus_q_y')}
    selected_upper={k:endpoints[1][k] for k in ('R_q_y','q_plus_q_y')}
    return {
      'schema':'historical_frame_gap_deferred_roles_v1',
      'status':'deferred_source_roles_not_numeric',
      'archive_raw_hashes':HASHES,'archived_receipts_passed':True,
      'same_source_family_metadata':True,
      'nested_C3_saved_fields_equal_source_archive':list(fields),
      'archive_input_substitution':substitution,
      'original_own_gap_field_template':{'archive':OWN,'pointer':'$.own_main_gap_end_fields.gap',
         'query':gap['query'],'radius':gap['radius'],'log_radius':gap['log_radius'],
         'domain':gap['domain'],'ordinary_Z_rows_are_source_refs':True},
      'original_gap_endpoint_roles':{'archive':C3,
         'pointer':'$.C5_pulse_C3_pressure_bindings.endpoints.gap',
         'lower':selected,'upper':selected_upper},
      'original_gap_m15_rows':original['pieces']['gap']['m15'],
      'original_global_m15_rows':original['pressure']['m15'],
      'adapter_global_m15_rows':frame['pressure']['m15'],
      'global_frame_role_matches_saved_comparison':True,
      'deferred_composition':{'operation':'simultaneous_saved_inlet_substitution',
         'input':'original gap radius, gap endpoints and gap m15 roles',
         'mapping_pointer':'archive_input_substitution','executed':False},
      'adapter_rewrites_global_pressure_not_gap_endpoints':True,
      'local_gap_decomposition_of_adapter_global_proved':False,
      'deferred_F_roles_equal_actual_rewritten_nodes_proved':False,
      'G_equals_F_proved':False,'P_equals_F_proved':False,
      'coordinates_and_units_between_sources_proved':False,
      'ordinary_Z1_transport_between_sources_proved':False,
      'source_values_computed':False,'source_graph_nodes_traversed':False,
      'old_operators_replayed':False,'ancestor_quadrature_calls':0,
      'old_error_terms_removed':False,'extra_ledger_terms_added':False,
      'independent_reference_offset_preserved':True,
      'actual_physical_correction_roots_preserved':[290917,290919,290921,290923],
      'admissible_for_numeric_use':False,'full_goal_complete':False}
