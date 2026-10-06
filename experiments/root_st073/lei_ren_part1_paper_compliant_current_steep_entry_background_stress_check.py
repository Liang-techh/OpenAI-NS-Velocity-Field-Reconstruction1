"""Check current entry full tensor and current angular-entry tensor join."""
import copy
import json
from pathlib import Path

from lei_ren_part1_paper_compliant_current_steep_entry_background_stress import (
    CurrentSteepEntryBackgroundStress,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,sha,pack,encode)
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import check_view,finite
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def check_interface(view):
    if view['current_common_tensor_contribution_count']!=65 or not view['actual_completed_tensor_and_divergence_remainder_common_source_bound']:
        raise ValueError('Actual completed tensor interface contribution layout incomplete')
    if not view['overlap_not_used_for_tensor_source_equality'] or any(view[k] for k in OPEN):
        raise ValueError('Actual tensor interface function proof or scope missing')
    rows=view['common_actual_angular_entry_tensor_rows']
    if len(rows)!=65:raise ValueError('All actual tensor/divergence/completion/remainder interface bounds required')
    nonzero=0
    for row in rows.values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Exact-zero common tensor row has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    if not nonzero:raise ValueError('Actual common tensor was replaced by the zero local support difference')
    return len(rows),nonzero


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentSteepEntryBackgroundStress(require_checked=False)
    field.assert_graph()
    omitted={'current_actual_steep_entry_tensor_views','current_actual_angular_entry_tensor_interface_bounds'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:
        raise ValueError('Actual current entry full moment/tensor/source manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Producer claims its own acceptance')
    if set(raw['current_actual_steep_entry_tensor_views'])!=set(VIEWS):raise ValueError('Current entry domain/boundary/fresh views incomplete')
    count=nonzero=0
    for name,args in VIEWS.items():
        view=field.entry(*args)
        if encode(pack(view))!=raw['current_actual_steep_entry_tensor_views'][name]:
            raise ValueError('Current actual entry tensor replay differs: '+name)
        n,m=check_view(field,view);count+=n;nonzero+=m
    interface=field.interface()
    if encode(pack(interface))!=raw['current_actual_angular_entry_tensor_interface_bounds']:
        raise ValueError('Actual current angular-entry common tensor bounds differ')
    n,m=check_interface(interface)
    fresh=field.interface(Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
    fn,fm=check_interface(fresh)
    rejected=[]
    for label,kwargs in (('outside_Z',dict(Z=2)),('outside_entry',dict(t=-1)),
            ('nonfinite_time',dict(log_tau='-inf')),('nonpositive_viscosity',dict(viscosity=0)),
            ('nonfinite_angle',dict(theta='inf'))):
        args=dict(Z='.2',t='.4',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.entry(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid current physical entry request accepted')
    clone=copy.copy(field);clone.steep=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('different_current_steep_owner')
    else:raise ValueError('Foreign current steep owner accepted')
    clone=copy.copy(field);clone.angular=copy.copy(field.angular);clone.angular.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_actual_angular_owner')
    else:raise ValueError('Unchecked actual angular tensor accepted')
    theorem_keys=('current_entry_transport_source_theorem','current_entry_full_moment_normalization_AST_theorem',
        'current_entry_shape_and_remaining_pressure_source_theorem','current_entry_full_A_KX_inertial_correlation_AST_theorem',
        'original_full_stress_baseline_cancellation_and_KR_units','current_angular_entry_actual_tensor_function_join_theorem',
        'original_entry_physical_tensor_remainder_theorem')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_entry_tensor_view_count=len(VIEWS),
        current_actual_entry_physical_tensor_divergence_remainder_rows_checked=count,
        current_actual_entry_nonzero_physical_contributions_checked=nonzero,
        actual_current_completed_tensor_adjacent_interface_count=1,
        actual_angular_entry_common_tensor_rows_checked=n,actual_angular_entry_nonzero_common_tensor_rows=m,
        fresh_sector_angular_entry_common_tensor_rows_checked=fn,fresh_sector_nonzero_common_tensor_rows=fm,
        fresh_sector_angular_entry_tensor_bounds=fresh,
        current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
        original_full_entry_moment_normalization_identity_count=len(field.moment_proof['original_full_entry_defect_AST_normalized_mixed4_identities']),
        full_entry_A_KX_and_inertial_correlation_identity_count=len(field.correlation['actual_full_A_KX_and_original_correlated_entry_stress_mixed3_identities']),
        current_angular_entry_primitive_function_trace_identity_count=60,
        source_theorems={k:raw[k] for k in theorem_keys},invalid_domains_and_current_owners_rejected=rejected,
        actual_full_current_histories_absolute_pressure_complete_future_retained=True,
        current_source_function_tensor_join_not_interval_overlap=True,
        direct_remaining_entry_pressure_is_same_original_integral_by_additivity=True,
        normalized_native_energy_consistency_not_used_as_functional_proof=True,
        source_domain=raw['source_domain'],
        scope='Actual current steep-entry full stress mixed3, completed regional physical tensor, divergence/leading axial-viscosity remainder mixed2 and one actual current angular-entry tensor function join with 65 common factored bounds, plus a fresh compact sector. Other actual chart tensors/joins, axis/full33 interfaces, global cone/lift/NS/temporal remainder/energy/resolved points and n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual steep-entry tensor and angular-entry completed tensor join PASS; global/temporal gates remain open',flush=True)
    return result


if __name__=='__main__':run()
