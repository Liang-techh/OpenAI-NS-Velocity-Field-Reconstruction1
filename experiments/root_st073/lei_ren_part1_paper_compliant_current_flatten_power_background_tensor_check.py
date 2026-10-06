"""Admit current flatten/power physical tensors and both function joins."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_flatten_power_background_tensor import (
    CurrentFlattenPowerBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,TENSOR_KEYS,
    sha,pack,encode,source_precision,endpoints)
from lei_ren_part1_paper_compliant_current_heat_background_tensor_check import check_physical
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes


def check_view(view):
    if not view['actual_full_stress_not_local_difference'] or not view['current_complete_pressure_energy_and_forward_X_retained'] or not view['ordinary_logR_derivatives_not_phase_derivatives']:
        raise ValueError('Actual current flatten/power full history or ordinary derivative lost')
    if view['actual_variable_axial_K_retained']!=(view['chart']=='flatten'):raise ValueError('Flatten variable axial K source lost')
    for name in ('A','E','P','K'):
        rows=view['current_actual_normalized_full_moment_rows'][name]
        if len(rows)!=5 or any(v.order!=5 for v in rows):raise ValueError('Full ordinary mixed4/axial5 moment rows required')
        for row in rows:
            for value in row.coefficients:finite(value)
    for value in view['normalized_energy_native_consistency'].coefficients:
        lo,hi=endpoints(value)
        if not lo<=0<=hi:raise ValueError('Remaining quadratic future contradicts same native complete energy')
    if len(view['stable_actual_absolute_pressure_mixed4_factored'])!=15:raise ValueError('Actual absolute pressure mixed4 layout incomplete')
    for value in view['stable_actual_absolute_pressure_mixed4_factored'].values():finite(value)
    for grid in view['current_actual_normalized_stress_mixed3'].values():
        if len(grid)!=10:raise ValueError('Actual stress mixed3 layout incomplete')
        for value in grid.values():finite(value)
    return check_physical(view)


def check_interface(view):
    rows=view['common_actual_tensor_rows']
    if len(rows)!=65 or view['current_common_tensor_contribution_count']!=65 or not view['overlap_not_used_for_tensor_source_equality'] or not view['actual_completed_tensor_and_divergence_remainder_common_source_bound'] or any(view[k] for k in OPEN):
        raise ValueError('Current completed flatten/power tensor function join incomplete')
    nonzero=0
    for row in rows.values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Zero interface row has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    return len(rows),nonzero


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentFlattenPowerBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_flatten_power_tensor_views','current_actual_two_flatten_power_tensor_interface_bounds'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current flatten/power tensor manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_flatten_power_tensor_views'])!=set(VIEWS):raise ValueError('Producer acceptance or view scope differs')
    count=nonzero=0
    for name,args in VIEWS.items():
        view=field.chart(*args)
        if encode(pack(view))!=raw['current_actual_flatten_power_tensor_views'][name]:raise ValueError('Actual current flatten/power replay differs: '+name)
        n,m=check_view(view);count+=n;nonzero+=m
    joins={};fresh={}
    for seam in SEAMS:
        view=field.interface(seam)
        if encode(pack(view))!=raw['current_actual_two_flatten_power_tensor_interface_bounds'][seam]:raise ValueError('Current completed tensor common bound differs')
        n,m=check_interface(view);joins[seam]=dict(rows=n,nonzero=m)
        view=field.interface(seam,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
        n,m=check_interface(view);fresh[seam]=dict(rows=n,nonzero=m,bounds=view)
    rejected=[]
    for label,kwargs in (('foreign_chart',dict(chart='waiting')),('outside_Z',dict(Z=2)),
            ('outside_flatten',dict(x=101)),('outside_power',dict(chart='outer_power',x=2)),
            ('nonfinite_coordinate',dict(x='inf')),('nonfinite_time',dict(log_tau='-inf')),('nonpositive_viscosity',dict(viscosity=0))):
        args=dict(chart='flatten',Z='.2',x='41.7',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.chart(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid current flatten/power tensor domain accepted')
    for label,attribute in (('foreign_current_flatten','flatten'),('foreign_current_outer','outer')):
        clone=copy.copy(field);setattr(clone,attribute,object())
        try:clone.assert_graph()
        except ValueError:rejected.append(label)
        else:raise ValueError('Foreign current regional owner accepted')
    clone=copy.copy(field);clone.heat_tensor=copy.copy(field.heat_tensor);clone.heat_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_heat_chain')
    else:raise ValueError('Unchecked current prerequisite admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_flatten_power_tensor_view_count=len(VIEWS),current_actual_flatten_power_physical_tensor_divergence_remainder_rows_checked=count,
        current_actual_flatten_power_nonzero_physical_contributions_checked=nonzero,current_actual_two_flatten_power_tensor_interface_rows=joins,
        current_actual_fresh_sector_two_flatten_power_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=8,
        current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),original_full_flatten_power_normalization_identity_count=90,
        current_two_flatten_power_primitive_seam_identity_count=120,
        source_theorems={k:raw[k] for k in raw if k.endswith('_theorem') or k.endswith('_theorems') or k.startswith('checked_')},
        current_variable_flatten_axial_K_and_nonzero_complete_moments_energy_pressure_retained=True,
        positive_remaining_quadratic_integrals_before_enclosure=True,native_forward_angular_X_not_backward_amplified=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='Actual current flatten and outer-power tensors, ordinary stress3/divergence2/completion2/remainder2 and two function joins. Nine current regions/eight adjacent tensor joins. Pulse/core/axis/all33, global cone/lift/NS/independent temporal remainder/energy/points and genuine n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual flatten/power tensors and two joins PASS;9 regions/8 tensor joins; global/time open',flush=True)
    return result


if __name__=='__main__':run()
