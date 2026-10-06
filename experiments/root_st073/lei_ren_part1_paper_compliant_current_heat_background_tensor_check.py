"""Admit actual current collar/heat tensors, two joins and exact exterior."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_heat_background_tensor import (
    CurrentHeatBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,TENSOR_KEYS,
    sha,pack,encode,source_precision,endpoints)
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes


def physical_rows(view):
    def flatten(value):
        if isinstance(value,dict) and 'signed_coefficient' in value:return [value]
        result=[]
        for v in (value.values() if isinstance(value,dict) else value):result.extend(flatten(v))
        return result
    return sum((flatten(view[key]) for key in TENSOR_KEYS),[])


def check_physical(view,exact_zero=False):
    if any(view[k] for k in OPEN):raise ValueError('Current regional heat scope exceeded')
    rows=physical_rows(view)
    if len(rows)!=65:raise ValueError('Actual tensor/divergence/remainder layout incomplete')
    nonzero=0
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Source zero has nonzero coefficient/bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Current exact scale logs lost')
    if exact_zero and nonzero:raise ValueError('Full Gamma physical tensor/remainder must be source-exact zero')
    return len(rows),nonzero


def check_view(view):
    if not view['actual_full_stress_not_local_difference'] or not view['actual_Cp_Dtheta_zero_from_checked_same_current_functions']:raise ValueError('Actual full closed-current heat source missing')
    for name in ('A','E','P','K'):
        rows=view['current_actual_normalized_full_moment_rows'][name]
        if len(rows)!=5 or any(v.order!=5 for v in rows):raise ValueError('Full moment mixed4/axial5 rows required')
        for row in rows:
            for value in row.coefficients:finite(value)
    for value in view['normalized_energy_native_consistency'].coefficients:
        lo,hi=endpoints(value)
        if not lo<=0<=hi:raise ValueError('Current native complete energy contradicts same remaining source')
    if len(view['stable_actual_absolute_pressure_mixed4_factored'])!=15:raise ValueError('Actual pressure mixed4 layout missing')
    for value in view['stable_actual_absolute_pressure_mixed4_factored'].values():finite(value)
    for name,grid in view['current_actual_normalized_stress_mixed3'].items():
        if len(grid)!=10:raise ValueError('Actual stress mixed3 layout missing')
        for value in grid.values():finite(value)
        if view['source_exact_Gamma_tensor_and_remainder_zero'] and name in ('theta','axial') and any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in grid.values()):raise ValueError('Full Gamma stress source zero lost')
    return check_physical(view,view['source_exact_Gamma_tensor_and_remainder_zero'])


def check_interface(view):
    rows=view['common_actual_tensor_rows']
    if len(rows)!=65 or view['current_common_tensor_contribution_count']!=65 or not view['overlap_not_used_for_tensor_source_equality'] or not view['actual_completed_tensor_and_divergence_remainder_common_source_bound'] or any(view[k] for k in OPEN):raise ValueError('Current heat tensor join incomplete')
    nonzero=0
    for row in rows.values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Zero interface has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    if view['source_exact_Gamma_zero_at_interface'] and nonzero:raise ValueError('Collar/exterior completed tensor join must be source zero')
    return len(rows),nonzero


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentHeatBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_heat_tensor_views','current_actual_two_heat_tensor_interface_bounds','current_unbounded_exact_Gamma_physical_identity'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current heat tensor manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_heat_tensor_views'])!=set(VIEWS):raise ValueError('Producer acceptance or view scope differs')
    count=nonzero=0
    for name,args in VIEWS.items():
        view=field.chart(*args)
        if encode(pack(view))!=raw['current_actual_heat_tensor_views'][name]:raise ValueError('Actual current heat tensor replay differs: '+name)
        n,m=check_view(view);count+=n;nonzero+=m
    joins={};fresh={}
    for seam in SEAMS:
        view=field.interface(seam)
        if encode(pack(view))!=raw['current_actual_two_heat_tensor_interface_bounds'][seam]:raise ValueError('Current actual heat common bounds differ')
        n,m=check_interface(view);joins[seam]=dict(rows=n,nonzero=m)
        view=field.interface(seam,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
        n,m=check_interface(view);fresh[seam]=dict(rows=n,nonzero=m,bounds=view)
    unbounded=field.unbounded_exterior()
    if encode(pack(unbounded))!=raw['current_unbounded_exact_Gamma_physical_identity']:raise ValueError('Unbounded current exterior function theorem differs')
    un,m=check_physical(unbounded,True)
    if not unbounded['original_unbounded_exterior_covered'] or not unbounded['actual_stress_tensor_divergence_remainder_and_momentum_exactly_zero']:raise ValueError('Full Gamma exterior function scope missing')
    rejected=[]
    for label,kwargs in (('foreign_chart',dict(chart='waiting')),('outside_Z',dict(Z=2)),('outside_collar',dict(t=4)),('nonfinite_offset',dict(t='inf')),('nonfinite_time',dict(log_tau='-inf')),('nonpositive_viscosity',dict(viscosity=0))):
        args=dict(chart='heat_collar',Z='.2',t='.4',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.chart(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid heat tensor domain accepted')
    clone=copy.copy(field);clone.heat=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_heat_owner')
    else:raise ValueError('Foreign current heat owner accepted')
    clone=copy.copy(field);clone.exterior=copy.copy(field.exterior);clone.exterior.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_closed_exterior')
    else:raise ValueError('Unchecked closed-current exterior accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_heat_tensor_view_count=len(VIEWS),current_actual_heat_physical_tensor_divergence_remainder_rows_checked=count,
        current_actual_heat_nonzero_physical_contributions_checked=nonzero,current_actual_two_heat_tensor_interface_rows=joins,
        current_actual_fresh_sector_two_heat_tensor_interface_rows=fresh,current_unbounded_exact_Gamma_zero_tensor_rows_checked=un,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=6,
        current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),original_full_collar_normalization_identity_count=45,
        current_two_heat_primitive_seam_identity_count=120,
        source_theorems={k:raw[k] for k in raw if k.endswith('_theorem') or k.endswith('_theorems') or k.startswith('checked_')},
        current_nonzero_complete_moments_energy_pressure_and_native_diagnostics_retained=True,
        original_small_defect_stress_cancellation_before_KR_scaling=True,source_exact_exterior_zero_from_current_five_histories_not_interval_overlap=True,
        invalid_domains_and_foreign_current_owner_rejected=rejected,
        scope='Actual current collar and full Gamma exterior physical tensors, two function joins, and exact regional exterior NS/remainder identity. Seven current tensor regions and six adjacent tensor joins; remaining charts, axis/full33, global cone/lift/NS/independent temporal remainder/energy/points and n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual collar/exterior tensor and two joins PASS; full-unbounded exterior regional NS identity; global/time open',flush=True)
    return result


if __name__=='__main__':run()
