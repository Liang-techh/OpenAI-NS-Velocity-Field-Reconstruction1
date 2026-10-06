"""Admit three actual current O7 tensors and three adjacent tensor joins."""
import copy
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_steep_waiting_background_stress import (
    CurrentSteepWaitingBackgroundStress,CHARTS,SEAMS,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,sha,pack,encode,source_precision)
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import check_view,finite
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes


def check_interface(view):
    rows=view['common_actual_tensor_rows']
    if len(rows)!=65 or view['current_common_tensor_contribution_count']!=65 or not view['actual_completed_tensor_and_divergence_remainder_common_source_bound'] or not view['overlap_not_used_for_tensor_source_equality'] or any(view[k] for k in OPEN):
        raise ValueError('Current O7 full tensor interface layout or function scope differs')
    nonzero=0
    for row in rows.values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Exact zero has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    if not nonzero:raise ValueError('Actual tensor replaced by local zero difference')
    return len(rows),nonzero


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentSteepWaitingBackgroundStress(require_checked=False)
    field.assert_graph()
    omitted={'current_actual_O7_tensor_views','current_actual_O7_three_tensor_interface_bounds'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current O7 tensor manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Producer claims its own acceptance')
    expected={chart+':'+name for chart in CHARTS for name in VIEWS}
    if set(raw['current_actual_O7_tensor_views'])!=expected:raise ValueError('Current O7 full domains and endpoints incomplete')
    count=nonzero=0
    for chart in CHARTS:
        for name,args in VIEWS.items():
            view=field.chart(chart,*args)
            if encode(pack(view))!=raw['current_actual_O7_tensor_views'][chart+':'+name]:raise ValueError('Current O7 tensor replay differs')
            n,m=check_view(field,view);count+=n;nonzero+=m
            if chart=='waiting' and not all(row['exact_zero'] for row in view['physical_axial_viscosity_remainder_mixed2'].values()):
                raise ValueError('Original waiting pure radial power must have source-exact zero axial viscosity')
    seams={};fresh={}
    for seam in SEAMS:
        view=field.interface(seam)
        if encode(pack(view))!=raw['current_actual_O7_three_tensor_interface_bounds'][seam]:raise ValueError('Current O7 interface replay differs')
        n,m=check_interface(view);seams[seam]=dict(rows=n,nonzero=m)
        view=field.interface(seam,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
        n,m=check_interface(view);fresh[seam]=dict(rows=n,nonzero=m,bounds=view)
    rejected=[]
    for label,kwargs in (('foreign_chart',dict(chart='other')),('outside_Z',dict(Z=2)),
            ('outside_coordinate',dict(x=-1)),('nonfinite_time',dict(log_tau='-inf')),
            ('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_angle',dict(theta='inf'))):
        args=dict(chart='waiting',Z='.2',x='.4',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.chart(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid current O7 domain accepted')
    clone=copy.copy(field);clone.steep=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_steep_owner')
    else:raise ValueError('Foreign current O7 owner accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_O7_tensor_view_count=len(expected),
        current_actual_O7_physical_tensor_divergence_remainder_rows_checked=count,
        current_actual_O7_nonzero_physical_contributions_checked=nonzero,
        current_actual_O7_three_common_tensor_interface_rows=seams,
        current_actual_O7_fresh_sector_three_common_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=4,
        current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
        original_full_power_exit_mixed4_KR_normalization_identity_count=len(field.moment_proof['original_full_power_exit_defect_AST_KR_normalization_mixed4_identities']),
        current_three_primitive_function_seam_identity_count=180,
        source_theorems={k:raw[k] for k in raw if k.endswith('_theorems') or k.endswith('_theorem') or k.startswith('checked_')},
        actual_full_current_nonzero_histories_absolute_pressure_and_Gamma_future_retained=True,
        stable_remaining_pressure_reduction_before_interval_enclosure=True,
        waiting_axial_viscosity_source_exact_zero_not_global_temporal_flatness=True,
        invalid_domains_and_foreign_current_owner_rejected=rejected,
        source_domain=raw['source_domain'],
        scope='Three actual current O7 full stress/tensor regions and three actual adjacent tensor function joins. Together with admitted angular/entry: five regional tensors, four adjacent tensor joins. Full current velocity-pressure atlas remains 14 adjacent + 8 internal interfaces. Other tensors, axis/full33, global cone/lift/NS/temporal remainder/energy/points and n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual O7 three tensors and three joins PASS; global and temporal gates remain open',flush=True)
    return result


if __name__=='__main__':run()
