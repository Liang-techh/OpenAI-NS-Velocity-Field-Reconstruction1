"""Scoped replay of current angular local stress/error difference bounds."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_angular_support_differences import (
    CurrentAngularSupportDifferences,exact_current_KR_source,EDGES,DISTANCES,GATES,OPEN,
    NAME,RECEIPT,HERE,PREFIX,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_pulse_interfaces_check import cloned_pulse_graph
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def difference_values(view):
    values=[]
    for key in ('local_beta_difference_rows','normalized_K_difference_rows','normalized_quadratic_K_difference_rows'):
        for jet in view[key]:values.extend(jet.coefficients)
    for rows in view['same_boundary_normalized_full_moment_difference_rows'].values():
        for jet in rows:values.extend(jet.coefficients)
    for grid in view['normalized_similarity_stress_difference_mixed3'].values():values.extend(grid.values())
    values.extend(view['normalized_axial_viscosity_remainder_difference_mixed2'].values())
    return values


def physical_values(view):
    return [row for grid in view['original_physical_stress_difference_mixed3_log_bounds'].values() for row in grid.values()]+list(view['original_physical_axial_viscosity_difference_mixed2_log_bounds'].values())


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentAngularSupportDifferences(require_checked=False)
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_local_angular_support_difference_views'}:
        raise ValueError('Current local angular source manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Angular difference producer claims admission')
    for key in ('original_angular_full_moment_transport_theorem','original_paper_stress_units','current_angular_local_difference_source_theorem'):
        if not all(raw[key]['identities'].values()):raise ValueError('Current full moment/source difference theorem missing')
    operator=raw['current_original_stress_and_viscosity_difference_operator_theorem']
    if len(operator['actual_AST_stress_mixed3_superposition_identities'])!=40 or not all(operator['actual_AST_stress_mixed3_superposition_identities'].values()):
        raise ValueError('Original stress mixed3 difference superposition incomplete')
    if len(operator['original_axial_viscosity_mixed2_superposition_identities'])!=6 or not all(operator['original_axial_viscosity_mixed2_superposition_identities'].values()):
        raise ValueError('Original axial viscosity mixed2 difference superposition incomplete')
    if len(raw['current_local_angular_support_difference_views'])!=12:raise ValueError('All four edges and three distances required')
    finite=zeros=physical=monotone=index=0
    for edge in EDGES:
        norms=[];logs=[]
        for h in DISTANCES:
            view=field.interface(edge,h)
            if encode(pack(view))!=raw['current_local_angular_support_difference_views'][index]:raise ValueError('Current local angular bound replay differs')
            index+=1;at_zero=endpoints(view['h'])==(mp.mpf(0),mp.mpf(0));norm=[];bounds=[]
            for value in difference_values(view):
                lo,hi=endpoints(value)
                if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite current normalized difference bound')
                finite+=1;norm.append(max(abs(lo),abs(hi)))
                if at_zero:
                    if lo or hi:raise ValueError('Nonzero exact angular endpoint difference')
                    zeros+=1
            rows=physical_values(view)
            if len(rows)!=26:raise ValueError('Angular stress3/error2 physical bound coverage incomplete')
            for row in rows:
                physical+=1
                if row['exact_zero']:
                    if row['log_absolute_upper'] is not None:raise ValueError('Zero physical difference has nonzero log bound')
                    bounds.append(mp.ninf)
                else:
                    if not all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])):raise ValueError('Nonfinite physical difference log bound')
                    bounds.append(endpoints(row['log_absolute_upper'])[1])
                if at_zero and not row['exact_zero']:raise ValueError('Nonzero physical endpoint difference')
            if not view['current_nonzero_actual_histories_and_absolute_pressure_preserved'] or any(view[k] for k in OPEN):
                raise ValueError('Current angular difference changed scope')
            norms.append(norm);logs.append(bounds)
        for groups in (norms,logs):
            for left,right in zip(groups,groups[1:]):
                if any(b>a for a,b in zip(left,right)):raise ValueError('Current difference majorant grew toward angular edge')
                monotone+=len(left)
        print('Current angular local stress/error replay PASS: '+edge['exact_edge'],flush=True)
    fresh=field.interface(EDGES[2],'.000001',Z='.631',log_tau='-2.7',viscosity='.8')
    if not all(row['exact_zero'] or all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])) for row in physical_values(fresh)):
        raise ValueError('Fresh current angular physical difference view nonfinite')
    rejected=[]
    for label in ('changed_current_steep_length','changed_current_waiting_normalization'):
        clone=cloned_pulse_graph(field.physical)
        if label=='changed_current_steep_length':clone.history.steep.Ts+=1
        else:clone.history.steep.logone+=1
        try:exact_current_KR_source(clone)
        except ValueError:rejected.append(label);continue
        raise ValueError('Wrong current exact angular normalization accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_source_identified_adjacent_interface_count=14,
        current_source_identified_internal_support_count=8,current_angular_support_edge_count=4,
        finite_normalized_local_difference_bounds_checked=finite,exact_endpoint_difference_zeros=zeros,
        original_physical_stress3_error2_difference_rows_checked=physical,monotone_difference_bound_comparisons=monotone,
        original_stress_mixed3_difference_operator_identity_count=40,original_axial_viscosity_mixed2_difference_operator_identity_count=6,
        current_exact_KR_source=raw['current_exact_KR_source'],
        wrong_current_normalizations_rejected=rejected,fresh_positive_time_viscosity_difference_view=fresh,
        source_common_KR_and_KR_squared_factors_never_materialized=True,
        actual_nonzero_A_E_P_K_histories_preserved_in_exact_operator_difference=True,
        local_interface_flatness_is_not_global_temporal_flatness=True,
        scope='Four current angular local stress3 and axial-viscosity error2 actual-reference difference bounds, same source boundary histories, exact logarithmic KR/KR^2 and original positive-time/viscosity physical operators. Combined quantitative interfaces, current global tensor/cone/NS/remainder/energy/points and actual temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current four angular local stress3/error2 flat differences PASS; global/time gates remain open',flush=True)
    return result


if __name__=='__main__':run()
