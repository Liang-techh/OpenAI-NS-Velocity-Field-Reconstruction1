"""Current four-support source transfer and local flat difference replay."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import (
    CurrentPulseSupportInterfaces,current_support_source_proof,EDGES,DISTANCES,
    GATES,OPEN,HERE,PREFIX,NAME,RECEIPT,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_pulse_interfaces_check import cloned_pulse_graph
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def difference_values(point):
    result=[]
    for jet in point['unscaled_beta_forcing_difference_rows']:result.extend(jet.coefficients)
    for jets in point['source_primitive_differences'].values():
        for jet in jets:result.extend(jet.coefficients)
    for name in ('similarity_stress_difference_mixed3','physical_stress_difference_mixed3_coefficients',
        'physical_three_component_error_difference_mixed2_coefficients'):
        for sectors in point[name].values():
            for grid in sectors.values():result.extend(grid.values())
    for name in ('similarity_velocity_difference_mixed4','current_native_velocity_pressure_difference_mixed4'):
        for grid in point[name].values():result.extend(grid.values())
    return result


def wrong_sources_rejected(field):
    rejected=[]
    for label in ('old_selected_pulse_provider','different_current_complete_future','different_flat_beta_owner'):
        p=cloned_pulse_graph(field.physical)
        if label=='old_selected_pulse_provider':p.dispatch.providers['pulse_end']=p.base.pulse
        elif label=='different_current_complete_future':p.heat.future=object()
        else:p.pulse.flat=copy.copy(p.pulse.flat);p.pulse.flat.normalization*=2
        try:current_support_source_proof(p)
        except ValueError:rejected.append(label);continue
        raise ValueError('Changed current support source accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),live_current_owners_unmodified=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseSupportInterfaces(require_checked=False)
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_whole_Z_support_difference_views'}:
        raise ValueError('Current pulse support manifest/source differs')
    proof=current_support_source_proof(field.physical)
    if encode(pack(proof))!=raw['current_pulse_support_source_proof']:raise ValueError('Current support defining-function theorem differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Support producer claims acceptance')
    ftc=proof['exact_one_sided_weighted_integral_theorem']
    if not ftc['passed'] or len(ftc['identities'])!=49 or not all(ftc['identities'].values()):raise ValueError('Complete one-sided FTC/flat beta source theorem missing')
    source=proof['current_source_endpoint_difference_theorem']
    if not source['passed'] or len(source['difference_mixed4_identities'])!=60 or not all(source['difference_mixed4_identities'].values()):
        raise ValueError('Current primitive/forcing endpoint difference theorem missing')
    if len(raw['current_whole_Z_support_difference_views'])!=len(EDGES)*len(DISTANCES):raise ValueError('Four complete source-edge difference ledgers required')
    finite=zeros=monotone=physical_count=0;bounds_by_edge={};index=0
    for edge in EDGES:
        norms=[];physical_norms=[]
        for h in DISTANCES:
            point=field.interface(edge,h)
            if encode(pack(point))!=raw['current_whole_Z_support_difference_views'][index]:raise ValueError('Current support bound replay differs')
            index+=1;at_zero=endpoints(point['h'])==(mp.mpf(0),mp.mpf(0));values=difference_values(point);norm=[]
            for value in values:
                lo,hi=endpoints(value)
                if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite current support difference bound')
                finite+=1;norm.append(max(abs(lo),abs(hi)))
                if at_zero:
                    if lo or hi:raise ValueError('Exact source-edge difference is nonzero')
                    zeros+=1
            pb=point['original_cartesian_spatial4_time1_difference_log_bounds'];logs=[]
            if len(pb)!=216:raise ValueError('Original physical spatial4/time1 difference coverage incomplete')
            for value in pb.values():
                if value['exact_zero']:
                    if value['log_absolute_upper'] is not None:raise ValueError('Zero physical difference has a nonzero log bound')
                    logs.append(mp.ninf)
                else:
                    lo,hi=endpoints(value['log_absolute_upper'])
                    if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite original physical difference log bound')
                    logs.append(hi)
                if at_zero and not value['exact_zero']:raise ValueError('Nonzero physical source-edge difference')
                physical_count+=1
            if not point['actual_nonzero_histories_and_current_complete_future_preserved'] or not point['h_zero_is_zero_difference_not_zero_actual_field']:
                raise ValueError('Support comparison dropped actual history')
            if any(point[k] for k in OPEN):raise ValueError('Local support bound overclaims global admission')
            norms.append(norm);physical_norms.append(logs)
        for rows in (norms,physical_norms):
            for a,b in zip(rows,rows[1:]):
                if any(right>left for left,right in zip(a,b)):raise ValueError('Current flat difference envelope grew toward the edge')
                monotone+=len(a)
        bounds_by_edge[edge['exact_edge']]=dict(difference_rows=len(norms[0]),physical_rows=216,source_distance_sectors=DISTANCES)
        print('Current pulse support replay PASS: '+edge['exact_edge'],flush=True)
    ledger=raw['current_internal_support_ledger']
    if len(ledger)!=8 or sum(row['current_full_field_interface_trace_admitted'] for row in ledger.values())!=4:
        raise ValueError('Current four pulse / four open angular edge split differs')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_source_identified_adjacent_interface_count=14,
        current_source_identified_internal_support_count=4,current_pulse_support_edge_count=4,
        one_sided_weighted_FTC_and_flat_beta_identity_count=49,
        source_endpoint_primitive_forcing_mixed4_difference_identity_count=60,
        current_instantiation_applies_to_all_four_exact_source_edges=True,
        finite_local_difference_bounds_checked=finite,exact_source_endpoint_difference_zeros=zeros,
        monotone_difference_bound_comparisons=monotone,original_physical_spatial4_time1_difference_rows_checked=physical_count,
        bounded_edge_sectors=bounds_by_edge,wrong_current_sources=wrong_sources_rejected(field),
        current_internal_support_ledger=ledger,remaining_internal_support_transfers=raw['remaining_internal_support_transfers'],
        source_identity_and_FTC_limits_precede_numeric_enclosures=True,
        inherited_nonzero_histories_not_reset=True,
        scope='Four current pulse-end internal support mixed4/physical spatial4-time1 traces, current local stress3/error2 flat difference coefficients and finite positive-time difference bounds. Four angular edges, quantitative all-interface/global tensor/NS/point/energy and true temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current four pulse support interfaces PASS; current coefficients/full histories, original local stress/error differences',flush=True)
    return result


if __name__=='__main__':run()
