"""Current angular support source theorem, endpoint replay and scoped admission."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_angular_support_interfaces import (
    CurrentAngularSupportInterfaces,current_angular_support_proof,EDGES,GATES,OPEN,
    HERE,PREFIX,NAME,RECEIPT,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_pulse_interfaces_check import cloned_pulse_graph
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def endpoint_consistency(view):
    exact=view['exact_source_endpoint_packet'];native=view['unchanged_native_coordinate_consistency_packet']
    count=0
    for label,grid in exact['physical_mixed_derivatives_total_order_le4'].items():
        if len(grid)!=15:raise ValueError('Complete angular mixed4 source grid required')
        for key,value in grid.items():
            lo,hi=endpoints(value)
            if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite exact angular endpoint bound')
            dl,dh=endpoints(value-native['physical_mixed_derivatives_total_order_le4'][label][key])
            if not dl<=0<=dh:raise ValueError('Unchanged native endpoint consistency does not enclose exact source limit')
            count+=1
    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for jet in exact['actual_angular_bump_y_derivatives'] for v in jet.coefficients):
        raise ValueError('Exact angular edge retained local beta input')
    F=exact['swirl_factor_one_plus_h_Taylor']
    if endpoints(F[0])!=(mp.mpf(1),mp.mpf(1)) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in F.coefficients[1:]):
        raise ValueError('Exact angular source factor is not one')
    if endpoints(exact['energy_Taylor'][0])[0]<=0:raise ValueError('Complete current endpoint future energy lost positivity')
    pb=view['original_physical_spatial4_time1_endpoint_log_bounds']
    if len(pb)!=216:raise ValueError('Original angular physical operator coverage incomplete')
    for row in pb.values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Zero source contribution has nonzero log bound')
        elif not all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])):
            raise ValueError('Nonfinite physical angular endpoint contribution')
    if not view['current_complete_future_and_absolute_pressure_retained'] or any(view[k] for k in OPEN):
        raise ValueError('Current angular support view changed scope')
    return count,len(pb)


def wrong_sources_rejected(field):
    rejected=[]
    for label in ('changed_current_full_weight','different_flat_normalization','old_angular_provider'):
        p=cloned_pulse_graph(field.physical)
        if label=='changed_current_full_weight':
            p.history.outer.weights=dict(p.history.outer.weights);p.history.outer.weights['D']*=2
        elif label=='different_flat_normalization':
            p.history.outer.flat=copy.copy(p.history.outer.flat)
            p.history.outer.flat.normalization*=2;p.history.outer.normalization=p.history.outer.flat.normalization
        else:p.dispatch.providers['outer_angular']=p.base.outer
        try:current_angular_support_proof(p)
        except ValueError:rejected.append(label);continue
        raise ValueError('Wrong current angular support source accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),live_current_graph_unmodified=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentAngularSupportInterfaces(require_checked=False)
    manifest={k:v for k,v in raw.items() if k!='current_whole_Z_angular_support_endpoint_views'}
    if encode(pack(field.manifest()))!=manifest:raise ValueError('Current angular support manifest differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Angular producer claims admission')
    proof=current_angular_support_proof(field.physical)
    if encode(pack(proof))!=raw['current_angular_support_source_proof']:raise ValueError('Current angular support defining-function proof differs')
    ftc=proof['weighted_integral_FTC_and_flat_beta_source_theorem']
    if not ftc['passed'] or len(ftc['identities'])!=49 or not all(ftc['identities'].values()):
        raise ValueError('Complete original weighted one-sided source theorem missing')
    bridge=proof['exact_normalization_and_five_weight_defining_source']
    if not bridge['passed'] or not all(bridge['current_preserved_source_graph'].values()) or not all(bridge['exact_full_past_future_five_weight_integral_identities'].values()):
        raise ValueError('Original full weights and flat beta do not share a defining source')
    theorem=proof['current_angular_endpoint_source_theorem']
    for key in ('current_primitive_mixed4_identities','current_velocity_pressure_mixed4_identities'):
        if len(theorem[key])!=240 or not all(theorem[key].values()):raise ValueError('Four angular source mixed4 grids incomplete')
    if len(raw['current_whole_Z_angular_support_endpoint_views'])!=4:raise ValueError('All four current angular edge views required')
    rows=physical=0
    for edge,saved in zip(EDGES,raw['current_whole_Z_angular_support_endpoint_views']):
        view=field.endpoint(edge)
        if encode(pack(view))!=saved:raise ValueError('Exact current angular source endpoint replay differs')
        n,m=endpoint_consistency(view);rows+=n;physical+=m
        print('Current angular support replay PASS: '+edge['exact_edge'],flush=True)
    fresh=[]
    for edge in EDGES:
        view=field.endpoint(edge,Z='.577',log_tau='-2.4',theta='.7')
        n,m=endpoint_consistency(view);rows+=n;physical+=m
        fresh.append(dict(exact_source_edge=edge['exact_edge'],Z=view['Z'],log_tau='-2.4',theta='.7',
            original_physical_spatial4_time1_endpoint_log_bounds=view['original_physical_spatial4_time1_endpoint_log_bounds']))
    ledger=raw['current_internal_support_ledger']
    if len(ledger)!=8 or not all(v['current_full_field_interface_trace_admitted'] for v in ledger.values()):
        raise ValueError('Combined current eight-support source ledger incomplete')
    mutations=wrong_sources_rejected(field)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_source_identified_adjacent_interface_count=14,
        current_source_identified_internal_support_count=8,current_angular_support_edge_count=4,
        one_sided_weighted_FTC_and_flat_beta_identity_count=49,
        current_angular_primitive_mixed4_source_identity_count=240,
        current_angular_velocity_pressure_mixed4_source_identity_count=240,
        current_native_endpoint_consistency_rows_checked=rows,
        original_physical_spatial4_time1_endpoint_rows_checked=physical,
        fresh_current_angular_source_endpoint_views=fresh,wrong_current_sources=mutations,
        current_internal_support_ledger=ledger,retained_current_pulse_support_count=4,
        same_defining_normalization_and_full_weight_source_consumed=True,
        source_function_one_sided_limits_precede_rounded_coordinate_consistency=True,
        inherited_X_absolute_P_and_complete_future_not_reset=True,
        all_current_angular_local_stress_error_difference_bounds_certified=False,
        scope='Four current angular internal support function/velocity-pressure mixed4 traces and original physical spatial4-time1 source maps. Combined 14 adjacent/8 internal source traces are identified; quantitative all-interface, full angular local stress/error bounds, global tensor/NS/remainder/energy/points and temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current four angular support source traces PASS; all eight internal source traces identified',flush=True)
    return result


if __name__=='__main__':run()
