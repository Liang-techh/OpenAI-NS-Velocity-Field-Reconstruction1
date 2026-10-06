"""Current source end/flatten physical traces, separate from all-interface gates."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_physical_interfaces import (
    CurrentPhysicalInterfaces,native_end_flatten_trace_proof,interface_ledger,
    functional_join_identities,supported_beta_edges,HERE,PREFIX,NAME,RECEIPT,GATES,OPEN,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def changed_terminal_sources_rejected(field):
    rejected=[]
    for label in ('old_flatten_provider','different_pressure_datum','different_U_constants'):
        p=copy.copy(field.physical);p.flatten=copy.copy(field.physical.flatten)
        p.flatten.inlet=copy.copy(p.flatten.inlet)
        if label=='old_flatten_provider':p.flatten=p.base.dispatch.provider('flatten')
        elif label=='different_pressure_datum':p.flatten.inlet.datum=object()
        else:p.flatten.inlet.constants=dict(p.flatten.inlet.constants)
        try:native_end_flatten_trace_proof(p)
        except ValueError:rejected.append(label);continue
        raise ValueError('Changed terminal source accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPhysicalInterfaces(require_checked=False)
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_end_flatten_views'}:
        raise ValueError('Current interface graph/source/scope differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Interface producer claims acceptance')
    proof=native_end_flatten_trace_proof(field.physical);canonical=functional_join_identities()
    if encode(pack(proof))!=raw['current_end_flatten_native_source_proof'] or encode(pack(canonical))!=raw['recomputed_arbitrary_shared_function_endpoint_theorem']:
        raise ValueError('Actual endpoint theorem differs')
    if len(proof['mixed4_identities'])!=60 or not all(proof['mixed4_identities'].values()) or len(canonical['identities'])!=155:
        raise ValueError('Native / canonical endpoint derivative identities omitted')
    if not all(proof['exact_positive_velocity_and_pressure_log_source_identities'].values()) or not proof['exact_pressure_source_bound_before_its_numeric_cap_enclosure']:
        raise ValueError('Exact pressure/velocity source was replaced by a cap')
    ledger=interface_ledger(field.physical)
    if len(ledger)!=14 or [k for k,v in ledger.items() if v['current_functional_mixed4_trace_admitted']]!=['end_flatten']:
        raise ValueError('Fourteen-seam ledger or restricted acceptance differs')
    edges=supported_beta_edges(field.physical)
    if len(edges)!=8 or any(v['current_full_field_interface_trace_admitted'] for v in edges.values()):
        raise ValueError('Interior beta-edge inventory exceeds source-transfer scope')
    c=field.ctx;native_rows=physical_rows=0;bounds={}
    if set(raw['current_end_flatten_views'])!={'whole_Z','fresh'}:raise ValueError('Whole-Z and fresh traces required')
    for name,args in {'whole_Z':([-1,1],'-1',None),'fresh':('.419','-2.31','.41')}.items():
        point=field.endpoint(*args)
        if encode(pack(point))!=raw['current_end_flatten_views'][name]:raise ValueError('Current endpoint call differs')
        if any(point[k] for k in OPEN):raise ValueError('Endpoint scope exceeds one trace')
        for rows in point['native_fixed_Ev0_velocity_and_Pstar_squared_pressure_difference_enclosures'].values():
            if len(rows)!=15:raise ValueError('Native mixed4 grid incomplete')
            for value in rows.values():
                lo,hi=endpoints(value)
                if not all(mp.isfinite(v) for v in (lo,hi)) or not lo<=0<=hi:
                    raise ValueError('Source-proved endpoint equality contradicts native enclosure')
                native_rows+=1
        left=point['left_physical_trace'];right=point['right_physical_trace']
        if endpoints(left['source_logR_enclosure'])!=endpoints(right['source_logR_enclosure']):
            raise ValueError('Physical trace radius source differs')
        view_bounds={}
        for category in ('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative'):
            def leaves(node,path=()):
                if 'exact_zero' in node:yield path,node
                else:
                    for key,value in node.items():yield from leaves(value,path+(key,))
            a=dict(leaves(left[category]));b=dict(leaves(right[category]))
            if set(a)!=set(b):raise ValueError('Physical moving-basis/time contributions differ')
            for key,row in a.items():
                other=b[key]
                if endpoints(row['physical_lambda_exponent'])!=endpoints(other['physical_lambda_exponent']):
                    raise ValueError('Physical trace lambda derivative power differs')
                if bool(row['exact_zero'])!=bool(other['exact_zero']):
                    raise ValueError('Source zero identity differs between traces')
                if row['exact_zero']:upper=None
                else:
                    values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other)]
                    if not all(mp.isfinite(v) for v in values):raise ValueError('Nonfinite physical trace upper')
                    upper=max(values)
                view_bounds[category+':'+':'.join(key)]=upper;physical_rows+=2
        if len(view_bounds)!=216:raise ValueError('Physical trace spatial4/time1 coverage incomplete')
        bounds[name]=view_bounds
    mutations=changed_terminal_sources_rejected(field)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_affected_interface_ledger=ledger,
        current_supported_beta_edge_inventory=edges,
        current_native_source_mixed4_identities_checked=60,
        recomputed_arbitrary_shared_function_endpoint_identities_checked=155,
        native_mixed4_difference_consistency_rows_checked=native_rows,
        physical_two_sided_trace_contributions_checked=physical_rows,
        finite_current_common_physical_trace_log_uppers=bounds,
        wrong_terminal_sources=mutations,
        source_identity_precedes_bounds_and_overlap=True,
        current_end_flatten_scope='All Z in [-1,1] and positive tau: source-identified ordinary mixed4 and original physical spatial4/fixed-x time1 traces at Rv; finite bounds reported at selected positive time sectors. Other seams, global tensor and temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current interfaces PASS:14 ledger entries,60 native identities,120 native checks,864 physical contributions',flush=True)
    return result


if __name__=='__main__':run()
