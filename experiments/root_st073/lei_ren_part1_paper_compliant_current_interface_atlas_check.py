"""Check common source composition and actual compact-time trace envelopes."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_interface_atlas import (
    CurrentInterfaceAtlas,time_sector_theorem,validate_rows,STEMS,GATES,OPEN,
    NAME,RECEIPT,HERE,PREFIX,sha,pack,encode,endpoints,restore_value)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def check_view(field,view):
    if not view['actual_full_trace_not_local_difference'] or not view['exact_source_equality_from_checked_function_theorem']:
        raise ValueError('Actual trace and its defining function theorem required')
    if not view['all_angles_covered'] or any(view[key] for key in OPEN):raise ValueError('Physical interface scope differs')
    rows=view['physical_contribution_log_bounds'];validate_rows(field.ctx,rows)
    totals=view['cartesian_velocity_pressure_total_log_bounds']
    if len(totals)!=144:raise ValueError('Complete four-field spatial4/time1 total layout required')
    for key,row in rows.items():
        total=totals[key.rsplit(':',1)[0]]
        if not row['exact_zero']:
            if total is None or endpoints(total)[1]<endpoints(row['log_absolute_upper'])[1]:
                raise ValueError('Cartesian triangle upper omits a physical contribution')
    if any(value is not None and not all(mp.isfinite(v) for v in endpoints(value)) for value in totals.values()):
        raise ValueError('Nonfinite actual Cartesian trace total upper')
    return len(rows),len(totals)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentInterfaceAtlas(require_checked=False)
    omitted={'current_22_compact_time_interface_bounds','current_four_pulse_full_endpoint_baseline_log_bounds'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:
        raise ValueError('Current common source manifest differs')
    if any(raw[key] for key in GATES+OPEN):raise ValueError('Atlas producer claims its own acceptance')
    if encode(pack(time_sector_theorem()))!=raw['compact_time_bound_source_theorem']:
        raise ValueError('Actual original positive-time source transfer theorem differs')
    # Each source receipt was independently replayed by its own checker. This
    # composition checks their hashes, shared owner and exact source identities;
    # it does not repeat historical quadratures to assert their same theorem.
    if len(raw['source_receipt_composition'])!=len(STEMS):raise ValueError('Every current interface layer must be composed')
    names=(*field.adjacent,*field.internal)
    if len(names)!=22 or set(raw['current_22_compact_time_interface_bounds'])!=set(names):
        raise ValueError('Fourteen adjacent and eight internal affected traces required')
    contributions=totals=0
    for name in names:
        view=field.bounds(name)
        if encode(pack(view))!=raw['current_22_compact_time_interface_bounds'][name]:
            raise ValueError('Actual current compact-time trace replay differs: '+name)
        n,m=check_view(field,view);contributions+=n;totals+=m
    # Re-enclose the actual full native pulse field at all four rational edges.
    # This is separate from the h=0 local difference, which is identically zero.
    pulse_names=[name for name,row in field.internal.items() if row['chart']=='pulse_end']
    if set(raw['current_four_pulse_full_endpoint_baseline_log_bounds'])!=set(pulse_names):
        raise ValueError('All four actual pulse endpoint bounds required')
    nonzero=0
    for name in pulse_names:
        field.baseline.pop(name,None)
        actual=field.baseline_rows(name)
        if encode(pack(actual))!=raw['current_four_pulse_full_endpoint_baseline_log_bounds'][name]:
            raise ValueError('Current native full endpoint bound differs: '+name)
        nonzero+=sum(not row['exact_zero'] for row in actual.values())
    if not nonzero:raise ValueError('Actual pulse endpoint field was replaced by the zero local difference')
    fresh={};fresh_rows=0
    for name in names:
        view=field.bounds(name,log_tau=('-5','-2'),Z=('-0.8','0.8'))
        n,m=check_view(field,view);fresh_rows+=n
        fresh[name]=dict(requested_log_tau=view['requested_log_tau'],requested_Z=view['requested_Z'],
            cartesian_velocity_pressure_total_log_bounds=view['cartesian_velocity_pressure_total_log_bounds'])
    # Reject nonphysical requests and a cross-layer owner mismatch without
    # modifying any live graph or launching new numerical history evaluations.
    rejected=[]
    for label,kwargs in (('nonfinite_time',dict(log_tau=('-inf','-1'))),('outside_source_Z',dict(Z=(-2,1)))):
        try:field.bounds(names[0],**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid physical sector accepted')
    clone=copy.copy(field);clone.owners=tuple(copy.copy(owner) for owner in field.owners)
    clone.owners[3].physical=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('different_interface_physical_owner')
    else:raise ValueError('A different current interface owner was accepted')
    for label,attribute,key in (('different_adjacent_ledger_owner','adjacent','left_owner'),
            ('different_internal_ledger_beta','internal','normalized_beta_source')):
        clone=copy.copy(field);ledger=copy.deepcopy(getattr(field,attribute))
        ledger[next(iter(ledger))][key]='different source recipe'
        setattr(clone,attribute,ledger)
        try:clone.assert_ledgers()
        except ValueError:rejected.append(label);continue
        raise ValueError('Receipt output ledger was accepted without its live source binding')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_source_identified_adjacent_interface_count=14,
        current_source_identified_internal_support_count=8,current_composed_checked_source_layer_count=6,
        original_physical_spatial4_time1_contributions_checked=contributions,
        actual_cartesian_four_field_spatial4_time1_total_bounds_checked=totals,
        fresh_compact_time_physical_contributions_checked=fresh_rows,
        actual_nonzero_pulse_full_endpoint_contributions_checked=nonzero,
        fresh_compact_time_sector_total_bounds=fresh,invalid_sectors_and_cross_layer_owner_rejected=rejected,
        current_adjacent_source_ledger=raw['current_adjacent_source_ledger'],
        current_internal_source_ledger=raw['current_internal_source_ledger'],
        current_live_ledger_source_proof=raw['current_live_ledger_source_proof'],
        retained_local_angular_stress_error_theorems_source_bound=True,
        source_receipt_composition=raw['source_receipt_composition'],
        source_domain=raw['source_domain'],quantitative_domain=raw['quantitative_domain'],
        uncovered_obligations=raw['uncovered_obligations'],
        native_interval_consistency_not_used_as_source_identity=True,
        actual_full_pulse_support_bounds_not_zero_difference=True,
        scope='One current graph composes 14 adjacent/8 internal function traces. Actual velocity/absolute-pressure Cartesian spatial4 and fixed-position time1 contribution/total upper bounds cover all angles, source Z=[-1,1] and requested finite compact log(tau) sectors. Other joins, actual completed stress, global tensor/NS/remainder/energy/points and n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current 22-interface source composition and actual compact-time physical bounds PASS; global/time gates open',flush=True)
    return result


if __name__=='__main__':run()
