"""Compare pressure budget source with the accepted coherent core snapshot.

Rebuilds only source/datum adapters, not radial core coefficients or collar
continuation. The accepted local snapshot remains untouched.
"""
import json,pickle,hashlib
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def fingerprint(schedule):
    values={str(k):str(v) for k,v in schedule.metadata()['inputs'].items()}
    return dict(inputs=values,sha256=hashlib.sha256(json.dumps(values,sort_keys=True).encode()).hexdigest())


def run():
    current=source_profile();adapter=ContinuousPreheatPressure(current,quadrature_order=192)
    print('Rebuilding accepted coherent source adapters only',flush=True)
    bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',logC='5e151',logPstar='14',delta='1e-200',
        continuous_pressure=True,coherent_waiting=True,complete_preheat_components=True,component_Z_jet_depth=3)
    accepted=bundle['complete_preheat_pressure_datum'];schedule=bundle['profile'].schedule
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache/second_axial_core_Z03.pkl'
    with cache.open('rb') as f:stored=pickle.load(f)
    previous=stored['coefficients']['pressure_datum_components']['components']
    with mp.workdps(300):
        a=adapter.taylor_components(0,center=0)['components'];b=accepted.taylor_components(0,center=0)['components']
        assert set(a)==set(b)==set(previous)
        rows={}
        for stage in a:
            current_mass=a[stage]['value_at_Z0'];accepted_mass=b[stage]['value_at_Z0'];cached=previous[stage]['value_at_Z0']
            rows[stage]=dict(current_mass=current_mass,accepted_mass=accepted_mass,
                accepted_matches_snapshot=(accepted_mass==cached),current_matches_accepted=(current_mass==accepted_mass),
                stored_mass_change=accepted_mass-current_mass)
        assert all(row['accepted_matches_snapshot'] for row in rows.values()),'Canonical rebuild differs from accepted cache'
        out=dict(current_schedule=fingerprint(current.schedule),accepted_schedule=fingerprint(schedule),stages=rows,
            changed_stages=[name for name,row in rows.items() if not row['current_matches_accepted']],
            accepted_matches_snapshot=True,coherent_waiting_length=str(schedule.waiting_length),
            accepted_waiting_match_receipt=bundle['waiting_match_receipt'],
            source_alignment_complete=False,old_pressure_budget_transfer_certified=False,
            accepted_core_snapshot_modified=False,global_field_modified=False)
        print('Matched snapshot; changed stages',out['changed_stages'],flush=True)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(out),indent=2)+'\n')
        print('Accepted source matches snapshot; changed pressure stages',out['changed_stages'],flush=True)


if __name__=='__main__':run()
