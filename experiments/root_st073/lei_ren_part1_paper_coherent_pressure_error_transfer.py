"""Transfer pressure enclosures onto the accepted coherent-waiting source.

The 11 unchanged stages transfer analytically because all schedule inputs
except waiting length are identical. The three terminal stages are rebuilt.
No pressure anchor or accepted field is changed.
"""
import json,copy
from decimal import Decimal
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
from lei_ren_part1_paper_continuous_angular_schedule import install_continuous_angular_schedule
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_variable_preheat_interval_integrals import VariablePreheatIntervalIntegrals
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval,fixed_beta_bounds
from lei_ren_part1_paper_uniform_pressure_high_derivatives import run as high_derivatives
from lei_ren_part1_paper_uniform_core_pressure_propagation import run as core_propagation


def exact_scalar(record):return mp.make_mpf(tuple(record['exact_mpf_tuple']))


def accepted_profile():
    receipt=json.loads(Path(__file__).with_name('lei_ren_part1_paper_coherent_pressure_source_alignment.json').read_text())
    assert receipt['accepted_matches_snapshot']
    current=receipt['current_schedule']['inputs'];chosen=receipt['accepted_schedule']['inputs']
    assert set(current)==set(chosen)
    for k in current:
        if k!='waiting_length':assert Decimal(current[k])==Decimal(chosen[k]),'Nonterminal source parameter changed: '+k
    assert set(receipt['changed_stages'])=={'waiting','heat_collar','exterior_power_tail'}
    schedule=PaperOuterSchedule(**chosen)
    profile=CorrectedSourceProfile(schedule=schedule,match_waiting=False,precision=260)
    install_continuous_angular_schedule(profile.schedule,precision=260,primitive_precision=100)
    profile.schedule._coherent_preheat_waiting=True
    return profile,receipt


def run():
    base=Path(__file__).parent;profile,alignment=accepted_profile()
    datum=ContinuousPreheatPressure(profile,quadrature_order=192);e=ScheduleEndpointEnclosures(profile.schedule);iv=e.iv
    fixed=json.loads((base/'lei_ren_part1_paper_uniform_fixed_beta_error.json').read_text())
    full=json.loads((base/'lei_ren_part1_paper_uniform_preheat_error_budget.json').read_text())
    with mp.workdps(e.precision+40):
        finite=datum.taylor_components(0,center=0)['components']
        for name,row in finite.items():
            assert row['value_at_Z0']==exact_scalar(alignment['stages'][name]['accepted_mass']),'Accepted source replay mismatch: '+name
        true_terminal=e.terminal_preheat_bounds()
        rebuilt={'waiting':VariablePreheatIntervalIntegrals(e).integrate_negative_stage('waiting')['normalized_mass_at_Z0_interval'],
            'heat_collar':true_terminal['heat_collar_mass_interval'],'exterior_power_tail':true_terminal['exterior_mass_interval']}
        for name,truth in rebuilt.items():
            mass=finite[name]['value_at_Z0'];error=truth-iv.mpf(mass)
            fixed['stages'][name]=encode(dict(beta=0,true_mass_interval=truth,finite_mass=mass,mass_error_interval=error,
                source_scope='rebuilt accepted coherent-waiting schedule',
                uniform_derivative_error_upper_bounds=dict(zip(('value','first_Z','second_Z'),[endpoints(v)[1] for v in fixed_beta_bounds(iv,error,0,'.8')]))))
        total=[iv.mpf(0) for _ in range(3)]
        for row in fixed['stages'].values():
            err=read_interval(iv,row['mass_error_interval']);bounds=fixed_beta_bounds(iv,err,row['beta'],'.8')
            total=[x+y for x,y in zip(total,bounds)]
        names=('value','first_Z','second_Z')
        fixed['normalized_uniform_derivative_error_upper_bounds']=encode(dict(zip(names,[endpoints(v)[1] for v in total])))
        fixed['normalized_uniform_weighted_C2_error_upper']=encode(endpoints(total[0]+total[1]+total[2]/2)[1])
        fixed['accepted_source_alignment_certified']=True;fixed['accepted_schedule_sha256']=alignment['accepted_schedule']['sha256']
        fixed_path=base/'lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json'
        fixed_path.write_text(json.dumps(fixed,indent=2)+'\n')
        flatten=[iv.mpf(exact_scalar(full['flatten_uniform_absolute_error_upper_bounds'][name])) for name in names]
        complete=[x+y for x,y in zip(total,flatten)]
        full['normalized_uniform_derivative_error_upper_bounds']=encode(dict(zip(names,[endpoints(v)[1] for v in complete])))
        full['normalized_uniform_weighted_C2_error_upper']=encode(endpoints(complete[0]+complete[1]+complete[2]/2)[1])
        full['fixed_beta_receipt']=fixed_path.name;full['accepted_source_alignment_certified']=True
        full['accepted_schedule_sha256']=alignment['accepted_schedule']['sha256']
        full['terminal_stages_rebuilt']=list(rebuilt)
        full_path=base/'lei_ren_part1_paper_coherent_uniform_preheat_error_budget.json'
        full_path.write_text(json.dumps(full,indent=2)+'\n')
        print('Accepted coherent pressure weighted C2 error',mp.nstr(endpoints(complete[0]+complete[1]+complete[2]/2)[1],16),flush=True)
    high_path=base/'lei_ren_part1_paper_coherent_uniform_pressure_high_derivatives.json'
    high_derivatives(fixed_receipt=fixed_path,full_receipt=full_path,output=high_path)
    core_propagation(profile=profile,pressure_receipt=high_path,
        output=base/'lei_ren_part1_paper_coherent_uniform_core_pressure_propagation.json',
        profile_origin='accepted coherent-waiting schedule, replayed against second_axial_core_Z03 snapshot')


if __name__=='__main__':run()
