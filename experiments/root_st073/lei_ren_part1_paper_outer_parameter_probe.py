"""Screen actual new preheat sources before expensive fresh-core generation.

Each trial rebuilds a conditional fourteen-stage pressure datum and nominal
waiting root. Point cone checks select candidates; they are not whole-stage
admissibility proofs. Radial placement is normalized using exact homogeneity.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_reference_endpoint_targets import accepted_parameters
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
from lei_ren_part1_paper_continuous_angular_schedule import install_continuous_angular_schedule
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals
from lei_ren_part1_paper_variable_preheat_interval_integrals import VariablePreheatIntervalIntegrals
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets,BETA2
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def trial(Md):
    with mp.workdps(340):
        lp=max(mp.mpf(14),mp.exp(mp.mpf(Md))+11)
        params=dict(accepted_parameters(),Md=Md,logPstar=mp.nstr(lp,280),waiting_length='0',
            logRref=mp.nstr(mp.log(110)+10*(mp.mpf('5e151')+lp),300))
        profile=CorrectedSourceProfile(schedule=PaperOuterSchedule(**params),match_waiting=True,precision=260)
        install_continuous_angular_schedule(profile.schedule,precision=260,primitive_precision=100)
        profile.schedule._coherent_preheat_waiting=True
        s=profile.schedule;e=ScheduleEndpointEnclosures(s,precision=300);c=e.iv
        v=VariablePreheatIntervalIntegrals(e);h=HighOrderPreheatIntegrals(e)
        refined=json.loads((HERE/'lei_ren_part1_paper_candidate_pressure_mass_refinement.json').read_bytes())
        masses=dict(reference_extension=c.mpf('2.5'),slope_transition_ref=
            read_interval(c,refined['stages']['slope_transition_ref']['refined_mass']))
        for stage in ('axial_turnoff','slope_transition_mu','power_buffer','pulse_reserved','power_buffer_rel',
                      'steep_transition_in','steep_power','steep_transition_out','waiting'):
            row=h.integrate_stage(stage,panels=64,order=12) if stage in (
                'slope_transition_mu','steep_transition_in','steep_transition_out') else v.integrate_negative_stage(stage)
            masses[stage]=row['normalized_mass_at_Z0_interval']
        flatten=c.mpf(e.stage_pressure_upper('z_flatten',axial_radius='.8')['mass_upper'])
        terminal=e.terminal_preheat_bounds()
        masses.update(heat_collar=terminal['heat_collar_mass_interval'],exterior_power_tail=terminal['exterior_mass_interval'])
        m2=sum((m for name,m in masses.items() if name in BETA2),c.mpf(0))
        m0=sum((m for name,m in masses.items() if name not in BETA2),c.mpf(0))
        inputs=s.metadata()['inputs'];sha=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
        datum=dict(stage_count_total=14,all_14_true_pressure_stages_included=True,new_schedule_sha256=sha,
            fixed_beta2_mass_upper_normalized=m2,fixed_beta0_mass_upper_normalized=m0,
            flatten_true_mass_upper_normalized=c.mpf([0,endpoints(flatten)[1]]),Pstar_squared=c.exp(2*c.mpf(inputs['logPstar'])))
        encoded=encode(datum)
        z=IntervalTaylor(c,[c.mpf(['.49','.51']),c.mpf(1)])
        A=(1+z*z).reciprocal()*c.exp(c.mpf(inputs['logPstar']))
        P0=IntervalTaylor(c,pressure_jets(c,z[0],1,encoded)['physical_pressure_coefficients'])
        # Outer inertial stresses /F are exactly linear in R. Compute their
        # coefficients at unit Rref and discard only the positive O(1/R)
        # shear contribution when screening the asymptotic cone. A failed
        # negative-direction condition already excludes every finite R.
        inlet=evaluate_transition(c,z,c.mpf(inputs['delta']),c.mpf(1),A,P0,1)
        inlet.update(pressure_schedule_Md=Md,pressure_datum_kind='rebuilt trial preheat source')
        field=IntervalOuterAxialTurnoff(c,z,c.mpf(inputs['delta']),inlet,Md,Md)
        samples=[]
        for phase in ('.25','.5','.75'):
            p=field.evaluate_phase(phase);F=p['F'][0];R=p['R']
            st=p['normalized_stress']['S_theta_over_F'];sz=p['normalized_stress']['S_z_over_F']
            nt=p['stress']['I_theta']/(F*R);nz=p['stress']['I_z']/(F*R)
            con=cone(c,st,sz,nt,nz)
            direction=st*nt+sz*nz
            samples.append(dict(phase=phase,inertial_direction_per_R=direction,
                normalized_inertial_cone=con,
                direction_obstructed_for_all_placements=endpoints(direction)[0]>0))
        result=dict(Md=Md,logPstar=inputs['logPstar'],new_schedule_inputs=inputs,new_schedule_sha256=sha,
            fourteen_stage_trial_pressure_datum=datum,pressure_stages=masses,flatten_mass_upper=flatten,
            sampled_inertial_cone=samples,
            sampled_all_inertial_cones_certified=all(r['normalized_inertial_cone']['admissible_cone_certified'] for r in samples),
            any_sample_direction_obstructed=any(r['direction_obstructed_for_all_placements'] for r in samples),
            only_local_axial_family=['.49','.51'],whole_stage_cone_certified=False,
            finite_radial_placement_shear_terms_certified=False,
            waiting_root_certified=False,unknown_paper_constants_verified=False,
            original_parameter_errors_enclosed=False,new_core_generated=False,temporal_recursion=False)
        print('Outer trial Md',Md,'logPstar',mp.nstr(lp,8),'sampled cones',
            [r['normalized_inertial_cone']['status'] for r in samples],flush=True)
        return result


def run():
    rows=[trial(md) for md in ('1.1','1.3','2','3')]
    result=dict(trials=rows,screen_before_new_core_generation=True,
        whole_outer_admissibility_certified=False,
        next_action='choose a source passing point screening, then certify complete axial/radial cells before any expensive new core',
        temporal_recursion=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_candidate_pressure_mass_refinement.json',
        'lei_ren_part1_paper_Md11_pressure_datum.py','lei_ren_part1_paper_interval_outer_slope_field.py',
        'lei_ren_part1_paper_interval_outer_axial_turnoff.py','lei_ren_part1_paper_interval_exit_stress_enclosure.py',
        'lei_ren_part1_paper_outer.py','lei_ren_part1_paper_corrected_profile.py',
        'lei_ren_part1_paper_high_order_preheat_integrals.py','lei_ren_part1_paper_variable_preheat_interval_integrals.py',
        'lei_ren_part1_paper_schedule_endpoint_enclosures.py','lei_ren_part1_paper_continuous_angular_schedule.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':run()
