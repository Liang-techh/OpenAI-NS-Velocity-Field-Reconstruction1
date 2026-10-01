"""Companion fourteen-stage preheat datum for Md=1.1, without old-chain edits.

Explicit Md/Td/Pstar inequalities are checked. Unknown paper smallness
constants and the waiting root/source parameter errors remain conditional.
This datum is not installed into the old certified core or moment inverse.
"""
import hashlib
import json
import time
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
from lei_ren_part1_paper_candidate_pressure_function import load_datum,q_jets
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
BETA2={'reference_extension','slope_transition_ref','axial_turnoff',
       'slope_transition_mu','power_buffer','pulse_reserved'}


def pressure_jets(c,z,order,datum):
    """Physical pressure jets from the new masses, with all-stage flatten bound."""
    from lei_ren_part1_paper_candidate_pressure_axis_jets import _interval_from_exact
    if not isinstance(order,int) or isinstance(order,bool) or order<0:
        raise ValueError('nonnegative integer Taylor order required')
    z=c.mpf(z)
    if endpoints(z)[0]<-1 or endpoints(z)[1]>1:raise ValueError('real axial family must lie in [-1,1]')
    if datum['stage_count_total']!=14 or not datum['all_14_true_pressure_stages_included']:
        raise ValueError('fourteen-stage pressure datum required')
    read=lambda name:_interval_from_exact(c,datum[name])
    m2=read('fixed_beta2_mass_upper_normalized');m0=read('fixed_beta0_mass_upper_normalized')
    flatten=read('flatten_true_mass_upper_normalized');scale=read('Pstar_squared')
    rho=c.mpf('.25');q_lower=(1-rho)**2
    coeff=q_jets(c,z,order);rows=[]
    for n,a in enumerate(coeff):
        if n==0:remainder=flatten
        else:
            b=flatten/q_lower**2/rho**n
            remainder=c.mpf([endpoints(-b)[0],endpoints(b)[1]])
        rows.append(-scale*(m2*a+(m0 if n==0 else 0)+remainder))
    return dict(center_Z=c.mpf(z),ordinary_Taylor_order=order,
        physical_pressure_coefficients=rows,new_schedule_sha256=datum['new_schedule_sha256'],
        all_14_stages_included=True,pressure_units='physical_P',
        original_parameter_errors_enclosed=False,new_core_generated=False)


def run():
    started=time.monotonic()
    _,old_hashes=load_datum()
    params=accepted_parameters();params=dict(params,Md='1.1',waiting_length='0')
    with mp.workdps(340):
        print('Constructing Md=1.1 schedule and recomputing nominal waiting root',flush=True)
        schedule=PaperOuterSchedule(**params)
        profile=CorrectedSourceProfile(schedule=schedule,match_waiting=True,precision=260)
        install_continuous_angular_schedule(profile.schedule,precision=260,primitive_precision=100)
        profile.schedule._coherent_preheat_waiting=True
        schedule=profile.schedule
        inputs=schedule.metadata()['inputs']
        digest=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
        e=ScheduleEndpointEnclosures(schedule,precision=300);c=e.iv
        md=c.mpf(inputs['Md']);td=c.exp(md)+10;lp=c.mpf(inputs['logPstar'])
        if endpoints(md)[0]<=1 or endpoints(lp-td)[0]<=0:
            raise ValueError('paper Md>1 / logPstar>Td failed')
        v=VariablePreheatIntervalIntegrals(e);h=HighOrderPreheatIntegrals(e)
        refined=json.loads((HERE/'lei_ren_part1_paper_candidate_pressure_mass_refinement.json').read_bytes())
        rows={}
        rows['reference_extension']=dict(beta=2,mass=c.mpf('2.5'),method='exact reference power primitive')
        # Eq.4.6 on y in [0,1] has no Md, waiting, mu or delta dependence.
        rows['slope_transition_ref']=dict(beta=2,
            mass=read_interval(c,refined['stages']['slope_transition_ref']['refined_mass']),
            method='unchanged Eq.4.6 directed high-order mass, hash-linked transfer')
        for stage in ('axial_turnoff','slope_transition_mu','power_buffer','pulse_reserved',
                      'power_buffer_rel','steep_transition_in','steep_power','steep_transition_out','waiting'):
            print('Enclosing new pressure mass:',stage,flush=True)
            result=h.integrate_stage(stage,panels=64,order=12) if stage in (
                'slope_transition_mu','steep_transition_in','steep_transition_out') else v.integrate_negative_stage(stage)
            rows[stage]=dict(beta=2 if stage in BETA2 else 0,
                mass=result['normalized_mass_at_Z0_interval'],integration_receipt=result)
        flatten=e.stage_pressure_upper('z_flatten',axial_radius='.8')
        upper=c.mpf(flatten['mass_upper'])
        rows['z_flatten']=dict(beta='variable [0,2]',mass=c.mpf([0,endpoints(upper)[1]]),
            method='complete positive mass and uniform axial derivative envelope',integration_receipt=flatten)
        terminal=e.terminal_preheat_bounds()
        for stage,key in (('heat_collar','heat_collar_mass_interval'),('exterior_power_tail','exterior_mass_interval')):
            rows[stage]=dict(beta=0,mass=terminal[key],integration_receipt=terminal)
        if len(rows)!=14:raise AssertionError('fourteen true pressure stages required')
        m2=sum((row['mass'] for row in rows.values() if row['beta']==2),c.mpf(0))
        m0=sum((row['mass'] for row in rows.values() if row['beta']==0),c.mpf(0))
        result=dict(new_schedule_inputs=inputs,new_schedule_sha256=digest,stage_count_total=14,
            stored_y_d=str(schedule.y_d),stored_mu=str(schedule.mu),
            stages=rows,pressure_units='physical_P',
            fixed_beta2_mass_upper_normalized=m2,fixed_beta0_mass_upper_normalized=m0,
            flatten_true_mass_upper_normalized=c.mpf([0,endpoints(upper)[1]]),Pstar_squared=c.exp(2*lp),
            explicit_Md_greater_than_one=True,explicit_logPstar_greater_than_Td=True,
            logPstar_minus_Td=lp-td,all_14_true_pressure_stages_included=True,
            finite_pressure_quadrature_approximation_used=False,pressure_parameter_order_truncated=False,
            waiting_length_recomputed=True,waiting_root_is_nominal_not_certified=True,
            unknown_paper_smallness_constants_verified=False,
            original_parameter_errors_enclosed=False,old_core_and_inverse_reused_as_new_source=False,
            new_core_generated=False,terminal_moment_matching_for_new_schedule_certified=False,
            full_admissibility_certified=False,temporal_recursion=False,elapsed_seconds=time.monotonic()-started)
        names=(Path(__file__).name,'lei_ren_part1_paper_outer.py','lei_ren_part1_paper_corrected_profile.py',
            'lei_ren_part1_paper_continuous_angular_schedule.py','lei_ren_part1_paper_schedule_endpoint_enclosures.py',
            'lei_ren_part1_paper_high_order_preheat_integrals.py','lei_ren_part1_paper_variable_preheat_interval_integrals.py',
            'lei_ren_part1_paper_interval_taylor.py','lei_ren_part1_paper_uniform_fixed_beta_error.py',
            'lei_ren_part1_paper_reference_endpoint_targets.py','lei_ren_part1_paper_candidate_pressure_mass_refinement.json')
        result['input_hashes']=dict(old_hashes)
        result['input_hashes'].update({n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
        encoded=encode(_pack(result))
        result['C3_local_pressure']=pressure_jets(c,c.mpf(['.49','.51']),3,encoded)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
        print('Md=1.1 fourteen-stage pressure complete; schedule',digest,'seconds',time.monotonic()-started,flush=True)
        return result

if __name__=='__main__':run()
