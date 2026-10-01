"""Restore O.2 pressure from the same remaining preheat masses.

This is an algebraic evaluation of P0+Mp, not an added/fitted pressure tail.
The original global atom is split at y; all post-y atoms are retained.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def tail_pressure(c,z,A,y,trial):
    schedule=PaperOuterSchedule(**trial['new_schedule_inputs'])
    yd=c.mpf(str(schedule.y_d))
    stages=trial['pressure_stages'];datum=trial['fourteen_stage_trial_pressure_datum']
    remaining=sum((read_interval(c,stages[k]) for k in ('slope_transition_mu','power_buffer','pulse_reserved')),c.mpf(0))
    remaining+=c.exp(c.mpf('.6'))*(c.exp(-y)-c.exp(-yd))/2
    m0=read_interval(c,datum['fixed_beta0_mass_upper_normalized'])
    flatten=read_interval(c,datum['flatten_true_mass_upper_normalized'])
    scale=read_interval(c,datum['Pstar_squared'])
    rho=c.mpf('.25');bound=flatten/((1-rho)**4*rho)
    remainder=IntervalTaylor(c,[flatten,c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])])
    return -A*A*remaining-(remainder+m0)*scale


def evaluate_trial(trial,center=('.49','.51'),cells=256):
    c=MPIntervalContext();c.dps=300
    with mp.workdps(340):
        z=IntervalTaylor(c,[c.mpf(list(center)),c.mpf(1)])
        inputs=trial['new_schedule_inputs'];A=(1+z*z).reciprocal()*c.exp(c.mpf(inputs['logPstar']))
        datum=trial['fourteen_stage_trial_pressure_datum']
        P0=IntervalTaylor(c,pressure_jets(c,z[0],1,datum)['physical_pressure_coefficients'])
        inlet=evaluate_transition(c,z,c.mpf(inputs['delta']),c.mpf(1),A,P0,1,cells)
        # Replace the coarse pressure increment by its existing directed
        # high-order source atom. It is the same defining integral.
        inlet['physical_moments']['p']=A*A*(c.mpf('2.5')+
            read_interval(c,trial['pressure_stages']['slope_transition_ref']))
        inlet.update(pressure_schedule_Md=trial['Md'])
        field=IntervalOuterAxialTurnoff(c,z,c.mpf(inputs['delta']),inlet,trial['Md'],trial['Md'])
        rows=[]
        for phase in ('.25','.5','.75'):
            p=field.evaluate_phase(phase,cells);P=tail_pressure(c,z,A,p['y'],trial)
            for k in (0,1):
                a,b=endpoints(P[k]);d,e=endpoints(p['P'][k])
                if max(a,d)>min(b,e):raise ValueError('same-source pressure evaluations disagree')
            stress=stress_values(c,z,c.mpf(inputs['delta']),p['R'],p['Utheta'],p['Utheta_y'],
                p['Uz'],p['Uz_y'],p['physical_moments'],P)
            F=p['F'][0];R=p['R'];st=c.mpf(-2);sz=2*p['Uz_y'][0]/p['Utheta'][0]
            nt=stress['stress']['I_theta']/(F*R);nz=stress['stress']['I_z']/(F*R)
            con=cone(c,st,sz,nt,nz)
            rows.append(dict(phase=phase,tail_pressure=P,normalized_inertial_stress=dict(theta=nt,z=nz),
                normalized_inertial_cone=con,same_source_pressure_overlap=True))
        print('Same-tail Md',trial['Md'],'axial family',center,'cones',[r['normalized_inertial_cone']['status'] for r in rows],flush=True)
        return dict(Md=trial['Md'],new_schedule_sha256=trial['new_schedule_sha256'],center_family=list(center),
            samples=rows,pressure_tail_is_same_preheat_datum=True,pressure_tail_fitted=False,
            sampled_all_inertial_cones_certified=all(r['normalized_inertial_cone']['admissible_cone_certified'] for r in rows),
            full_phase_or_finite_placement_cone_certified=False,original_parameter_errors_enclosed=False)


def run():
    name='lei_ren_part1_paper_outer_parameter_probe.json';raw=json.loads((HERE/name).read_bytes())
    results=[]
    for trial in raw['trials']:
        if trial['Md'] not in ('2','3'):continue
        results.append(evaluate_trial(trial))
        results.append(evaluate_trial(trial,center=('.5','.5')))
    result=dict(trials=results,stage='same-source analytic pressure restoration for O.2',
        whole_outer_admissibility_certified=False,temporal_recursion=False,
        input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
            (name,Path(__file__).name,'lei_ren_part1_paper_interval_outer_axial_turnoff.py',
             'lei_ren_part1_paper_interval_outer_slope_field.py','lei_ren_part1_paper_Md11_pressure_datum.py',
             'lei_ren_part1_paper_interval_repaired_reference_field.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':run()
