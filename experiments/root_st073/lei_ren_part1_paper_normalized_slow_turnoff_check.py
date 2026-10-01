"""Compare normalized propagation bounds against existing actual source data."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from fractions import Fraction
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_normalized_slow_turnoff import certificate,sigma_derivative_upper
from lei_ren_part1_paper_large_Md_screen import weighted_cutoff_integrals
from lei_ren_part1_paper_outer_tail_pressure_probe import tail_pressure
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    name='lei_ren_part1_paper_large_Md_screen.json';raw=json.loads((HERE/name).read_bytes())
    rows=[];c=MPIntervalContext();c.dps=300
    with mp.workdps(340):
        for trial in raw['trials']:
            bound=certificate(trial['Md'],center=('.5','.5'),phase_cells=16)
            z=IntervalTaylor(c,[c.mpf('.5'),c.mpf(1)]);inputs=trial['new_schedule_inputs']
            A=(1+z*z).reciprocal()*c.exp(c.mpf(inputs['logPstar']))
            P0=IntervalTaylor(c,pressure_jets(c,z[0],1,trial['fourteen_stage_trial_pressure_datum'])['physical_pressure_coefficients'])
            inlet=evaluate_transition(c,z,c.mpf(inputs['delta']),c.mpf(1),A,P0,1,512)
            inlet['physical_moments']['p']=A*A*(c.mpf('2.5')+read_interval(c,trial['pressure_stages']['slope_transition_ref']))
            inlet.update(pressure_schedule_Md=trial['Md'])
            field=IntervalOuterAxialTurnoff(c,z,c.mpf(inputs['delta']),inlet,trial['Md'],trial['Md'])
            for phase in ('.25','.5','.75'):
                q=c.mpf(phase);md=c.mpf(trial['Md']);y=c.exp(md*q)
                sig,ds=sigma_value_derivative(c,q)
                p=field._packet(y,y-1,1-sig,-ds/(md*y),weighted_cutoff_integrals(c,trial['Md'],phase,512),Fraction(phase))
                P=tail_pressure(c,z,A,p['y'],trial)
                stress=stress_values(c,z,c.mpf(inputs['delta']),p['R'],p['Utheta'],p['Utheta_y'],p['Uz'],p['Uz_y'],p['physical_moments'],P)
                u=p['Utheta'][0];R=p['R'];L=1-c.mpf(inputs['delta'])*z[0]**2
                Q=stress['stress']['I_theta']*L/(p['F'][0]*R)
                J=stress['stress']['I_z']*L/(c.sqrt(R/2)*u*u)
                Pi=(2*(1+c.mpf(inputs['delta']))*z[0]*P[0]-(1-z[0]**2)*P[1])/(u*u)
                assert endpoints(Q)[0]>=endpoints(bound['Q_uniform_lower'])[0]
                # Older pressure atom boxes may extend beyond the sharper
                # analytic envelope. Require consistency of the enclosures;
                # do not pretend their full widths are exact values.
                pib=endpoints(bound['pressure_Pi_absolute_upper'])[1]
                assert max(endpoints(Pi)[0],-pib)<=min(endpoints(Pi)[1],mp.mpf(0))
                invy=1/p['y']
                jb=read_interval(c,encode(bound['J_radial_derivative_absolute_upper']))
                j0=read_interval(c,encode(bound['J_inlet_absolute_upper']))
                jy_bound=jb+(j0-jb)*invy
                assert max(abs(v) for v in endpoints(J/p['y']))<=endpoints(jy_bound)[1]
                rows.append(dict(Md=trial['Md'],phase=phase,Q=Q,J_over_y=J/p['y'],Pi=Pi,
                    Q_and_J_bounds_contain_actual_source_boxes=True,
                    pressure_box_overlaps_sharper_analytic_envelope=True))
        # Evaluate an independent scalar derivative formula inside all cells,
        # including the two flat endpoint cells. These are smoke checks;
        # the monotonic-product proof supplies endpoint inclusion.
        count=0
        for i in range(512):
            left=c.mpf(i)/512;right=c.mpf(i+1)/512;bound=sigma_derivative_upper(c,left,right)
            for offset in ('.1','.5','.9'):
                q=mp.mpf(i)/512+mp.mpf(offset)/512
                phase=1/q**2-1/(1-q)**2
                sig=1/(1+mp.exp(phase))
                derivative=sig*(1-sig)*(2/q**3+2/(1-q)**3)
                assert derivative<=endpoints(bound)[1]
                count+=1
    result=dict(actual_source_comparisons=rows,independent_derivative_point_checks=count,
        all_passed=True,point_checks_are_not_uniform_proof=True,
        input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
            (name,Path(__file__).name,'lei_ren_part1_paper_normalized_slow_turnoff.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Normalized common-source and derivative checks PASS',flush=True)
    return result


if __name__=='__main__':run()
