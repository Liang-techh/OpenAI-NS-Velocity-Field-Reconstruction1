"""Coherent larger-Md preheat source screening before fresh core construction."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_reference_endpoint_targets import accepted_parameters
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_waiting_length import incoming_angular_ratio,solve_waiting_length,apply_waiting_root
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


def build_trial(Md):
    with mp.workdps(340):
        lp=max(mp.mpf(14),mp.exp(mp.mpf(Md))+11)
        delta=min(mp.mpf('1e-200'),mp.exp(-4*lp-30))
        params=dict(accepted_parameters(),Md=Md,delta=mp.nstr(delta,280),logPstar=mp.nstr(lp,280),waiting_length='0',
            logRref=mp.nstr(mp.log(110)+10*(mp.mpf('5e151')+lp),300))
        base=PaperOuterSchedule(**params)
        incoming=incoming_angular_ratio(base)
        s=apply_waiting_root(base,solve_waiting_length(base,incoming['incoming_X']))
        # Pressure atoms need the schedule and waiting root, not the separate
        # finite-float bump coefficient factory (which underflows at Md>=6).
        install_continuous_angular_schedule(s,precision=260,primitive_precision=100)
        s._coherent_preheat_waiting=True
        e=ScheduleEndpointEnclosures(s,precision=300);c=e.iv
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
        result=dict(Md=Md,logPstar=inputs['logPstar'],new_schedule_inputs=inputs,new_schedule_sha256=sha,
            fourteen_stage_trial_pressure_datum=datum,pressure_stages=masses,flatten_mass_upper=flatten,
            delta_below_c_delta_mu=endpoints(c.mpf(inputs['c_delta'])*c.mpf(str(s.mu))-c.mpf(inputs['delta']))[0]>0,
            unknown_paper_constants_verified=False,waiting_root_certified=False,
            original_parameter_errors_enclosed=False,new_core_generated=False)
        return encode(_pack(result))


def weighted_cutoff_integrals(c, Md, phase, cells):
    """Integrate exp(s) exactly; enclose only the monotone cutoff on cells."""
    from fractions import Fraction
    from lei_ren_part1_paper_interval_outer_slope_field import fraction_box
    from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
    md_exact=Fraction(Md);q=Fraction(phase)
    if md_exact<=0 or not 0<=q<=1 or not isinstance(cells,int) or cells<1:
        raise ValueError('positive Md, phase in [0,1], positive integer cells required')
    md=fraction_box(c,md_exact)
    mass=c.mpf(0);square=c.mpf(0)
    for i in range(cells):
        left=fraction_box(c,q*i/cells);right=fraction_box(c,q*(i+1)/cells)
        lo=1-sigma_value_derivative(c,right)[0];hi=1-sigma_value_derivative(c,left)[0]
        cutoff=c.mpf([endpoints(lo)[0],endpoints(hi)[1]])
        weight=c.exp(c.exp(md*right))-c.exp(c.exp(md*left))
        mass+=weight*cutoff;square+=weight*cutoff**2
    return dict(B_mass=mass,B_squared_mass=square)


def screen(trial, center=('.49','.51'), cells=512):
    from mpmath.ctx_iv import MPIntervalContext
    from lei_ren_part1_paper_outer_tail_pressure_probe import tail_pressure
    from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
    c=MPIntervalContext();c.dps=300
    with mp.workdps(340):
        inputs=trial['new_schedule_inputs'];z=IntervalTaylor(c,[c.mpf(list(center)),c.mpf(1)])
        A=(1+z*z).reciprocal()*c.exp(c.mpf(inputs['logPstar']))
        P0=IntervalTaylor(c,pressure_jets(c,z[0],1,trial['fourteen_stage_trial_pressure_datum'])['physical_pressure_coefficients'])
        inlet=evaluate_transition(c,z,c.mpf(inputs['delta']),c.mpf(1),A,P0,1,cells)
        inlet['physical_moments']['p']=A*A*(c.mpf('2.5')+read_interval(c,trial['pressure_stages']['slope_transition_ref']))
        inlet.update(pressure_schedule_Md=trial['Md'])
        field=IntervalOuterAxialTurnoff(c,z,c.mpf(inputs['delta']),inlet,trial['Md'],trial['Md'])
        rows=[]
        for phase in ('.25','.5','.75'):
            from fractions import Fraction
            from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
            from lei_ren_part1_paper_interval_outer_slope_field import fraction_box
            q=fraction_box(c,Fraction(phase));md=fraction_box(c,Fraction(trial['Md']))
            y=c.exp(md*q);sig,ds=sigma_value_derivative(c,q)
            p=field._packet(y,y-1,1-sig,-ds/(md*y),weighted_cutoff_integrals(c,trial['Md'],phase,cells),Fraction(phase))
            P=tail_pressure(c,z,A,p['y'],trial)
            stress=stress_values(c,z,c.mpf(inputs['delta']),p['R'],p['Utheta'],p['Utheta_y'],p['Uz'],p['Uz_y'],p['physical_moments'],P)
            F=p['F'][0];R=p['R'];b=-2*p['Uz_y'][0]/p['Utheta'][0]
            x=stress['stress']['I_theta']/(F*R);nz=stress['stress']['I_z']/(F*R)
            B=b*b;Q=b*nz;D=2*x+Q
            # Exact algebraic cancellation of Q^2. Never obtain kappa-2
            # by subtracting two nearly equal rounded numbers.
            bracket=(8-B*B/2)*x+(8+2*B)*Q
            margin=x*bracket
            certified=endpoints(b)[0]>0 and endpoints(x)[0]>0 and endpoints(bracket)[0]>0 and endpoints(D)[0]>0
            obstructed=endpoints(D)[1]<0 or endpoints(margin)[1]<0
            status='sampled asymptotic strong cone certified' if certified else ('ruled out' if obstructed else 'unresolved')
            rows.append(dict(phase=phase,b_squared=B,Ntheta=x,b_Nz=Q,negative_dot_per_R=D,strong_margin_per_R_squared=margin,
                status=status,sampled_asymptotic_strong_cone_certified=certified,
                asymptotic_direction_ruled_out=endpoints(D)[1]<0,
                direction_obstructed_for_all_positive_placements=endpoints(D)[1]<0,
                asymptotic_margin_ruled_out=endpoints(margin)[1]<0,
                actual_finite_placement_certified=False))
        print('Large Md',trial['Md'],'family',center,[r['status'] for r in rows],flush=True)
        return encode(_pack(dict(center_family=list(center),samples=rows,
            all_sampled_asymptotic_cones_certified=all(r['sampled_asymptotic_strong_cone_certified'] for r in rows))))


def run():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--Md',nargs='+',default=['4','5','6'])
    parser.add_argument('--reuse-trials',action='store_true',help='Reuse stored fourteen-stage masses; rerun all cone samples')
    args=parser.parse_args();out=Path(__file__).with_suffix('.json')
    saved=json.loads(out.read_bytes()) if args.reuse_trials and out.exists() else None
    dependencies=set(json.loads((HERE/'lei_ren_part1_paper_outer_parameter_probe.json').read_bytes())['input_hashes'])
    dependencies.update([Path(__file__).name,'lei_ren_part1_paper_outer_tail_pressure_probe.py',
        'lei_ren_part1_paper_waiting_length.py','lei_ren_part1_paper_interval_long_reshape_field.py',
        'lei_ren_part1_paper_interval_repaired_reference_field.py'])
    result=dict(trials=[],whole_outer_admissibility_certified=False,temporal_recursion=False,
        input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sorted(dependencies)})
    for md in args.Md:
        existing=next((t for t in saved['trials'] if t['Md']==md),None) if saved else None
        trial=existing if existing is not None else build_trial(md)
        if existing is not None:
            trial.setdefault('datum_generation_source_sha256',saved['input_hashes'][Path(__file__).name])
        else:
            trial['datum_generation_source_sha256']=result['input_hashes'][Path(__file__).name]
        if not trial['delta_below_c_delta_mu']:raise ValueError('delta hierarchy violated')
        trial['screens']=[screen(trial),screen(trial,center=('.5','.5'))]
        result['trials'].append(trial)
        out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':run()

