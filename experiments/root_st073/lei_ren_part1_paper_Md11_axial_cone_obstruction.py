"""Radius-independent direction obstruction in the new Md1.1 turnoff.

On this outer segment I_theta/F=R*Ntheta and I_z/F=R*Nz exactly.
Hence (S/F).(T/F)/R = st*Ntheta+sz*Nz+(st^2+sz^2)/R.
If the first two terms have positive lower bound, the required negative
direction fails for every positive radial placement, not just a huge number.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    datum=json.loads((HERE/'lei_ren_part1_paper_Md11_pressure_datum.json').read_bytes())
    c=MPIntervalContext();c.dps=100
    with mp.workdps(140):
        z=IntervalTaylor(c,[c.mpf(['.49','.51']),c.mpf(1)])
        A=(1+z*z).reciprocal()*c.exp(14)
        P0=IntervalTaylor(c,pressure_jets(c,z[0],1,datum)['physical_pressure_coefficients'])
        rows=[]
        for placement in ('1','10','1e20'):
            inlet=evaluate_transition(c,z,c.mpf('1e-200'),c.mpf(placement),A,P0,1)
            inlet.update(pressure_schedule_Md='1.1',pressure_datum_kind='new analytic preheat datum')
            field=IntervalOuterAxialTurnoff(c,z,c.mpf('1e-200'),inlet,'1.1','1.1')
            p=field.evaluate_phase('.5');F=p['F'][0];R=p['R']
            st=p['normalized_stress']['S_theta_over_F'];sz=p['normalized_stress']['S_z_over_F']
            Ntheta=p['stress']['I_theta']/(F*R);Nz=p['stress']['I_z']/(F*R)
            direction=st*Ntheta+sz*Nz
            if endpoints(direction)[0]<=0:raise ValueError('radius-independent direction obstruction not proved')
            rows.append(dict(reference_radius=placement,inertial_direction_per_R=direction,
                nonnegative_shear_square_term=(st*st+sz*sz)/R,
                total_direction_per_R=direction+(st*st+sz*sz)/R))
        lower=max(endpoints(r['inertial_direction_per_R'])[0] for r in rows)
        upper=min(endpoints(r['inertial_direction_per_R'])[1] for r in rows)
        if lower>upper:raise ValueError('radius-normalized diagnostic inconsistent')
        result=dict(new_schedule_sha256=datum['new_schedule_sha256'],Md='1.1',phase='.5',
            axial_family=['.49','.51'],radius_independent_inertial_direction_positive=True,
            necessary_negative_direction_condition_ruled_out_for_this_candidate=True,
            obstruction_is_not_only_interval_noncertificate=True,
            original_parameter_errors_enclosed=False,whole_axis=False,
            normalized_placement_consistency_overlap=c.mpf([lower,upper]),diagnostics=rows,
            remedy_requires_new_parameter_or_construction_data_not_more_radial_orders=True,
            temporal_recursion=False)
        names=(Path(__file__).name,'lei_ren_part1_paper_Md11_pressure_datum.json',
            'lei_ren_part1_paper_Md11_pressure_datum.py','lei_ren_part1_paper_interval_outer_slope_field.py',
            'lei_ren_part1_paper_interval_outer_axial_turnoff.py')
        result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Md1.1 axial direction obstruction proved over local family; normalized direction interval',
            mp.nstr(lower,16),mp.nstr(upper,16),flush=True)
        return result

if __name__=='__main__':run()
