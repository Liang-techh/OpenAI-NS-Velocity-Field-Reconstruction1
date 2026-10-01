"""Whole-cell cone atlas for axial restoration and continuation Rz..Rm."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_axial_restore_field import IntervalAxialRestoreField
from lei_ren_part1_paper_interval_long_reshape_field import positive_decay_integral,sigma_value_derivative
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def evaluate_cell(field,left,right):
    a,b=Fraction(left),Fraction(right)
    if not 0<=a<=b<=2:raise ValueError('restoration cell outside [0,2]')
    c=field.calc.ctx
    with mp.workdps(field.calc.precision+60):
        lo=c.mpf(a.numerator)/a.denominator;hi=c.mpf(b.numerator)/b.denominator
        t=c.mpf([endpoints(lo)[0],endpoints(hi)[1]])
        if a>=1:cutoff,ds=c.mpf(0),c.mpf(0)
        else:
            cutoff=alpha_box(c,t+1,c.mpf(1))
            _,ds=sigma_value_derivative(c,t)
        # Unnormalized primitives are nonnegative and monotone; their saved
        # full-support upper bounds enclose every partial integral.
        K={}
        for name,rate in (('K1','1'),('K16','1.6'),('K2','1')):
            primitive=field.full[name] if a>=1 else c.mpf([0,endpoints(field.full[name])[1]])
            K[name]=primitive*c.exp(-c.mpf(rate)*t)
        R=field.start['R']*c.exp(t);u=field.start['Utheta']*c.exp(t/10)
        uy=u*c.mpf('.1');V=field.z*4+field.g*cutoff;Vy=-field.g*ds
        theta=u*(c.sqrt(2)*R**c.mpf('1.5'))*positive_decay_integral(c,c.mpf('1.6'),t)
        pressure=u*u/2*positive_decay_integral(c,c.mpf('.2'),t)
        swirl=u*u*(R/2)*positive_decay_integral(c,c.mpf('1.2'),t)
        deltaR=field.start['R']*(c.exp(t)-1)
        mass=field.z*(4*deltaR)+field.g*(R*K['K1'])
        mixed=field.z*theta*4+field.g*u*(c.sqrt(2)*R**c.mpf('1.5')*K['K16'])
        axial=field.z*field.z*(16*deltaR)+field.z*field.g*(8*R*K['K1'])+field.g*field.g*(R*K['K2'])
        m0=field.start['physical_moments']
        moments=dict(z=m0['z']+mass,theta=m0['theta']+theta,theta_z=m0['theta_z']+mixed,
            z_theta=m0['z_theta']+axial-swirl,p=m0['p']+pressure)
        P=field.P0+moments['p'];root=c.sqrt(2*R);F=u/root
        stress=stress_values(c,field.z,field.calc.delta,R,u,uy,V,Vy,moments,P)
        st=-c.mpf('.8');sz=2*Vy[0]/u[0]
        Itheta=stress['stress']['I_theta'];Iz=stress['stress']['I_z']
        normalized=dict(S_theta_over_F=st,S_z_over_F=sz,T_theta_over_F=Itheta/F[0]+st,T_z_over_F=Iz/F[0]+sz)
        certificate=cone(c,st,sz,normalized['T_theta_over_F'],normalized['T_z_over_F'])
        return dict(left=str(a),right=str(b),normalized_stress=normalized,cone=certificate,
            Utheta=u,Uz=V,Uz_y=Vy,P=P,physical_moments=moments,
            full_partial_primitive_envelopes=True,same_axis_pressure=True)

def run():
    field=IntervalAxialRestoreField();cells=[(Fraction(i,32),Fraction(i+1,32)) for i in range(64)]
    records=[]
    for i,(a,b) in enumerate(cells):
        data=evaluate_cell(field,a,b)
        records.append({k:data[k] for k in ('left','right','normalized_stress','cone')})
        if (i+1)%16==0:print('Restoration full cells:',i+1,'/',len(cells),flush=True)
    complete=all(r['cone'].get('relaxed_cone_certified',False) for r in records)
    coverage=cells[0][0]==0 and cells[-1][1]==2 and all(a[1]==b[0] for a,b in zip(cells,cells[1:]))
    unresolved=[dict(left=r['left'],right=r['right'],status=r['cone']['status']) for r in records
        if not r['cone'].get('relaxed_cone_certified',False)]
    result=dict(center_family=field.calc.center_family,state_sha256=field.calc.state_hash,
        phase_domain=['0','2'],coverage_complete=coverage,cell_count=len(cells),records=records,
        entire_Rz_to_Rm_relaxed_cone_certified=complete and coverage,unresolved_cells=unresolved,
        whole_cells_not_point_sampling=True,strong_admissibility_claimed=False,whole_axis=False,
        original_parameter_remainders_enclosed=False,heat_exterior_matched=False,temporal_recursion=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_axial_restore_field.py',
        'lei_ren_part1_paper_interval_long_reshape_field.py','lei_ren_part1_paper_interval_functional_defects.json',
        'lei_ren_part1_paper_interval_repaired_reference_field.py','lei_ren_part1_paper_interval_exit_stress_enclosure.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Whole Rz-to-Rm relaxed cone:',result['entire_Rz_to_Rm_relaxed_cone_certified'],
        '; unresolved cells:',len(unresolved),flush=True)
    return result

if __name__=='__main__':run()
