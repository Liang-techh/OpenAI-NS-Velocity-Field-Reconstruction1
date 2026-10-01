"""Full-cell relaxed-cone atlas for the repaired reference annulus.

Uses fixed-parameter C1 implicit controls and directed cumulative primitive
envelopes. A cell certificate encloses every radial and axial coordinate
in that cell, rather than inferring coverage from point samples.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_repaired_reference_field import IntervalRepairedReferenceField,stress_values
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def radial_box(c,left,right):
    a=Fraction(left);b=Fraction(right)
    if not 1<=a<=b<=2:raise ValueError('cell must lie within reference annulus [1,2]')
    lo=c.mpf(a.numerator)/a.denominator;hi=c.mpf(b.numerator)/b.denominator
    return c.mpf([endpoints(lo)[0],endpoints(hi)[1]])

def evaluate_cell(field,cell_map,left,right):
    """Enclose the same physical field as the point adapter over a full cell."""
    c=field.calc.ctx;z=field.z;delta=field.calc.delta
    with mp.workdps(field.calc.precision+60):
        xx=radial_box(c,left,right);Rm=field.Rm;Am=field.Am
        h=field.h;d=field.d
        rows=cell_map.apply(h,(Am*Am).reciprocal(),left,right)
        bumps=cell_map.bump_values_on(left,right);beta=bumps['values'];beta_x=bumps['derivatives']
        zero=h[0]*0
        f=sum((h[i+2]*beta[i] for i in range(3)),zero)
        fx=sum((h[i+2]*beta_x[i] for i in range(3)),zero)
        g=h[0]*beta[0]+h[1]*beta[2];gx=h[0]*beta_x[0]+h[1]*beta_x[2]
        power=xx**c.mpf('.1');R=Rm*xx;root=c.sqrt(2*R)
        u=Am*(f+power);uy=Am*(fx*xx+power*c.mpf('.1'))
        uz=z*4+g;uzy=gx*xx;F=u/root
        ref=dict(z=z*(4*R),theta=Am*(c.mpf(5)/8*c.sqrt(2)*R**c.mpf('1.5')*power))
        ref['theta_z']=z*ref['theta']*4
        ref['z_theta']=z*z*(16*R)-Am*Am*(c.mpf(5)/12*R*power*power)
        ref['p']=Am*Am*(c.mpf('2.5')*power*power)
        centered=[d[i]+rows[i] for i in range(5)]
        scale=Am*(c.sqrt(2)*Rm**c.mpf('1.5'))
        inc=dict(z=centered[0]*Rm,theta=centered[2]*scale,
            theta_z=(centered[1]+z*centered[2]*4)*scale,
            z_theta=centered[3]*Am*Am*Rm+z*centered[0]*(8*Rm),p=centered[4]*Am*Am)
        moments={key:ref[key]+inc[key] for key in ref};P=field.P0+moments['p']
        zv=z[0];L=1-delta*zv*zv;mz=moments['z']
        Ur=(2*zv*R*uz[0]-(1-delta)*zv*mz[0]-(1-zv*zv)*mz[1])/(L*root)
        Ur_R=((1+delta)*zv*uz[0]+2*zv*uzy[0]-(1-zv*zv)*uz[1])/(L*root)-Ur/(2*R)
        stress=stress_values(c,z,delta,R,u,uy,uz,uzy,moments,P)
        return dict(left=str(Fraction(left)),right=str(Fraction(right)),x=xx,R=R,F=F,
            Utheta=u,Uz=uz,Utheta_y=uy,Uz_y=uzy,P=P,Ur_value=Ur,Ur_R_value=Ur_R,
            physical_moments=moments,partial_correction_rows=rows,
            radial_bump_envelopes=bumps,**stress)

def partition(subdivisions=4):
    """Include every support edge and centre, plus all outside gaps."""
    if subdivisions<1:raise ValueError('subdivisions must be positive')
    knots=[Fraction(1),Fraction(49,40),Fraction(5,4),Fraction(51,40),
        Fraction(59,40),Fraction(3,2),Fraction(61,40),Fraction(69,40),
        Fraction(7,4),Fraction(71,40),Fraction(2)]
    grid=[]
    for a,b in zip(knots,knots[1:]):
        grid.extend(a+(b-a)*j/subdivisions for j in range(subdivisions))
    grid.append(Fraction(2))
    return list(zip(grid,grid[1:]))

def run():
    from lei_ren_part1_paper_interval_five_bump_cell_map import IntervalFiveBumpCellMap
    field=IntervalRepairedReferenceField();cell_map=IntervalFiveBumpCellMap(field.map)
    cells=partition();records=[]
    for index,(left,right) in enumerate(cells):
        data=evaluate_cell(field,cell_map,left,right)
        # Preserve signs/margins and whole-cell field bounds without repeating
        # the large source moment packet forty times.
        records.append({k:data[k] for k in ('left','right','x','Utheta','Uz','Utheta_y','Uz_y',
            'normalized_stress','cone')})
        if (index+1)%10==0:print('Full radial cells processed:',index+1,'/',len(cells),flush=True)
    coverage=(cells[0][0]==1 and cells[-1][1]==2 and all(a[1]==b[0] for a,b in zip(cells,cells[1:])))
    all_relaxed=all(r['cone'].get('relaxed_cone_certified',False) for r in records)
    unresolved=[dict(left=r['left'],right=r['right'],status=r['cone']['status'])
        for r in records if not r['cone'].get('relaxed_cone_certified',False)]
    result=dict(center_family=field.calc.center_family,state_sha256=field.calc.state_hash,
        radial_domain=['1','2'],coverage_complete=coverage,cell_count=len(cells),
        entire_reference_annulus_relaxed_cone_certified=coverage and all_relaxed,
        whole_cells_not_point_sampling=True,controls_projected_to_midpoints=False,
        same_pressure_datum=True,directed_cumulative_primitive_envelopes=True,
        records=records,unresolved_cells=unresolved,whole_axis=False,strong_admissibility_claimed=False,
        R110_to_Rm_connecting_field_installed=False,heat_exterior_matched=False,
        temporal_recursion=False,original_parameter_remainders_enclosed=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_five_bump_cell_map.py',
        'lei_ren_part1_paper_interval_partial_five_bump_map.py',
        'lei_ren_part1_paper_interval_repaired_reference_field.py',field.inverse_name,field.defects_name,
        'lei_ren_part1_paper_bump_integral_enclosures_check.json','lei_ren_part1_paper_interval_exit_stress_enclosure.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Whole reference annulus relaxed cone:',result['entire_reference_annulus_relaxed_cone_certified'],
        '; unresolved cells:',len(unresolved),flush=True)
    return result

if __name__=='__main__':run()
