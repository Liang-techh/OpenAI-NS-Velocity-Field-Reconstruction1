"""Regression from known restoration point packets to full radial cells.

Analytic positive primitive bounds prove cell coverage; this comparison checks
that the new physical cell adapter preserves the existing fields and moments.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_axial_restore_field import IntervalAxialRestoreField
from lei_ren_part1_paper_interval_axial_restore_cone_atlas import evaluate_cell
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent

def run():
    field=IntervalAxialRestoreField();count=0
    with mp.workdps(field.calc.precision+60):
        for phase in (Fraction(0),Fraction(1,2),Fraction(1),Fraction(2)):
            point=field.evaluate_phase(phase)
            left=phase if phase<2 else phase-Fraction(1,32)
            right=left+Fraction(1,32)
            cell=evaluate_cell(field,left,right)
            values=[]
            for key in ('Utheta','Uz','Uz_y','P'):
                values.extend((key,k,point[key][k],cell[key][k]) for k in (0,1))
            for key in ('z','theta','theta_z','z_theta','p'):
                values.extend((key,k,point['physical_moments'][key][k],cell['physical_moments'][key][k]) for k in (0,1))
            for key in ('S_theta_over_F','S_z_over_F','T_theta_over_F','T_z_over_F'):
                values.append((key,0,point['normalized_stress'][key],cell['normalized_stress'][key]))
            for key,k,a,b in values:
                al,ah=endpoints(a);bl,bh=endpoints(b)
                if not bl<=al<=ah<=bh:raise AssertionError(('cell excludes existing point box',str(phase),key,k))
                count+=1
    result=dict(existing_point_packet_coefficients_contained=count,phase_samples=['0','1/2','1','2'],
        production_cell_proof_uses_analytic_envelopes=True,independent_quadrature_repeated=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_axial_restore_cone_atlas.py',
        'lei_ren_part1_paper_interval_axial_restore_field.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Restoration atlas regression passed:',count,'point-packet coefficients contained',flush=True)
    return result

if __name__=='__main__':run()
