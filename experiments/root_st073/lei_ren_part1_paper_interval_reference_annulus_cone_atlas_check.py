"""Point reduction and enclosing-cell regression for the physical adapter.

Independent primitive/bump fixtures live in the cell-map receipt; independent
physical quadrature lives in the preceding repaired-field fixture. This check
detects drift while translating that field to a full radial box.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_repaired_reference_field import IntervalRepairedReferenceField
from lei_ren_part1_paper_interval_five_bump_cell_map import IntervalFiveBumpCellMap
from lei_ren_part1_paper_interval_reference_annulus_cone_atlas import evaluate_cell,partition
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent

def run():
    field=IntervalRepairedReferenceField();cells=IntervalFiveBumpCellMap(field.map)
    count=0;contained=0
    with mp.workdps(field.calc.precision+60):
        point=field.evaluate_x('1.5');degenerate=evaluate_cell(field,cells,'1.5','1.5')
        enclosing=evaluate_cell(field,cells,'1.475','1.525')
        values=[]
        for key in ('F','Utheta','Uz','Utheta_y','Uz_y','P'):
            values.extend((key,k,point[key][k],degenerate[key][k],enclosing[key][k]) for k in (0,1))
        for key in ('z','theta','theta_z','z_theta','p'):
            values.extend((key,k,point['physical_moments'][key][k],degenerate['physical_moments'][key][k],enclosing['physical_moments'][key][k]) for k in (0,1))
        for key in ('I_theta','I_z','S_theta','S_z','T_theta','T_z'):
            values.append((key,0,point['stress'][key],degenerate['stress'][key],enclosing['stress'][key]))
        for key in ('Ur_value','Ur_R_value'):
            values.append((key,0,point[key],degenerate[key],enclosing[key]))
        for key,k,a,b,outer in values:
            al,ah=endpoints(a);bl,bh=endpoints(b);lo,hi=endpoints(outer)
            if not(al==bl and ah==bh):raise AssertionError(('degenerate physical adapter differs',key,k))
            count+=1
            if not(lo<=al<=ah<=hi):raise AssertionError(('cell excludes point enclosure',key,k))
            contained+=1
        grid=partition()
        if not(grid[0][0]==1 and grid[-1][1]==2 and all(a[1]==b[0] for a,b in zip(grid,grid[1:]))):
            raise AssertionError('partition coverage gap')
    report=dict(point_reduction_exact_coefficients=count,point_enclosures_contained_in_full_cell=contained,
        partition_gap_free=True,independent_quadrature_repeated=False,
        full_cell_certificate_from_interval_bounds_not_this_regression=True)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_reference_annulus_cone_atlas.py',
        'lei_ren_part1_paper_interval_five_bump_cell_map.py','lei_ren_part1_paper_interval_repaired_reference_field.py')
    report['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Physical cell adapter checks:',count,'exact point reductions;',contained,'cell containments',flush=True)
    return report

if __name__=='__main__':run()
