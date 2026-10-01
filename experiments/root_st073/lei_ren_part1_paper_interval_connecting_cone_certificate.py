"""Compose gap-free local relaxed-cone certificates from R110 through Rm."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from lei_ren_part1_paper_interval_long_reshape_field import T_EXACT,RZ_LOG_EXACT,RM_LOG_EXACT

HERE=Path(__file__).parent

def load(name):
    raw=json.loads((HERE/name).read_text())
    for path,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('certificate dependency changed: '+path)
    return raw

def run():
    angular_name='lei_ren_part1_paper_interval_long_reshape_angular_cone.json'
    restore_name='lei_ren_part1_paper_interval_axial_restore_cone_atlas.json'
    angular=load(angular_name);restore=load(restore_name)
    if angular['state_sha256']!=restore['state_sha256'] or angular['center_family']!=restore['center_family']:
        raise ValueError('connecting certificates do not share the same core/axial family')
    shape=[(Fraction(r['left']),Fraction(r['right'])) for r in angular['cells']]
    restoration=[(Fraction(r['left']),Fraction(r['right'])) for r in restore['records']]
    coverage=(shape==[(Fraction(0),T_EXACT),(T_EXACT,RZ_LOG_EXACT)]
        and restoration[0][0]==0 and restoration[-1][1]==2
        and all(a[1]==b[0] for a,b in zip(restoration,restoration[1:]))
        and RZ_LOG_EXACT+2==RM_LOG_EXACT)
    checks=angular['cells']+restore['records']
    complete=coverage and all(r['cone'].get('relaxed_cone_certified',False)
        and r['cone']['branch']=='kappa<=2' for r in checks)
    result=dict(center_family=angular['center_family'],state_sha256=angular['state_sha256'],
        profile_domain='R110 <= R <= Rm',log_radius_domain=['0',str(RM_LOG_EXACT)],
        gap_free_piecewise_coverage=coverage,whole_cell_count=len(checks),
        full_R110_to_Rm_relaxed_cone_certified=complete,
        fixed_parameter_source_family_only=True,strong_admissibility_claimed=False,
        higher_smoothness_or_full_terminal_matching_claimed=False,whole_axis=False,
        original_parameter_remainders_enclosed=False,heat_exterior_matched=False,
        temporal_recursion=False)
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
        for n in (Path(__file__).name,angular_name,restore_name,'lei_ren_part1_paper_interval_long_reshape_field.py')}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Composed R110-to-Rm relaxed cone:',complete,';',len(checks),'whole cells',flush=True)
    return result

if __name__=='__main__':run()
