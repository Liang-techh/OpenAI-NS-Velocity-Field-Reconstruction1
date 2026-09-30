"""Same-source actual reference background joined with a finite bump response."""
import json
import pickle
from types import SimpleNamespace
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_switches_check import build_provider
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response
from lei_ren_part1_paper_five_bump_field import evaluate_correction,join
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'reference_defect_Z03.pkl'
    derived_cache_used=bool(resume and cache.exists())
    if derived_cache_used:
        # Only this command's internally written local cache is supported.
        with cache.open('rb') as stream:snapshot=pickle.load(stream)
        if snapshot['version']!=1 or snapshot['Z']!='.3' or snapshot['x']!='1.25':raise ValueError('Unsupported local reference cache')
        baseline=snapshot['baseline'];data=snapshot['data']
        defects=SimpleNamespace(work_precision=snapshot['precision'],Rm=snapshot['Rm'])
        reference=SimpleNamespace(delta=snapshot['delta'])
        print('resuming locally cached actual reference and defect inputs',flush=True)
    else:
        source=build_provider(resume=resume)
        switches=PressureWidthExitSwitches(source,steps=8)
        defects=CenteredComponentDefects(switches,A='1e150',logC='5e151',logPstar='14',order=32,window=24,restore_order=64)
        reference=ReferenceDefectBackground(defects,delta=switches.comparison.delta)
        print('actual common-source reference background and bump join',flush=True)
        with mp.workdps(defects.work_precision):
            baseline=reference.evaluate_x('1.25','.3');data=reference.defect_data('.3')
        cache.parent.mkdir(parents=True,exist_ok=True)
        snapshot=dict(version=1,Z='.3',x='1.25',baseline=baseline,data=data,
                      precision=defects.work_precision,Rm=defects.Rm,delta=reference.delta)
        with cache.open('wb') as stream:pickle.dump(snapshot,stream)
        print('local actual reference source snapshot saved',flush=True)
    with mp.workdps(defects.work_precision):
        amplitude=AxialDual(data['Am'],data['Am_Z'])
        moment_map=FiveBumpMomentMap(precision=defects.work_precision,order=96)
        response=build_response(moment_map,amplitude,degree=3)
        directions=[data['d_dual'][i] for i in range(1,6)]
        terms=response.evaluate_terms(directions)
        h=[sum(rows[i] for rows in terms.values()) for i in range(5)]
        correction=evaluate_correction(moment_map,h,amplitude,defects.Rm,'.3','1.25',delta=reference.delta)
        result=join(baseline,correction,delta=reference.delta)
        assert result['P0']==baseline['P0'] and result['P0_Z']==baseline['P0_Z']
        assert result['stress'] is not None,result.get('stress_error')
        assert 'five_bump_correction' in result['moment_parts']
        assert len(baseline['moment_parts'])>=69
        def encode(v):
            if hasattr(v,'value'):v=v.value
            return {str(k):signed_log(a,120) for k,a in v.atoms.items()}
        report=dict(Z='.3',x='1.25',local_R100_source_cache_used=resume,derived_reference_cache_used=derived_cache_used,
          source_scope='Same live common R110 source and complete centered defect integrals; cached R100 mode restricted to Z=.3',
          baseline_fields={k:encode(baseline[k]) for k in ['F','Uz','P','P0','Ur']},
          corrected_fields={k:encode(result[k]) for k in ['F','Uz','P','P0','Ur']},
          baseline_moments={k:encode(v) for k,v in baseline['moments'].items()},
          corrected_moments={k:encode(v) for k,v in result['moments'].items()},
          P0_preserved=True,source_part_count=baseline.get('source_part_count',len(baseline['moment_parts'])),
          baseline_part_count=len(baseline['moment_parts']),
          corrected_part_labels=list(result['moment_parts']),
          finite_defect_response_degree=3,coefficient_monomial_count=len(terms),
          correction_terms_separately_available=True,public_materialization_can_lose_tiny_terms=True,
          functional_closure=False,uniform_Z_certified=False,cone_certified=False,
          exact_heat_certified=False,source_error_enclosed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('same-source actual reference/bump join saved',flush=True)
        return report

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
