"""Actual common R110 source assembled into five centered finite-ring defects."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_switches_check import build_provider
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects
    provider=build_provider(resume=resume)
    switches=PressureWidthExitSwitches(provider,steps=8)
    print('assembling actual five centered defect rows',flush=True)
    defects=CenteredComponentDefects(switches,A='1e150',logC='5e151',logPstar='14',order=32,window=24,restore_order=64)
    result=defects.evaluate('.3')
    with mp.workdps(defects.work_precision):
        def encode(a):
            if not hasattr(a,'atoms') and hasattr(a,'value'):a=a.value
            return {str(k):signed_log(v,120) for k,v in a.atoms.items()}
        d=result['d'];dZ=result['d_Z']
        assert len(d)==5 and len(dZ)==5
        assert all(a.component(0,0)!=0 for a in d.values())
        assert any(a.component(1,0)!=0 for a in d.values())
        assert any(a.component(0,1)!=0 for a in d.values())
        report=dict(Z='.3',metadata=result.get('metadata',{}),local_source_cache_used=resume,
            defects={k:encode(v) for k,v in d.items()},defects_Z={k:encode(v) for k,v in dZ.items()},
            parts={label:{k:encode(v) for k,v in rows.items()} for label,rows in result['parts'].items()},
            parts_Z={label:{k:encode(v) for k,v in rows.items()} for label,rows in result['parts_Z'].items()},
            all_five_baseline_defects_nonzero=True,source_scope='Same actual R110 core/pressure/moments and automatic first Z tangents',
            source_error_enclosed=False,smallness_certified=False,functional_terminal_moments_closed=False,
            five_bump_repair_installed=False,cone_certified=False,temporal_recursion_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,default=lambda v:mp.nstr(v,100) if isinstance(v,mp.mpf) else str(v))+'\n',encoding='utf-8')
        print('actual five centered defect receipt saved',flush=True)
        return report

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
