"""Actual component chain through axial restoration and reference continuation."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_switches_check import build_provider
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_pressure_width_long_reshape import PressureWidthLongReshape
from lei_ren_part1_paper_pressure_width_reference_extension import PressureWidthReferenceExtension
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    provider=build_provider(resume=resume)
    switches=PressureWidthExitSwitches(provider,steps=8)
    reshape=PressureWidthLongReshape(switches,A='1e150',logC='5e151',order=24,tail_digits=60)
    from lei_ren_part1_paper_pressure_width_axial_restore import PressureWidthAxialRestore
    with mp.workdps(reshape.work_precision):
        reference=PressureWidthReferenceExtension(reshape,logRref=mp.log(110)+10*(mp.mpf('5e151')+14))
        restore=PressureWidthAxialRestore(reference,order=48)
        print('computing common actual axial restoration entry',flush=True)
        initial=reference.evaluate_log_offset(reference.extension_length,'.3')
        encode=lambda v:{str(k):signed_log(a,60) for k,a in v.atoms.items()}
        error=lambda a,b:max((abs(a.component(*k)-b.component(*k))/max(1,abs(b.component(*k)))
            for k in set(a.atoms)|set(b.atoms)),default=mp.mpf(0))
        start=restore.evaluate_phase(0,'.3')
        boundary={k:error(start[k],initial[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur')}
        assert max(boundary.values())<mp.mpf('1e-200')
        rows=[]
        for phase in ('.5','1','3'):
            print('actual axial restoration phase='+phase,flush=True)
            v=restore.evaluate_phase(phase,'.3')
            assert error(v['P0'],initial['P0'])==0
            for part in ('inner_seed','reshape','reference'):
                assert all(error(v['moment_parts'][part][k],initial['moment_parts'][part][k])==0 for k in initial['moments'])
            if phase!='.5':
                assert error(v['Uz'],reference._base_jet('1.2'))<mp.mpf('1e-200')
                assert error(v['Uz_Z'],reference._base_jet(4))<mp.mpf('1e-200')
            rows.append(dict(phase=phase,fields={k:encode(v[k]) for k in ('F','F_Z','Utheta','Utheta_Z','Uz','Uz_Z','Uz_y','P','P_Z','Ur')},
                moment_parts={part:{k:encode(a) for k,a in values.items()} for part,values in v['moment_parts'].items()},
                moment_parts_Z={part:{k:encode(a) for k,a in values.items()} for part,values in v['moment_parts_Z'].items()}))
        report=dict(Z='.3',metadata=restore.metadata(),local_source_cache_used=resume,
            boundary_scaled_errors={k:mp.nstr(v,40) for k,v in boundary.items()},rows=rows,
            inherited_parts_preserved=True,axis_pressure_preserved=True,
            functional_terminal_moments_closed=False,finite_energy_certified=False,
            cone_certified=False,temporal_recursion_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual axial restoration receipt saved',flush=True)
        return report

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
