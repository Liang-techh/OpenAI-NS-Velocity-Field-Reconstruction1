"""Common actual-source reference continuation to the axial restoration entry."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_switches_check import build_provider
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_pressure_width_long_reshape import PressureWidthLongReshape
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    provider=build_provider(resume=resume)
    switches=PressureWidthExitSwitches(provider,steps=8)
    reshape=PressureWidthLongReshape(switches,A='1e150',logC='5e151',order=24,tail_digits=60)
    from lei_ren_part1_paper_pressure_width_reference_extension import PressureWidthReferenceExtension
    with mp.workdps(reshape.work_precision):
        logRref=mp.log(110)+10*(mp.mpf('5e151')+14)
        extension=PressureWidthReferenceExtension(reshape,logRref=logRref)
        encode=lambda v:{str(k):signed_log(a,60) for k,a in v.atoms.items()}
        error=lambda a,b:max((abs(a.component(*k)-b.component(*k))/max(1,abs(b.component(*k)))
            for k in set(a.atoms)|set(b.atoms)),default=mp.mpf(0))
        print('computing common actual reshape endpoint',flush=True)
        initial=extension._endpoint(extension._z_key('.3'))
        start=extension.evaluate_log_offset(0,'.3')
        boundary={k:error(start[k],initial[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur')}
        assert max(boundary.values())<mp.mpf('1e-200')
        rows=[]
        for label,x in [('one_log_unit',mp.mpf(1)),('axial_restore_entry',extension.extension_length)]:
            print('reference continuation '+label,flush=True)
            v=extension.evaluate_log_offset(x,'.3')
            assert v['Uz_Z'].component(1,0)!=0
            assert all(error(v['moment_parts']['inner_seed'][k],initial['moment_seeds'][k])==0 for k in initial['moment_seeds'])
            assert all(error(v['moment_parts']['reshape'][k],initial['moment_increments'][k])==0 for k in initial['moment_increments'])
            assert error(v['P0'],initial['P0'])==0
            rows.append(dict(label=label,x=mp.nstr(x,60),logR=mp.nstr(v['logR'],60),
                fields={k:encode(v[k]) for k in ('F','F_Z','Utheta','Utheta_Z','Uz','Uz_Z','P','P_Z','Ur')},
                moment_parts={part:{k:encode(a) for k,a in values.items()} for part,values in v['moment_parts'].items()},
                moment_parts_Z={part:{k:encode(a) for k,a in values.items()} for part,values in v['moment_parts_Z'].items()},
                raw_quadratic_parts={part:{k:encode(a) for k,a in values.items()} for part,values in v['raw_quadratic_parts'].items()}))
        report=dict(Z='.3',metadata=extension.metadata(),local_source_cache_used=resume,
            boundary_scaled_errors={k:mp.nstr(v,40) for k,v in boundary.items()},rows=rows,
            inherited_parts_preserved=True,axis_pressure_preserved=True,
            functional_terminal_moments_closed=False,axial_restoration_installed=False,
            finite_energy_certified=False,cone_certified=False,temporal_recursion_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('reference continuation receipt saved',flush=True)
        return report

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
