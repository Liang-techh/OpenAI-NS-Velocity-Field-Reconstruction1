"""Actual common-source reshape with separately retained R110 moment seeds."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_switches_check import build_provider
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    provider=build_provider(resume=resume)
    switches=PressureWidthExitSwitches(provider,steps=8)
    print('propagating common R110 input for long reshape',flush=True)
    initial=switches.evaluate_R(110,'.3')
    from lei_ren_part1_paper_pressure_width_long_reshape import PressureWidthLongReshape
    reshape=PressureWidthLongReshape(switches,A='1e150',logC='5e151',order=24,tail_digits=60)
    with mp.workdps(reshape.work_precision):
        encode=lambda v:{str(k):signed_log(a,70) for k,a in v.atoms.items()}
        error=lambda a,b:max((abs(a.component(*k)-b.component(*k))/max(1,abs(b.component(*k)))
            for k in set(a.atoms)|set(b.atoms)),default=mp.mpf(0))
        start=reshape.evaluate_phase(0,'.3')
        boundary_errors={k:error(start[k],initial[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur')}
        assert max(boundary_errors.values())<mp.mpf('1e-200'),boundary_errors
        rows=[]
        for phase in ('.5','1'):
            print(f'actual normalized reshape phase={phase}',flush=True)
            end=reshape.evaluate_phase(phase,'.3')
            seeds=end['moment_seeds']; increments=end['moment_increments']
            assert all(error(seeds[k],initial['moments'][k])==0 for k in seeds)
            assert seeds['z'].component(1,0)!=0
            assert increments['theta'].component(0,0)!=0
            assert end['Uz_Z'].component(1,0)!=0
            rows.append(dict(phase=phase,log_radius_offset=mp.nstr(end['log_radius_offset'],60),
                fields={k:encode(end[k]) for k in ('F','F_Z','Utheta','Utheta_Z','Uz','Uz_Z','P','P_Z','Ur')},
                moment_seeds={k:encode(v) for k,v in seeds.items()},
                moment_seeds_Z={k:encode(v) for k,v in end['moment_seeds_Z'].items()},
                moment_increments={k:encode(v) for k,v in increments.items()},
                moment_increments_Z={k:encode(v) for k,v in end['moment_increments_Z'].items()},
                normalized_integrals={k:encode(v) for k,v in end['normalized_integrals'].items()},
                normalized_tail_bounds={k:mp.nstr(v,60) for k,v in end['normalized_tail_bounds'].items()},
                seed_components_retained=True))
        z=mp.mpf('.3'); target=-2*z/(1+z*z)
        end_slope=end['Utheta_Z']/end['Utheta']
        slope_error=abs(end_slope.component(0,0)-target)
        assert slope_error<mp.mpf('1e-200'),slope_error
        expected_log=reshape.T/10-reshape.logC-mp.log(1+z*z)
        actual_log=end['Utheta'].log()
        log_error=abs(actual_log.component(0,0)-expected_log)
        assert log_error<mp.mpf('1e-200'),log_error
        report=dict(Z='.3',metadata=reshape.metadata(),local_source_cache_used=resume,
            boundary_scaled_errors={k:mp.nstr(v,40) for k,v in boundary_errors.items()},rows=rows,
            endpoint_reference_log_error=mp.nstr(log_error,40),endpoint_reference_Z_slope_error=mp.nstr(slope_error,40),
            endpoint_seeds_preserved_separately=True,quadrature_error_enclosed=False,
            reference_extension_installed=False,functional_terminal_moments_closed=False,
            finite_energy_certified=False,cone_certified=False,temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('actual reshape receipt saved with separate inherited seeds',flush=True)
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
