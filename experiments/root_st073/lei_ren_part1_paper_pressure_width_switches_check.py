"""Actual-source two-switch field and five-moment propagation to R=110."""
import json
import pickle
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_pressure_width_axial_comparison import AxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_axial_bridge import AxialPressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_continuation import PressureWidthExitContinuation
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False):
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'pressure_width_R100_source.pkl'
    if resume:
        # Only the local cache written by this command is supported. It is not
        # an import format for downloaded or user-supplied pickle files.
        with cache.open('rb') as stream:data=pickle.load(stream)
        if data['version']!=1:raise ValueError('Unsupported R100 cache version')
        axis=SimpleNamespace(precision=260,Lambda=data['Lambda'],delta=data['delta'])
        def coefficients(Z):
            if mp.mpf(str(Z))!=mp.mpf('.3'):raise ValueError('Cached source supports Z=.3 only')
            return data['coefficients']
        component=SimpleNamespace(axis=axis,precision=260,coefficients=coefficients)
        bundle=dict(precision=260,axis=axis,component_pressure_core=component)
        print('resuming locally cached actual R100 source',flush=True)
    else:
        print('building actual source for component shear switches',flush=True)
        bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
            logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
            coherent_waiting=True,complete_preheat_components=True,component_Z_jet_depth=2)
    with mp.workdps(bundle['precision']):
        comparison=AxialPressureWidthComparison(bundle,h_b=mp.exp(-100-100*mp.mpf('1e152')),
            pressure_order=9,width_order=2,transition_steps=4)
        bridge=AxialPressureWidthExitBridge(comparison,steps=8)
        provider=PressureWidthExitContinuation(bridge)
        if resume:
            def cached_functions(Z):
                if mp.mpf(str(Z))!=mp.mpf('.3'):raise ValueError('Cached source supports Z=.3 only')
                return data['functions']
            provider.functions=cached_functions
        print('propagating common collar and R100 data',flush=True)
        start=provider.evaluate_R(100,'.3')
        if not resume:
            cache.parent.mkdir(parents=True,exist_ok=True)
            data=dict(version=1,Lambda=comparison.Lambda,delta=comparison.delta,
                coefficients=bundle['component_pressure_core'].coefficients('.3'),
                functions=provider.functions('.3'))
            with cache.open('wb') as stream:pickle.dump(data,stream)
            print('local R100 source cache saved',flush=True)
        from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
        switches=PressureWidthExitSwitches(provider,steps=8)
        boundary=switches.evaluate_R(100,'.3')
        encode=lambda v:{str(k):signed_log(a,70) for k,a in v.atoms.items()}
        error=lambda a,b:max((abs(a.component(*k)-b.component(*k))/max(1,abs(b.component(*k)))
            for k in set(a.atoms)|set(b.atoms)),default=mp.mpf(0))
        boundary_errors={k:error(boundary[k],start[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur')}
        boundary_errors.update({k:error(boundary['moments'][k],start['moments'][k]) for k in start['moments']})
        assert max(boundary_errors.values())<mp.mpf('1e-200'),boundary_errors
        print('integrating both tiny switches and exact constant-power segment',flush=True)
        end=switches.evaluate_R(110,'.3')
        ratio=end['F']/start['F']
        expected=mp.mpf('1.1')**mp.mpf('-.4')
        assert abs(ratio.component(0,0)-expected)<mp.mpf('1e-200')
        assert end['Uz_Z'].component(1,0)!=0
        assert end['g_y'].component(0,0)==mp.mpf('-.4')
        assert not end['Uz_y'].atoms
        raw=end['raw_quadratic_integrals']
        quadratic_error=error(raw['axial']-raw['swirl'],end['moments']['z_theta'])
        assert quadratic_error<mp.mpf('1e-200'),quadratic_error
        # First-Z local input diagnostics for v2 equations 9.28--9.31.
        # Full C2 bounds and uniform cone conditions are separate requirements.
        z=mp.mpf('.3')
        B=end['F'].log()+mp.log(220)/2+mp.mpf('5e151')+mp.log(1+z*z)
        BZ=end['F_Z']/end['F']+2*z/(1+z*z)
        inputs=dict(D_R110=end['D'],V_minus_4Z=end['Uz']-4*z,V_Z_minus_4=end['Uz_Z']-4,
            reshape_B=B,reshape_B_Z=BZ)
        report=dict(Z='.3',R_start=100,R_end=110,metadata=switches.metadata(),local_source_cache_used=resume,
            boundary_max_scaled_errors={k:mp.nstr(v,40) for k,v in boundary_errors.items()},
            F_ratio_components=encode(ratio),leading_F_ratio=mp.nstr(expected,50),
            fields={k:encode(end[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur','g_y','Uz_y')},
            moments={k:encode(v) for k,v in end['moments'].items()},
            moments_Z={k:encode(v) for k,v in end['moments_Z'].items()},
            quadratic_moment_identity_error=mp.nstr(quadratic_error,40),
            local_reshape_input_components={k:encode(v) for k,v in inputs.items()},
            uniform_C2_source_bounds_verified=False,
            complete_preheat_pressure_tail_retained=True,constant_power_segment_verified=True,
            long_reshape_installed=False,functional_terminal_moments_closed=False,
            finite_energy_certified=False,cone_certified=False,temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('actual component switches R110 receipt saved',flush=True)
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    run(resume=parser.parse_args().resume)
