"""Resolved-family and actual-source checks of the joint component exit ODE."""
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison
from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
from lei_ren_part1_paper_component_pressure_core import ComponentPressureCore,PressurePolynomial
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_core_adapter import CorePolynomial,build_source_core
from lei_ren_part1_paper_exit_comparison import Section923Comparison
from lei_ren_part1_paper_exit_bridge import ExitBridge
from lei_ren_part1_paper_axial_correction import signed_log


def resolved_fixture():
    with mp.workdps(100):
        f=[mp.mpf(2),mp.mpf('.3')]+[mp.mpf(0)]*3
        u=[mp.mpf('.4'),mp.mpf('.2')]+[mp.mpf(0)]*3
        p=[mp.mpf(-1),mp.mpf('.2')]+[mp.mpf(0)]*3
        tail=[mp.mpf('1e-6'),mp.mpf('3e-6'),mp.mpf('-2e-6')]+[mp.mpf(0)]*2
        kw=dict(F0_Z_taylor=f,U0_Z_taylor=u,radial_degree=3,precision=100)
        formal=core_coefficients('.3','.01',P0_Z_taylor=[PressurePolynomial({0:a,1:b}) for a,b in zip(p,tail)],scalar_converter=PressurePolynomial,**kw)
        axis=SimpleNamespace(precision=100,Lambda=mp.mpf(10),delta=mp.mpf('.01'))
        core=ComponentPressureCore(axis,SimpleNamespace(precision=100),3)
        core.coefficients=lambda Z:formal
        base=dict(axis=axis,Lambda=axis.Lambda,precision=100,radial_degree=3,component_pressure_core=core)
        hb=mp.mpf('1e-6')
        comparison=PressureWidthComparison(base,h_b=hb,pressure_order=9,width_order=2,transition_steps=4)
        bridge=PressureWidthExitBridge(comparison,steps=8)
        maximum=mp.mpf(0)
        for parameter in (0,1):
            scalar=core_coefficients('.3','.01',P0_Z_taylor=[a+parameter*b for a,b in zip(p,tail)],**kw)
            scalar_core=CorePolynomial(lambda Z:scalar,Lambda=10,delta='.01',precision=100)
            scalar_comparison=Section923Comparison(dict(base,core=scalar_core),h_b=hb,transition_steps=4)
            scalar_bridge=ExitBridge(scalar_comparison,epsilon=hb,steps=8)
            for s in (0,2):
                a=bridge.evaluate(s,'.3');b=scalar_bridge.evaluate(s*hb,'.3')
                for name in ('F','Uz','P','log_F_over_Fa','Uz_R','F_R_over_F'):
                    error=abs(a[name].evaluate(pressure=parameter,width=1)-b[name])/max(1,abs(b[name]))
                    maximum=max(maximum,error)
                for name,jet in a['moments'].items():
                    error=abs(jet.evaluate(pressure=parameter,width=1)-b['moments'][name])/max(1,abs(b['moments'][name]))
                    maximum=max(maximum,error)
        assert maximum<mp.mpf('1e-14'),maximum
        return dict(maximum_scaled_difference=mp.nstr(maximum,40),pressure_order=9,width_order=2,
            pressure_parameter_remainder_enclosed=False,width_parameter_remainder_enclosed=False)


def run(pressure_order=9):
    print('resolved width/pressure fixture',flush=True)
    fixture=resolved_fixture();print(fixture,flush=True)
    print('building actual source width-aware exit bridge',flush=True)
    bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
        logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
        coherent_waiting=True,complete_preheat_components=True)
    with mp.workdps(bundle['precision']):
        hb=mp.exp(-100-100*mp.mpf('1e152'))
        comparison=PressureWidthComparison(bundle,h_b=hb,pressure_order=pressure_order,width_order=2,transition_steps=4)
        bridge=PressureWidthExitBridge(comparison,steps=8)
        start=bridge.evaluate(0,'.3');print('actual joint bridge integrating to s=2',flush=True)
        end=bridge.evaluate(2,'.3')
        assert end['log_F_over_Fa'].component(0,1)!=0
        assert end['Uz'].component(1,0)!=0
        assert end['P'].component(1,0)!=0
        increments={k:end['moments'][k]-start['moments'][k] for k in end['moments']}
        assert all(j.component(0,1)!=0 for j in increments.values())
        encode=lambda jet:{str(k):signed_log(v,70) for k,v in jet.atoms.items()}
        report=dict(Z='.3',s_end=2,actual_h_b=signed_log(hb,80),resolved_fixture=fixture,
            metadata=bridge.metadata(),comparison_metadata=comparison.metadata(),
            field_components={k:encode(end[k]) for k in ('F','Uz','P','log_F_over_Fa','F_R_over_F','Uz_R')},
            moment_increments={k:encode(v) for k,v in increments.items()},
            pressure_tail_and_width_atoms_retained=True,all_five_first_width_moment_increments_nonzero=True,
            pressure_parameter_remainder_enclosed=False,width_parameter_remainder_enclosed=False,
            ODE_error_enclosed=False,Z_tangents_installed=False,global_matching_installed=False,
            functional_terminal_moments_closed=False,stress_cone_certified=False,
            finite_energy_certified=False,temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('actual width-aware bridge completed with nonzero five moment increments',flush=True)
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--fixture-only',action='store_true')
    if parser.parse_args().fixture_only:print(resolved_fixture())
    else:run()
