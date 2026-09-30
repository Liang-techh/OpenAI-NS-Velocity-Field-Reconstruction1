"""Actual source continuation retaining positive-epsilon width atoms."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_pressure_width_axial_comparison import AxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_axial_bridge import AxialPressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_continuation import PressureWidthExitContinuation
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    print('building actual source for post-collar continuation',flush=True)
    bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
        logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
        coherent_waiting=True,complete_preheat_components=True,component_Z_jet_depth=2)
    with mp.workdps(bundle['precision']):
        comparison=AxialPressureWidthComparison(bundle,h_b=mp.exp(-100-100*mp.mpf('1e152')),
            pressure_order=9,width_order=2,transition_steps=4)
        bridge=AxialPressureWidthExitBridge(comparison,steps=8)
        continuation=PressureWidthExitContinuation(bridge)
        print('integrating differentiated collar once',flush=True)
        functions=continuation.functions('.3'); start=functions['start']
        encode=lambda v:{str(k):signed_log(a,70) for k,a in v.atoms.items()}
        rows=[]
        for radius in (1,10,100):
            end=continuation.evaluate_R(radius,'.3')
            increments={k:end[k]-start[k].value for k in ('F','Uz')}
            assert all(v.component(0,1)!=0 for v in increments.values())
            moment_increments={k:end['moments'][k]-start['moments'][k].value for k in end['moments']}
            assert end['Uz_Z'].component(1,0)!=0
            rows.append(dict(R=radius,field_increments={k:encode(v) for k,v in increments.items()},
                moments={k:encode(v) for k,v in moment_increments.items()},
                fields={k:encode(end[k]) for k in ('F','F_Z','Uz','Uz_Z','P','P_Z','Ur')},
                shears={k:encode(end['stress'][k]) for k in ('S_theta','S_z')},
                positive_epsilon_field_changes_retained=True))
            print(f'analytic continuation R={radius} complete',flush=True)
        # Check the finite exponential driver representation at an unused radius.
        y=mp.log(3); R=functions['Rb']*3
        direct=continuation.auxiliary.evaluate(R,'.3')
        driver_errors={}
        for name,a,b in [('D',functions['D'].evaluate(y),direct['D']),
                         ('J',functions['J'].evaluate(y),(R/2).sqrt()*direct['I_z'])]:
            keys=set(a.value.atoms)|set(b.value.atoms)
            driver_errors[name]=max((abs(a.value.component(*k)-b.value.component(*k))/
                max(1,abs(b.value.component(*k))) for k in keys),default=mp.mpf(0))
        assert max(driver_errors.values())<mp.mpf('1e-200'),driver_errors
        report=dict(Z='.3',metadata=continuation.metadata(),rows=rows,
            independent_driver_representation_errors={k:mp.nstr(v,40) for k,v in driver_errors.items()},
            exact_within_finite_ring=True,untruncated_ODE_error_enclosed=False,
            switching_reshape_installed=False,global_matching_complete=False,
            temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('actual positive-epsilon continuation receipt saved',flush=True)
    return report


if __name__=='__main__':run()
