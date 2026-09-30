"""Actual-source auxiliary exit comparison with retained pressure components."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_exit_comparison import Section923Comparison
from lei_ren_part1_paper_component_exit_comparison import ComponentExitComparison,fixture
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    print('building actual-source component exit comparison',flush=True)
    bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
        logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
        coherent_waiting=True,complete_preheat_components=True)
    with mp.workdps(bundle['precision']):
        hb=mp.exp(-100-100*mp.mpf('1e152'))
        component=ComponentExitComparison(bundle,h_b=hb,pressure_order=9,transition_steps=8)
        scalar=Section923Comparison(bundle,h_b=hb,transition_steps=8)
        points=[mp.mpf(0),mp.mpf('1.5')*hb,2*hb,3*hb,mp.mpf('.02')]
        samples=[];prefix_errors=[]
        for y in points:
            print('evaluating y/h_b or frozen branch',mp.nstr(y/hb,12),flush=True)
            a=component.evaluate(y,'.3');b=scalar.evaluate(y,'.3')
            prefix={}
            for name in ('F','Uz','P','PZ','D','E','T_theta','T_z'):
                base=a[name].component(0);reference=b[name]
                error=abs(base/reference-1) if reference else abs(base)
                prefix[name]=mp.nstr(error,40);prefix_errors.append(error)
            assert a['P'].component(1)!=0
            assert a['Uz'].component(1)!=0
            samples.append(dict(y=signed_log(y,60),y_over_hb=signed_log(y/hb,60),
                components={name:{str(k):signed_log(v,60) for k,v in a[name].atoms.items()}
                    for name in ('F','Uz','P','PZ','D','E','T_theta','T_z')},
                moments={name:{str(k):signed_log(v,60) for k,v in jet.atoms.items()}
                    for name,jet in a['moments'].items()},prefix_relative_replay=prefix))
        maximum=max(prefix_errors)
        assert maximum<mp.mpf('1e-150'),maximum
        report=dict(Z='.3',actual_h_b=signed_log(hb,80),samples=samples,
            zero_pressure_parameter_scalar_replay_max=mp.nstr(maximum,40),fixture=fixture(),
            metadata=component.metadata(),actual_pressure_tail_retained_through_auxiliary_exit=True,
            actual_exit_width_increments_may_be_below_atom_precision=True,
            pressure_parameter_remainder_enclosed=False,ODE_error_enclosed=False,
            prescribed_exit_bridge_installed=False,annular_matching_installed=False,
            five_terminal_moments_closed=False,pressure_compatibility_certified=False,
            finite_energy_certified=False,stress_cone_certified=False,temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('actual auxiliary exit replay complete; prefix error',report['zero_pressure_parameter_scalar_replay_max'],flush=True)
    return report


if __name__=='__main__':run()
