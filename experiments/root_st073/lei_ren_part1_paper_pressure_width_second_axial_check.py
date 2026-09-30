"""Fresh actual third-depth core source for the second-Z comparison.

An independent local cache is written; the old first-Z R100 cache is untouched.
This receipt covers the inner comparison endpoint, not the R110 defect chain.
"""
import json
import pickle
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_pressure_width_second_axial_comparison import SecondAxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_second_axial_bridge import SecondAxialPressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_second_axial_continuation import SecondAxialPressureWidthExitContinuation
from lei_ren_part1_paper_pressure_width_second_axial_switches import SecondAxialPressureWidthExitSwitches
from lei_ren_part1_paper_axial_correction import signed_log


def run(resume=False,post_collar=False,switches=False):
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'second_axial_core_Z03.pkl'
    if resume:
        with cache.open('rb') as stream:data=pickle.load(stream)
        if data['version']!=1 or data['Z']!='.3' or data['Z_depth']!=3:
            raise ValueError('Unsupported local second-Z core snapshot')
        axis=SimpleNamespace(precision=260,Lambda=data['Lambda'],delta=data['delta'])
        def coefficients(z):
            if mp.mpf(str(z))!=mp.mpf('.3'):raise ValueError('Snapshot supports Z=.3 only')
            return data['coefficients']
        component=SimpleNamespace(axis=axis,precision=260,coefficients=coefficients)
        bundle=dict(precision=260,axis=axis,component_pressure_core=component)
    else:
        print('building actual source with third-depth axial core data',flush=True)
        bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
            logC='5e151',logPstar='14',delta='1e-200',continuous_pressure=True,
            coherent_waiting=True,complete_preheat_components=True,component_Z_jet_depth=3)
        coefficients=bundle['component_pressure_core'].coefficients('.3')
        axis=bundle['axis'];cache.parent.mkdir(parents=True,exist_ok=True)
        with cache.open('wb') as stream:
            pickle.dump(dict(version=1,Z='.3',Z_depth=3,Lambda=axis.Lambda,
                             delta=axis.delta,coefficients=coefficients),stream)
        print('fresh actual third-depth source cached independently',flush=True)
    with mp.workdps(260):
        comparison=SecondAxialPressureWidthComparison(bundle,
            h_b=mp.exp(-100-100*mp.mpf('1e152')),pressure_order=9,width_order=2)
        value=comparison.evaluate(0,'.3')
        def encode_jet(jet):return {str(k):signed_log(v,70) for k,v in jet.atoms.items()}
        def encode(v):return dict(value=encode_jet(v.value),first_Z=encode_jet(v.tangent),second_Z=encode_jet(v.second))
        report=dict(Z='.3',s=0,source_Z_depth=3,source_radial_degree=18,precision=260,
            pressure_order=9,width_order=2,resumed_new_second_Z_cache=resume,
            fields={k:encode(value[k]) for k in ('F','Uz','P','I_theta','I_z','D','E')},
            moments={k:encode(v) for k,v in value['moments'].items()},
            pressure_datum=encode(comparison.core_coefficients('.3')['P'][0][0]),
            source_scope='actual same-preheat core at the inner comparison endpoint only',
            downstream_R110_second_Z_installed=False,defect_second_Z_installed=False,
            source_error_enclosed=False,uniform_Z_certified=False,
            endpoint_Z_extension_certified=False,temporal_recursion=False)
        if post_collar or switches:
            print('propagating actual second-Z collar and analytic continuation to R100',flush=True)
            bridge=SecondAxialPressureWidthExitBridge(comparison,steps=8)
            continuation=SecondAxialPressureWidthExitContinuation(bridge)
            extended=cache.with_name('second_axial_R100_source_Z03.pkl')
            if resume and switches and extended.exists():
                with extended.open('rb') as stream:existing=pickle.load(stream)
                if existing['version']!=1 or existing['Z']!='.3' or existing['Z_depth']!=3:
                    raise ValueError('Unsupported local second-Z R100 snapshot')
                def cached_functions(z):
                    if mp.mpf(str(z))!=mp.mpf('.3'):raise ValueError('R100 snapshot supports Z=.3 only')
                    return existing['functions']
                continuation.functions=cached_functions
            outer=continuation.evaluate_R(100,'.3')
            report['post_collar_R100']=dict(
                fields={k:encode(v) for k,v in outer['second_jet_fields'].items()},
                moments={k:encode(v) for k,v in outer['second_jet_moments'].items()},
                raw_quadratic={k:encode(v) for k,v in outer['second_jet_raw_quadratic_integrals'].items()},
                Ur_Z=encode_jet(outer['Ur_Z']),P0_ZZ=encode_jet(outer['P0_ZZ']),
                second_Z_available=True,third_Z_not_available=True,
                finite_ring_remainder_enclosed=False,collar_RK_error_enclosed=False)
            with extended.open('wb') as stream:
                pickle.dump(dict(version=1,Z='.3',Z_depth=3,Lambda=comparison.Lambda,
                    delta=comparison.delta,coefficients=bundle['component_pressure_core'].coefficients('.3'),
                    functions=continuation.functions('.3')),stream)
            print('actual second-Z R100 source saved in independent cache',flush=True)
            report['source_scope']='Actual same-preheat inner core, collar and analytic R100 continuation at Z=.3'
            if switches:
                print('propagating actual second-Z R100--R110 switches',flush=True)
                switch=SecondAxialPressureWidthExitSwitches(continuation,steps=8)
                boundary=switch.evaluate_R(100,'.3')
                errors=[]
                def compare(a,b):
                    errors.extend(abs(a.component(*key)-b.component(*key))/max(1,abs(b.component(*key)))
                                  for key in set(a.atoms)|set(b.atoms))
                for name in ('F','Uz','P'):
                    for suffix in ('','_Z','_ZZ'):compare(boundary[name+suffix],outer[name+suffix])
                for name in outer['moments']:
                    for key in ('moments','moments_Z','moments_ZZ'):compare(boundary[key][name],outer[key][name])
                assert max(errors)<mp.mpf('1e-200'), 'Second-Z switch inlet mismatch'
                end=switch.evaluate_R(110,'.3')
                report['R110']=dict(fields={k:encode(v) for k,v in end['second_jet_fields'].items()},
                    moments={k:encode(v) for k,v in end['second_jet_moments'].items()},
                    raw_quadratic={k:encode(v) for k,v in end['second_jet_raw_quadratic_integrals'].items()},
                    P0_ZZ=encode_jet(end['P0_ZZ']),Ur_Z=encode_jet(end['Ur_Z']),
                    inlet_scaled_error=mp.nstr(max(errors),70),second_Z_available=True,
                    switch_RK_error_enclosed=False)
                report['downstream_R110_second_Z_installed']=True
                report['source_scope']='Actual same-preheat inner core, collar, analytic R100 continuation and R110 switches at Z=.3'
                with cache.with_name('second_axial_R110_source_Z03.pkl').open('wb') as stream:
                    pickle.dump(dict(version=1,Z='.3',Z_depth=3,Lambda=comparison.Lambda,
                        delta=comparison.delta,source=end),stream)
                print('actual second-Z R110 source saved',flush=True)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual second-Z inner comparison receipt saved',flush=True)
        return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true')
    parser.add_argument('--post-collar',action='store_true')
    parser.add_argument('--switches',action='store_true');args=parser.parse_args()
    run(resume=args.resume,post_collar=args.post_collar,switches=args.switches)
