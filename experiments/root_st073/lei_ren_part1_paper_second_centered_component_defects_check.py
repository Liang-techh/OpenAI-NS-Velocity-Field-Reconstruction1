"""Actual common R110 source: second-Z centered rows and inverse increments."""
import json
import pickle
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_second_centered_component_defects import SecondCenteredComponentDefects
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'second_axial_R110_source_Z03.pkl'
    with cache.open('rb') as stream:snapshot=pickle.load(stream)
    if snapshot['version']!=1 or snapshot['Z']!='.3' or snapshot['Z_depth']!=3:
        raise ValueError('Unsupported internal actual second-Z source snapshot')
    source=snapshot['source']
    def evaluate_R(R,Z):
        if mp.mpf(R)!=110 or mp.mpf(Z)!=mp.mpf('.3'):
            raise ValueError('Source snapshot supplies R110/Z=.3 only')
        return source
    switches=SimpleNamespace(precision=260,pressure_order=9,width_order=2,
        comparison=SimpleNamespace(delta=snapshot['delta']),evaluate_R=evaluate_R,
        shared_parameters=dict(logC='5e151',logPstar='14'))
    defects=SecondCenteredComponentDefects(switches,A='1e150',logC='5e151',logPstar='14',
                                          order=32,window=24,restore_order=64)
    print('assembling actual C2 centered defects including N+2 flat sensitivities',flush=True)
    with mp.workdps(defects.work_precision):
        data=defects.evaluate('.3')
        assert data['P0']==source['P0'] and data['P0_ZZ']==source['P0_ZZ']
        assert len(data['parts_ZZ'])==69
        for row in range(1,6):
            total=sum((p.get(row,0) for p in data['parts_ZZ'].values()),defects._jet(0))
            assert not (total-data['d_ZZ'][row]).atoms
        moment_map=FiveBumpMomentMap(precision=defects.work_precision,order=96)
        amplitude=AxialSecondJet(data['Am'],data['Am_Z'],data['Am_ZZ'])
        inverse=iterate_inverse(moment_map,[data['d_second_jet'][i] for i in range(1,6)],amplitude,steps=10)
        def encode_jet(v):return {str(k):signed_log(a,70) for k,a in v.atoms.items()}
        def encode(v):return dict(value=encode_jet(v.value),first_Z=encode_jet(v.tangent),second_Z=encode_jet(v.second))
        report=dict(Z='.3',precision=defects.work_precision,pressure_order=9,width_order=2,
            source_scope='fresh third-depth core and same actual C2 collar/R100/R110 source',
            d={str(i):encode(v) for i,v in data['d_second_jet'].items()},
            parts={label:{str(i):encode(v) for i,v in rows.items()} for label,rows in data['parts_second_jet'].items()},
            flat_kernel_cache_orders={k:max(rec['derivative_cache']) for k,rec in data['flat_records'].items()},
            second_Z_part_count=len(data['parts_ZZ']),P0=encode(data['P0_dual']),P0_ZZ_preserved=True,
            nonlinear_inverse_updates=10,inverse_increments=[[encode(v) for v in values] for values in inverse['increments']],
            inverse_terminal_residual=[encode(v) for v in inverse['terminal_residual']],
            inverse_max_nominal_second_Z_residual=signed_log(max(abs(v.second.evaluate()) for v in inverse['terminal_residual']),70),
            source_error_enclosed=False,flat_quadrature_enclosed=False,
            uniform_Z_certified=False,infinite_inverse_error_enclosed=False,
            functional_closure=False,temporal_recursion=False)
        output_cache=cache.with_name('second_centered_defects_Z03.pkl')
        with output_cache.open('wb') as stream:
            pickle.dump(dict(version=1,Z='.3',precision=defects.work_precision,Rm=defects.Rm,
                             delta=snapshot['delta'],data=data,inverse=inverse),stream)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual second-Z defect source and inverse receipt saved',flush=True)
        return report


if __name__=='__main__':run()
