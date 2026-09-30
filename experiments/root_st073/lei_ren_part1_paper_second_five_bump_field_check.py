"""Actual C2 source/inverse joined to the same reference similarity field."""
import json
import pickle
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_second_five_bump_field import evaluate_second_bump_field
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'second_centered_defects_Z03.pkl'
    with cache.open('rb') as stream:snapshot=pickle.load(stream)
    if snapshot['version']!=1 or snapshot['Z']!='.3':raise ValueError('Unsupported actual C2 snapshot')
    with mp.workdps(snapshot['precision']):
        moment_map=FiveBumpMomentMap(precision=snapshot['precision'],order=96)
        data=snapshot['data'];inverse=snapshot['inverse']
        result=evaluate_second_bump_field(moment_map,inverse,data,Rm=snapshot['Rm'],
                                         Z='.3',x='1.25',delta=snapshot['delta'])
        datum=result['second_jet_fields']['P0']
        assert datum.second==data['P0_ZZ'] and datum.value==data['P0']
        a=result['Ur'];b=result['stress']['U_r']
        ur_error=max((abs(a.component(*k)-b.component(*k))/max(1,abs(b.component(*k)))
                     for k in set(a.atoms)|set(b.atoms)),default=mp.mpf(0))
        assert ur_error<mp.mpf('1e-400'), 'Shared radial velocity recovery mismatch'
        def encode_jet(v):return {str(k):signed_log(a,70) for k,a in v.atoms.items()}
        def encode(v):return dict(value=encode_jet(v.value),first_Z=encode_jet(v.tangent),second_Z=encode_jet(v.second))
        report=dict(Z='.3',x='1.25',precision=snapshot['precision'],pressure_order=9,width_order=2,
            fields={k:encode(v) for k,v in result['second_jet_fields'].items()},
            moments={k:encode(v) for k,v in result['second_jet_moments'].items()},
            moment_parts={label:{k:encode(v) for k,v in rows.items()}
                          for label,rows in result['second_jet_moment_parts'].items()},
            Ur=encode_jet(result['Ur']),Ur_Z=encode_jet(result['Ur_Z']),Ur_ZZ_available=False,
            shared_stress_Ur_scaled_error=mp.nstr(ur_error,70),
            source_part_count=result['source_part_count'],physical_part_count=len(result['second_jet_moment_parts']),
            P0_second_preserved=True,finite_inverse_updates=10,
            source_scope='actual same-preheat C2 R110 and centered source at Z=.3 only',
            uniform_Z_certified=False,functional_closure=False,source_error_enclosed=False,
            infinite_inverse_error_enclosed=False,cone_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual C2 reference/inverse field joined; source parts',result['source_part_count'],flush=True)


if __name__=='__main__':run()
