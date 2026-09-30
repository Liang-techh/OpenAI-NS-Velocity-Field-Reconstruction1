"""Actual serialized P9/W2 correction field, retaining each defect monomial."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_axial_correction import from_signed_log,signed_log
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response
from lei_ren_part1_paper_five_bump_response_field import evaluate_response_correction


def run():
    src=Path(__file__).with_name('lei_ren_part1_paper_centered_component_defects_check.json')
    data=json.loads(src.read_text())
    with mp.workdps(260):
        def jet(rows):
            return PressureWidthJet({tuple(int(i) for i in k.strip('()').split(',')):from_signed_log(v) for k,v in rows.items()},pressure_order=9,width_order=2)
        d=[AxialDual(jet(data['defects'][str(i)]),jet(data['defects_Z'][str(i)])) for i in range(1,6)]
        z=mp.mpf('.3');am=mp.exp(mp.mpf('13.4'))/(1+z*z)
        amplitude=AxialDual(am,-2*z*am/(1+z*z),pressure_order=9,width_order=2)
        moment_map=FiveBumpMomentMap(precision=260,order=96)
        response=build_response(moment_map,amplitude,degree=3)
        rm=110*mp.exp(10*(mp.mpf('5e151')+14)-6)
        print('actual termwise velocity and full quadratic cumulative moments',flush=True)
        result=evaluate_response_correction(response,d,Rm=rm,Z=z,x='1.25',delta='1e-200')
        def encode(v):
            if isinstance(v,AxialDual):
                return dict(baseline=signed_log(v.value.component(0,0),120),
                    baseline_Z=signed_log(v.tangent.component(0,0),120),
                    value_atom_count=len(v.value.atoms),Z_atom_count=len(v.tangent.atoms))
            return dict(baseline=signed_log(v.component(0,0) if hasattr(v,'component') else v,120),
                        value_atom_count=len(v.atoms) if hasattr(v,'atoms') else 1)
        velocity={str(k):[encode(v) for v in rows] for k,rows in result['velocity_terms'].items()}
        moments={str(k):[encode(v) for v in rows] for k,rows in result['moment_terms'].items()}
        radial={str(k):encode(v) for k,v in result['radial_velocity_terms'].items()}
        for key in [(0,0,1,0,0),(0,0,0,0,1)]:
            assert any(v.value.component(0,0)!=0 for v in result['velocity_terms'][key])
        high=sum(any(v.value.component(0,0)!=0 for v in rows) for key,rows in result['moment_terms'].items() if sum(key)>3)
        assert high>0
        report=dict(x='1.25',source=src.name,source_scope='Committed 120-digit actual aggregate defect jets and first-Z data',
          pressure_order=9,width_order=2,metadata=result['metadata'],
          receipt_scope='Per-monomial baseline and first-Z logs, plus finite pressure/width atom counts; calculations retain the full declared ring',
          velocity_terms=velocity,moment_terms=moments,radial_velocity_terms=radial,
          nonzero_higher_quadratic_moment_terms=high,tiny_d3_d5_velocity_responses_preserved=True,
          full_baseline_join_installed=False,functional_closure=False,cone_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('actual field correction receipt saved',flush=True)
        return report

if __name__=='__main__':run()
