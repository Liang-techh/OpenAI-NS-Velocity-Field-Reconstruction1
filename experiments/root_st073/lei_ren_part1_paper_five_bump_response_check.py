"""Actual serialized centered defects: separated formal bump responses."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_axial_correction import from_signed_log,signed_log
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response


def run():
    source_path=Path(__file__).with_name('lei_ren_part1_paper_centered_component_defects_check.json')
    source=json.loads(source_path.read_text())
    with mp.workdps(260):
        def jet(rows):
            atoms={tuple(int(i) for i in key.strip('()').split(',')):from_signed_log(v) for key,v in rows.items()}
            return PressureWidthJet(atoms,pressure_order=9,width_order=2)
        d=[AxialDual(jet(source['defects'][str(i)]),jet(source['defects_Z'][str(i)])) for i in range(1,6)]
        z=mp.mpf('.3');am=mp.exp(mp.mpf('13.4'))/(1+z*z)
        amplitude=AxialDual(am,-2*z*am/(1+z*z),pressure_order=9,width_order=2)
        moment_map=FiveBumpMomentMap(precision=260,order=96)
        response=build_response(moment_map,amplitude,degree=3)
        terms=response.evaluate_terms(d)
        def encode(v):
            return {'value':{str(k):signed_log(x,120) for k,x in v.value.atoms.items()},
                    'Z':{str(k):signed_log(x,120) for k,x in v.tangent.atoms.items()}}
        encoded={str(k):[encode(v) for v in rows] for k,rows in terms.items()}
        unit3=(0,0,1,0,0);unit5=(0,0,0,0,1)
        assert unit3 in terms and unit5 in terms
        assert any(v.value.component(0,0)!=0 for v in terms[unit3])
        assert any(v.value.component(0,0)!=0 for v in terms[unit5])
        # Preserve original source labels independently in the linear response.
        source_linear={}
        zero=PressureWidthJet(0,pressure_order=9,width_order=2)
        for label,rows in source['parts'].items():
            direction=[AxialDual(jet(rows[str(i)]) if str(i) in rows else zero,
                       jet(source['parts_Z'][label][str(i)]) if str(i) in source['parts_Z'][label] else zero) for i in range(1,6)]
            source_linear[label]=[encode(-v) for v in moment_map.linear_inverse(direction)]
        report=dict(source=str(source_path.name),source_scope='Committed 120-digit aggregate centered defect jets and first-Z data; original source labels retained in linear response only',
          formal_defect_degree=3,pressure_order=9,width_order=2,
          separated_monomial_responses=encoded,linear_source_responses=source_linear,
          tiny_d3_and_d5_responses_nonzero=True,
          temporal_recursion=False,functional_closure=False,actual_corrected_field_installed=False,
          nonlinear_source_stage_products_expanded=False,formal_degree_remainder_enclosed=False,
          quadrature_enclosed=False,cone_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual serialized five-defect response saved; tiny d3/d5 retained',flush=True)
        return report

if __name__=='__main__':run()
