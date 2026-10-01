"""Reproducible independent physical-moment replay of persisted controls."""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_candidate_bump_integral_check import integrate_rows

HERE = Path(__file__).resolve().parent
INPUT = 'lei_ren_part1_paper_candidate_five_bump_inverse.json'


def decode(packet):
    def jet(rows):
        return PressureWidthJet({(r['pressure_order'],r['width_order']):
            mp.make_mpf(tuple(r['exact_mpf_tuple'])) for r in rows},pressure_order=9,width_order=2)
    return AxialDual(jet(packet['value']),jet(packet['first_Z']))


def maximum(vector,part):
    return max((abs(v) for dual in vector for v in getattr(dual,part).atoms.values()),default=mp.mpf(0))


def run():
    raw = (HERE/INPUT).read_bytes()
    data = json.loads(raw)
    for name,digest in data['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Candidate inverse input changed: '+name)
    with mp.workdps(180):
        h = tuple(decode(v) for v in data['materialized_controls'])
        amplitude = decode(data['amplitude'])
        defects = tuple(decode(v) for v in data['source_defects'])
        geometry = FiveBumpMomentMap(precision=180,order=data['bump_quadrature_order'])
        rows = integrate_rows(geometry,h,amplitude,mp.mpf(data['center']),order=128,precision=180)
        residual = tuple(rows[i]+defects[i] for i in range(5))
        def packed(dual):
            return {part:[dict(pressure_order=p,width_order=w,exact_mpf_tuple=list(v._mpf_))
                for (p,w),v in sorted(getattr(dual,part).atoms.items())] for part in ('value','tangent')}
        report = dict(source_inverse_sha256=hashlib.sha256(raw).hexdigest(),center=data['center'],
                      independent_quadrature_order=128,precision=180,
                      production_quadrature_order=data['bump_quadrature_order'],Rm_chart='1',
                      physical_scaling_Rm_three_fixture_checked_in_helper=True,
                      rows=[packed(v) for v in rows],rows_plus_source_defects=[packed(v) for v in residual],
                      maximum_value_residual=mp.nstr(maximum(residual,'value'),80),
                      maximum_first_Z_residual=mp.nstr(maximum(residual,'tangent'),80),
                      independent_physical_polarized_integration_used=True,
                      production_map_apply_or_cached_weights_used_for_integration=False,
                      beta_normalization_reconstructed_at_precision=180,
                      actual_source_profile_replayed=False,
                      reference_profile='paper analytic power law on bump supports',
                      source_defects_from_actual_candidate=True,
                      quadrature_error_enclosed=False,functional_closure=False,
                      input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in (INPUT,'lei_ren_part1_paper_candidate_bump_integral_check.py',
                                       'lei_ren_part1_paper_five_bump_map.py')},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Independent candidate physical-moment replay:',report['maximum_value_residual'],
              'first Z',report['maximum_first_Z_residual'],flush=True)
        return report


if __name__ == '__main__':
    run()
