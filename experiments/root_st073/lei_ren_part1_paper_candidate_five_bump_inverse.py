"""Same-candidate finite five-bump inverse and repaired reference field.

The finite inverse preserves incremental responses and first axial tangents.
No point calculation is promoted to an axial functional closure certificate.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_component_defects import CandidateEndpointSource
from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse_field import FiveBumpInverseField

HERE = Path(__file__).resolve().parent


def atoms(jet):
    return [dict(pressure_order=p,width_order=w,exact_mpf_tuple=list(v._mpf_))
            for (p,w),v in sorted(jet.atoms.items())]


def dual_packet(dual):
    return dict(value=atoms(dual.value),first_Z=atoms(dual.tangent))


def max_atom(vector):
    return max((abs(v) for item in vector for jet in (item.value,item.tangent)
                for v in jet.atoms.values()),default=mp.mpf(0))


def build_problem():
    source = CandidateEndpointSource()
    centered = CenteredComponentDefects(source,order=16,restore_order=32)
    reference = ReferenceDefectBackground(centered,delta='1e-200')
    moment_map = FiveBumpMomentMap(precision=centered.work_precision,order=64)
    adapter = FiveBumpInverseField(reference,moment_map,steps=10)
    return source,centered,reference,moment_map,adapter


def run():
    source,centered,reference,moment_map,adapter = build_problem()
    with mp.workdps(moment_map.precision):
        z = source.center
        amplitude,inverse = adapter.inverse_data(z)
        baseline = reference.evaluate_x('1.25',z)
        field = adapter.evaluate_x('1.25',z)
        if field['P0'] != baseline['P0'] or field['P0_Z'] != baseline['P0_Z']:
            raise AssertionError('Repair changed the prescribed axis pressure datum')
        if field['stress'] is None:
            raise AssertionError('Corrected field stress failed: '+str(field.get('stress_error')))
        history = []
        for row in inverse['iterations']:
            error = max_atom(row['residual'])
            print('Candidate five-bump update',row['update'],'maximum retained value/tangent residual',
                  mp.nstr(error,12),flush=True)
            history.append(dict(update=row['update'],residual_max_atom=mp.nstr(error,70),
                                residual=[dual_packet(v) for v in row['residual']]))
        data = reference.defect_data(z)
        # Keep cancellation-resistant hierarchy and materialized algebra
        # separate: short finite precision sums cannot retain every tiny term.
        direct = moment_map.apply(inverse['h'],amplitude)
        replay = [direct[i]+data['d_dual'][i+1] for i in range(5)]
        if max_atom(inverse['terminal_residual']) >= mp.mpf('1e-180'):
            raise AssertionError('Finite inverse did not sufficiently reduce the candidate defect')
        inputs = ('lei_ren_part1_paper_candidate_transition_analytic.json',
                  'lei_ren_part1_paper_candidate_component_defects.py',
                  'lei_ren_part1_paper_centered_component_defects.py',
                  'lei_ren_part1_paper_five_moment_reference_background.py',
                  'lei_ren_part1_paper_five_bump_inverse.py',
                  'lei_ren_part1_paper_five_bump_inverse_field.py',
                  'lei_ren_part1_paper_five_bump_map.py','lei_ren_part1_paper_five_bump_field.py')
        report = dict(center='.3',Lambda='1e120',actual_candidate_endpoint_used=True,
                      core_state_sha256=source.receipt['state_sha256'],
                      precision=moment_map.precision,bump_quadrature_order=64,nonlinear_updates=10,
                      coefficient_order=list(moment_map.coefficient_order),
                      coefficient_increments=[[dual_packet(v) for v in step] for step in inverse['increments']],
                      materialized_controls=[dual_packet(v) for v in inverse['h']],
                      amplitude=dual_packet(amplitude),
                      source_defects=[dual_packet(data['d_dual'][i]) for i in range(1,6)],
                      iterations=history,
                      polarization_terminal_residual=[dual_packet(v) for v in inverse['terminal_residual']],
                      materialized_same_map_residual=[dual_packet(v) for v in replay],
                      field_sample_x='1.25',
                      corrected_field={name:atoms(field[name]) for name in ('F','Uz','P','P0','Ur')},
                      P0_preserved=True,analytic_first_axial_tangents_used=True,
                      temporal_recursion=False,exact_infinite_inverse=False,
                      independent_integral_validation_available=False,
                      quadrature_error_enclosed=False,source_error_enclosed=False,
                      whole_axis_inverse=False,uniform_cone_certified=False,
                      terminal_functional_five_moment_closure=False,
                      input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Actual candidate finite five-bump inverse and corrected field generated; independent integral and functional closure still open',flush=True)
        return report


if __name__ == '__main__':
    run()
