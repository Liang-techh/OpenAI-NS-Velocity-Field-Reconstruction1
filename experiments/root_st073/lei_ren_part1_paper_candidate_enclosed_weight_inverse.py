"""Re-solve the actual candidate using midpoints of directed bump weights.

The finite center solve is followed by a directed integral residual check.
Weight-midpoint arithmetic is not itself an interval inverse or a source
error certificate. The original finite-Gauss receipt remains unchanged.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_bump_integral_replay import decode
from lei_ren_part1_paper_candidate_five_bump_inverse import dual_packet, max_atom
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_candidate_bump_enclosed_residual import run as enclose_residual

HERE = Path(__file__).resolve().parent
ORIGINAL = 'lei_ren_part1_paper_candidate_five_bump_inverse.json'
WEIGHTS = 'lei_ren_part1_paper_bump_integral_enclosures_check.json'


class DirectedWeightMidpointMap(FiveBumpMomentMap):
    """Full-support algebra only; no inherited partial-field quadrature."""

    def __init__(self, receipt, precision=473):
        self.precision = precision
        with mp.workdps(precision):
            self.centers = tuple(mp.mpf(v) for v in ('1.25', '1.5', '1.75'))
            table = {}
            for row in receipt['weight_records'].values():
                bound = row['weight_interval']
                lo = mp.make_mpf(tuple(bound['lower_exact_mpf_tuple']))
                hi = mp.make_mpf(tuple(bound['upper_exact_mpf_tuple']))
                table[Fraction(row['center']), Fraction(row['power']), row['multiplicity']] = (lo+hi)/2
            def weight(i, power, k=1):
                return table[Fraction(str(self.centers[i])), Fraction(power), k]
            self.linear_matrix = mp.matrix(5, 5)
            # Exact normalized mass, rather than a midpoint of a dependent ratio.
            self.linear_matrix[0, 0] = self.linear_matrix[0, 1] = 1
            for j, i in enumerate((0, 2)):
                self.linear_matrix[1, j] = weight(i, '.6')
            for i in range(3):
                self.linear_matrix[2, i+2] = weight(i, '.5')
                self.linear_matrix[3, i+2] = -weight(i, '.1')
                self.linear_matrix[4, i+2] = weight(i, '-.9')
            self.inverse_matrix = mp.inverse(self.linear_matrix)
            fg, gg, ff, ff_over_x = mp.matrix(2, 3), mp.matrix(2), mp.matrix(3), mp.matrix(3)
            for j, i in enumerate((0, 2)):
                fg[j, i] = weight(i, '.5', 2)
                gg[j, j] = weight(i, '0', 2)
            for i in range(3):
                ff[i, i] = weight(i, '0', 2)
                ff_over_x[i, i] = weight(i, '-1', 2)
            self.quadratic_weights = dict(fg_sqrt=fg, gg=gg, ff=ff, ff_over_x=ff_over_x)

    def _weights_for_interval(self, left, right):
        raise NotImplementedError('This map supplies full support algebra only')


def run():
    original = json.loads((HERE/ORIGINAL).read_bytes())
    weights = json.loads((HERE/WEIGHTS).read_bytes())
    for name, digest in original['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Original source input changed: '+name)
    for name, digest in (
        ('lei_ren_part1_paper_bump_integral_enclosures.py', weights['module_sha256']),
        ('lei_ren_part1_paper_bump_integral_enclosures_check.py', weights['check_source_sha256']),
    ):
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Directed weight input changed: '+name)
    moment_map = DirectedWeightMidpointMap(weights)
    with mp.workdps(moment_map.precision):
        amplitude = decode(original['amplitude'])
        defects = tuple(decode(row) for row in original['source_defects'])
        inverse = iterate_inverse(moment_map, defects, amplitude, steps=10)
        old = tuple(decode(row) for row in original['materialized_controls'])
        if max_atom(inverse['terminal_residual']) >= mp.mpf('1e-180'):
            raise ArithmeticError('Refined finite inverse failed its arithmetic check')
        inputs = (ORIGINAL, WEIGHTS, 'lei_ren_part1_paper_five_bump_map.py',
                  'lei_ren_part1_paper_five_bump_inverse.py',
                  'lei_ren_part1_paper_candidate_bump_integral_replay.py',
                  'lei_ren_part1_paper_candidate_five_bump_inverse.py',
                  'lei_ren_part1_paper_candidate_bump_enclosed_residual.py')
        report = dict(center=original['center'], Lambda=original['Lambda'],
                      amplitude=original['amplitude'], source_defects=original['source_defects'],
                      materialized_controls=[dual_packet(row) for row in inverse['h']],
                      changes_from_original=[dual_packet(inverse['h'][i]-old[i]) for i in range(5)],
                      coefficient_increments=[[dual_packet(row) for row in increment]
                                              for increment in inverse['increments']],
                      maximum_arithmetic_polarization_residual=mp.nstr(max_atom(inverse['terminal_residual']),70),
                      directed_weight_midpoints_used=True, precision=moment_map.precision,
                      nonlinear_updates=10, original_source_targets_preserved=True,
                      original_inverse_receipt_unchanged=True, full_support_map_only=True,
                      corrected_field_installed=False, rigorous_inverse_enclosed=False,
                      source_error_enclosed=False, functional_terminal_closure=False,
                      temporal_recursion=False, cone_certified=False,
                      input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                    for name in inputs},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        destination = Path(__file__).with_suffix('.json')
        destination.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print('Directed-weight midpoint candidate inverse generated:',
              report['maximum_arithmetic_polarization_residual'], flush=True)
    enclosed = enclose_residual(destination.name,
                               Path(__file__).with_name('lei_ren_part1_paper_candidate_bump_enclosed_refined_residual.json'))
    if not (enclosed['all_value_residuals_include_zero'] and enclosed['all_first_Z_residuals_include_zero']):
        raise ArithmeticError('Refined controls still exclude zero in directed integral residuals')
    return report, enclosed


if __name__ == '__main__':
    run()
