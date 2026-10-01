"""Install refined candidate controls with directed partial bump integrals.

Field values are midpoint approximations; partial integral bounds are retained
separately. Neither source errors nor the finite formal-ring remainder are
enclosed. The source intentionally accepts only its persisted axial center.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_bump_partial_enclosures import PaperBumpPartialEnclosures
from lei_ren_part1_paper_candidate_enclosed_weight_inverse import DirectedWeightMidpointMap
from lei_ren_part1_paper_candidate_bump_integral_replay import decode
from lei_ren_part1_paper_candidate_component_defects import CandidateEndpointSource
from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_five_bump_field import evaluate_correction, join
from lei_ren_part1_paper_candidate_repaired_cone_scan import cone
from lei_ren_part1_paper_candidate_five_bump_inverse import atoms

HERE = Path(__file__).resolve().parent
INVERSE = 'lei_ren_part1_paper_candidate_enclosed_weight_inverse.json'
WEIGHTS = 'lei_ren_part1_paper_bump_integral_enclosures_check.json'


def midpoint(interval):
    lo, hi = (mp.make_mpf(v) for v in interval._mpi_)
    return (lo + hi) / 2


def interval_packet(interval):
    lo, hi = (mp.make_mpf(v) for v in interval._mpi_)
    return dict(lower_exact_mpf_tuple=list(lo._mpf_), upper_exact_mpf_tuple=list(hi._mpf_),
                lower=mp.nstr(lo, 70), upper=mp.nstr(hi, 70), width=mp.nstr(hi-lo, 40))


class RefinedPartialBumpMap(DirectedWeightMidpointMap):
    """Full and partial weights from one directed paper-family backend."""

    def __init__(self, receipt, precision=473):
        super().__init__(receipt, precision)
        self.partial_backend = PaperBumpPartialEnclosures(precision=80, order=12, tolerance='1e-20')
        with mp.workdps(precision):
            self.radius = mp.mpf(1) / 40
            self.beta_normalization = midpoint(self.partial_backend.normalization)
        self._partial_cache = {}
        self.partial_bounds = {}

    def beta(self, x, center):
        with mp.workdps(self.precision):
            t = (mp.mpf(x) - center) / self.radius
            return mp.mpf(0) if abs(t) >= 1 else mp.exp(-1 / (1-t*t)) / (self.radius*self.beta_normalization)

    def _weights_for_interval(self, left, right):
        if mp.mpf(left) != 1:
            raise ValueError('Candidate cumulative map requires left=1')
        x = min(mp.mpf(right), mp.mpf(2))
        key = mp.nstr(x, self.precision)
        if key in self._partial_cache:
            return self._partial_cache[key]
        if not 1 <= x <= 2:
            raise ValueError('Require cumulative endpoint in [1,2]')
        bounds = {}
        def weight(i, power, k=1):
            interval = self.partial_backend.weight(str(self.centers[i]), power,
                                                    Fraction(key), multiplicity=k)
            bounds[(i, power, k)] = interval
            return midpoint(interval)
        matrix = mp.matrix(5, 5)
        for j, i in enumerate((0, 2)):
            matrix[0, j] = weight(i, '0')
            if x >= self.centers[i] + self.radius:
                matrix[0, j] = 1  # Exact normalized full support mass.
            matrix[1, j] = weight(i, '.6')
        for i in range(3):
            matrix[2, i+2] = weight(i, '.5')
            matrix[3, i+2] = -weight(i, '.1')
            matrix[4, i+2] = weight(i, '-.9')
        fg, gg, ff, ff_over_x = mp.matrix(2, 3), mp.matrix(2), mp.matrix(3), mp.matrix(3)
        for j, i in enumerate((0, 2)):
            fg[j, i] = weight(i, '.5', 2)
            gg[j, j] = weight(i, '0', 2)
        for i in range(3):
            ff[i, i] = weight(i, '0', 2)
            ff_over_x[i, i] = weight(i, '-1', 2)
        weights = dict(fg_sqrt=fg, gg=gg, ff=ff, ff_over_x=ff_over_x)
        self.partial_bounds[key] = bounds
        self._partial_cache[key] = matrix, weights
        return matrix, weights


class RefinedCandidateBumpField:
    """Callable actual-candidate matching field at its supplied axial center."""

    def __init__(self):
        self.inverse = json.loads((HERE/INVERSE).read_bytes())
        self.weights = json.loads((HERE/WEIGHTS).read_bytes())
        for name, digest in self.inverse['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Refined inverse input changed: '+name)
        self.source = CandidateEndpointSource()
        self.centered = CenteredComponentDefects(self.source, order=16, restore_order=32)
        self.reference = ReferenceDefectBackground(self.centered, delta='1e-200')
        self.precision = self.centered.work_precision
        self.moment_map = RefinedPartialBumpMap(self.weights, precision=self.precision)
        with mp.workdps(self.precision):
            self.controls = tuple(decode(v) for v in self.inverse['materialized_controls'])
            self.amplitude = decode(self.inverse['amplitude'])
        self._cache = {}

    def evaluate_x(self, x, Z='.3'):
        with mp.workdps(self.precision):
            x, z = mp.mpf(x), mp.mpf(Z)
            if abs(z-self.source.center) > mp.mpf('1e-255'):
                raise ValueError('Candidate source only supplies its axial center')
            key = mp.nstr(x, self.precision)
            if key not in self._cache:
                baseline = self.reference.evaluate_x(x, z)
                correction = evaluate_correction(self.moment_map, self.controls, self.amplitude,
                                                 self.reference.Rm, z, x, delta=self.reference.delta)
                field = join(baseline, correction, delta=self.reference.delta)
                if field['P0'] != baseline['P0'] or field['P0_Z'] != baseline['P0_Z']:
                    raise AssertionError('Refined field changed prescribed pressure datum')
                field.update(refined_controls_installed=True, partial_integral_bounds_available=True,
                             values_are_interval_midpoint_approximations=True,
                             normalization_interval=self.moment_map.partial_backend.normalization,
                             nominal_stress_only=True, source_error_enclosed=False,
                             uniform_cone_certified=False, terminal_functional_closure=False,
                             temporal_recursion=False)
                self._cache[key] = field
            return self._cache[key]


def run():
    candidate = RefinedCandidateBumpField()
    rows = []
    with mp.workdps(candidate.precision):
        for x in ('1', '1.225', '1.25', '1.275', '1.475', '1.5', '1.525',
                  '1.725', '1.75', '1.775', '2'):
            field = candidate.evaluate_x(x)
            if field['stress'] is None:
                raise ArithmeticError(field.get('stress_error'))
            row = dict(x=x, pointwise_cone=cone(field),
                       F=mp.nstr(field['F'].evaluate(), 50),
                       Uz=mp.nstr(field['Uz'].evaluate(), 50),
                       logR=mp.nstr(field['logR'], 70),
                       P0_preserved=field['P0_preserved'],
                       field_atoms={name: atoms(field[name]) for name in
                                    ('F', 'F_Z', 'Utheta', 'Utheta_Z', 'Uz', 'Uz_Z', 'Ur', 'P', 'P_Z', 'P0', 'P0_Z')},
                       moment_atoms={name: atoms(value) for name, value in field['moments'].items()},
                       moment_Z_atoms={name: atoms(value) for name, value in field['moments_Z'].items()},
                       partial_weight_bounds={str(k): interval_packet(v) for k, v in
                                             candidate.moment_map.partial_bounds[mp.nstr(mp.mpf(x),candidate.precision)].items()})
            rows.append(row)
            print('Refined partial-integral candidate field x', x,
                  'relaxed cone', row['pointwise_cone']['relaxed_passed'], flush=True)
        endpoint = candidate.evaluate_x('2')
        amplitude = candidate.amplitude
        partial = candidate.moment_map.partial(candidate.controls, amplitude, 1, 2)
        full = candidate.moment_map.apply(candidate.controls, amplitude)
        if any((partial[i]-full[i]).value.atoms or (partial[i]-full[i]).tangent.atoms for i in range(5)):
            raise AssertionError('Partial endpoint and full refined map disagree')
        inputs = (INVERSE, WEIGHTS, 'lei_ren_part1_paper_bump_partial_enclosures.py',
                  'lei_ren_part1_paper_bump_integral_enclosures.py',
                  'lei_ren_part1_paper_candidate_enclosed_weight_inverse.py',
                  'lei_ren_part1_paper_five_bump_field.py',
                  'lei_ren_part1_paper_five_moment_reference_background.py',
                  'lei_ren_part1_paper_candidate_component_defects.py')
        report = dict(center='.3', Lambda='1e120', rows=rows,
                      refined_controls_installed=True, P0_preserved=True,
                      partial_endpoint_matches_full_map=True,
                      all_sampled_relaxed_cone_passed=all(r['pointwise_cone']['relaxed_passed'] for r in rows),
                      actual_partial_weight_bounds_saved=True,
                      normalization=interval_packet(candidate.moment_map.partial_backend.normalization),
                      whole_axis=False, source_error_enclosed=False,
                      formal_ring_remainder_enclosed=False, uniform_cone_certified=False,
                      temporal_recursion=False, independent_Cartesian_residual_validated=False,
                      field_values_are_midpoint_approximations=True,
                      input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        return report


if __name__ == '__main__':
    run()
