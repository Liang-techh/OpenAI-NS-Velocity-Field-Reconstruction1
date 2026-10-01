"""Directed bump-integral residual of fixed candidate controls and targets.

Only the integral weights are enclosed. Source targets and coefficient data
are fixed nominal inputs, not enclosures of the reconstructed NS source.
"""
import hashlib
import json
from pathlib import Path
from fractions import Fraction

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

HERE = Path(__file__).resolve().parent
INPUT = 'lei_ren_part1_paper_candidate_five_bump_inverse.json'
WEIGHTS = 'lei_ren_part1_paper_bump_integral_enclosures_check.json'


def run(input_name=INPUT, output=None):
    raw = (HERE / input_name).read_bytes()
    source = json.loads(raw)
    for name, digest in source['input_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Candidate inverse dependency changed: ' + name)
    weight_data = json.loads((HERE / WEIGHTS).read_bytes())
    for name, digest in (
        ('lei_ren_part1_paper_bump_integral_enclosures.py', weight_data['module_sha256']),
        ('lei_ren_part1_paper_bump_integral_enclosures_check.py', weight_data['check_source_sha256']),
    ):
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Directed weight source changed: ' + name)
    if not weight_data['directed_integrals_certified']:
        raise ValueError('Weight receipt does not attest directed integrals')
    ctx = MPIntervalContext()
    ctx.dps = weight_data['settings']['precision']
    def restore_interval(packet):
        return ctx.mpf([mp.make_mpf(tuple(packet['lower_exact_mpf_tuple'])),
                        mp.make_mpf(tuple(packet['upper_exact_mpf_tuple']))])
    def fixed_atoms(atoms):
        return sum((ctx.mpf(mp.make_mpf(tuple(row['exact_mpf_tuple'])))
                    for row in atoms), ctx.mpf(0))
    def dual(packet):
        return fixed_atoms(packet['value']), fixed_atoms(packet['first_Z'])
    h = [dual(packet) for packet in source['materialized_controls']]
    amplitude, amplitude_Z = dual(source['amplitude'])
    defects = [dual(packet) for packet in source['source_defects']]
    centers = ('1.25', '1.5', '1.75')
    weights = {}
    for record in weight_data['weight_records'].values():
        weights[Fraction(record['center']), Fraction(record['power']), record['multiplicity']] = \
            restore_interval(record['weight_interval'])
    def weight(i, p, k=1):
        return weights[Fraction(centers[i]), Fraction(p), k]
    c = [h[0][0], h[1][0]]
    cz = [h[0][1], h[1][1]]
    xi = [packet[0] for packet in h[2:]]
    xiz = [packet[1] for packet in h[2:]]
    axial = (0, 2)
    zero = ctx.mpf(0)
    def total(terms):
        return sum(terms, zero)
    rows = [c[0] + c[1],
            total(weight(i, '.6') * c[j] + weight(i, '.5', 2) * c[j] * xi[i]
                  for j, i in enumerate(axial)),
            total(weight(i, '.5') * xi[i] for i in range(3)),
            -total(weight(i, '.1') * xi[i] for i in range(3))
            + total(weight(i, '0', 2) * c[j]**2 for j, i in enumerate(axial)) / amplitude**2
            - total(weight(i, '0', 2) * xi[i]**2 for i in range(3)) / 2,
            total(weight(i, '-.9') * xi[i] + weight(i, '-1', 2) * xi[i]**2 / 2
                  for i in range(3))]
    gg = total(weight(i, '0', 2) * c[j]**2 for j, i in enumerate(axial))
    rows_Z = [cz[0] + cz[1],
              total(weight(i, '.6') * cz[j]
                    + weight(i, '.5', 2) * (cz[j] * xi[i] + c[j] * xiz[i])
                    for j, i in enumerate(axial)),
              total(weight(i, '.5') * xiz[i] for i in range(3)),
              -total(weight(i, '.1') * xiz[i] for i in range(3))
              + total(2 * weight(i, '0', 2) * c[j] * cz[j]
                      for j, i in enumerate(axial)) / amplitude**2
              - 2 * gg * amplitude_Z / amplitude**3
              - total(weight(i, '0', 2) * xi[i] * xiz[i] for i in range(3)),
              total(weight(i, '-.9') * xiz[i] + weight(i, '-1', 2) * xi[i] * xiz[i]
                    for i in range(3))]
    def packet(interval):
        lo, hi = (mp.make_mpf(v) for v in interval._mpi_)
        with mp.workdps(90):
            return dict(lower=mp.nstr(lo, 85), upper=mp.nstr(hi, 85),
                        width=mp.nstr(hi-lo, 85), contains_zero=bool(lo <= 0 <= hi))
    residual = [rows[i] + defects[i][0] for i in range(5)]
    residual_Z = [rows_Z[i] + defects[i][1] for i in range(5)]
    inputs = (input_name, WEIGHTS, 'lei_ren_part1_paper_bump_integral_enclosures.py',
              'lei_ren_part1_paper_bump_integral_enclosures_check.py',
              'lei_ren_part1_paper_interval_taylor.py')
    report = dict(center=source['center'], input_inverse=input_name, finite_controls_fixed=True,
                  nominal_source_targets_fixed=True,
                  actual_exponential_bump_integrals_enclosed=True,
                  value_residuals=[packet(v) for v in residual],
                  first_Z_residuals=[packet(v) for v in residual_Z],
                  all_value_residuals_include_zero=all(packet(v)['contains_zero'] for v in residual),
                  all_first_Z_residuals_include_zero=all(packet(v)['contains_zero'] for v in residual_Z),
                  normalization=packet(restore_interval(weight_data['normalization'])),
                  directed_weight_settings=weight_data['settings'],
                  persisted_directed_weights_reused=True,
                  weight_enclosures={str(key): packet(v) for key, v in weights.items()},
                  source_error_enclosed=False, coefficient_ring_remainder_enclosed=False,
                  functional_terminal_closure=False, temporal_recursion=False,
                  input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                for name in inputs},
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    destination = Path(output) if output is not None else Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('Directed candidate integral residual generated; zero included in every value row:',
          report['all_value_residuals_include_zero'], flush=True)
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=INPUT)
    parser.add_argument('--output')
    args = parser.parse_args()
    run(args.input, args.output)
