"""Same-datum Lambda120 core seeds and finite rows at new axial centers.

Every center or interval recomputes the exact amplitude primitive, analytic
pressure and coupled coefficient rows. Frozen Z=.3 tensors are never read.
A center interval encloses a family of local series, not an extrapolation.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_exact_amplitude import exact_amplitude, _encode
from lei_ren_part1_paper_candidate_pressure_function import load_datum, pressure_jets, ACCEPTED_SHA
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_factored_core_positivity import squared_axis_rows
from lei_ren_part1_paper_functional_core_recursion import coupled_rows, _source_hashes
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).resolve().parent
PARAMETERS = dict(j='1e-14', Lambda='1e120', logC='5e151', delta='1e-200', sigma_denominator='500')


class GeneralCenterCoreFactory:
    def __init__(self, precision=260):
        self.precision = precision
        self.ctx = MPIntervalContext()
        self.ctx.dps = precision
        self.datum, self.pressure_hashes = load_datum()

    def build(self, center, degree=8, required_depth=3):
        ctx = self.ctx
        length = degree + required_depth + 1
        with mp.workdps(self.precision+40):
            z = ctx.mpf(center)
            lo, hi = endpoints(z)
            if lo < -1 or hi > 1:
                raise ValueError('Axial center interval must lie in [-1,1]')
            amplitude = exact_amplitude(ctx, Lambda=PARAMETERS['Lambda'], Z=z)
            axis = uniform_axis_jets(ctx, radius=1, j=PARAMETERS['j'],
                                     Lambda=PARAMETERS['Lambda'], logC=PARAMETERS['logC'],
                                     delta=PARAMETERS['delta'], length=length, axial_interval=z)
            pressure = pressure_jets(ctx, z, order=length-1, datum=self.datum)
            if pressure['accepted_schedule_sha256'] != ACCEPTED_SHA:
                raise ValueError('Accepted pressure datum changed')
            gradient = list(axis['gradient_coefficients'])
            f0 = amplitude['F0_interval']
            coarse_lo, coarse_hi = endpoints(axis['F0_interval'])
            exact_lo, exact_hi = endpoints(f0)
            if not coarse_lo <= exact_lo <= exact_hi <= coarse_hi:
                raise ArithmeticError('Exact amplitude outside uniform analytic amplitude enclosure')
            lam = ctx.mpf(PARAMETERS['Lambda'])
            fixed = dict(ell_Z_taylor=[-lam*g for g in gradient],
                         S_Z_taylor=squared_axis_rows(ctx, gradient, lam, f0, length),
                         U0_Z_taylor=list(axis['U0']),
                         P0_Z_taylor=list(pressure['physical_pressure_coefficients']),
                         gradient_coefficients=gradient, F0_interval=f0)
            rows = coupled_rows(ctx, fixed, z, degree, PARAMETERS['delta'],
                                required_depth=required_depth)
            if rows['axis_jet_length'] != length or any(len(row) != length-n
                                                       for rowset in (rows['A'], rows['Uz'], rows['P'])
                                                       for n, row in enumerate(rowset)):
                raise ArithmeticError('Axial depth was lost in the coupled recursion')
            return dict(requested_center=center, center_Z=z, fixed=fixed, core_rows=rows,
                        amplitude={name: amplitude[name] for name in
                                   ('G_interval', 'logF0_interval', 'F0_interval', 'checks')},
                        pressure=pressure, exact_amplitude_within_coarse_enclosure=True,
                        same_accepted_pressure_datum=True,
                        requested_center_recomputed=True, old_center_tensor_loaded=False,
                        source_parameter_errors_enclosed=False, infinite_radial_remainder_enclosed=False,
                        comparison_transition_generated=False, terminal_functional_closure=False,
                        temporal_recursion=False, full_field=False)


def contains_family(family, point):
    count = 0
    for name in ('A', 'Uz', 'P'):
        for wide_row, narrow_row in zip(family['core_rows'][name], point['core_rows'][name]):
            for wide, narrow in zip(wide_row, narrow_row):
                lo, hi = endpoints(wide)
                pl, ph = endpoints(narrow)
                if not lo <= pl <= ph <= hi:
                    raise ArithmeticError('New scalar center escaped interval family: '+name)
                count += 1
    return count


def run():
    factory = GeneralCenterCoreFactory()
    cases = {}
    for label, center in (('new_scalar_Z05', '.5'), ('local_family_Z049_Z051', ['.49', '.51'])):
        cases[label] = factory.build(center)
        print('Same-source new-center core generated:', label, 'radial degree8, retained axial depth3', flush=True)
    count = contains_family(cases['local_family_Z049_Z051'], cases['new_scalar_Z05'])
    inputs = ('lei_ren_part1_paper_candidate_exact_amplitude.py',
              'lei_ren_part1_paper_candidate_pressure_function.py',
              'lei_ren_part1_paper_uniform_axis_jets.py',
              'lei_ren_part1_paper_factored_core_positivity.py')
    hashes = _source_hashes()
    hashes.update({name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in inputs})
    hashes.update(factory.pressure_hashes)
    with mp.workdps(factory.precision+40):
        report = dict(parameters=PARAMETERS, precision=factory.precision,
                      radial_degree=8, retained_axial_depth=3,
                      cases=cases, scalar_coefficient_intervals_contained_in_local_family=count,
                      accepted_schedule_sha256=ACCEPTED_SHA,
                      all_fourteen_accepted_pressure_stages_retained=True,
                      old_center_tensor_used=False, whole_axis_core_generated=False,
                      comparison_transition_generated=False, temporal_recursion=False,
                      infinite_radial_remainder_enclosed=False,
                      input_hashes=hashes,
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(report), indent=2)+'\n', encoding='utf-8')
    print('New-center/local-family core receipt generated; coefficient inclusion count', count, flush=True)
    return report


if __name__ == '__main__':
    run()
