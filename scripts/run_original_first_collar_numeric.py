"""Execute the original first-collar numerical approximation and saved receipt."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_numeric import OriginalFirstCollarNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY = ROOT/'experiments/root_st073'
NAME = 'lei_ren_part1_paper_compliant_current_original_first_collar_numeric.json'


def split_identities():
    u, e, f, p, h, m, k, a, b, cc = sy.symbols('U0 epsilon phi psi H Mpsi Kpsi Apsi B C')
    dh, dm, dk, da, db, dc = 2*f-2*h, p-m, 2*f*p-2*k, p*p-a, f*f-2*b, f*f-cc
    V, M, K, A = u+e*p, u+e*m, u*h+e*k, u*u+2*e*u*m+e*e*a
    equations = (dh-(2*f-2*h), e*dm-(V-M), u*dh+e*dk-(2*f*V-2*K),
                 2*e*u*dm+e*e*da-(V*V-A), db-(f*f-2*b), dc-(f*f-cc))
    for row in equations:
        assert sy.expand(row) == 0
    return len(equations)


def quadrature_inlets(field, Z, degree):
    """Independent Gauss integration of products, not coefficient primitives."""
    c = field.ctx
    packet = field.core.coefficients(Z, degree=degree, depth=1)
    rows = packet['point']
    phi = [row[0] for row in rows['A']]
    psi = [c.mpf(0)]+[rows['Uz'][n][0]/field.core.epsilon for n in range(1, degree+1)]
    nodes, weights = c.gauss_quadrature(degree+2, 'legendre')
    values = {name:c.mpf(0) for name in ('H', 'Mpsi', 'Kpsi', 'Apsi', 'B', 'C')}
    for x, w in zip(nodes, weights):
        r, weight = 2*(x+1), 2*w
        f, p = c.polyval(list(reversed(phi)), r), c.polyval(list(reversed(psi)), r)
        integrands = dict(H=r*f/8, Mpsi=p/4, Kpsi=r*f*p/8, Apsi=p*p/4, B=r*f*f/16, C=f*f/4)
        for name, integrand in integrands.items():
            values[name] += weight*integrand
    expected = field.inlet(Z, degree=degree)['point']
    for name, value in values.items():
        assert abs(value-expected[name]) <= c.mpf('1e-60')*max(1, abs(value), abs(expected[name])), name
    return len(values)


def run():
    started = time.perf_counter()
    field = OriginalFirstCollarNumeric(precision=80)
    samples, refinements, maxima = [], [], dict(Phi=field.ctx.mpf(0), Psi=field.ctx.mpf(0))
    quadrature_checks = split_checks = 0
    splits = split_identities()
    for Z in ('-.3', '0', '.3', '.17'):
        quadrature_checks += quadrature_inlets(field, Z, 22)
        for phase in ('0', '.25', '.5', '1'):
            fine = field.evaluate(Z, phase, degree=22)
            coarse = field.evaluate(Z, phase, degree=16)
            for name, value in fine['normalized_fields_and_own_histories'].items():
                for key in ('inlet', 'increment_per_original_hb'):
                    p, q = value[key], coarse['normalized_fields_and_own_histories'][name][key]
                    error = field.ctx.mpf(p['absolute_error_upper'])+field.ctx.mpf(q['absolute_error_upper'])
                    delta = abs(field.ctx.mpf(p['approximation'])-field.ctx.mpf(q['approximation']))
                    assert delta <= error, (Z, phase, name, key)
                    split_checks += 1
            if phase == '0':
                for row in fine['normalized_fields_and_own_histories'].values():
                    inc = row['increment_per_original_hb']
                    assert field.ctx.mpf(inc['approximation']) == 0
                    assert inc['source_enclosure_exact_mpf_tuples'] == [[0, 0, 0, 0], [0, 0, 0, 0]]
            else:
                # This is a source-controlled collar, whose leading angular
                # increment differs from analytic core continuation by wc/s.
                ratio = field.ctx.mpf(fine['original_flat_complement_mass'])/field.ctx.mpf(phase)
                assert 0 < ratio < 1
                for name in maxima:
                    row = fine['normalized_fields_and_own_histories'][name]['increment_per_original_hb']
                    value = field.ctx.mpf(row['approximation'])
                    if value and Z != '0':
                        relative = field.ctx.mpf(row['absolute_error_upper'])/abs(value)
                        maxima[name] = max(maxima[name], relative)
                        assert relative < field.ctx.mpf('1e-4'), (Z, phase, name, relative)
            samples.append(fine)
            refinements.append(dict(Z=Z, phase=phase, directed_coarse_fine_comparisons=16))
        print('Computed original actual first collar, inlet through s=1: Z='+Z, flush=True)
    rejected = 0
    for args in (('.3', '-.1'), ('.3', '1.01'), ('1', '.5'), ('nan', '.5'), ('.3', .5)):
        try:
            field.evaluate(*args)
        except (ValueError, ArithmeticError):
            rejected += 1
        else:
            raise AssertionError(('unsupported input accepted', args))
    code = dict(field.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,
        checks=dict(passed=True, exact_original_six_history_split_identities=splits,
            independent_Gauss_core_inlet_integrals=quadrature_checks,
            original_first_collar_sample_packets=len(samples), original_axial_points=4,
            directed_degree16_degree22_comparisons=split_checks, exact_zero_phase_included=True,
            actual_first_exit_included=True, offgrid_original_inlet_included=True,
            maximum_nonaxis_field_increment_relative_error={k:field.ctx.nstr(v, 35) for k,v in maxima.items()},
            original_width_product_log_proofs=len(field.scale_proofs), rejected_unsupported_inputs=rejected,
            actual_source_used_instead_of_comparison_histories=True,
            all_pressure_swirl_hydro_drive_sources_enclosed=True, source_midpoint_selected=False),
        samples=samples, radial_refinement_comparisons=refinements,
        normalized_source_product_log_proofs=field.scale_proofs,
        input_hashes=field.hashes, implementation_hashes=code,
        limitations=dict(original_leading_first_collar_values_with_directed_absolute_errors=True,
            nominal_first_width_coefficients_only=True, exact_original_point_oracle=False,
            tiny_cancellation_in_normalized_Kpsi_increment_not_resolved=True,
            original_hb_F0_materialized=False, axial_derivative_numeric_transport=False,
            physical_collar_velocity=False, n1_collar_transport=False, global_moment_repair=False,
            global_NS_field=False, temporal_n_ge_2_recursion=False),
        execution_seconds=round(time.perf_counter()-started, 3))


def inspect(path):
    report = json.loads(path.read_text(encoding='utf-8'))
    assert report['checks']['passed']
    for name, digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest, name
    for name, digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-report', action='store_true')
    parser.add_argument('--inspect-report', action='store_true')
    args = parser.parse_args()
    path = DIRECTORY/NAME
    result = inspect(path) if args.inspect_report else run()
    if args.write_report:
        path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result['checks'], indent=2), flush=True)
