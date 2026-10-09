"""Audit the live source coverage, original equations and frame independence."""
import gzip
import json
import math
import time
from pathlib import Path
from unittest.mock import patch

import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source as current

ep, read = current.ep, current.read


def forbidden(*args, **kwargs):
    raise AssertionError('An old source constructor or saved frame owner was invoked')


def symbolic_histories():
    y = sy.symbols('y', real=True)
    phi, V = sy.symbols('phi V', real=True)
    rates = dict(H=2, M=1, K=2, A=1, B=2, C=1)
    targets = dict(H=phi, M=V, K=phi*V, A=V*V, B=phi*phi/2, C=phi*phi)
    drives = dict(H=2*phi, M=V, K=2*phi*V, A=V*V, B=phi*phi, C=phi*phi)
    for key, rate in rates.items():
        inlet = sy.symbols(key+'0')
        value = sy.exp(-rate*y)*inlet+(1-sy.exp(-rate*y))*targets[key]
        assert sy.simplify(sy.diff(value, y)+rate*value-drives[key]) == 0
        assert sy.simplify(value.subs(y, 0)-inlet) == 0
    # Independent radial normalization of each of the six moments.
    rho = sy.symbols('rho', nonnegative=True)
    p, v = 1+rho/7, 2-rho/9
    polynomials = dict(H=(p, 1, 8), M=(v, 0, 4), K=(p*v, 1, 8),
                       A=(v*v, 0, 4), B=(p*p, 1, 16), C=(p*p, 0, 4))
    c = current.macro.MPIntervalContext()
    c.dps = 100
    phi_rows = [[c.mpf(1)], [c.mpf(1)/7]]
    V_rows = [[c.mpf(2)], [-c.mpf(1)/9]]
    actual = current.bridge.atom_module.finite_atom_coefficients(c, phi_rows, V_rows, 0)
    for key, (polynomial, power, normalization) in polynomials.items():
        exact = sy.integrate(polynomial*rho**power, (rho, 0, 4))/normalization
        value = c.mpf(int(exact.p))/int(exact.q)
        lo, hi = ep(actual[key][0])
        vl, vh = ep(value)
        assert lo <= vh and vl <= hi
    return dict(original_all_six_own_ODEs_and_core_inlet_identities=True,
                independent_exact_radial_atom_normalizations=6)


def same_source(actual, saved, counts):
    """Require all emitted source interval branches to match a live replay."""
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and set(actual) == set(saved)
        if 'lower_exact_mpf_tuple' in actual:
            assert actual['lower_exact_mpf_tuple'] == saved['lower_exact_mpf_tuple']
            assert actual['upper_exact_mpf_tuple'] == saved['upper_exact_mpf_tuple']
            counts['exact_directed_interval_rows'] += 1
            return
        for key in actual:
            same_source(actual[key], saved[key], counts)
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(actual) == len(saved)
        for left, right in zip(actual, saved):
            same_source(left, right, counts)
    else:
        assert actual == saved


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] is True
    assert all(saved[key] is False for key in current.OPEN)
    assert saved['exact_Z_domain'] == ['-1', '1']
    assert saved['exact_Z_partition'] == [list(cell) for cell in current.CELLS]
    assert len(saved['source_cells']) == len(current.CELLS)
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Leading source input changed: '+name
    guards = [
        patch.object(current.macro.OriginalBridgeMacroFunctions, 'owner', forbidden),
        patch.object(current.bridge.BridgeComparison, '__init__', forbidden),
        patch.object(current.bridge.BridgeCoreAtoms, '__init__', forbidden),
        patch.object(current.bridge.CompliantActualBridgeIntegrals, '__init__', forbidden),
        patch.object(current.rebuild_module.CompliantCoreCoefficientRebuild, '__init__', forbidden),
        patch.object(current.amplitude_module.CompliantAnchoredAxisAmplitude, '__init__', forbidden),
        patch.object(current.core_module.CompliantCorePhysicalField, '__init__', forbidden),
    ]
    counts = dict(exact_directed_interval_rows=0, whole_Z_cells=0,
                  fresh_radial_coefficients=0, rejected_domains=0)
    for guard in guards:
        guard.start()
    try:
        with mp.workdps(540):
            provider = current.OriginalWholeZBridgeSource()
            assert provider.source == saved['implicit_source_sha256']
            assert provider.family == saved['actual_five_defect_family_sha256']
            assert provider.datum == saved['datum_enclosure_sha256']
            # Remove the old labels and coefficient state before computation.
            # Successful production now cannot depend on any of those values.
            provider.admission.records['actual_bridge_integrals'].pop('packets')
            provider.admission.records['actual_bridge_integrals'].pop('whole_axis_R100')
            provider.admission.records['comparison_point_integrals'].pop('comparison_point_packets')
            provider.admission.records['anchored_axis_amplitude'].pop('anchored_amplitude_packets')
            provider.admission.records['anchored_axis_amplitude'].pop('physical_core_value_packets')
            c = provider.c
            previous_right = mp.mpf(-1)
            for ends, expected in zip(current.CELLS, saved['source_cells']):
                Z = provider.cell(ends)
                assert ep(Z)[0] == previous_right
                previous_right = ep(Z)[1]
                actual = provider.evaluate(ends)
                same_source(current.bridge._encode(current.serialized(actual)), expected, counts)
                flow, proof = provider.owner(ends)
                assert flow.c is c and flow.logs[0]._mpi_ == provider.logh._mpi_
                assert ep(flow.weight)[0] > 0
                assert flow.sources_set and ep(flow.logs[2])[1] < 0
                assert ep(proof['positive_nonlinear_swirl_bound_retained'])[0] > 0
                assert ep(proof['actual_anchored_amplitude']['G'])[0] >= 0
                assert ep(actual['live_macro_R100']['geometry']['R1'].coefficient) == (100, 100)
                for key in current.OPEN:
                    assert actual[key] is False
                for row in actual['live_macro_R100']['phi']+actual['live_macro_R100']['V']:
                    assert row.scale.bases is flow.logs and row.ledger is flow.ledger
                    assert all(mp.isfinite(x) for x in ep(row.coefficient))
                core = proof['fresh_actual_core']
                data = provider.actual.prepare(Z)['data']
                assert set(core['actual_own_six_moments_axial5']) == set(current.bridge.RATES)
                assert ep(core['actual_phi_axial5'][0])[0] > 0
                for name in current.bridge.RATES:
                    for n, value in enumerate(core['actual_own_six_moments_axial5'][name]):
                        assert value._mpi_ == data['moments'][name][n]._mpi_
                for name in ('actual_delta_phi_axial5', 'actual_delta_V_axial5'):
                    assert all(ep(row) == (0, 0) for row in core[name])
                radial = provider.rebuild.rebuild(Z, 24, 6)
                assert radial['old_finite_coefficient_rows_read'] is False
                assert radial['temporal_recursion'] is False
                assert len(radial['rows']['A']) >= 25
                counts['fresh_radial_coefficients'] += 25*7
                counts['whole_Z_cells'] += 1
            assert previous_right == 1
            for ends in (('-1.5','0'), ('0','1.5'), ('.5','-.5'), ('0',)):
                try:
                    provider.owner(ends)
                except ValueError:
                    counts['rejected_domains'] += 1
                else:
                    raise AssertionError('Invalid original source domain accepted')
    finally:
        for guard in reversed(guards):
            guard.stop()
    symbolic = symbolic_histories()
    result = dict(all_passed=True, **{current.GATE: True},
        implicit_source_sha256=provider.source,
        actual_five_defect_family_sha256=provider.family,
        datum_enclosure_sha256=provider.datum,
        full_real_Z_partition_checked=True,
        old_source_constructors_and_frame_owners_rejected_during_replay=True,
        saved_frame_and_archived_whole_axis_packets_removed_before_replay=True,
        leading_functions_are_not_current_finite_N_correction_histories=True,
        symbolic_source= symbolic, replay_counts=counts,
        **dict.fromkeys(current.OPEN, False),
        input_hashes={**provider.hashes, current.NAME: current.sha(current.NAME),
                     Path(__file__).name: current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    print('Live whole-Z leading bridge source checks passed', flush=True)
    return result


if __name__ == '__main__':
    run()
