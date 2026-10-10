"""Independent rational P0 derivatives and both live numeric source owners."""
import copy
import gzip
import json
import math
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_pressure_binding as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError):
        return True
    raise AssertionError('Invalid current pressure source admitted')


def polynomial_at(ctx, polynomial, z):
    result = ctx.mpf(0)
    for coefficient in s.Poly(polynomial, s.Symbol('Z', real=True)).all_coeffs():
        p, q = coefficient.as_numer_denom()
        result = result*z+ctx.mpf(int(p))/int(q)
    return result


def rational_derivative(ctx, order, z):
    symbol = s.Symbol('Z', real=True)
    # Independently differentiate; no use of the production q_jets recurrence.
    expression = s.cancel(s.diff((1+symbol**2)**-2, symbol, order)/math.factorial(order))
    numerator, denominator = expression.as_numer_denom()
    return polynomial_at(ctx, numerator, z)/polynomial_at(ctx, denominator, z)


@source_precision
def run(owner, deliveries):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family'] == owner.family_record and all(owner.assert_graph().values())
    # Prior dependency hashes belong to the already accepted live source. New
    # blobs are checked here; publication audits unchanged Git blob identities.
    assert all(raw['input_hashes'][name] == digest for name, digest in owner.before.hashes.items())
    for name in (current.signed.NAME, current.signed.RECEIPT, Path(current.__file__).name):
        assert current.sha(name) == raw['input_hashes'][name]
    c = MPIntervalContext()
    c.dps = 400
    findings = {}
    row_count = 0
    for name, delivery in deliveries.items():
        result = owner.source_binding(delivery)
        assert current.correlated.report(result) == raw['actual_P0_source_bindings'][name]
        _, _, request = owner.amplitude.before._validate_delivery(delivery)
        # The exact rational coordinate is independently enclosed at 400 digits.
        z = c.mpf(request.Z.numerator)/request.Z.denominator
        checks = []
        for datum_name, datum in owner.datums.items():
            source = result['actual_live_owners'][datum_name]
            m0, m2 = c.mpf(datum.m0), c.mpf(datum.m2)
            for order, row in enumerate(source['normalized_P0_Taylor_rows']):
                q = rational_derivative(c, order, z)
                flat = c.mpf(source['flatten_Taylor_enclosures'][order])
                reference = -(m2*q+(m0 if order == 0 else 0)+flat)
                assert current.amplitude.contains(row, reference), (name, datum_name, order)
                assert current.amplitude.contains(source['normalized_q_Taylor_rows'][order], q)
                assert result['actual_inlet_P0'][order]._mpi_ == owner.inlet.ctx.mpf(current.ends(row))._mpi_
                if request.Z == 0 and order % 2:
                    assert current.ends(row) == current.ends(flat) == (0, 0)
                elif order == 0:
                    assert current.ends(flat)[1] > 0
                else:
                    assert current.ends(flat)[0] < 0 < current.ends(flat)[1]
                row_count += 1
            assert len(datum.stages) == 14 and current.ends(source['positive_late_tail'])[0] > 0
            checks.append(dict(owner=datum_name, six_independent_rational_derivatives_included=True,
                mass_uncertainty_and_flatten_Cauchy_error_retained=True, positive_late_tail_retained=True))
        for n in range(6):
            expected = result['actual_inlet_P0'][n]+result['separate_inlet_Mp'][n]
            assert expected._mpi_ == result['actual_inlet_pressure'][n]._mpi_
        lp = c.exp(c.mpf(40))+11
        assert current.amplitude.contains(result['directed_pressure_scale_log'], 2*lp)
        assert result['physical_Pstar_squared_not_expanded'] and not result['full_physical_accuracy']
        packet = owner.velocity_pressure(delivery)
        assert set(packet['values']) == {'u', 'v', 'w', 'p'}
        assert packet['current_P0_numeric_inclusion_pending'] and not packet['full_certified_physical_accuracy']
        assert packet['values']['p']['current_P0_numeric_inclusion_pending']
        assert not packet['values']['p']['ordinary_pressure_materialization_not_a_current_P0_binding']
        assert not any(packet[key] for key in current.OPEN)
        findings[name] = dict(exact_Z=str(request.Z), owners=checks,
            pressure_P0_plus_Mp_not_double_counted=True, original_exact_scale_log_included=True,
            actual_supported_uvw_pressure_interface_preserved=True)
    # One actual coordinate box checks both owners and parity behavior off a
    # point. This is an axial P0 source box, not a uniform 3D field certificate.
    box = owner.normalized_jets(['.370', '.372'])
    for exact in ('.370', '.371', '.372'):
        point = owner.normalized_jets(exact)
        for datum_name in owner.datums:
            for outer, inner in zip(box['actual_live_owners'][datum_name]['normalized_P0_Taylor_rows'],
                                    point['actual_live_owners'][datum_name]['normalized_P0_Taylor_rows']):
                assert current.amplitude.contains(outer, inner)
    bad = dict(copied_delivery=rejected(lambda: owner.source_binding(copy.copy(next(iter(deliveries.values()))))),
        caller_pressure_oracle=rejected(lambda: owner.normalized_jets({'P0': [0, 0]})),
        out_of_support_Z=rejected(lambda: owner.normalized_jets(['-2', '0'])))
    datum = owner.datums['selected']
    saved = datum.m2
    try:
        datum.m2 = datum.ctx.mpf(0)
        bad['mutated_selected_mass'] = rejected(owner.assert_graph)
    finally:
        datum.m2 = saved
    saved = owner.selected.datum
    try:
        owner.selected.datum = owner.inlet.datum
        bad['substituted_selected_datum'] = rejected(owner.assert_graph)
    finally:
        owner.selected.datum = saved
    saved = owner.scale_log
    try:
        owner.scale_log = owner.graph.zero
        bad['substituted_pressure_scale'] = rejected(owner.assert_graph)
    finally:
        owner.scale_log = saved
    assert all(bad.values()) and all(owner.assert_graph().values())
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family_record,
        actual_live_P0_source_rows_independently_included=row_count,
        independent_directed_precision=400, source_findings=findings,
        actual_axial_coordinate_box_both_owner_rows_enclosed=True,
        actual_analytic_P0_six_Taylor_rows_and_q_recurrence_identified=True,
        distinct_inlet_and_selected_datum_owners_preserved=True,
        fourteen_atoms_flatten_Cauchy_errors_and_positive_late_tail_retained=True,
        original_pressure_scale_log_bound_without_expansion=True,
        source_graph_assertions=owner.assert_graph(), invalid_inputs_rejected=bad,
        P0_source_inclusion_does_not_certify_absolute_physical_accuracy=True,
        **dict.fromkeys(current.GATES, True), **dict.fromkeys(current.OPEN, False),
        input_hashes=hashes, execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt), indent=2)+'\n', encoding='utf8', newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_P0_BINDING both live owners, six rows, original scale and positive tails', flush=True)
    return receipt
