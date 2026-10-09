"""Same-source all-Z micro/macro joins, positive proofs and five-ODE units."""
import gzip
import json
import math
from pathlib import Path
import time
from unittest.mock import patch

import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_micro_macro as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source

ep = current.ep


def forbidden(*args, **kwargs):
    raise AssertionError('Old constructor or frame owner invoked by live interval adapter')


def original_five_ODE_units():
    y = sy.symbols('y', real=True)
    F0, S = sy.symbols('F0 S', positive=True)
    R = sy.exp(y)
    phi, V = sy.Function('phi')(y), sy.Function('V')(y)
    old = {key: sy.Function(key)(y) for key in current.moments.RATES}
    E = sy.sqrt(2*R)*F0*phi/S
    H, M, K, A, B, C = (old[key] for key in ('H','M','K','A','B','C'))
    deriv = {sy.diff(H,y):2*phi-2*H, sy.diff(M,y):V-M,
             sy.diff(K,y):2*phi*V-2*K, sy.diff(A,y):V**2-A,
             sy.diff(B,y):phi**2-2*B, sy.diff(C,y):phi**2-C}
    five = dict(m=M, h=sy.sqrt(R)*F0*H/(S*sy.sqrt(2)),
                k=sy.sqrt(R)*F0*K/(S*sy.sqrt(2)),
                e=A/S**2-R*F0**2*B/S**2, p=R*F0**2*C/S**2)
    rhs = dict(m=V-five['m'], h=E-sy.Rational(3,2)*five['h'],
               k=E*V-sy.Rational(3,2)*five['k'],
               e=V**2/S**2-five['e']-E**2/2, p=E**2/2)
    for key, value in five.items():
        assert sy.simplify((sy.diff(value,y)-rhs[key]).subs(deriv)) == 0
    P0 = sy.symbols('independent_P0')
    assert sy.simplify(sy.diff(P0+five['p'],y).subs(deriv)-E**2/2) == 0
    return dict(passed=True, original_all_five_physical_histories_y_ODEs=True,
                analytic_P0_independent_of_cumulative_pressure=True)


def finite_rows(row, counts):
    if isinstance(row, dict):
        if 'lower_exact_mpf_tuple' in row:
            lo, hi = ep(current.read(counts['context'], row))
            assert mp.isfinite(lo) and mp.isfinite(hi) and lo <= hi
            counts['saved_finite_directed_interval_rows'] += 1
        else:
            for value in row.values():
                finite_rows(value, counts)
    elif isinstance(row, list):
        for value in row:
            finite_rows(value, counts)


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and all(saved[key] is False for key in current.OPEN)
    assert saved['exact_Z_partition'] == [list(cell) for cell in current.source.CELLS]
    assert saved['exact_Z_domain'] == ['-1','1']
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Current source dependency changed: '+name
    guards = [patch.object(current.micro.HydratedComparison, '__init__', forbidden),
        patch.object(current.micro.OriginalMicroFunctions, '__init__', forbidden),
        patch.object(current.macro.OriginalMacroFiniteN, '__init__', forbidden),
        patch.object(current.first.FirstSwitchFunctions, '__init__', forbidden),
        patch.object(current.endpoint.OriginalR100Endpoint, '__init__', forbidden),
        patch.object(current.moments.OriginalBridgeMacroMoments, '__init__', forbidden),
        patch.object(current.source.macro.OriginalBridgeMacroFunctions, 'owner', forbidden)]
    counts = dict(whole_Z_cells=0, saved_finite_directed_interval_rows=0,
                  exact_directed_interval_rows=0, live_replayed_radial_cells=0,
                  live_seam_overlap_rows=0, original_positive_Dbar_calls=0,
                  typed_source_rejections=0)
    for guard in guards:
        guard.start()
    try:
        with mp.workdps(540):
            owner = current.WholeZMicroMacro()
            assert owner.identity == saved['source_family']
            assert ep(owner.logDlower)[0] == ep(-owner.c.ln(2)-owner.logK)[0]
            # A live query must remain independent of every saved source frame.
            owner.source.admission.records['actual_bridge_integrals'].pop('packets')
            owner.source.admission.records['actual_bridge_integrals'].pop('whole_axis_R100')
            owner.source.admission.records['comparison_point_integrals'].pop('comparison_point_packets')
            owner.source.admission.records['anchored_axis_amplitude'].pop('anchored_amplitude_packets')
            owner.source.admission.records['anchored_axis_amplitude'].pop('physical_core_value_packets')
            counts['context'] = owner.c
            right = mp.mpf(-1)
            for ends, packet in zip(current.source.CELLS, saved['source_cells']):
                z = owner.source.cell(ends)
                assert ep(z)[0] == right
                right = ep(z)[1]
                assert packet['exact_Z_cell'] == list(ends)
                assert len(packet['whole_radial_source_cells']) == 6
                finite_rows(packet, counts)
                op = owner.owner(ends)
                assert op.micro.flow is op.flow and op.macro.flow is op.flow
                assert op.series.flow is op.flow and op.micro.decoder.flow is op.flow
                assert op.reference.P0 is op.axis['pressure_axis_over_Pstar_squared_axial5']
                seam = owner.seam(ends)
                assert seam['actual_micro_histories_are_the_live_macro_incoming']
                counts['live_seam_overlap_rows'] += seam['overlap_rows']
                for index, query in ((0, current.MICRO_CELLS[0]),
                                     (4, ('frozen_macro', *current.MACRO_CELLS[0]))):
                    chart, left, end = query
                    result = owner.generic(ends, chart, left, end)
                    same_source(current.source.bridge._encode(current.source.serialized(
                        current.exported(result))), packet['whole_radial_source_cells'][index], counts)
                    background = result['background']
                    proof = background['actual_positive_Dbar_proof']
                    assert proof['original_analytic_theorem'] is owner.positive_theorem
                    assert proof['strict_positive_complete_source_enclosure']
                    assert ep(proof['source_row'].coefficient)[0] > 0
                    assert proof['source_log_lower']._mpi_ == owner.logDlower._mpi_
                    assert ep(background['correlated_a_axial5'][0].coefficient)[0] > 0
                    generic = result['full_signed_generic_source']
                    assert generic['common_original_P0_axial5'] is op.reference.P0
                    assert generic['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
                    assert generic['actual_positive_E_C_and_generic_cone_admission_certified'] is False
                    assert set(generic['common_own_five_histories_axial5']) == {'m','h','k','e','p'}
                    assert len(generic['common_radial_Q_axial4']) == 5
                    for row in generic['common_velocity_E_axial5']:
                        assert row.scale.bases is op.flow.logs and row.ledger is op.flow.ledger
                    counts['live_replayed_radial_cells'] += 1
                    counts['original_positive_Dbar_calls'] += 1
                assert ep(op.axis['radius'].coefficient) == (100,100)
                assert op.axis['analytic_P0_kept_separate'] and op.axis['actual_moments_not_comparison']
                counts['whole_Z_cells'] += 1
                print('Whole-Z live consumer checks: '+str(ends),flush=True)
            assert right == 1
            op = owner.owner(current.source.CELLS[0])
            # A crossing-zero cover is sharpened by the bound, not by choosing
            # its endpoint. A nonpositive cover contradicting the theorem fails.
            crossing = op.flow.scalar(owner.c.mpf([-1,4]))
            proof = owner.positive_Dbar(op.flow,crossing,'actual_micro_comparison_Dbar')
            assert proof['original_raw_source_row'] is crossing
            assert ep(proof['source_row'].coefficient) == (1,1)
            assert ep(proof['source_row'].scale.offset)[0] == ep(owner.logDlower)[0]
            bad = [(op.flow.scalar(-1),'actual_micro_comparison_Dbar'),
                   (op.flow.scalar(1),'arbitrary_other_denominator'),
                   (current.source.macro.prior.ScaledEnclosure(crossing.scale,
                    crossing.coefficient,{}),'actual_micro_comparison_Dbar')]
            for row, role in bad:
                try:
                    owner.positive_Dbar(op.flow,row,role)
                except (ValueError, ArithmeticError):
                    counts['typed_source_rejections'] += 1
                else:
                    raise AssertionError('Wrong Dbar source/theorem accepted')
    finally:
        for guard in reversed(guards):
            guard.stop()
    counts.pop('context')
    result = dict(all_passed=True, **{current.GATE:True}, source_family=owner.identity,
        source_function_bindings=owner.bindings, independent_own_history_units=original_five_ODE_units(),
        replay_counts=counts, full_Z_source_cell_coverage_and_saved_finiteness_checked=True,
        two_live_radial_cells_and_live_micro_macro_seam_replayed_per_Z_cell=True,
        leading_history_and_signed_generic_source_not_finite_N_correction=True,
        old_constructor_and_saved_frame_owner_calls_rejected=True,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
                     Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z live micro/macro source and history checks passed',flush=True)
    return result


if __name__ == '__main__':
    run()
