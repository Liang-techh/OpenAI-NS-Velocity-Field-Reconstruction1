"""Focused checker for the pre-``R_p`` C1 input archive."""
import gzip
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_preRp_C1_inputs as current


def _read_raw():
    return json.loads(gzip.decompress(
        (current.HERE / current.NAME).read_bytes()))


def _independent_compact_fixture():
    """Check the differentiated radial-free identity at nonunit ``R_c``."""
    z, x = sy.symbols("z x", real=True)
    delta = sy.Rational(1, 5)
    Rc = sy.Rational(37, 10)
    M = 1 + z / 5 + z ** 2 / 7
    m = M * x
    U = 2 * M * x                    # m_y = U - m, since x*m_x = m
    d = 1 - z ** 2
    L = 1 - delta * z ** 2
    m_z = sy.diff(m, z)
    q = sy.cancel((2 * z * U - (1 - delta) * z * m - d * m_z) / L)
    q_y = x * sy.diff(q, x)
    S = ((q + (1 - delta) * z * sy.diff(q, z) / 2 + q_y) / L
         + q * (q / 2 + q_y)
         + U * (-2 * z * q + d * sy.diff(q, z) - 2 * z * q_y) / L)
    I0 = sy.integrate(S * Rc, (x, 1, 2))
    I1 = sy.integrate(S * Rc ** 2 * x, (x, 1, 2))
    def B(expr):
        return sy.expand(expr.subs(x, 2) - expr.subs(x, 1))
    moments0 = {
        "q_Z": sy.integrate(sy.diff(q, z) * Rc, (x, 1, 2)),
        "U_q_Z": sy.integrate(U * sy.diff(q, z) * Rc, (x, 1, 2)),
        "U_squared": sy.integrate(U ** 2 * Rc, (x, 1, 2)),
        "U_U_Z": sy.integrate(U * sy.diff(U, z) * Rc, (x, 1, 2)),
    }
    moments1 = {
        "q": sy.integrate(q * Rc ** 2 * x, (x, 1, 2)),
        "q_Z": sy.integrate(sy.diff(q, z) * Rc ** 2 * x, (x, 1, 2)),
        "q_squared": sy.integrate(q ** 2 * Rc ** 2 * x, (x, 1, 2)),
        "U_q_Z": sy.integrate(U * sy.diff(q, z) * Rc ** 2 * x, (x, 1, 2)),
        "U_squared": sy.integrate(U ** 2 * Rc ** 2 * x, (x, 1, 2)),
        "U_U_Z": sy.integrate(U * sy.diff(U, z) * Rc ** 2 * x, (x, 1, 2)),
    }
    boundaries = dict(
        R_q=B(Rc * x * q),
        R_q_squared=B(Rc * x * q ** 2),
        R_U_squared=B(Rc * x * U ** 2),
        R_squared_q=B(Rc ** 2 * x ** 2 * q),
        R_squared_q_squared=B(Rc ** 2 * x ** 2 * q ** 2),
        R_squared_U_squared=B(Rc ** 2 * x ** 2 * U ** 2),
    )
    a = 1 - delta
    S0 = (boundaries["R_q"] / L + a * z * moments0["q_Z"] / (2 * L)
          + boundaries["R_q_squared"] / 2
          + d * moments0["U_q_Z"] / L
          - 2 * z ** 2 * boundaries["R_U_squared"] / L ** 2
          - 2 * delta * z ** 2 * moments0["U_squared"] / L ** 2
          + 2 * z * d * moments0["U_U_Z"] / L ** 2)
    S1 = (boundaries["R_squared_q"] / L - moments1["q"] / L
          + a * z * moments1["q_Z"] / (2 * L)
          + boundaries["R_squared_q_squared"] / 2
          - moments1["q_squared"] / 2
          + d * moments1["U_q_Z"] / L
          - 2 * z ** 2 * boundaries["R_squared_U_squared"] / L ** 2
          + 2 * a * z ** 2 * moments1["U_squared"] / L ** 2
          + 2 * z * d * moments1["U_U_Z"] / L ** 2)
    residual0 = sy.factor(sy.together(I0 - S0))
    residual1 = sy.factor(sy.together(I1 - S1))
    derivative0 = sy.factor(sy.together(sy.diff(I0, z) - sy.diff(S0, z)))
    derivative1 = sy.factor(sy.together(sy.diff(I1, z) - sy.diff(S1, z)))
    if any(value != 0 for value in (residual0, residual1, derivative0, derivative1)):
        raise AssertionError("independent compact radial-free identity failed")
    return dict(
        passed=True, Rc="37/10", delta="1/5", residual_m13="0",
        residual_m15="0", derivative_residual_m13="0",
        derivative_residual_m15="0", nonunit_endpoint_fixture=True)


def _assert_scratch_and_baseline(owner):
    # Re-read the retained checked physical source directly.  The exploratory
    # scratch atlas is optional and is deliberately not a publication input.
    fields = owner.fields
    fields.assert_graph()
    bd = current.rg.Bounds(owner.c)
    import math
    for key, view in owner.baseline.items():
        ends = tuple(view["exact_Z_cell"])
        chart = view["chart"]
        partition = tuple(fields.partitions[chart])
        source_view = current.physical.source_rows(
            fields, ends, chart, partition[0], partition[-1])
        packet = source_view["source_packet"]
        assert packet["actual_phase_Z_exact_zero"] is True
        assert source_view["Pstar"] is fields.modified.ranges.phase.outer.owner.owner(ends).Pstar
        V = current.rg.JetBound(*(
            bd.row(source_view["V"][k] * math.factorial(k))
            for k in range(4)))
        M = current.rg.JetBound(*(
            bd.row(source_view["m"][k] * math.factorial(k))
            for k in range(4)))
        Pstar = bd.row(source_view["Pstar"])
        U = bd.scaled(V, Pstar)
        M = bd.scaled(M, Pstar)
        assert view["exact_native_endpoints"] == [partition[0], partition[-1]]
        assert view["source_cell_count"] == len(partition) - 1
        assert view["physical_m_Z3_magnitude"] == [item.record() for item in current.rg.rows(M)]
        assert view["physical_U_Z3_magnitude"] == [item.record() for item in current.rg.rows(U)]
        assert view["Pstar_record"] == Pstar.record()
    fields.assert_graph()


def _assert_scale_and_compact(owner):
    for power, key in ((1, "m13"), (2, "m15")):
        expected = current._record(
            current.preRp.Rc_unit_log(owner.c, owner.pulse, power))
        assert current._record(owner.scale_ratios[key]) == expected
    assert len(owner.compact) == 4
    for key, value in owner.compact.items():
        assert value["baseline_first_Z_derivative_exact_zero"]
        assert value["physical_common_scale_Pstar_applied_once"]
        assert value["selected_N_rows_source_bound"]
        assert value["selected_N_division_not_repeated"]
        assert value["finite_endpoint_q_and_q_y_derivatives_retained"]
        assert value["no_differentiation_of_old_C0_caps"]
        for name in ("m13", "m15"):
            rec = value["normalized_first_Z_upper"][name]
            cap = current._raw_bound(owner.c, rec)
            assert cap.log is not None
            assert mp.isfinite(current.ep(cap.log)[1])
        inv = value["endpoint_derivative_bounds"]["inverse_Rc"]
        assert current.ep(current._raw_bound(owner.c, inv).log)[1] < 0


def _assert_raw_live(owner, raw):
    assert raw["source_family"] == owner.family
    assert raw["baseline"] == current._record(owner.baseline)
    assert raw["compact"] == current._record(owner.compact)
    assert raw["scale_ratios"] == current._record(owner.scale_ratios)
    assert raw["quiet_source"] == owner.quiet_source
    assert raw["absolute_Rc_source"]["native_frame_nodes"] == owner.radius_binding["expected_nodes"]
    assert raw["absolute_Rc_source"]["log_Rc"] == current._record(owner.log_Rc)
    for name, digest in raw["input_hashes"].items():
        assert current.sha(name) == digest, name


@current.source_precision
def run(owner, packets=None):
    """Check a live producer owner and write the local receipt."""
    began = mp.mpf("0")
    owner.assert_graph()
    raw = _read_raw()
    _assert_raw_live(owner, raw)
    _assert_scratch_and_baseline(owner)
    _assert_scale_and_compact(owner)
    fixture = _independent_compact_fixture()
    if packets is None:
        packets = owner.evaluate_all()
    assert set(packets) == {"|".join(cell) for cell in current.CELLS}
    for key, packet in packets.items():
        view = owner.report(packet)
        assert tuple(view["exact_Z_cell"]) == tuple(key.split("|"))
        assert view[current.GATE] is False
        assert view["complete_pressure_C1_not_installed"]
    before = mp.mp.dps
    with mp.workdps(15):
        packet = owner.evaluate(current.CELLS[0])
        owner.report(packet)
    assert mp.mp.dps == before
    inputs = dict(owner.hashes)
    inputs[current.NAME] = current.sha(current.NAME)
    inputs[Path(__file__).name] = current.sha(Path(__file__).name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            live_4x17_baseline_source_rows=True,
            scratch_layout_and_live_baseline_match=True,
            factorial_conversion_once=True,
            physical_Pstar_once_and_Z_independent=True,
            correlated_Rc_unit_logs=True,
            compact_differentiated_finite_N_bounds=True,
            finite_endpoint_q_and_q_y_derivatives_retained=True,
            quiet_Z1_exact_zero_record_retained=True,
            independent_nonunit_compact_identity=fixture,
            default_precision_restored=True,
            no_radius_or_N_materialization=True,
            no_full_pressure_C1_or_global_claim=True),
        exact_Z_cells=[list(cell) for cell in current.CELLS],
        input_hashes=inputs,
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        execution_seconds=0,
    )
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(result, separators=(",", ":")) + "\n",
        encoding="utf8", newline="\n")
    print("PASS_PRE_RP_C1_INPUTS", len(owner.baseline), len(owner.compact),
          flush=True)
    return result


if __name__ == "__main__":
    raise SystemExit("run(owner, packets) requires a retained owner")
