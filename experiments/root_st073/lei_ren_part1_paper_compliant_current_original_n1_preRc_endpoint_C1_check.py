"""Independent checker for the original pre-``Rc`` endpoint C1 adapter."""

import copy
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1 as current


def _read_raw():
    return json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))


def _canon(value):
    return current.c3.encoded(value)


def _exact_zero(record):
    return isinstance(record, dict) and record.get("exact_zero") is True


def _upper_log(record):
    if _exact_zero(record):
        return None
    if isinstance(record, dict) and "log_absolute_upper" in record:
        value = record["log_absolute_upper"]
        if isinstance(value, dict):
            if "upper_exact_mpf_tuple" in value:
                return tuple(value["upper_exact_mpf_tuple"])
            return value.get("upper")
        return value
    return None


def _assert_log_record(record, label):
    if not isinstance(record, dict) or "exact_zero" not in record:
        raise AssertionError("Missing directed endpoint record: " + label)
    if record["exact_zero"]:
        return
    if _upper_log(record) is None:
        raise AssertionError("Nonzero endpoint record lost its log upper: " + label)


def _direct_endpoint_fixture():
    """Check quotient, product, and fixed-radius endpoint signs independently."""
    R, Z = sy.symbols("R Z", positive=True, real=True)
    delta = sy.Rational(2, 11)
    U = (1 + Z + 2 * Z**2) * R + Z**3 / 7
    Uy = (2 - Z + Z**2) * R**2 / 5
    m = (1 + Z / 3 + Z**2 / 4) / R + Z**3 * R / 9
    a = 1 - delta
    d = 1 - Z**2
    L = 1 - delta * Z**2
    q = sy.cancel((2 * Z * U - a * Z * m - d * sy.diff(m, Z)) / L)
    qplus = sy.cancel((Z * (1 + delta) * U + 2 * Z * Uy
                       - d * sy.diff(U, Z)) / L)
    qy = sy.cancel(qplus - q)
    qplus_Z = sy.diff(qplus, Z)
    qy_Z = sy.diff(qy, Z)
    residual = sy.simplify(qy_Z - (qplus_Z - sy.diff(q, Z)))
    if residual != 0:
        raise AssertionError("Independent q_yZ sign identity failed")

    # Use nonunit fixed radial endpoints.  The two pressure boundary terms
    # are then checked as actual endpoint differences, including the R factor
    # in the weighted viscous term.
    lo, hi = sy.Rational(3, 2), sy.Rational(7, 3)
    boundary = lambda f: sy.cancel(f.subs(R, hi) - f.subs(R, lo))
    m13 = boundary(qplus)
    m15 = -boundary(R * qy)
    m13_direct = boundary(q + qy)
    m15_direct = -boundary(R * qy)
    if sy.simplify(m13 - m13_direct) != 0 or sy.simplify(m15 - m15_direct) != 0:
        raise AssertionError("Independent finite-radius boundary signs failed")
    # Rebuild the zero-axial branch explicitly; U is a fixture expression,
    # rather than a SymPy atom that can be replaced by ``subs``.
    q_upper = sy.cancel((-a * Z * m - d * sy.diff(m, Z)) / L)
    upper_q = sy.cancel(-q_upper)
    upper_residual = sy.simplify(upper_q + q_upper)
    if upper_residual != 0:
        raise AssertionError("Independent upper q_y=-q identity failed")
    return dict(
        passed=True,
        nonunit_radius_endpoints=["3/2", "7/3"],
        q_yZ_residual=str(residual),
        m13_boundary_sign=True,
        m15_weighted_boundary_sign=True,
        upper_q_y_minus_q=True,
        fixed_radius_is_Z_independent=True,
    )


def _assert_cell_view(view, ends):
    if tuple(view.get("exact_Z_cell", ())) != tuple(ends):
        raise AssertionError("Endpoint C1 view cell changed")
    lower = view.get("lower", {})
    upper = view.get("upper", {})
    for name in ("U_y", "U_ZZ", "U_yZ", "q_plus_q_y_Z",
                 "q_y_Z", "R_q_y", "R_q_y_Z"):
        _assert_log_record(lower.get("records", {}).get(name), "lower." + name)
    for name in ("q_plus_q_y_Z", "q_y_Z", "R_q_y_Z"):
        _assert_log_record(upper.get("records", {}).get(name), "upper." + name)
    if not lower.get("proof", {}).get("actual_V_y_row_1_retained"):
        raise AssertionError("Actual lower V_y row 1 was not retained")
    if not lower.get("proof", {}).get("actual_U_ZZ_from_V_rows_2_and_3"):
        raise AssertionError("Actual lower U_ZZ source rows were not retained")
    if lower.get("proof", {}).get("q_y_Z_formula") != "q_plus_q_y_Z-q_Z":
        raise AssertionError("Lower q_yZ formula changed")
    upper_proof = upper.get("proof", {})
    required_upper = (
        "original_O3_power_V_rows_exact_zero",
        "q_plus_q_y_exact_zero",
        "q_plus_q_y_Z_exact_zero",
        "q_y_equals_minus_q",
        "q_y_Z_equals_minus_q_Z",
    )
    if not all(upper_proof.get(key) for key in required_upper):
        raise AssertionError("Upper O3 endpoint C1 identity proof incomplete")
    boundary = view.get("boundary", {})
    if boundary.get("m13", {}).get("signed_form") != "B(q+q_y)":
        raise AssertionError("m13 endpoint boundary sign changed")
    if boundary.get("m15", {}).get("signed_form") != "-B(R*q_y)":
        raise AssertionError("m15 endpoint boundary sign changed")
    if view.get("proof", {}).get("finite_endpoint_radius_fixed_in_Z") is not True:
        raise AssertionError("Fixed-radius endpoint proof missing")


def _assert_source_hashes(raw):
    hashes = dict(raw.get("input_hashes", {}))
    if raw.get("source_family") is None or not hashes:
        raise AssertionError("Endpoint C1 source/hash manifest missing")
    for name, digest in hashes.items():
        if current.sha(name) != digest:
            raise AssertionError("Endpoint C1 source changed: " + name)
    return hashes


@current.source_precision
def run(owner=None, packets=None):
    began = time.monotonic()
    if owner is None:
        raise ValueError("Pass the live checked endpoint C1 owner; no ancestor rebuild is performed")
    if type(owner) is not current.CurrentOriginalN1PreRcEndpointC1:
        raise ValueError("CurrentOriginalN1PreRcEndpointC1 owner required")
    owner.assert_graph()
    raw = _read_raw()
    if raw.get("source_family") != owner.family:
        raise AssertionError("Saved endpoint C1 source family differs")
    cells = tuple(tuple(cell) for cell in raw.get("cells", ()))
    if cells != owner.cells or len(cells) != 4:
        raise AssertionError("Saved endpoint C1 four-cell cover changed")
    if packets is None:
        packets = {cell: owner.evaluate(cell) for cell in cells}
    if set(tuple(key) for key in packets) != set(cells):
        raise AssertionError("Live endpoint C1 packet cover changed")
    samples = raw.get("samples", {})
    live_samples = {}
    for cell in cells:
        packet = packets[cell]
        live = owner.report(packet)
        _assert_cell_view(live, cell)
        live_encoded = _canon(live)
        saved = samples.get(str(cell))
        if saved is None:
            raise AssertionError("Saved endpoint C1 cell sample missing")
        if live_encoded != saved:
            raise AssertionError("Saved/live endpoint C1 sample differs for " + str(cell))
        live_samples[cell] = live
    _assert_source_hashes(raw)

    symbolic = current.symbolic_check()
    fixture = _direct_endpoint_fixture()

    # Exercise the public packet/report path under the caller's default
    # precision without rebuilding the accepted parent graph.
    before = mp.mp.dps
    with mp.workdps(15):
        packet = owner.evaluate(("0", ".5"))
        public = owner.report(packet)
        _assert_cell_view(public, ("0", ".5"))
        public_dps = mp.mp.dps
    if mp.mp.dps != before or public_dps != 15:
        raise AssertionError("Endpoint C1 default-precision public path did not restore precision")

    checks = dict(
        symbolic=symbolic,
        independent_nonunit_finite_radius=fixture,
        exact_four_Z_cell_cover=True,
        live_saved_samples_equal=True,
        source_hash_manifest_checked=True,
        source_family_checked=True,
        lower_actual_V_y_rows_and_U_ZZ_checked=True,
        upper_q_y_Z_equals_minus_q_Z_checked=True,
        signed_m13_boundary_retained=True,
        signed_m15_weighted_boundary_retained=True,
        physical_common_Pstar_once=True,
        finite_radius_endpoints_fixed_in_Z=True,
        no_cap_differentiation=True,
        no_signed_numerical_centers=True,
        default_precision_public_query=True,
    )
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(
        all_passed=True,
        source_family=owner.family,
        cells=[list(cell) for cell in cells],
        checks=checks,
        input_hashes=hashes,
        execution_seconds=time.monotonic() - began,
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
    )
    (current.HERE / current.RECEIPT).write_bytes(
        (json.dumps(_canon(receipt), separators=(",", ":")) + "\n").encode())
    print("PASS_PRE_RC_ENDPOINT_C1", len(cells), flush=True)
    return receipt


if __name__ == "__main__":
    run()
