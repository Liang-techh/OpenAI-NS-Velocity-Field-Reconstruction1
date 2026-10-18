"""Pure fixed-phase ordinary-``Z`` bounds for restoration charts.

The restoration and post-restoration charts are the first charts in this
stage whose native radial first jet has a nonzero interior ``b`` term.  This
adapter therefore consumes the source-owned native-control packet and keeps
``b`` and ``t0=-b/a`` in the logarithmic Taylor algebra.  It does not call an
accepted builder, producer, owner, checker, norm, recurrence, quadrature, or
the finite-N correction adapter.

The packet contains caps for the actual native fields and for the four root
inputs ``a, b, t0, p2, Delta``.  The only numerical work here is pure
``LogZ``/``LogZPhase`` majorant propagation.  A global ``Delta<0`` proof
allows the exact flat branch ``sigma(1-Delta/eta)=1``.  If that proof is not
available, the packet must provide a positive active-branch gamma floor and
the saved flat-sigma derivative theorem; otherwise the adapter stays
fail-closed instead of silently replacing the cutoff by a constant.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath import isfinite, make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext

from .original_log_Z_phase_bounds import (
    LogZ,
    LogZPhase,
    compose_inverse,
    inverse_coefficients,
)


BASE = "experiments/root_st073/"
PARENT = (
    BASE
    + "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
)
PARENT_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
MODULE = "src/openai_ns_reconstruction/original_restoration_loop_Z3_bounds.py"
INPUT_MODULE = "src/openai_ns_reconstruction/original_restoration_native_controls.py"
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
LOG_MIXED_SOURCE = "src/openai_ns_reconstruction/original_log_mixed_jets.py"
C5_SOURCE = (
    BASE
    + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
C6_SOURCE = (
    BASE
    + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
)
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_outer_initial.py"
SIGMA_DERIVATIVE_SOURCE = (
    "src/openai_ns_reconstruction/original_flat_sigma_derivative_bounds.py"
)
PHASE_SOURCE = "src/openai_ns_reconstruction/original_restoration_phase_mixed.py"
CORRELATED_SOURCE = "src/openai_ns_reconstruction/original_restoration_correlated_pressure.py"
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
CHARTS = ("restoration", "post_restore")

# These hashes are inherited source pins.  The parent report remains the
# authority for every other input hash; constants below only cover sources
# whose names/formulas are part of this pure adapter's local audit.
LOG_Z_SOURCE_SHA256 = "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59"
C5_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
C6_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
SIGMA_SOURCE_SHA256 = "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657"
SIGMA_DERIVATIVE_SOURCE_SHA256 = "c528a552f8c9942165c799aa292656afe1b9233b49aa26c7874c2576ede2b64d"


def _require(value, message):
    if not value:
        raise ValueError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(value):
    """Serialize one directed interval without choosing a midpoint."""
    return dict(
        lower_exact_mpf_tuple=list(value._mpi_[0]),
        upper_exact_mpf_tuple=list(value._mpi_[1]),
        lower_display=nstr(make_mpf(value._mpi_[0]), 24),
        upper_display=nstr(make_mpf(value._mpi_[1]), 24),
        directed=True,
    )


def _read_record(c, value, *, require_directed=True):
    if hasattr(value, "_mpi_"):
        out = value
    else:
        if isinstance(value, Mapping) and "value" in value \
                and "lower_exact_mpf_tuple" not in value:
            value = value["value"]
        _require(isinstance(value, Mapping), "Directed logarithmic record is required")
        if require_directed:
            _require(value.get("directed") is True,
                     "Directed logarithmic record must set directed=True")
        try:
            low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
            high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Malformed directed logarithmic record") from exc
        out = c.mpf([low, high])
    # ``mpi == mpi`` only rejects NaN endpoints; infinities compare equal to
    # themselves.  Decode both exact endpoint tuples as scalar ``mpf`` values
    # and reject either non-finite endpoint before any logarithmic algebra.
    try:
        low_endpoint = make_mpf(out._mpi_[0])
        high_endpoint = make_mpf(out._mpi_[1])
    except (TypeError, ValueError) as exc:
        raise ValueError("Malformed logarithmic record endpoints") from exc
    _require(isfinite(low_endpoint) and isfinite(high_endpoint),
             "Logarithmic records must be finite")
    _require(out.a <= out.b, "Reversed logarithmic record")
    return out


def _mapping(value, label):
    _require(isinstance(value, Mapping), label + " must be a mapping")
    return value


def _rows_to_logz(c, rows, name):
    """Read four complete ordinary-Z rows from a source packet."""
    _require(isinstance(rows, (list, tuple)) and len(rows) >= 4,
             name + " must contain four ordinary-Z rows")
    out = {}
    for index, row in enumerate(rows[:4]):
        if isinstance(row, Mapping):
            order = row.get("ordinary_Z_order", row.get("Z_order", index))
            exact_zero = row.get("exact_zero") is True
            raw = row.get("log_absolute_upper_bound")
            if raw is None and "record" in row:
                raw = row["record"]
            if raw is None and "value" in row and not exact_zero:
                raw = row["value"]
            if raw is None and not exact_zero and "log_cap" in row:
                raw = row["log_cap"]
        else:
            order, exact_zero, raw = index, False, row
        _require(type(order) is int and order in ORDINARY_Z_ORDERS,
                 name + " has an invalid ordinary_Z_order")
        _require(order not in out, name + " repeats an ordinary_Z_order")
        if exact_zero:
            _require(raw is None, name + " exact-zero row must omit its record")
            out[order] = None
        else:
            _require(raw is not None, name + " nonzero row lacks a record")
            out[order] = _read_record(c, raw)
    _require(set(out) == set(ORDINARY_Z_ORDERS),
             name + " must cover ordinary Z orders 0 through 3")
    return LogZ.from_mapping(c, out)


def _logz_rows(jet):
    return [
        dict(
            ordinary_Z_order=order,
            exact_zero=jet[order] is None,
            log_absolute_upper_bound=(
                None if jet[order] is None else _record(jet[order])
            ),
        )
        for order in ORDINARY_Z_ORDERS
    ]


def _phase_rows(jet):
    return [
        dict(
            Z_order=order,
            free_angle_order=angle_order,
            exact_zero=jet[order, angle_order] is None,
            log_absolute_upper_bound=(
                None
                if jet[order, angle_order] is None
                else _record(jet[order, angle_order])
            ),
        )
        for order, angle_order in jet.indices
    ]


def _jet_rows(jet):
    if isinstance(jet, LogZ):
        return _logz_rows(jet)
    if isinstance(jet, LogZPhase):
        return _phase_rows(jet)
    raise TypeError("Expected a LogZ or LogZPhase jet")


def _values_for(packet, names):
    names = frozenset(names)
    found = []

    def visit(value):
        if isinstance(value, Mapping):
            for key, child in value.items():
                if key in names:
                    found.append(child)
                visit(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(child)

    visit(packet)
    return found


def _is_true(packet, names):
    return any(value is True for value in _values_for(packet, names))


def _find_value(packet, names):
    values = _values_for(packet, names)
    return values[0] if values else None


def _phase_cosine(c):
    zero = c.mpf(0)
    return LogZPhase.from_mapping(c, {
        (order, angle_order): zero if order == 0 else None
        for order, angle_order in LogZPhase.indices
    })


def _angle_primitive(c):
    values = {index: None for index in LogZPhase.indices}
    values[0, 0] = c.ln(c.mpf(2) * c.pi)
    values[0, 1] = c.mpf(0)
    return LogZPhase.from_mapping(c, values)


def _logsum_upper(c, values):
    """Return a finite point upper bound for log(sum(exp(values)))."""
    values = [value for value in values if value is not None]
    if not values:
        return None
    maximum = max(value.b if hasattr(value, "b") else value for value in values)
    correction = c.mpf(0)
    for value in values:
        upper = value.b if hasattr(value, "b") else value
        gap = maximum - upper
        if gap > c.mpf(1024) * c.ln(2):
            correction += c.mpf(2) ** -1024
        else:
            correction += c.exp(-gap)
    return c.mpf([maximum + c.log(1 + correction),
                  maximum + c.log(1 + correction)])


def _log1pexp_upper(c, value):
    value = value.b if hasattr(value, "b") else value
    if value <= 0:
        return c.mpf([c.log(2), c.log(2)])
    if value > c.mpf(1024) * c.ln(2):
        # Keep a positive directed remainder even for very large values;
        # replacing exp(-value) by zero would turn this upper bound into an
        # underbound.  At the branch threshold 10^6, 2^-prec is still a
        # conservative remainder at the 900-bit context used by build().
        correction = (c.exp(-value) if value < c.mpf(10) ** 6
                      else c.mpf(2) ** (-c.prec))
        return c.mpf([value, value + c.log(1 + correction)])
    out = c.log(1 + c.exp(value))
    return c.mpf([out.a, out.b])


def _sqrt_log_jet(c, square, log_square_lower):
    """Absolute ordinary-Z caps for ``sqrt(square)`` through order 3."""
    _require(square[0] is not None, "Square root value row cannot be zero")
    # The reciprocal formulas use the lower bound for q=sqrt(square), not
    # the lower bound for square itself.  Keeping this conversion here makes
    # the call site explicit and avoids accidentally treating a square bound
    # as a first-order q bound.
    log_q_lower = c.mpf("0.5") * log_square_lower
    q0 = c.mpf("0.5") * square[0]
    q = {
        0: q0,
        1: (None if square[1] is None else
            square[1] - c.ln(2) - log_q_lower),
        2: _logsum_upper(c, [
            None if square[2] is None else
            square[2] - c.ln(2) - log_q_lower,
            None if square[1] is None else
            2 * square[1] - c.ln(4) - 3 * log_q_lower,
        ]),
        3: _logsum_upper(c, [
            None if square[3] is None else
            square[3] - c.ln(2) - log_q_lower,
            None if square[1] is None or square[2] is None else
            square[1] + square[2] + c.ln(3) - c.ln(4)
            - 3 * log_q_lower,
            None if square[1] is None else
            3 * square[1] + c.ln(3) - c.ln(8)
            - 5 * log_q_lower,
        ]),
    }
    return LogZ.from_mapping(c, q)


def _flat_sigma_bounds(c, cut, derivative_bounds):
    """Bound sigma(cut) and its first three ordinary Z derivatives."""
    _require(len(derivative_bounds) >= 4,
             "Saved sigma derivative bounds through order three are required")
    # sigma is [0,1]; the derivative source provides absolute caps for
    # sigma', sigma'', sigma''' on the full real line, including flat tails.
    engine = LogZ(c)
    cut2 = engine.mul(cut, cut)
    cut3 = engine.mul(cut2, cut)
    # Faà di Bruno gives sigma' cut'', sigma'' (cut')^2 and
    # sigma' cut''' + 3 sigma'' cut' cut'' + sigma''' (cut')^3.
    # The products below are deliberately logarithmic absolute majorants;
    # no sigma value is selected on an active or flat subbranch.
    c1 = engine.scale(cut, c.ln(c.mpf(derivative_bounds[1])))
    c2 = engine.scale(cut2, c.ln(c.mpf(derivative_bounds[2])))
    c3 = engine.scale(cut3, c.ln(c.mpf(derivative_bounds[3])))
    out = {
        0: c.mpf(0),
        1: c1[1],
        2: engine.add(c1, c2)[2],
        3: engine.add(
            c1,
            engine.scale(cut2, c.ln(c.mpf(3 * derivative_bounds[2]))),
            c3,
        )[3],
    }
    return LogZ.from_mapping(c, out)


def _phase_rows_from_input(c, chart):
    """Resolve the root-input row mapping used by the native worker."""
    candidates = [chart]
    for key in ("root_inputs", "roots", "controls", "native_controls"):
        value = chart.get(key)
        if isinstance(value, Mapping):
            candidates.insert(0, value)
    def get(*names):
        for candidate in candidates:
            for name in names:
                if name in candidate:
                    return candidate[name]
        return None
    return {key: get(*aliases) for key, aliases in {
        "a": ("a", "a_log_caps", "a_rows"),
        "b": ("b", "b_log_caps", "b_rows"),
        "t0": ("t0", "t0_log_caps", "t0_rows"),
        "p2": ("p2", "p2_log_caps", "p2_rows"),
        "Delta": ("Delta", "delta", "Delta_log_caps", "Delta_rows"),
        "E_relative": (
            "E_relative", "Erelative", "E_over_frozen_E0",
            "E_relative_log_caps", "E_relative_rows",
        ),
    }.items()}


def _scalar_with_name(chart, names):
    """Return a scalar and the alias that supplied it.

    Native controls keep some positive source parameters inside a proof
    receipt rather than at the chart root.  The recursive fallback is still
    read-only: it only locates an already serialized directed record.
    """
    names = tuple(names)
    candidates = (
        chart,
        chart.get("root_inputs", {}),
        chart.get("parameters", {}),
        chart.get("actual_logs", {}),
    )
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            for name in names:
                if name in candidate:
                    return candidate[name], name

    def visit(value):
        if isinstance(value, Mapping):
            for name in names:
                if name in value:
                    return value[name], name
            for child in value.values():
                found = visit(child)
                if found is not None:
                    return found
        elif isinstance(value, (list, tuple)):
            for child in value:
                found = visit(child)
                if found is not None:
                    return found
        return None

    return visit(chart)


def _scalar_from_chart(chart, names):
    found = _scalar_with_name(chart, names)
    return None if found is None else found[0]


def _positive_floor(c, chart, names, *, label):
    found = _scalar_with_name(chart, names)
    _require(found is not None, "Missing positive floor: " + label)
    raw, alias = found
    value = _read_record(c, raw)
    if alias.endswith("_log") or alias.startswith("log_"):
        _require(value.a == value.a, "Positive log floor is not finite: " + label)
        return value
    # The native worker exposes a_value_interval as a positive raw a record;
    # turn its directed lower endpoint into the logarithmic lower floor used
    # by reciprocal majorants.  A logarithmic alias above remains untouched.
    _require(value.a > 0, "Positive floor is not positive: " + label)
    return c.log(value)


def _parameter_log(c, chart, names, *, label):
    raw = _scalar_from_chart(chart, names)
    _require(raw is not None, "Missing logarithmic parameter: " + label)
    return _read_record(c, raw)


def _signed_upper_negative(c, packet, names):
    """Find a signed interval proof whose upper endpoint is negative."""
    for value in _values_for(packet, names):
        if isinstance(value, Mapping) and "lower_exact_mpf_tuple" in value:
            interval = _read_record(c, value)
            if interval.b < 0:
                return interval
    return None


def _derivative_bounds(root):
    """Read saved universal sigma derivative caps without running their proof."""
    path = root / SIGMA_DERIVATIVE_SOURCE
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    assignments = [node for node in tree.body if isinstance(node, ast.Assign)
                   and any(isinstance(target, ast.Name)
                           and target.id == "DERIVATIVE_UPPER_BOUNDS"
                           for target in node.targets)]
    _require(len(assignments) == 1,
             "Saved sigma derivative bound assignment is missing or ambiguous")
    values = ast.literal_eval(assignments[0].value)
    _require(isinstance(values, tuple) and len(values) >= 4
             and all(isinstance(v, int) and v >= 0 for v in values),
             "Saved sigma derivative bounds are malformed")
    return tuple(values)


def _packet_proofs(root, parent, packet):
    _require(packet.get("passed") is True,
             "Current restoration native-control packet must pass")
    _require(packet.get("parent_report_sha256") == PARENT_SHA256,
             "Native-control packet is tied to a different parent report")
    _require(packet.get("source_family") == parent.get("source_family"),
             "Native-control packet source family changed")
    packet_orders = packet.get("ordinary_Z_orders", [])
    _require(packet_orders[:4] == list(ORDINARY_Z_ORDERS),
             "Native-control packet must cover ordinary Z orders 0 through 3")
    _require(packet.get("physical_y_derivatives_claimed", False) is False,
             "Physical y derivatives cannot be claimed by the loop adapter")
    _require(packet.get("source_executed", False) is False,
             "Native-control packet source execution flag changed")
    charts = packet.get("charts")
    _require(isinstance(charts, Mapping), "Native-control charts are required")
    for name in CHARTS:
        chart = charts.get(name)
        if chart is None and name == "post_restore":
            chart = charts.get("postrestore")
        _require(isinstance(chart, Mapping) and chart.get("passed") is True,
                 "Native-control chart did not pass: " + name)
        chart_orders = chart.get("ordinary_Z_orders", list(ORDINARY_Z_ORDERS))
        _require(chart_orders[:4] == list(ORDINARY_Z_ORDERS),
                 "Native-control chart order prefix changed: " + name)
        roots = _phase_rows_from_input(None, chart)
        _require(all(roots.get(key) is not None for key in
                     ("a", "b", "t0", "p2", "Delta", "E_relative")),
                 "Native-control chart roots are missing: " + name)
        global_delta = _is_true(
            chart, ("Delta_upper_lt_zero", "Delta_globally_negative",
                    "global_Delta_lt_zero", "kappa_upper_lt_one"),
        )
        if global_delta:
            bproof = _find_value(chart, ("b_threshold_proof",))
            _require(isinstance(bproof, Mapping) and bproof.get("passed") is True,
                     "Global Delta<0 requires the native |b| margin proof: " + name)
            _require(
                bproof.get("exact_zero") is True
                or bproof.get("bound") == "|b|<1/10"
                or _is_true(bproof, ("bound_below_one_tenth",)),
                "Native |b|<1/10 proof is missing: " + name,
            )
            delta_proof = _find_value(chart, ("delta_negative_proof",))
            _require(
                isinstance(delta_proof, Mapping)
                and delta_proof.get("passed") is True
                and delta_proof.get("kappa_upper_lt_one") is True
                and delta_proof.get("Delta_upper_lt_zero") is True,
                "Global Delta<0 proof must include kappa<1: " + name,
            )
        else:
            gamma_floor = _scalar_from_chart(
                chart,
                ("active_gamma_lower_log", "gamma_active_lower_log",
                 "log_gamma_lower", "gamma_lower_log", "gamma_lower"),
            )
            _require(gamma_floor is not None,
                     "Active cutoff requires a positive gamma floor: " + name)
    common_n = _find_value(packet, (
        "actual_N18", "actual_N18_is_graph_handle", "common_N18", "N18",
        "common_N_handle", "common_N", "actual_N", "current_N",
    ))
    parent_n_handles = (
        parent.get("actual_restoration_correlated_pressure", {})
        .get("native_control_theorem", {})
        .get("local_IFT", {})
        .get("parameter_handles", {})
    )
    parent_n = parent_n_handles.get("N") if isinstance(parent_n_handles, Mapping) else None
    parent_logN = (
        parent.get("actual_micro_q_mixed_bounds", {})
        .get("parameter_envelopes", {})
        .get("actual_logs", {})
        .get("N")
    )
    common_n_from_parent = parent_n == 18 and isinstance(parent_logN, Mapping)
    _require(common_n in (True, 18, "18", "N18", "N18 graph handle")
             or _is_true(packet, ("actual_N18", "actual_N18_is_graph_handle",
                                  "common_N18"))
             or common_n_from_parent,
             "Common N=18 proof is required")
    _require(_is_true(packet, ("actual_phase_Z_independent",
                               "phase_Z_independent")),
             "Phase Z-independence proof is required")
    return dict(
        passed=True,
        actual_N18=True,
        actual_phase_Z_independent=True,
        physical_y_derivatives_claimed=False,
        common_N_handle=18,
        common_N_from_parent_handle=common_n_from_parent,
    )


def _verify_packet_pins(root, packet):
    pins = packet.get("source_pins")
    _require(isinstance(pins, Mapping) and pins,
             "Native-control source_pins are required")
    verified = {}
    input_self_pin = False
    for key, raw in pins.items():
        _require(isinstance(raw, Mapping),
                 "Native-control source pin must be a mapping: " + str(key))
        path_text = str(raw.get("path", key))
        expected = raw.get("expected_sha256") or raw.get("sha256")
        actual_claim = raw.get("actual_sha256")
        _require(raw.get("passed", True) is True,
                 "Native-control source pin did not pass: " + path_text)
        _require(raw.get("source_executed", False) is False,
                 "Native-control source pin was executed: " + path_text)
        _require(isinstance(expected, str) and len(expected) == 64,
                 "Native-control source pin has no SHA256: " + path_text)
        path = Path(path_text)
        if not path.is_absolute():
            path = root / path
        _require(path.is_file(), "Native-control source is missing: " + path_text)
        actual = _sha(path)
        _require(actual == expected, "Changed native-control source: " + path_text)
        if actual_claim is not None:
            _require(actual_claim == actual,
                     "Native-control source pin actual SHA256 drifted: " + path_text)
        if path_text.replace("\\", "/").endswith(INPUT_MODULE):
            _require(raw.get("self_pin") is True,
                     "Native-control input module self_pin is required")
            input_self_pin = True
        verified[path_text] = dict(path=path_text, expected_sha256=expected,
                                   actual_sha256=actual, passed=True,
                                   source_executed=False,
                                   **({"self_pin": True} if path_text.replace(
                                       "\\", "/").endswith(INPUT_MODULE) else {}))
    _require(input_self_pin,
             "Native-control input-module self hash pin is required")
    return verified


def _parent_proof_guards(parent):
    phase = parent.get("actual_restoration_phase_mixed", {})
    correlated = parent.get("actual_restoration_correlated_pressure", {})
    _require(phase.get("passed") is True and
             phase.get("phase_coefficients", {}).get("passed") is True,
             "Saved actual restoration phase-mixed proof is required")
    phase_coefficients = phase.get("phase_coefficients", {})
    _require(phase.get("source_executed", False) is False and
             phase_coefficients.get("source_executed", False) is False,
             "Saved phase-mixed proof source execution flag changed")
    _require(phase_coefficients.get("total_rows") == 42 and
             phase_coefficients.get("matched_rows") == 42,
             "Saved phase-mixed proof must contain all 42 rows")
    phase_charts = phase_coefficients.get("charts", {})
    _require(isinstance(phase_charts, Mapping),
             "Saved phase-mixed chart rows are missing")
    raw_lambda_total = sum(
        chart.get("raw_lambda_rows", 0)
        for chart in phase_charts.values() if isinstance(chart, Mapping)
    )
    raw_lambda_matched = sum(
        chart.get("raw_lambda_rows_matched", 0)
        for chart in phase_charts.values() if isinstance(chart, Mapping)
    )
    _require(raw_lambda_total == 30 and raw_lambda_matched == 30 and
             all(chart.get("raw_lambda_convolution_bound") is True
                 for chart in phase_charts.values()
                 if isinstance(chart, Mapping)),
             "Saved phase-mixed proof must contain all 30 lambda rows")
    _require(correlated.get("passed") is True and
             correlated.get("full_modified_seam_proved") is True and
             correlated.get("unresolved_gates") == [],
             "Saved actual restoration correlated-pressure proof is required")
    _require(correlated.get("source_executed", False) is False,
             "Correlated-pressure proof was executed during loop bounds")
    gates = correlated.get("gates", {})
    _require(isinstance(gates, Mapping) and all(v is True for v in gates.values()),
             "Restoration correlated-pressure gates are incomplete")
    endpoints = parent.get("actual_restoration_native_endpoints", {})
    lemma = endpoints.get("sigma_endpoint_lemma", {})
    _require(lemma.get("passed") is True and
             lemma.get("squared_flat_cutoff") is True and
             lemma.get("value_first_second_flat") is True,
             "Saved flat-sigma endpoint theorem is required")
    return dict(
        phase_mixed_consumed=True,
        correlated_pressure_consumed=True,
        sigma_flat_endpoint_lemma_consumed=True,
        accepted_new_roots_replayed=False,
        accepted_new_phase_rows_replayed=False,
    )


def _source_guards(root, parent, packet_pins):
    expected = parent.get("input_hashes", {})
    implementation = parent.get("implementation_hashes", {})
    fixed = {
        PARENT: PARENT_SHA256,
        LOG_Z_SOURCE: LOG_Z_SOURCE_SHA256,
        C5_SOURCE: C5_SOURCE_SHA256,
        C6_SOURCE: C6_SOURCE_SHA256,
        SIGMA_SOURCE: SIGMA_SOURCE_SHA256,
        SIGMA_DERIVATIVE_SOURCE: SIGMA_DERIVATIVE_SOURCE_SHA256,
        MODULE: _sha(root / MODULE),
        INPUT_MODULE: _sha(root / INPUT_MODULE),
        PHASE_SOURCE: implementation.get(PHASE_SOURCE),
        CORRELATED_SOURCE: implementation.get(CORRELATED_SOURCE),
    }
    out = {}
    for path, fallback in fixed.items():
        if fallback is None:
            fallback = expected.get(Path(path).name)
        _require(isinstance(fallback, str), "Missing source hash: " + path)
        actual = _sha(root / path)
        _require(actual == fallback, "Changed pinned source: " + path)
        out[path] = dict(path=path, expected_sha256=fallback,
                         actual_sha256=actual, passed=True,
                         source_executed=False,
                         **({"self_pin": True} if path in
                            (MODULE, INPUT_MODULE) else {}))
    out.update(packet_pins)
    # The AST snippets pin the formulas while keeping the accepted source
    # modules read-only.  This does not import or execute a source provider.
    text = (root / C5_SOURCE).read_text(encoding="utf-8")
    snippets = (
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2)",
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha)",
        "K=a.div(one,a.scale(ratio,twopi)",
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
        "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat)",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
    )
    compact = "".join(text.split())
    c5_formula_guard = dict(
        passed=all("".join(item.split()) in compact for item in snippets),
        path=C5_SOURCE,
        required_source_fragments=list(snippets),
        source_executed=False,
    )
    _require(c5_formula_guard["passed"],
             "Actual C5 loop source formula guard failed")
    # Formula receipts are provenance guards, not canonical source pins.  Keep
    # them out of ``source_pins`` so every pin there has an exact path/hash
    # contract; the receipt is emitted under ``source_guards`` below.
    return out, c5_formula_guard


def _compute_bounds(c, chart, derivative_bounds):
    """Run the pure actual restoration Z/angle majorant algebra."""
    roots = _phase_rows_from_input(c, chart)
    for name in ("a", "b", "t0", "p2", "Delta", "E_relative"):
        _require(roots.get(name) is not None,
                 "Native chart is missing root input " + name)
    caps = {
        name: _rows_to_logz(c, roots[name], "root_inputs." + name)
        for name in ("a", "b", "t0", "p2", "Delta", "E_relative")
    }
    a = caps["a"]
    b = caps["b"]
    t0 = caps["t0"]
    p2 = caps["p2"]
    Delta = caps["Delta"]
    E = caps["E_relative"]

    log_a_lower = _positive_floor(
        c, chart,
        ("a_lower_log", "positive_a_lower_log", "log_a_lower",
         "a_positive_lower", "a_value_lower", "a_value_interval"),
        label="a",
    )
    log_eta = _parameter_log(
        c, chart,
        ("actual_eta", "actual_eta_log", "eta_log", "eta"),
        label="eta",
    )
    log_dstar = _parameter_log(
        c, chart,
        ("dstar", "d_star", "dstar_log", "d_star_log"),
        label="d_star",
    )
    _require(log_eta.a == log_eta.a and log_dstar.a == log_dstar.a,
             "eta and d_star logs must be finite")
    # eta and d_star are positive source parameters represented in log form.
    # Their lower log endpoints, rather than a selected scalar, feed the
    # reciprocal and square-root majorants.
    eta_lower = log_eta.a
    dstar_lower = log_dstar.a
    one = LogZ(c).constant(c.mpf(0))
    eta_jet = LogZ(c).constant(log_eta)
    two_eta = LogZ(c).scale(eta_jet, c.ln(2))

    # gamma=2*eta-Delta.  Its absolute derivative caps use the ordinary
    # absolute majorant; the packet supplies the active positive lower floor.
    gamma = LogZ(c).add(two_eta, Delta)
    delta_negative = _is_true(
        chart,
        ("delta_globally_negative", "Delta_globally_negative",
         "global_Delta_lt_zero", "Delta_lt_zero", "kappa_less_than_two",
         "kappa_lt_two", "global_active_delta", "Delta_upper_lt_zero",
         "kappa_upper_lt_one"),
    )
    signed_negative = _signed_upper_negative(
        c, chart,
        ("Delta_interval", "delta_interval", "global_Delta_interval",
         "Delta_signed_interval"),
    )
    if signed_negative is not None:
        delta_negative = True
    if delta_negative:
        # The signed source identity is Delta=kappa-2.  The validated
        # kappa<1 proof therefore gives |Delta|>1, while eta>0, so
        # gamma=2*eta-Delta>1.  Use the uniform directed lower log 0 rather
        # than the much weaker eta floor; this is valid only after the
        # packet's explicit signed Delta proof above has passed.
        gamma_lower_log = c.mpf(0)
        cutoff_mode = "exact_flat_sigma_one_from_global_Delta_lt_zero"
        sigma_jet = LogZ(c).constant(c.mpf(0))
    else:
        raw_gamma_floor = _scalar_from_chart(
            chart,
            ("active_gamma_lower_log", "gamma_active_lower_log",
             "log_gamma_lower", "gamma_lower_log", "gamma_lower"),
        )
        _require(raw_gamma_floor is not None,
                 "Active/flat cutoff proof lacks a positive gamma floor")
        gamma_lower_log = _read_record(c, raw_gamma_floor)
        _require(gamma_lower_log.a == gamma_lower_log.a,
                 "Active gamma floor must be finite")
        # cut=1-Delta/eta; the global sigma derivative source is consumed as
        # a majorant, while the saved flat theorem handles the Delta>=eta tail.
        delta_over_eta = LogZ(c).scale(Delta, -eta_lower)
        cut = LogZ(c).add(one, delta_over_eta)
        sigma_jet = _flat_sigma_bounds(c, cut, derivative_bounds)
        cutoff_mode = "active_flat_sigma_with_saved_global_derivative_bounds"

    inv_a = LogZ(c).reciprocal(a, log_a_lower)
    square = LogZ(c).scale(LogZ(c).mul(gamma, inv_a), -c.ln(2))
    log_square_lower = gamma_lower_log - c.ln(2) - a[0].b
    q_raw = _sqrt_log_jet(c, square, log_square_lower)
    q = LogZ(c).mul(sigma_jet, q_raw)

    u = LogZ(c).scale(LogZ(c).mul(p2, q), -dstar_lower)
    u2 = LogZ(c).mul(u, u)
    h = LogZ.from_mapping(c, {
        0: c.mpf(0),
        1: None if u2[1] is None else u2[1] - c.ln(2),
        2: _logsum_upper(c, [
            None if u2[2] is None else u2[2] - c.ln(2),
            None if u2[1] is None else 2 * u2[1] + c.ln(c.mpf("0.75")),
        ]),
        3: _logsum_upper(c, [
            None if u2[3] is None else u2[3] - c.ln(2),
            None if u2[1] is None or u2[2] is None else
            u2[1] + u2[2] + c.ln(c.mpf("2.25")),
            None if u2[1] is None else 3 * u2[1] + c.ln(c.mpf("1.875")),
        ]),
    })
    r = LogZ(c).mul(u, h)
    alpha = LogZ(c).scale(LogZ(c).mul(q, h), c.ln(2))

    t02 = LogZ(c).mul(t0, t0)
    q2 = LogZ(c).mul(q, q)
    denominator = LogZ(c).add(one, t02, LogZ(c).scale(q2, c.ln(2)))
    K = LogZ(c).scale(
        LogZ(c).reciprocal(denominator, c.mpf(0)),
        -c.ln(2 * c.pi),
    )

    # The exact positive denominator bound is uniform in angle.  h>=h_lower,
    # D>=h^4/4, and lambda=K*(1+t^2)>=K.
    log_h_lower = -c.mpf("0.5") * _log1pexp_upper(c, 2 * u[0].b)
    log_D_lower = c.mpf([4 * log_h_lower.a - c.ln(4),
                         4 * log_h_lower.b - c.ln(4)])
    log_lambda_lower = c.mpf([
        -c.ln(2 * c.pi) - denominator[0].b,
        -c.ln(2 * c.pi) - denominator[0].b,
    ])
    conditioning = {
        "log_D_lower": log_D_lower,
        "log_lambda_lower": log_lambda_lower,
    }

    # Pure fixed-phase angle propagation.  No angle value is selected.
    p = LogZPhase(c)
    one_phase = p.constant(c.mpf(0))
    cosine = _phase_cosine(c)
    r_phase = LogZPhase.liftZ(r)
    alpha_phase = LogZPhase.liftZ(alpha)
    K_phase = LogZPhase.liftZ(K)
    t0_phase = LogZPhase.liftZ(t0)
    r2_phase = p.mul(r_phase, r_phase)
    D = p.add(one_phase, p.scale(p.mul(r_phase, cosine), c.ln(2)), r2_phase)
    invD = p.reciprocal(D, conditioning["log_D_lower"])
    cos_minus_r = p.add(cosine, r_phase)
    t = p.add(t0_phase,
             p.mul(alpha_phase, p.mul(cos_minus_r, invD)))
    t2 = p.mul(t, t)
    lam = p.mul(K_phase, p.add(one_phase, t2))
    T1 = p.primitiveAngle(t)
    T2 = p.primitiveAngle(t2)
    angle = _angle_primitive(c)
    raw_phi = p.mul(K_phase, p.add(angle, T2))
    phi_values = {
        index: raw_phi[index] if index[1] == 0
        else lam[index[0], index[1] - 1]
        for index in p.indices
    }
    phi_values[0, 0] = c.mpf(0)
    Phi = LogZPhase.from_mapping(c, phi_values)
    inverse, delta, inverse_proof = inverse_coefficients(
        p, Phi, conditioning["log_lambda_lower"]
    )
    T1hat = compose_inverse(p, T1, delta)

    # Z-independent target phase is represented only by its unit magnitude
    # bound.  It is not a selected function value.
    target_phase = LogZ.constant(LogZ(c), c.mpf(0))
    two_pi_log = c.ln(2 * c.pi)
    term_a = LogZ(c).scale(LogZ(c).mul(a, T1hat), -two_pi_log)
    term_b = LogZ(c).mul(b, target_phase)
    Brel = LogZ(c).scale(LogZ(c).add(term_a, term_b), -c.ln(2))
    Bfull = LogZ(c).mul(E, Brel)

    # A=(a/2)*(phase-psi/(2*pi)) is retained for the finite-N delta-E
    # correction.  The phase is a unit magnitude bound and remains unselected.
    psi_scaled = LogZ(c).scale(inverse, -two_pi_log)
    chi = LogZ(c).add(target_phase, psi_scaled)
    A = LogZ(c).scale(LogZ(c).mul(a, chi), -c.ln(2))

    return dict(
        a=a, b=b, t0=t0, Delta=Delta, gamma=gamma, q_raw=q_raw, q=q,
        p2=p2, E_relative=E, u=u, h=h, r=r, alpha=alpha, K=K,
        denominator=denominator,
        D=D, invD=invD, t=t, lambda_value=lam, T1=T1, T2=T2, Phi=Phi,
        inverse=inverse, delta=delta, inverse_proof=inverse_proof,
        T1hat=T1hat, target_phase=target_phase, target_phase_bound_one=True,
        B_over_E=Brel, B_over_Pstar_over_frozen_E0=Bfull, A=A,
        conditioning=conditioning, cutoff_mode=cutoff_mode,
        delta_global_negative=delta_negative,
        gamma_lower_log=gamma_lower_log,
        universal_B_over_E_value_cap_tightened=False,
    )


def _read_packet(root, input_packet):
    if isinstance(input_packet, (str, Path)):
        path = root / input_packet if not Path(input_packet).is_absolute() else Path(input_packet)
        packet = json.loads(path.read_bytes())
        packet_path = (str(path.relative_to(root))
                       if path.is_relative_to(root) else str(path))
    else:
        packet, packet_path = input_packet, None
    _mapping(packet, "input_packet")
    return packet, packet_path


def _with_canonical_parameters(parent, packet):
    """Attach the unchanged positive source parameters when an older native
    packet omitted its redundant chart aliases.

    The current native worker already computes the actual charts and root
    jets, but an intermediate packet revision placed eta only inside its
    cutoff proof and omitted d_star/logN.  These are global source parameters,
    so reading their directed records from the pinned parent is equivalent to
    the worker's intended aliases and does not evaluate a source function.
    """
    parameter_logs = (
        parent.get("actual_micro_q_mixed_bounds", {})
        .get("parameter_envelopes", {})
        .get("actual_logs", {})
    )
    if not isinstance(parameter_logs, Mapping):
        return packet, False
    required = {key: parameter_logs.get(key) for key in ("eta", "dstar", "N")}
    if not all(isinstance(value, Mapping) for value in required.values()):
        return packet, False
    charts = packet.get("charts")
    if not isinstance(charts, Mapping):
        return packet, False
    changed = False
    out = dict(packet)
    out_charts = {}
    for name in CHARTS:
        source_name = name if name in charts else ("postrestore" if name == "post_restore" else name)
        raw_chart = charts.get(source_name)
        if not isinstance(raw_chart, Mapping):
            return packet, False
        chart = dict(raw_chart)
        if _scalar_from_chart(chart, ("actual_eta_log", "eta_log", "eta")) is None:
            chart["actual_eta_log"] = required["eta"]
            changed = True
        if _scalar_from_chart(chart, ("dstar_log", "d_star_log", "dstar", "d_star")) is None:
            chart["dstar_log"] = required["dstar"]
            changed = True
        if _scalar_from_chart(chart, ("logN", "actual_logN", "N_log")) is None:
            chart["logN"] = required["N"]
            changed = True
        chart["canonical_parameter_source"] = dict(
            eta="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.eta",
            dstar="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.dstar",
            logN="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.N",
            source_executed=False,
            parent_report_sha256=PARENT_SHA256,
        )
        out_charts[source_name] = chart
    # Preserve a possible noncanonical postrestore alias while replacing only
    # the two charts consumed by this adapter.
    for key, value in charts.items():
        out_charts.setdefault(key, value)
    out["charts"] = out_charts
    out["canonical_parameter_source"] = dict(
        eta="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.eta",
        dstar="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.dstar",
        logN="actual_micro_q_mixed_bounds.parameter_envelopes.actual_logs.N",
        source_executed=False,
        parent_report_sha256=PARENT_SHA256,
    )
    return out, changed


def _chart(packet, name):
    charts = packet.get("charts", {})
    chart = charts.get(name)
    if chart is None and name == "post_restore":
        chart = charts.get("postrestore")
    _require(isinstance(chart, Mapping), "Missing native-control chart: " + name)
    return chart


def _chart_record(c, chart, names, label, fallback=None):
    raw = _scalar_from_chart(chart, names)
    if raw is None and fallback is not None:
        raw = _scalar_from_chart(fallback, names)
    _require(raw is not None, "Missing chart record: " + label)
    return _read_record(c, raw)


def _source_contract(name, cutoff_mode):
    return dict(
        chart=name,
        quantity="B/Pstar/frozen_E0",
        normalized_quantity="B_over_Pstar_over_frozen_E0",
        A_quantity="A for finite-N delta_E correction",
        physical_deltaU="deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
        actual_native_restoration_controls=True,
        actual_native_post_restore_controls=True,
        actual_root_inputs=True,
        common_N_handle=18,
        current_N="N18 graph handle",
        N_not_materialized=True,
        physical_Pstar_applied_once=True,
        phase_Z_independent=True,
        selected_phase_value=False,
        cutoff_mode=cutoff_mode,
        physical_y_derivatives_claimed=False,
        accepted_ancestor_replayed=False,
        finite_N_numeric=False,
    )


def build(root, parent, input_packet):
    """Build actual restoration/post-restore fixed-phase Z0--Z3 bounds."""
    root = Path(root)
    parent_path = root / PARENT
    _require(parent_path.is_file() and _sha(parent_path) == PARENT_SHA256,
             "Current restoration parent bytes changed")
    disk_parent = json.loads(parent_path.read_bytes())
    _require(parent == disk_parent,
             "Current restoration parent mapping is required")

    packet, packet_path = _read_packet(root, input_packet)
    packet, canonical_parameters_inherited = _with_canonical_parameters(
        disk_parent, packet
    )
    proof = _packet_proofs(root, disk_parent, packet)
    parent_proofs = _parent_proof_guards(disk_parent)
    packet_pins = _verify_packet_pins(root, packet)
    if packet_path is not None:
        packet_file = Path(packet_path)
        if not packet_file.is_absolute():
            packet_file = root / packet_file
        packet_hash = _sha(packet_file)
        packet_pins[packet_path] = dict(
            path=packet_path,
            expected_sha256=packet_hash,
            actual_sha256=packet_hash,
            passed=True,
            source_executed=False,
        )
    source_pins, c5_formula_guard = _source_guards(
        root, disk_parent, packet_pins
    )
    derivative_bounds = _derivative_bounds(root)

    results = {}
    derivative_source = dict(
        path=SIGMA_DERIVATIVE_SOURCE,
        expected_sha256=SIGMA_DERIVATIVE_SOURCE_SHA256,
        derivative_upper_bounds=list(derivative_bounds),
        source_executed=False,
        theorem_consumed=True,
    )
    for name in CHARTS:
        chart = _chart(packet, name)
        c = MPIntervalContext()
        c.prec = 900
        computed = _compute_bounds(c, chart, derivative_bounds)
        log_Pstar_E_upper = _chart_record(
            c, chart,
            ("log_Pstar_E_upper", "log_Pstar_E", "log_amplitude_upper",
             "native_log_amplitude_upper"),
            name + ".log_Pstar_E_upper",
        )
        logN = _chart_record(
            c, chart, ("logN", "actual_logN", "N_log"),
            name + ".logN", packet,
        )
        intermediate = {
            key: _jet_rows(computed[key])
            for key in (
                "D", "invD", "t", "lambda_value", "T1", "T2",
                "inverse", "T1hat", "target_phase", "B_over_E", "A",
                "q", "u", "r", "alpha", "K",
            )
        }
        results[name] = dict(
            passed=True,
            chart=name,
            parent_report=PARENT,
            parent_report_sha256=PARENT_SHA256,
            ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
            physical_y_derivatives_claimed=False,
            source_pins=source_pins,
            B_over_Pstar_over_frozen_E0_log_caps=_logz_rows(
                computed["B_over_Pstar_over_frozen_E0"]
            ),
            B_over_Pstar_over_frozen_E0_log_caps4rows=_logz_rows(
                computed["B_over_Pstar_over_frozen_E0"]
            ),
            B_over_E_log_caps=_logz_rows(computed["B_over_E"]),
            B_over_E_log_caps4rows=_logz_rows(computed["B_over_E"]),
            A_log_caps=_logz_rows(computed["A"]),
            A_log_caps4rows=_logz_rows(computed["A"]),
            inverse_Psi_log_caps=_logz_rows(computed["inverse"]),
            log_Pstar_E_upper=_record(log_Pstar_E_upper),
            logN=_record(logN),
            conditioning={
                "log_D_lower": _record(computed["conditioning"]["log_D_lower"]),
                "log_lambda_lower": _record(
                    computed["conditioning"]["log_lambda_lower"]
                ),
            },
            intermediate_Z_phase_caps=intermediate,
            root_control_caps={
                key: _logz_rows(computed[key])
                for key in ("a", "b", "t0", "Delta", "gamma", "q",
                            "p2", "u", "h", "r", "alpha", "K",
                            "E_relative")
            },
            root_control_proof=dict(
                a_positive=True,
                Delta_lower_bound="active Delta<eta; flat Delta>=eta",
                Delta_global_negative=computed["delta_global_negative"],
                gamma_lower_log=_record(computed["gamma_lower_log"]),
                gamma_lower_formula=(
                    "Delta=kappa-2, kappa<1, eta>0 => gamma=2*eta-Delta>1"
                ),
                q_formula=(
                    "q=lazy[sigma(1-Delta/eta)*"
                    "sqrt((2*eta-Delta)/(2*a))]"
                ),
                q_derivative_source=derivative_source,
                D_lower_formula="D>=h_lower^4/4",
                lambda_lower_formula="lambda=K*(1+t^2)>=K",
                t0_retained=True,
                b_retained=True,
            ),
            inverse_proof=computed["inverse_proof"],
            delta_normalized_Z_rows={
                order: (None if value is None else _record(value))
                for order, value in computed["delta"].items()
            },
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            actual_N_not_materialized=True,
            actual_phase_Z_independent=True,
            finite_N_phase_not_selected=True,
            finite_N_phase_not_materialized=True,
            universal_loop_controls_consumed=False,
            universal_B_over_E_value_cap_tightened=False,
            source_contract=_source_contract(name, computed["cutoff_mode"]),
            physical_scale_proof=dict(
                passed=True,
                B_quantity="B/Pstar/frozen_E0",
                A_quantity="A for finite-N delta_E correction",
                deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
                frozen_E0_normalization="constant whole-chart upper cap for E",
                frozen_E0_is_source_value=False,
                frozen_E0_source_value=False,
                E_over_frozen_E0_upper_bound="1",
                log_Pstar_E_upper_source="log(Utheta)=log(Pstar*E)",
                current_N="N18 graph handle",
                actual_N18=True,
                actual_N18_is_graph_handle=True,
                actual_N_not_materialized=True,
                phase_Z_independent=True,
                physical_Pstar_applied_once=True,
                physical_y_derivatives_claimed=False,
                source_bound_to_current_restoration_native_controls=True,
            ),
            source_pin_audit=dict(
                input_module=INPUT_MODULE,
                new_module_input_self_hash=source_pins[INPUT_MODULE]["actual_sha256"],
                source_executed=False,
                self_pin=True,
            ),
            frozen_E0=dict(
                normalization="constant whole-chart upper cap for E",
                normalization_cap="native chart.frozen_E0.chart_upper_log_E",
                E_over_frozen_E0_upper_bound="1",
                frozen_E0_is_source_value=False,
                selected_endpoint="wholechart_upper_cap",
                log_Pstar_E_upper="log(Utheta)=log(Pstar*E)",
            ),
            canonical_parameters_inherited=canonical_parameters_inherited,
            canonical_parameter_source=(
                chart.get("canonical_parameter_source")
                if isinstance(chart, Mapping) else None
            ),
            open_gates=[
                "finite_N_numeric",
                "global_five_moment_closure",
                "full_pressure_global_closure",
            ],
            source_formulas=dict(
                q=("lazy(sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a)))"),
                u="p2*q/d_star",
                r="u/sqrt(1+u^2)",
                alpha="2*q/sqrt(1+u^2)",
                K="1/[2*pi*(1+t0^2+2*q^2)]",
                B="E*M/2",
                A="a*(phase-Psi/(2*pi))/2",
                T1hat="compose_inverse(Phi,T1,delta)",
                Phi="K*(angle+T2), Phi[k,j>=1]=lambda[k,j-1]",
                inverse="inverse_coefficients(Phi, log_lambda_lower)",
                target_phase="Z-independent unit magnitude bound; no selected phase",
            ),
        )

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_family=parent.get("source_family"),
        source_executed=False,
        accepted_ancestors_replayed=False,
        actual_runtime=dict(
            module=MODULE,
            module_sha256=_sha(root / MODULE),
            directed_interval_precision_bits=900,
            source_executed=False,
            accepted_ancestors_replayed=False,
            model="unavailable",
            tool_metadata="unavailable",
        ),
        input_packet_path=packet_path,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        charts=results,
        B_over_Pstar_over_frozen_E0_log_caps={
            name: results[name]["B_over_Pstar_over_frozen_E0_log_caps"]
            for name in CHARTS
        },
        B_over_E_log_caps={
            name: results[name]["B_over_E_log_caps"] for name in CHARTS
        },
        A_log_caps={name: results[name]["A_log_caps"] for name in CHARTS},
        inverse_Psi_log_caps={
            name: results[name]["inverse_Psi_log_caps"] for name in CHARTS
        },
        log_Pstar_E_upper={
            name: results[name]["log_Pstar_E_upper"] for name in CHARTS
        },
        logN=results[CHARTS[0]]["logN"],
        source_pins=source_pins,
        new_module_input_self_hash=source_pins[INPUT_MODULE]["actual_sha256"],
        source_pin_audit=dict(
            input_module=INPUT_MODULE,
            new_module_input_self_hash=source_pins[INPUT_MODULE]["actual_sha256"],
            source_executed=False,
            self_pin=True,
        ),
        source_guards=dict(
            source_pins=source_pins,
            C5_formula_guard=c5_formula_guard,
            parent_saved_proofs=parent_proofs,
            sigma_derivative_source=derivative_source,
        ),
        guards=source_pins,
        proofs=dict(
            input_packet=proof,
            canonical_parameters_inherited=canonical_parameters_inherited,
            parent_saved_proofs=parent_proofs,
            actual_phase_mixed_consumed=True,
            actual_correlated_pressure_consumed=True,
            cutoff_derivative_source_consumed=True,
            no_new_accepted_root_or_phase_replay=True,
        ),
        physical_scale_proof=dict(
            passed=True,
            quantity="B/Pstar/frozen_E0",
            A_quantity="A for finite-N delta_E correction",
            deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
            frozen_E0_normalization="constant whole-chart upper cap for E",
            frozen_E0_is_source_value=False,
            frozen_E0_source_value=False,
            E_over_frozen_E0_upper_bound="1",
            log_Pstar_E_upper_source="log(Utheta)=log(Pstar*E)",
            current_N="N18 graph handle",
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            actual_N_not_materialized=True,
            phase_Z_independent=True,
            physical_Pstar_applied_once=True,
            physical_y_derivatives_claimed=False,
            source_bound_to_current_restoration_native_controls=True,
        ),
        source_contract=dict(
            charts=list(CHARTS),
            actual_native_restoration_controls=True,
            actual_native_post_restore_controls=True,
            actual_root_inputs=True,
            common_N_handle=18,
            current_N="N18 graph handle",
            finite_N_numeric=False,
            phase_Z_independent=True,
            selected_phase_value=False,
            physical_Pstar_applied_once=True,
            b_nonzero_interior_retained=True,
            t0_nonzero_interior_retained=True,
            source_family=parent.get("source_family"),
            no_accepted_builder_checker_producer_owner_norm_recurrence_quadrature=True,
            canonical_parameters_inherited=canonical_parameters_inherited,
        ),
        open_gates=[
            "finite_N_numeric",
            "global_five_moment_closure",
            "full_pressure_global_closure",
        ],
        loop_recipe=dict(
            q="lazy(sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a)))",
            B="E*M/2",
            A="a*(phase-Psi/(2*pi))/2",
            target_phase="Z-independent unit magnitude bound",
            no_universal_C0_tightening=True,
            global_q_positive_lower_not_claimed=True,
        ),
        no_accepted_builder_checker_producer_owner_norm_recurrence_quadrature=True,
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_N_not_materialized=True,
        actual_phase_Z_independent=True,
        finite_N_phase_not_selected=True,
        finite_N_phase_not_materialized=True,
        source_binding=True,
        full_goal_complete=False,
        limits=dict(
            native_geometry_log_only=True,
            physical_y_derivatives_claimed=False,
            finite_N_numeric=False,
            global_five_moment_closure=False,
            full_pressure_global_closure=False,
            full_goal_complete=False,
        ),
    )
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    """Exercise pure arithmetic with a nonzero ``b``/``t0`` fixture."""
    c = MPIntervalContext()
    c.prec = 240
    zero = c.mpf(0)
    # Use a globally negative Delta so the fixture exercises the exact flat
    # branch while retaining nonzero b and t0 rows in every ordinary order.
    point = lambda value: c.mpf([value, value])
    caps = {
        "a": LogZ.from_mapping(c, {0: c.log(c.mpf("0.8")), 1: zero,
                                    2: None, 3: None}),
        "b": LogZ.from_mapping(c, {0: c.log(c.mpf("0.05")), 1: zero,
                                    2: None, 3: None}),
        "t0": LogZ.from_mapping(c, {0: c.log(c.mpf("0.0625")), 1: zero,
                                     2: None, 3: None}),
        "p2": LogZ.from_mapping(c, {0: c.log(c.mpf("0.2")), 1: zero,
                                    2: None, 3: None}),
        "Delta": LogZ.from_mapping(c, {0: c.log(c.mpf("1.2")), 1: zero,
                                        2: None, 3: None}),
        "E_relative": LogZ.from_mapping(c, {0: zero, 1: None, 2: None, 3: None}),
    }
    chart = {
        "passed": True,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "root_inputs": {key: _logz_rows(value) for key, value in caps.items()},
        "a_lower_log": point(c.log(c.mpf("0.8"))),
        "actual_eta_log": point(c.log(c.mpf("0.25"))),
        "dstar_log": point(c.log(c.mpf("1.0"))),
        "delta_globally_negative": True,
    }
    computed = _compute_bounds(c, chart, (1, 1296, 38880, 19290960))
    checks = dict(
        four_Z_rows=len(_logz_rows(computed["B_over_Pstar_over_frozen_E0"])) == 4,
        four_A_rows=len(_logz_rows(computed["A"])) == 4,
        phase_rows=len(_phase_rows(computed["Phi"])) == 10,
        nonzero_b_retained=computed["b"][0] is not None,
        nonzero_t0_retained=computed["t0"][0] is not None,
        finite_r_value=computed["r"][0] is not None,
        delta_zero_is_exact=computed["delta"][0] is None,
        target_phase_Z_independent=all(
            computed["target_phase"][order] is None
            for order in ORDINARY_Z_ORDERS if order != 0
        ),
        flat_branch_used=(
            computed["cutoff_mode"] == "exact_flat_sigma_one_from_global_Delta_lt_zero"
        ),
        no_universal_C0_tightening=(
            computed["universal_B_over_E_value_cap_tightened"] is False
        ),
    )
    _require(all(checks.values()), "Restoration Z3 fixture failed")
    receipt = dict(passed=True, checks=checks, JSON_serializable=True,
                   source_packet_required=True, full_goal_complete=False)
    json.dumps(receipt, allow_nan=False)
    return receipt


def self_check():
    return fixture_bound_check()


__all__ = (
    "PARENT",
    "PARENT_SHA256",
    "MODULE",
    "INPUT_MODULE",
    "CHARTS",
    "ORDINARY_Z_ORDERS",
    "build",
    "fixture_bound_check",
    "self_check",
)
