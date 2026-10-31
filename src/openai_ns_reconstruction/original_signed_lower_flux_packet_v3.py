"""Live signed lower-flux evaluator with the physical ``V_y`` normalization.

This module consumes an already-built live ``source_rows`` view.  It never
opens a source packet, calls a producer, walks a graph, or turns a source cap
into a point value.  ``V`` and ``m`` are the normalized generic rows, so their
physical rows are ``Pstar*V`` and ``Pstar*m``.  In contrast,
``actual_background_source.V_y`` is already the raw physical ``U_y`` source;
its first two rows are used directly.  Applying ``Pstar`` to those rows again
would be the archived endpoint's duplicate-scale defect.

The source certificates below are caller-supplied evidence tied to the exact
live packet and target objects.  A certificate check is reported as a
conditional admission, never as an independently established fact merely
because a caller supplied ``True``.  The two error layers are carried
verbatim, and their presence is reported dynamically.

The physical source expression is

    L = 1 - delta*Z**2
    Ny = 2*Z*Uy0 - (1-delta)*Z*A - (1-Z**2)*Az
    RQy = R*Ny/L

with its ordinary first ``Z`` derivative retained.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_signed_lower_flux_packet_v3"
COMPONENT = "original_signed_lower_flux_packet"
SOURCE_QUERY_CHART = "first_micro"
SOURCE_QUERY_ENDPOINT = ("sc_half",)


class SignedLowerFluxError(ValueError):
    """The supplied live source view or source conditions cannot be admitted."""


OriginalSignedLowerFluxError = SignedLowerFluxError


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SignedLowerFluxError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    return value


def _sequence(value: Any, label: str) -> Sequence[Any]:
    _require(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
        f"{label} must be a sequence",
    )
    return value


def _attr(value: Any, name: str, label: str) -> Any:
    try:
        result = getattr(value, name)
    except AttributeError as exc:
        raise SignedLowerFluxError(f"{label}.{name} is required") from exc
    return result


def _same_identity(value: Any, ctx: Any, flow: Any, label: str) -> None:
    """Require the exact live context, scale basis, and source ledger."""

    _require(_attr(value, "ctx", label) is ctx, f"{label} uses a different ctx")
    scale = _attr(value, "scale", label)
    _require(
        _attr(scale, "bases", label + ".scale") is _attr(flow, "logs", "flow"),
        f"{label} uses a different flow.logs basis",
    )
    _require(
        _attr(value, "ledger", label) is _attr(flow, "ledger", "flow"),
        f"{label} uses a different flow.ledger",
    )


def _source_rows(view: Any) -> Mapping[str, Any]:
    return _mapping(view, "source_rows")


def _source_packet(view: Mapping[str, Any]) -> Mapping[str, Any]:
    packet = view.get("source_packet")
    _require(packet is not None, "source_rows.source_packet is required")
    packet = _mapping(packet, "source_rows.source_packet")
    _require(packet.get("source_identity") is not None,
             "source_packet.source_identity is required")
    return packet


def _packet_background(packet: Mapping[str, Any]) -> tuple[Mapping[str, Any], Sequence[Any]]:
    background = _mapping(packet.get("actual_background_source"),
                          "source_packet.actual_background_source")
    values = _sequence(background.get("V_y"),
                       "source_packet.actual_background_source.V_y")
    _require(len(values) >= 2, "actual_background_source.V_y rows 0 and 1 are required")
    return background, values


def _canonical_endpoint(value: Any, label: str) -> tuple[Any, ...]:
    values = _sequence(value, label)
    _require(len(values) == 1 and values[0] == SOURCE_QUERY_ENDPOINT[0],
             f"{label} must be the exact ('sc_half',) source endpoint")
    return tuple(values)


def _live_source_links(view: Mapping[str, Any], packet: Mapping[str, Any],
                       V: Sequence[Any], m: Sequence[Any],
                       vy: Sequence[Any], op: Any) -> None:
    """Bind live rows to exact packet leaves, not merely a family name."""

    generic = _mapping(packet.get("original_generic_source"),
                       "source_packet.original_generic_source")
    _require(generic.get("common_velocity_V_axial5") is V,
             "source_rows.V is not the packet's common velocity rows")
    histories = _mapping(generic.get("common_own_five_histories_axial5"),
                         "source_packet.original_generic_source.common_own_five_histories_axial5")
    _require(histories.get("m") is m,
             "source_rows.m is not the packet's original m history")
    _require(generic.get("original_physical_radius") is view.get("R"),
             "source_rows.R is not the packet's original physical radius")
    actual_background = _mapping(packet.get("actual_background_source"),
                                  "source_packet.actual_background_source")
    _require(actual_background.get("V_y") is vy,
             "source_rows.V_y is not packet.actual_background_source.V_y")

    if hasattr(op, "P0"):
        _require("exact_common_P0_axial5" in packet,
                 "source_packet exact_common_P0_axial5 is required")
        _require(packet.get("exact_common_P0_axial5") is op.P0,
                 "source_packet exact_common_P0_axial5 is not op.P0")


def _validate_query(packet: Mapping[str, Any], assumptions: Mapping[str, Any],
                    Z: Any) -> dict[str, Any]:
    _require(packet.get("chart") == SOURCE_QUERY_CHART,
             "source packet must be queried on first_micro")
    left = _canonical_endpoint(packet.get("exact_left"),
                               "source_packet.exact_left")
    right = _canonical_endpoint(packet.get("exact_right"),
                                "source_packet.exact_right")
    _require(left == right == SOURCE_QUERY_ENDPOINT,
             "source packet must be the exact sc_half lower query")
    _require(packet.get("actual_phase_Z_exact_zero") is True,
             "source packet does not certify its phase is Z independent")
    _require("exact_Z_cell" in packet,
             "source packet exact_Z_cell is required")
    candidate_N = packet.get("candidate_N")
    _require(candidate_N is not None, "source packet candidate_N is required")

    cert = _certificate(assumptions, ("query", "query_certificate", "source_query"),
                        packet, "source query")
    checks = _mapping(cert.get("checks"), "source query certificate.checks")
    for name in ("packet_query_exact", "coordinate_exact",
                 "candidate_N_is_source_owned", "Z_matches_exact_Z_cell"):
        _require(checks.get(name) is True,
                 f"source query certificate lacks {name}")
    _require(cert.get("Z") is Z,
             "source query certificate is not tied to this target Z interval")
    details = cert.get("details")
    if details is not None:
        details = _mapping(details, "source query certificate.details")
        if "chart" in details:
            _require(details.get("chart") == SOURCE_QUERY_CHART,
                     "source query certificate chart changed")
        if "exact_left" in details:
            _require(_canonical_endpoint(details.get("exact_left"),
                                         "source query certificate.details.exact_left") == left,
                     "source query certificate left coordinate changed")
        if "exact_right" in details:
            _require(_canonical_endpoint(details.get("exact_right"),
                                         "source query certificate.details.exact_right") == right,
                     "source query certificate right coordinate changed")
        if "candidate_N" in details:
            _require(details.get("candidate_N") == candidate_N,
                     "source query certificate candidate_N changed")
        if "exact_Z_cell" in details:
            _require(details.get("exact_Z_cell") == packet.get("exact_Z_cell"),
                     "source query certificate exact_Z_cell changed")
    return {"chart": SOURCE_QUERY_CHART, "exact_left": list(left),
            "exact_right": list(right), "candidate_N_present": True}


def _certificate(assumptions: Mapping[str, Any], names: Sequence[str],
                 packet: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    cert = next((assumptions.get(name) for name in names
                 if isinstance(assumptions.get(name), Mapping)), None)
    _require(cert is not None, f"{label} certificate is required")
    _require(cert.get("source_packet") is packet,
             f"{label} certificate is not tied to this live source_packet")
    source = cert.get("source")
    _require(isinstance(source, str) and bool(source),
             f"{label} certificate needs a non-empty source path")
    return cert


def _validate_flat_collar_and_history(packet: Mapping[str, Any],
                                      m: Sequence[Any],
                                      assumptions: Mapping[str, Any]) -> dict[str, Any]:
    support = packet.get("source_owned_modification_support")
    _require(support in ("unchanged_left", "original_flat_left_collar"),
             "source packet is not on the source-owned flat collar")
    _require(packet.get("source_owned_collar_q_C0_Z_exact_zero") is True,
             "source-owned flat-collar C0 flag is missing")
    _require(packet.get("original_active_inverse_not_evaluated_on_source_owned_flat_collar") is True,
             "source-owned flat collar did not retain the original branch")
    initial = _mapping(packet.get("actual_source_owned_initial_condition"),
                       "source_packet.actual_source_owned_initial_condition")
    for name in (
        "exact_zero_correction_follows_from_defining_Duhamel_lower_bound",
        "actual_background_incoming_histories_remain_nonzero",
        "q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar",
        "no_saved_frame_or_N1024_initial_correction_used",
    ):
        _require(initial.get(name) is True,
                 f"source-owned collar/history certificate lacks {name}")
    _require(packet.get("original_background_nonzero_histories_remain_separate", True) is True,
             "original nonzero background history was reset or collapsed")

    # A zero row can be a valid exact source value.  It is intentionally not
    # evidence for, or against, historical preservation.  Only object identity
    # and the explicit source-owned certificate establish that provenance.
    zero_flags = [getattr(row, "zero") for row in m if hasattr(row, "zero")]
    cert = _certificate(assumptions, ("history", "history_certificate"),
                        packet, "original m history")
    checks = _mapping(cert.get("checks"), "original m history certificate.checks")
    for name in ("packet_m_rows_exact", "original_history_not_zeroed",
                 "source_owned_flat_collar"):
        _require(checks.get(name) is True,
                 f"original m history certificate lacks {name}")
    return {
        "support": support,
        "source_owned_flat_collar": True,
        "original_m_history_not_zeroed": True,
        "zero_rows_are_not_provenance_evidence": True,
        "m_rows_exact_zero_count": sum(flag is True for flag in zero_flags),
    }


def _validate_z_independence(assumptions: Mapping[str, Any], packet: Mapping[str, Any],
                             Z: Any, delta: Any) -> dict[str, Any]:
    cert = _certificate(assumptions,
                        ("z_independence", "z_independence_certificate", "geometry"),
                        packet, "Z-independence")
    checks = _mapping(cert.get("checks"), "Z-independence certificate.checks")
    _require(cert.get("Z") is Z,
             "Z-independence certificate is not tied to this target Z interval")
    _require(cert.get("delta") is delta,
             "Z-independence certificate is not tied to this delta")
    for name in ("R_Z_independent", "Pstar_Z_independent", "delta_Z_independent",
                 "full_target_Z_interval"):
        _require(checks.get(name) is True,
                 f"Z-independence certificate lacks {name}")
    return {
        "status": "conditionally_admitted_from_caller_certificate",
        "certificate_checks": tuple(
            name for name in ("R_Z_independent", "Pstar_Z_independent",
                              "delta_Z_independent", "full_target_Z_interval")
        ),
    }


def _lower_bound(value: Any) -> Any:
    """Read a lower endpoint without evaluating a source enclosure."""

    for name in ("lo", "lower", "a"):
        if hasattr(value, name):
            endpoint = getattr(value, name)
            if callable(endpoint):
                endpoint = endpoint()
            return endpoint
    if isinstance(value, (int, float)):
        return value
    raise SignedLowerFluxError("cannot inspect the requested Z interval lower endpoint")


def _strict_positive(value: Any, label: str) -> Any:
    lower = _lower_bound(value)
    try:
        ok = bool(lower > 0)
    except Exception as exc:
        raise SignedLowerFluxError(f"cannot prove strict positivity of {label}") from exc
    _require(ok, f"{label} is not strictly positive on the full requested interval")
    return lower


def _validate_L_certificate(assumptions: Mapping[str, Any], packet: Mapping[str, Any],
                            Z: Any, delta: Any, L: Any) -> dict[str, Any]:
    lower = _strict_positive(L, "L=1-delta*Z^2")
    cert = _certificate(assumptions,
                        ("L_positive", "L_positive_certificate", "domain"),
                        packet, "L positivity")
    _require(cert.get("Z") is Z and cert.get("delta") is delta,
             "L positivity certificate is not tied to this Z/delta input")
    checks = _mapping(cert.get("checks"), "L positivity certificate.checks")
    _require(checks.get("L_strictly_positive") is True,
             "L positivity certificate lacks L_strictly_positive")
    if "lower_bound" in cert:
        try:
            _require(cert.get("lower_bound") <= lower,
                     "L positivity certificate lower bound is not conservative")
        except TypeError as exc:
            raise SignedLowerFluxError("L positivity lower bound is not comparable") from exc
    return {"strictly_positive": True, "lower_bound": lower}


def _row(value: Any, index: int, ctx: Any, flow: Any, label: str) -> Any:
    _same_identity(value, ctx, flow, f"{label}[{index}]")
    return value


def _validate_live_inputs(source_rows: Any, op: Any, ctx: Any, Z: Any, delta: Any,
                          assumptions: Any) -> tuple[Mapping[str, Any], Any,
                                                      Mapping[str, Any], Sequence[Any],
                                                      Sequence[Any], Sequence[Any]]:
    view = _source_rows(source_rows)
    _require(op is not None, "op is required")
    flow = _attr(op, "flow", "op")
    _require(flow is not None, "op.flow is required")
    _require(_attr(flow, "c", "flow") is ctx, "flow.c is not the supplied ctx")
    if hasattr(op, "ctx"):
        _require(op.ctx is ctx, "op.ctx is not the supplied ctx")
    if hasattr(op, "c"):
        _require(op.c is ctx, "op.c is not the supplied ctx")
    if hasattr(Z, "ctx"):
        _require(Z.ctx is ctx, "Z interval uses a different ctx")
    if hasattr(delta, "ctx"):
        _require(delta.ctx is ctx, "delta uses a different ctx")
    _mapping(assumptions, "assumptions")
    packet = _source_packet(view)
    _background, vy = _packet_background(packet)
    V = _sequence(view.get("V"), "source_rows.V")
    m = _sequence(view.get("m"), "source_rows.m")
    _require(len(V) >= 3, "source_rows.V rows 0 through 2 are required")
    _require(len(m) >= 3, "source_rows.m rows 0 through 2 are required")
    for index, value in enumerate(V):
        _row(value, index, ctx, flow, "source_rows.V")
    for index, value in enumerate(m):
        _row(value, index, ctx, flow, "source_rows.m")
    for index, value in enumerate(vy):
        _row(value, index, ctx, flow,
             "source_packet.actual_background_source.V_y")
    for name in ("R", "Pstar"):
        _same_identity(view.get(name), ctx, flow, "source_rows." + name)
    _require(hasattr(op, "P0"), "op.P0 is required for exact source query binding")
    _require(hasattr(op, "Pstar"), "op.Pstar is required for exact source scale binding")
    _require(view.get("Pstar") is op.Pstar,
             "source_rows.Pstar is not the exact op.Pstar object")
    _live_source_links(view, packet, V, m, vy, op)
    _validate_query(packet, _mapping(assumptions, "assumptions"), Z)
    _validate_z_independence(_mapping(assumptions, "assumptions"), packet, Z, delta)
    _validate_flat_collar_and_history(packet, m, _mapping(assumptions, "assumptions"))
    return view, flow, packet, V, m, vy


def _error_report(view: Mapping[str, Any], packet: Mapping[str, Any]) -> dict[str, Any]:
    """Describe both opaque error layers without inventing propagation claims."""

    source_errors = view.get("errors")
    packet_errors = packet.get("errors")
    source_present = source_errors is not None
    packet_present = packet_errors is not None
    return {
        "source_rows_present": source_present,
        "source_packet_present": packet_present,
        "both_layers_present": source_present and packet_present,
        "source_rows": source_errors,
        "source_packet": packet_errors,
    }


def evaluate(source_rows: Mapping[str, Any], op: Any, ctx: Any, Z: Any, delta: Any,
             *, assumptions: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Evaluate physical lower ``R*Q_y`` and ordinary ``Z1``.

    ``V`` and ``m`` are normalized generic source rows and receive ``Pstar``.
    ``actual_background_source.V_y`` is already physical, so ``Uy0`` and
    ``Uy1`` are the raw rows themselves.  This distinction is retained by the v3 repair.
    """

    view, flow, packet, V, m, vy = _validate_live_inputs(
        source_rows, op, ctx, Z, delta, assumptions,
    )
    assumptions_map = _mapping(assumptions, "assumptions")

    Pstar = view["Pstar"]
    R = view["R"]
    U0, U1, U2 = Pstar * V[0], Pstar * V[1], Pstar * V[2] * 2
    M0, M1, M2 = Pstar * m[0], Pstar * m[1], Pstar * m[2] * 2
    # ``vy`` is actual_background_source.V_y and is already physical Uy.
    # Do not multiply these rows by Pstar a second time.
    Uy0, Uy1 = vy[0], vy[1]
    A, Az, Azz = U0 - M0, U1 - M1, U2 - M2

    one = ctx.mpf(1) if hasattr(ctx, "mpf") else 1
    Z2 = Z * Z
    one_minus_Z2 = one - Z2
    one_minus_delta = one - delta
    L = one - delta * Z2
    L_info = _validate_L_certificate(assumptions_map, packet, Z, delta, L)
    try:
        inverse_L = one / L
    except Exception as exc:
        raise SignedLowerFluxError("positive numeric reciprocal for L is unavailable") from exc

    two_Z = Z * 2
    two_delta_Z = delta * Z * 2
    Ny = Uy0 * two_Z - A * (one_minus_delta * Z) - Az * one_minus_Z2
    NyZ = (Uy0 * 2 + Uy1 * two_Z
           - (A + Az * Z) * one_minus_delta
           + Az * two_Z - Azz * one_minus_Z2)
    RQy = R * Ny * inverse_L
    dRQy = R * (NyZ * L + Ny * two_delta_Z) * inverse_L * inverse_L
    errors = _error_report(view, packet)

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_source_enclosure",
        "source_enclosures": {
            "U": (U0, U1, U2),
            "m": (M0, M1, M2),
            "U_y": (Uy0, Uy1),
            "A": A,
            "A_Z": Az,
            "A_ZZ": Azz,
            "N_y": Ny,
            "N_y_Z": NyZ,
            "R_Q_y": RQy,
            "dZ_R_Q_y": dRQy,
            "ordinary_Z1": dRQy,
        },
        "RQy": RQy,
        "R_Q_y": RQy,
        "ordinary_Z1": dRQy,
        "dZ_RQy": dRQy,
        "denominator": L,
        "proof": {
            "same_live_flow_context_basis_ledger": True,
            "exact_source_packet_and_query": True,
            "source_owned_flat_collar": True,
            "original_m_history_not_zeroed": True,
            "zero_m_rows_not_used_as_history_evidence": True,
            "physical_Uy_rows_used_without_Pstar_rescale": True,
            "normalized_V_and_m_receive_Pstar_once": True,
            "R_Pstar_delta_Z_independence": "conditional_on_caller_certificate",
            "L_strictly_positive_on_full_Z": True,
            "Taylor_second_row_factorial_converted_once": True,
            "signed_physical_R_times_Qy": True,
            "ordinary_Z1_retained": True,
            "source_enclosures_not_point_values": True,
            "candidate_N_source_certificate": True,
            "native_or_corrected_C1_identity_proved": False,
            "graph_node_identity_proved": False,
        },
        "scope": {
            "conditional_on_source_certificates": True,
            "target_Z_is_full_requested_interval": True,
            "R_Pstar_delta_Z_independence": "conditional_source_certificate",
            "L_positive_lower_bound": L_info["lower_bound"],
            "query": {
                "chart": SOURCE_QUERY_CHART,
                "exact_left": list(SOURCE_QUERY_ENDPOINT),
                "exact_right": list(SOURCE_QUERY_ENDPOINT),
                "exact_Z_cell": packet.get("exact_Z_cell"),
                "Z_packet_coordinate_match": "conditional_source_certificate",
            },
            "candidate_N_not_materialized": True,
            "error_layers": errors,
            "declared_error_layers_retained": errors["both_layers_present"],
            "all_inherited_errors_verified_preserved": False,
            "retained_error_terms": {
                "source_rows": errors["source_rows"],
                "source_packet": errors["source_packet"],
            },
            "native_RQy_identity_proved": False,
            "corrected_RQy_identity_proved": False,
        },
    }


def evaluate_signed_lower_flux(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Descriptive alias for :func:`evaluate`."""

    return evaluate(*args, **kwargs)


def consume(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return evaluate(*args, **kwargs)


__all__ = [
    "SCHEMA", "COMPONENT", "SignedLowerFluxError", "OriginalSignedLowerFluxError",
    "evaluate", "evaluate_signed_lower_flux", "consume",
]

