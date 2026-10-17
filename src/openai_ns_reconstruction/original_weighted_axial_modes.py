"""Hydrate the signed axial ordinary modes from accepted source packets.

This module is intentionally a saved-data adapter.  It does not construct an
ancestor source, recompute the accepted mass, or materialize the positive
swirl amplitude.  The pressure rows are already macro-enclosed.  The two
factored hydro sectors and the three swirl sectors are evaluated here with
ordinary order-two arithmetic and then receive their own pure macro and
positive-amplitude enclosure steps.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .original_leading_viscous_kernel import decode_pair, encode_pair
from .original_signed_C1_source import macro_nonlinear_enclosure
from .original_source_jet_kernel import ordinary_product, ordinary_quotient
from .original_swirl_source_point import scaled_atom_enclosure


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _exact_interval(ctx: Any, iv: Any, value: Any) -> Any:
    """Hydrate an interval from live arithmetic or a saved exact packet."""

    if hasattr(value, "_mpi_"):
        return value
    if hasattr(value, "_mpf_"):
        return iv.mpf(ctx.make_mpf(value._mpf_))
    if isinstance(value, Mapping):
        endpoints = value.get("source_enclosure_exact_mpf_tuples")
        if endpoints is None:
            endpoints = value.get("exact_mpf_tuples")
        if endpoints is not None:
            if not _is_sequence(endpoints) or len(endpoints) != 2:
                raise ValueError("A saved interval requires two exact endpoints")
            return iv.mpf(
                [ctx.make_mpf(tuple(endpoints[0])), ctx.make_mpf(tuple(endpoints[1]))]
            )
        lower_tuple = value.get("lower_exact_mpf_tuple")
        upper_tuple = value.get("upper_exact_mpf_tuple")
        if lower_tuple is not None or upper_tuple is not None:
            if lower_tuple is None or upper_tuple is None:
                raise ValueError("Both saved lower and upper exact endpoints are required")
            return iv.mpf(
                [ctx.make_mpf(tuple(lower_tuple)), ctx.make_mpf(tuple(upper_tuple))]
            )
        if "exact_mpf_tuple" in value:
            point = ctx.make_mpf(tuple(value["exact_mpf_tuple"]))
            return iv.mpf(point)
        if "lower" in value and "upper" in value:
            return iv.mpf([ctx.mpf(value["lower"]), ctx.mpf(value["upper"])])
        raise ValueError("Unsupported saved scalar or interval packet")
    if _is_sequence(value):
        if len(value) != 2:
            raise ValueError("An interval sequence requires two endpoints")
        if all(_is_sequence(v) and len(v) == 4 for v in value):
            return iv.mpf([ctx.make_mpf(tuple(value[0])), ctx.make_mpf(tuple(value[1]))])
        return iv.mpf([ctx.mpf(value[0]), ctx.mpf(value[1])])
    return iv.mpf(value)


def _serial(value: Any) -> Any:
    """Return JSON-safe exact endpoint evidence without selecting an interval value."""

    if hasattr(value, "_mpi_"):
        return {"directed_exact_mpf_tuples": [list(v) for v in value._mpi_]}
    if hasattr(value, "_mpf_"):
        return {"exact_mpf_tuple": list(value._mpf_)}
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_serial(item) for item in value]
    return value


def _rows(report: Mapping[str, Any], *paths: tuple[str, ...]) -> list[Mapping[str, Any]]:
    for path in paths:
        current: Any = report
        for key in path:
            if not isinstance(current, Mapping) or key not in current:
                current = None
                break
            current = current[key]
        if _is_sequence(current):
            return list(current)
    raise ValueError("Saved source rows are missing")


def _find_row(rows: Sequence[Mapping[str, Any]], ztext: str) -> tuple[int, Mapping[str, Any]]:
    for index, row in enumerate(rows):
        if row.get("Z") == ztext or str(row.get("Z")) == ztext:
            return index, row
    raise ValueError(f"Saved source row for Z={ztext!r} is missing")


def _decode_jet(ctx: Any, iv: Any, rows: Any, *, name: str) -> tuple[list[Any], list[Any]]:
    if not _is_sequence(rows) or len(rows) != 3:
        raise ValueError(f"{name} must retain ordinary Z0/Z1/Z2 rows")
    pairs = [decode_pair(ctx, iv, row) for row in rows]
    return [pair[0] for pair in pairs], [pair[1] for pair in pairs]


def _factored_entries(
    ctx: Any,
    iv: Any,
    row: Mapping[str, Any],
    key: str,
    expected: set[int],
) -> dict[int, tuple[list[Any], list[Any]]]:
    entries = row.get(key)
    if not _is_sequence(entries):
        raise ValueError(f"{key} must retain every saved epsilon sector")
    result: dict[int, tuple[list[Any], list[Any]]] = {}
    for entry in entries:
        if not isinstance(entry, Mapping) or "epsilon_degree" not in entry:
            raise ValueError(f"Malformed {key} epsilon sector")
        degree = entry["epsilon_degree"]
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise ValueError(f"Malformed {key} epsilon degree")
        if degree in result:
            raise ValueError(f"Duplicate {key} epsilon degree {degree}")
        result[degree] = _decode_jet(ctx, iv, entry.get("ordinary_rows"), name=f"{key}[{degree}]")
    if set(result) != expected:
        raise ValueError(f"{key} omitted or added epsilon sectors: {sorted(result)}")
    return result


def _evaluate_sectors(
    ctx: Any,
    entries: Mapping[int, tuple[list[Any], list[Any]]],
    epsilon: Any,
    *,
    truth: bool,
    factor: Any = 1,
    exponent_shift: int = 0,
) -> list[Any]:
    """Evaluate a saved epsilon polynomial in either point or interval context."""

    zero = epsilon * 0
    result = [zero, zero, zero]
    for degree, (point_rows, truth_rows) in entries.items():
        rows = truth_rows if truth else point_rows
        scalar = factor * epsilon ** (degree + exponent_shift)
        result = [value + scalar * rows[order] for order, value in enumerate(result)]
    return result


def _negative(row: Sequence[Any]) -> list[Any]:
    return [-value for value in row]


def _saved_q(
    ctx: Any,
    iv: Any,
    pressure_report: Mapping[str, Any],
    pressure_index: int,
) -> tuple[list[Any], list[Any], Mapping[str, Any]]:
    source_jet = pressure_report.get("source_jet", pressure_report)
    comparison_q = source_jet.get("comparison_q")
    if not isinstance(comparison_q, Mapping):
        raise ValueError("Saved pressure comparison-q proof is missing")
    proofs = comparison_q.get("correlated_shared_phi_proof")
    if not _is_sequence(proofs) or pressure_index >= len(proofs):
        raise ValueError("Saved correlated comparison-q proof is missing for this Z")
    proof = proofs[pressure_index]
    if not isinstance(proof, Mapping) or proof.get("passed") is not True:
        raise ValueError("Saved correlated comparison-q proof did not pass")
    if proof.get("shared_Phi0_cancellation_is_algebraic") is not True:
        raise ValueError("The shared-Phi q identity is required")
    orders = proof.get("orders")
    if not _is_sequence(orders) or len(orders) != 3:
        raise ValueError("Saved q proof must retain ordinary orders 0 through 2")
    caps: list[Any] = []
    for order, item in enumerate(orders):
        if not isinstance(item, Mapping) or item.get("ordinary_Z_order") != order:
            raise ValueError("Saved q proof ordinary orders are not contiguous")
        cap = item.get("representable_absolute_error_cap")
        if cap is None:
            raise ValueError("Saved q proof error cap is missing")
        caps.append(cap)
    point = [ctx.mpf(1), ctx.mpf(0), ctx.mpf(0)]
    truth = [
        iv.mpf(1 if order == 0 else 0) + iv.mpf([-1, 1]) * iv.mpf(cap)
        for order, cap in enumerate(caps)
    ]
    return point, truth, proof


def _saved_rho2(ctx: Any, iv: Any, hydro_row: Mapping[str, Any]) -> tuple[Any, Any, Mapping[str, Any]]:
    finite_width = hydro_row.get("finite_width_comparison")
    if not isinstance(finite_width, Mapping):
        raise ValueError("Saved hydro finite-width comparison is missing")
    proofs = finite_width.get("rho_error_proof")
    if not _is_sequence(proofs) or not proofs:
        raise ValueError("Saved rho proof is missing")
    proof = proofs[0]
    if not isinstance(proof, Mapping) or proof.get("passed") is not True:
        raise ValueError("Saved rho proof did not pass")
    cap = proof.get("representable_absolute_error_cap")
    if cap is None:
        raise ValueError("Saved rho proof error cap is missing")
    return ctx.mpf(4), iv.mpf(4) + iv.mpf([0, 1]) * iv.mpf(cap), proof


def _epsilon_from_sample(ctx: Any, iv: Any, sample: Mapping[str, Any]) -> tuple[Any, Any, Any, Any, Mapping[str, Any]]:
    if not isinstance(sample, Mapping):
        raise ValueError("The phase-zero first-collar sample is required")
    if "phase" in sample and str(sample["phase"]) != "0":
        raise ValueError("hydrate_axial_modes requires the saved phase-zero sample")
    scales = sample.get("source_scales", sample)
    if not isinstance(scales, Mapping) or "log_epsilon" not in scales:
        raise ValueError("The saved phase-zero log-epsilon source scale is required")
    log_point, log_truth = decode_pair(ctx, iv, scales["log_epsilon"])
    epsilon_point = ctx.exp(log_point)
    epsilon_truth = iv.exp(log_truth)
    return log_point, log_truth, epsilon_point, epsilon_truth, scales["log_epsilon"]


def _macro_truth(
    ctx: Any,
    iv: Any,
    truth: Sequence[Any],
    weight: Any,
    hcap: Any,
) -> tuple[list[Any], Mapping[str, Any]]:
    result, proof = macro_nonlinear_enclosure(
        ctx,
        iv,
        list(truth),
        weight,
        hcap,
        macro_norm_bound_proved=True,
    )
    return result, proof


def hydrate_axial_modes(
    ctx: Any,
    iv: Any,
    ztext: str,
    pressure_report: Mapping[str, Any],
    hydro_report: Mapping[str, Any],
    swirl_report: Mapping[str, Any],
    bridge: Mapping[str, Any],
    weight: Any,
    hcap: Any,
    *,
    sample: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Hydrate pressure, hydro, and swirl axial modes at one saved ``Z``.

    ``sample`` is the phase-zero first-collar source record, supplied by the
    caller so this adapter never searches or replays an earlier owner.  Its
    exact directed log-epsilon pair is exponentiated separately in the point
    and interval contexts.  The return contains live point/interval rows for
    the weighted primitive and JSON-safe evidence for every new enclosure.
    """

    if sample is None:
        raise ValueError("hydrate_axial_modes requires sample=phase-zero first-collar record")
    weight_iv = _exact_interval(ctx, iv, weight)
    hcap_iv = _exact_interval(ctx, iv, hcap)
    logepsilon, logepsilon_truth, epsilon, epsilon_truth, epsilon_source = _epsilon_from_sample(
        ctx, iv, sample
    )

    pressure_rows = _rows(pressure_report, ("source_jet", "rows"), ("rows",))
    hydro_rows = _rows(hydro_report, ("hydro_source_rows",), ("rows",))
    swirl_rows = _rows(swirl_report, ("swirl_source_rows",), ("rows",))
    pressure_index, pressure_row = _find_row(pressure_rows, ztext)
    _, hydro_row = _find_row(hydro_rows, ztext)
    _, swirl_row = _find_row(swirl_rows, ztext)

    pressure_point, pressure_truth = _decode_jet(
        ctx, iv, pressure_row.get("source_rows"), name="pressure B0 source"
    )
    qpoint, qtruth, qproof = _saved_q(ctx, iv, pressure_report, pressure_index)
    rho_point, rho_truth, rhoproof = _saved_rho2(ctx, iv, hydro_row)

    target_entries = _factored_entries(
        ctx, iv, hydro_row, "epsilon_factored_target_mode", {0, 1, 2, 3}
    )
    difference_entries = _factored_entries(
        ctx, iv, hydro_row, "epsilon_factored_moment_difference_mode", {1, 2, 3}
    )
    denominator_entries = _factored_entries(
        ctx, iv, hydro_row, "epsilon_factored_denominator", {0, 1}
    )

    denominator_point = _evaluate_sectors(
        ctx, denominator_entries, epsilon, truth=False, exponent_shift=0
    )
    denominator_truth = _evaluate_sectors(
        iv, denominator_entries, epsilon_truth, truth=True, exponent_shift=0
    )

    target_numerator_point = _evaluate_sectors(
        ctx, target_entries, epsilon, truth=False, factor=100, exponent_shift=-1
    )
    target_numerator_truth = _evaluate_sectors(
        iv, target_entries, epsilon_truth, truth=True, factor=100, exponent_shift=-1
    )
    target_point = ordinary_product(
        qpoint, ordinary_quotient(_negative(target_numerator_point), denominator_point)
    )
    target_truth = ordinary_product(
        qtruth, ordinary_quotient(_negative(target_numerator_truth), denominator_truth)
    )
    target_truth, target_macro_proof = _macro_truth(
        ctx, iv, target_truth, weight_iv, hcap_iv
    )

    difference_numerator_point = _evaluate_sectors(
        ctx,
        difference_entries,
        epsilon,
        truth=False,
        factor=rho_point,
        exponent_shift=0,
    )
    difference_numerator_truth = _evaluate_sectors(
        iv,
        difference_entries,
        epsilon_truth,
        truth=True,
        factor=rho_truth,
        exponent_shift=0,
    )
    difference_point = ordinary_product(
        qpoint,
        ordinary_quotient(_negative(difference_numerator_point), denominator_point),
    )
    difference_truth = ordinary_product(
        qtruth,
        ordinary_quotient(_negative(difference_numerator_truth), denominator_truth),
    )
    difference_truth, difference_macro_proof = _macro_truth(
        ctx, iv, difference_truth, weight_iv, hcap_iv
    )

    modes: list[dict[str, Any]] = [
        {
            "sector": "pressure",
            "mode_j": 0,
            "mode": "B0",
            "rate": 2,
            "point": pressure_point,
            "truth": pressure_truth,
            "serialized_rows": [
                encode_pair(ctx, iv, point, truth)
                for point, truth in zip(pressure_point, pressure_truth)
            ],
            "identity": "saved pressure B0 source row, including its macro enclosure",
        },
        {
            "sector": "hydro",
            "mode_j": 0,
            "mode": "B0_target",
            "rate": 2,
            "point": target_point,
            "truth": target_truth,
            "serialized_rows": [
                encode_pair(ctx, iv, point, truth)
                for point, truth in zip(target_point, target_truth)
            ],
            "identity": "-q*sum_n 100*epsilon^(n-1)*Ntarget_n / D(epsilon), with exp(ell) enclosed",
            "macro_proof": _serial(target_macro_proof),
        },
        {
            "sector": "hydro",
            "mode_j": 1,
            "mode": "B1_moment_difference",
            "rate": 1,
            "point": difference_point,
            "truth": difference_truth,
            "serialized_rows": [
                encode_pair(ctx, iv, point, truth)
                for point, truth in zip(difference_point, difference_truth)
            ],
            "identity": "-q*sum_n rho2*epsilon^n*Ndifference_n / D(epsilon), with exp(ell) enclosed",
            "macro_proof": _serial(difference_macro_proof),
        },
    ]

    pieces = swirl_row.get("mode_delta_pieces")
    if not _is_sequence(pieces):
        raise ValueError("Saved swirl mode-delta pieces are missing")
    grouped: dict[int, dict[int, tuple[list[Any], list[Any]]]] = {}
    for piece in pieces:
        if not isinstance(piece, Mapping):
            raise ValueError("Malformed saved swirl mode-delta piece")
        mode_j = piece.get("mode_j")
        delta_degree = piece.get("delta_degree")
        if (
            isinstance(mode_j, bool)
            or not isinstance(mode_j, int)
            or mode_j not in (0, 1, 2)
            or isinstance(delta_degree, bool)
            or not isinstance(delta_degree, int)
            or delta_degree not in (0, 1)
        ):
            raise ValueError("Saved swirl mode-delta indices are outside the original sectors")
        mode_group = grouped.setdefault(mode_j, {})
        if delta_degree in mode_group:
            raise ValueError("Duplicate saved swirl mode-delta sector")
        mode_group[delta_degree] = _decode_jet(
            ctx,
            iv,
            piece.get("relative_coefficient_ordinary_rows"),
            name=f"swirl relative coefficient ({mode_j},{delta_degree})",
        )
    if set(grouped) != {0, 1, 2} or any(set(group) != {0, 1} for group in grouped.values()):
        raise ValueError("Saved swirl mode-delta pieces omitted an original radial sector")

    scale_packet = bridge.get("swirl_scale_log_enclosure")
    if scale_packet is None:
        raise ValueError("The bridge swirl log-amplitude enclosure is required")
    swirl_log_scale = _exact_interval(ctx, iv, scale_packet)
    relative_swirl_modes: list[dict[str, Any]] = []
    swirl_evidence: list[dict[str, Any]] = []
    for mode_j in (0, 1, 2):
        point_rows: list[Any] = []
        truth_rows: list[Any] = []
        for delta_degree in (0, 1):
            point_piece, truth_piece = grouped[mode_j][delta_degree]
            if not point_rows:
                point_rows = list(point_piece)
                truth_rows = list(truth_piece)
            else:
                point_rows = [left + right for left, right in zip(point_rows, point_piece)]
                truth_rows = [left + right for left, right in zip(truth_rows, truth_piece)]
        relative_truth, relative_macro_proof = _macro_truth(
            ctx, iv, truth_rows, weight_iv, hcap_iv
        )
        absolute_point, absolute_truth, atom_proof = scaled_atom_enclosure(
            ctx,
            iv,
            relative_truth,
            swirl_log_scale,
            source_proved=True,
        )
        rate = 3 - mode_j
        relative_swirl_modes.append(
            {
                "mode_j": mode_j,
                "mode": f"relative_B{mode_j}",
                "rate": rate,
                "point": point_rows,
                "truth": relative_truth,
                "serialized_rows": [
                    encode_pair(ctx, iv, point, truth)
                    for point, truth in zip(point_rows, relative_truth)
                ],
                "identity": "grouped signed relative coefficient sectors; actual positive F0^2 remains factored",
            }
        )
        modes.append(
            {
                "sector": "swirl",
                "mode_j": mode_j,
                "mode": f"B{mode_j}",
                "rate": rate,
                "point": absolute_point,
                "truth": absolute_truth,
                "serialized_rows": [
                    encode_pair(ctx, iv, point, truth)
                    for point, truth in zip(absolute_point, absolute_truth)
                ],
                "identity": "scaled positive F0^2 atom applied to the signed relative mode; zero is a controlled proxy",
                "relative_point": point_rows,
                "relative_truth": relative_truth,
                "macro_proof": _serial(relative_macro_proof),
                "scaled_atom_proof": _serial(atom_proof),
            }
        )
        swirl_evidence.append(
            {
                "mode_j": mode_j,
                "rate": rate,
                "relative_rows": [
                    encode_pair(ctx, iv, point, truth)
                    for point, truth in zip(point_rows, relative_truth)
                ],
                "absolute_rows": [
                    encode_pair(ctx, iv, point, truth)
                    for point, truth in zip(absolute_point, absolute_truth)
                ],
                "macro_proof": _serial(relative_macro_proof),
                "scaled_atom_proof": _serial(atom_proof),
                "identity": "F0^2 is positive and retained as a source atom; no logarithmic amplitude is exponentiated",
            }
        )

    saved_atom = swirl_row.get("absolute_source_atom")
    if not isinstance(saved_atom, Mapping) or saved_atom.get("actual_amplitude_positive") is not True:
        raise ValueError("The saved positive swirl source atom is missing")

    evidence = {
        "passed": True,
        "Z": ztext,
        "ordinary_Z_orders": [0, 1, 2],
        "source_identity": {
            "pressure": "saved pressure B0 ordinary source already macro-enclosed",
            "hydro_target": "-q*sum_n 100*epsilon^(n-1)*Ntarget_n / D(epsilon)",
            "hydro_moment_difference": "-q*sum_n rho2*epsilon^n*Ndifference_n / D(epsilon)",
            "swirl": "group saved relative coefficient delta sectors by mode_j and retain F0^2 as a positive atom",
        },
        "epsilon_source": _serial(epsilon_source),
        "epsilon_point": _serial(epsilon),
        "epsilon_truth": _serial(epsilon_truth),
        "log_epsilon_point": _serial(logepsilon),
        "log_epsilon_truth": _serial(logepsilon_truth),
        "q_proof": _serial(qproof),
        "rho2_point": _serial(rho_point),
        "rho2_truth": _serial(rho_truth),
        "rho2_proof": _serial(rhoproof),
        "denominator_sectors": _serial(
            {
                str(degree): {
                    "point": point_rows,
                    "truth": truth_rows,
                }
                for degree, (point_rows, truth_rows) in denominator_entries.items()
            }
        ),
        "hydro_mode_proofs": {
            "target": _serial(target_macro_proof),
            "moment_difference": _serial(difference_macro_proof),
        },
        "swirl_mode_proofs": swirl_evidence,
        "swirl_scale_log_enclosure": _serial(swirl_log_scale),
        "all_original_sectors_retained": True,
        "source_width_selected": False,
        "positive_swirl_amplitude_materialized": False,
        "absolute_zero_is_controlled_proxy_only": True,
    }
    return {
        "Z": ztext,
        "ordinary_Z_orders": [0, 1, 2],
        "modes": modes,
        "relative_swirl_modes": relative_swirl_modes,
        "absolute_swirl_atom": saved_atom,
        "evidence": evidence,
    }


__all__ = ["hydrate_axial_modes"]
