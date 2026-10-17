"""Pure factored leaves for the original weighted viscous primitive.

The weighted source identity used here is

``J = sum_j B_j K_(q-j+1) = M0*sum_j B_j + h*R``

and hence

``V = R100*epsilon*h**2*(M0*C_delta(sum_j B_j) + h*C_delta(R))``.

This module only performs the ordinary source-row algebra.  It does not load
an accepted packet, integrate a source, materialize a width from a logarithm,
select an error cap, or replace a zero source by a synthetic row.  ``h`` is
kept as a symbolic name in the returned factored records.  A caller that has
an independently evaluated finite-width source can evaluate the two returned
width powers itself.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from mpmath.libmp import fzero, mpf_cmp


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _is_scalar(value: Any) -> bool:
    return not _is_sequence(value) and not isinstance(value, Mapping)


def _row(value: Any, *, name: str) -> list[Any]:
    """Validate one full ordinary row through ``Z2``.

    A full source row is required even though ``C_delta`` returns only its
    ordinary ``Z0/Z1`` coefficients.  In particular, the third source row is
    never inferred or zero padded.
    """

    if not _is_sequence(value) or len(value) != 3:
        raise ValueError(f"{name} must contain ordinary Z0/Z1/Z2 rows")
    if any(item is None for item in value):
        raise ValueError(f"{name} contains a missing ordinary derivative row")
    if any(not _is_scalar(item) for item in value):
        raise ValueError(f"{name} contains a nested ordinary row")
    return list(value)


def _sum_rows(rows: Sequence[Sequence[Any]], *, name: str) -> list[Any]:
    if not rows:
        raise ValueError(f"{name} must contain at least one source row")
    normalized = [_row(row, name=f"{name}[{index}]") for index, row in enumerate(rows)]
    return [
        sum((row[order] for row in normalized), normalized[0][order] * 0)
        for order in range(3)
    ]


def _extract_candidate(value: Any, *, name: str, aliases: tuple[str, ...]) -> Any:
    """Find source rows in a small raw-packet mapping.

    The accepted raw 46/raw 47 packet names have changed between reports.  The
    parser accepts only structural aliases and never decodes a saved positive
    enclosure or chooses a midpoint as a signed source value.
    """

    if isinstance(value, Mapping):
        for key in aliases:
            if key in value:
                return value[key]
        for key in (
            "raw46",
            "raw47",
            "mainpacket",
            "signedquery",
            "packet",
            "query",
            "source",
            "source_packet",
        ):
            if key in value:
                return _extract_candidate(value[key], name=name, aliases=aliases)
        raise ValueError(f"{name} is missing full ordinary source rows")
    return value


def _source_rows(value: Any, *, name: str, aliases: tuple[str, ...]) -> tuple[list[Any], list[list[Any]]]:
    """Return ``(sum_j row_j, row_j list)`` from a raw packet.

    A direct three-scalar row is one already-summed signed source.  A list of
    ordinary rows is interpreted as the explicit modal ``B_j`` decomposition
    and is summed componentwise.  Both forms retain every input row in the
    returned leaf for provenance.
    """

    candidate = _extract_candidate(value, name=name, aliases=aliases)
    if isinstance(candidate, Mapping):
        # One extra structural wrapper is useful for records such as
        # {"ordinary_rows": {"rows": [...]}}.
        candidate = _extract_candidate(candidate, name=name, aliases=aliases)
    if not _is_sequence(candidate):
        raise ValueError(f"{name} must contain full ordinary source rows")
    if len(candidate) == 3 and all(_is_scalar(item) for item in candidate):
        row = _row(candidate, name=name)
        return row, [row]
    if not candidate:
        raise ValueError(f"{name} must contain at least one full ordinary source row")
    rows: list[list[Any]] = []
    for index, item in enumerate(candidate):
        if isinstance(item, Mapping):
            item = _extract_candidate(
                item,
                name=f"{name}[{index}]",
                aliases=("ordinary_rows", "rows", "source_rows", "coefficients"),
            )
        rows.append(_row(item, name=f"{name}[{index}]"))
    return _sum_rows(rows, name=name), rows


def _scalar(value: Any, *, name: str) -> Any:
    if value is None:
        raise ValueError(f"{name} is required")
    if not _is_scalar(value):
        raise ValueError(f"{name} must be a Z-independent scalar")
    return value


def _scale(row: Sequence[Any], scalar: Any) -> list[Any]:
    if not _is_sequence(row) or not row:
        raise ValueError("ordinary coefficient row is required")
    if any(not _is_scalar(value) for value in row):
        raise ValueError("ordinary coefficient row contains a nested value")
    return [value * scalar for value in row]


def ordinary_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply two full ordinary rows through ``Z2``.

    This tiny helper keeps the ordinary-derivative convention explicit for
    callers assembling a modal ``B_j`` product before passing its sum to the
    leaf builder.  The ``2 * left[1] * right[1]`` term is required because the
    rows store derivatives, rather than Taylor coefficients.
    """

    a = _row(left, name="left ordinary row")
    b = _row(right, name="right ordinary row")
    return [
        a[0] * b[0],
        a[1] * b[0] + a[0] * b[1],
        a[2] * b[0] + 2 * a[1] * b[1] + a[0] * b[2],
    ]


def _interval_endpoints(value: Any) -> tuple[Any, Any] | None:
    if hasattr(value, "_mpi_"):
        return value._mpi_
    return None


def _check_nonzero_denominator(value: Any) -> None:
    endpoints = _interval_endpoints(value)
    if endpoints is not None:
        lower, upper = endpoints
        # ``ivmpf._mpi_`` stores raw four-tuples, not ordinary mpf objects.
        # Use libmp comparisons so this check remains exact and directed.
        if mpf_cmp(upper, lower) < 0:
            raise ValueError("C_delta denominator has reversed directed endpoints")
        if mpf_cmp(lower, fzero) <= 0 <= mpf_cmp(upper, fzero):
            raise ValueError("C_delta denominator enclosure crosses zero")
    elif value == 0:
        raise ValueError("C_delta denominator is zero")


def c_delta(z: Any, delta: Any, source_rows: Sequence[Any]) -> list[Any]:
    """Return ordinary ``C_delta(B)`` through ``Z1`` from ``B/B_Z/B_ZZ``.

    The input rows are ordinary derivatives (not Taylor coefficients), so the
    product derivative in the ``Z1`` slot carries its usual factor two only
    where it arises from an ordinary product.  The formula is

    ``C_delta(B) = ((1-delta)*Z*B + (1-Z**2)*B_Z)/(1-delta*Z**2)``.
    """

    b, bz, bzz = _row(source_rows, name="source_rows")
    z = _scalar(z, name="z")
    delta = _scalar(delta, name="delta")
    denominator = 1 - delta * z * z
    _check_nonzero_denominator(denominator)
    numerator = (1 - delta) * z * b + (1 - z * z) * bz
    numerator_z = (1 - delta) * b - (1 + delta) * z * bz + (1 - z * z) * bzz
    return [
        numerator / denominator,
        numerator_z / denominator
        + 2 * delta * z * numerator / (denominator * denominator),
    ]


# Clear spelling for callers that prefer the displayed mathematical name.
apply_c_delta = c_delta
C_delta = c_delta


def _source_guard(source_proved: Any, shared_source_proved: Any) -> None:
    if source_proved is not True:
        raise ValueError("an explicit same-source weighted-primitive proof is required")
    if shared_source_proved is not None and shared_source_proved is not True:
        raise ValueError("the exact shared source proof is required")


def _mass_guard(mass: Any, mass_Z_independent: Any) -> Any:
    if mass_Z_independent is not True:
        raise ValueError("the supplied M0 mass must be proved Z-independent")
    mass = _scalar(mass, name="mass")
    endpoints = _interval_endpoints(mass)
    if endpoints is not None:
        lower, upper = endpoints
        if mpf_cmp(upper, lower) < 0 or mpf_cmp(lower, fzero) < 0:
            raise ValueError("mass must be a nonnegative directed scalar")
    elif mass < 0:
        raise ValueError("mass must be nonnegative")
    return mass


def _width_symbol(value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("width_symbol must be a nonempty symbolic name")
    return value.strip()


def _leaf(
    *,
    name: str,
    raw_number: int,
    source_rows: list[Any],
    source_modes: list[list[Any]],
    z: Any,
    delta: Any,
    mass: Any,
    common_scale: Any,
    width_symbol: str,
) -> dict[str, Any]:
    c_rows = c_delta(z, delta, source_rows)
    # raw 46 is M0*C(B), while raw 47 is C(R), the coefficient of one
    # additional width.  Keep both the unscaled and physical coefficient rows.
    unscaled = _scale(c_rows, mass) if raw_number == 46 else list(c_rows)
    physical = _scale(unscaled, common_scale)
    power = 2 if raw_number == 46 else 3
    if raw_number == 46:
        formula = f"{common_scale!r}*{width_symbol}^2*mass*C_delta(sum_j B_j)"
        role = "old_raw46_mainpacket"
    else:
        formula = f"{common_scale!r}*{width_symbol}^3*C_delta(R)"
        role = "new_raw47_signedquery"
    return {
        "kind": "factored_signed_leaf",
        "role": role,
        "raw_number": raw_number,
        "signed": True,
        "ordinary_source_rows": list(source_rows),
        "ordinary_source_modes": [list(row) for row in source_modes],
        "source_ordinary_orders": [0, 1, 2],
        "ordinary_C0_C1": list(c_rows),
        "coefficient_rows": list(physical),
        "unscaled_coefficient_rows": list(unscaled),
        "ordinary_coefficient_orders": [0, 1],
        "width_factor": {"symbol": width_symbol, "power": power},
        "width_power": power,
        "physical_scale": common_scale,
        "physical_scale_factors": ["R100", "epsilon"],
        "mass_factor": mass if raw_number == 46 else None,
        "formula": formula,
        "source_proved": True,
        "mass_Z_independent": True,
        "proof": {
            "same_source_identity": True,
            "source_proved": True,
            "mass_Z_independent": True,
            "raw46_is_M0_C_of_sum_B": True,
            "raw47_is_C_of_signed_R": True,
            "full_source_row_through_Z2_required": True,
        },
    }


def _combined(
    raw46: Mapping[str, Any],
    raw47: Mapping[str, Any],
    *,
    common_scale: Any,
    width_symbol: str,
) -> dict[str, Any]:
    # Store the polynomial in h as separate terms.  This is intentionally a
    # symbolic representation: no numerical width is multiplied into a row.
    main = list(raw46["unscaled_coefficient_rows"])
    remainder = list(raw47["unscaled_coefficient_rows"])
    return {
        "kind": "combined_factored_jet",
        "signed": True,
        "width_symbol": width_symbol,
        "outer_width_power": 2,
        "physical_scale": common_scale,
        "physical_scale_factors": ["R100", "epsilon"],
        "ordinary_orders": [0, 1],
        "terms": [
            {
                "width_power_inside_factor": 0,
                "coefficient_rows": main,
                "source": "raw46",
            },
            {
                "width_power_inside_factor": 1,
                "coefficient_rows": remainder,
                "source": "raw47",
            },
        ],
        "ordinary_C0_C1": {
            "constant": main,
            "linear_in_width": remainder,
        },
        "formula": (
            f"{common_scale!r}*{width_symbol}^2*"
            "(M0*C_delta(sum_j B_j)+"
            f"{width_symbol}*C_delta(R))"
        ),
        "source_proved": True,
        "mass_Z_independent": True,
        "proof": {
            "same_source_identity": True,
            "source_proved": True,
            "mass_Z_independent": True,
            "full_source_row_through_Z2_required": True,
        },
        "width_selected": False,
        "width_materialized": False,
    }


def build_factored_leaf_records(
    raw46: Any = None,
    raw47: Any = None,
    z: Any = None,
    delta: Any = None,
    mass: Any = None,
    *,
    raw46_mainpacket: Any = None,
    raw47_signedquery: Any = None,
    R100: Any = 100,
    epsilon: Any = 1,
    width_symbol: str = "h",
    source_proved: Any = None,
    shared_source_proved: Any = None,
    mass_Z_independent: Any = None,
    include_combined_factored_jet: bool = True,
    optionalcombinedfactoredjet: Any = None,
) -> dict[str, Any]:
    """Build raw 46/raw 47 signed leaves and an optional combined jet.

    ``raw46`` and ``raw47`` may be full ``[Z0,Z1,Z2]`` rows, explicit lists of
    modal rows, or mappings containing those rows under ``ordinary_rows``,
    ``source_rows``, ``B_rows``/``R_rows`` and the corresponding packet names.
    The returned coefficients have ordinary order ``Z0/Z1`` while retaining
    the full source rows for the proof that generated them.
    """

    if raw46 is None:
        raw46 = raw46_mainpacket
    elif raw46_mainpacket is not None:
        raise ValueError("supply only one raw46/mainpacket input")
    if raw47 is None:
        raw47 = raw47_signedquery
    elif raw47_signedquery is not None:
        raise ValueError("supply only one raw47/signedquery input")
    if raw46 is None or raw47 is None:
        raise ValueError("both raw46 mainpacket and raw47 signedquery are required")
    if z is None or delta is None:
        raise ValueError("z and delta are required")
    _source_guard(source_proved, shared_source_proved)
    mass = _mass_guard(mass, mass_Z_independent)
    z = _scalar(z, name="z")
    delta = _scalar(delta, name="delta")
    R100 = _scalar(R100, name="R100")
    epsilon = _scalar(epsilon, name="epsilon")
    width_symbol = _width_symbol(width_symbol)
    common_scale = R100 * epsilon
    b_sum, b_modes = _source_rows(
        raw46,
        name="raw46 mainpacket",
        aliases=(
            "ordinary_rows",
            "rows",
            "source_rows",
            "B_rows",
            "raw46_rows",
            "main_rows",
            "mainpacket",
        ),
    )
    r_sum, r_modes = _source_rows(
        raw47,
        name="raw47 signedquery",
        aliases=(
            "ordinary_rows",
            "rows",
            "source_rows",
            "R_rows",
            "raw47_rows",
            "remainder_rows",
            "signedquery",
        ),
    )
    raw46_leaf = _leaf(
        name="raw46",
        raw_number=46,
        source_rows=b_sum,
        source_modes=b_modes,
        z=z,
        delta=delta,
        mass=mass,
        common_scale=common_scale,
        width_symbol=width_symbol,
    )
    raw47_leaf = _leaf(
        name="raw47",
        raw_number=47,
        source_rows=r_sum,
        source_modes=r_modes,
        z=z,
        delta=delta,
        mass=mass,
        common_scale=common_scale,
        width_symbol=width_symbol,
    )
    if optionalcombinedfactoredjet is not None:
        include_combined_factored_jet = bool(optionalcombinedfactoredjet)
    result: dict[str, Any] = {
        "passed": True,
        "raw46": raw46_leaf,
        "raw47": raw47_leaf,
        "raw46_leaf": raw46_leaf,
        "raw47_leaf": raw47_leaf,
        "raw46_mainpacket": raw46_leaf,
        "raw47_signedquery": raw47_leaf,
        "source_proved": True,
        # ``source_proved`` is the required exact same-source certificate;
        # the optional spelling is accepted only as an additional guard.
        "shared_source_proved": True,
        "mass_Z_independent": True,
        "width_symbol": width_symbol,
        "width_selected": False,
        "width_materialized": False,
        "source_rows_preserved_through_Z2": True,
        "proof": {
            "same_source_identity": True,
            "source_proved": True,
            "mass_Z_independent": True,
            "full_source_row_through_Z2_required": True,
            "width_is_symbolic": True,
        },
        "old_raw46_retained": True,
        "new_raw47_retained": True,
        "combined_factored_jet": None,
        "optional_combined_factored_jet": None,
        "optionalcombinedfactoredjet": None,
    }
    if include_combined_factored_jet:
        combined = _combined(
            raw46_leaf,
            raw47_leaf,
            common_scale=common_scale,
            width_symbol=width_symbol,
        )
        result["combined_factored_jet"] = combined
        result["optional_combined_factored_jet"] = combined
        result["optionalcombinedfactoredjet"] = combined
    return result


def combine_coefficients(
    z: Any = None,
    delta: Any = None,
    raw46: Any = None,
    raw47: Any = None,
    mass: Any = None,
    *,
    raw46_mainpacket: Any = None,
    raw47_signedquery: Any = None,
    R100: Any = 100,
    epsilon: Any = 1,
    width_symbol: str = "h",
    source_proved: Any = None,
    shared_source_proved: Any = None,
    mass_Z_independent: Any = None,
    include_combined_factored_jet: bool = True,
    optionalcombinedfactoredjet: Any = None,
) -> dict[str, Any]:
    """Alias-shaped entry point for combining raw 46 and signed raw 47."""

    return build_factored_leaf_records(
        raw46,
        raw47,
        z,
        delta,
        mass,
        raw46_mainpacket=raw46_mainpacket,
        raw47_signedquery=raw47_signedquery,
        R100=R100,
        epsilon=epsilon,
        width_symbol=width_symbol,
        source_proved=source_proved,
        shared_source_proved=shared_source_proved,
        mass_Z_independent=mass_Z_independent,
        include_combined_factored_jet=include_combined_factored_jet,
        optionalcombinedfactoredjet=optionalcombinedfactoredjet,
    )


combinecoefficients = combine_coefficients
combine_factored_coefficients = combine_coefficients
factored_leaf_records = build_factored_leaf_records


def combine_factored_rows(
    raw46: Any,
    raw47: Any,
    z: Any,
    delta: Any,
    mass: Any,
    **kwargs: Any,
) -> dict[str, Any]:
    """Input-order convenience wrapper for callers starting from the leaves."""

    return combine_coefficients(z, delta, raw46, raw47, mass, **kwargs)


__all__ = [
    "c_delta",
    "C_delta",
    "apply_c_delta",
    "ordinary_product",
    "build_factored_leaf_records",
    "factored_leaf_records",
    "combine_coefficients",
    "combinecoefficients",
    "combine_factored_coefficients",
    "combine_factored_rows",
]
