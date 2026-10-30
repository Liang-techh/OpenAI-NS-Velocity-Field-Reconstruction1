"""Scaled signed finite-``N`` density arithmetic.

This module is the normalized bridge used by the finite-``N`` source ports.
It accepts ordinary derivatives of the normalized baseline fields ``V`` and
``q`` and of a same-source correction ``dV``/``dm``.  The correction pressure
is recovered from the source constraint before the density is evaluated::

    L = 1 - delta*Z**2,       D = 1 - Z**2
    dq = (2*Z*dV - Z*(1-delta)*dm - D*dm_Z)/L
    dq_y = (Z*(1+delta)*dV + 2*Z*dV_y - D*dV_Z)/L - dq

The physical fields have a formal common factor ``Pstar``.  Consequently
the signed density difference has the exact sector form

    Delta v = Pstar*v1 + Pstar**2*v2,
    Delta F = Pstar*F1 + Pstar**2*F2.

Only the normalized sectors are returned.  ``Pstar`` and any finite-``N``
power are intentionally not represented in this module.  All arithmetic is
delegated to :mod:`original_finite_N_density_kernel`, so exact ``Fraction``
values, SymPy expressions, scalar-expression graphs, point mpmath values,
and directed ``mpmath.iv`` values use the same ordinary-derivative rules.

Rows are ordinary derivative values, never Taylor coefficients.  A row
``[f, f_Z, f_ZZ]`` therefore stores the derivatives themselves; no factorial
conversion or automatic zero padding is performed.  ``dm`` needs four rows
because the C2 pressure recovery differentiates ``dm_Z`` twice.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from numbers import Number
from typing import Any, Optional

from .original_finite_N_density_kernel import (
    OrdinaryJet,
    PYTHON_OPERATIONS,
    ScalarOperations,
    endpoint_difference_terms,
    full_endpoint_density_difference,
    full_volume_density_difference,
    ordinary_jet,
    physical_fields,
    volume_difference_terms,
)


class ScaledDensityJetError(ValueError):
    """Raised when a normalized ordinary derivative row is incomplete."""


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _contains_nested_ordinary_jet(value: Any) -> bool:
    """Return whether a row container contains an ``OrdinaryJet`` component.

    A row is a container of scalar components.  An ``OrdinaryJet`` inside a
    mapping or sequence is therefore a malformed nested row, not a source
    adapter to borrow while selecting the outer row's adapter.
    """

    if isinstance(value, Mapping):
        return any(_contains_nested_ordinary_jet(item) or isinstance(item, OrdinaryJet)
                   for item in value.values())
    if _is_sequence(value):
        return any(_contains_nested_ordinary_jet(item) or isinstance(item, OrdinaryJet)
                   for item in value)
    return False


def _reject_nested_ordinary_jet(value: Any, name: str) -> None:
    if not isinstance(value, OrdinaryJet) and _contains_nested_ordinary_jet(value):
        raise ScaledDensityJetError(
            f"{name} contains a nested OrdinaryJet; supply one outer row jet"
        )


def _row_operations(value: Any, operations: Optional[ScalarOperations]) -> Optional[ScalarOperations]:
    """Validate and preserve an outer row jet's explicit adapter identity."""

    if isinstance(value, OrdinaryJet):
        if value.operations is not None and operations is not None and value.operations is not operations:
            raise TypeError(
                "ordinary row uses a different scalar operation adapter; "
                "explicit adapters cannot be relabeled"
            )
        return value.operations if value.operations is not None else operations
    return operations


def _literal_or_source(value: Any, operations: Optional[ScalarOperations]) -> Any:
    """Lift numeric literals while preserving opaque source scalar identity."""

    if operations is not None and isinstance(value, Number):
        return operations.constant(value)
    return value


def _zero_scalar(operations: Optional[ScalarOperations]) -> Any:
    return 0 if operations is None else operations.constant(0)


def _one_scalar(operations: Optional[ScalarOperations]) -> Any:
    return 1 if operations is None else operations.constant(1)


def _row(value: Any, order: int, name: str, *, operations: Optional[ScalarOperations] = None) -> list[Any]:
    """Validate one explicit ordinary derivative row.

    ``order`` is the highest derivative required by the formula.  Exact row
    length is enforced so a Taylor prefix, an omitted derivative, or a hidden
    zero-padding convention cannot enter the source algebra accidentally.
    ``OrdinaryJet`` and value/dz/dzz mappings are accepted only where they
    contain all requested components; a fourth derivative still must be
    supplied as a sequence for ``dm``.
    """

    required = order + 1
    if value is None:
        raise ScaledDensityJetError(f"{name} is required; missing ordinary derivatives")
    _reject_nested_ordinary_jet(value, name)
    _row_operations(value, operations)
    if isinstance(value, OrdinaryJet):
        if order > 2:
            raise ScaledDensityJetError(
                f"{name} needs ordinary Z0..Z{order}; an OrdinaryJet has no explicit Z{order}"
            )
        entries = [value.value, value.dz, value.dzz]
        if any(item is None for item in entries[:required]):
            raise ScaledDensityJetError(f"{name} contains a missing ordinary derivative")
        return entries[:required]
    if isinstance(value, Mapping):
        # Mappings are deliberately narrow.  A generic coefficient mapping
        # is not treated as a Taylor row, and cannot supply Z3 implicitly.
        keys = ("value", "val", "f")
        raw = next((value[key] for key in keys if key in value), None)
        if raw is None:
            raise ScaledDensityJetError(f"{name} mapping has no ordinary value component")
        dz = next((value[key] for key in ("dz", "z", "Z", "first", "d_dZ") if key in value), None)
        dzz = next((value[key] for key in ("dzz", "zz", "ZZ", "second", "d2_dZ2") if key in value), None)
        if order > 2:
            raise ScaledDensityJetError(
                f"{name} needs an explicit ordinary Z3 row; mappings cannot hide Z3"
            )
        if order >= 1 and dz is None:
            raise ScaledDensityJetError(f"{name} mapping is missing ordinary Z1")
        if order >= 2 and dzz is None:
            raise ScaledDensityJetError(f"{name} mapping is missing ordinary Z2")
        entries = [_literal_or_source(item, operations) for item in (raw, dz, dzz)]
        if any(item is None for item in entries[:required]):
            raise ScaledDensityJetError(f"{name} contains a missing ordinary derivative")
        return entries[:required]
    if not _is_sequence(value):
        raise ScaledDensityJetError(
            f"{name} must be an ordinary derivative sequence of length {required}"
        )
    if len(value) != required:
        raise ScaledDensityJetError(
            f"{name} must contain exactly ordinary Z0..Z{order}; no zero padding or Taylor coefficients"
        )
    result = [_literal_or_source(item, operations) for item in value]
    if any(item is None for item in result):
        raise ScaledDensityJetError(f"{name} contains a missing ordinary derivative")
    return result


def explicit_shift_jet(
    rows: Sequence[Any],
    *,
    name: str = "ordinary_rows",
    operations: Optional[ScalarOperations] = None,
) -> OrdinaryJet:
    """Return the explicit ordinary derivative jet ``[f_Z, f_ZZ, f_ZZZ]``.

    The input must contain exactly four ordinary rows.  This helper is kept
    separate to make the source requirement visible at the pressure recovery
    call site; it is not a Taylor shift and does not multiply by factorials.
    """

    row_operations = _row_operations(rows, operations)
    checked = _row(rows, 3, name, operations=row_operations)
    return ordinary_jet(checked[1], checked[2], checked[3], operations=row_operations)


# The compact spelling is useful to callers mirroring the source notation.
explicitshiftjet = explicit_shift_jet


def _as_jet(
    rows: Any,
    order: int,
    name: str,
    *,
    operations: Optional[ScalarOperations] = None,
) -> OrdinaryJet:
    row_operations = _row_operations(rows, operations)
    checked = _row(rows, order, name, operations=row_operations)
    if order == 0:
        return ordinary_jet(checked[0], operations=row_operations)
    if order == 1:
        return ordinary_jet(checked[0], checked[1], operations=row_operations)
    return ordinary_jet(checked[0], checked[1], checked[2], operations=row_operations)


def _coordinate_jet(
    value: Any,
    *,
    operations: Optional[ScalarOperations] = None,
) -> OrdinaryJet:
    """Coerce the independent ``Z`` coordinate to an explicit C2 jet.

    A scalar is accepted as shorthand for ``[Z, 1, 0]``.  A supplied row or
    ``OrdinaryJet`` is checked as an explicit ordinary coordinate jet.
    """

    if isinstance(value, OrdinaryJet):
        if value.dzz is None:
            raise ScaledDensityJetError("Z_jet must provide ordinary Z0, Z1 and Z2")
        _row_operations(value, operations)
        return value if operations is None else ordinary_jet(
            value.value, value.dz, value.dzz, operations=operations
        )
    if _is_sequence(value) or isinstance(value, Mapping):
        return _as_jet(value, 2, "Z_jet", operations=operations)
    # Z is the independent variable.  Its ordinary derivatives are therefore
    # fixed explicitly even when a caller supplies its scalar value only.
    return ordinary_jet(
        _literal_or_source(value, operations),
        _one_scalar(operations),
        _zero_scalar(operations),
        operations=operations,
    )


def _constant_jet(value: Any, operations: Optional[ScalarOperations]) -> OrdinaryJet:
    return ordinary_jet(
        _literal_or_source(value, operations),
        _zero_scalar(operations),
        _zero_scalar(operations),
        operations=operations,
    )


def _contains_interval(value: Any) -> bool:
    """Recognize an ``mpmath.iv`` scalar without importing it eagerly."""

    if isinstance(value, OrdinaryJet):
        return any(_contains_interval(item) for item in (value.value, value.dz, value.dzz))
    if isinstance(value, Mapping):
        return any(_contains_interval(item) for item in value.values())
    if _is_sequence(value):
        return any(_contains_interval(item) for item in value)
    cls = type(value)
    return cls.__name__ == "ivmpf" or cls.__module__ == "mpmath.ctx_iv"


def _interval_operations() -> ScalarOperations:
    """Use interval constants instead of Fraction constants for ``mpmath.iv``.

    The pinned kernel deliberately keeps integer-only Python arithmetic exact;
    mpmath's interval context does not accept a Fraction produced by an exact
    ``1/4`` operation on every supported version.  This adapter preserves the
    same ordinary jet algebra while converting only numeric constants to
    directed intervals.
    """

    import mpmath

    iv = mpmath.iv

    def cast(value: Any) -> Any:
        if isinstance(value, Fraction):
            return iv.mpf(str(value))
        if isinstance(value, int) and not isinstance(value, bool):
            return iv.mpf(str(value))
        return value

    return ScalarOperations(
        lambda a, b: cast(a) + cast(b),
        lambda a, b: cast(a) - cast(b),
        lambda a, b: cast(a) * cast(b),
        lambda a, b: cast(a) / cast(b),
        lambda a: -cast(a),
        cast,
    )


def _select_operations(
    operations: Optional[ScalarOperations], *values: Any
) -> ScalarOperations:
    explicit: list[ScalarOperations] = []
    for value in values:
        _reject_nested_ordinary_jet(value, "scaled density input")
        if not isinstance(value, OrdinaryJet) or value.operations is None:
            continue
        adapter = value.operations
        if operations is not None and adapter is not operations:
            raise TypeError(
                "ordinary jet uses an explicit scalar operation adapter "
                "different from the requested adapter"
            )
        if not any(adapter is previous for previous in explicit):
            explicit.append(adapter)
    if operations is not None:
        return operations
    if len(explicit) > 1:
        raise TypeError(
            "ordinary jets use multiple explicit scalar operation adapters"
        )
    if explicit:
        return explicit[0]
    if any(_contains_interval(value) for value in values):
        return _interval_operations()
    return PYTHON_OPERATIONS


def _c1(jet: OrdinaryJet) -> OrdinaryJet:
    """Drop an unneeded second row while preserving the scalar adapter."""

    return ordinary_jet(jet.value, jet.dz, operations=jet.operations)


def _sum_jets(jets: Sequence[OrdinaryJet]) -> OrdinaryJet:
    if not jets:
        raise ValueError("at least one ordinary jet is required")
    result = jets[0]
    for jet in jets[1:]:
        result = result + jet
    return result


def _scaled_recovery(
    baseline_V_rows: Any,
    baseline_q_rows: Any,
    delta_V_rows: Any,
    delta_m_rows: Any,
    delta_V_y_rows: Any,
    Z_jet: Any,
    delta: Any,
    *,
    operations: Optional[ScalarOperations] = None,
) -> dict[str, Any]:
    """Build normalized fields and pressure recovery before density sectors."""

    # Select the scalar adapter before validating rows so interval constants
    # are directed throughout the pinned kernel rather than becoming exact
    # Fractions in an intermediate denominator.
    ops = _select_operations(
        operations,
        baseline_V_rows,
        baseline_q_rows,
        delta_V_rows,
        delta_m_rows,
        delta_V_y_rows,
        Z_jet,
        delta,
    )
    # Validate every row before doing any arithmetic.  In particular, dm's
    # fourth row must not be inferred from a lower-order record.
    V = _as_jet(baseline_V_rows, 2, "baseline_V_rows", operations=ops)
    q = _as_jet(baseline_q_rows, 2, "baseline_q_rows", operations=ops)
    dV = _as_jet(delta_V_rows, 2, "delta_V_rows", operations=ops)
    dm_rows = _row(delta_m_rows, 3, "delta_m_rows", operations=ops)
    dVy = _as_jet(delta_V_y_rows, 1, "delta_V_y_rows", operations=ops)
    dM = ordinary_jet(dm_rows[0], dm_rows[1], dm_rows[2], operations=ops)
    # dm_rows was already validated and lifted; do not lift it a second time.
    dM_Z = ordinary_jet(dm_rows[1], dm_rows[2], dm_rows[3], operations=ops)
    z = _coordinate_jet(Z_jet, operations=ops)
    de = _constant_jet(delta, ops)
    one = ordinary_jet(ops.constant(1), ops.constant(0), ops.constant(0), operations=ops)
    two = ordinary_jet(ops.constant(2), ops.constant(0), ops.constant(0), operations=ops)
    D = one - z * z
    L = one - de * z * z

    # Pressure recovery is performed once from the same-source correction.
    # The derivative rows below are ordinary derivatives, not coefficient rows.
    dV_Z = ordinary_jet(dV.dz, dV.dzz, operations=ops)
    dQ = (two * z * dV - z * (one - de) * dM - D * dM_Z) / L
    dQ_y = (z * (one + de) * dV + two * z * dVy - D * dV_Z) / L - dQ

    baseline = physical_fields(V, q, operations=ops)
    correction = physical_fields(dV, dQ, operations=ops)
    return {
        "Z": z,
        "delta": de,
        "D": D,
        "L": L,
        "baseline_V": V,
        "baseline_q": q,
        "delta_V": dV,
        "delta_m": dM,
        "delta_m_Z": dM_Z,
        "delta_V_y": dVy,
        "delta_V_Z": dV_Z,
        "delta_q": dQ,
        "delta_q_y": dQ_y,
        "baseline_fields": baseline,
        "correction_fields": correction,
    }


def correlated_scaled_density(
    baseline_V_rows: Any,
    baseline_q_rows: Any,
    delta_V_rows: Any,
    delta_m_rows: Any,
    delta_V_y_rows: Any,
    Z_jet: Any,
    delta: Any,
    *,
    operations: Optional[ScalarOperations] = None,
) -> dict[str, Any]:
    """Return normalized C2/C1 recovery and signed density sectors.

    The baseline and correction are assembled as correlated normalized source
    fields and then passed to the pinned six-term density kernel.  The first
    two volume terms form ``v1`` (formal ``Pstar`` degree one); the remaining
    four form ``v2`` (formal degree two).  Endpoint terms are split in the
    same way: one linear term ``F1`` and two quadratic terms ``F2``.

    No independent full/background boxes are subtracted, and no physical
    ``Pstar`` or ``N`` exponent is constructed.
    """

    state = _scaled_recovery(
        baseline_V_rows,
        baseline_q_rows,
        delta_V_rows,
        delta_m_rows,
        delta_V_y_rows,
        Z_jet,
        delta,
        operations=operations,
    )
    base = state["baseline_fields"]
    corr = state["correction_fields"]
    z = state["Z"]
    de = state["delta"]

    # These are the six exact factored source terms from the pinned kernel.
    selected_operations = z.operations
    volume_terms = volume_difference_terms(
        base, corr, z, de, operations=selected_operations
    )
    endpoint_terms = endpoint_difference_terms(
        base, corr, z, de, operations=selected_operations
    )
    volume_terms_c1 = {name: _c1(value) for name, value in volume_terms.items()}
    endpoint_terms_c1 = {name: _c1(value) for name, value in endpoint_terms.items()}

    v1 = _sum_jets(
        [volume_terms_c1["linear_q_over_2L"], volume_terms_c1["linear_z_qz_over_4L"]]
    )
    v2 = _sum_jets(
        [
            volume_terms_c1["quadratic_q_squared_over_4"],
            volume_terms_c1["correlated_D_u_qz_over_2L"],
            volume_terms_c1["correlated_z2_u_squared_over_L2"],
            volume_terms_c1["correlated_z_D_u_uz_over_L2"],
        ]
    )
    F1 = endpoint_terms_c1["linear_q_over_2L"]
    F2 = _sum_jets(
        [
            endpoint_terms_c1["quadratic_q_squared_over_4"],
            endpoint_terms_c1["quadratic_minus_z2_u_squared_over_L2"],
        ]
    )

    delta_v = v1 + v2
    delta_F = F1 + F2
    # Return C1 sectors as the public normalized density interface while
    # retaining all six exact factored terms for source inspection.
    result = dict(state)
    result.update(
        {
            "delta_q_C2": state["delta_q"],
            "delta_q_y_C1": _c1(state["delta_q_y"]),
            "delta_qy_C1": _c1(state["delta_q_y"]),
            "delta_Q": state["delta_q"],
            "delta_Q_y": _c1(state["delta_q_y"]),
            "delta_Qy": _c1(state["delta_q_y"]),
            "normalized_delta_Q": state["delta_q"],
            "normalized_delta_Q_y": _c1(state["delta_q_y"]),
            "normalized_delta_Qy": _c1(state["delta_q_y"]),
            "normalized_deltaQ": state["delta_q"],
            "normalized_deltaQy": _c1(state["delta_q_y"]),
            "normalized_source_fields": {
                "baseline": base,
                "correction": corr,
                "V": state["baseline_V"],
                "q": state["baseline_q"],
                "dV": state["delta_V"],
                "dm": state["delta_m"],
                "dV_y": state["delta_V_y"],
                "delta_q": state["delta_q"],
                "delta_q_y": _c1(state["delta_q_y"]),
            },
            "volume_terms": volume_terms_c1,
            "endpoint_terms": endpoint_terms_c1,
            "volume_sectors": {"linear": v1, "quadratic": v2},
            "endpoint_sectors": {"linear": F1, "quadratic": F2},
            "normalized_C1_linear": {"v": v1, "F": F1},
            "normalized_C1_quadratic": {"v": v2, "F": F2},
            "normalized_C1_sectors": {
                "linear": {"v": v1, "F": F1},
                "quadratic": {"v": v2, "F": F2},
            },
            "v1": v1,
            "v2": v2,
            "F1": F1,
            "F2": F2,
            "normalized_delta_v": delta_v,
            "normalized_delta_F": delta_F,
            "delta_v": delta_v,
            "delta_F": delta_F,
            "physical_sector_formula": {
                "delta_v": "Pstar*v1 + Pstar**2*v2",
                "delta_F": "Pstar*F1 + Pstar**2*F2",
                "Pstar_materialized": False,
                "N_physical_exponent_materialized": False,
            },
            "same_source_correction": True,
            "independent_full_background_subtraction": False,
        }
    )
    return result


def evaluate(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Assemble normalized source fields and all signed output sectors."""

    return correlated_scaled_density(*args, **kwargs)


def _full_scaled_difference(
    normalized: dict[str, Any],
    pstar: Any,
    *,
    endpoint: bool = False,
) -> OrdinaryJet:
    """Synthetic control helper: materialize ``Pstar`` only in a test."""

    base = normalized["baseline_fields"]
    corr = normalized["correction_fields"]
    physical_base = physical_fields(pstar * base.U, pstar * base.Q)
    physical_corr = physical_fields(pstar * corr.U, pstar * corr.Q)
    if endpoint:
        return full_endpoint_density_difference(
            physical_base,
            physical_corr,
            normalized["Z"],
            normalized["delta"],
        )
    return full_volume_density_difference(
        physical_base,
        physical_corr,
        normalized["Z"],
        normalized["delta"],
    )


def _interval_contains(interval: Any, point: Any) -> bool:
    try:
        return point in interval
    except Exception:
        # mpmath.iv membership can be conservative for a plain Fraction;
        # converting only the synthetic control point is safe here.
        try:
            import mpmath

            scalar = mpmath.iv.mpf(str(point)) if isinstance(point, Fraction) else mpmath.iv.mpf(point)
            return scalar in interval
        except Exception:
            return False


def run_exact_controls() -> dict[str, Any]:
    """Run synthetic exact rational, symbolic, and signed-interval controls.

    The controls use new algebraic rows only.  They do not load or evaluate
    any archived source function, cap, enclosure, fixture, recurrence, or
    quadrature.  A passing control therefore validates the arithmetic shape,
    not a numerical claim about the actual C5 source.
    """

    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover - repository normally has SymPy
        return {"passed": False, "sympy_available": False, "error": str(exc)}

    z, de, pstar = sp.symbols("Z delta Pstar")
    V = [z**2 + 2 * z + 3, 2 * z + 2, 2]
    q = [z**3 - z + 1, 3 * z**2 - 1, 6 * z]
    dV = [z**2 - 2 * z + 1, 2 * z - 2, 2]
    dm = [z**3 + z**2 - z + 2, 3 * z**2 + 2 * z - 1, 6 * z + 2, 6]
    dVy = [z + 1, 1]
    normalized = correlated_scaled_density(V, q, dV, dm, dVy, ordinary_jet(z, 1, 0), de)
    full_v = _full_scaled_difference(normalized, pstar)
    full_F = _full_scaled_difference(normalized, pstar, endpoint=True)
    split_v = sp.simplify(full_v.value - pstar * normalized["v1"].value - pstar**2 * normalized["v2"].value) == 0
    split_F = sp.simplify(full_F.value - pstar * normalized["F1"].value - pstar**2 * normalized["F2"].value) == 0
    c1_chain = all(
        sp.simplify(jet.dz - sp.diff(jet.value, z)) == 0
        for jet in (normalized["delta_q_C2"], normalized["v1"], normalized["v2"], normalized["F1"], normalized["F2"])
    )

    # Same-source linearity of the recovered pressure correction, with exact
    # rational rows and a nonzero baseline.  This is separate from the full
    # density split so a hidden independent-box subtraction cannot pass.
    rz = ordinary_jet(Fraction(1, 3), 1, 0)
    rd = Fraction(1, 5)
    rV = [Fraction(2), Fraction(3), Fraction(5)]
    rq = [Fraction(7, 4), Fraction(2), Fraction(1)]
    rVy = [Fraction(1, 2), Fraction(1, 3)]
    dV1 = [Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)]
    dM1 = [Fraction(2, 5), Fraction(1, 3), Fraction(1, 4), Fraction(1, 7)]
    dV2 = [Fraction(3, 5), Fraction(1, 2), Fraction(2, 7)]
    dM2 = [Fraction(1, 6), Fraction(2, 9), Fraction(3, 10), Fraction(4, 11)]
    one = correlated_scaled_density(rV, rq, dV1, dM1, rVy, rz, rd)
    two = correlated_scaled_density(rV, rq, dV2, dM2, rVy, rz, rd)
    both = correlated_scaled_density(
        rV,
        rq,
        [a + b for a, b in zip(dV1, dV2)],
        [a + b for a, b in zip(dM1, dM2)],
        rVy,
        rz,
        rd,
    )
    q_linear = all(
        both["delta_q_C2"].__getattribute__(component)
        == one["delta_q_C2"].__getattribute__(component) + two["delta_q_C2"].__getattribute__(component)
        for component in ("value", "dz", "dzz")
    )

    incomplete_mapping_rejected = False
    try:
        correlated_scaled_density(
            {"value": rV[0], "dzz": rV[2]},
            rq,
            dV1,
            dM1,
            rVy,
            rz,
            rd,
        )
    except ScaledDensityJetError:
        incomplete_mapping_rejected = True
    truncated_dm_rejected = False
    try:
        correlated_scaled_density(rV, rq, dV1, dM1[:3], rVy, rz, rd)
    except ScaledDensityJetError:
        truncated_dm_rejected = True

    interval_control = False
    interval_nonpoint = False
    try:
        import mpmath

        iv = mpmath.iv
        iz = iv.mpf([mpmath.mpf("0.28"), mpmath.mpf("0.32")])
        ide = iv.mpf([mpmath.mpf("0.18"), mpmath.mpf("0.22")])
        # ``iv`` does not accept Fraction objects directly on all mpmath
        # versions; its decimal-string parser does, without routing through
        # a binary float.
        irows = lambda row: [iv.mpf(str(x)) if isinstance(x, Fraction) else iv.mpf(x) for x in row]
        interval = correlated_scaled_density(
            irows(rV), irows(rq), irows(dV1), irows(dM1), irows(rVy),
            [iz, iv.mpf(1), iv.mpf(0)], ide,
        )
        point = correlated_scaled_density(
            rV, rq, dV1, dM1, rVy, rz, rd,
        )
        interval_control = _interval_contains(interval["delta_q_C2"].value, point["delta_q_C2"].value)
        interval_control = interval_control and _interval_contains(interval["v1"].value, point["v1"].value)
        interval_control = interval_control and _interval_contains(interval["F2"].value, point["F2"].value)
        point_v1_interval = iv.mpf(str(point["v1"].value)) if isinstance(point["v1"].value, Fraction) else iv.mpf(point["v1"].value)
        interval_nonpoint = str(interval["v1"].value) != str(point_v1_interval)
    except Exception:
        interval_control = False
        interval_nonpoint = False

    passed = all(
        (
            split_v,
            split_F,
            c1_chain,
            q_linear,
            incomplete_mapping_rejected,
            truncated_dm_rejected,
            interval_control,
            interval_nonpoint,
        )
    )
    return {
        "passed": passed,
        "sympy_available": True,
        "symbolic_Pstar_volume_split": split_v,
        "symbolic_Pstar_endpoint_split": split_F,
        "C1_chain_rule": c1_chain,
        "same_source_delta_q_linearity": q_linear,
        "incomplete_mapping_rejected": incomplete_mapping_rejected,
        "truncated_dm_rejected": truncated_dm_rejected,
        "signed_interval_function_enclosure": interval_control,
        "interval_function_enclosure_control": interval_control,
        "bound_vs_function_enclosure_distinguished": interval_control and interval_nonpoint,
        "interval_is_not_point_value": interval_nonpoint,
        "bounds_are_not_function_values": True,
        "synthetic_controls_only": True,
        "actual_sources_numeric_values_claimed": False,
        "physical_Pstar_materialized_in_kernel": False,
        "physical_N_exponent_materialized_in_kernel": False,
    }


__all__ = [
    "ScaledDensityJetError",
    "explicit_shift_jet",
    "explicitshiftjet",
    "correlated_scaled_density",
    "evaluate",
    "run_exact_controls",
]
