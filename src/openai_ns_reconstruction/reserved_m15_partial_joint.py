"""Pure caller-owned interval arithmetic for the partial joint leaf.

The caller supplies interval objects and their arithmetic implementation.  No
endpoint arithmetic, file access, value selection, or source reconstruction is
performed here.

For each branch the caller-established identity is

    C - T = own + r2 * J - (1 - r2) * T,

where ``T = -P`` and ``J = -inner - lower - pre``.  ``target_C0_C1`` and
``joint_C0_C1`` are ordered ``(C0, C1)``.  Each source cell supplies its own
corresponding pair, and the caller's hull is applied separately across cells
for C0 and C1.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, Callable


IntervalObject = Any
CallerHull = Callable[[Iterable[IntervalObject]], IntervalObject]


def _pair(value: Any, *, name: str) -> tuple[IntervalObject, IntervalObject]:
    """Read a caller-supplied two-branch pair without inspecting endpoints."""

    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a length-2 (C0, C1) pair")
    try:
        size = len(value)
        first = value[0]
        second = value[1]
    except (TypeError, IndexError, KeyError) as exc:
        raise TypeError(f"{name} must be a length-2 (C0, C1) pair") from exc
    if size != 2:
        raise ValueError(f"{name} must be a length-2 (C0, C1) pair")
    return first, second


def _cell_pair(cell: Mapping[str, Any]) -> tuple[Any, tuple[Any, Any]]:
    """Extract one source cell without modifying the caller's mapping."""

    try:
        cell_id = cell["cell_id"]
        own_pair = _pair(cell["own_C0_C1"], name="source cell own_C0_C1")
    except KeyError as exc:
        raise KeyError("each source cell needs cell_id and own_C0_C1") from exc
    return cell_id, own_pair


def _terms(
    *,
    r2: IntervalObject,
    own: IntervalObject,
    target: IntervalObject,
    joint: IntervalObject,
) -> dict[str, Any]:
    """Keep the three formula terms and their sum using interval operators."""

    # These operators are intentionally dispatched to the caller's interval
    # object.  No endpoint is extracted, converted, or rounded here.
    r2_times_joint = r2 * joint
    minus_complement_r2_times_target = -((1 - r2) * target)
    residual = own + r2_times_joint + minus_complement_r2_times_target
    return {
        "own": own,
        "r2_times_joint": r2_times_joint,
        "minus_complement_r2_times_target": minus_complement_r2_times_target,
        "terms": (
            own,
            r2_times_joint,
            minus_complement_r2_times_target,
        ),
        "residual": residual,
    }


def assemble_partial_joint(
    *,
    r2: Any,
    source_cells: Iterable[Mapping[str, Any]],
    target_C0_C1: Any,
    joint_C0_C1: Any,
    hull: CallerHull,
    source_identity: Any,
    fixed_unit: bool,
    r2_z_constant: bool,
    error_ledger: Any,
) -> dict[str, Any]:
    """Assemble per-cell partial residuals and separate C0/C1 hulls.

    API and shapes:

    * ``target_C0_C1`` and ``joint_C0_C1`` are caller interval pairs ordered
      ``(C0, C1)``;
    * each item of ``source_cells`` is a mapping with ``cell_id`` and
      ``own_C0_C1``, where the latter is the corresponding own pair;
    * ``r2`` is one caller interval shared by both branches;
    * ``hull`` receives an iterable of residual interval objects.  It is called
      once for all source-cell C0 residuals and once for all source-cell C1
      residuals.  Those two hulls are never combined.

    The caller must establish the common-source relation
    ``F = inner + lower + pre + P``, ``T = -P``, and
    ``J = -inner - lower - pre``.  ``source_identity`` identifies that source;
    ``fixed_unit`` and ``r2_z_constant`` are required proof markers.  For the
    ordinary Z1 branch uses the same value of ``r2`` because its Z derivative
    is zero.  ``error_ledger`` is returned unchanged.
    """

    if source_identity is None or (
        isinstance(source_identity, str) and not source_identity
    ):
        raise ValueError("source_identity is required")
    if fixed_unit is not True:
        raise ValueError("fixed_unit proof marker must be True")
    if r2_z_constant is not True:
        raise ValueError("r2_z_constant proof marker must be True")
    if not callable(hull):
        raise TypeError("hull must be the caller's interval hull callable")

    target_c0, target_c1 = _pair(target_C0_C1, name="target_C0_C1")
    joint_c0, joint_c1 = _pair(joint_C0_C1, name="joint_C0_C1")
    cell_results: list[dict[str, Any]] = []
    c0_residuals: list[IntervalObject] = []
    c1_residuals: list[IntervalObject] = []

    for source_cell in source_cells:
        cell_id, (own_c0, own_c1) = _cell_pair(source_cell)
        if any(item['cell_id'] == cell_id for item in cell_results):
            raise ValueError(f"duplicate source cell id: {cell_id!r}")

        c0 = _terms(
            r2=r2,
            own=own_c0,
            target=target_c0,
            joint=joint_c0,
        )
        c1 = _terms(
            r2=r2,
            own=own_c1,
            target=target_c1,
            joint=joint_c1,
        )
        cell_results.append({"cell_id": cell_id, "own_C0_C1": (own_c0, own_c1),
            "joint_term_C0_C1": (c0['r2_times_joint'], c1['r2_times_joint']),
            "target_term_C0_C1": (c0['minus_complement_r2_times_target'], c1['minus_complement_r2_times_target']),
            "residual_C0_C1": (c0['residual'], c1['residual'])})
        c0_residuals.append(c0["residual"])
        c1_residuals.append(c1["residual"])

    if not cell_results:
        raise ValueError("source_cells must contain at least one source cell")

    hull_c0 = hull(c0_residuals)
    hull_c1 = hull(c1_residuals)

    return {
        "cell_results": cell_results,
        "hull_C0_C1": (hull_c0, hull_c1),
        "error_ledger": error_ledger,
        "source_identity": source_identity,
        "fixed_unit": fixed_unit,
        "r2_z_constant": r2_z_constant,
    }


__all__ = ["assemble_partial_joint"]
