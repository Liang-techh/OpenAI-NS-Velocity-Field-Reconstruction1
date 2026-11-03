"""Independent synthetic fixture for the live-axis potential adapter.

The runner passes the adapter module once to :func:`run_fixture`.  This file
does not import the implementation or any repository numerical helper.  It
creates fresh mock weighted/physical/raw objects and returns a JSON-safe
report covering ownership, ordinary Z1 extraction, operator signs, shared
snapshot use, caller-frame declaration guards, graph replacement refusal,
public snapshot mutation isolation, bound-coordinate override refusal,
continuous coefficient mappings, and vector shape rejection.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Callable


__all__ = ["run_fixture", "run"]


class _Scalar:
    """Tiny tagged scalar used to detect zero-type preservation."""

    __slots__ = ("value", "tag")

    def __init__(self, value: int, tag: str = "synthetic") -> None:
        self.value = value
        self.tag = tag

    def __mul__(self, other: Any) -> "_Scalar":
        if isinstance(other, _Scalar):
            return _Scalar(self.value * other.value, self.tag)
        return _Scalar(self.value * other, self.tag)

    __rmul__ = __mul__

    def __add__(self, other: Any) -> "_Scalar":
        if isinstance(other, _Scalar):
            return _Scalar(self.value + other.value, self.tag)
        return _Scalar(self.value + other, self.tag)

    def __sub__(self, other: Any) -> "_Scalar":
        if isinstance(other, _Scalar):
            return _Scalar(self.value - other.value, self.tag)
        return _Scalar(self.value - other, self.tag)

    def __neg__(self) -> "_Scalar":
        return _Scalar(-self.value, self.tag)

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, _Scalar):
            return self.value == other.value and self.tag == other.tag
        if other == 0:
            return self.value == 0
        return False

    def __repr__(self) -> str:
        return f"Scalar({self.value},{self.tag})"


class _Chi:
    def __init__(self, coefficients: Any) -> None:
        self.coefficients = (
            dict(coefficients) if isinstance(coefficients, Mapping)
            else list(coefficients)
        )


class _AxisRecord:
    def __init__(self, chi: _Chi) -> None:
        self.chi = chi


class _Physical:
    def __init__(self, rows: Mapping[Any, _AxisRecord]) -> None:
        self.rows = dict(rows)
        self.calls: list[Any] = []

    def axis_inputs(self, Z: Any) -> _AxisRecord:
        self.calls.append(Z)
        return self.rows[Z]


class _Box:
    def __init__(self, **values: Any) -> None:
        self.__dict__.update(values)


class _Weighted:
    def __init__(self, physical: _Physical) -> None:
        self.physical = physical
        self.inputs = _Box(n1=_Box(source=_Box(core=physical)))


class _RawVolterra:
    def __init__(self) -> None:
        self.calls = 0

    def __getattr__(self, name: str) -> Any:
        raise AssertionError(f"fixture raw owner was unexpectedly invoked: {name}")


def _jsonable(value: Any) -> Any:
    """Convert fixture values to JSON-safe primitives without importing JSON."""

    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, _Scalar):
        return {"value": value.value, "tag": value.tag}
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_jsonable(item) for item in value]
    if isinstance(value, BaseException):
        return {"type": type(value).__name__, "message": str(value)}
    return repr(value)


def _resolve_factory(adapter_module: Any) -> Callable[..., Any]:
    if callable(adapter_module):
        return adapter_module
    for name in (
        "build",
        "build_live_axis_potential_adapter",
        "make_adapter",
        "construct",
    ):
        candidate = getattr(adapter_module, name, None)
        if callable(candidate):
            return candidate
    raise TypeError("adapter_module must expose a callable build entry point")


def _capture(
    checks: dict[str, bool],
    cases: dict[str, dict[str, Any]],
    name: str,
    passed: bool,
    *,
    observed: Any = None,
    expected: Any = None,
    error: Any = None,
) -> None:
    checks[name] = bool(passed)
    cases[name] = {
        "passed": bool(passed),
        "observed": _jsonable(observed),
        "expected": _jsonable(expected),
        "error": _jsonable(error),
    }


def run_fixture(adapter_module: Any) -> dict[str, Any]:
    """Run all synthetic checks once and return a JSON-serializable report."""

    factory = _resolve_factory(adapter_module)
    scalar = lambda value: _Scalar(value)
    rows = {
        "z0": _AxisRecord(_Chi([scalar(3), scalar(5), scalar(7), scalar(-1)])),
        "z1": _AxisRecord(_Chi([scalar(-4), scalar(6), scalar(9)])),
    }
    physical = _Physical(rows)
    weighted = _Weighted(physical)
    raw = _RawVolterra()
    proof = {"frame": "weighted_axis", "declared": True}
    adapter = factory(weighted, raw, frame_proof=proof)

    checks: dict[str, bool] = {}
    cases: dict[str, dict[str, Any]] = {}
    vector = [scalar(2), scalar(11), scalar(13), scalar(17), scalar(19), scalar(23)]

    calls_before_materialize = len(physical.calls)
    pair = adapter.materialize("z0")
    calls_after_materialize = len(physical.calls)
    A = pair["A"]
    A_Z = pair["A_Z"]
    output_A = A(vector)
    output_A_Z = A_Z(vector)
    calls_after_operators = len(physical.calls)
    jets = pair["axis_jets"]
    snapshot_before_public_mutation = A.snapshot
    history_before_public_mutation = list(jets["coefficient_history"])
    public_mutation_output = None
    try:
        jets["chi"] = scalar(99)
        jets["chi_Z"] = scalar(101)
        jets["coefficients"][0] = scalar(77)
        jets["coefficient_history"][1] = scalar(88)
        public_mutation_output = A(vector)
    except Exception as exc:
        public_mutation_output = exc

    bound_Z_error: Any = None
    try:
        A(vector, "z1")
    except Exception as exc:
        bound_Z_error = exc

    _capture(
        checks,
        cases,
        "operator_signature",
        callable(adapter.A)
        and callable(adapter.A_Z)
        and callable(A)
        and callable(A_Z)
        and A.owner is adapter
        and A_Z.owner is adapter,
        observed={
            "root_A_callable": callable(adapter.A),
            "root_A_Z_callable": callable(adapter.A_Z),
            "bound_A_owner": A.owner is adapter,
            "bound_A_Z_owner": A_Z.owner is adapter,
        },
        expected="two callable operators owned by one adapter",
    )
    _capture(
        checks,
        cases,
        "materialize_reads_one_shared_axis_snapshot",
        calls_after_materialize - calls_before_materialize == 1
        and calls_after_operators == calls_after_materialize
        and pair["axis_jets"] is pair["jets"]
        and A.snapshot is A_Z.snapshot,
        observed={
            "calls_before_materialize": calls_before_materialize,
            "calls_after_materialize": calls_after_materialize,
            "calls_after_operators": calls_after_operators,
            "shared_snapshot": A.snapshot is A_Z.snapshot,
        },
        expected="one axis_inputs read shared by A and A_Z",
    )
    _capture(
        checks,
        cases,
        "public_jets_mutation_does_not_change_snapshot",
        not isinstance(public_mutation_output, BaseException)
        and A.snapshot is snapshot_before_public_mutation
        and A_Z.snapshot is snapshot_before_public_mutation
        and public_mutation_output[4].value == -12
        and A.snapshot[0] is snapshot_before_public_mutation[0]
        and A.snapshot[1] is snapshot_before_public_mutation[1],
        observed={
            "public_chi": jets["chi"],
            "public_chi_Z": jets["chi_Z"],
            "snapshot_chi": A.snapshot[0],
            "snapshot_chi_Z": A.snapshot[1],
            "output_after_public_mutation": public_mutation_output,
        },
        expected="detached public jets leave the shared immutable pair unchanged",
    )
    _capture(
        checks,
        cases,
        "materialized_operator_rejects_coordinate_override",
        bound_Z_error is not None
        and type(bound_Z_error).__name__ == "TypeError",
        observed=bound_Z_error,
        expected="TypeError",
    )
    _capture(
        checks,
        cases,
        "ordinary_Z1_without_factorial",
        snapshot_before_public_mutation[1] is rows["z0"].chi.coefficients[1]
        and history_before_public_mutation[1] is rows["z0"].chi.coefficients[1]
        and jets["factorial_denominator"] == 1
        and jets["ordinary_Z1_factorial_divided"] is False,
        observed={
            "chi_Z": snapshot_before_public_mutation[1],
            "ordinary_Z1": snapshot_before_public_mutation[1],
            "factorial_denominator": jets["factorial_denominator"],
            "history": history_before_public_mutation,
        },
        expected={"chi_Z": 5, "factorial_denominator": 1},
    )
    _capture(
        checks,
        cases,
        "owner_identity_and_exact_raw_handle",
        adapter.physical is physical
        and adapter.axis_owner is physical
        and adapter.raw_volterra_owner is raw
        and adapter.source_flags["same_live_axis_owner_proved"] is True,
        observed={
            "physical_is_exact": adapter.physical is physical,
            "axis_owner_is_exact": adapter.axis_owner is physical,
            "raw_is_exact": adapter.raw_volterra_owner is raw,
            "same_live_axis_owner_proved": adapter.source_flags[
                "same_live_axis_owner_proved"
            ],
        },
        expected="exact physical and raw handles retained",
    )
    _capture(
        checks,
        cases,
        "operator_sign_and_fifth_component",
        output_A[4].value == -12
        and output_A_Z[4].value == -20
        and all(item.value == 0 for index, item in enumerate(output_A) if index != 4)
        and all(item.value == 0 for index, item in enumerate(output_A_Z) if index != 4),
        observed={"A": output_A, "A_Z": output_A_Z},
        expected={"A_fifth": -12, "A_Z_fifth": -20, "other_components": 0},
    )
    _capture(
        checks,
        cases,
        "compatible_zero_type",
        all(isinstance(item, _Scalar) for item in output_A)
        and all(isinstance(item, _Scalar) for item in output_A_Z),
        observed={"A_types": [type(item).__name__ for item in output_A]},
        expected="all six output values retain synthetic scalar type",
    )
    _capture(
        checks,
        cases,
        "source_flags_do_not_claim_raw_merge_or_same_source",
        adapter.source_flags["raw_weighted_auto_merged"] is False
        and adapter.source_flags["same_source_proof_established"] is False
        and adapter.source_flags["real_same_source_proof_established"] is False
        and adapter.source_flags["same_initial_datum_proved"] is False
        and adapter.source_flags["full_Kchi_proved"] is False
        and adapter.source_flags["real_ns_verified"] is False,
        observed={
            key: adapter.source_flags[key]
            for key in (
                "raw_weighted_auto_merged",
                "same_source_proof_established",
                "real_same_source_proof_established",
                "same_initial_datum_proved",
                "full_Kchi_proved",
                "real_ns_verified",
            )
        },
        expected="all unproved source identities remain false",
    )
    _capture(
        checks,
        cases,
        "caller_frame_declaration_is_not_actual_proof",
        adapter.source_flags["caller_frame_declared"] is True
        and adapter.source_flags["explicit_frame_proved"] is False
        and adapter.frame_proven is False,
        observed={
            "caller_frame_declared": adapter.source_flags["caller_frame_declared"],
            "explicit_frame_proved": adapter.source_flags["explicit_frame_proved"],
            "frame_proven_property": adapter.frame_proven,
        },
        expected={"caller_frame_declared": True, "explicit_frame_proved": False},
    )

    replacement = _Physical(rows)
    weighted.physical = replacement
    replacement_error: Any = None
    try:
        A(vector)
    except Exception as exc:
        replacement_error = exc
    _capture(
        checks,
        cases,
        "graph_replacement_refused",
        replacement_error is not None
        and type(replacement_error).__name__ == "GraphReplacementError",
        observed=replacement_error,
        expected="GraphReplacementError",
    )

    # Build a fresh graph for the remaining negative checks after replacement.
    physical_negative = _Physical(rows)
    weighted_negative = _Weighted(physical_negative)
    unproved = factory(weighted_negative, _RawVolterra())
    unproved_operator = unproved.materialize("z0")["A"]
    frame_error: Any = None
    try:
        unproved_operator(vector)
    except Exception as exc:
        frame_error = exc
    _capture(
        checks,
        cases,
        "unproved_frame_fails_closed",
        frame_error is not None and type(frame_error).__name__ == "FrameProofError",
        observed=frame_error,
        expected="FrameProofError",
    )

    declared_frame_adapter = factory(weighted_negative, _RawVolterra(), frame_proof=proof)
    raw_frame_error: Any = None
    try:
        declared_frame_adapter.materialize("z0")["A"](vector, frame="raw_volterra")
    except Exception as exc:
        raw_frame_error = exc
    _capture(
        checks,
        cases,
        "raw_frame_cannot_authorize_operator",
        raw_frame_error is not None and type(raw_frame_error).__name__ == "FrameProofError",
        observed=raw_frame_error,
        expected="FrameProofError",
    )

    noncanonical_frame_error: Any = None
    try:
        factory(
            _Weighted(_Physical(rows)),
            _RawVolterra(),
            frame_proof={"frame": "caller_weighted_axis", "declared": True},
        )
    except Exception as exc:
        noncanonical_frame_error = exc
    _capture(
        checks,
        cases,
        "noncanonical_frame_declaration_rejected",
        noncanonical_frame_error is not None
        and type(noncanonical_frame_error).__name__ == "FrameProofError",
        observed=noncanonical_frame_error,
        expected="FrameProofError",
    )

    mapping_physical = _Physical(
        {"bad": _AxisRecord(_Chi({0: scalar(3), 2: scalar(5)}))}
    )
    mapping_error: Any = None
    try:
        factory(_Weighted(mapping_physical), _RawVolterra(), frame_proof=proof).axis_jets("bad")
    except Exception as exc:
        mapping_error = exc
    _capture(
        checks,
        cases,
        "coefficient_mapping_gap_rejected",
        mapping_error is not None
        and type(mapping_error).__name__ == "AxisInterfaceError",
        observed=mapping_error,
        expected="AxisInterfaceError",
    )

    vector_error: Any = None
    try:
        A(vector[:2])
    except Exception as exc:
        vector_error = exc
    _capture(
        checks,
        cases,
        "vector_length_rejected",
        vector_error is not None and type(vector_error).__name__ == "VectorShapeError",
        observed=vector_error,
        expected="VectorShapeError",
    )

    same_raw_error: Any = None
    try:
        factory(_Weighted(physical_negative), physical_negative, frame_proof=proof)
    except Exception as exc:
        same_raw_error = exc
    _capture(
        checks,
        cases,
        "raw_owner_cannot_be_physical_owner",
        same_raw_error is not None and type(same_raw_error).__name__ == "RawOwnerError",
        observed=same_raw_error,
        expected="RawOwnerError",
    )

    flags = dict(unproved.source_flags)
    report = {
        "schema": "n1_live_axis_potential_adapter_fixture_v1",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "cases": cases,
        "source_flags": _jsonable(flags),
        "axis_calls": _jsonable(physical.calls),
        "raw_owner_calls": raw.calls,
        "synthetic_only": True,
        "real_objects_used": False,
        "old_bridge_used": False,
        "returned_data_complete": True,
    }
    return _jsonable(report)


run = run_fixture

