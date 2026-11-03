"""Live weighted-axis source adapter for the n=1 potential operators.

This module is deliberately a small boundary object.  The caller supplies a
live weighted propagator and a separately owned raw Volterra handle.  The
adapter retains the exact physical axis owner, reads its ``axis_inputs(Z)``
interface at execution time, and materializes only the two source operators
needed by the n=1 potential row::

    A(v)   = (0, 0, 0, 0, -2*chi*v[0],   0)
    A_Z(v) = (0, 0, 0, 0, -2*chi_Z*v[0], 0)

The fifth component is index four in the six-component vector.  ``chi_Z`` is
the ordinary coefficient at index one of ``chi.coefficients``; it is not
divided by a factorial.  All other components are made as zero-like values
from the input's first component.

No old constructor, graph assertion, evaluator, bridge, numerical helper, or
raw Volterra operation is called here.  The raw handle is retained only as a
provenance boundary.  A source-equality proof, an actual frame proof, a common
initial datum, a complete Kchi proof, and a real Navier--Stokes verification
remain explicitly unclaimed.  A caller may declare the exact symbolic frame
label ``weighted_axis``; that declaration is not evidence that the runtime
object is in that frame.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


__all__ = [
    "SCHEMA",
    "WEIGHTED_AXIS_FRAME",
    "LiveAxisPotentialAdapter",
    "LiveAxisPotentialError",
    "SourceOwnershipError",
    "RawOwnerError",
    "AxisInterfaceError",
    "GraphReplacementError",
    "FrameProofError",
    "VectorShapeError",
    "build",
    "build_live_axis_potential_adapter",
    "make_adapter",
    "construct",
]


SCHEMA = "n1_live_axis_potential_adapter_v1"
WEIGHTED_AXIS_FRAME = "weighted_axis"

_MISSING = object()
_UNBOUND = object()


class LiveAxisPotentialError(ValueError):
    """Base error for a live-axis adapter boundary failure."""


class SourceOwnershipError(LiveAxisPotentialError):
    """The weighted physical object is not the declared live axis owner."""


class RawOwnerError(LiveAxisPotentialError):
    """The raw Volterra handle is missing or is not independent."""


class AxisInterfaceError(LiveAxisPotentialError):
    """The weighted physical object does not expose the required axis API."""


class GraphReplacementError(LiveAxisPotentialError):
    """The live weighted graph replaced an owner captured by this adapter."""


class FrameProofError(LiveAxisPotentialError):
    """An operator call lacks the required weighted-axis declaration."""


class VectorShapeError(LiveAxisPotentialError):
    """An operator was called with something other than a six-vector."""


def _member(value: Any, name: str, label: str) -> Any:
    """Read one named field without invoking any legacy source operation."""

    if isinstance(value, Mapping):
        try:
            if name in value:
                return value[name]
        except (TypeError, ValueError):
            pass
    try:
        return getattr(value, name)
    except AttributeError as exc:
        raise AxisInterfaceError(f"{label} is missing {name!r}") from exc


def _try_member(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        try:
            if name in value:
                return value[name]
        except (TypeError, ValueError):
            pass
    try:
        return getattr(value, name)
    except AttributeError:
        return _MISSING


def _read_live_owner(weighted: Any) -> tuple[Any, Any]:
    """Read physical and ``inputs.n1.source.core`` handles by identity."""

    physical = _member(weighted, "physical", "weighted propagator")
    inputs = _member(weighted, "inputs", "weighted propagator")
    n1 = _member(inputs, "n1", "weighted.inputs")
    source = _member(n1, "source", "weighted.inputs.n1")
    core = _member(source, "core", "weighted.inputs.n1.source")
    return physical, core


def _is_vector(value: Any) -> bool:
    if isinstance(value, (str, bytes, bytearray, Mapping)):
        return False
    try:
        return len(value) == 6 and hasattr(value, "__getitem__")
    except (TypeError, ValueError, AttributeError):
        return False


def _zero_like(value: Any) -> Any:
    """Make a zero using exact-zero multiplication, or fail closed."""

    try:
        return value * 0
    except Exception:
        right_error: Exception | None = None
    try:
        return 0 * value
    except Exception as exc:
        right_error = exc
    raise AxisInterfaceError(
        "vector first component cannot produce a compatible exact-zero product"
    ) from right_error


def _negative_two_times(value: Any) -> Any:
    """Form ``-2*value`` without converting or selecting a point value."""

    try:
        return -2 * value
    except Exception:
        # A scalar that implements multiplication only on the right is still
        # compatible with the requested exact sign and scale.
        return -(value * 2)


def _same_nonzero(value: Any) -> bool:
    """Return a conservative nonzero test for history bookkeeping only."""

    try:
        result = value != 0
    except Exception:
        return True
    if result is NotImplemented:
        return True
    try:
        return bool(result)
    except Exception:
        return True


def _coefficient_history(coefficients: Any, label: str) -> list[Any]:
    """Copy the coefficient container's references without normalizing them."""

    if isinstance(coefficients, Mapping):
        # Integer and integer-string spellings are accepted, but the mapping
        # must describe exactly 0, 1, ..., n.  Rejecting gaps and aliases
        # prevents an accidental map order from changing the ordinary Z1 row.
        pairs: dict[int, Any] = {}
        for key, value in coefficients.items():
            if not isinstance(key, (int, str)) or isinstance(key, bool):
                raise AxisInterfaceError(
                    f"{label} mapping keys must be integer orders"
                )
            try:
                order = int(key)
            except (TypeError, ValueError) as exc:
                raise AxisInterfaceError(
                    f"{label} mapping keys must be ordinary integer orders"
                ) from exc
            if order < 0 or order in pairs:
                raise AxisInterfaceError(
                    f"{label} mapping keys must be unique nonnegative orders"
                )
            pairs[order] = value
        expected = list(range(len(pairs)))
        actual = sorted(pairs)
        if actual != expected:
            raise AxisInterfaceError(
                f"{label} mapping keys must be the continuous sequence 0..n"
            )
        return [pairs[order] for order in expected]
    if (
        isinstance(coefficients, Sequence)
        and not isinstance(coefficients, (str, bytes, bytearray))
    ):
        return list(coefficients)
    if not isinstance(coefficients, (str, bytes, bytearray, Mapping)):
        try:
            return list(coefficients)
        except TypeError:
            pass
    raise AxisInterfaceError(f"{label} must be a coefficient sequence or mapping")


def _vector_like(original: Any, values: list[Any]) -> Any:
    """Retain ordinary list/tuple vector containers where it is unambiguous."""

    if isinstance(original, tuple):
        return tuple(values)
    if isinstance(original, list):
        return values
    try:
        return type(original)(values)
    except Exception:
        return tuple(values)


def _frame_declaration_record(frame_proof: Any) -> tuple[bool, str | None]:
    """Parse a caller declaration without calling it actual frame evidence."""

    if frame_proof is None:
        return False, None
    if isinstance(frame_proof, Mapping):
        frame = frame_proof.get("frame", frame_proof.get("name"))
        declared = (
            frame_proof.get("declared") is True
            or frame_proof.get("caller_declared") is True
            or frame_proof.get("proven") is True
            or frame_proof.get("explicit") is True
        )
    elif isinstance(frame_proof, str):
        # A bare exact label is only a caller declaration, never proof.
        frame = frame_proof
        declared = True
    else:
        return False, None
    if not isinstance(frame, str) or not frame:
        return False, None
    if frame != WEIGHTED_AXIS_FRAME:
        raise FrameProofError(
            "only the exact weighted_axis frame declaration is accepted"
        )
    return declared, frame


class _OwnedOperator:
    """Callable operator carrying a stable owner and optional bound Z."""

    __slots__ = ("owner", "name", "bound_Z", "_snapshot")

    def __init__(
        self,
        owner: "LiveAxisPotentialAdapter",
        name: str,
        bound_Z: Any = _UNBOUND,
        snapshot: Any = _UNBOUND,
    ):
        self.owner = owner
        self.name = name
        self.bound_Z = bound_Z
        self._snapshot = snapshot

    @property
    def snapshot(self) -> Any:
        return self._snapshot

    def __repr__(self) -> str:
        suffix = "" if self.bound_Z is _UNBOUND else f"(Z={self.bound_Z!r})"
        return f"<{self.name}{suffix} owner={SCHEMA}>"

    def __call__(
        self,
        vector: Any,
        Z: Any = _UNBOUND,
        *,
        frame: Any = _UNBOUND,
    ) -> Any:
        # Unbound operators support ``op(vector, Z)`` and the unambiguous
        # ``op(Z, vector)`` form.  A materialized operator is fixed at its
        # captured Z and therefore rejects any second coordinate argument.
        if self.bound_Z is _UNBOUND:
            if Z is _UNBOUND:
                raise TypeError(f"{self.name} requires Z when not materialized")
            if not _is_vector(vector) and _is_vector(Z):
                vector, Z = Z, vector
            selected_Z = Z
        else:
            if Z is not _UNBOUND:
                raise TypeError(
                    f"{self.name} is materialized at Z and rejects a second Z"
                )
            selected_Z = self.bound_Z

        if not _is_vector(vector):
            raise VectorShapeError(f"{self.name} requires a six-component vector")
        self.owner._require_frame(frame)
        self.owner._assert_graph_stable()
        if self._snapshot is _UNBOUND:
            jets = self.owner.axis_jets(selected_Z)
            coefficient = jets["chi"] if self.name == "A" else jets["chi_Z"]
        else:
            coefficient = self._snapshot[0] if self.name == "A" else self._snapshot[1]
        first = vector[0]
        zero_values = [_zero_like(first) for _ in range(6)]
        scaled = _negative_two_times(coefficient)
        zero_values[4] = scaled * first
        return _vector_like(vector, zero_values)


class LiveAxisPotentialAdapter:
    """Identity-preserving adapter around one live weighted physical axis."""

    __slots__ = (
        "_weighted",
        "_physical",
        "_raw_volterra_owner",
        "_caller_frame_declared",
        "_frame_name",
        "_source_flags",
        "A",
        "A_Z",
    )

    def __init__(
        self,
        weighted: Any,
        raw_volterra_owner: Any,
        frame_proof: Any = None,
    ) -> None:
        if weighted is None:
            raise SourceOwnershipError("weighted propagator is required")
        if raw_volterra_owner is None:
            raise RawOwnerError("an independent raw Volterra owner is required")

        physical, core = _read_live_owner(weighted)
        if physical is not core:
            raise SourceOwnershipError(
                "weighted.physical is not weighted.inputs.n1.source.core"
            )
        if raw_volterra_owner is physical or raw_volterra_owner is weighted:
            raise RawOwnerError(
                "raw Volterra owner must remain independent from weighted physical graph"
            )

        caller_frame_declared, frame_name = _frame_declaration_record(frame_proof)
        self._weighted = weighted
        self._physical = physical
        self._raw_volterra_owner = raw_volterra_owner
        self._caller_frame_declared = caller_frame_declared
        self._frame_name = frame_name
        axis_inputs = _try_member(physical, "axis_inputs")
        self._source_flags = {
            "schema": SCHEMA,
            "interface_available": callable(axis_inputs),
            "weighted_axis_interface_available": callable(axis_inputs),
            "chi_coefficients_interface_available": False,
            "ordinary_Z1_available": False,
            "nonzero_history_preserved": False,
            "same_live_axis_owner_proved": True,
            "physical_is_weighted_inputs_n1_source_core": True,
            "same_source_proof_established": False,
            "real_same_source_proof_established": False,
            "same_source_proof_not_established": True,
            "real_same_source_proof_not_established": True,
            "source_identity_proved": False,
            "raw_handle_stored": True,
            "raw_owner_is_physical": False,
            "raw_physical_identity_proved": False,
            "raw_physical_same_object_claimed": False,
            "raw_weighted_auto_merged": False,
            "raw_weighted_same_frame_proved": False,
            # A declaration is not runtime proof of an actual frame.
            "explicit_frame_proved": False,
            "caller_frame_declared": caller_frame_declared,
            "frame_name": frame_name,
            "same_initial_datum_proved": False,
            "full_Kchi_proved": False,
            "real_ns_verified": False,
            "old_bridge_called": False,
            "graph_replacement_detected": False,
            "operators_materialized": False,
        }
        # Root operators intentionally remain unbound to Z.  They still carry
        # this one adapter as owner; callers may bind them with ``Z`` or use
        # the pair returned by ``materialize(Z)``.
        self.A = _OwnedOperator(self, "A")
        self.A_Z = _OwnedOperator(self, "A_Z")

    @property
    def weighted(self) -> Any:
        return self._weighted

    @property
    def physical(self) -> Any:
        return self._physical

    @property
    def physical_owner(self) -> Any:
        return self._physical

    @property
    def axis_owner(self) -> Any:
        return self._physical

    @property
    def raw_volterra_owner(self) -> Any:
        return self._raw_volterra_owner

    @property
    def frame_name(self) -> str | None:
        return self._frame_name

    @property
    def frame_proven(self) -> bool:
        # Compatibility spelling retained without ever representing evidence.
        return False

    @property
    def caller_frame_declared(self) -> bool:
        return self._caller_frame_declared

    @property
    def source_flags(self) -> dict[str, Any]:
        flags = dict(self._source_flags)
        try:
            self._assert_graph_stable(update=False)
            flags["graph_stable"] = True
        except GraphReplacementError:
            flags["graph_stable"] = False
        return flags

    def _assert_graph_stable(self, *, update: bool = True) -> None:
        try:
            current_physical, current_core = _read_live_owner(self._weighted)
        except Exception as exc:
            if update:
                self._source_flags["graph_replacement_detected"] = True
            raise GraphReplacementError(
                "weighted live graph no longer exposes its captured axis owner"
            ) from exc
        if current_physical is not self._physical or current_core is not self._physical:
            if update:
                self._source_flags["graph_replacement_detected"] = True
            raise GraphReplacementError(
                "weighted graph replaced the exact physical axis owner"
            )

    def _require_frame(self, frame: Any) -> None:
        if not self._caller_frame_declared or self._frame_name != WEIGHTED_AXIS_FRAME:
            raise FrameProofError(
                "caller declaration for the exact weighted_axis frame is required;"
                " actual frame proof is not established"
            )
        selected = self._frame_name if frame is _UNBOUND else frame
        if selected != WEIGHTED_AXIS_FRAME:
            raise FrameProofError(
                "operator calls are limited to the exact weighted_axis frame label"
            )

    def axis_jets(self, Z: Any) -> dict[str, Any]:
        """Read the current physical axis coefficients at ``Z``.

        The first coefficient is ``chi`` and the second is the ordinary Z1
        coefficient ``chi_Z``.  The whole coefficient history is retained by
        reference, including any later nonzero terms.  No factorial
        normalization is performed.
        """

        self._assert_graph_stable()
        axis_inputs = _try_member(self._physical, "axis_inputs")
        if not callable(axis_inputs):
            raise AxisInterfaceError("physical.axis_inputs must be callable")
        axis_record = axis_inputs(Z)
        chi_record = _member(axis_record, "chi", "physical.axis_inputs(Z)")
        coefficients = _member(chi_record, "coefficients", "axis chi record")
        history = _coefficient_history(coefficients, "axis chi coefficients")
        if len(history) < 2:
            raise AxisInterfaceError(
                "axis chi coefficients require chi and ordinary Z1 entries"
            )
        chi = history[0]
        chi_Z = history[1]
        nonzero_history = [item for item in history if _same_nonzero(item)]
        self._source_flags.update(
            {
                "chi_coefficients_interface_available": True,
                "ordinary_Z1_available": True,
                "nonzero_history_preserved": True,
            }
        )
        return {
            "schema": SCHEMA,
            "Z": Z,
            "chi": chi,
            "chi_Z": chi_Z,
            "Z1": chi_Z,
            "chi_Z1": chi_Z,
            "ordinary_Z1": chi_Z,
            "ordinary_Z_order": 1,
            "factorial_denominator": 1,
            "ordinary_Z1_factorial_divided": False,
            "ordinary_Z1_is_raw_coefficient": True,
            "coefficients": history,
            "coefficient_history": history,
            "nonzero_history": nonzero_history,
            "history_preserved": True,
            "source_flags": self.source_flags,
        }

    def materialize(self, Z: Any) -> dict[str, Any]:
        """Bind ``A`` and ``A_Z`` to one shared live ``axis_jets(Z)`` snapshot.

        The returned operators still check the caller declaration and graph
        identity on every application.  The unbound ``self.A`` and
        ``self.A_Z`` remain dynamic convenience operators; this shared pair is
        the recommended path when the two source rows must stay correlated.
        """

        jets = self.axis_jets(Z)
        self._source_flags["operators_materialized"] = True
        if isinstance(jets.get("source_flags"), dict):
            jets["source_flags"]["operators_materialized"] = True
        # Operators retain only this immutable pair.  The public jet packet is
        # a detached container so callers cannot rewrite the coefficients used
        # by either materialized operator through the returned mapping.
        snapshot = (jets["chi"], jets["chi_Z"])
        public_jets = dict(jets)
        for key in ("coefficients", "coefficient_history", "nonzero_history"):
            value = public_jets.get(key)
            if isinstance(value, list):
                public_jets[key] = list(value)
        if isinstance(public_jets.get("source_flags"), dict):
            public_jets["source_flags"] = dict(public_jets["source_flags"])
        return {
            "schema": SCHEMA,
            "Z": Z,
            "axis_jets": public_jets,
            "jets": public_jets,
            "A": _OwnedOperator(self, "A", Z, snapshot),
            "A_Z": _OwnedOperator(self, "A_Z", Z, snapshot),
            "owner": self,
            "source_flags": self.source_flags,
        }

    def materialize_operators(self, Z: Any) -> tuple[Any, Any]:
        pair = self.materialize(Z)
        return pair["A"], pair["A_Z"]

    def operators(self, Z: Any) -> tuple[Any, Any]:
        return self.materialize_operators(Z)

    def as_record(self) -> dict[str, Any]:
        """Return metadata without pretending opaque handles are serializable."""

        return {
            "schema": SCHEMA,
            "physical_owner": self._physical,
            "raw_volterra_owner": self._raw_volterra_owner,
            "A": self.A,
            "A_Z": self.A_Z,
            "source_flags": self.source_flags,
        }

    def __getitem__(self, key: str) -> Any:
        if key == "A":
            return self.A
        if key == "A_Z":
            return self.A_Z
        if key in {"physical", "physical_owner", "axis_owner"}:
            return self._physical
        if key in {"raw", "raw_volterra_owner"}:
            return self._raw_volterra_owner
        if key == "source_flags":
            return self.source_flags
        if key == "schema":
            return SCHEMA
        raise KeyError(key)


def build(
    weighted: Any,
    raw_volterra_owner: Any,
    frame_proof: Any = None,
    *,
    frame: Any = _UNBOUND,
) -> LiveAxisPotentialAdapter:
    """Construct an identity-preserving live-axis adapter.

    ``frame_proof`` is only a caller declaration, such as
    ``{"frame": "weighted_axis", "declared": True}``; it is never runtime
    proof of the actual frame.  Operator calls are permitted only with the
    exact ``weighted_axis`` label.  The optional ``frame=`` spelling records
    the same caller declaration and does not establish actual frame evidence.
    """

    if frame is not _UNBOUND:
        if frame_proof is not None:
            raise FrameProofError("provide either frame_proof or frame, not both")
        frame_proof = {"frame": frame, "caller_declared": True}
    return LiveAxisPotentialAdapter(weighted, raw_volterra_owner, frame_proof)


build_live_axis_potential_adapter = build
make_adapter = build
construct = build

