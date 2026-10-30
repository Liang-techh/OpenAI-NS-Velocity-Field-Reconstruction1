"""Source-compatible signed finite-N density algebra.

This module contains the ordinary (not Taylor-coefficient) algebra used by the
physical ``m15`` volume density.  The source formula is copied from
``original_macro_m15.half_S1_coefficients`` and is deliberately evaluated on
the correlated baseline and correction before any interval or error boxing::

    L = 1 - delta*Z**2
    D = 1 - Z**2
    v = -Q/(2*L) + Z*(1-delta)*Qz/(4*L) - Q**2/4
        + D*U*Qz/(2*L) + Z**2*(1-delta)*U**2/L**2
        + Z*D*U*Uz/L**2

The source endpoint density is also exported independently::

    F = Q/(2*L) + Q**2/4 - Z**2*U**2/L**2

The physical m15 endpoint remains ``R**2*F - R*Q_y``; the preceding
``-R*Q_y`` pressure inlet is deliberately not folded into F.  The bounded
controls verify the source FTC identity ``S/2 = d_y F + 2*F + v`` under the
source constraint for ``Q_y``.

The symbols ``U``, ``Q``, ``Uz`` and ``Qz`` are source fields, while ``Z`` and
``delta`` are the source geometry fields.  ``delta`` below therefore means the
geometric parameter; a finite-N correction is represented by a separate
``PhysicalFields`` object.  Every jet component is an ordinary derivative
value.  A three-component sequence means ``(value, d/dZ, d^2/dZ^2)`` and is
never interpreted as a Taylor coefficient vector.

The arithmetic accepts Python numbers, exact ``Fraction`` values, SymPy
expressions, and the small exact ``ScalarExpression`` graph used by the C5
source binding.  A caller with a scalar graph that does not overload Python
operators may pass ``ScalarOperations`` callbacks.

No finite-N absolute-error cap is a source field here.  Provenance is checked
only by a caller-supplied guard; when no guard is supplied the report marks the
check as absent rather than manufacturing an identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Number
from collections.abc import Mapping, Sequence
from typing import Any, Callable, Optional


class SourceProvenanceError(ValueError):
    """Raised when a caller's source identity guard rejects an input field."""


class JetFormatError(TypeError):
    """Raised when a field is not an explicit ordinary derivative jet."""


class ScalarOperations:
    """Adapter for scalar expression graphs with no Python operator overloads.

    The callbacks use two scalar arguments except ``constant`` and ``neg``.
    ``constant`` is optional; if omitted, an integer/Fraction is passed through.
    The adapter is intentionally tiny so a C5 graph can provide its existing
    operation constructors without changing this kernel.
    """

    def __init__(self, add: Callable[[Any, Any], Any],
                 sub: Callable[[Any, Any], Any],
                 mul: Callable[[Any, Any], Any],
                 div: Callable[[Any, Any], Any],
                 neg: Callable[[Any], Any] = lambda x: -x,
                 constant: Optional[Callable[[Any], Any]] = None):
        self.add = add
        self.sub = sub
        self.mul = mul
        self.div = div
        self.neg = neg
        self.constant = constant or (lambda x: x)


def _python_div(a: Any, b: Any) -> Any:
    """Keep integer/Fraction controls exact while leaving symbolic scalars alone."""

    if isinstance(a, (int, Fraction)) and isinstance(b, (int, Fraction)):
        return Fraction(a, b)
    return a / b


class _PythonOperations(ScalarOperations):
    def __init__(self):
        super().__init__(lambda a, b: a + b, lambda a, b: a - b,
                         lambda a, b: a * b, _python_div,
                         lambda a: -a, lambda x: x)


PYTHON_OPERATIONS = _PythonOperations()
# ``None`` on an OrdinaryJet means unbound; this singleton is an explicit
# Python arithmetic adapter and is equivalent only by object identity.


def _ops_or_default(operations: Optional[ScalarOperations]) -> ScalarOperations:
    return PYTHON_OPERATIONS if operations is None else operations


@dataclass(frozen=True)
class OrdinaryJet:
    """A scalar and its ordinary first and second Z derivatives.

    ``value`` is f(Z), ``dz`` is f'(Z), and ``dzz`` is f''(Z).  ``dzz`` may be
    ``None`` when the caller only needs C1 output.  In particular, the jet is
    not a truncated Taylor series and no factorial or power-of-Z normalization
    is applied.
    """

    value: Any
    dz: Any = 0
    dzz: Any = None
    operations: Optional[ScalarOperations] = None

    @property
    def first(self) -> Any:
        """Alias for the ordinary first derivative."""

        return self.dz

    @property
    def second(self) -> Any:
        """Alias for the ordinary second derivative."""

        return self.dzz

    @property
    def z(self) -> Any:
        """Alias for the ordinary first derivative."""

        return self.dz

    @property
    def zz(self) -> Any:
        """Alias for the ordinary second derivative."""

        return self.dzz

    @classmethod
    def constant(cls, value: Any,
                 operations: Optional[ScalarOperations] = None) -> "OrdinaryJet":
        ops = _ops_or_default(operations)
        return cls(value, ops.constant(0), ops.constant(0), ops)

    def _ops(self) -> ScalarOperations:
        return _ops_or_default(self.operations)

    def _other(self, other: Any) -> "OrdinaryJet":
        if isinstance(other, OrdinaryJet):
            return other
        if isinstance(other, Number):
            ops = self._ops()
            return OrdinaryJet(ops.constant(other), ops.constant(0),
                               ops.constant(0), ops)
        return OrdinaryJet.constant(other, self._ops())

    def _binary(self, other: Any, operation: str) -> "OrdinaryJet":
        b = self._other(other)
        # If one operand carries an explicit graph adapter, use it for both
        # operands.  Two different explicit adapters are almost certainly a
        # graph-namespace mistake and must not be silently mixed.
        if self.operations is not None and b.operations is not None:
            if self.operations is not b.operations:
                raise TypeError("ordinary jets use different scalar operation adapters")
            ops = self.operations
        else:
            ops = self.operations or b.operations or PYTHON_OPERATIONS
        # Preserve deferred literal conversion and reject adapter relabeling.
        a = _normalize(self, ops)
        b = _normalize(b, ops)
        if operation == "add":
            value = ops.add(a.value, b.value)
            dz = ops.add(a.dz, b.dz)
            dzz = None if a.dzz is None or b.dzz is None else ops.add(a.dzz, b.dzz)
        elif operation == "sub":
            value = ops.sub(a.value, b.value)
            dz = ops.sub(a.dz, b.dz)
            dzz = None if a.dzz is None or b.dzz is None else ops.sub(a.dzz, b.dzz)
        elif operation == "mul":
            value = ops.mul(a.value, b.value)
            dz = ops.add(ops.mul(a.dz, b.value), ops.mul(a.value, b.dz))
            if a.dzz is None or b.dzz is None:
                dzz = None
            else:
                dzz = ops.add(ops.add(ops.mul(a.dzz, b.value),
                                     ops.mul(ops.constant(2), ops.mul(a.dz, b.dz))),
                              ops.mul(a.value, b.dzz))
        elif operation == "div":
            return a * b.inverse()
        else:
            raise ValueError("unknown ordinary-jet operation: " + operation)
        return OrdinaryJet(value, dz, dzz, ops)

    def inverse(self) -> "OrdinaryJet":
        ops = self._ops()
        one = ops.constant(1)
        value = ops.div(one, self.value)
        denominator2 = ops.mul(self.value, self.value)
        dz = ops.neg(ops.div(self.dz, denominator2))
        if self.dzz is None:
            dzz = None
        else:
            denominator3 = ops.mul(denominator2, self.value)
            first = ops.div(ops.mul(ops.constant(2), ops.mul(self.dz, self.dz)), denominator3)
            second = ops.div(self.dzz, denominator2)
            dzz = ops.sub(first, second)
        return OrdinaryJet(value, dz, dzz, ops)

    def __add__(self, other: Any) -> "OrdinaryJet":
        return self._binary(other, "add")

    def __radd__(self, other: Any) -> "OrdinaryJet":
        return self._other(other)._binary(self, "add")

    def __sub__(self, other: Any) -> "OrdinaryJet":
        return self._binary(other, "sub")

    def __rsub__(self, other: Any) -> "OrdinaryJet":
        return self._other(other)._binary(self, "sub")

    def __mul__(self, other: Any) -> "OrdinaryJet":
        return self._binary(other, "mul")

    def __rmul__(self, other: Any) -> "OrdinaryJet":
        return self._other(other)._binary(self, "mul")

    def __truediv__(self, other: Any) -> "OrdinaryJet":
        return self._binary(other, "div")

    def __rtruediv__(self, other: Any) -> "OrdinaryJet":
        return self._other(other)._binary(self, "div")

    def __neg__(self) -> "OrdinaryJet":
        ops = self._ops()
        return OrdinaryJet(ops.neg(self.value), ops.neg(self.dz),
                           None if self.dzz is None else ops.neg(self.dzz), ops)

    def __pow__(self, power: int) -> "OrdinaryJet":
        if type(power) is not int:
            raise TypeError("ordinary jets only support exact integer powers")
        if power == 0:
            return OrdinaryJet.constant(self._ops().constant(1), self._ops())
        if power < 0:
            return (self ** (-power)).inverse()
        out = OrdinaryJet.constant(self._ops().constant(1), self._ops())
        for _ in range(power):
            out = out * self
        return out


def _numeric_constant_jet(value: Number,
                          operations: Optional[ScalarOperations]) -> OrdinaryJet:
    """Construct a numeric literal with adapter-owned value and zero rows."""

    ops = _ops_or_default(operations)
    return OrdinaryJet(ops.constant(value), ops.constant(0),
                       ops.constant(0), ops)


def _lift_numeric_component(value: Any,
                            operations: Optional[ScalarOperations]) -> Any:
    """Lift only a numeric literal; preserve opaque source handles verbatim."""

    if operations is not None and isinstance(value, Number):
        return operations.constant(value)
    return value


@dataclass(frozen=True)
class _UnboundMixedRow(OrdinaryJet):
    """Unbound mapping/sequence row with per-component literal provenance."""

    numeric_components: tuple[bool, bool, bool] = (False, False, False)


@dataclass(frozen=True)
class _UnboundScalarLiteral(OrdinaryJet):
    """Private marker for a bare scalar awaiting graph-constant conversion.

    Ordinary jets constructed explicitly with ``ordinary_jet`` and no adapter
    remain ordinary unbound rows.  Only ``coerce_ordinary_jet`` uses this
    marker, so later adapter binding can distinguish a literal that needs
    ``operations.constant`` from an opaque source value whose rows must be
    preserved verbatim.
    """


def ordinary_jet(value: Any, dz: Any = 0, dzz: Any = None,
                 *, operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Construct an ordinary derivative jet explicitly.

    This constructor is the preferred adapter for C5 graph nodes.  ``dz`` and
    ``dzz`` are derivative values, not Taylor coefficients.
    """

    if isinstance(value, OrdinaryJet):
        if dz != 0 or dzz is not None or operations is not None:
            raise JetFormatError("an existing OrdinaryJet cannot receive extra components")
        return value
    return OrdinaryJet(value, dz, dzz, operations)


def _mapping_component(data: Mapping[str, Any], names: tuple[str, ...], default: Any = None) -> Any:
    for name in names:
        if name in data:
            return data[name]
    return default


def coerce_ordinary_jet(value: Any, *, operations: Optional[ScalarOperations] = None,
                        name: str = "field") -> OrdinaryJet:
    """Convert a source scalar or explicit ordinary-derivative record to a jet.

    Accepted records use ``value`` / ``dz`` / ``dzz`` (or ``z`` / ``zz``), a
    three-element sequence ``(value, dz, dzz)``, or an existing
    :class:`OrdinaryJet`.  A scalar is a constant source field.  Arbitrary
    sequences are rejected so a Taylor polynomial cannot be mistaken for a
    derivative jet.
    """

    if isinstance(value, OrdinaryJet):
        if isinstance(value, (_UnboundScalarLiteral, _UnboundMixedRow)):
            if operations is None or value.operations is operations:
                return value
            return _normalize(value, operations)
        if value.operations is not None and operations is not None and value.operations is not operations:
            raise TypeError(
                name + " uses a different scalar operation adapter; "
                "an explicit adapter cannot be rebound"
            )
        if operations is None or value.operations is operations:
            return value
        return OrdinaryJet(value.value, value.dz, value.dzz, operations)
    if isinstance(value, Mapping):
        raw = _mapping_component(value, ("value", "val", "f"), None)
        if raw is None:
            raise JetFormatError(name + " mapping has no ordinary value component")
        dz = _mapping_component(value, ("dz", "z", "Z", "first", "d_dZ"), 0)
        dzz = _mapping_component(value, ("dzz", "zz", "ZZ", "second", "d2_dZ2"), None)
        components = (raw, dz, dzz)
        if operations is None:
            return _UnboundMixedRow(
                raw, dz, dzz, None,
                tuple(item is not None and isinstance(item, Number)
                      for item in components),
            )
        return OrdinaryJet(
            *(_lift_numeric_component(item, operations) for item in components),
            operations,
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            components = (value[0], value[1], None)
            if operations is None:
                return _UnboundMixedRow(
                    *components, None,
                    tuple(item is not None and isinstance(item, Number)
                          for item in components),
                )
            return OrdinaryJet(
                *(_lift_numeric_component(item, operations) for item in components),
                operations,
            )
        if len(value) == 3:
            components = (value[0], value[1], value[2])
            if operations is None:
                return _UnboundMixedRow(
                    *components, None,
                    tuple(item is not None and isinstance(item, Number)
                          for item in components),
                )
            return OrdinaryJet(
                *(_lift_numeric_component(item, operations) for item in components),
                operations,
            )
        raise JetFormatError(name + " sequence must be (value, dz) or (value, dz, dzz)")
    # A bare scalar with no requested adapter is unbound.  In particular, do
    # not route through ``OrdinaryJet.constant`` in that case: that constructor
    # stores the singleton ``PYTHON_OPERATIONS`` even when its argument is
    # omitted, which would make an ordinary scalar look like an explicit graph
    # adapter and block safe inheritance by a unique source adapter later on.
    # For an explicit adapter, lift numeric literals and zero rows through
    # its constant callback, while preserving opaque source scalars.
    if operations is None:
        # Only numeric literals may be converted with the adapter constant
        # callback. An opaque source scalar retains its identity and rows.
        if isinstance(value, Number):
            return _UnboundScalarLiteral(value, 0, 0, None)
        return OrdinaryJet(value, 0, 0, None)
    if isinstance(value, Number):
        return OrdinaryJet(operations.constant(value), operations.constant(0),
                           operations.constant(0), operations)
    # Preserve opaque source scalars; only their zero rows are constants.
    return OrdinaryJet.constant(value, operations)


def _derivative_jet(jet: OrdinaryJet) -> OrdinaryJet:
    """Return the C1 jet for the ordinary derivative of ``jet``."""

    if isinstance(jet, _UnboundScalarLiteral):
        return _UnboundMixedRow(
            jet.dz, jet.dzz, None, None, (True, True, False)
        )
    if isinstance(jet, _UnboundMixedRow):
        return _UnboundMixedRow(
            jet.dz, jet.dzz, None, None,
            (jet.numeric_components[1], jet.numeric_components[2], False),
        )
    return OrdinaryJet(jet.dz, jet.dzz, None, jet.operations)


def _lookup(data: Any, names: tuple[str, ...], label: str) -> Any:
    if isinstance(data, Mapping):
        result = _mapping_component(data, names, None)
        if result is not None:
            return result
    for name in names:
        if hasattr(data, name):
            return getattr(data, name)
    raise KeyError(label + " source field is missing")


@dataclass(frozen=True)
class PhysicalFields:
    """Correlated physical U/Q source fields and their ordinary Z derivatives."""

    U: OrdinaryJet
    Q: OrdinaryJet
    Uz: OrdinaryJet
    Qz: OrdinaryJet

    @classmethod
    def from_fields(cls, U: Any, Q: Any, *, Uz: Any = None, Qz: Any = None,
                    operations: Optional[ScalarOperations] = None) -> "PhysicalFields":
        u = coerce_ordinary_jet(U, operations=operations, name="U")
        q = coerce_ordinary_jet(Q, operations=operations, name="Q")
        uz = (_derivative_jet(u) if Uz is None else
              coerce_ordinary_jet(Uz, operations=operations, name="Uz"))
        qz = (_derivative_jet(q) if Qz is None else
              coerce_ordinary_jet(Qz, operations=operations, name="Qz"))
        return cls(u, q, uz, qz)

    @classmethod
    def from_mapping(cls, fields: Any, *, operations: Optional[ScalarOperations] = None,
                     label: str = "fields") -> "PhysicalFields":
        U = _lookup(fields, ("U", "u", "dU", "delta_U", "correction_U", "velocity"), label + ".U")
        Q = _lookup(fields, ("Q", "q", "dQ", "delta_Q", "correction_Q", "pressure"), label + ".Q")
        Uz = None
        Qz = None
        if isinstance(fields, Mapping):
            Uz = _mapping_component(fields, ("Uz", "U_z", "u_z", "dUz", "delta_Uz"), None)
            Qz = _mapping_component(fields, ("Qz", "Q_z", "q_z", "dQz", "delta_Qz"), None)
        else:
            for names, target in ((("Uz", "U_z", "u_z", "dUz", "delta_Uz"), "Uz"),
                                  (("Qz", "Q_z", "q_z", "dQz", "delta_Qz"), "Qz")):
                val = next((getattr(fields, n) for n in names if hasattr(fields, n)), None)
                if target == "Uz":
                    Uz = val
                else:
                    Qz = val
        return cls.from_fields(U, Q, Uz=Uz, Qz=Qz, operations=operations)

    def __add__(self, other: "PhysicalFields") -> "PhysicalFields":
        if not isinstance(other, PhysicalFields):
            return NotImplemented
        return PhysicalFields(self.U + other.U, self.Q + other.Q,
                              self.Uz + other.Uz, self.Qz + other.Qz)


def physical_fields(U: Any, Q: Any, *, Uz: Any = None, Qz: Any = None,
                    operations: Optional[ScalarOperations] = None) -> PhysicalFields:
    """Convenience constructor for correlated U/Q source fields."""

    return PhysicalFields.from_fields(U, Q, Uz=Uz, Qz=Qz, operations=operations)


def _coerce_state(state: Any, *, operations: Optional[ScalarOperations] = None,
                  label: str = "fields") -> PhysicalFields:
    if isinstance(state, PhysicalFields):
        if operations is None:
            return state
        return PhysicalFields.from_fields(state.U, state.Q, Uz=state.Uz, Qz=state.Qz,
                                          operations=operations)
    return PhysicalFields.from_mapping(state, operations=operations, label=label)


def _provenance_call(guard: Any, role: str, value: Any) -> Any:
    """Call a user guard while accommodating the two common callback shapes."""

    if guard is None:
        return None
    if hasattr(guard, "assert_source"):
        result = guard.assert_source(role, value)
    elif hasattr(guard, "check"):
        result = guard.check(role, value)
    elif callable(guard):
        try:
            result = guard(role, value)
        except TypeError:
            result = guard(value, role=role)
    else:
        raise TypeError("source_provenance_guard must be callable or expose check/assert_source")
    if result is False:
        raise SourceProvenanceError("source provenance guard rejected " + role)
    return result


def _check_provenance(guard: Any, baseline: Any, correction: Any, Z: Any, delta: Any) -> dict[str, Any]:
    checked: dict[str, Any] = {}
    if guard is None:
        return {"provided": False, "checked_roles": (), "results": {}}
    for role, value in (("baseline", baseline), ("correction", correction),
                        ("Z", Z), ("delta", delta)):
        checked[role] = _provenance_call(guard, role, value)
    return {"provided": True, "checked_roles": tuple(checked), "results": checked}


def _common_operations(*jets: OrdinaryJet,
                       operations: Optional[ScalarOperations] = None) -> ScalarOperations:
    """Resolve one adapter without silently merging graph namespaces.

    ``None`` on a jet means that its scalar components are unbound and may
    inherit the unique explicit adapter in this operation.  In contrast,
    ``PYTHON_OPERATIONS`` is an explicit adapter identity, just like any
    caller-supplied ``ScalarOperations`` object; it is compatible only with
    itself (or with an unbound jet), never with a different graph adapter.
    """

    explicit = []
    for jet in jets:
        adapter = jet.operations
        if adapter is None:
            continue
        if operations is not None and adapter is not operations:
            raise TypeError(
                "ordinary jets use an explicit scalar operation adapter "
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
    return explicit[0] if explicit else PYTHON_OPERATIONS


def _normalize(jet: OrdinaryJet, operations: Optional[ScalarOperations]) -> OrdinaryJet:
    """Bind only an unbound jet, preserving explicit adapter identity.

    Explicit adapters are never relabeled.  ``None`` is the unbound marker;
    the singleton ``PYTHON_OPERATIONS`` remains a real explicit identity and
    therefore cannot be replaced by a graph adapter.  The private literal
    markers are the exceptions to ordinary row-preserving inheritance: bare
    numeric literals and numeric components of mapping/sequence rows are
    converted through the target adapter's ``constant`` callback, while
    opaque components remain unchanged.
    """

    if jet.operations is operations:
        return jet
    if isinstance(jet, _UnboundScalarLiteral) and jet.operations is None:
        if operations is None:
            return jet
        return OrdinaryJet(
            operations.constant(jet.value),
            operations.constant(0),
            operations.constant(0),
            operations,
        )
    if isinstance(jet, _UnboundMixedRow) and jet.operations is None:
        if operations is None:
            return jet
        rows = (jet.value, jet.dz, jet.dzz)
        bound = tuple(
            operations.constant(item) if is_numeric and item is not None else item
            for item, is_numeric in zip(rows, jet.numeric_components)
        )
        return OrdinaryJet(bound[0], bound[1], bound[2], operations)
    if jet.operations is not None:
        if operations is None:
            raise TypeError(
                "cannot normalize an explicit scalar operation adapter to an unbound jet"
            )
        raise TypeError(
            "ordinary jet uses a different scalar operation adapter; "
            "explicit adapters cannot be relabeled"
        )
    return OrdinaryJet(jet.value, jet.dz, jet.dzz, operations)


def _volume_terms(fields: PhysicalFields, Z: Any, delta: Any,
                  *, operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    z_raw = coerce_ordinary_jet(Z, operations=operations, name="Z")
    de_raw = coerce_ordinary_jet(delta, operations=operations, name="delta")
    ops = _common_operations(fields.U, fields.Q, fields.Uz, fields.Qz,
                             z_raw, de_raw, operations=operations)
    U, Q, Uz, Qz = (_normalize(j, ops) for j in (fields.U, fields.Q, fields.Uz, fields.Qz))
    z = _normalize(z_raw, ops)
    de = _normalize(de_raw, ops)
    one = OrdinaryJet.constant(ops.constant(1), ops)
    two = OrdinaryJet.constant(ops.constant(2), ops)
    four = OrdinaryJet.constant(ops.constant(4), ops)
    L = one - de * z * z
    D = one - z * z
    L2 = L * L
    return {
        "q_over_2L": -(Q / (two * L)),
        "z_qz_over_4L": z * (one - de) * Qz / (four * L),
        "minus_q_squared_over_4": -(Q * Q / four),
        "D_u_qz_over_2L": D * U * Qz / (two * L),
        "z2_one_minus_delta_u2_over_L2": z * z * (one - de) * U * U / L2,
        "z_D_u_uz_over_L2": z * D * U * Uz / L2,
    }


def volume_density(state: Any, Z: Any, delta: Any, *,
                   operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Evaluate the source volume density v on one correlated physical state."""

    fields = _coerce_state(state, operations=operations, label="state")
    terms = _volume_terms(fields, Z, delta, operations=operations)
    ops = _common_operations(fields.U, fields.Q, operations=operations)
    out = _numeric_constant_jet(0, ops)
    for term in terms.values():
        out = out + term
    return out


def volume_density_value(U: Any, Q: Any, Uz: Any, Qz: Any,
                         Z: Any, delta: Any, *,
                         operations: Optional[ScalarOperations] = None) -> Any:
    """Evaluate only v itself, retaining the exact scalar operation graph."""

    fields = PhysicalFields.from_fields(U, Q, Uz=Uz, Qz=Qz, operations=operations)
    return volume_density(fields, Z, delta, operations=operations).value


def _endpoint_terms(fields: PhysicalFields, Z: Any, delta: Any, *,
                    operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    """Return exact source boundary-density terms for F."""

    z_raw = coerce_ordinary_jet(Z, operations=operations, name="Z")
    de_raw = coerce_ordinary_jet(delta, operations=operations, name="delta")
    ops = _common_operations(fields.U, fields.Q, fields.Uz, fields.Qz,
                             z_raw, de_raw, operations=operations)
    U, Q = (_normalize(j, ops) for j in (fields.U, fields.Q))
    z = _normalize(z_raw, ops)
    de = _normalize(de_raw, ops)
    one = OrdinaryJet.constant(ops.constant(1), ops)
    two = OrdinaryJet.constant(ops.constant(2), ops)
    four = OrdinaryJet.constant(ops.constant(4), ops)
    L = one - de * z * z
    L2 = L * L
    return {
        "q_over_2L": Q / (two * L),
        "q_squared_over_4": Q * Q / four,
        "minus_z2_u_squared_over_L2": -(z * z * U * U / L2),
    }


def endpoint_density(state: Any, Z: Any, delta: Any, *,
                     operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Evaluate the source endpoint density ``F``.

    The physical m15 endpoint uses ``R**2 * F - R * Q_y``.  The preceding
    pressure inlet term ``-R*Q_y`` is intentionally kept explicit by the
    caller; this function supplies only the source ``F`` polynomial.
    """

    fields = _coerce_state(state, operations=operations, label="state")
    terms = _endpoint_terms(fields, Z, delta, operations=operations)
    ops = _common_operations(fields.U, fields.Q, operations=operations)
    out = _numeric_constant_jet(0, ops)
    for term in terms.values():
        out = out + term
    return out


def endpoint_density_value(U: Any, Q: Any, Z: Any, delta: Any, *,
                           Uz: Any = None, Qz: Any = None,
                           operations: Optional[ScalarOperations] = None) -> Any:
    """Evaluate only F itself, retaining the exact scalar operation graph."""

    fields = PhysicalFields.from_fields(U, Q, Uz=Uz, Qz=Qz, operations=operations)
    return endpoint_density(fields, Z, delta, operations=operations).value


def _difference_terms(base: PhysicalFields, correction: PhysicalFields,
                      Z: Any, delta: Any, *,
                      operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    z_raw = coerce_ordinary_jet(Z, operations=operations, name="Z")
    de_raw = coerce_ordinary_jet(delta, operations=operations, name="delta")
    ops = _common_operations(base.U, base.Q, base.Uz, base.Qz,
                             correction.U, correction.Q, correction.Uz, correction.Qz,
                             z_raw, de_raw,
                             operations=operations)
    B = PhysicalFields(*(_normalize(j, ops) for j in (base.U, base.Q, base.Uz, base.Qz)))
    C = PhysicalFields(*(_normalize(j, ops) for j in (correction.U, correction.Q,
                                                      correction.Uz, correction.Qz)))
    z = _normalize(z_raw, ops)
    de = _normalize(de_raw, ops)
    one = OrdinaryJet.constant(ops.constant(1), ops)
    two = OrdinaryJet.constant(ops.constant(2), ops)
    four = OrdinaryJet.constant(ops.constant(4), ops)
    L = one - de * z * z
    D = one - z * z
    L2 = L * L
    dU, dQ, dUz, dQz = C.U, C.Q, C.Uz, C.Qz
    return {
        "linear_q_over_2L": -(dQ / (two * L)),
        "linear_z_qz_over_4L": z * (one - de) * dQz / (four * L),
        "quadratic_q_squared_over_4": -((two * B.Q * dQ + dQ * dQ) / four),
        "correlated_D_u_qz_over_2L": D * (B.U * dQz + B.Qz * dU + dU * dQz) / (two * L),
        "correlated_z2_u_squared_over_L2": z * z * (one - de) *
            (two * B.U * dU + dU * dU) / L2,
        "correlated_z_D_u_uz_over_L2": z * D *
            (B.U * dUz + B.Uz * dU + dU * dUz) / L2,
    }


def _endpoint_difference_terms(base: PhysicalFields, correction: PhysicalFields,
                               Z: Any, delta: Any, *,
                               operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    """Return exact factored terms for ΔF, preserving both quadratic terms."""

    z_raw = coerce_ordinary_jet(Z, operations=operations, name="Z")
    de_raw = coerce_ordinary_jet(delta, operations=operations, name="delta")
    ops = _common_operations(base.U, base.Q, base.Uz, base.Qz,
                             correction.U, correction.Q, correction.Uz, correction.Qz,
                             z_raw, de_raw, operations=operations)
    B = PhysicalFields(*(_normalize(j, ops) for j in (base.U, base.Q, base.Uz, base.Qz)))
    C = PhysicalFields(*(_normalize(j, ops) for j in (correction.U, correction.Q,
                                                      correction.Uz, correction.Qz)))
    z = _normalize(z_raw, ops)
    de = _normalize(de_raw, ops)
    one = OrdinaryJet.constant(ops.constant(1), ops)
    two = OrdinaryJet.constant(ops.constant(2), ops)
    four = OrdinaryJet.constant(ops.constant(4), ops)
    L = one - de * z * z
    L2 = L * L
    dU, dQ = C.U, C.Q
    return {
        "linear_q_over_2L": dQ / (two * L),
        "quadratic_q_squared_over_4": (two * B.Q * dQ + dQ * dQ) / four,
        "quadratic_minus_z2_u_squared_over_L2":
            -(z * z * (two * B.U * dU + dU * dU) / L2),
    }


def factored_endpoint_density_difference(baseline: Any, correction: Any,
                                         Z: Any, delta: Any, *,
                                         source_provenance_guard: Any = None,
                                         operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Return exact factored ``ΔF = F(base + correction) - F(base)``."""

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    terms = _endpoint_difference_terms(base, corr, Z, delta, operations=operations)
    ops = _common_operations(base.U, base.Q, corr.U, corr.Q,
                             operations=operations)
    out = _numeric_constant_jet(0, ops)
    for term in terms.values():
        out = out + term
    return out


def full_endpoint_density_difference(baseline: Any, correction: Any,
                                     Z: Any, delta: Any, *,
                                     source_provenance_guard: Any = None,
                                     operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Evaluate ΔF independently from the full endpoint polynomial."""

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    return endpoint_density(base + corr, Z, delta, operations=operations) - endpoint_density(
        base, Z, delta, operations=operations)


def endpoint_density_difference(baseline: Any, correction: Any,
                                Z: Any, delta: Any, *,
                                source_provenance_guard: Any = None,
                                operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Public signed endpoint-difference entry point; returns factored ΔF."""

    return factored_endpoint_density_difference(
        baseline, correction, Z, delta,
        source_provenance_guard=source_provenance_guard,
        operations=operations)


def endpoint_difference_terms(baseline: Any, correction: Any, Z: Any, delta: Any, *,
                              source_provenance_guard: Any = None,
                              operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    """Return individual exact ΔF terms for source-equivalence inspection."""

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    return _endpoint_difference_terms(base, corr, Z, delta, operations=operations)


def factored_volume_density_difference(baseline: Any, correction: Any,
                                       Z: Any, delta: Any, *,
                                       source_provenance_guard: Any = None,
                                       operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Return the exact factored Δv = v(base + correction) - v(base).

    The products U*dQz, Qz*dU, dU*dQz, 2*U*dU, dU*dU, U*dUz, Uz*dU and
    dU*dUz are retained separately in the source algebra.  This matters when
    the baseline is nonzero or a supposedly fast-mean correction is zero.
    """

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    terms = _difference_terms(base, corr, Z, delta, operations=operations)
    ops = _common_operations(base.U, base.Q, corr.U, corr.Q,
                             operations=operations)
    out = _numeric_constant_jet(0, ops)
    for term in terms.values():
        out = out + term
    return out


def full_volume_density_difference(baseline: Any, correction: Any,
                                   Z: Any, delta: Any, *,
                                   source_provenance_guard: Any = None,
                                   operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Evaluate Δv independently from the full source density expression."""

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    total = base + corr
    return volume_density(total, Z, delta, operations=operations) - volume_density(
        base, Z, delta, operations=operations)


def volume_density_difference(baseline: Any, correction: Any,
                              Z: Any, delta: Any, *,
                              source_provenance_guard: Any = None,
                              operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    """Public signed-difference entry point; returns the factored C1 jet."""

    return factored_volume_density_difference(
        baseline, correction, Z, delta,
        source_provenance_guard=source_provenance_guard,
        operations=operations)


def volume_difference_terms(baseline: Any, correction: Any, Z: Any, delta: Any, *,
                            source_provenance_guard: Any = None,
                            operations: Optional[ScalarOperations] = None) -> dict[str, OrdinaryJet]:
    """Return individual exact Δv terms for source-equivalence inspection."""

    _check_provenance(source_provenance_guard, baseline, correction, Z, delta)
    base = _coerce_state(baseline, operations=operations, label="baseline")
    corr = _coerce_state(correction, operations=operations, label="correction")
    return _difference_terms(base, corr, Z, delta, operations=operations)


SOURCE_EVIDENCE = {
    "volume_source": {
        "path": "src/openai_ns_reconstruction/original_macro_m15.py",
        "lines": {"half_S1_coefficients": [47, 53], "densities": [56, 59],
                  "correction_densities": [274, 284], "finite_errors": [286, 304]},
        "ast_targets": ["half_S1_coefficients", "densities",
                        "OriginalMacroM15.correction_densities"],
        "source_variables": {"L": "1-delta*Z^2", "D": "1-Z^2",
                             "Uz": "dU/dZ", "Qz": "dQ/dZ"},
        "source_formula": "-Q/(2L)+Z*(1-delta)*Qz/(4L)-Q^2/4+D*U*Qz/(2L)+Z^2*(1-delta)*U^2/L^2+Z*D*U*Uz/L^2",
        "caps_are_values": False,
    },
    "endpoint_density_source": {
        "path": "src/openai_ns_reconstruction/original_macro_m15.py",
        "lines": {"half_S1_coefficients": [47, 53]},
        "ast_targets": ["half_S1_coefficients endpoint coefficients bc"],
        "source_formula": "F=Q/(2L)+Q^2/4-Z^2*U^2/L^2",
        "m15_endpoint_combination": "R^2*F-R*Q_y",
        "preceding_pressure_inlet_is_explicit": True,
    },
    "post_axial_quiet_sources": {
        "paths": ["src/openai_ns_reconstruction/original_post_axial_m15.py",
                   "src/openai_ns_reconstruction/original_compact_quiet_m15.py"],
        "lines": {"post_axial_evaluate": [43, 93],
                  "compact_quiet_finite_errors": [45, 60],
                  "compact_quiet_evaluate": [62, 109]},
        "ast_targets": ["OriginalPostAxialM15.evaluate",
                        "OriginalCompactQuietM15.finite_errors",
                        "OriginalCompactQuietM15.evaluate"],
        "finite_N_error_logs_are_enclosures": True,
        "signed_centers_selected": False,
    },
    "pressure_ftc_source": {
        "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
        "lines": {"physical_pressure_fields": [193, 202],
                  "reduced_pressure_S": [205, 213]},
        "ast_targets": ["physical_pressure_fields", "reduced_pressure_S"],
        "source_formula": "S=div(Q+(1-delta)*Z*Qz/2+Qy)+Q*(Q/2+Qy)+div(U*(-2*Z*Q+D*Qz-2*Z*Qy))",
        "same_L": "1-delta*Z^2",
        "verified_identity": "S/2=d_y F+2*F+v under L*(Q+Qy)=Z*(1+delta)*U+2*Z*Uy-D*Uz",
    },
    "compact_actual_C5": {
        "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py",
        "lines": {"profiles": [91, 101], "density": [104, 107],
                  "background_history_modified_fields": [196, 215],
                  "physical_and_pressure": [216, 218]},
        "ast_targets": ["profiles.radial_y_rows", "density", "self.history",
                        "self.modified_fields", "self.radial", "self.physical", "self.S",
                        "self.pressure"],
        "numeric_signed_centers_available": False,
    },
    "rm_actual_C5": {
        "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py",
        "lines": {"pressure_continuation": [65, 78],
                  "finite_history_modified_fields": [95, 111],
                  "radial_physical_pressure": [112, 119]},
        "ast_targets": ["five_signed_density_increments", "self.history",
                        "self.chart.modified_fields", "self.radial", "self.physical", "self.S",
                        "self.pressure"],
        "numeric_signed_centers_available": False,
    },
    "outer_pressure_geometry": {
        "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
        "lines": {"pressure_geometry": [79, 90], "chart_physical": [211, 215]},
        "ast_targets": ["pressure_geometry", "chart.physical", "chart.S", "chart.pressure"],
        "physical_m15_sign": "+1/2 cumulative S1 minus boundary R*Q_y difference",
    },
}


def source_equivalence_evidence() -> dict[str, Any]:
    """Return source paths and AST targets used by this kernel."""

    return SOURCE_EVIDENCE.copy()


def remaining_numerical_blocked_handles() -> dict[str, Any]:
    """Describe actual graph handles still lacking numeric signed centers."""

    return {
        "numeric_signed_centers_available": False,
        "blocked_reason": "Current C5 sources retain exact graph functions and finite-N signed increments, but no accepted numeric signed center values.",
        "compact": {
            "baseline_handles": ["background.fields", "physical.U", "physical.Q", "physical.Q_y"],
            "correction_handles": ["history", "modified_fields - background.fields", "physical(modified) - physical(background)"],
            "correlated_delta_rule": "Build delta U and delta Q from the same modified physical graph; do not independently box history components.",
            "profiles": ["profiles.radial_y_rows", "density", "history", "modified_fields", "radial", "physical", "S", "pressure"],
        },
        "rm_patch": {
            "baseline_handles": ["chart.source_chart.fields", "physical.U", "physical.Q", "physical.Q_y"],
            "correction_handles": ["five_signed_density_increments", "history", "chart.modified_fields", "physical(modified) - physical(background)"],
            "correlated_delta_rule": "Use the same finite-N signed increment history for U and Q before forming delta v.",
            "profiles": ["finite_N", "history", "chart.modified_fields", "radial", "physical", "S", "pressure"],
        },
        "macro_caps": {
            "path": "src/openai_ns_reconstruction/original_macro_m15.py",
            "field": "correction_densities / finite_errors",
            "status": "absolute error enclosures only; never a signed source value",
        },
    }


def _sum_jets(terms: dict[str, OrdinaryJet], operations: Optional[ScalarOperations] = None) -> OrdinaryJet:
    if not terms:
        return _numeric_constant_jet(0, operations)
    first = next(iter(terms.values()))
    ops = _common_operations(first, operations=operations)
    out = _numeric_constant_jet(0, ops)
    for term in terms.values():
        out = out + term
    return out


def run_exact_controls() -> dict[str, Any]:
    """Run bounded exact rational/SymPy controls for the signed kernel.

    These controls intentionally evaluate independently supplied rational or
    symbolic fields only; no archived cap or source constructor is executed.
    """

    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover - repository normally has SymPy
        return {"passed": False, "sympy_available": False, "error": str(exc)}

    z, de = sp.symbols("Z delta")
    u, q, uz, qz = sp.symbols("U Q Uz Qz")
    du, dq, duz, dqz = sp.symbols("dU dQ dUz dQz")
    zjet = ordinary_jet(z, 1, 0)
    base = physical_fields(ordinary_jet(u, uz, sp.Symbol("Uzz")),
                           ordinary_jet(q, qz, sp.Symbol("Qzz")))
    corr = physical_fields(ordinary_jet(du, duz, sp.Symbol("dUzz")),
                           ordinary_jet(dq, dqz, sp.Symbol("dQzz")))
    full = full_volume_density_difference(base, corr, zjet, de)
    factored = factored_volume_density_difference(base, corr, zjet, de)
    value_identity = sp.simplify(full.value - factored.value) == 0
    derivative_identity = sp.simplify(full.dz - factored.dz) == 0
    endpoint_full = full_endpoint_density_difference(base, corr, zjet, de)
    endpoint_factored = factored_endpoint_density_difference(base, corr, zjet, de)
    endpoint_value_identity = sp.simplify(endpoint_full.value - endpoint_factored.value) == 0
    endpoint_derivative_identity = sp.simplify(endpoint_full.dz - endpoint_factored.dz) == 0
    # A separate polynomial fixture supplies the ordinary derivatives of U and
    # Q.  Differentiating the symbolic-name fixture above would incorrectly
    # treat Uz and Qz as constants and would therefore test a different graph.
    U_poly = z**2 + 3*z + 2
    Q_poly = z**3 - z + 1
    dU_poly = 2*z - 1
    dQ_poly = z**2 + 2
    poly_base = physical_fields(ordinary_jet(U_poly, sp.diff(U_poly, z),
                                              sp.diff(U_poly, z, 2)),
                                ordinary_jet(Q_poly, sp.diff(Q_poly, z),
                                              sp.diff(Q_poly, z, 2)))
    poly_corr = physical_fields(ordinary_jet(dU_poly, sp.diff(dU_poly, z),
                                              sp.diff(dU_poly, z, 2)),
                                ordinary_jet(dQ_poly, sp.diff(dQ_poly, z),
                                              sp.diff(dQ_poly, z, 2)))
    poly_full = full_volume_density_difference(poly_base, poly_corr, zjet, de)
    direct_derivative = sp.simplify(poly_full.dz - sp.diff(poly_full.value, z)) == 0

    # Independent source FTC check.  Qy is eliminated using the actual source
    # constraint L*(Q + Qy) = Z*(1 + delta)*U + 2*Z*Uy - D*Uz.
    Uy = sp.symbols("Uy")
    Ls = 1 - de * z**2
    Ds = 1 - z**2
    Qy = (z * (1 + de) * u + 2 * z * Uy - Ds * uz) / Ls - q
    F_source = q / (2 * Ls) + q**2 / 4 - z**2 * u**2 / Ls**2
    v_source = (-q / (2 * Ls) + z * (1 - de) * qz / (4 * Ls) - q**2 / 4
                + Ds * u * qz / (2 * Ls) + z**2 * (1 - de) * u**2 / Ls**2
                + z * Ds * u * uz / Ls**2)
    dy_F_source = Qy / (2 * Ls) + q * Qy / 2 - 2 * z**2 * u * Uy / Ls**2
    S_source = ((q + (1 - de) * z * qz / 2 + Qy) / Ls
                + q * (q / 2 + Qy)
                + u * (-2 * z * q + Ds * qz - 2 * z * Qy) / Ls)
    pressure_ftc_identity = sp.factor(S_source / 2 - (dy_F_source + 2 * F_source + v_source)) == 0

    # Exact rational sign and quadratic control.  The fast-mean correction is
    # zero (dU = dUz = 0), while p0 / baseline U is nonzero and dQ^2 survives.
    rational = dict(Z=Fraction(1, 3), delta=Fraction(1, 5), U=Fraction(2),
                    Q=Fraction(1, 7), Uz=Fraction(3, 5), Qz=Fraction(4, 9),
                    dU=Fraction(0), dQ=Fraction(2, 11), dUz=Fraction(0),
                    dQz=Fraction(0))
    rb = physical_fields((rational["U"], rational["Uz"], Fraction(0)),
                         (rational["Q"], rational["Qz"], Fraction(0)))
    rc = physical_fields((rational["dU"], rational["dUz"], Fraction(0)),
                         (rational["dQ"], rational["dQz"], Fraction(0)))
    rfull = full_volume_density_difference(rb, rc, rational["Z"], rational["delta"])
    rterms = volume_difference_terms(rb, rc, rational["Z"], rational["delta"])
    rfactored = _sum_jets(rterms).value
    rational_identity = rfull.value == rfactored
    linear = rterms["linear_q_over_2L"] + rterms["linear_z_qz_over_4L"]
    linear = linear + rterms["correlated_D_u_qz_over_2L"]
    linear = linear + rterms["correlated_z2_u_squared_over_L2"]
    linear = linear + rterms["correlated_z_D_u_uz_over_L2"]
    linear_only_gap = rfull.value - linear.value
    nonzero_quadratic_p0 = rational["U"] != 0 and rterms["quadratic_q_squared_over_4"].value != 0
    zero_fast_mean_counterexample = (rational["dU"] == 0 and rational["dUz"] == 0 and
                                     linear_only_gap != 0)

    # Provenance callback control: rejection is observable and no value is
    # accepted as a source function without the caller's identity decision.
    roles: list[str] = []
    guarded = lambda role, value: roles.append(role) or True
    factored_volume_density_difference(rb, rc, rational["Z"], rational["delta"],
                                       source_provenance_guard=guarded)
    provenance_control = roles == ["baseline", "correction", "Z", "delta"]
    rejected = False
    try:
        factored_volume_density_difference(rb, rc, rational["Z"], rational["delta"],
                                           source_provenance_guard=lambda role, value: False)
    except SourceProvenanceError:
        rejected = True

    passed = all((value_identity, derivative_identity, endpoint_value_identity,
                  endpoint_derivative_identity, direct_derivative,
                  pressure_ftc_identity,
                  rational_identity, nonzero_quadratic_p0,
                  zero_fast_mean_counterexample, provenance_control, rejected))
    return {
        "passed": passed,
        "sympy_available": True,
        "value_identity": value_identity,
        "derivative_identity": derivative_identity,
        "endpoint_value_identity": endpoint_value_identity,
        "endpoint_derivative_identity": endpoint_derivative_identity,
        "direct_derivative_identity": direct_derivative,
        "pressure_ftc_identity": pressure_ftc_identity,
        "rational_identity": rational_identity,
        "nonzero_quadratic_p0": nonzero_quadratic_p0,
        "zero_fast_mean_counterexample": zero_fast_mean_counterexample,
        "quadratic_gap": linear_only_gap,
        "provenance_guard_control": provenance_control and rejected,
        "source_caps_selected": False,
        "numeric_signed_centers_used": False,
    }


# Short aliases used by parent integration code.
signed_density_difference = volume_density_difference
factored_signed_density_difference = factored_volume_density_difference
signed_endpoint_density_difference = endpoint_density_difference
factored_signed_endpoint_density_difference = factored_endpoint_density_difference


__all__ = [
    "JetFormatError", "SourceProvenanceError", "ScalarOperations",
    "OrdinaryJet", "PhysicalFields", "ordinary_jet", "coerce_ordinary_jet",
    "physical_fields", "volume_density", "volume_density_value",
    "endpoint_density", "endpoint_density_value",
    "factored_volume_density_difference", "full_volume_density_difference",
    "volume_density_difference", "signed_density_difference",
    "factored_signed_density_difference", "volume_difference_terms",
    "factored_endpoint_density_difference", "full_endpoint_density_difference",
    "endpoint_density_difference", "signed_endpoint_density_difference",
    "factored_signed_endpoint_density_difference", "endpoint_difference_terms",
    "source_equivalence_evidence", "remaining_numerical_blocked_handles",
    "run_exact_controls", "SOURCE_EVIDENCE",
]
