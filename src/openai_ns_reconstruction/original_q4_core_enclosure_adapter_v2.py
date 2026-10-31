"""Caller-owned adapter for the canonical Q4 core component enclosure.

The original Q4 evaluator calls an oracle as ``component(record, rho, Z)``.
The canonical core source, however, returns a *record* containing a finite
polynomial enclosure, an infinite radial tail bound, and the complete
directed enclosure.  This module is the narrow boundary between those two
interfaces.  It never imports a producer, invokes an old provider, evaluates
a source graph, takes an interval midpoint, or turns a tail/cap into a value.

Only a caller-supplied provider is used, and it is called lazily after all
record-side admission checks pass.  A strict caller-supplied interval
validator must accept the provider's ``enclosure`` object.  That callback is
treated as caller evidence only; the receipt deliberately remains
uncertified and does not call the callback a proof.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Number
from typing import Any, Callable


SCHEMA = "original_q4_core_enclosure_adapter_v2"
SUPPORTED_QUANTITIES = ("Phi", "Psi", "mean_Psi")

_MISSING = object()
_PARTIAL_FIELDS = (
    "finite_polynomial_enclosure",
    "infinite_radial_tail_upper",
    "enclosure",
)
_ORDINARY_FLAG_NAMES = (
    "all_rows_are_ordinary_derivatives",
    "genuine_ordinary_derivatives",
)
_NONLINEAR_FLAG_NAMES = (
    "actual_twenty_term_nonlinear_solution_enclosed",
    "actual_twenty_term_solution_enclosed",
)
_CENTER_FLAG_NAMES = (
    "numerical_center_not_selected",
    "center_not_selected",
)


class EnclosureAdapterError(ValueError):
    """Fail-closed error raised for a malformed or unadmitted component."""

    def __init__(self, code: str, detail: Any = None):
        self.code = code
        self.detail = detail
        message = code if detail is None else f"{code}: {detail}"
        super().__init__(message)


@dataclass(frozen=True, slots=True)
class ComponentEnclosure:
    """A tagged enclosure; deliberately has no scalar arithmetic or conversion.

    An interval-aware evaluator must explicitly consume this type and preserve
    its enclosure semantics through arithmetic and quadrature.  The current
    scalar evaluator alone does not provide that contract.
    """
    interval: Any
    kind: str = "enclosure"
    actual_value_selected: bool = False
    certified_by_adapter: bool = False

    def __float__(self) -> float:
        raise TypeError("an enclosure is not a selected point value")

    def __bool__(self) -> bool:
        raise TypeError("an enclosure has no scalar truth value")


def _same(left: Any, right: Any) -> bool:
    """Compare caller metadata without coercing references or containers."""

    if left is right:
        return True
    try:
        result = left == right
    except Exception:
        return False
    if type(result) is bool:
        return result
    try:
        return bool(result)
    except Exception:
        return False


def _nonempty(value: Any) -> bool:
    if value is _MISSING or value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (bytes, bytearray)):
        return bool(value)
    if isinstance(value, Mapping):
        return bool(value)
    return True


def _integer(value: Any, label: str) -> int:
    if type(value) is not int:
        raise EnclosureAdapterError("invalid_integer", label)
    return value


def _first_present(record: Mapping[str, Any], names: tuple[str, ...]) -> Any:
    """Read a field from the emitted row and its opaque source row.

    ``original_q4_scoped_ir_evaluator_v1._source_record`` preserves the
    source row and selected direct node fields.  Looking at those two explicit
    containers is sufficient and avoids treating an integer graph reference
    or an arbitrary object as a source lookup.
    """

    containers: list[Mapping[str, Any]] = [record]
    nested = record.get("source_row", _MISSING)
    if isinstance(nested, Mapping) and nested is not record:
        containers.append(nested)
    found: list[Any] = []
    for container in containers:
        for name in names:
            if name in container:
                found.append(container[name])
    if not found:
        return _MISSING
    value = found[0]
    if any(not _same(value, item) for item in found[1:]):
        raise EnclosureAdapterError("record_field_mismatch", "/".join(names))
    return value


def _required_field(record: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    value = _first_present(record, names)
    if value is _MISSING:
        raise EnclosureAdapterError("missing_" + label)
    return value


def _optional_field(record: Mapping[str, Any], names: tuple[str, ...]) -> Any:
    return _first_present(record, names)


def _require_exact_flag(result: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    values = [result[name] for name in names if name in result]
    if not values or any(type(value) is not bool or value is not True for value in values):
        raise EnclosureAdapterError("missing_or_false_" + label)


def _descriptor_from_record(record: Mapping[str, Any]) -> Any:
    # v4 may expose ``provider`` as a mapping of all provider roles while the
    # canonical source row separately has the specific ``source_provider``.
    # They are different roles and must not be compared as aliases.
    exact = _first_present(record, ("source_provider",))
    if exact is not _MISSING:
        return exact
    return _required_field(
        record,
        (
            "provider_descriptor",
            "source_provider_descriptor",
            "source_provider_description",
            "provider_description",
            "provider",
        ),
        "source_provider_descriptor",
    )


def _evidence_from_record(record: Mapping[str, Any]) -> Any:
    return _first_present(
        record,
        (
            "source_evidence",
            "source_provider_evidence",
            "provider_evidence",
        ),
    )


def _resolve_alias(
    primary: Any,
    aliases: tuple[Any, ...],
    *,
    label: str,
    missing: Any = _MISSING,
) -> Any:
    # ``None`` is the public default for optional callable arguments.  When a
    # named alias is supplied, that default is absent rather than a conflicting
    # second value.  Explicit non-None values are still compared exactly.
    values = [value for value in (primary,) + aliases if value is not _MISSING and value is not missing]
    if not values:
        return missing
    first = values[0]
    if any(not _same(first, value) for value in values[1:]):
        raise EnclosureAdapterError("constructor_argument_mismatch", label)
    return first


class CoreEnclosureAdapter:
    """Adapt a lazy canonical ``component_enclosure`` provider.

    Parameters are intentionally caller-owned.  ``provider`` must implement
    the canonical call shape::

        provider(quantity, rho, Z, radial_order=i, axial_order=k,
                 degree=degree, depth=depth)

    ``provider_descriptor`` is the exact descriptor expected in each source
    record.  ``source_evidence`` is optional additional exact metadata; when
    supplied, the record must carry the same evidence.  A validator is
    required at component-call time and must return the boolean ``True`` for
    the complete enclosure object.  It is not a proof and is reported as an
    unverified caller claim.
    """

    def __init__(
        self,
        provider: Callable[..., Any] | None = None,
        source_family: Any = _MISSING,
        provider_descriptor: Any = _MISSING,
        degree: int = 6,
        depth: int = 3,
        interval_validator: Callable[[Any], Any] | None = None,
        source_evidence: Any = _MISSING,
        *,
        provider_callable: Any = _MISSING,
        expected_provider_descriptor: Any = _MISSING,
        expected_source_provider: Any = _MISSING,
        source_provider: Any = _MISSING,
        source_provider_description: Any = _MISSING,
        expected_source_provider_description: Any = _MISSING,
        expected_source_evidence: Any = _MISSING,
        validate_interval: Any = _MISSING,
        strict_interval_validator: Any = _MISSING,
        expected_degree: Any = _MISSING,
        expected_depth: Any = _MISSING,
    ) -> None:
        self.provider = _resolve_alias(
            provider,
            (provider_callable,),
            label="provider",
            missing=None,
        )
        self.source_family = source_family
        self.expected_provider_descriptor = _resolve_alias(
            provider_descriptor,
            (
                expected_provider_descriptor,
                expected_source_provider,
                source_provider,
                source_provider_description,
                expected_source_provider_description,
            ),
            label="provider_descriptor",
        )
        self.expected_source_evidence = _resolve_alias(
            source_evidence,
            (expected_source_evidence,),
            label="source_evidence",
        )
        self.degree = _resolve_alias(
            degree,
            (expected_degree,),
            label="degree",
        )
        self.depth = _resolve_alias(
            depth,
            (expected_depth,),
            label="depth",
        )
        self.interval_validator = _resolve_alias(
            interval_validator,
            (validate_interval, strict_interval_validator),
            label="interval_validator",
            missing=None,
        )
        self.last_receipt: dict[str, Any] = self._receipt_base()

    def _receipt_base(self) -> dict[str, Any]:
        descriptor = (
            None
            if self.expected_provider_descriptor is _MISSING
            else self.expected_provider_descriptor
        )
        evidence = None if self.expected_source_evidence is _MISSING else self.expected_source_evidence
        family = None if self.source_family is _MISSING else self.source_family
        return {
            "schema": SCHEMA,
            "kind": "enclosure",
            "status": "uninitialized",
            "global_status": "uncertified",
            "certified": False,
            "source_binding": {
                "source_family": family,
                "provider_descriptor": descriptor,
                "source_evidence": evidence,
                "independently_verified": False,
            },
            "source_binding_independently_verified": False,
            "errors": "not independently verified",
            "errors_independently_verified": False,
            "caller_claims": {
                "source_family_supplied": self.source_family is not _MISSING,
                "provider_descriptor_supplied": self.expected_provider_descriptor is not _MISSING,
                "source_evidence_supplied": self.expected_source_evidence is not _MISSING,
                "provider_callable_injected": callable(self.provider),
                "degree": self.degree,
                "depth": self.depth,
                "interval_validator_supplied": callable(self.interval_validator),
                "interval_validator_is_proof": False,
            },
            "actual_value_selected": False,
            "enclosure_is_point_value": False,
            "finite_polynomial_is_output": False,
            "tail_cap_is_output": False,
            "midpoint_selected": False,
        }

    def _reject(self, error: EnclosureAdapterError, request: Mapping[str, Any] | None = None) -> None:
        receipt = self._receipt_base()
        receipt["status"] = "rejected"
        if request is not None:
            receipt["request"] = dict(request)
        receipt["error"] = {"code": error.code}
        if error.detail is not None:
            receipt["error"]["detail"] = error.detail
        self.last_receipt = receipt

    def _accept(
        self,
        request: Mapping[str, Any],
        result: Mapping[str, Any],
    ) -> None:
        receipt = self._receipt_base()
        receipt["status"] = "accepted_enclosure"
        receipt["request"] = dict(request)
        receipt["returned"] = {
            "kind": "enclosure",
            "quantity": request["quantity"],
            "radial_order": request["radial_order"],
            "axial_order": request["axial_order"],
            "degree": request["degree"],
            "depth": request["depth"],
            "complete_record_received": True,
            "interval_validation_callback_returned_true": True,
            "interval_validation_is_proof": False,
            "actual_value_selected": False,
            "result_metadata_flags_retained": True,
        }
        receipt["provider_result_metadata"] = {
            name: result[name]
            for name in (
                "all_rows_are_ordinary_derivatives",
                "genuine_ordinary_derivatives",
                "actual_twenty_term_nonlinear_solution_enclosed",
                "actual_twenty_term_solution_enclosed",
                "numerical_center_not_selected",
                "center_not_selected",
            )
            if name in result
        }
        self.last_receipt = receipt

    def _read_request(self, record: Any) -> dict[str, Any]:
        if not isinstance(record, Mapping):
            raise EnclosureAdapterError("record_not_mapping")
        quantity = _required_field(record, ("quantity", "source_quantity"), "quantity")
        if not isinstance(quantity, str) or quantity not in SUPPORTED_QUANTITIES:
            raise EnclosureAdapterError("unsupported_quantity", quantity)
        radial_order = _required_field(
            record,
            ("radial_order", "radial_derivative_order", "radial_index", "i"),
            "radial_order",
        )
        axial_order = _required_field(
            record,
            (
                "ordinary_slow_Z_derivative_order",
                "ordinary_Z_derivative_order",
                "axial_order",
                "Z_derivative_order",
                "Z_derivative",
                "k",
            ),
            "ordinary_slow_Z_derivative_order",
        )
        radial_order = _integer(radial_order, "radial_order")
        axial_order = _integer(axial_order, "ordinary_slow_Z_derivative_order")
        if radial_order < 0 or axial_order < 0:
            raise EnclosureAdapterError("negative_derivative_order")
        if radial_order + axial_order > 5:
            raise EnclosureAdapterError("derivative_order_sum_exceeds_five")
        if type(self.degree) is not int or self.degree < 6:
            raise EnclosureAdapterError("degree_below_six")
        if type(self.depth) is not int or self.depth < 0:
            raise EnclosureAdapterError("invalid_depth")
        if axial_order > self.depth:
            raise EnclosureAdapterError("axial_order_exceeds_depth")
        if self.source_family is _MISSING or not _nonempty(self.source_family):
            raise EnclosureAdapterError("missing_source_family_binding")
        record_family = _required_field(record, ("source_family", "family"), "source_family")
        if not _same(record_family, self.source_family):
            raise EnclosureAdapterError("source_family_mismatch")
        if self.expected_provider_descriptor is _MISSING or not _nonempty(self.expected_provider_descriptor):
            raise EnclosureAdapterError("missing_expected_provider_descriptor")
        descriptor = _descriptor_from_record(record)
        if not _same(descriptor, self.expected_provider_descriptor):
            raise EnclosureAdapterError("provider_descriptor_mismatch")
        if self.expected_source_evidence is not _MISSING:
            evidence = _evidence_from_record(record)
            if evidence is _MISSING or not _same(evidence, self.expected_source_evidence):
                raise EnclosureAdapterError("source_evidence_mismatch")

        record_degree = _optional_field(record, ("degree", "radial_degree"))
        if record_degree is not _MISSING and record_degree != self.degree:
            raise EnclosureAdapterError("record_degree_mismatch")
        record_depth = _optional_field(record, ("depth", "axial_depth"))
        if record_depth is not _MISSING and record_depth != self.depth:
            raise EnclosureAdapterError("record_depth_mismatch")
        return {
            "quantity": quantity,
            "radial_order": radial_order,
            "axial_order": axial_order,
            "degree": self.degree,
            "depth": self.depth,
        }

    def _validate_result(
        self,
        result: Any,
        request: Mapping[str, Any],
    ) -> Any:
        if not isinstance(result, Mapping):
            raise EnclosureAdapterError("provider_result_not_mapping")
        for field in _PARTIAL_FIELDS:
            if field not in result or result[field] is None:
                raise EnclosureAdapterError("missing_provider_field", field)
        finite = result["finite_polynomial_enclosure"]
        tail = result["infinite_radial_tail_upper"]
        enclosure = result["enclosure"]
        if enclosure is finite or enclosure is tail:
            raise EnclosureAdapterError("enclosure_aliases_partial_record")
        if isinstance(enclosure, (str, bytes, bytearray, bool)) or isinstance(enclosure, Number):
            raise EnclosureAdapterError("enclosure_is_not_interval_object")

        for field, expected in (
            ("quantity", request["quantity"]),
            ("radial_order", request["radial_order"]),
            ("axial_order", request["axial_order"]),
            ("source_family", self.source_family),
            ("radial_degree", request["degree"]),
            ("axial_depth", request["depth"]),
        ):
            if field not in result or not _same(result[field], expected):
                raise EnclosureAdapterError("provider_metadata_mismatch", field)
        _require_exact_flag(result, _ORDINARY_FLAG_NAMES, "genuine_ordinary_derivatives")
        _require_exact_flag(result, _NONLINEAR_FLAG_NAMES, "actual_twenty_term_nonlinear_solution_enclosed")
        _require_exact_flag(result, _CENTER_FLAG_NAMES, "numerical_center_not_selected")

        if not callable(self.interval_validator):
            raise EnclosureAdapterError("strict_interval_validator_required")
        try:
            valid = self.interval_validator(enclosure)
        except Exception as exc:
            raise EnclosureAdapterError("interval_validator_error", type(exc).__name__) from exc
        if type(valid) is not bool or valid is not True:
            raise EnclosureAdapterError("interval_validator_did_not_accept_enclosure")
        return enclosure

    def component(self, record: Any, rho: Any, Z: Any) -> Any:
        """Return the complete directed enclosure for one Q4 component node."""

        request: dict[str, Any] | None = None
        try:
            request = self._read_request(record)
            if not callable(self.provider):
                raise EnclosureAdapterError("provider_callable_required")
            if not callable(self.interval_validator):
                # Reject before invoking a potentially expensive source
                # provider; an unvalidated enclosure cannot be accepted.
                raise EnclosureAdapterError("strict_interval_validator_required")
            # The call shape intentionally mirrors canonical component_enclosure.
            # No old provider or source module is discovered implicitly.
            result = self.provider(
                request["quantity"],
                rho,
                Z,
                radial_order=request["radial_order"],
                axial_order=request["axial_order"],
                degree=request["degree"],
                depth=request["depth"],
            )
            enclosure = self._validate_result(result, request)
            self._accept(request, result)
            return ComponentEnclosure(enclosure)
        except EnclosureAdapterError as error:
            self._reject(error, request)
            raise
        except Exception as exc:
            error = EnclosureAdapterError("provider_error", type(exc).__name__)
            self._reject(error, request)
            raise error from exc

    __call__ = component


Q4CoreEnclosureAdapter = CoreEnclosureAdapter
ComponentEnclosureAdapter = CoreEnclosureAdapter


def make_adapter(
    provider: Callable[..., Any] | None = None,
    source_family: Any = _MISSING,
    provider_descriptor: Any = _MISSING,
    degree: int = 6,
    depth: int = 3,
    interval_validator: Callable[[Any], Any] | None = None,
    source_evidence: Any = _MISSING,
    **kwargs: Any,
) -> CoreEnclosureAdapter:
    """Factory alias convenient for injecting the adapter into Q4 tests."""

    return CoreEnclosureAdapter(
        provider,
        source_family,
        provider_descriptor,
        degree,
        depth,
        interval_validator,
        source_evidence,
        **kwargs,
    )


build_adapter = make_adapter
make_component_adapter = make_adapter


__all__ = [
    "SCHEMA",
    "SUPPORTED_QUANTITIES",
    "EnclosureAdapterError",
    "ComponentEnclosure",
    "CoreEnclosureAdapter",
    "Q4CoreEnclosureAdapter",
    "ComponentEnclosureAdapter",
    "make_adapter",
    "build_adapter",
    "make_component_adapter",
]

