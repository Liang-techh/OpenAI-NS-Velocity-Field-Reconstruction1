"""Independent saved-coefficient core value and ordinary-jet kernel.

This is an *unexecuted* staging leaf.  It consumes one literal capsule that a
root runner freezes before invocation; it does not read files, import a source
provider, call a checker, or replay a recurrence.  The capsule supplies the
accepted ``/recomputed_center_packets/0`` packet, the authoritative common
fixed-point receipt, and the already-certified majorants/tail rows.

The numerical method in this file is deliberately new.  A saved coefficient
is kept as a directed interval.  For a polynomial

``p(rho) = sum(a[n] * rho**n, n=0..N)``

the exact rational change of variables ``x=rho/4`` gives Bernstein
coefficients

``b[j] = sum(c[n] * C(j,n)/C(N,n), n=0..j)``, ``c[n]=4**n*a[n]``.

The interval Bernstein coefficients are then evaluated by interval
de Casteljau at the four scaled radii ``1/16, 1/4, 1, 3``.  No interval
midpoint or endpoint is selected as a coefficient or a function value.

The capsule contract is intentionally explicit.  The root supplies:

* ``schema`` equal to :data:`CAPSULE_SCHEMA`;
* ``source_identity`` with the three accepted family/source/datum hashes;
* an authoritative common-check mapping whose admitted Banach, same-source,
  coefficient, and tail fields are true.  Raw pre-receipt false flags are
  preserved as metadata and are not consulted as gates;
* ``saved_packet`` (or an equivalent nested mapping) containing packet ``0``
  at exact ``Z=0``, degree ``24``, depth ``5``, rows ``A``, ``Uz``, ``P``, and
  fixed ``U0_Z_taylor``/``P0_Z_taylor``;
* selected degree-24 ``geometric_tail_rows`` at reference radius ``4.1``;
* positive interval majorants ``Phi_ball_norm_upper``,
  ``Psi_ball_norm_upper``, ``epsilon_core``, and the saved
  ``mixed_pressure_increment_Xh_norm_upper``;
* an accepted ``delta`` interval for the local ``Q_Z`` recovery functional;
* optional preservation/source-formula records.  They are copied as data only.

The component output remains a directed enclosure packet.  It is not a point
physical-field evaluation, does not materialize ``N``, and includes the new
G0-amplitude transport subcertificate merged before the once-only entry.  The local ``Q``,
``Q_rho``, and ``Q_Z`` packets are linear functionals of the same saved Uz
coefficients and retain their finite and tail parts separately.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_saved_nonlinear_core_value_jet_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_saved_nonlinear_core_value_jet_source_capsule_v1"
COMPONENT = "native_O3_saved_nonlinear_core_value_jet"
EXPECTED_FROZEN_CAPSULE_SHA256 = "012511ca881a0207e04b7f494ce8ca275d7ea8a0a04d2914a824b5f846d64e0e"

INTERVAL_DPS = 220
RADIAL_DEGREE = 24
AXIAL_DEPTH = 5
TAIL_REFERENCE_RADIUS = "4.1"
RHO_POINTS = ("1/16", "1/4", "1", "3")
REQUESTED_JETS = ((0, 0), (1, 0), (0, 1), (0, 2), (1, 1))

# The values are literal contracts for the root runner.  In particular, all
# old/provider/checker counters stay zero and the one overall numeric call is
# mapped to the single entry event.  Entry is recorded before interval math.
SUCCESS_COUNTS = {'N_materialization_calls': 0, 'ancestor_calls': 0, 'old_math_calls': 0, 'old_provider_calls': 0, 'provider_calls': 0, 'source_value_reducer_calls': 0, 'checker_calls': 0, 'selftest_calls': 0, 'source_family_calls': 1, 'new_math_family_calls': 4, 'new_bernstein_family_calls': 1, 'new_decasteljau_family_calls': 1, 'new_tail_dilation_family_calls': 1, 'numeric_evaluation_calls': 1, 'actual_core_value_jet_enclosure_calls': 1, 'new_mathematical_proofs': 1, 'new_g0_transport_family_calls': 1, 'new_g0_subproof_packets': 1}

SUCCESS_FLAGS = {'actual_core_value_jet_enclosure_evaluated': True, 'numeric_evaluation': True, 'source_only': True, 'metadata_only': False, 'same_actual_unique_analytic_core': True, 'finite_rows_consumed_as_intervals': True, 'ordinary_jets_computed': True, 'local_core_recovery_computed': True, 'shared_source_linear_functionals': True, 'mean_minus_axial_exact_axis_zero': True, 'axis_increment_join_structural': True, 'bounds_are_error_bounds_not_function_values': True, 'source_Z_jets_admitted': False, 'actual_source_Z_jets_admitted': False, 'actual_parameter_admission': False, 'physical_layer_completed': False, 'physical_field_materialized': False, 'full_point_physical_field_evaluation': False, 'full_goal': False, 'full_goal_complete': False, 'N_materialized': False, 'old_math_reentered': False, 'old_provider_called': False, 'source_value_point_selected': False, 'source_parameter_point_selected': False, 'nominal_anchor_selected': False, 'optional_new_function_hook_called': True, 'all_old_caps_and_errors_preserved': True, 'global_errors_preserved': True, 'actual_same_source_G0_enclosure': True, 'nonanchor_log_amplitude_enclosed': True}

CALLBACK_EVENT_COUNTS = {'family_entered': 1, 'bernstein_entered': 1, 'bernstein_complete': 1, 'decasteljau_entered': 1, 'decasteljau_complete': 1, 'tail_dilation_entered': 1, 'tail_dilation_complete': 1, 'proof_bundle': 1, 'family_completed': 1, 'g0_transport_entered': 1, 'g0_proof_bundle': 1, 'g0_transport_complete': 1}

COUNT_EVENT_MAP = {'source_family_calls': ['family_entered'], 'new_math_family_calls': ['bernstein_entered', 'decasteljau_entered', 'tail_dilation_entered', 'g0_transport_entered'], 'new_bernstein_family_calls': ['bernstein_entered'], 'new_decasteljau_family_calls': ['decasteljau_entered'], 'new_tail_dilation_family_calls': ['tail_dilation_entered'], 'numeric_evaluation_calls': ['family_entered'], 'actual_core_value_jet_enclosure_calls': ['family_entered'], 'new_mathematical_proofs': ['proof_bundle'], 'new_g0_transport_family_calls': ['g0_transport_entered'], 'new_g0_subproof_packets': ['g0_proof_bundle']}

EXPECTED_SOURCE_IDENTITY = {
    "actual_five_defect_family_sha256":
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256":
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256":
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

EXPECTED_HASH_PINS = {
    "addendum_sha256":
        "642c266ae551671a8798f2c7fdb73966b9ef9e5b046423eeaa91cb8ea470f5df",
    "common_fixed_point_checker_source_sha256":
        "f6fea64369400c6669eddda14840c88eedaf4c91a7fb91a0d78b47f7c36f9e0c",
    "coefficient_rebuild_json_sha256":
        "f4811f93a387913d3e3f2778e1b3f5e44adc3fd7ab5639c9e8fe2e4e858bb372",
    "coefficient_rebuild_source_sha256":
        "d37bfa6fe89b8955e17b38a00e1c8ae53b6725bb5a73c4648067737ce3be990f",
    "coefficient_rebuild_checker_source_sha256":
        "0aed6d045d79c53ec73c383f7839b81c4f40e219ba6d7d539f0479b0ce753c45",
    "coefficient_rebuild_checker_json_sha256":
        "02790017aeddc309624dcff46f2d4a2cd9d557075cc0e88096ff0c63cd9b65d0",
    "common_fixed_point_check_json_sha256":
        "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672",
    "raw_common_fixed_point_json_sha256":
        "c96ece0110427b141f870771d07e1a2db4a7b85ae0a4bd0b3ad7d8ff38bf8286",
}


class KernelInputError(ValueError):
    """Raised when the frozen literal capsule violates this narrow contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _authenticate_frozen_capsule(capsule: Mapping[str, Any]) -> dict[str, Any]:
    """Admit only this root-authenticated, source-bound frozen data packet."""
    raw = (json.dumps(capsule, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    actual = hashlib.sha256(raw).hexdigest()
    if actual != EXPECTED_FROZEN_CAPSULE_SHA256:
        raise KernelInputError("capsule is not the authenticated fixed-source packet")
    return {"canonical_capsule_sha256": actual, "all_source_joins_are_fixed_data": True}


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _json_value(value: Any) -> Any:
    """Keep source records JSON-shaped without evaluating opaque values."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if _is_mapping(value):
        return {str(key): _json_value(child) for key, child in value.items()}
    if _is_sequence(value):
        return [_json_value(child) for child in value]
    return str(value)


def _emit(
    emit: Callable[[str, Mapping[str, Any]], Any] | None,
    events: list[dict[str, Any]],
    name: str,
    payload: Mapping[str, Any] | None = None,
) -> None:
    if name not in CALLBACK_EVENT_COUNTS:
        raise KernelInputError("event outside frozen callback contract " + name)
    row: dict[str, Any] = {
        "event": name,
        "component": COMPONENT,
        "source_only": True,
        "numeric_evaluation": name in {
            "bernstein_complete",
            "decasteljau_complete",
            "tail_dilation_complete",
            "family_completed",
        },
    }
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def _new_interval_context() -> Any:
    """Create the private directed context used by this new method only."""

    # Deliberately local: no project module, old numerical helper, provider,
    # recurrence, checker, or source reducer is imported by this staging leaf.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _raw_endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower_endpoint(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, False)


def _upper_endpoint(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, True)


def _interval_from_record(context: Any, value: Any, label: str) -> Any:
    """Rehydrate a saved directed interval without selecting a representative."""

    if _is_mapping(value):
        # A few saved producers wrap interval records under one of these names.
        if not ("lower" in value or "lower_exact_mpf_tuple" in value):
            for key in (
                "interval",
                "enclosure",
                "record",
                "data",
                "value",
                "bound",
                "upper_bound",
            ):
                if key in value:
                    return _interval_from_record(context, value[key], label)
        lower_tuple = value.get("lower_exact_mpf_tuple")
        upper_tuple = value.get("upper_exact_mpf_tuple")
        if _is_sequence(lower_tuple) and _is_sequence(upper_tuple):
            if len(lower_tuple) != 4 or len(upper_tuple) != 4:
                raise KernelInputError(label + " MPF tuple must have length four")
            try:
                left = tuple(int(item) for item in lower_tuple)
                right = tuple(int(item) for item in upper_tuple)
            except (TypeError, ValueError) as error:
                raise KernelInputError(label + " MPF tuple is malformed") from error
            try:
                return context.make_mpf((left, right))
            except Exception as error:  # pragma: no cover - defensive boundary
                raise KernelInputError(label + " MPF tuple cannot be rehydrated") from error
        if "lower" in value and "upper" in value:
            try:
                return context.mpf((str(value["lower"]), str(value["upper"])))
            except Exception as error:  # pragma: no cover - defensive boundary
                raise KernelInputError(label + " endpoint strings are malformed") from error
        raise KernelInputError(label + " has no directed endpoints")
    if _is_sequence(value):
        if len(value) != 2:
            raise KernelInputError(label + " endpoint sequence must have length two")
        try:
            return context.mpf((str(value[0]), str(value[1])))
        except Exception as error:  # pragma: no cover - defensive boundary
            raise KernelInputError(label + " endpoint sequence is malformed") from error
    if value is None:
        raise KernelInputError(label + " is missing")
    try:
        return context.mpf(str(value))
    except Exception as error:  # pragma: no cover - defensive boundary
        raise KernelInputError(label + " is not interval-shaped") from error


def _display_endpoint(context: Any, value: Any) -> str:
    text = str(context.nstr(context.convert(value), n=INTERVAL_DPS)).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
        if "," in text:
            text = text.split(",", 1)[0].strip()
    return text


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    interval = context.convert(value)
    lower_raw, upper_raw = interval._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    width = upper - lower
    return {
        "lower": _display_endpoint(context, lower),
        "upper": _display_endpoint(context, upper),
        "width": _display_endpoint(context, width),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
    }


def _literal_zero(value: Any) -> bool:
    """Structural zero check used before numerical interval work."""

    if _is_mapping(value):
        if "lower" in value and "upper" in value:
            return _literal_zero(value["lower"]) and _literal_zero(value["upper"])
        if "value" in value:
            return _literal_zero(value["value"])
        return False
    if isinstance(value, bool) or value is None:
        return False
    text = str(value).strip().lower()
    if text.startswith("+"):
        text = text[1:]
    if text.startswith("-"):
        text = text[1:]
    if "e" in text:
        mantissa, exponent = text.split("e", 1)
        return _literal_zero(mantissa) and exponent.lstrip("+-").isdigit()
    return bool(text) and all(char in "0." for char in text) and "." in text or text == "0"


def _find_value(container: Any, names: Sequence[str], depth: int = 5) -> Any:
    """Find a named literal record through shallow capsule wrappers."""

    if depth < 0:
        return None
    if _is_mapping(container):
        for name in names:
            if name in container:
                return container[name]
        for value in container.values():
            found = _find_value(value, names, depth - 1)
            if found is not None:
                return found
    elif _is_sequence(container):
        for value in container:
            found = _find_value(value, names, depth - 1)
            if found is not None:
                return found
    return None


def _find_mapping(container: Any, names: Sequence[str], depth: int = 5) -> Mapping[str, Any] | None:
    found = _find_value(container, names, depth)
    return found if _is_mapping(found) else None


def _find_path(container: Any, path: Sequence[str]) -> Any:
    current = container
    for key in path:
        if _is_mapping(current) and key in current:
            current = current[key]
        else:
            return None
    return current


def _require_admitted(mapping: Mapping[str, Any], keys: Sequence[str], label: str) -> None:
    for key in keys:
        if key not in mapping or mapping[key] is not True:
            raise KernelInputError(label + "." + key + " is missing or not admitted")


def _require_false(mapping: Mapping[str, Any], keys: Sequence[str], label: str) -> None:
    for key in keys:
        if key in mapping and mapping[key] is not False:
            raise KernelInputError(label + "." + key + " must remain false")


def _validate_admission(capsule: Mapping[str, Any]) -> dict[str, Any]:
    common = _find_mapping(
        capsule,
        (
            "authoritative_common_fixed_point_check",
            "common_fixed_point_check",
            "accepted_common_fixed_point_check",
            "common_check",
            "accepted_common_check",
            "authoritative_common_receipt",
            "common_receipt",
        ),
    )
    if common is None:
        raise KernelInputError("authoritative common fixed-point check is missing")
    if not any(key in common for key in ("all_passed", "passed")):
        for nested_name in ("receipt", "admission", "accepted_receipt"):
            nested = common.get(nested_name)
            if _is_mapping(nested):
                common = nested
                break

    if common.get("all_passed") is not True and common.get("passed") is not True:
        raise KernelInputError("authoritative common fixed-point check did not pass")
    banach = _find_mapping(
        common,
        ("actual_Banach_and_common_coefficient_proof", "banach_proof"),
        depth=3,
    )
    if banach is None:
        raise KernelInputError("authoritative Banach coefficient proof is missing")
    _require_admitted(
        banach,
        (
            "passed",
            "unique_actual_analytic_fixed_point_for_each_fixed_implicit_source",
            "all_original_finite_rows_enclose_coefficients_of_this_unique_solution",
            "no_interval_overlap_or_fitted_source_equality_argument",
            "same_source_analytic_functions_not_just_shared_receipt_ids",
        ),
        "common.banach",
    )
    tail = _find_mapping(
        common,
        ("same_solution_all_order_tail_proof", "all_order_tail_proof"),
        depth=3,
    )
    if tail is None:
        raise KernelInputError("authoritative all-order tail proof is missing")
    _require_admitted(
        tail,
        (
            "passed",
            "same_actual_fixed_point_coefficient_norm_used",
            "finite_and_infinite_parts_of_same_analytic_solution",
            "ordinary_derivative_factorials_and_radial_units_preserved",
        ),
        "common.tail",
    )
    for key in (
        "current_core_actual_operator_Banach_fixed_point_admitted",
        "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point",
        "current_core_actual_operator_contraction_and_tail_certified",
        "all_scoped_checks_passed",
        "all_passed",
    ):
        if key not in common or common[key] is not True:
            raise KernelInputError("common." + key + " is missing or not admitted")

    rebuild = _find_mapping(
        capsule,
        (
            "coefficient_rebuild_check",
            "accepted_coefficient_rebuild_check",
            "rebuild_check",
            "packet_check",
            "coefficient_rebuild_receipt",
            "accepted_rebuild_receipt",
            "rebuild_receipt",
        ),
    )
    if rebuild is None:
        raise KernelInputError("coefficient rebuild check is missing")
    if not any(key in rebuild for key in ("all_passed", "passed", "evidence_fields")):
        for nested_name in ("receipt", "admission", "accepted_receipt"):
            nested = rebuild.get(nested_name)
            if _is_mapping(nested):
                rebuild = nested
                break
    if rebuild.get("all_passed") is not True and rebuild.get("passed") is not True:
        evidence = rebuild.get("evidence_fields")
        if not _is_mapping(evidence) or evidence.get("all_passed") is not True:
            raise KernelInputError("coefficient rebuild check did not pass")
    evidence = rebuild.get("evidence_fields", rebuild)
    if not _is_mapping(evidence):
        raise KernelInputError("coefficient rebuild evidence fields are missing")
    # The rebuild receipt is the second half of the authoritative admission
    # chain.  Its pre-receipt old-row flag is required to remain false; the
    # accepted saved packet must explicitly carry a certified infinite tail.
    _require_false(
        evidence,
        (
            "old_finite_coefficient_rows_read",
            "full_point_physical_field_evaluation",
            "temporal_scope",
        ),
        "rebuild",
    )
    if evidence.get("old_finite_coefficient_rows_read") is not False:
        raise KernelInputError("rebuild must explicitly certify old finite rows were not read")
    if evidence.get("infinite_radial_tails_bound") is not True:
        raise KernelInputError("rebuild must explicitly certify infinite radial tails")

    # Raw report flags intentionally are not read as gates.  The addendum
    # records them as expected pre-receipt values; only the authoritative
    # receipt above controls admission.
    return {
        "authoritative_common_check": _json_value(common),
        "coefficient_rebuild_check": _json_value(rebuild) if rebuild is not None else None,
        "raw_pending_flags_ignored_after_authoritative_receipt": True,
        "same_actual_unique_analytic_core": True,
        "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point": True,
    }


def _validate_source_identity(capsule: Mapping[str, Any]) -> dict[str, Any]:
    identity = _find_mapping(
        capsule,
        ("source_identity", "family_identity", "source_binding_identity"),
    )
    if identity is None:
        identity = {
            key: capsule.get(key) for key in EXPECTED_SOURCE_IDENTITY if key in capsule
        }
    if not identity:
        raise KernelInputError("source-family identity record is missing")
    for key, expected in EXPECTED_SOURCE_IDENTITY.items():
        if identity.get(key) != expected:
            raise KernelInputError("source-family identity mismatch at " + key)

    # Hash pins are data joins, never file reads.  Accept nested aliases used
    # by root capsules but reject a conflicting pin when one is supplied.
    supplied_hashes = _find_mapping(
        capsule, ("hash_pins", "source_hashes", "input_hashes", "artifact_hashes")
    )
    if supplied_hashes is not None:
        aliases = {
            "addendum_sha256": ("addendum_sha256", "route_addendum_sha256"),
            "common_fixed_point_checker_source_sha256": (
                "common_fixed_point_checker_source_sha256",
                "common_check_source_sha256",
                "current_core_common_fixed_point_source_sha256",
                "common_checker_source",
            ),
            "coefficient_rebuild_json_sha256": (
                "coefficient_rebuild_json_sha256",
                "rebuild_json_sha256",
                "core_coefficient_rebuild_json_sha256",
                "rebuild_json",
            ),
            "coefficient_rebuild_source_sha256": (
                "coefficient_rebuild_source_sha256",
                "rebuild_source_sha256",
                "core_coefficient_rebuild_source_sha256",
                "rebuild_source",
            ),
            "coefficient_rebuild_checker_source_sha256": (
                "coefficient_rebuild_checker_source_sha256",
                "rebuild_checker_source_sha256",
                "rebuild_checker_source",
            ),
            "coefficient_rebuild_checker_json_sha256": (
                "coefficient_rebuild_checker_json_sha256",
                "rebuild_checker_json_sha256",
                "rebuild_checker_json",
            ),
            "common_fixed_point_check_json_sha256": (
                "common_fixed_point_check_json_sha256",
                "common_check_json_sha256",
                "current_core_common_fixed_point_check_json_sha256",
                "common_check_json",
            ),
            "raw_common_fixed_point_json_sha256": (
                "raw_common_fixed_point_json_sha256",
                "raw_common_report_json_sha256",
                "current_core_common_fixed_point_json_sha256",
                "raw_common_report_json",
            ),
        }
        for canonical, names in aliases.items():
            value = next((supplied_hashes[name] for name in names if name in supplied_hashes), None)
            if value is not None and value != EXPECTED_HASH_PINS[canonical]:
                raise KernelInputError("source hash mismatch at " + canonical)

    return {
        "accepted_hashes": dict(EXPECTED_SOURCE_IDENTITY),
        "required_artifact_hash_pins": dict(EXPECTED_HASH_PINS),
        "supplied_artifact_hash_pins": _json_value(supplied_hashes),
        "capsule_identity": _json_value(identity),
        "formal_source_definition": _json_value(
            _find_value(capsule, ("formal_source_definition", "formal_source_identity"))
        ),
        "ordinary_order_definition": _json_value(
            _find_value(capsule, ("ordinary_order_definition", "ordinary_order_contract"))
        ),
        "raw_jlogC_definition": _json_value(
            _find_value(capsule, ("raw_jlogC_definition", "raw_jlogC_defs", "raw_logC_definition"))
        ),
    }


def _packet_from_capsule(capsule: Mapping[str, Any]) -> Mapping[str, Any]:
    direct = _find_value(
        capsule,
        (
            "saved_packet",
            "saved_coefficient_packet",
            "coefficient_packet",
            "packet_zero",
            "packet_0",
            "saved_packet_zero",
        ),
    )
    if _is_mapping(direct):
        # A direct packet has rows/fixed; a rebuild object is unwrapped below.
        if "rows" in direct and "fixed" in direct:
            return direct
        nested = _find_path(direct, ("recomputed_center_packets", "0"))
        if _is_mapping(nested):
            return nested
    for root_name in (
        "saved_coefficient_data",
        "coefficient_rebuild_json",
        "rebuild_json",
        "coefficient_rebuild",
    ):
        root = capsule.get(root_name)
        if _is_mapping(root):
            packets = root.get("recomputed_center_packets")
            if _is_mapping(packets):
                packet = packets.get("0")
                if _is_mapping(packet):
                    return packet
    packets = _find_value(capsule, ("recomputed_center_packets",))
    if _is_mapping(packets) and _is_mapping(packets.get("0")):
        return packets["0"]
    raise KernelInputError("saved recomputed_center_packets/0 is missing")


def _validate_packet(packet: Mapping[str, Any]) -> dict[str, Any]:
    z = packet.get("Z")
    if not _is_mapping(z) or not _literal_zero(z.get("lower")) or not _literal_zero(z.get("upper")):
        raise KernelInputError("saved packet center must be exact Z=0")
    if packet.get("radial_degree") != RADIAL_DEGREE:
        raise KernelInputError("saved packet radial degree must be 24")
    if packet.get("axial_depth") != AXIAL_DEPTH:
        raise KernelInputError("saved packet axial depth must be 5")
    rows = _as_mapping(packet.get("rows"), "saved packet rows")
    fixed = _as_mapping(packet.get("fixed"), "saved packet fixed data")
    for field in ("A", "Uz", "P"):
        table = rows.get(field)
        if not _is_sequence(table) or len(table) != RADIAL_DEGREE + 1:
            raise KernelInputError("saved rows." + field + " must have 25 rows")
        for n, row in enumerate(table):
            if not _is_sequence(row) or len(row) != 30 - n:
                raise KernelInputError(
                    "saved rows." + field + "[" + str(n) + "] must have "
                    + str(30 - n)
                    + " axial entries"
                )
            for k, coefficient in enumerate(row):
                if not _is_mapping(coefficient) and not _is_sequence(coefficient):
                    # Scalar literals are accepted as exact point intervals in
                    # small synthetic capsules; saved production rows are maps.
                    if coefficient is None:
                        raise KernelInputError(
                            "saved rows." + field + "[" + str(n) + "][" + str(k) + "] missing"
                        )
    for field in ("U0_Z_taylor", "P0_Z_taylor"):
        values = fixed.get(field)
        if not _is_sequence(values) or len(values) < AXIAL_DEPTH + 1:
            raise KernelInputError("saved fixed." + field + " is too short")
    return {
        "packet": packet,
        "rows": rows,
        "fixed": fixed,
        "row_units": _json_value(packet.get("row_units", {})),
        "packet_identity": {
            "Z_lower": "0.0",
            "Z_upper": "0.0",
            "radial_degree": RADIAL_DEGREE,
            "axial_depth": AXIAL_DEPTH,
            "rows": ["A[n][k]", "Uz[n][k]", "P[n][k]"],
        },
    }


def _tail_rows_from_capsule(capsule: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    candidate = _find_value(
        capsule,
        (
            "geometric_tail_rows",
            "selected_geometric_tail_rows",
            "tail_rows",
        ),
    )
    if _is_mapping(candidate):
        candidate = list(candidate.values())
    if not _is_sequence(candidate):
        raise KernelInputError("accepted geometric_tail_rows are missing")
    rows = [row for row in candidate if _is_mapping(row)]
    if len(rows) == 0:
        raise KernelInputError("accepted geometric_tail_rows are empty")
    return rows


def _select_tail_rows(capsule: Mapping[str, Any]) -> dict[tuple[int, int], Mapping[str, Any]]:
    rows = _tail_rows_from_capsule(capsule)
    selected: dict[tuple[int, int], Mapping[str, Any]] = {}
    for row in rows:
        if row.get("degree") != RADIAL_DEGREE:
            continue
        try:
            i = int(row.get("radial_order"))
            k = int(row.get("axial_order"))
        except (TypeError, ValueError) as error:
            raise KernelInputError("tail row derivative orders are malformed") from error
        if (i, k) in REQUESTED_JETS:
            if (i, k) in selected:
                raise KernelInputError("duplicate saved tail row for (" + str(i) + "," + str(k) + ")")
            if "tail_per_Xh_norm" not in row:
                raise KernelInputError("tail row lacks tail_per_Xh_norm")
            selected[(i, k)] = row
    missing = [jet for jet in REQUESTED_JETS if jet not in selected]
    if missing:
        raise KernelInputError("saved degree-24 tail rows missing " + repr(missing))
    return selected


def _majorant_value(capsule: Mapping[str, Any], names: Sequence[str], label: str) -> Any:
    value = _find_value(capsule, names)
    if value is None:
        raise KernelInputError(label + " is missing from accepted capsule")
    return value


def _validate_majorants(capsule: Mapping[str, Any]) -> dict[str, Any]:
    phi = _majorant_value(
        capsule,
        (
            "Phi_ball_norm_upper",
            "fullPhi_ball_norm",
            "full_Phi_ball_norm",
            "Phi_ball_upper",
            "fullPhi_ball_norm_upper",
        ),
        "Phi_ball_norm_upper",
    )
    psi = _majorant_value(
        capsule,
        (
            "Psi_ball_norm_upper",
            "fullPsi_ball_norm",
            "full_Psi_ball_norm",
            "Psi_ball_upper",
            "fullPsi_ball_norm_upper",
        ),
        "Psi_ball_norm_upper",
    )
    epsilon = _majorant_value(
        capsule,
        (
            "epsilon_core",
            "accepted_epsilon_core",
            "core_epsilon",
            "epsilon_core_upper",
        ),
        "epsilon_core",
    )
    pressure = _majorant_value(
        capsule,
        (
            "mixed_pressure_increment_Xh_norm_upper",
            "pressure_increment_Xh_norm_upper",
            "saved_pressure_increment_norm",
        ),
        "mixed_pressure_increment_Xh_norm_upper",
    )
    delta = _majorant_value(
        capsule,
        (
            "delta",
            "accepted_delta",
            "delta_parameter",
            "core_delta",
            "delta_interval",
            "accepted_delta_interval",
        ),
        "delta",
    )
    return {
        "Phi_ball_norm_upper": phi,
        "Psi_ball_norm_upper": psi,
        "epsilon_core": epsilon,
        "mixed_pressure_increment_Xh_norm_upper": pressure,
        "delta": delta,
        "majorant_source": _json_value(
            _find_mapping(capsule, ("majorants", "accepted_majorants"))
        ),
    }


def _interval_zero(context: Any) -> Any:
    return context.mpf("0")


def _interval_one(context: Any) -> Any:
    return context.mpf("1")


def _positive_error_interval(context: Any, upper_value: Any) -> Any:
    upper = _upper_endpoint(context, upper_value)
    zero_raw = _interval_zero(context)._mpi_[0]
    return context.make_mpf((zero_raw, upper._mpi_[1]))


def _symmetric_error(context: Any, upper_value: Any) -> Any:
    upper = _upper_endpoint(context, upper_value)
    lower = -upper
    return context.make_mpf((lower._mpi_[0], upper._mpi_[1]))


def _factorial(value: int) -> int:
    out = 1
    for item in range(2, value + 1):
        out *= item
    return out


def _falling(value: int, order: int) -> int:
    out = 1
    for item in range(order):
        out *= value - item
    return out


def _binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    out = 1
    for item in range(1, k + 1):
        out = (out * (n - k + item)) // item
    return out


def _rational(context: Any, numerator: int, denominator: int) -> Any:
    if denominator == 0:
        raise KernelInputError("zero denominator in exact rational transform")
    return context.mpf(str(numerator)) / context.mpf(str(denominator))


def _monomial_to_bernstein(
    context: Any,
    coefficients: Sequence[Any],
    degree: int,
) -> list[Any]:
    """Apply the exact-rational monomial-to-Bernstein transform on [0,4]."""

    if len(coefficients) != degree + 1:
        raise KernelInputError("monomial coefficient length does not match degree")
    scaled = [
        coefficients[n] * context.mpf(str(4 ** n)) for n in range(degree + 1)
    ]
    bernstein: list[Any] = []
    for j in range(degree + 1):
        value = _interval_zero(context)
        for n in range(j + 1):
            value += scaled[n] * _rational(
                context, _binomial(j, n), _binomial(degree, n)
            )
        bernstein.append(value)
    return bernstein


def _decasteljau(context: Any, coefficients: Sequence[Any], x: Any) -> Any:
    """Evaluate a Bernstein interval by de Casteljau, preserving intervals."""

    values = list(coefficients)
    while len(values) > 1:
        values = [
            (context.one - x) * values[j] + x * values[j + 1]
            for j in range(len(values) - 1)
        ]
    return values[0] if values else _interval_zero(context)


def _coefficient(
    context: Any,
    rows: Mapping[str, Any],
    field: str,
    n: int,
    k: int,
    average: bool,
    cache: dict[tuple[str, int, int], Any],
) -> Any:
    key = (field, n, k)
    if key not in cache:
        table = rows[field]
        cache[key] = _interval_from_record(
            context, table[n][k], "rows." + field + "[" + str(n) + "][" + str(k) + "]"
        )
    value = cache[key]
    # M divides the original Uz coefficient before radial differentiation.
    return value / context.mpf(str(n + 1)) if average else value


def _finite_bernstein_result(
    context: Any,
    rows: Mapping[str, Any],
    field: str,
    i: int,
    k: int,
    start_n: int,
    average: bool,
    rho: Any,
    coefficient_cache: dict[tuple[str, int, int], Any],
    weights: Sequence[Any] | None = None,
    axis_override: Any | None = None,
    functional_label: str | None = None,
) -> dict[str, Any]:
    degree = RADIAL_DEGREE - i
    coefficients: list[Any] = []
    radial_factorials = _factorial(k)
    if weights is not None and len(weights) != RADIAL_DEGREE + 1:
        raise KernelInputError("linear-functional weight table must have 25 entries")
    for power in range(degree + 1):
        n = power + i
        if n < start_n:
            coefficients.append(_interval_zero(context))
            continue
        if n == 0 and axis_override is not None:
            # The n=0 source is the saved axis packet, not a second copy of
            # the row-zero box.  This keeps shared-axis structure explicit.
            value = axis_override
        else:
            value = _coefficient(context, rows, field, n, k, average, coefficient_cache)
        if weights is not None:
            value *= weights[n]
        value *= context.mpf(str(_falling(n, i)))
        value *= context.mpf(str(radial_factorials))
        coefficients.append(value)
    bernstein = _monomial_to_bernstein(context, coefficients, degree)
    x = rho / context.mpf("4")
    value = _decasteljau(context, bernstein, x)
    return {
        "field": field,
        "radial_order": i,
        "axial_order": k,
        "start_radial_index": start_n,
        "linear_functional": functional_label,
        "source_weights_applied_before_radial_factor": weights is not None,
        "axis_source_override": axis_override is not None,
        "degree": degree,
        "domain": ["0", "4"],
        "coordinate": "rho=Lambda*R",
        "x": _interval_record(context, x),
        "monomial_coefficients": [_interval_record(context, item) for item in coefficients],
        "bernstein_coefficients": [_interval_record(context, item) for item in bernstein],
        "decasteljau_value": _interval_record(context, value),
        "finite_polynomial": _interval_record(context, value),
        "transform": {
            "kind": "exact_rational_monomial_to_Bernstein_on_[0,4]",
            "formula": "b[j]=sum_{n=0}^j (4^n a[n]) C(j,n)/C(degree,n)",
            "rational_coefficients_exact": True,
            "interval_coefficients_preserved": True,
        },
        "decasteljau": {
            "kind": "interval_deCasteljau",
            "convex_blend_parameter": "x=rho/4",
            "interval_value_not_midpoint": True,
        },
    }


def _axis_value(
    context: Any,
    fixed: Mapping[str, Any],
    key: str,
    k: int,
) -> Any:
    values = fixed.get(key)
    if not _is_sequence(values) or k >= len(values):
        raise KernelInputError("fixed." + key + " lacks axial coefficient " + str(k))
    value = _interval_from_record(context, values[k], "fixed." + key + "[" + str(k) + "]")
    return value * context.mpf(str(_factorial(k)))


def _tail_bound(
    context: Any,
    row: Mapping[str, Any],
    rho: Any,
    component_norm_upper: Any,
) -> dict[str, Any]:
    """Dilate an accepted degree-24 tail row from radius 4.1.

    This uses only the saved positive tail row.  It does not recompute a
    recurrence or an analytic tail helper.  For the first omitted radial power
    ``n0=25``, positive-coefficient dilation gives
    ``T(rho)<=T(4.1)*(rho/4.1)**(25-i)``.
    """

    tail_norm = _interval_from_record(
        context, row["tail_per_Xh_norm"], "tail_per_Xh_norm"
    )
    tail_norm_upper = _upper_endpoint(context, tail_norm)
    norm_upper = _upper_endpoint(context, component_norm_upper)
    reference = context.mpf(TAIL_REFERENCE_RADIUS)
    ratio = rho / reference
    i = int(row["radial_order"])
    exponent = RADIAL_DEGREE + 1 - i
    dilation = ratio ** exponent
    bound = tail_norm_upper * norm_upper * dilation
    error = _positive_error_interval(context, bound)
    return {
        "degree": RADIAL_DEGREE,
        "radial_order": i,
        "axial_order": int(row["axial_order"]),
        "reference_radius": TAIL_REFERENCE_RADIUS,
        "tail_per_Xh_norm_saved": _interval_record(context, tail_norm),
        "component_norm_upper_saved": _interval_record(context, component_norm_upper),
        "component_norm_upper_endpoint": _display_endpoint(context, norm_upper),
        "dilation_exponent": exponent,
        "dilation_factor": _interval_record(context, dilation),
        "bound_upper": _interval_record(context, error),
        "error_interval": _interval_record(context, error),
        "is_error_bound_not_function_value": True,
        "positive_coefficient_dilation_theorem": {
            "omitted_first_radial_index": RADIAL_DEGREE + 1,
            "inequality": "rho^(n-i)<= (rho/4.1)^(25-i)*4.1^(n-i), 0<rho<=4.1",
            "saved_reference_row_consumed": True,
            "old_tail_helper_recomputed": False,
        },
    }


def _field_norms(context: Any, majorants: Mapping[str, Any]) -> dict[str, Any]:
    phi = _interval_from_record(
        context, majorants["Phi_ball_norm_upper"], "Phi_ball_norm_upper"
    )
    psi = _interval_from_record(
        context, majorants["Psi_ball_norm_upper"], "Psi_ball_norm_upper"
    )
    epsilon = _interval_from_record(context, majorants["epsilon_core"], "epsilon_core")
    pressure = _interval_from_record(
        context,
        majorants["mixed_pressure_increment_Xh_norm_upper"],
        "mixed_pressure_increment_Xh_norm_upper",
    )
    # Only upper endpoints are used for nonnegative error caps.  The lower
    # endpoint is retained in the source record, never promoted to a value.
    epsilon_upper = _upper_endpoint(context, epsilon)
    psi_upper = _upper_endpoint(context, psi)
    uz_norm_upper = epsilon_upper * psi_upper
    return {
        "Phi": phi,
        "Psi": psi,
        "Uz": uz_norm_upper,
        "Mz_over_R": uz_norm_upper,
        "epsilon_core": epsilon,
        "P": pressure,
        "Uz_M_error_norm": uz_norm_upper,
        "norm_semantics": {
            "Phi": "upper endpoint of accepted Phi_ball_norm_upper used only as error cap",
            "Uz": "upper endpoint epsilon_core times Psi_ball_norm_upper used only as error cap",
            "Mz_over_R": "same epsilon_core*Psi_ball_norm cap; averaging is no larger",
            "P": "saved mixed_pressure_increment_Xh_norm_upper used only as error cap",
        },
    }


def _interval_abs_upper(context: Any, value: Any) -> Any:
    """Return a directed upper cap for ``abs(value)`` without a midpoint."""

    lower = _lower_endpoint(context, value)
    upper = _upper_endpoint(context, value)
    lower_point = context.mpf(_display_endpoint(context, lower))
    upper_point = context.mpf(_display_endpoint(context, upper))
    zero = context.mpf("0")
    if upper_point <= zero:
        return -lower
    if lower_point >= zero:
        return upper
    negative_lower = -lower
    return negative_lower if negative_lower >= upper else upper


def _qz_weights(context: Any, delta: Any) -> list[Any]:
    """Weights for the shared-source functional ``2 U-(1-delta) M``."""

    one = context.mpf("1")
    return [
        context.mpf("2") - (one - delta) / context.mpf(str(n + 1))
        for n in range(RADIAL_DEGREE + 1)
    ]


def _mean_minus_axial_weights(context: Any) -> list[Any]:
    """Weights for ``M-U``; n=0 is exactly the shared axis zero."""

    return [
        context.mpf(str(-n)) / context.mpf(str(n + 1))
        for n in range(RADIAL_DEGREE + 1)
    ]


def _build_negated_mean_cell(
    context: Any,
    base: Mapping[str, Any],
) -> dict[str, Any]:
    """Recover ``-M``-derived cells while retaining the saved tail cap."""

    finite_base = _interval_from_record(
        context, base["finite_total"], "mean finite total"
    )
    finite = -finite_base
    tail = base["certified_tail"]
    tail_cap = _interval_from_record(context, tail["error_interval"], "mean tail")
    enclosure = finite + _symmetric_error(context, tail_cap)
    base_axis = _interval_from_record(context, base["axis"]["interval"], "mean axis")
    return {
        "finite_source": _json_value(base["increment"]["finite_bernstein"]),
        "axis": {
            "active": bool(base["axis"].get("active")),
            "interval": _interval_record(context, -base_axis),
            "source": "negative of the same saved U0/M axis packet",
            "not_subtracted_from_increment": True,
        },
        "finite_total": _interval_record(context, finite),
        "certified_tail": _json_value(tail),
        "tail_error_cap": _interval_record(context, tail_cap),
        "enclosure": _interval_record(context, enclosure),
        "point_value": False,
        "error_bound_not_function_value": True,
    }


def _build_qz_cell(
    context: Any,
    packet_data: Mapping[str, Any],
    norm_data: Mapping[str, Any],
    tail_rows: Mapping[tuple[int, int], Mapping[str, Any]],
    delta: Any,
    rho_text: str,
    coefficient_cache: dict[tuple[str, int, int], Any],
    mzz_cell: Mapping[str, Any],
) -> dict[str, Any]:
    """Build Q_Z from one shared-source weighted polynomial and M_ZZ.

    The weighted source polynomial is formed before interval subtraction, so
    the U/M axis is not recovered by subtracting two unrelated boxes.
    """

    rho = context.mpf(rho_text)
    weights = _qz_weights(context, delta)
    axis = _axis_value(context, packet_data["fixed"], "U0_Z_taylor", 0)
    finite_weighted = _finite_bernstein_result(
        context,
        packet_data["rows"],
        "Uz",
        0,
        0,
        0,
        False,
        rho,
        coefficient_cache,
        weights=weights,
        axis_override=axis,
        functional_label="2U-(1-delta)M; weight 2-(1-delta)/(n+1)",
    )
    weighted_value = _interval_from_record(
        context, finite_weighted["decasteljau_value"], "QZ weighted finite"
    )
    mzz_value = _interval_from_record(context, mzz_cell["finite_total"], "MZZ finite")
    finite_total = weighted_value - mzz_value

    # The saved positive tail cap is dilated separately for the weighted U
    # tail and the M_ZZ tail, then combined as a valid error bound.
    delta_deviation = _interval_abs_upper(context, context.mpf("1") - delta)
    weighted_factor = context.mpf("2") + delta_deviation / context.mpf(
        str(RADIAL_DEGREE + 2)
    )
    weighted_norm = norm_data["Uz_M_error_norm"] * weighted_factor
    weighted_tail = _tail_bound(
        context, tail_rows[(0, 0)], rho, weighted_norm
    )
    mzz_tail = mzz_cell["certified_tail"]
    weighted_cap = _interval_from_record(
        context, weighted_tail["error_interval"], "QZ weighted tail"
    )
    mzz_cap = _interval_from_record(context, mzz_tail["error_interval"], "MZZ tail")
    combined_cap = weighted_cap + mzz_cap
    enclosure = finite_total + _symmetric_error(context, combined_cap)
    return {
        "definition": "Q_Z = 2U-(1-delta)M-M_ZZ at Z=0",
        "source_weight": "2-(1-delta)/(n+1)",
        "delta_interval": _interval_record(context, delta),
        "shared_source_axis": {
            "U0_interval": _interval_record(context, axis),
            "weight_n0": _interval_record(context, weights[0]),
            "MZZ_axis_interval": _json_value(mzz_cell["axis"]["interval"]),
            "axis_not_recovered_by_box_subtraction": True,
        },
        "weighted_finite_bernstein": finite_weighted,
        "subtracted_MZZ_finite": _interval_record(context, mzz_value),
        "finite_total": _interval_record(context, finite_total),
        "tail_components": {
            "weighted_U": weighted_tail,
            "MZZ": _json_value(mzz_tail),
            "combined_error_cap": _interval_record(context, combined_cap),
            "weighted_tail_multiplier_upper": _interval_record(context, weighted_factor),
        },
        "certified_tail": {
            "error_interval": _interval_record(context, combined_cap),
            "tail_is_error_bound": True,
            "components_are_combined_by_valid_triangle_bound": True,
        },
        "enclosure": _interval_record(context, enclosure),
        "point_value": False,
        "error_bound_not_function_value": True,
    }


def _build_mean_minus_axial_cell(
    context: Any,
    packet_data: Mapping[str, Any],
    norm_data: Mapping[str, Any],
    tail_rows: Mapping[tuple[int, int], Mapping[str, Any]],
    i: int,
    k: int,
    rho_text: str,
    coefficient_cache: dict[tuple[str, int, int], Any],
    weights: Sequence[Any],
) -> dict[str, Any]:
    """Evaluate ``M-U`` with shared-source weights and exact zero axis."""

    rho = context.mpf(rho_text)
    axis = _axis_value(context, packet_data["fixed"], "U0_Z_taylor", k)
    finite = _finite_bernstein_result(
        context,
        packet_data["rows"],
        "Uz",
        i,
        k,
        0,
        False,
        rho,
        coefficient_cache,
        weights=weights,
        axis_override=axis,
        functional_label="M-U; weight -n/(n+1)",
    )
    finite_value = _interval_from_record(
        context, finite["decasteljau_value"], "M-U finite"
    )
    # Use the requested conservative 2 * Uz tail cap for this correlated
    # difference, while preserving the exact shared axis zero in metadata.
    tail = _tail_bound(
        context,
        tail_rows[(i, k)],
        rho,
        norm_data["Uz_M_error_norm"] * context.mpf("2"),
    )
    tail_cap = _interval_from_record(context, tail["error_interval"], "M-U tail")
    enclosure = finite_value + _symmetric_error(context, tail_cap)
    return {
        "definition": "M-U with weights -n/(n+1) on the original Uz coefficients",
        "axis": {
            "interval": _interval_record(context, _interval_zero(context)),
            "exact_zero_from_shared_U0": True,
            "source_U0_interval": _interval_record(context, axis),
        },
        "finite_bernstein": finite,
        "finite_total": _interval_record(context, finite_value),
        "certified_tail": tail,
        "tail_bound_statement": "tail <= 2 * UzTail",
        "enclosure": _interval_record(context, enclosure),
        "point_value": False,
        "error_bound_not_function_value": True,
    }


def _jet_key(i: int, k: int) -> str:
    return "radial_" + str(i) + "_axial_" + str(k)


def _build_jet_cell(
    context: Any,
    packet_data: Mapping[str, Any],
    norm_data: Mapping[str, Any],
    tail_rows: Mapping[tuple[int, int], Mapping[str, Any]],
    field: str,
    i: int,
    k: int,
    rho_text: str,
    coefficient_cache: dict[tuple[str, int, int], Any],
) -> dict[str, Any]:
    rho = context.mpf(rho_text)
    rows = packet_data["rows"]
    fixed = packet_data["fixed"]
    if field == "Phi":
        row_field = "A"
        start_n = 0
        average = False
        axis_key = None
        component_norm = norm_data["Phi"]
    elif field == "Uz":
        row_field = "Uz"
        start_n = 1
        average = False
        axis_key = "U0_Z_taylor"
        component_norm = norm_data["Uz_M_error_norm"]
    elif field == "Mz_over_R":
        row_field = "Uz"
        start_n = 1
        average = True
        axis_key = "U0_Z_taylor"
        component_norm = norm_data["Uz_M_error_norm"]
    elif field == "P_scaled":
        row_field = "P"
        start_n = 1
        average = False
        axis_key = "P0_Z_taylor"
        component_norm = norm_data["P"]
    else:  # pragma: no cover - fixed internal dispatch
        raise KernelInputError("unknown jet field " + field)

    finite = _finite_bernstein_result(
        context,
        rows,
        row_field,
        i,
        k,
        start_n,
        average,
        rho,
        coefficient_cache,
    )
    if axis_key is None:
        axis = _interval_zero(context)
        axis_source = "rows.A[0][k] remains inside finite Bernstein polynomial"
        axis_active = False
        total_finite = _interval_from_record(
            context, finite["decasteljau_value"], "finite Bernstein result"
        )
    elif i == 0:
        axis = _axis_value(context, fixed, axis_key, k)
        axis_source = "fixed." + axis_key + "[k] times k!"
        axis_active = True
        increment = _interval_from_record(
            context, finite["decasteljau_value"], "finite increment Bernstein result"
        )
        total_finite = axis + increment
    else:
        axis = _interval_zero(context)
        axis_source = "axis affine/pressure seed has zero radial derivative"
        axis_active = False
        total_finite = _interval_from_record(
            context, finite["decasteljau_value"], "finite increment Bernstein result"
        )
    tail = _tail_bound(context, tail_rows[(i, k)], rho, component_norm)
    tail_interval = _interval_from_record(context, tail["error_interval"], "tail error")
    total = total_finite + _symmetric_error(context, tail_interval)
    return {
        "field": field,
        "radial_order": i,
        "axial_order": k,
        "rho": rho_text,
        "Z": "0",
        "ordinary_derivative_definition": {
            "radial_coordinate": "rho=Lambda*R",
            "axial_coordinate": "Z",
            "axial_order_is_Taylor_index": True,
            "axial_factor": "k!",
            "axial_factor_integer": _factorial(k),
            "radial_factor": "n!/(n-i)!",
            "ordinary_not_temporal": True,
        },
        "axis": {
            "active": axis_active,
            "source": axis_source,
            "interval": _interval_record(context, axis),
            "function_definition": (
                "U0 axis Taylor seed" if field in {"Uz", "Mz_over_R"}
                else "epsilon_core*physical pressure axis Taylor seed"
                if field == "P_scaled"
                else "Phi row-zero coefficient"
            ),
            "not_subtracted_from_increment": True,
        },
        "increment": {
            "finite_bernstein": finite,
            "start_radial_index": start_n,
            "average_Uz_before_radial_derivative": average,
            "interval": _interval_record(
                context,
                _interval_from_record(context, finite["decasteljau_value"], "finite increment"),
            ),
        },
        "finite_total": _interval_record(context, total_finite),
        "certified_tail": tail,
        "enclosure": _interval_record(context, total),
        "point_value": False,
        "function_value": False,
        "error_bound_not_function_value": True,
    }


def derive(
    capsule: Mapping[str, Any],
    emit: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Consume the frozen saved packet once and return directed jet enclosures.

    The root runner owns source freezing, entry persistence, invocation, and
    complete-result saving.  This function intentionally raises on malformed
    input after preserving every callback event already emitted; it never
    resets counters or manufactures a rejected zero result.
    """

    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")
    events: list[dict[str, Any]] = []
    # This is the durable failure prefix and is emitted before any source
    # validation, interval-context creation, or numerical conversion.
    _emit(
        emit,
        events,
        "family_entered",
        {
            "capsule_schema": capsule.get("schema") if _is_mapping(capsule) else None,
            "entry_before_interval_math": True,
            "overall_numeric_evaluation_call": 1,
            "source_only": True,
            "N_materialized": False,
        },
    )
    try:
        if not _is_mapping(capsule):
            raise KernelInputError("capsule must be a mapping")
        if capsule.get("schema") != CAPSULE_SCHEMA:
            raise KernelInputError("unexpected saved nonlinear core capsule schema")
        status = capsule.get("status")
        if _is_mapping(status):
            if status.get("metadata_only") not in (None, True):
                raise KernelInputError("input source capsule must be metadata-only before entry")
            for key in ("executed", "numeric_evaluation", "new_recipe_executed"):
                if key in status and status[key] not in (False, None, 0):
                    raise KernelInputError("input source capsule is already executed: " + key)

        frozen_input_authentication = _authenticate_frozen_capsule(capsule)
        admission = _validate_admission(capsule)
        identity = _validate_source_identity(capsule)
        identity["frozen_input_authentication"] = frozen_input_authentication
        packet = _validate_packet(_packet_from_capsule(capsule))
        selected_tail_rows = _select_tail_rows(capsule)
        majorants = _validate_majorants(capsule)

        _emit(
            emit,
            events,
            "bernstein_entered",
            {
                "method": "exact_rational_monomial_to_Bernstein_on_[0,4]",
                "degree": RADIAL_DEGREE,
                "saved_rows_as_directed_intervals": True,
                "midpoint_coefficients_selected": False,
                "ordinary_jets": [list(item) for item in REQUESTED_JETS],
            },
        )
        context = _new_interval_context()
        norm_data = _field_norms(context, majorants)
        delta = _interval_from_record(context, majorants["delta"], "delta")
        amplitude_Z0 = bound_core_amplitude_at_Z0(
            context,
            capsule["accepted_parameter_payload"],
            capsule["accepted_anchor"],
            capsule["saved_Z0_gauge"],
            emit=lambda name, row: _emit(emit, events, name, row),
        )
        coefficient_cache: dict[tuple[str, int, int], Any] = {}

        _emit(
            emit,
            events,
            "decasteljau_entered",
            {
                "radii": list(RHO_POINTS),
                "x_values": ["1/64", "1/16", "1/4", "3/4"],
                "interval_deCasteljau": True,
            },
        )
        _emit(
            emit,
            events,
            "tail_dilation_entered",
            {
                "method": "saved_degree24_positive_coefficient_tail_dilation",
                "reference_radius": TAIL_REFERENCE_RADIUS,
                "saved_tail_rows": len(selected_tail_rows),
                "old_tail_helper_called": False,
            },
        )

        finite_bernstein: dict[str, dict[str, dict[str, Any]]] = {}
        decasteljau: dict[str, dict[str, dict[str, Any]]] = {}
        ordinary_jets: dict[str, dict[str, dict[str, Any]]] = {}
        fields = ("Phi", "Uz", "Mz_over_R", "P_scaled")
        for field in fields:
            finite_bernstein[field] = {}
            decasteljau[field] = {}
            ordinary_jets[field] = {}
            for i, k in REQUESTED_JETS:
                key = _jet_key(i, k)
                finite_bernstein[field][key] = {}
                decasteljau[field][key] = {}
                ordinary_jets[field][key] = {}
                for rho_text in RHO_POINTS:
                    cell = _build_jet_cell(
                        context,
                        packet,
                        norm_data,
                        selected_tail_rows,
                        field,
                        i,
                        k,
                        rho_text,
                        coefficient_cache,
                    )
                    finite_bernstein[field][key][rho_text] = cell["increment"]["finite_bernstein"]
                    decasteljau[field][key][rho_text] = {
                        "finite": cell["increment"]["finite_bernstein"]["decasteljau"],
                        "axis": cell["axis"],
                        "finite_total": cell["finite_total"],
                    }
                    ordinary_jets[field][key][rho_text] = cell

        # P is kept under the explicit scaled name internally and also under
        # the short field name expected by callers: it is epsilon_core times
        # physical pressure, never an independently reconstructed pressure.
        finite_bernstein["P"] = finite_bernstein["P_scaled"]
        decasteljau["P"] = decasteljau["P_scaled"]
        ordinary_jets["P"] = ordinary_jets["P_scaled"]

        # Local core recovery uses the same saved Uz coefficients.  Q and
        # Q_rho are negatives of the already-built M_Z and M_rhoZ cells, so
        # their finite and certified-tail packets stay visibly separate.
        core_recovery: dict[str, Any] = {
            "Q": {
                "definition": "Q=-M_Z at Z=0",
                "source": "same Uz coefficients through M=Uz_n/(n+1)",
                "values": {},
            },
            "Q_rho": {
                "definition": "Q_rho=-M_rhoZ at Z=0",
                "source": "same Uz coefficients through M=Uz_n/(n+1)",
                "values": {},
            },
            "Q_Z": {
                "definition": "Q_Z=2U-(1-delta)M-M_ZZ at Z=0",
                "source": "shared-source weighted coefficients, then M_ZZ source data",
                "values": {},
            },
            "shared_source_axis_structure": True,
            "uncorrelated_box_subtraction": False,
        }
        qz_mzz_key = _jet_key(0, 2)
        for rho_text in RHO_POINTS:
            q_base = ordinary_jets["Mz_over_R"][_jet_key(0, 1)][rho_text]
            qrho_base = ordinary_jets["Mz_over_R"][_jet_key(1, 1)][rho_text]
            q_cell = _build_negated_mean_cell(context, q_base)
            qrho_cell = _build_negated_mean_cell(context, qrho_base)
            qz_cell = _build_qz_cell(
                context,
                packet,
                norm_data,
                selected_tail_rows,
                delta,
                rho_text,
                coefficient_cache,
                ordinary_jets["Mz_over_R"][qz_mzz_key][rho_text],
            )
            for cell in (q_cell, qrho_cell, qz_cell):
                cell["rho"] = rho_text
                cell["Z"] = "0"
            core_recovery["Q"]["values"][rho_text] = q_cell
            core_recovery["Q_rho"]["values"][rho_text] = qrho_cell
            core_recovery["Q_Z"]["values"][rho_text] = qz_cell

        # Optional correlated source functional: M-U has an exact zero axis,
        # and its tail is deliberately retained as the conservative 2*Uz cap.
        mean_minus_axial: dict[str, Any] = {
            "definition": "M-U with weights -n/(n+1) on original Uz coefficients",
            "axis_exact_zero": True,
            "tail_bound": "tail <= 2 * UzTail",
            "values": {},
        }
        mean_weights = _mean_minus_axial_weights(context)
        for i, k in REQUESTED_JETS:
            key = _jet_key(i, k)
            mean_minus_axial["values"][key] = {}
            for rho_text in RHO_POINTS:
                mean_minus_axial["values"][key][rho_text] = _build_mean_minus_axial_cell(
                    context,
                    packet,
                    norm_data,
                    selected_tail_rows,
                    i,
                    k,
                    rho_text,
                    coefficient_cache,
                    mean_weights,
                )

        _emit(
            emit,
            events,
            "decasteljau_complete",
            {
                "fields": list(fields),
                "jet_count": len(REQUESTED_JETS),
                "radius_count": len(RHO_POINTS),
                "local_core_recovery": ["Q", "Q_rho", "Q_Z"],
                "mean_minus_axial": True,
                "numeric_evaluation": True,
                "point_values_selected": False,
            },
        )
        _emit(
            emit,
            events,
            "bernstein_complete",
            {
                "domain": ["0", "4"],
                "radii": list(RHO_POINTS),
                "rational_transform_only": True,
                "saved_rows_consumed": len(coefficient_cache),
                "coefficient_interval_cache_entries": len(coefficient_cache),
                "midpoint_coefficients_selected": False,
            },
        )
        # Tail arithmetic is performed while constructing cells above; this
        # entry/completion pair records that one accepted tail family was used.
        _emit(
            emit,
            events,
            "tail_dilation_complete",
            {
                "requested_jets": [list(item) for item in REQUESTED_JETS],
                "tail_is_error_not_value": True,
                "positive_upper_caps_retained": True,
                "norms": {
                    "Phi": _interval_record(context, norm_data["Phi"]),
                    "Uz_M": _interval_record(context, norm_data["Uz_M_error_norm"]),
                    "P": _interval_record(context, norm_data["P"]),
                },
            },
        )

        proof = {
            "kind": "saved_interval_Bernstein_deCasteljau_with_positive_tail_dilation_v1",
            "source_identity": identity,
            "authoritative_admission": admission,
            "packet_identity": packet["packet_identity"],
            "ordinary_jet_definition": {
                "requested": [list(item) for item in REQUESTED_JETS],
                "axial_index_factor": "k!",
                "radial_factor": "n!/(n-i)!",
                "average": "M_n=Uz_n/(n+1) before radial differentiation",
                "axis": "U0 and pressure axis are joined structurally, never by subtracting boxes",
                "pressure": "P_scaled=epsilon_core*physical_P; row zero is the saved pressure axis",
            },
            "local_core_recovery": {
                "Q": "-M_Z",
                "Q_rho": "-M_rhoZ",
                "Q_Z": "2U-(1-delta)M-M_ZZ",
                "Q_Z_weight": "2-(1-delta)/(n+1) before radial differentiation",
                "shared_source_axis_structure": True,
                "tail_errors_combined_by_triangle_bound": True,
            },
            "coordinate_definition": {
                "rho": "Lambda*R (scaled radial coordinate)",
                "Z": "exactly 0",
                "rho_not_physical_r": True,
            },
            "tail_theorem": {
                "reference_degree": RADIAL_DEGREE,
                "first_omitted_index": RADIAL_DEGREE + 1,
                "reference_radius": TAIL_REFERENCE_RADIUS,
                "source": "accepted common-check geometric_tail_rows degree 24",
                "dilation": "T(rho)<=T(4.1)*(rho/4.1)^(25-i)",
                "old_formula_recomputed": False,
            },
            "G0_amplitude_subcertificate": {
                "method": "new monotone log-H transport",
                "called_once": True,
                "same_source_parameter_and_Z0_gauge_DATA": True,
                "included_in_one_overall_proof_bundle": True,
                "old_gauge_derivatives_recomputed": False,
            },
            "runtime_model_metadata": None,
            "no_old_math_or_provider_calls": True,
        }
        _emit(
            emit,
            events,
            "proof_bundle",
            {
                "proof_kind": proof["kind"],
                "proof_count": 1,
                "same_actual_unique_analytic_core": True,
                "finite_and_tail_parts_same_solution": True,
                "point_value": False,
            },
        )
        _emit(
            emit,
            events,
            "family_completed",
            {
                "proof_count": 1,
                "numeric_evaluation": True,
                "metadata_only": False,
                "source_only": True,
                "full_goal_complete": False,
            },
        )

        preserved = {
            "old_caps_and_errors_unchanged": True,
            "global_errors_retained": True,
            "old_signed_m15_and_error_terms_unchanged": _json_value(
                capsule.get("old_signed_m15_and_error_terms_unchanged")
            ),
            "cumulative_212_terms_unchanged": True,
            "preservation_records": _json_value(
                capsule.get("preservation", capsule.get("preservation_metadata", {}))
            ),
        }
        payload = {
            "component": COMPONENT,
            "source_only": True,
            "numeric_evaluation": True,
            "metadata_only": False,
            "coordinate": "rho=Lambda*R; Z=0",
            "core_amplitude_Z0": amplitude_Z0,
            "source_identity": identity,
            "accepted_function_definition_records": {
                "formal_source_definition": identity["formal_source_definition"],
                "ordinary_order_definition": identity["ordinary_order_definition"],
                "raw_jlogC_definition": identity["raw_jlogC_definition"],
            },
            "admission": admission,
            "packet_identity": packet["packet_identity"],
            "row_units": packet["row_units"],
            "majorant_records": {
                key: _json_value(value) for key, value in majorants.items()
            },
            "norm_records": {
                "Phi_ball_norm_upper": _interval_record(context, norm_data["Phi"]),
                "Psi_ball_norm_upper": _interval_record(context, norm_data["Psi"]),
                "epsilon_core": _interval_record(context, norm_data["epsilon_core"]),
                "epsilon_core_times_Psi_ball_norm_upper": _interval_record(
                    context, norm_data["Uz_M_error_norm"]
                ),
                "delta": _interval_record(context, delta),
                "mixed_pressure_increment_Xh_norm_upper": _interval_record(
                    context, norm_data["P"]
                ),
                "all_are_error_caps_not_values": True,
            },
            "finite_bernstein": finite_bernstein,
            "decasteljau": decasteljau,
            "certified_tail": {
                "reference_radius": TAIL_REFERENCE_RADIUS,
                "degree": RADIAL_DEGREE,
                "selected_rows": {
                    _jet_key(i, k): _json_value(row)
                    for (i, k), row in selected_tail_rows.items()
                },
                "component_norms": {
                    "Phi": _interval_record(context, norm_data["Phi"]),
                    "Uz": _interval_record(context, norm_data["Uz_M_error_norm"]),
                    "Mz_over_R": _interval_record(context, norm_data["Uz_M_error_norm"]),
                    "P_scaled": _interval_record(context, norm_data["P"]),
                },
                "tail_is_error_bound": True,
                "function_values_not_tail_corrected": True,
            },
            "ordinary_jets": ordinary_jets,
            "local_core_recovery": core_recovery,
            "mean_minus_axial": mean_minus_axial,
            "axis_packets": {
                "U0": "fixed.U0_Z_taylor[k] times k!; active only at radial order 0",
                "pressure_axis": "fixed.P0_Z_taylor[k] times k!; active only at radial order 0",
                "pressure_units": "epsilon_core*physical_P",
                "structural_join": True,
                "uncorrelated_box_subtraction": False,
            },
            "exact_function_definition": {
                "Phi": "saved A radial coefficient enclosure of the admitted unique analytic core",
                "Uz": "saved physical Uz radial coefficient enclosure",
                "Mz_over_R": "radial average of Uz, Uz coefficient divided by n+1 before D_rho",
                "P_scaled": "epsilon_core*physical_P with fixed P0 axis row",
                "P": "alias of P_scaled; epsilon_core*physical_P",
                "Q": "-M_Z local core recovery from shared Uz source",
                "Q_rho": "-M_rhoZ local core recovery from shared Uz source",
                "Q_Z": "2U-(1-delta)M-M_ZZ local core recovery",
                "mean_minus_axial": "M-U with weights -n/(n+1), exact shared-axis zero",
                "ordinary_axial_derivative": "k! times saved Taylor coefficient at exact Z=0",
                "radial_derivative": "falling factorial n!/(n-i)! in rho",
            },
            "preservation": preserved,
            "proof_packet": proof,
            "events": events,
            "runtime_model_metadata": None,
            "no_old_math_or_provider_calls": True,
            "no_physical_huge_field_materialized": True,
            "no_NS_claim": True,
        }
        return {
            "schema": SCHEMA,
            "component": COMPONENT,
            "passed": True,
            "metadata_only": False,
            "source_only": True,
            "numeric_evaluation": True,
            "runtime_model_metadata": None,
            "physical_layer_completed": False,
            "full_goal_complete": False,
            "proof_count": 1,
            "actual_counts": dict(SUCCESS_COUNTS),
            "new_flags": dict(SUCCESS_FLAGS),
            "flags": dict(SUCCESS_FLAGS),
            "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
            "count_event_map": {
                key: list(value) for key, value in COUNT_EVENT_MAP.items()
            },
            "payload": payload,
        }
    except KernelInputError:
        # Preserve the callback entry prefix and the root runner's traceback.
        # Never reset counters, reject to zero, or fabricate a failed result.
        raise


__all__ = [
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "CALLBACK_EVENT_COUNTS",
    "EXPECTED_HASH_PINS",
    "EXPECTED_SOURCE_IDENTITY",
    "INTERVAL_DPS",
    "KernelInputError",
    "RADIAL_DEGREE",
    "REQUESTED_JETS",
    "RHO_POINTS",
    "SCHEMA",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]

# New unexecuted G0 transport implementation merged before entry.
"""Unexecuted source-only O3 core amplitude transport staging.

This module deliberately has no imports and no module-level evaluation.  The
root runner supplies the already accepted interval context and the three
read-only JSON packets.  The one new numerical operation is the monotone
log-H change of variables for the same-source quantity ``G(0)``.  In
particular, this file does not select an anchor midpoint, solve an inverse
map, call a source provider, or materialize an astronomical exponential.

The source definitions carried by the packet are

    H(s) = -4*s**3 - j*s**2 + (9-delta)*s/2 + j,
    L(s) = 1-delta*s**2,
    sigma = j/500,
    g(s) = L(s)*H(s)/(H(s)**2+sigma**2).

For the accepted unique root ``a`` of ``H`` in ``[-j,0]``, ``H'`` is
positive on the saved bracket.  With ``w=H/sigma`` and
``u=log(1+w**2)/2``, the exact endpoint values are ``w(a)=0`` and
``w(0)=500``.  Hence the new enclosure is the directed interval product

    G(0) in (L([a,0])/H'([a,0])) * log(250001)/2.

The log-amplitude at the actual source coordinate is then kept in logarithmic
form as ``logF0(0)=-logC-Lambda*G(0)``.  The ordinary Z1/Z2 rows are composed
from the saved accepted Z0 gauge rows; their gauge formulas are not replayed.
"""


G0_SCHEMA = "native_O3_core_G0_log_transport_result_v1"
G0_COMPONENT = "native_O3_core_G0_log_transport"
G0_CAPSULE_SCHEMA = "native_O3_anchor_parameter_axial_jet_kernel_result_v1"

# Exact source constants.  The interval supplied by the caller is used for
# enclosure arithmetic; this literal remains the exact construction identity.
G0_J_NUMERATOR = 1
G0_J_DENOMINATOR = 5721363714816000000000
G0_SIGMA_DENOMINATOR = 500
G0_DISPLAY_DPS = 220
G0_SAVED_RESULT_PATH = "tmp/native_O3_anchor_parameter_axial_jet/result.json"
G0_SAVED_RESULT_SHA256 = "51daffcc06d39b683d67a392bcfe375921c1571a68c1f3968040a4a49ea36b29"
G0_CORE_RHO_VALUES = ("1/16", "1/4", "1", "3")

# These are literal audit contracts, not mutable counters.  The function
# returns a fresh copy so a runner cannot reset or mutate this module state.
G0_STATIC_COUNT_CONTRACT = {
    "N_materialization_calls": 0,
    "failed_family_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "checker_calls": 0,
    "source_function_provider_calls": 0,
    "new_monotone_log_H_transport_family_calls": 1,
    "new_core_amplitude_family_calls": 1,
    "new_mathematical_proofs": 1,
    "numeric_evaluation_calls": 1,
}

G0_STATIC_FLAGS = {
    "source_only": True,
    "actual_same_source_G0_enclosure": True,
    "actual_parameter_value_selected": False,
    "anchor_midpoint_selected": False,
    "inverse_root_or_inverse_transport_solved": False,
    "quadrature_family_used": False,
    "pole_family_used": False,
    "Rouche_family_used": False,
    "old_kernel_called": False,
    "old_provider_called": False,
    "checker_called": False,
    "F0_value_materialized": False,
    "strict_positive_F0_definition_retained": True,
    "full_goal_complete": False,
    "physical_layer_completed": False,
    "actual_parameter_admission": False,
    "full_nonlinear_core_values_evaluated": False,
    "physical_NS_claim": False,
    "temporal_recursion_claim": False,
    "global_five_moment_closure_claim": False,
    "outer_epsilon_used": False,
}


def g0_is_mapping(value):
    return hasattr(value, "get") and hasattr(value, "keys")


def g0_is_sequence(value):
    return isinstance(value, (list, tuple))


def g0_require(mapping, key, label):
    if not g0_is_mapping(mapping):
        raise ValueError(label + " must be a JSON mapping")
    if key not in mapping:
        raise ValueError(label + " is missing " + str(key))
    value = mapping[key]
    if value is None:
        raise ValueError(label + " has null " + str(key))
    return value


def g0_optional(mapping, key, default=None):
    if g0_is_mapping(mapping) and key in mapping:
        return mapping[key]
    return default


def g0_unwrap_parameters(value):
    """Accept the direct parameter payload or its saved result wrappers."""

    if not g0_is_mapping(value):
        raise ValueError("accepted_parameter_payload must be a JSON mapping")
    current = value
    if "result" in current and g0_is_mapping(current["result"]):
        current = current["result"]
    if "payload" in current and g0_is_mapping(current["payload"]):
        current = current["payload"]
    if "parameters" in current and g0_is_mapping(current["parameters"]):
        current = current["parameters"]
    return current


def g0_unwrap_anchor(value):
    """Return the saved root record without choosing an endpoint."""

    if not g0_is_mapping(value):
        raise ValueError("accepted_anchor must be a JSON mapping")
    current = value
    if "result" in current and g0_is_mapping(current["result"]):
        current = current["result"]
    if "payload" in current and g0_is_mapping(current["payload"]):
        current = current["payload"]
    if "anchor" in current and g0_is_mapping(current["anchor"]):
        nested = current["anchor"]
        if "lower" not in nested and "lower_exact_mpf_tuple" not in nested:
            if "anchor" in nested and g0_is_mapping(nested["anchor"]):
                nested = nested["anchor"]
        current = nested
    return current


def g0_unwrap_saved_Z0(value):
    """Return the accepted ordinary gauge rows at the coordinate Z=0."""

    if not g0_is_mapping(value):
        raise ValueError("saved_Z0_gauge must be a JSON mapping")
    current = value
    if "result" in current and g0_is_mapping(current["result"]):
        current = current["result"]
    if "payload" in current and g0_is_mapping(current["payload"]):
        current = current["payload"]
    if "axial_gauge_jets" in current and g0_is_mapping(current["axial_gauge_jets"]):
        current = current["axial_gauge_jets"]
    if "jets" in current and g0_is_mapping(current["jets"]):
        current = current["jets"]
    if "Z_zero" in current and g0_is_mapping(current["Z_zero"]):
        current = current["Z_zero"]
    return current


def g0_safe_json(value):
    """Copy only JSON-shaped provenance; opaque values become text."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if g0_is_mapping(value):
        return {str(key): g0_safe_json(child) for key, child in value.items()}
    if g0_is_sequence(value):
        return [g0_safe_json(child) for child in value]
    return str(value)


def g0_endpoint(context, value, upper=False):
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def g0_lower(context, value):
    return g0_endpoint(context, value, False)


def g0_upper(context, value):
    return g0_endpoint(context, value, True)


def g0_display(context, value):
    text = str(context.nstr(context.convert(value), n=G0_DISPLAY_DPS)).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
        if "," in text:
            text = text.split(",", 1)[0].strip()
    return text


def g0_interval_record(context, value):
    """Serialize a directed interval with display strings and raw MPF tuples."""

    interval = context.convert(value)
    lower_raw, upper_raw = interval._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    width = upper - lower
    return {
        "lower": g0_display(context, lower),
        "upper": g0_display(context, upper),
        "width": g0_display(context, width),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
    }


def g0_interval_from_record(context, value, label):
    """Rehydrate a saved interval without promoting a representative value."""

    if hasattr(value, "_mpi_"):
        return context.convert(value)
    if g0_is_mapping(value):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            lower_raw = value.get("lower_exact_mpf_tuple")
            upper_raw = value.get("upper_exact_mpf_tuple")
            if not g0_is_sequence(lower_raw) or not g0_is_sequence(upper_raw):
                raise ValueError(label + " has malformed raw MPF endpoint tuples")
            if len(lower_raw) != 4 or len(upper_raw) != 4:
                raise ValueError(label + " raw MPF endpoint tuples must have length four")
            left = tuple(int(item) for item in lower_raw)
            right = tuple(int(item) for item in upper_raw)
            return context.make_mpf((left, right))
        if "lower" in value and "upper" in value:
            return context.mpf((str(value["lower"]), str(value["upper"])))
        for nested_key in ("interval", "enclosure", "record", "root_interval", "data"):
            if nested_key in value:
                return g0_interval_from_record(context, value[nested_key], label)
        for scalar_key in ("value", "point", "selected", "scalar"):
            if scalar_key in value:
                return g0_interval_from_record(context, value[scalar_key], label)
        raise ValueError(label + " has no interval endpoints")
    if g0_is_sequence(value):
        if len(value) != 2:
            raise ValueError(label + " endpoint sequence must have length two")
        return context.mpf((str(value[0]), str(value[1])))
    if value is None:
        raise ValueError(label + " is missing")
    return context.mpf(str(value))


def g0_exact_dyadic(context, record, label):
    sign = int(g0_require(record, "sign", label))
    mantissa = int(g0_require(record, "mantissa", label))
    exponent = int(g0_require(record, "binary_exponent", label))
    bitcount = int(g0_require(record, "bitcount", label))
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(label + " has invalid dyadic fields")
    raw = (sign, mantissa, exponent, bitcount)
    return context.make_mpf((raw, raw)), raw


def g0_emit_event(emit, events, name, payload=None):
    row = {
        "event": name,
        "component": G0_COMPONENT,
        "source_only": True,
        "numeric_evaluation": name == "g0_transport_complete",
    }
    if payload is not None:
        row.update(g0_safe_json(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def g0_context_log(context):
    if hasattr(context, "log"):
        return context.log
    if hasattr(context, "ln"):
        return context.ln
    raise ValueError("the supplied interval context has no logarithm operation")


def g0_require_zero_coordinate(context, saved):
    coordinate = g0_optional(saved, "Z")
    if coordinate is None:
        return True
    interval = g0_interval_from_record(context, coordinate, "saved Z0 coordinate")
    zero = context.mpf("0")
    if not bool(g0_lower(context, interval) == zero):
        raise ValueError("saved gauge coordinate is not the accepted Z0 row")
    if not bool(g0_upper(context, interval) == zero):
        raise ValueError("saved gauge coordinate is not the accepted Z0 row")
    return True


def bound_core_amplitude_at_Z0(
    context,
    accepted_parameter_payload,
    accepted_anchor,
    saved_Z0_gauge,
    emit=None,
):
    """Bound the same-source core amplitude and ordinary Z0 jets.

    ``context`` is caller-owned and must be the root runner's directed
    interval context (220 decimal digits in the accepted run).  Every input is
    data-only.  Errors intentionally propagate to the root runner; no local
    counter reset or fallback approximation is permitted.
    """

    events = []
    g0_emit_event(
        emit,
        events,
        "g0_transport_entered",
        {
            "hard_entry": True,
            "method": "monotone_log_H_transport",
            "context_precision_required_dps": 220,
            "source_result_sha256": G0_SAVED_RESULT_SHA256,
        },
    )

    params = g0_unwrap_parameters(accepted_parameter_payload)
    anchor_record = g0_unwrap_anchor(accepted_anchor)
    saved = g0_unwrap_saved_Z0(saved_Z0_gauge)

    # Bind the exact recipe metadata before any interval arithmetic.  The
    # intervals themselves remain enclosures; none is promoted to a source
    # parameter representative.
    j_exact = g0_require(params, "j_exact", "accepted parameters")
    j_num = int(g0_require(j_exact, "numerator", "accepted j_exact"))
    j_den = int(g0_require(j_exact, "denominator", "accepted j_exact"))
    # The fraction is a provenance contract only.  The accepted j enclosure
    # below is the sole arithmetic input; this branch does not re-evaluate the
    # prior minimum/selection recipe.
    if j_num != G0_J_NUMERATOR or j_den != G0_J_DENOMINATOR:
        raise ValueError("accepted j is not the required exact source fraction")
    one = context.mpf("1")
    two = context.mpf("2")
    j_enclosure_record = g0_require(params, "j_enclosure", "accepted parameters")
    sigma_record = g0_require(params, "sigma", "accepted parameters")
    # Both records are accepted directed data; sigma is not reconstructed from
    # j in this consumer.
    j = g0_interval_from_record(context, j_enclosure_record, "accepted j_enclosure")
    sigma = g0_interval_from_record(context, sigma_record, "accepted sigma")

    delta_record = g0_require(params, "delta", "accepted parameters")
    logP_record = g0_require(params, "logP", "accepted parameters")
    logLambda_record = g0_require(params, "logLambda", "accepted parameters")
    Lambda_record = g0_require(params, "Lambda", "accepted parameters")
    logC_dyadic = g0_require(params, "logC_exact_dyadic", "accepted parameters")
    delta = g0_interval_from_record(context, delta_record, "accepted delta")
    logP = g0_interval_from_record(context, logP_record, "accepted logP")
    logLambda = g0_interval_from_record(context, logLambda_record, "accepted logLambda")
    Lambda = g0_interval_from_record(context, Lambda_record, "accepted Lambda")
    logC, logC_raw = g0_exact_dyadic(context, logC_dyadic, "accepted logC_exact_dyadic")
    anchor = g0_interval_from_record(context, anchor_record, "accepted anchor interval")
    anchor_raw = context.convert(anchor)._mpi_
    zero = context.mpf("0")
    zero_raw = context.convert(zero)._mpi_[0]
    if not bool(g0_upper(context, anchor) <= zero):
        raise ValueError("accepted unique root bracket must lie at or below Z0")
    source_bracket = context.make_mpf((anchor_raw[0], zero_raw))

    # Fresh interval formula for H' and L over the saved [a,0] source
    # bracket.  This is a monotonicity certificate, not a root or inverse
    # solve.  The source bracket is a superset of the exact a-to-zero path.
    s = source_bracket
    s_squared = s * s
    half = one / two
    Hprime = -context.mpf("12") * s_squared - context.mpf("2") * j * s
    Hprime = Hprime + (context.mpf("9") - delta) * half
    L = one - delta * s_squared
    if not bool(g0_lower(context, Hprime) > zero):
        raise ValueError("fresh Hprime interval does not certify positivity")
    if not bool(g0_lower(context, L) > zero):
        raise ValueError("fresh L interval does not certify positivity")
    weight = L / Hprime
    if not bool(g0_lower(context, weight) > zero):
        raise ValueError("fresh L/Hprime transport weight is not positive")

    # At w=500 the transformed endpoint is exactly log(1+500^2)/2.  No
    # inverse s(u), root iteration, or numerical quadrature is introduced.
    log_function = g0_context_log(context)
    u_extent = log_function(context.mpf("250001")) / two
    G0 = weight * u_extent
    if not bool(g0_lower(context, u_extent) > zero):
        raise ValueError("log-H endpoint length is not positive")
    if not bool(g0_lower(context, G0) > zero):
        raise ValueError("new same-source G0 interval is not strictly positive")

    logF0_Z0 = -logC - Lambda * G0
    if not bool(g0_upper(context, logF0_Z0) < zero):
        raise ValueError("new logF0(Z0) interval is not strictly negative")

    # The accepted gauge rows are the only inputs for g(0) and g'(0).  Do not
    # reconstruct these rows from H/L or from the root bracket.
    g0_require_zero_coordinate(context, saved)
    saved_g = g0_optional(saved, "g")
    if saved_g is None:
        saved_g = g0_optional(saved, "g_at_Z0")
    if saved_g is None:
        saved_g = g0_optional(saved, "g_zero")
    if saved_g is None:
        raise ValueError("saved accepted Z0 gauge is missing g(Z0) data")
    saved_gprime = g0_optional(saved, "ordinary_DZ_g")
    if saved_gprime is None:
        saved_gprime = g0_optional(saved, "g_prime")
    if saved_gprime is None:
        saved_gprime = g0_optional(saved, "g_prime_at_Z0")
    if saved_gprime is None:
        saved_gprime = g0_optional(saved, "DZ_g")
    if saved_gprime is None:
        raise ValueError("saved accepted Z0 gauge is missing ordinary g-prime data")
    g_at_Z0 = g0_interval_from_record(context, saved_g, "saved g(Z0)")
    gprime_at_Z0 = g0_interval_from_record(context, saved_gprime, "saved g-prime(Z0)")
    saved_DZ_logF0 = g0_interval_from_record(
        context,
        g0_require(saved, "ordinary_DZ_logF0", "saved accepted Z0 gauge"),
        "saved ordinary_DZ_logF0(Z0)",
    )
    saved_DZZ_logF0 = g0_interval_from_record(
        context,
        g0_require(saved, "ordinary_DZZ_logF0", "saved accepted Z0 gauge"),
        "saved ordinary_DZZ_logF0(Z0)",
    )
    saved_F0_Z_over_F0 = g0_interval_from_record(
        context,
        g0_require(saved, "F0_Z_over_F0", "saved accepted Z0 gauge"),
        "saved F0_Z_over_F0(Z0)",
    )
    saved_F0_ZZ_over_F0 = g0_interval_from_record(
        context,
        g0_require(saved, "F0_ZZ_over_F0", "saved accepted Z0 gauge"),
        "saved F0_ZZ_over_F0(Z0)",
    )
    # These four ordinary derivative/relative-amplitude rows are accepted
    # numerical data.  They are transported verbatim; no gauge formula or
    # Bell postcomposition is replayed in this new G0 family.
    ell1 = saved_DZ_logF0
    ell2 = saved_DZZ_logF0
    F0_Z_over_F0 = saved_F0_Z_over_F0
    F0_ZZ_over_F0 = saved_F0_ZZ_over_F0

    # Four saved core radial labels all use this one Z0 log-amplitude.  They
    # are labels only; no rho-dependent amplitude calculation is performed.
    rho_log_rows = []
    for rho in G0_CORE_RHO_VALUES:
        rho_log_rows.append(
            {
                "rho": rho,
                "coordinate_namespace": "core_rho",
                "logF0_Z0": g0_interval_record(context, logF0_Z0),
                "reuse": "the same function of Z0; no rho-dependent recomputation",
            }
        )

    proof_packet = {
        "method": "monotone_log_H_transport",
        "variable_change": [
            "w(s)=H(s)/sigma",
            "u(s)=log(1+w(s)^2)/2",
            "du=H(s)*Hprime(s)/(H(s)^2+sigma^2) ds",
            "g(s) ds=(L(s)/Hprime(s)) du",
            "H(a)=0 and H(0)=j with sigma=j/500 give w(a)=0 and w(0)=500 exactly",
            "u(0)-u(a)=log(1+500^2)/2=log(250001)/2",
            "G0=integral_0^(log(250001)/2) L(s(u))/Hprime(s(u)) du",
        ],
        "Hprime_formula": "-12*s^2-2*j*s+(9-delta)/2",
        "L_formula": "1-delta*s^2",
        "sigma_formula": "j/500",
        "j_sigma_exact_ratio": "j/sigma=500 exactly at the same exact j",
        "w_endpoint_values": {"w_at_anchor": "0 exactly", "w_at_Z0": "500 exactly"},
        "G0_variation_scale": "narrow same-source variation O(j^2); no midpoint is selected",
        "inverse_map_used": False,
        "root_solver_used": False,
        "quadrature_used": False,
        "Hprime_positive_on_saved_source_bracket": True,
        "positive_bounds": {
            "u_extent_gt_0": True,
            "transport_weight_gt_0": True,
            "G0_gt_0": True,
            "logF0_Z0_lt_0": True,
        },
        "Hprime_bracket": g0_interval_record(context, Hprime),
        "L_bracket": g0_interval_record(context, L),
        "transport_weight_L_over_Hprime": g0_interval_record(context, weight),
        "u_extent": g0_interval_record(context, u_extent),
        "G0": g0_interval_record(context, G0),
        "anchor_bracket_used_as_domain_only": True,
        "anchor_value_selected": False,
    }
    g0_emit_event(emit, events, "g0_proof_bundle", proof_packet)

    result = {
        "schema": G0_SCHEMA,
        "component": G0_COMPONENT,
        "source_only": True,
        "passed": True,
        "numeric_evaluation": True,
        "proof_count": 1,
        "actual_source_value_enclosure": True,
        "actual_source_Z0_value_enclosure": True,
        "ordinary_Z0_gauge_consumed": True,
        "conditional_operator_only": False,
        "full_goal_complete": False,
        "physical_layer_completed": False,
        "events": events,
        "hard_entry_event": "g0_transport_entered",
        "actual_counts": dict(G0_STATIC_COUNT_CONTRACT),
        "flags": dict(G0_STATIC_FLAGS),
        "source_data": {
            "saved_result_path": G0_SAVED_RESULT_PATH,
            "saved_result_sha256": G0_SAVED_RESULT_SHA256,
            "saved_Z0_gauge_data_id": G0_SAVED_RESULT_PATH + "#result.payload.axial_gauge_jets.jets.Z_zero",
            "accepted_parameter_payload_used_as_readonly": True,
            "accepted_anchor_used_as_readonly": True,
            "raw_mpf_tuples_are_certificates_only": True,
            "decimal_endpoint_displays_are_not_selected_values": True,
        },
        "exact_parameter_recipes": {
            "j": "1/5721363714816000000000",
            "j_exact_metadata": g0_safe_json(j_exact),
            "logP": "exp(40)+11",
            "delta": "exp(-4*logP-30)",
            "logLambda": "4*logP+1000",
            "Lambda": "exp(logLambda)",
            "epsilon_core": "exp(-logLambda)",
            "logC": "the selected singleton dyadic supplied by accepted_parameter_payload",
            "logC_selected_dyadic_raw_mpf_tuple": list(logC_raw),
            "Lambda_interval_source": "accepted_parameter_payload.Lambda; exp(logLambda) is not materialized here",
        },
        "parameter_intervals": {
            "j": g0_interval_record(context, j),
            "sigma": g0_interval_record(context, sigma),
            "logP": g0_interval_record(context, logP),
            "delta": g0_interval_record(context, delta),
            "logLambda": g0_interval_record(context, logLambda),
            "Lambda": g0_interval_record(context, Lambda),
            "selected_logC": g0_interval_record(context, logC),
            "core_epsilon_namespace": "core_epsilon_only; no outer_epsilon is used",
        },
        "accepted_anchor": {
            "root_interval": g0_interval_record(context, anchor),
            "source_bracket": g0_interval_record(context, source_bracket),
            "definition": "the unique exact real root H(a)=0 in [-j,0]",
            "root_interval_is_not_a_selected_scalar": True,
        },
        "transport": proof_packet,
        "amplitude": {
            "Z0_coordinate": "Z=0",
            "G0": g0_interval_record(context, G0),
            "logF0_Z0": g0_interval_record(context, logF0_Z0),
            "logF0_formula": "-logC-Lambda*G0",
            "F0_exact_definition": "exp(logF0_Z0)",
            "strictly_positive_F0": True,
            "F0_value_materialized": False,
            "astronomical_exponential_not_materialized": True,
        },
        "ordinary_Z0_gauge": {
            "saved_data_id": G0_SAVED_RESULT_PATH + "#result.payload.axial_gauge_jets.jets.Z_zero",
            "saved_g_Z0": g0_interval_record(context, g_at_Z0),
            "saved_g_prime_Z0": g0_interval_record(context, gprime_at_Z0),
            "saved_ordinary_DZ_logF0": g0_interval_record(context, saved_DZ_logF0),
            "saved_ordinary_DZZ_logF0": g0_interval_record(context, saved_DZZ_logF0),
            "saved_F0_Z_over_F0": g0_interval_record(context, saved_F0_Z_over_F0),
            "saved_F0_ZZ_over_F0": g0_interval_record(context, saved_F0_ZZ_over_F0),
            "saved_rows_reused_without_gauge_formula_replay": True,
            "saved_relative_rows_reused_without_Bell_postcomposition": True,
            "logF0_Z1_formula": "saved ordinary_DZ_logF0 at Z0",
            "logF0_Z2_formula": "saved ordinary_DZZ_logF0 at Z0",
            "logF0_Z1": g0_interval_record(context, ell1),
            "logF0_Z2": g0_interval_record(context, ell2),
            "F0_Z_over_F0": g0_interval_record(context, F0_Z_over_F0),
            "F0_ZZ_over_F0": g0_interval_record(context, F0_ZZ_over_F0),
            "F0_ZZ_over_F0_formula": "saved F0_ZZ_over_F0 at Z0",
            "relative_rows_recomputed": False,
            "coordinate_namespace": {
                "Z0": "value at coordinate Z=0",
                "Z1": "first ordinary derivative at coordinate Z=0",
                "Z2": "second ordinary derivative at coordinate Z=0",
                "Z1_is_not_coordinate_Z_equals_1": True,
                "Z2_is_not_coordinate_Z_equals_2": True,
            },
        },
        "logF0_by_core_rho": rho_log_rows,
        "core_rho_values": list(G0_CORE_RHO_VALUES),
        "core_rho_reuse_rule": "all four rows reuse exactly the Z0 log-amplitude interval",
        "negative_scope": {
            "physical_NS": False,
            "temporal_n_recursion": False,
            "global_five_moment_closure": False,
            "outer_layer_epsilon": False,
        },
    }
    g0_emit_event(
        emit,
        events,
        "g0_transport_complete",
        {
            "actual_source_value_enclosure": True,
            "G0_gt_0": True,
            "logF0_Z0_lt_0": True,
            "F0_value_materialized": False,
        },
    )
    return result
