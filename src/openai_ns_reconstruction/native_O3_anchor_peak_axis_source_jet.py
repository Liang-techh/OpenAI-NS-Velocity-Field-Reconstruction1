"""Unexecuted native O3 anchor-axis source jet.

This leaf is deliberately independent of the entered O3 kernels. It accepts
one frozen data["capsule"] wrapper and, when the root invokes it, evaluates
only the exact-root local source scale

    Z = a + s*eta,       s = sigma*sqrt(epsilon_core),

using the correlated identity H(a)=0. It does not import a project module,
old G0/log-H transport, root solver, recurrence, provider, quadrature, or
checker. The absolute F0 exponential is never materialized. The capsule
canonical digest is computed from the capsule body; a supplied digest is not
used as a substitute for ROOT_CAPSULE_CANONICAL.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_anchor_peak_axis_source_jet_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_anchor_peak_axis_source_jet_source_capsule_v1"
COMPONENT = "native_O3_anchor_peak_axis_source_jet"
INTERVAL_DPS = 280
SOURCE_PRECISION_BITS = 867

# The root freezes this literal after constructing the capsule. An empty
# value is intentionally a hard rejection, not a permissive placeholder.
ROOT_CAPSULE_CANONICAL = "506ada9f2916ac74b9af523bf3667bbaf435301082d3f7ce4ffab26496b277ec"

ANCHOR_SOURCE_SHA256 = (
    "51daffcc06d39b683d67a392bcfe375921c1571a68c1f3968040a4a49ea36b29"
)
PRESSURE_SOURCE_SHA256 = (
    "909ad64901823b68beae8ef5f95abdfd0f90fc9de509460928b7283ebb7b9d99"
)
SOURCE_IDS = {
    "anchor_parameter_result_sha256": ANCHOR_SOURCE_SHA256,
    "pressure_result_sha256": PRESSURE_SOURCE_SHA256,
}

ETA_TEXTS = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
FAMILY_NAMES = (
    "anchor_polynomial",
    "scaled_gauge_integral",
    "relative_profile_jets",
    "axis_pressure_envelope",
    "peak_width_comparison",
)
CALLBACK_EVENTS = (
    "family_entered",
    "anchor_polynomial_entered",
    "anchor_polynomial_complete",
    "scaled_gauge_integral_entered",
    "scaled_gauge_integral_complete",
    "relative_profile_jets_entered",
    "relative_profile_jets_complete",
    "axis_pressure_envelope_entered",
    "axis_pressure_envelope_complete",
    "peak_width_comparison_entered",
    "peak_width_comparison_complete",
    "proof_bundle",
    "family_completed",
)
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "anchor_polynomial_entered": 1,
    "anchor_polynomial_complete": 1,
    "scaled_gauge_integral_entered": 1,
    "scaled_gauge_integral_complete": 1,
    "relative_profile_jets_entered": 1,
    "relative_profile_jets_complete": 1,
    "axis_pressure_envelope_entered": 1,
    "axis_pressure_envelope_complete": 1,
    "peak_width_comparison_entered": 1,
    "peak_width_comparison_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}
MANAGER_JOURNAL_EVENT_COUNT = 15

# This is one derive/source entry with five new mathematical families, not
# five independent source/provider calls.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 5,
    "anchor_polynomial_family_calls": 1,
    "scaled_gauge_integral_family_calls": 1,
    "relative_profile_jets_family_calls": 1,
    "axis_pressure_envelope_family_calls": 1,
    "peak_width_comparison_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
}
SUCCESS_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "local_anchor_axis_ports": True,
    "anchor_root_correlated": True,
    "relative_profile_log_packet": True,
    "relative_profile_exp_safe_only": True,
    "axis_pressure_envelope_source_defined": True,
    "axis_pressure_envelope_value_admitted": True,
    "nonlinear_off_axis_core_values_admitted": False,
    "full_physical": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
    "N_materialized": False,
    "F0_absolute_exp_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "anchor_midpoint_selected": False,
    "Gaussian_replacement_used": False,
}
COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "anchor_polynomial_entered",
        "scaled_gauge_integral_entered",
        "relative_profile_jets_entered",
        "axis_pressure_envelope_entered",
        "peak_width_comparison_entered",
    ],
    "anchor_polynomial_family_calls": ["anchor_polynomial_entered"],
    "scaled_gauge_integral_family_calls": ["scaled_gauge_integral_entered"],
    "relative_profile_jets_family_calls": ["relative_profile_jets_entered"],
    "axis_pressure_envelope_family_calls": ["axis_pressure_envelope_entered"],
    "peak_width_comparison_family_calls": ["peak_width_comparison_entered"],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}

# Static literal contracts for the runner and for the persisted metadata.
EXPECTED_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 5,
    "anchor_polynomial_family_calls": 1,
    "scaled_gauge_integral_family_calls": 1,
    "relative_profile_jets_family_calls": 1,
    "axis_pressure_envelope_family_calls": 1,
    "peak_width_comparison_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "source_value_reducer_calls": 0,
    "F0_absolute_exp_calls": 0,
}
EXPECTED_FLAGS = {
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "local_anchor_axis_ports": True,
    "anchor_root_correlated": True,
    "relative_profile_log_packet": True,
    "relative_profile_exp_safe_only": True,
    "axis_pressure_envelope_source_defined": True,
    "axis_pressure_envelope_value_admitted": True,
    "nonlinear_off_axis_core_values_admitted": False,
    "full_physical": False,
    "physical_layer_completed": False,
    "global_five_moment": False,
    "time_n_recursion": False,
    "full_goal": False,
    "full_goal_complete": False,
    "N_materialized": False,
    "F0_absolute_exp_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "provider_called": False,
    "checker_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "anchor_midpoint_selected": False,
    "Gaussian_replacement_used": False,
}
API = {
    "name": "derive",
    "signature": "derive(data, emit=None)",
    "capsule_selector": "data['capsule']",
    "result_schema": "native_O3_anchor_peak_axis_source_jet_kernel_result_v1",
    "source_only": True,
    "numeric": True,
    "metadata_only": False,
}
API_METADATA = {
    "function": "derive",
    "signature": "derive(data, emit=None)",
    "input": "data['capsule']",
    "metadata_only": False,
    "source_only": True,
    "numeric_evaluation": True,
    "runtime_model_metadata": None,
}

class KernelInputError(ValueError):
    """Raised for a missing or unauthenticated source capsule field."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _json_value(value: Any) -> Any:
    """Copy JSON-shaped metadata without interpreting opaque source objects."""

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
    row: dict[str, Any] = {
        "event": name,
        "component": COMPONENT,
        "source_only": True,
        "numeric_evaluation": name.endswith("_complete")
        or name in {"family_entered", "proof_bundle", "family_completed"},
    }
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def _canonical_bytes(capsule: Mapping[str, Any]) -> bytes:
    # The frozen capsule has no canonical field. Hash the entire deterministic
    # JSON object; the root keeps the resulting digest in the code constant.
    body = _json_value(capsule)
    return (
        json.dumps(
            body,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        .encode("utf-8")
    )


def _valid_sha(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdefABCDEF" for char in value)
    )


def _lookup(mapping: Mapping[str, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _exact_mpf_from_tuple(context: Any, value: Any, label: str) -> Any:
    if not _is_sequence(value) or len(value) != 4:
        raise KernelInputError(label + " must be a four-word MPF tuple")
    try:
        raw = tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise KernelInputError(label + " contains a non-integer MPF word") from exc
    return context.make_mpf(raw)


def _interval_value(context: Any, value: Any, label: str) -> Any:
    """Read an exact scalar or directed interval without choosing a point."""

    if _is_mapping(value):
        if "interval" in value:
            return _interval_value(context, value["interval"], label + ".interval")
        if "enclosure" in value:
            return _interval_value(context, value["enclosure"], label + ".enclosure")
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lower = _exact_mpf_from_tuple(
                context, value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple"
            )
            upper = _exact_mpf_from_tuple(
                context, value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple"
            )
            return context.mpf([lower, upper])
        if "lower" in value and "upper" in value:
            return context.mpf([str(value["lower"]), str(value["upper"])])
        if "numerator" in value and "denominator" in value:
            denominator = int(value["denominator"])
            if denominator == 0:
                raise KernelInputError(label + " has zero denominator")
            return context.mpf(int(value["numerator"])) / context.mpf(denominator)
        if "value" in value:
            return _interval_value(context, value["value"], label + ".value")
        if "raw" in value:
            return _interval_value(context, value["raw"], label + ".raw")
        raise KernelInputError(label + " has no scalar or interval record")
    if isinstance(value, bool) or value is None:
        raise KernelInputError(label + " must not be boolean or null")
    if isinstance(value, (str, int, float)):
        return context.mpf(str(value))
    raise KernelInputError(label + " is not a scalar or interval record")


def _raw_endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, False)


def _upper(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, True)


def _display_endpoint(context: Any, value: Any) -> str:
    text = str(context.nstr(context.convert(value), n=INTERVAL_DPS)).strip()
    if text.startswith("[") and text.endswith("]") and "," in text:
        text = text[1:-1].split(",", 1)[0].strip()
    return text


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    interval = context.convert(value)
    lower_raw, upper_raw = interval._mpi_
    lower = context.make_mpf((lower_raw, lower_raw))
    upper = context.make_mpf((upper_raw, upper_raw))
    return {
        "lower": _display_endpoint(context, lower),
        "upper": _display_endpoint(context, upper),
        "width": _display_endpoint(context, upper - lower),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
        "point_value": False,
    }


def _strict_positive(context: Any, value: Any) -> bool:
    return bool(_lower(context, value) > context.mpf(0))


def _exact_eta(context: Any, text: str) -> Any:
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return context.mpf(int(numerator)) / context.mpf(int(denominator))
    return context.mpf(text)


def _eta_interval(context: Any, eta: Any) -> Any:
    lower = _lower(context, eta)
    upper = _upper(context, eta)
    if bool(upper < context.mpf(0)):
        return context.mpf([lower, context.mpf(0)])
    if bool(lower > context.mpf(0)):
        return context.mpf([context.mpf(0), upper])
    return context.mpf([context.mpf(0), context.mpf(0)])


def _lookup_parameter(
    parameters: Mapping[str, Any],
    names: Sequence[str],
) -> Any:
    raw = _lookup(parameters, names)
    if raw is None:
        intervals = parameters.get("parameter_intervals")
        if _is_mapping(intervals):
            raw = _lookup(intervals, names)
    return raw


def _parameter_value(
    context: Any,
    parameters: Mapping[str, Any],
    names: Sequence[str],
    label: str,
) -> tuple[Any, Any]:
    raw = _lookup_parameter(parameters, names)
    if raw is None:
        raise KernelInputError("missing source parameter " + label)
    return raw, _interval_value(context, raw, "parameters." + label)


def _pressure_selectors(pressure: Mapping[str, Any]) -> Any:
    selectors = _lookup(
        pressure,
        (
            "selectors",
            "selector",
            "source_selector",
            "pressure_source_ordinary_Z_jets",
        ),
    )
    if selectors is None and "zero_row" in pressure and "anchor_row" in pressure:
        # These are literal capsule-owned DATA paths, not a synthetic value
        # selector and not a source-point substitution.
        selectors = {
            "zero_row": "pressure_source.zero_row",
            "zero_pressure": "pressure_source.zero_row.pressure",
            "zero_D_ZZ_pressure": "pressure_source.zero_row.D_ZZ_pressure",
            "anchor_row": "pressure_source.anchor_row",
            "anchor_pressure": "pressure_source.anchor_row.pressure",
            "anchor_D_Z_pressure": "pressure_source.anchor_row.D_Z_pressure",
            "anchor_D_ZZ_pressure": "pressure_source.anchor_row.D_ZZ_pressure",
        }
    return selectors


def _validate_capsule(
    capsule: Any,
) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]:
    capsule = _as_mapping(capsule, "data.capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected anchor-peak source capsule schema")
    if capsule.get("component") != COMPONENT:
        raise KernelInputError("source capsule component does not match")
    if not _valid_sha(ROOT_CAPSULE_CANONICAL) or not ROOT_CAPSULE_CANONICAL.strip("0"):
        raise KernelInputError("ROOT_CAPSULE_CANONICAL is not frozen")
    actual = hashlib.sha256(_canonical_bytes(capsule)).hexdigest()
    if actual.lower() != ROOT_CAPSULE_CANONICAL.lower():
        raise KernelInputError("source capsule canonical digest mismatch")
    grid = _as_mapping(capsule.get("grid"), "capsule.grid")
    eta_grid = grid.get("eta")
    tau_grid = grid.get("physical_tau")
    if not _is_sequence(eta_grid) or list(eta_grid) != list(ETA_TEXTS):
        raise KernelInputError("capsule eta grid does not match the frozen source grid")
    if not _is_sequence(tau_grid) or list(tau_grid) != ["1", "1/4", "1/16", "1/64"]:
        raise KernelInputError("capsule physical_tau grid does not match the frozen source grid")
    guards = _as_mapping(capsule.get("guards"), "capsule.guards")
    required_true_guards = (
        "H_anchor_identity_retained",
        "Lambda_epsilon_same_source_exact_identity",
        "relative_profile_is_not_absolute_amplitude",
        "reserved_contribution_once_retained",
        "all_212_terms_retained",
        "all_global_errors_and_caps_retained",
        "six_signed_channels_retained",
        "source_input_is_saved_DATA",
    )
    for guard_name in required_true_guards:
        if guards.get(guard_name) is not True:
            raise KernelInputError("required source guard is not true: " + guard_name)
    required_false_guards = (
        "full_goal_complete",
        "full_physical_field_materialized",
        "nonlinear_off_axis_core_values_admitted",
        "physical_layer_completed",
        "root_midpoint_selected",
    )
    for guard_name in required_false_guards:
        if guards.get(guard_name) is not False:
            raise KernelInputError("forbidden source guard is not false: " + guard_name)
    for guard_name in (
        "N_materialization_calls",
        "absolute_F0_exponential_materialization_calls",
        "old_math_calls",
        "old_provider_calls",
        "old_quadrature_calls",
        "old_recurrence_calls",
    ):
        if guards.get(guard_name) != 0:
            raise KernelInputError("source guard count is not zero: " + guard_name)
    source_ids = _as_mapping(capsule.get("source_ids"), "capsule.source_ids")
    for name, expected in SOURCE_IDS.items():
        if source_ids.get(name) != expected:
            raise KernelInputError("source id mismatch: " + name)

    parameters = _as_mapping(capsule.get("parameters"), "capsule.parameters")
    anchor = _as_mapping(capsule.get("anchor"), "capsule.anchor")
    pressure = _as_mapping(capsule.get("pressure_source"), "capsule.pressure_source")
    exact_identities = _as_mapping(
        capsule.get("exact_identities", {}),
        "capsule.exact_identities",
    )
    identity = _lookup(
        exact_identities,
        ("Lambda_epsilon_equals_one", "Lambda_epsilon_exact_identity"),
    )
    if identity is None:
        identity = _lookup(
            parameters,
            ("Lambda_epsilon_exact_identity", "Lambda_epsilon_equals_one"),
        )
    if identity is None:
        identity = _lookup(
            capsule,
            ("Lambda_epsilon_exact_identity", "Lambda_epsilon_equals_one"),
        )
    if identity is None:
        raise KernelInputError("exact Lambda*epsilon=1 source identity is missing")

    definition = str(
        _lookup(anchor, ("exact_definition", "definition", "root_definition")) or ""
    )
    if "H(a)=0" not in definition.replace(" ", ""):
        raise KernelInputError("anchor exact H(a)=0 definition is missing")
    if bool(anchor.get("root_interval_midpoint_selected", False)):
        raise KernelInputError("anchor midpoint selection is forbidden")
    if bool(anchor.get("nominal_anchor_selected", False)):
        raise KernelInputError("nominal anchor selection is forbidden")

    selectors = _pressure_selectors(pressure)
    if selectors is None or selectors == {} or selectors == "":
        raise KernelInputError("pressure source selectors are required")
    return capsule, parameters, anchor, pressure


def _new_interval_context() -> Any:
    # Deliberately local: importing an interval library is permitted only when
    # the root chooses to invoke this new leaf. No repository module is read.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _parameter_packet(
    context: Any,
    parameters: Mapping[str, Any],
) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for label, names in (
        ("j", ("j", "j_enclosure", "j_exact")),
        ("delta", ("delta", "delta_enclosure")),
        ("epsilon_core", ("epsilon_core", "epsilon", "epsilon_enclosure")),
        ("Lambda", ("Lambda", "Lambda_enclosure")),
        ("logC", ("logC", "logC_enclosure", "logCstar", "selected_logCstar")),
        ("sigma", ("sigma", "sigma_enclosure")),
    ):
        raw, interval = _parameter_value(context, parameters, names, label)
        rows[label] = {"raw": _json_value(raw), "interval": interval}
    j = rows["j"]["interval"]
    delta = rows["delta"]["interval"]
    epsilon = rows["epsilon_core"]["interval"]
    Lambda = rows["Lambda"]["interval"]
    sigma = rows["sigma"]["interval"]
    if not _strict_positive(context, j):
        raise KernelInputError("j must be strictly positive")
    if not _strict_positive(context, epsilon):
        raise KernelInputError("epsilon_core must be strictly positive")
    if not _strict_positive(context, Lambda):
        raise KernelInputError("Lambda must be strictly positive")
    if not _strict_positive(context, delta):
        raise KernelInputError("delta must be strictly positive")
    if not bool(_upper(context, delta) < context.mpf(1) / 200):
        raise KernelInputError("delta must be below 1/200")
    if not _strict_positive(context, sigma):
        raise KernelInputError("sigma=j/500 must be strictly positive")
    scale = sigma * context.sqrt(epsilon)
    if not _strict_positive(context, scale):
        raise KernelInputError("source peak scale must be strictly positive")
    rows["sigma"]["definition"] = (
        "saved sigma DATA; same-source exact recipe sigma=j/500"
    )
    rows["scale_s"] = {
        "raw": "sigma*sqrt(epsilon_core)",
        "interval": scale,
        "definition": "s=sigma*sqrt(epsilon_core)=sigma/sqrt(Lambda)",
    }
    rows["Lambda_epsilon_identity"] = _lookup(
        parameters,
        ("Lambda_epsilon_exact_identity", "Lambda_epsilon_equals_one"),
    )
    rows["parameter_point_selection"] = False
    return rows


def _anchor_packet(context: Any, anchor: Mapping[str, Any]) -> dict[str, Any]:
    raw = _lookup(
        anchor, ("interval", "root_interval", "anchor_interval", "anchor")
    )
    if raw is None:
        raw = anchor
    value = _interval_value(context, raw, "anchor")
    return {
        "interval": value,
        "record": _interval_record(context, value),
        "exact_definition": str(
            _lookup(anchor, ("exact_definition", "definition", "root_definition"))
            or ""
        ),
        "root_tag": "exact implicit root a of H(a)=0",
        "midpoint_selected": False,
        "H_at_anchor_used_as_exact_identity": True,
    }


def _source_polynomial(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
    eta: Any,
) -> dict[str, Any]:
    a = anchor["interval"]
    j = packet["j"]["interval"]
    delta = packet["delta"]["interval"]
    scale = packet["scale_s"]["interval"]
    h = scale * eta
    H1 = (context.mpf(9) - delta) / 2 - 12 * a**2 - 2 * j * a
    L0 = 1 - delta * a**2
    K = H1 - (12 * a + j) * h - 4 * h**2
    Lh = L0 - 2 * delta * a * h - delta * h**2
    z = a + h
    return {
        "eta": eta,
        "h": h,
        "Z": z,
        "H1": H1,
        "L0": L0,
        "K": K,
        "Lh": Lh,
        "H_factorization": "H(a+h)=h*K(h), exact H(a)=0",
        "L_factorization": "L(a+h)=L0-2*delta*a*h-delta*h^2",
    }


def _B_range(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
    xi_range: Any,
) -> dict[str, Any]:
    a = anchor["interval"]
    j = packet["j"]["interval"]
    delta = packet["delta"]["interval"]
    epsilon = packet["epsilon_core"]["interval"]
    scale = packet["scale_s"]["interval"]
    H1 = (context.mpf(9) - delta) / 2 - 12 * a**2 - 2 * j * a
    L0 = 1 - delta * a**2
    h = scale * xi_range
    K = H1 - (12 * a + j) * h - 4 * h**2
    Lh = L0 - 2 * delta * a * h - delta * h**2
    den = 1 + epsilon * xi_range**2 * K**2
    B = Lh * K / den
    if not _strict_positive(context, K):
        raise KernelInputError("K is not strictly positive on the signed eta segment")
    if not _strict_positive(context, Lh):
        raise KernelInputError("L is not strictly positive on the signed eta segment")
    if not _strict_positive(context, den):
        raise KernelInputError("scaled denominator is not strictly positive")
    return {
        "xi_range": xi_range,
        "h_range": h,
        "K_range": K,
        "L_range": Lh,
        "denominator_range": den,
        "B": B,
        "positive_weight_theorem": (
            "xi=eta*t gives xi dxi=eta^2*t dt; B is enclosed over the entire "
            "signed segment, so D is eta^2/2 times the directed B enclosure."
        ),
    }


def _D_for_eta(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
    eta: Any,
) -> tuple[Any, dict[str, Any]]:
    xi_range = _eta_interval(context, eta)
    bpacket = _B_range(context, packet, anchor, xi_range)
    D = (eta**2 / 2) * bpacket["B"]
    return D, bpacket


def _safe_relative_profile(context: Any, D: Any) -> dict[str, Any]:
    zero = context.mpf(0)
    ten = context.mpf(10)
    lower = _lower(context, D)
    upper = _upper(context, D)
    if bool(lower < zero) or bool(upper > ten):
        raise KernelInputError(
            "relative D interval is not certified inside [0,10]"
        )
    return {
        "status": "safe_relative_only",
        "D": _interval_record(context, D),
        "relative_profile": _interval_record(context, context.exp(-D)),
        "definition": "F0(a+s*eta)/F0(a)=exp(-D(eta))",
        "exp_materialized": True,
        "absolute_F0_exp_materialized": False,
    }


def _relative_jet(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
    eta: Any,
    D: Any,
) -> dict[str, Any]:
    source = _source_polynomial(context, packet, anchor, eta)
    a = anchor["interval"]
    j = packet["j"]["interval"]
    delta = packet["delta"]["interval"]
    epsilon = packet["epsilon_core"]["interval"]
    scale = packet["scale_s"]["interval"]
    h = source["h"]
    K = source["K"]
    Lh = source["Lh"]
    den = 1 + epsilon * eta**2 * K**2
    K_eta = -(12 * a + j) * scale - 8 * h * scale
    L_eta = -2 * delta * (a + h) * scale
    den_eta = (
        2 * epsilon * eta * K**2
        + 2 * epsilon * eta**2 * K * K_eta
    )
    B = Lh * K / den
    B_eta = (
        (L_eta * K + Lh * K_eta) / den
        - Lh * K * den_eta / den**2
    )
    D_eta = eta * B
    D_etaeta = B + eta * B_eta
    numerator_F0ZZ = D_eta**2 - D_etaeta
    DZ_logF0 = -D_eta / scale
    DZZ_logF0 = -D_etaeta / scale**2
    return {
        "eta": _interval_record(context, eta),
        "Z_formal": "a+s*eta",
        "Z_interval": _interval_record(context, source["Z"]),
        "D": _interval_record(context, D),
        "B_at_eta": _interval_record(context, B),
        "B_eta": _interval_record(context, B_eta),
        "D_eta": _interval_record(context, D_eta),
        "D_etaeta": _interval_record(context, D_etaeta),
        "ordinary_DZ_logF0": _interval_record(context, DZ_logF0),
        "ordinary_DZZ_logF0": _interval_record(context, DZZ_logF0),
        "F0Z_over_F0": _interval_record(context, DZ_logF0),
        "F0ZZ_over_F0": _interval_record(context, numerator_F0ZZ / scale**2),
        "F0ZZ_numerator_before_scale": _interval_record(context, numerator_F0ZZ),
        "eta_zero_first_order_exact": bool(
            _lower(context, eta) == context.mpf(0)
            and _upper(context, eta) == context.mpf(0)
        ),
        "absolute_logF0": {
            "definition": "logF0(a+s*eta)=-logC-D(eta)",
            "record": _interval_record(
                context, -packet["logC"]["interval"] - D
            ),
        },
        "absolute_F0": {
            "definition": "F0(a+s*eta)=exp(-logC-D(eta))",
            "strictly_positive_formal_atom": True,
            "materialized": False,
        },
    }


def _axis_source_rows(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
    eta: Any,
) -> dict[str, Any]:
    source = _source_polynomial(context, packet, anchor, eta)
    a = anchor["interval"]
    j = packet["j"]["interval"]
    delta = packet["delta"]["interval"]
    h = source["h"]
    z = source["Z"]
    uz0_at_exact_anchor = -((1 - delta) / 2) * a / (1 - a**2)
    uz0 = uz0_at_exact_anchor + 4 * h
    d = 1 - z**2
    L = 1 - delta * z**2
    q0 = ((1 + delta) * z * uz0 - 4 * d) / L
    d_uz0 = h * source["K"] - (1 - delta) * (a + h) / 2
    d_uz0_z = (
        source["H1"]
        - 2 * (12 * a + j) * h
        - 12 * h**2
        - (1 - delta) / 2
    )
    d_uz0_zz = -2 * j - 24 * (a + h)
    return {
        "eta": _interval_record(context, eta),
        "Z": _interval_record(context, source["Z"]),
        "Uz0_at_exact_anchor": _interval_record(context, uz0_at_exact_anchor),
        "Uz0": _interval_record(context, uz0),
        "Uz0_Z": _interval_record(context, context.mpf(4)),
        "Uz0_ZZ": _interval_record(context, context.mpf(0)),
        "q0": _interval_record(context, q0),
        "d_times_Uz0": _interval_record(context, d_uz0),
        "d_times_Uz0_Z": _interval_record(context, d_uz0_z),
        "d_times_Uz0_ZZ": _interval_record(context, d_uz0_zz),
        "definitions": {
            "Uz0": "Uz0(Z)=4Z+j",
            "Uz0_at_exact_anchor": (
                "U0(a)=-((1-delta)/2)*a/(1-a^2) from exact H(a)=0"
            ),
            "Uz0_anchor_affine_identity": (
                "U0(a)+4h=4(a+h)+j; same H(a)=0 correlation and "
                "H=(1-Z^2)*U0+(1-delta)*Z/2"
            ),
            "q0": "q0=[(1+delta)*Z*Uz0-(1-Z^2)*4]/L",
            "d_times_Uz0": "d*Uz0=H(Z)-(1-delta)*Z/2",
            "exact_root_factor": "H(a+h)=h*K(h)",
        },
        "point_value": False,
    }


def _axis_binding(capsule: Mapping[str, Any]) -> Mapping[str, Any]:
    axis = capsule.get("axis_source")
    if axis is None:
        axis = capsule.get("source_bindings")
    if axis is None:
        axis = capsule.get("source_recipes")
    axis = _as_mapping(axis, "capsule.axis_source")
    uz_definition = _lookup(
        axis,
        ("Uz0", "U0", "source_Uz0", "axis_Uz0_definition"),
    )
    h_definition = _lookup(axis, ("H", "source_H", "H_definition"))
    if (
        uz_definition is None
        or "4" not in str(uz_definition)
        or "j" not in str(uz_definition)
    ):
        raise KernelInputError("exact Uz0=4Z+j axis source binding is missing")
    h_text = str(h_definition or "")
    expanded_h = "-4" in h_text
    factored_h = "1-Z" in h_text and ("U0" in h_text or "Uz0" in h_text)
    if h_definition is None or not (expanded_h or factored_h):
        raise KernelInputError("exact cubic H axis source binding is missing")
    return {
        "source_bindings": _json_value(axis),
        "Uz0_definition": str(uz_definition),
        "H_definition": str(h_definition),
        "source_axis_port_admitted": True,
    }


def _pressure_record(
    context: Any,
    pressure: Mapping[str, Any],
    names: Sequence[str],
    label: str,
) -> tuple[Any, Any]:
    raw = _lookup(pressure, names)
    if raw is None:
        raise KernelInputError("missing pressure source record " + label)
    return raw, _interval_value(context, raw, "pressure_source." + label)


def _pressure_slot(
    context: Any,
    pressure: Mapping[str, Any],
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the actual normalized P0/P0Z/P0ZZ function enclosures.

    The capsule supplies P0(0), P0,ZZ(0), the physical log-scale record and
    the stage/cap DATA.  The beta support is [0,2].  No missing pressure row
    is replaced by zero and no bound is promoted to a point value.
    """

    selectors = _pressure_selectors(pressure)
    if selectors is None or selectors == {} or selectors == "":
        raise KernelInputError("pressure source selectors are required")
    source_definition = pressure.get("source_definition")
    if source_definition is None:
        source_definition = pressure.get("proof_pressure_source")
    if source_definition is None or source_definition == "":
        raise KernelInputError("pressure source definition is required")

    beta_support = pressure.get("beta_support")
    if _is_mapping(beta_support):
        beta_low = beta_support.get("lower")
        beta_high = beta_support.get("upper")
    elif _is_sequence(beta_support) and len(beta_support) == 2:
        beta_low, beta_high = beta_support
    else:
        raise KernelInputError("pressure beta_support [0,2] is required")
    if str(beta_low) not in {"0", "0.0"} or str(beta_high) not in {"2", "2.0"}:
        raise KernelInputError("pressure beta_support must be exactly [0,2]")

    stage_value = pressure.get("source_stage_metadata")
    if stage_value is None:
        stage_value = pressure.get("stage_counts")
    if stage_value is None:
        stage_value = pressure.get("stage_partition")
    stage = _as_mapping(stage_value, "pressure_source.source_stage_metadata")
    stage_partition_value = pressure.get("stage_partition")
    if stage_partition_value is None:
        raise KernelInputError("pressure stage_partition DATA is required")
    stage_partition = (
        _as_mapping(stage_partition_value, "pressure_source.stage_partition")
        if stage_partition_value is not None
        else None
    )
    if pressure.get("all_prior_errors_and_joint_late_bound_retained") is not True:
        raise KernelInputError(
            "pressure all_prior_errors_and_joint_late_bound_retained DATA is required"
        )
    same_family_admission = pressure.get("same_family_admission")
    if same_family_admission is None:
        raise KernelInputError("pressure same_family_admission DATA is required")
    parameter_data_consumed = pressure.get("parameter_data_consumed")
    if parameter_data_consumed is None:
        raise KernelInputError("pressure parameter_data_consumed DATA is required")
    uniform_cell = pressure.get("uniform_cell")
    if uniform_cell is None:
        raise KernelInputError("pressure uniform_cell DATA is required")

    def stage_count(names: Sequence[str], label: str) -> int:
        value = _lookup(stage, names)
        if value is None:
            raise KernelInputError("missing pressure stage count " + label)
        try:
            return int(value)
        except (TypeError, ValueError) as exc:
            raise KernelInputError("invalid pressure stage count " + label) from exc

    stage_total = stage_count(
        ("total", "stage_count", "count", "source_stage_count"), "total"
    )
    stage_early = stage_count(
        ("early", "early_count", "early_stage_count"), "early"
    )
    stage_late = stage_count(
        ("late", "late_count", "late_stage_count"), "late"
    )
    if (stage_total, stage_early, stage_late) != (14, 3, 11):
        raise KernelInputError("pressure stage metadata must be 14/3/11")
    late_cap_names = (
        "old_late_caps",
        "old_latecaps",
        "late_caps",
        "old_late_cap_data",
        "late_caps_data",
    )
    late_caps = _lookup(pressure, late_cap_names)
    if late_caps is None:
        late_caps = _lookup(stage, late_cap_names)
    if late_caps is None and stage_partition is not None:
        late_caps = _lookup(stage_partition, late_cap_names)
    if late_caps is None and stage_partition is not None:
        late_caps = stage_partition.get("late_stage_definitions")
    if late_caps is None:
        late_caps = stage.get("late_stage_definitions")
    if late_caps is None:
        raise KernelInputError("old late-cap DATA must be retained")

    zero_row = pressure.get("zero_row")
    anchor_row = pressure.get("anchor_row")
    if zero_row is None or anchor_row is None:
        raise KernelInputError("pressure zero_row and anchor_row DATA are required")
    if zero_row is not None:
        zero_row = _as_mapping(zero_row, "pressure_source.zero_row")
    if anchor_row is not None:
        anchor_row = _as_mapping(anchor_row, "pressure_source.anchor_row")

    def row_lookup(row: Mapping[str, Any] | None, names: Sequence[str]) -> Any:
        return _lookup(row, names) if row is not None else None

    p0_raw = row_lookup(zero_row, ("pressure", "P0_norm", "P0_norm_Z0"))
    if p0_raw is None:
        p0_raw = _lookup(pressure, ("P0_norm_Z0", "P0_norm_at_Z0"))
    if p0_raw is None:
        raise KernelInputError("missing pressure zero-row P0_norm record")
    p0_zero = _interval_value(context, p0_raw, "pressure_source.P0_norm_Z0")

    dzz_raw = row_lookup(
        zero_row,
        ("D_ZZ_pressure", "P0_norm_DZZ", "M1_norm_Z0_DZZ"),
    )
    if dzz_raw is None:
        dzz_raw = _lookup(
            pressure,
            ("M1_norm_Z0_DZZ", "P0_norm_DZZ_Z0", "P0_norm_Z0_DZZ"),
        )
    if dzz_raw is None and anchor_row is not None:
        dzz_raw = row_lookup(
            anchor_row, ("D_ZZ_pressure", "P0_norm_DZZ", "M1_norm_Z0_DZZ")
        )
    if dzz_raw is None:
        raise KernelInputError("missing pressure zero-row D_ZZ record")
    p0_dzz_zero = _interval_value(
        context, dzz_raw, "pressure_source.M1_norm_Z0_DZZ"
    )

    anchor_p_raw = row_lookup(anchor_row, ("pressure", "P0_norm", "P0"))
    anchor_dz_raw = row_lookup(
        anchor_row, ("D_Z_pressure", "P0_norm_Z", "P0_Z")
    )
    anchor_dzz_raw = row_lookup(
        anchor_row, ("D_ZZ_pressure", "P0_norm_ZZ", "P0_ZZ")
    )
    if anchor_p_raw is None or anchor_dz_raw is None or anchor_dzz_raw is None:
        raise KernelInputError(
            "anchor pressure row must contain pressure, D_Z_pressure and D_ZZ_pressure"
        )
    anchor_p = _interval_value(context, anchor_p_raw, "pressure_source.anchor.pressure")
    anchor_dz = _interval_value(context, anchor_dz_raw, "pressure_source.anchor.D_Z")
    anchor_dzz = _interval_value(
        context, anchor_dzz_raw, "pressure_source.anchor.D_ZZ"
    )
    logscale_raw, logscale = _pressure_record(
        context,
        pressure,
        ("physical_log_scale_Pstar_squared", "physical_log_scale"),
        "physical_log_scale_Pstar_squared",
    )
    epsilon_raw, epsilon_source = _pressure_record(
        context,
        pressure,
        ("epsilon_core", "epsilon"),
        "epsilon_core",
    )
    core_seed_scale_raw, core_seed_scale = _pressure_record(
        context,
        pressure,
        ("core_seed_scale", "epsilon_Pstar_squared"),
        "core_seed_scale",
    )
    core_seed_log_scale_raw, core_seed_log_scale = _pressure_record(
        context,
        pressure,
        ("core_seed_log_scale", "log_core_seed_scale"),
        "core_seed_log_scale",
    )
    if not _strict_positive(context, epsilon_source):
        raise KernelInputError("pressure epsilon_core must be positive")
    if not _strict_positive(context, core_seed_scale):
        raise KernelInputError("pressure core_seed_scale must be positive")

    B0 = -p0_zero
    B1 = p0_dzz_zero / 2
    if not _strict_positive(context, B0):
        raise KernelInputError("pressure B0=-P0_norm(0) must be positive")
    if not _strict_positive(context, B1):
        raise KernelInputError("pressure B1=P0_norm_DZZ(0)/2 must be positive")

    rows: dict[str, Any] = {}
    for eta_text in ETA_TEXTS:
        eta = _exact_eta(context, eta_text)
        anchor_reused = eta_text == "0"
        if anchor_reused:
            # eta=0 is the saved anchor DATA row, not a newly evaluated
            # beta-kernel envelope.  Keep its complete ordinary pressure jets.
            z = anchor["interval"]
            actual_p0 = anchor_p
            p0_z = anchor_dz
            p0_zz = anchor_dzz
            I1 = None
            I2 = None
            positive_kernel_evaluated = False
        else:
            source = _source_polynomial(context, packet, anchor, eta)
            z = source["Z"]
            x = z**2
            q = 1 + x
            q_lower = _lower(context, q)
            q_upper = _upper(context, q)
            if not bool(q_lower > context.mpf(0)):
                raise KernelInputError("pressure q=1+Z^2 is not positive")
            p0_lower = -_upper(context, B0)
            p0_upper = -_lower(context, B0) / q_upper**2
            actual_p0 = context.mpf([p0_lower, p0_upper])
            I1 = context.mpf(
                [
                    _lower(context, B1) / q_upper**3,
                    _upper(context, B1) / q_lower,
                ]
            )
            I2 = context.mpf(
                [
                    _lower(context, B1) / q_upper**4,
                    3 * _upper(context, B1) / q_lower**2,
                ]
            )
            p0_z = 2 * z * I1
            p0_zz = 2 * I1 - 4 * x * I2
            positive_kernel_evaluated = True
        x = z**2
        q = 1 + x
        scaled_p0 = core_seed_scale * actual_p0
        scaled_p0_z = core_seed_scale * p0_z
        scaled_p0_zz = core_seed_scale * p0_zz
        rows[eta_text] = {
            "eta": _interval_record(context, eta),
            "Z": _interval_record(context, z),
            "x=Z^2": _interval_record(context, x),
            "q=1+x": _interval_record(context, q),
            "P0_norm": _interval_record(context, actual_p0),
            "I1": _interval_record(context, I1) if I1 is not None else None,
            "P0_norm_Z": _interval_record(context, p0_z),
            "I2": _interval_record(context, I2) if I2 is not None else None,
            "P0_norm_ZZ": _interval_record(context, p0_zz),
            "definitions": (
                {
                    "P0_norm": "saved anchor_row.pressure DATA at eta=0",
                    "P0_norm_Z": "saved anchor_row.D_Z_pressure DATA at eta=0",
                    "P0_norm_ZZ": "saved anchor_row.D_ZZ_pressure DATA at eta=0",
                    "I1": "not newly evaluated at eta=0",
                    "I2": "not newly evaluated at eta=0",
                }
                if anchor_reused
                else {
                    "P0_norm": "-integral m*q^(-beta)",
                    "I1": "integral beta*m*q^(-(beta+1))",
                    "P0_norm_Z": "2*Z*I1",
                    "I2": "integral beta*(beta+1)*m*q^(-(beta+2))",
                    "P0_norm_ZZ": "2*I1-4*x*I2",
                }
            ),
            "positive_kernel_envelope_evaluated": positive_kernel_evaluated,
            "actual_function_enclosure": True,
            "point_value": False,
            "anchor_row_reused": anchor_reused,
        }
        rows[eta_text]["scaled_core_pressure"] = _interval_record(
            context, scaled_p0
        )
        rows[eta_text]["scaled_core_pressure_Z"] = _interval_record(
            context, scaled_p0_z
        )
        rows[eta_text]["scaled_core_pressure_ZZ"] = _interval_record(
            context, scaled_p0_zz
        )
    return {
        "status": "actual_normalized_pressure_enclosure",
        "selectors": _json_value(selectors),
        "source_definition": _json_value(source_definition),
        "beta_support": [0, 2],
        "source_stage_metadata": {
            "total": stage_total,
            "early": stage_early,
            "late": stage_late,
        },
        "stage_partition": _json_value(stage_partition),
        "old_late_caps": _json_value(late_caps),
        "old_late_caps_retained": True,
        "all_prior_errors_and_joint_late_bound_retained": True,
        "same_family_admission": _json_value(same_family_admission),
        "parameter_data_consumed": _json_value(parameter_data_consumed),
        "uniform_cell": _json_value(uniform_cell),
        "zero_row": _json_value(zero_row),
        "anchor_row_raw": _json_value(anchor_row),
        "P0_norm_Z0": {
            "raw": _json_value(p0_raw),
            "interval": _interval_record(context, p0_zero),
        },
        "M1_norm_Z0_DZZ": {
            "raw": _json_value(dzz_raw),
            "interval": _interval_record(context, p0_dzz_zero),
        },
        "B0": _interval_record(context, B0),
        "B1": _interval_record(context, B1),
        "anchor_row": {
            "pressure": _interval_record(context, anchor_p),
            "D_Z_pressure": _interval_record(context, anchor_dz),
            "D_ZZ_pressure": _interval_record(context, anchor_dzz),
            "direct_data_reused_at_eta_zero": True,
        },
        "core_seed_scale": {
            "raw": _json_value(core_seed_scale_raw),
            "interval": _interval_record(context, core_seed_scale),
            "definition": "epsilon_core*Pstar^2 saved DATA",
            "saved_DATA_required": True,
            "formal_only_if_missing": False,
        },
        "core_seed_log_scale": {
            "raw": _json_value(core_seed_log_scale_raw),
            "interval": _interval_record(context, core_seed_log_scale),
            "formal_only": True,
        },
        "physical_log_scale_Pstar_squared": {
            "raw": _json_value(logscale_raw),
            "interval": _interval_record(context, logscale),
            "formal_only": True,
            "old_exp_parameter_family_reentered": False,
        },
        "epsilon_core": {
            "raw": _json_value(epsilon_raw),
            "interval": _interval_record(context, epsilon_source),
        },
        "rows": rows,
        "axis_envelope_present": True,
        "axis_envelope_value_admitted": True,
        "accepted": True,
        "missing_is_not_zero": True,
        "point_value": False,
    }


def _peak_width(
    context: Any,
    packet: Mapping[str, Any],
    anchor: Mapping[str, Any],
) -> dict[str, Any]:
    minus_one = _exact_eta(context, "-1")
    plus_one = _exact_eta(context, "1")
    all_range = context.mpf([-1, 1])
    all_b = _B_range(context, packet, anchor, all_range)
    if not _strict_positive(context, all_b["B"]):
        raise KernelInputError("B_all is not strictly positive")
    D_minus, minus_b = _D_for_eta(context, packet, anchor, minus_one)
    D_plus, plus_b = _D_for_eta(context, packet, anchor, plus_one)
    log_two = context.log(2)
    if not bool(_lower(context, D_minus) > log_two):
        raise KernelInputError("D(-1)>log(2) source certificate failed")
    if not bool(_lower(context, D_plus) > log_two):
        raise KernelInputError("D(+1)>log(2) source certificate failed")
    eta_lower = context.sqrt(2 * log_two / _upper(context, all_b["B"]))
    eta_upper = context.sqrt(2 * log_two / _lower(context, all_b["B"]))
    if not bool(_upper(context, eta_upper) <= context.mpf(1)):
        raise KernelInputError("half-height width does not fit the certified eta cell")
    offset = packet["scale_s"]["interval"] * context.mpf(
        [_lower(context, eta_lower), _upper(context, eta_upper)]
    )
    return {
        "B_all_on_signed_eta_segment": _interval_record(context, all_b["B"]),
        "D_minus_one": _interval_record(context, D_minus),
        "D_plus_one": _interval_record(context, D_plus),
        "log_two": _interval_record(context, log_two),
        "eta_halfheight_lower": _interval_record(context, eta_lower),
        "eta_halfheight_upper": _interval_record(context, eta_upper),
        "offset_width": _interval_record(context, offset),
        "center": "exact implicit anchor a; no midpoint selected",
        "comparison_theorem": (
            "B_all>0 makes D monotone away from eta=0; "
            "sqrt(2 log(2)/B_all.upper) <= |eta_half| <= "
            "sqrt(2 log(2)/B_all.lower)."
        ),
        "width_kind": "F0-factor half-height only",
        "not_full_profile_width": True,
        "not_material_trajectory": True,
        "not_physical_swirl_maximum": True,
        "signed_segment_records": {
            "minus": minus_b["positive_weight_theorem"],
            "plus": plus_b["positive_weight_theorem"],
        },
    }


def derive(
    data: Mapping[str, Any],
    emit: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Evaluate only the new exact-anchor source-axis leaf.

    The root supplies {"capsule": frozen_capsule}. The function is not
    called by this worker turn; the root owns freeze, one invocation, journaling
    and result persistence.
    """

    if not _is_mapping(data) or "capsule" not in data:
        raise KernelInputError("derive requires data['capsule']")
    if emit is not None and not callable(emit):
        raise KernelInputError("emit must be callable or None")

    capsule, parameters, anchor_raw, pressure = _validate_capsule(data["capsule"])
    events: list[dict[str, Any]] = []

    _emit(
        emit,
        events,
        "family_entered",
        {
            "capsule_schema": CAPSULE_SCHEMA,
            "precision_dps": INTERVAL_DPS,
            "precision_bits": SOURCE_PRECISION_BITS,
            "source_ids": _json_value(capsule["source_ids"]),
            "exact_root_correlated": True,
            "old_math_or_provider_calls": 0,
        },
    )

    _emit(emit, events, "anchor_polynomial_entered", {"family": "anchor_polynomial"})
    context = _new_interval_context()
    packet = _parameter_packet(context, parameters)
    if packet["Lambda_epsilon_identity"] is None:
        packet["Lambda_epsilon_identity"] = _lookup(
            capsule,
            ("Lambda_epsilon_exact_identity", "Lambda_epsilon_equals_one"),
        )
    anchor = _anchor_packet(context, anchor_raw)
    axis_binding = _axis_binding(capsule)
    anchor_rows: dict[str, Any] = {}
    for eta_text in ETA_TEXTS:
        eta = _exact_eta(context, eta_text)
        raw_row = _source_polynomial(context, packet, anchor, eta)
        anchor_rows[eta_text] = {
            key: _interval_record(context, value)
            if key in {"eta", "h", "Z", "H1", "L0", "K", "Lh"}
            else value
            for key, value in raw_row.items()
        }
    _emit(
        emit,
        events,
        "anchor_polynomial_complete",
        {
            "family": "anchor_polynomial",
            "eta_count": len(ETA_TEXTS),
            "exact_H_anchor_identity": True,
        },
    )

    _emit(
        emit,
        events,
        "scaled_gauge_integral_entered",
        {"family": "scaled_gauge_integral"},
    )
    integral_rows: dict[str, Any] = {}
    for eta_text in ETA_TEXTS:
        eta = _exact_eta(context, eta_text)
        D, bpacket = _D_for_eta(context, packet, anchor, eta)
        integral_rows[eta_text] = {
            "eta": _interval_record(context, eta),
            "D": _interval_record(context, D),
            "B_signed_segment": _interval_record(context, bpacket["B"]),
            "segment": _interval_record(context, bpacket["xi_range"]),
            "definition": (
                "D(eta)=integral_0^eta xi*B(xi)d xi; "
                "B is evaluated over the whole signed segment"
            ),
            "positive_weight_transform": bpacket["positive_weight_theorem"],
            "G_anchor_zero": True,
            "source_root_correlation": "H(a)=0 exact formal identity",
            "Gaussian_replacement": False,
        }
    _emit(
        emit,
        events,
        "scaled_gauge_integral_complete",
        {
            "family": "scaled_gauge_integral",
            "eta_count": len(ETA_TEXTS),
            "old_G0_transport_calls": 0,
            "old_quadrature_calls": 0,
        },
    )

    _emit(
        emit,
        events,
        "relative_profile_jets_entered",
        {"family": "relative_profile_jets"},
    )
    jet_rows: dict[str, Any] = {}
    axis_rows: dict[str, Any] = {}
    for eta_text in ETA_TEXTS:
        eta = _exact_eta(context, eta_text)
        D_interval = _interval_value(
            context, integral_rows[eta_text]["D"], "integral_rows.D"
        )
        jet_rows[eta_text] = _relative_jet(
            context, packet, anchor, eta, D_interval
        )
        axis_rows[eta_text] = _axis_source_rows(
            context, packet, anchor, eta
        )
        jet_rows[eta_text]["relative_profile"] = _safe_relative_profile(
            context, D_interval
        )
    _emit(
        emit,
        events,
        "relative_profile_jets_complete",
        {
            "family": "relative_profile_jets",
            "ordinary_Z1_definition": "D_Z(logF0)=-D_eta/s",
            "absolute_F0_exp_materialized": False,
        },
    )

    _emit(
        emit,
        events,
        "axis_pressure_envelope_entered",
        {"family": "axis_pressure_envelope"},
    )
    pressure_slot = _pressure_slot(context, pressure, packet, anchor)
    if pressure_slot["axis_envelope_value_admitted"] is not True:
        raise KernelInputError("axis pressure envelope DATA was not admitted")
    actual_flags = dict(SUCCESS_FLAGS)
    pressure_row = {
        "family": "axis_pressure_envelope",
        "pressure_source": pressure_slot,
        "definition": "P=P0+integral F_c^2 dR; P_scaled=epsilon_core*P",
        "missing_pressure_is_not_zero": True,
        "accepted_as_nonlinear_off_axis_core": False,
    }
    _emit(
        emit,
        events,
        "axis_pressure_envelope_complete",
        {
            "family": "axis_pressure_envelope",
            "selector_bound": True,
            "value_admitted": pressure_slot["axis_envelope_value_admitted"],
        },
    )

    _emit(
        emit,
        events,
        "peak_width_comparison_entered",
        {"family": "peak_width_comparison"},
    )
    width = _peak_width(context, packet, anchor)
    _emit(
        emit,
        events,
        "peak_width_comparison_complete",
        {
            "family": "peak_width_comparison",
            "width_kind": "F0-factor half-height",
            "center_exact_anchor": True,
        },
    )

    proof = {
        "kind": "exact_root_scaled_source_interval_v1",
        "proof_count": 1,
        "root_identity": "H(a)=0 is retained symbolically; no root midpoint",
        "integral_identity": (
            "Lambda*(G(a+s*eta)-G(a))="
            "integral_0^eta xi*L(a+s*xi)*K(s*xi)/"
            "(1+epsilon_core*xi^2*K(s*xi)^2) dxi"
        ),
        "positive_weight_change_of_variable": True,
        "B_all_strictly_positive": True,
        "D_plus_minus_one_strictly_above_log_two": True,
        "absolute_F0_exp_materialized": False,
        "nonlinear_off_axis_core_values_admitted": False,
        "pressure_missing_not_defaulted": True,
    }
    _emit(emit, events, "proof_bundle", proof)
    _emit(
        emit,
        events,
        "family_completed",
        {
            "proof_count": 1,
            "source_only": True,
            "numeric_evaluation": True,
            "full_goal": False,
        },
    )

    payload = {
        "component": COMPONENT,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "source_capsule_schema": CAPSULE_SCHEMA,
        "runtime_model_metadata": None,
        "api": dict(API),
        "api_metadata": dict(API_METADATA),
        "canonical_gate": {
            "algorithm": "sha256(canonical_json(entire_capsule))",
            "actual_sha256": hashlib.sha256(_canonical_bytes(capsule)).hexdigest(),
            "root_constant_used": ROOT_CAPSULE_CANONICAL,
            "canonical_field_in_capsule": False,
        },
        "source_ids": _json_value(capsule["source_ids"]),
        "parameter_data": {
            key: {
                "raw": value["raw"],
                "interval": _interval_record(context, value["interval"])
                if _is_mapping(value) and "interval" in value
                else value,
            }
            for key, value in packet.items()
            if key not in {"Lambda_epsilon_identity", "parameter_point_selection"}
        },
        "parameter_identities": {
            "Lambda_epsilon_exact_identity": packet["Lambda_epsilon_identity"],
            "sigma": "saved sigma DATA; same-source exact recipe sigma=j/500",
            "scale": "sigma*sqrt(epsilon_core)=sigma/sqrt(Lambda)",
            "parameter_point_selected": False,
        },
        "anchor": {
            "record": anchor["record"],
            "exact_definition": anchor["exact_definition"],
            "root_tag": anchor["root_tag"],
            "midpoint_selected": False,
            "H_at_anchor_exact_identity": True,
        },
        "coordinate": {
            "eta_values": list(ETA_TEXTS),
            "Z": "exact formal a+s*eta",
            "s": "sigma*sqrt(epsilon_core)",
            "ordinary_Z_is_not_eta": True,
            "root_interval_midpoint_used": False,
        },
        "anchor_polynomial": {
            "rows": anchor_rows,
            "definitions": {
                "Uz0": "Uz0=4Z+j (affine source axis)",
                "H_factorized": "H=(1-Z^2)*Uz0+(1-delta)*Z/2",
                "H1": "(9-delta)/2-12*a^2-2*j*a",
                "K": "H1-(12*a+j)*h-4*h^2",
                "Lh": "L0-2*delta*a*h-delta*h^2",
                "h": "s*eta",
            },
            "exact_H_anchor_identity": True,
        },
        "scaled_gauge_integral": {
            "rows": integral_rows,
            "D_definition": "D(eta)=Lambda*G(a+s*eta), G(a)=0",
            "B_definition": (
                "B(xi)=L(a+s*xi)*K(s*xi)/"
                "(1+epsilon_core*xi^2*K(s*xi)^2)"
            ),
            "signed_segment_directed_enclosure": True,
            "positive_weight_change_of_variable": "xi=eta*t",
            "old_G0_transport_called": False,
            "old_quad_called": False,
        },
        "relative_profile_jets": {
            "rows": jet_rows,
            "axis_rows": axis_rows,
            "relative_profile_definition": "F0(a+s*eta)/F0(a)=exp(-D(eta))",
            "absolute_logF0_definition": "logF0=-logC-D",
            "absolute_F0_definition": "F0=exp(-logC-D)",
            "absolute_F0_exp_materialized": False,
            "ordinary_Z1_not_coordinate_Z1": True,
        },
        "axis_source": axis_binding,
        "axis_pressure_envelope": pressure_row,
        "peak_width_comparison": width,
        "physical_width_slot": {
            "status": "definition_only_not_evaluated",
            "tau_values": ["1", "1/4", "1/16", "1/64"],
            "map": "lambda=sqrt(tau/(1-Z^2)); kappa=lambda^(1-delta); z=kappa*Z",
            "dz_dZ": "kappa*L/(1-Z^2)",
            "center": "exact root a",
            "width_is_not_material_trajectory": True,
        },
        "proof_packet": proof,
        "events": events,
        "actual_counts_contract": dict(SUCCESS_COUNTS),
        "actual_flags_contract": dict(actual_flags),
        "expected_counts": dict(EXPECTED_COUNTS),
        "expected_flags": dict(EXPECTED_FLAGS),
        "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
        "callback_event_count": len(CALLBACK_EVENTS),
        "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "count_event_map": {
            key: list(value) for key, value in COUNT_EVENT_MAP.items()
        },
        "preservation": {
            "all_old_caps_and_errors_preserved": True,
            "nonlinear_off_axis_core_values_admitted": False,
            "global_five_moment": False,
            "time_n_recursion": False,
        },
    }
    result_out = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "proof_count": 1,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "runtime_model_metadata": None,
        "api": dict(API),
        "api_metadata": dict(API_METADATA),
        "actual_counts": dict(SUCCESS_COUNTS),
        "new_flags": dict(actual_flags),
        "expected_counts": dict(EXPECTED_COUNTS),
        "expected_flags": dict(EXPECTED_FLAGS),
        "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
        "callback_event_count": len(CALLBACK_EVENTS),
        "manager_journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "journal_event_count": MANAGER_JOURNAL_EVENT_COUNT,
        "count_event_map": {
            key: list(value) for key, value in COUNT_EVENT_MAP.items()
        },
        "physical_layer_completed": False,
        "full_goal": False,
        "full_goal_complete": False,
        "payload": payload,
    }
    json.dumps(
        _json_value(result_out),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return result_out


__all__ = [
    "API",
    "API_METADATA",
    "ANCHOR_SOURCE_SHA256",
    "CALLBACK_EVENT_COUNTS",
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "EXPECTED_COUNTS",
    "EXPECTED_FLAGS",
    "INTERVAL_DPS",
    "KernelInputError",
    "MANAGER_JOURNAL_EVENT_COUNT",
    "PRESSURE_SOURCE_SHA256",
    "ROOT_CAPSULE_CANONICAL",
    "SCHEMA",
    "SOURCE_IDS",
    "SOURCE_PRECISION_BITS",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]
