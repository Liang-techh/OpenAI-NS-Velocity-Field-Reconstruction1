"""Saved-core physical geometry binding leaf.

This module is intentionally a fresh consumer of a frozen ``core_payload``.
It does not import a repository provider, recurrence, checker, or previous
kernel.  The root runner supplies the saved result data and decides when this
leaf is frozen and invoked.  All finite core quantities are directed interval
enclosures; the enormous positive amplitude ``F0=exp(logF0_Z0)`` remains a
formal atom and is never exponentiated or materialized here.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_core_physical_source_binding_result_v1"
CAPSULE_SCHEMA = "native_O3_core_physical_source_binding_source_capsule_v1"
COMPONENT = "native_O3_core_physical_source_binding"
SOURCE_RESULT_SHA256 = (
    "d5553374bc5ce9919e50892c02e09ada1e8d1c18ea34bb2894dbcf722424618f"
)
CANONICAL_CAPSULE_SHA256 = (
    "012511ca881a0207e04b7f494ce8ca275d7ea8a0a04d2914a824b5f846d64e0e"
)
SOURCE_CAPSULE_CANONICAL_SHA256 = (
    "c4f0dd60a7856a489ad1ddbc3f119da33d59036001eb8ce53ac2d8830eb4788e"
)
INTERVAL_DPS = 240
RHO_POINTS = ("1/16", "1/4", "1", "3")
TAU_POINTS = ("1", "1/4", "1/16", "1/64")
JET_KEYS = (
    "radial_0_axial_0",
    "radial_1_axial_0",
    "radial_0_axial_1",
    "radial_0_axial_2",
    "radial_1_axial_1",
)
THETA_CASES = {
    "0": {"cos": "1", "sin": "0"},
    "pi/2": {"cos": "0", "sin": "1"},
}

SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 4,
    "physical_geometry_family_calls": 1,
    "physical_values_family_calls": 1,
    "physical_jets_family_calls": 1,
    "common_source_growth_family_calls": 1,
    "derive_calls": 1,
    "numeric_evaluation_calls": 1,
    "physical_core_packets_evaluated_calls": 1,
    "new_mathematical_proofs": 1,
}

SUCCESS_FLAGS = {
    "physical_core_packets_evaluated": True,
    "physical_geometry_evaluated": True,
    "physical_values_evaluated": True,
    "physical_jets_evaluated": True,
    "common_source_growth_evaluated": True,
    "source_only": True,
    "metadata_only": False,
    "numeric_evaluation": True,
    "full_physical_field_materialized": False,
    "physical_field_materialized": False,
    "full_physical_NS_evaluation": False,
    "full_goal": False,
    "full_goal_complete": False,
    "N_materialized": False,
    "F0_value_materialized": False,
    "F0_exp_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "checker_called": False,
    "selftest_called": False,
    "trajectory_computed": False,
    "cumulative_winding_computed": False,
    "formal_streamline_pitch_only": True,
    "old_caps_and_errors_preserved": True,
    "global_errors_preserved": True,
}

CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "physical_geometry_entered": 1,
    "physical_geometry_complete": 1,
    "physical_values_entered": 1,
    "physical_values_complete": 1,
    "physical_jets_entered": 1,
    "physical_jets_complete": 1,
    "common_source_growth_entered": 1,
    "common_source_growth_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}

COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": [
        "physical_geometry_entered",
        "physical_values_entered",
        "physical_jets_entered",
        "common_source_growth_entered",
    ],
    "physical_geometry_family_calls": ["physical_geometry_entered"],
    "physical_values_family_calls": ["physical_values_entered"],
    "physical_jets_family_calls": ["physical_jets_entered"],
    "common_source_growth_family_calls": ["common_source_growth_entered"],
    "derive_calls": ["family_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "physical_core_packets_evaluated_calls": ["family_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}

EXPECTED_COUNTS = SUCCESS_COUNTS
EXPECTED_FLAGS = SUCCESS_FLAGS
API = {
    "name": "derive",
    "signature": "derive(data, emit=None)",
    "capsule_selector": "data['capsule']",
    "result_schema": SCHEMA,
    "source_only": True,
    "numeric": True,
    "metadata_only": False,
}


class BindingInputError(ValueError):
    """Raised when a frozen saved-core binding capsule is malformed."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _json_value(value: Any) -> Any:
    """Copy only JSON-shaped source metadata; never evaluate opaque data."""

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
        raise BindingInputError("event outside frozen callback contract " + name)
    row: dict[str, Any] = {
        "event": name,
        "component": COMPONENT,
        "source_only": True,
        "numeric_evaluation": name.endswith("complete") or name in {
            "family_completed",
        },
    }
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def _new_interval_context() -> Any:
    # Local standard-library-adjacent interval context only.  No project
    # module or previous numerical family is imported.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _endpoint_text(context: Any, value: Any) -> str:
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
    return {
        "lower": _endpoint_text(context, lower),
        "upper": _endpoint_text(context, upper),
        "width": _endpoint_text(context, upper - lower),
        "lower_exact_mpf_tuple": list(lower_raw),
        "upper_exact_mpf_tuple": list(upper_raw),
    }


def _interval_from_record(context: Any, value: Any, label: str) -> Any:
    if _is_mapping(value):
        low_tuple = value.get("lower_exact_mpf_tuple")
        high_tuple = value.get("upper_exact_mpf_tuple")
        if _is_sequence(low_tuple) and _is_sequence(high_tuple):
            if len(low_tuple) != 4 or len(high_tuple) != 4:
                raise BindingInputError(label + " MPF tuple must have length four")
            try:
                low = tuple(int(item) for item in low_tuple)
                high = tuple(int(item) for item in high_tuple)
                return context.make_mpf((low, high))
            except Exception as error:
                raise BindingInputError(label + " MPF tuple is malformed") from error
        if "lower" in value and "upper" in value:
            try:
                return context.mpf((str(value["lower"]), str(value["upper"])))
            except Exception as error:
                raise BindingInputError(label + " endpoints are malformed") from error
        for key in ("interval", "enclosure", "record", "value", "bound"):
            if key in value:
                return _interval_from_record(context, value[key], label)
        raise BindingInputError(label + " has no interval endpoints")
    if _is_sequence(value):
        if len(value) != 2:
            raise BindingInputError(label + " endpoint sequence must have length two")
        try:
            return context.mpf((str(value[0]), str(value[1])))
        except Exception as error:
            raise BindingInputError(label + " endpoint sequence is malformed") from error
    if value is None:
        raise BindingInputError(label + " is missing")
    try:
        return context.mpf(str(value))
    except Exception as error:
        raise BindingInputError(label + " is not interval-shaped") from error


def _exact_fraction(context: Any, text: str) -> Any:
    if "/" not in text:
        return context.mpf(text)
    numerator, denominator = text.split("/", 1)
    return context.mpf(numerator) / context.mpf(denominator)


def _lookup(mapping: Mapping[str, Any], names: Sequence[str], label: str) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    raise BindingInputError(label + " is missing")


def _path(root: Mapping[str, Any], keys: Sequence[str], label: str) -> Any:
    current: Any = root
    for key in keys:
        if not _is_mapping(current) or key not in current:
            raise BindingInputError(label + " is missing at " + key)
        current = current[key]
    return current


def _payload_from_data(data: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    if "capsule" not in data:
        raise BindingInputError("data['capsule'] is required")
    capsule: Any = data["capsule"]
    if not _is_mapping(capsule):
        raise BindingInputError("data.capsule must be a mapping")
    if "core_payload" not in capsule:
        raise BindingInputError("capsule.core_payload is required")
    payload = capsule["core_payload"]
    if not _is_mapping(payload):
        raise BindingInputError("capsule.core_payload must be a mapping")
    return capsule, payload


def _validate_source(capsule: Mapping[str, Any], payload: Mapping[str, Any]) -> dict[str, Any]:
    canonical_bytes = json.dumps(
        capsule,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    computed_capsule_hash = hashlib.sha256(canonical_bytes).hexdigest()
    if computed_capsule_hash != SOURCE_CAPSULE_CANONICAL_SHA256:
        raise BindingInputError("canonical source capsule content hash mismatch")
    source_result = capsule.get("source_result")
    if not _is_mapping(source_result):
        raise BindingInputError("capsule.source_result is required")
    if source_result.get("sha256") != SOURCE_RESULT_SHA256:
        raise BindingInputError("capsule.source_result.sha256 mismatch")
    if source_result.get("readonly_DATA") is not True:
        raise BindingInputError("saved source result is not marked readonly DATA")
    if source_result.get("ancestor_reexecution") is not False:
        raise BindingInputError("saved source result ancestor_reexecution must be false")
    if source_result.get("root_selector") != "result.payload":
        raise BindingInputError("saved source result root_selector mismatch")
    identity = payload.get("source_identity", {})
    if not _is_mapping(identity):
        raise BindingInputError("core source_identity is missing")
    hashes = identity.get("accepted_hashes", identity)
    if not _is_mapping(hashes):
        raise BindingInputError("accepted source hashes are missing")
    required = {
        "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
        "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
        "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
    }
    for key, value in required.items():
        if hashes.get(key) != value:
            raise BindingInputError("core source identity mismatch at " + key)
    supplied_pins = identity.get("supplied_artifact_hash_pins")
    if supplied_pins is not None and not _is_mapping(supplied_pins):
        raise BindingInputError("supplied artifact hash pins are malformed")
    frozen_authentication = identity.get("frozen_input_authentication")
    if not _is_mapping(frozen_authentication):
        raise BindingInputError("frozen input authentication is missing")
    core_capsule_hash = frozen_authentication.get("canonical_capsule_sha256")
    if core_capsule_hash != CANONICAL_CAPSULE_SHA256:
        raise BindingInputError("canonical capsule hash is missing or mismatched")
    return {
        "source_result_sha256": SOURCE_RESULT_SHA256,
        "canonical_capsule_sha256": core_capsule_hash,
        "source_capsule_canonical_sha256": computed_capsule_hash,
        "accepted_hashes": _json_value(hashes),
        "required_artifact_hash_pins": _json_value(
            identity.get("required_artifact_hash_pins")
        ),
    }


def _validate_capsule_guards(capsule: Mapping[str, Any]) -> dict[str, Any]:
    guards = capsule.get("guards", {})
    if not _is_mapping(guards):
        raise BindingInputError("capsule guards are malformed")
    required_true = (
        "all_212_terms_retained",
        "all_global_errors_retained",
        "finite_polynomial_is_not_complete_value",
        "reserved_contribution_once_retained",
        "six_signed_channels_retained",
        "tail_is_bound_not_value",
    )
    for key in required_true:
        if guards.get(key) is not True:
            raise BindingInputError("capsule guard is not admitted: " + key)
    for key in (
        "N_materialization_calls",
        "F0_exp_materialization_calls",
        "newton_or_root_calls",
        "old_family_calls",
    ):
        if key not in guards or guards[key] != 0:
            raise BindingInputError("capsule guard count is missing/nonzero: " + key)
    for key in ("full_goal_complete", "global_five_moment_closure", "physical_layer_completed"):
        if key not in guards or guards[key] is not False:
            raise BindingInputError("capsule guard is missing/unexpectedly complete: " + key)
    return _json_value(guards)


def _validate_grid(capsule: Mapping[str, Any]) -> dict[str, Any]:
    grid = capsule.get("grid", {})
    if not _is_mapping(grid):
        raise BindingInputError("capsule grid is malformed")
    if grid.get("Z") != "0":
        raise BindingInputError("capsule grid Z must be exact string '0'")
    for key, expected in (
        ("rho", list(RHO_POINTS)),
        ("tau", list(TAU_POINTS)),
        ("theta", ["0", "pi/2"]),
    ):
        value = grid.get(key)
        if not _is_sequence(value) or list(value) != expected:
            raise BindingInputError("capsule grid mismatch at " + key)
    return _json_value(grid)


def _core_cell(
    context: Any,
    payload: Mapping[str, Any],
    field: str,
    jet: str,
    rho: str,
) -> tuple[Any, Mapping[str, Any]]:
    jets = _path(payload, ("ordinary_jets", field), "ordinary_jets." + field)
    jet_map = _path(jets, (jet,), "ordinary_jets." + field + "." + jet)
    cell = _path(jet_map, (rho,), "ordinary_jets." + field + "." + jet + "." + rho)
    if not _is_mapping(cell):
        raise BindingInputError("saved ordinary jet cell is not a mapping")
    enclosure = cell.get("enclosure")
    if enclosure is None:
        raise BindingInputError("ordinary jet cell lacks complete enclosure")
    return _interval_from_record(context, enclosure, "ordinary jet enclosure"), cell


def _recovery_cell(
    context: Any,
    payload: Mapping[str, Any],
    name: str,
    rho: str,
) -> tuple[Any, Mapping[str, Any]]:
    root = _path(payload, ("local_core_recovery", name, "values"), name)
    cell = _path(root, (rho,), name + ".values." + rho)
    if not _is_mapping(cell) or "enclosure" not in cell:
        raise BindingInputError(name + " saved enclosure is missing")
    return _interval_from_record(context, cell["enclosure"], name + ".enclosure"), cell


def _parameter_interval(
    context: Any,
    payload: Mapping[str, Any],
    name: str,
) -> Any:
    majorants = payload.get("majorant_records", {})
    if _is_mapping(majorants) and name in majorants:
        return _interval_from_record(context, majorants[name], "majorant " + name)
    amplitude = _path(payload, ("core_amplitude_Z0", "parameter_intervals"), "parameter_intervals")
    return _interval_from_record(context, _lookup(amplitude, (name,), name), name)


def _sign_status(context: Any, value: Any) -> str:
    lower = context.mpf(_endpoint_text(context, _endpoint(context, value)))
    upper = context.mpf(_endpoint_text(context, _endpoint(context, value, True)))
    zero = context.mpf("0")
    if lower > zero:
        return "positive"
    if upper < zero:
        return "negative"
    if lower == zero and upper == zero:
        return "zero"
    return "uncertain"


def _formal_F0_product(
    context: Any,
    coefficient: Any,
    logF0: Any,
    label: str,
) -> dict[str, Any]:
    sign = _sign_status(context, coefficient)
    logabs = None
    if sign in {"positive", "negative"}:
        magnitude = coefficient if sign == "positive" else -coefficient
        logabs = logF0 + context.log(magnitude)
    return {
        "kind": "formal_F0_product",
        "label": label,
        "coefficient_interval": _interval_record(context, coefficient),
        "coefficient_sign": sign,
        "logF0_interval": _interval_record(context, logF0),
        "logabs_interval": _interval_record(context, logabs) if logabs is not None else None,
        "exact_F0_atom_ref": "result.payload.core_amplitude_Z0.amplitude.logF0_Z0",
        "exactF0atomref": "result.payload.core_amplitude_Z0.amplitude.logF0_Z0",
        "exact_definition": "F0=exp(logF0_Z0), product=F0*coefficient",
        "strictly_positive_F0_formal": True,
        "exp_materialized": False,
        "zero": False if sign in {"positive", "negative"} else (True if sign == "zero" else None),
        "F0_atom_zero": False,
        "coefficient_zero_possible": sign == "uncertain",
        "materialized_function_value": False,
        "exact_source_function_defined": True,
    }


def _source_records(context: Any, payload: Mapping[str, Any], rho: str) -> dict[str, Any]:
    values: dict[str, Any] = {}
    selectors: dict[str, str] = {}
    ordinary = payload.get("ordinary_jets", {})
    p_field = "P_scaled" if _is_mapping(ordinary) and _is_mapping(ordinary.get("P_scaled")) else "P"
    for field, jet in (
        ("Phi", "radial_0_axial_0"),
        ("Uz", "radial_0_axial_0"),
        ("Mz_over_R", "radial_0_axial_0"),
        (p_field, "radial_0_axial_0"),
    ):
        field_value, _ = _core_cell(context, payload, field, jet, rho)
        output_field = "P_scaled" if field == p_field and field == "P" else field
        values[output_field] = field_value
        selectors[output_field] = (
            "result.payload.ordinary_jets."
            + field
            + "."
            + jet
            + "."
            + rho
            + ".enclosure"
        )
    values["U_rho"], _ = _core_cell(
        context, payload, "Uz", "radial_1_axial_0", rho
    )
    selectors["U_rho"] = "result.payload.ordinary_jets.Uz.radial_1_axial_0." + rho + ".enclosure"
    values["U_Z"], _ = _core_cell(
        context, payload, "Uz", "radial_0_axial_1", rho
    )
    selectors["U_Z"] = "result.payload.ordinary_jets.Uz.radial_0_axial_1." + rho + ".enclosure"
    values["A_rho"], _ = _core_cell(
        context, payload, "Phi", "radial_1_axial_0", rho
    )
    selectors["A_rho"] = "result.payload.ordinary_jets.Phi.radial_1_axial_0." + rho + ".enclosure"
    values["A_Z"], _ = _core_cell(
        context, payload, "Phi", "radial_0_axial_1", rho
    )
    selectors["A_Z"] = "result.payload.ordinary_jets.Phi.radial_0_axial_1." + rho + ".enclosure"
    for name in ("Q", "Q_rho", "Q_Z"):
        values[name], _ = _recovery_cell(context, payload, name, rho)
        selectors[name] = "result.payload.local_core_recovery." + name + ".values." + rho + ".enclosure"
    return {
        "values": values,
        "source_selectors": selectors,
        "all_values_from_complete_enclosures": True,
        "tail_values_not_reused_as_function_values": True,
        "Q_scope": "saved Z=0 local core recovery only; no global Q=-M_Z extension",
    }


def _source_public_record(context: Any, source: Mapping[str, Any]) -> dict[str, Any]:
    """Serialize internal interval objects without selecting endpoints."""

    return {
        "values": {
            key: _interval_record(context, value)
            for key, value in source["values"].items()
        },
        "source_selectors": _json_value(source["source_selectors"]),
        "all_values_from_complete_enclosures": True,
        "tail_values_not_reused_as_function_values": True,
        "Q_scope": source["Q_scope"],
    }


def _geometry_packet(
    context: Any,
    rho_text: str,
    tau_text: str,
    epsilon: Any,
    delta: Any,
) -> dict[str, Any]:
    rho = _exact_fraction(context, rho_text)
    tau = _exact_fraction(context, tau_text)
    lam = context.sqrt(tau)
    one = context.mpf("1")
    a = one - delta
    kappa = lam ** a
    R = rho * epsilon
    radius_r = lam * context.sqrt(context.mpf("2") * R)
    return {
        "rho_input": rho_text,
        "tau_input": tau_text,
        "rho": _interval_record(context, rho),
        "tau": _interval_record(context, tau),
        "lambda": _interval_record(context, lam),
        "a_one_minus_delta": _interval_record(context, a),
        "kappa_lambda_to_a": _interval_record(context, kappa),
        "epsilon_core": _interval_record(context, epsilon),
        "R": _interval_record(context, R),
        "r": _interval_record(context, radius_r),
        "z": "0",
        "lambda_source": "sqrt(tau) exact fraction; no Lambda construction",
        "R_source": "rho*epsilon_core saved-data scaling; no old provider",
        "coordinate": "rho=Lambda*R, with R=rho*epsilon_core",
    }


def _instantaneous_streamline_slope(
    context: Any,
    source: Mapping[str, Any],
    geometry: Mapping[str, Any],
    logF0: Any,
) -> dict[str, Any]:
    """Formal local direction only; no material trajectory is integrated."""

    lam = _interval_from_record(context, geometry["lambda"], "lambda")
    values = source["values"]
    u_sign = _sign_status(context, values["Uz"])
    if u_sign != "positive":
        return {
            "kind": "local_instantaneous_direction",
            "definition": "dtheta/dz=ut/(r*uz)=F0*Phi/(lambda*U)",
            "status": "future_U_positive_uncertified",
            "ordinary_U_sign": u_sign,
            "ordinary_U_positive_required": True,
            "coefficient_interval": None,
            "logabs_interval": None,
            "exact_F0_atom_ref": "result.payload.core_amplitude_Z0.amplitude.logF0_Z0",
            "exp_materialized": False,
            "materialized_function_value": False,
            "exact_source_function_defined": True,
            "material_trajectory": False,
            "cumulative_winding": False,
        }
    coefficient = values["Phi"] / (lam * values["Uz"])
    result = _formal_F0_product(context, coefficient, logF0, "dtheta/dz")
    result.update(
        {
            "definition": "dtheta/dz=ut/(r*uz)=F0*Phi/(lambda*U)",
            "kind": "local_instantaneous_direction",
            "material_trajectory": False,
            "cumulative_winding": False,
            "ordinary_U_positive_required": True,
        }
    )
    return result


def _physical_values(
    context: Any,
    source: Mapping[str, Any],
    geometry: Mapping[str, Any],
    epsilon: Any,
    delta: Any,
    logF0: Any,
) -> dict[str, Any]:
    rho = _interval_from_record(context, geometry["rho"], "rho")
    lam = _interval_from_record(context, geometry["lambda"], "lambda")
    kappa = _interval_from_record(context, geometry["kappa_lambda_to_a"], "kappa")
    R = _interval_from_record(context, geometry["R"], "R")
    values = source["values"]
    one = context.mpf("1")
    ur = lam ** (-one) * context.sqrt(R / context.mpf("2")) * values["Q"]
    uz = lam ** (-one - delta) * values["Uz"]
    pressure = lam ** (-context.mpf("2") - context.mpf("2") * delta) * values["P_scaled"] / epsilon
    psi = kappa * R * values["Mz_over_R"]
    ut_coefficient = lam ** (-one - delta) * context.sqrt(context.mpf("2") * R) * values["Phi"]
    ut = _formal_F0_product(context, ut_coefficient, logF0, "ut")
    streamline_slope = _instantaneous_streamline_slope(
        context, source, geometry, logF0
    )
    return {
        "ur": _interval_record(context, ur),
        "uz": _interval_record(context, uz),
        "p": _interval_record(context, pressure),
        "psi": _interval_record(context, psi),
        "ut": ut,
        "instantaneous_streamline_slope": streamline_slope,
        "definitions": {
            "ur": "lambda^-1*sqrt(R/2)*Q",
            "uz": "lambda^(-1-delta)*U",
            "p": "lambda^(-2-2delta)*P_scaled/epsilon_core",
            "psi": "kappa*R*M",
            "ut": "F0*[lambda^(-1-delta)*sqrt(2R)*Phi]",
            "dtheta_dz": "F0*Phi/(lambda*U); local instantaneous direction only",
        },
        "Q_scope": "saved Z=0 Q only; global Q=-[d*M_Z+(1-delta)Z*M-2Z*U]/L is not evaluated here",
        "F0_exp_materialized": False,
        "F0_zero": False,
    }


def _cylindrical_jets(
    context: Any,
    source: Mapping[str, Any],
    geometry: Mapping[str, Any],
    epsilon: Any,
    delta: Any,
    logF0: Any,
    logF0_Z: Any,
) -> dict[str, Any]:
    rho = _interval_from_record(context, geometry["rho"], "rho")
    lam = _interval_from_record(context, geometry["lambda"], "lambda")
    R = _interval_from_record(context, geometry["R"], "R")
    values = source["values"]
    one = context.mpf("1")
    two = context.mpf("2")
    ur_r = lam ** (-two) * (values["Q"] / two + rho * values["Q_rho"])
    ur_over_r = lam ** (-two) * values["Q"] / two
    ur_z = lam ** (-two + delta) * context.sqrt(R / two) * values["Q_Z"]
    uz_r = lam ** (-two - delta) * context.sqrt(two * R) / epsilon * values["U_rho"]
    uz_z = lam ** (-two) * values["U_Z"]
    ut_r_coeff = lam ** (-two - delta) * (values["Phi"] + two * rho * values["A_rho"])
    ut_over_r_coeff = lam ** (-two - delta) * values["Phi"]
    ut_z_coeff = lam ** (-two) * context.sqrt(two * R) * (
        values["A_Z"] + logF0_Z * values["Phi"]
    )
    omega_theta = ur_z - uz_r
    omega_r = _formal_F0_product(context, -ut_z_coeff, logF0, "omega_r=-ut_z")
    omega_z_coeff = two * lam ** (-two - delta) * (
        values["Phi"] + rho * values["A_rho"]
    )
    omega_z = _formal_F0_product(context, omega_z_coeff, logF0, "omega_z")
    divergence = ur_r + ur_over_r + uz_z
    return {
        "ur_r": _interval_record(context, ur_r),
        "ur_over_r": _interval_record(context, ur_over_r),
        "ur_z": _interval_record(context, ur_z),
        "uz_r": _interval_record(context, uz_r),
        "uz_z": _interval_record(context, uz_z),
        "ut_r": _formal_F0_product(context, ut_r_coeff, logF0, "ut_r"),
        "ut_over_r": _formal_F0_product(context, ut_over_r_coeff, logF0, "ut_over_r"),
        "ut_z": _formal_F0_product(context, ut_z_coeff, logF0, "ut_z"),
        "omega_theta": _interval_record(context, omega_theta),
        "omega_r": omega_r,
        "omega_z": omega_z,
        "local_angular_velocity": _formal_F0_product(
            context, ut_over_r_coeff, logF0, "ut/r"
        ),
        "divergence_diagnostic": {
            "interval": _interval_record(context, divergence),
            "definition": "ur_r+ur_over_r+uz_z",
            "diagnostic_only": True,
            "not_a_divergence_proof": True,
            "not_removed_by_contains_zero": True,
        },
        "definitions": {
            "ur_r": "lambda^-2*(Q/2+rho*Qrho)",
            "ur_over_r": "lambda^-2*Q/2",
            "ur_z": "lambda^(-2+delta)*sqrt(R/2)*QZ",
            "uz_r": "lambda^(-2-delta)*sqrt(2R)/epsilon_core*U_rho",
            "uz_z": "lambda^-2*U_Z",
            "ut_r": "F0*lambda^(-2-delta)*(Phi+2rho*Phi_rho)",
            "ut_over_r": "F0*lambda^(-2-delta)*Phi",
            "ut_z": "F0*lambda^-2*sqrt(2R)*(Phi_Z+logF0_Z*Phi)",
            "local_angular_velocity": "ut/r=F0*lambda^(-2-delta)*Phi",
            "omega_theta": "ur_z-uz_r",
            "omega_r": "-ut_z",
            "omega_z": "F0*2lambda^(-2-delta)*(Phi+rho*Phi_rho)",
        },
        "Q_scope": "saved Z=0 Q/Qrho/QZ only; no full-Z symbolic pullback",
    }


def _cartesian_values(
    context: Any,
    values: Mapping[str, Any],
) -> dict[str, Any]:
    # ``_physical_values`` already publishes ``ur`` as a complete interval
    # record.  Reusing that record avoids feeding a JSON mapping back through
    # the interval context (and avoids selecting an endpoint).
    ur_record = _json_value(values["ur"])
    formal_ut = values["ut"]
    # Rebuild the negative formal product from its saved interval operands so
    # that coefficient sign, zero semantics, and exact F0 provenance agree
    # with the actual Cartesian component at theta=pi/2.
    coefficient = _interval_from_record(
        context, formal_ut["coefficient_interval"], "ut coefficient"
    )
    logF0 = _interval_from_record(
        context, formal_ut["logF0_interval"], "ut logF0"
    )
    negative_ut = _formal_F0_product(
        context, -coefficient, logF0, "Cartesian u at theta=pi/2"
    )
    negative_ut["rotation_multiplier"] = -1
    negative_ut["rotation_already_applied_to_coefficient"] = True
    negative_ut["rotation_source_expression"] = "-ut_coefficient_before_rotation"
    return {
        "theta_0": {
            "u": ur_record,
            "v": formal_ut,
            "w": values["uz"],
            "rotation": {"cos": "1", "sin": "0"},
        },
        "theta_pi_over_2": {
            "u": negative_ut,
            "v": ur_record,
            "w": values["uz"],
            "rotation": {"cos": "0", "sin": "1"},
            "rotation_multiplier": -1,
            "rotation_already_applied_to_coefficient": True,
        },
        "any_theta_definition": {
            "u": "ur*cos(theta)-ut*sin(theta)",
            "v": "ur*sin(theta)+ut*cos(theta)",
            "w": "uz",
            "theta_domain": "symbolic exact rotation; no approximate cos/sin evaluation",
        },
        "mixed_numeric_and_formal_ut_preserved": True,
        "ut_not_set_to_zero": True,
    }


def _growth_factor(
    context: Any,
    tau_a: str,
    tau_b: str,
    delta: Any,
    geometry_a: Mapping[str, Any],
    geometry_b: Mapping[str, Any],
    source_nonzero: Mapping[str, bool],
) -> dict[str, Any]:
    ta = _exact_fraction(context, tau_a)
    tb = _exact_fraction(context, tau_b)
    tau_ratio = tb / ta
    lam_a = _interval_from_record(context, geometry_a["lambda"], "lambda_a")
    lam_b = _interval_from_record(context, geometry_b["lambda"], "lambda_b")
    kap_a = _interval_from_record(context, geometry_a["kappa_lambda_to_a"], "kappa_a")
    kap_b = _interval_from_record(context, geometry_b["kappa_lambda_to_a"], "kappa_b")
    lam_ratio = lam_b / lam_a
    kap_ratio = kap_b / kap_a
    one = context.mpf("1")
    two = context.mpf("2")

    def positive_excess_certificate() -> dict[str, Any]:
        """Enclose exp(x)-1 stably for x=-(delta/2)log(tau_b/tau_a)."""

        x = -(delta / two) * context.log(tau_ratio)
        x_lower = _endpoint(context, x)
        x_upper = _endpoint(context, x, True)
        zero = context.mpf("0")
        one_local = context.mpf("1")
        x_lower_point = context.mpf(_endpoint_text(context, x_lower))
        x_upper_point = context.mpf(_endpoint_text(context, x_upper))
        if x_lower_point <= zero or x_upper_point >= one_local:
            return {
                "status": "positive_excess_assumption_not_certified",
                "x": _interval_record(context, x),
                "one_plus_positive_excess": None,
                "theorem": "requires delta>0 and 0<x<1",
                "delta_strict_positive_assumption": True,
                "delta_source": "saved delta interval",
            }
        quotient = x / (one_local - x)
        upper_excess = _endpoint(context, quotient, True)
        excess = context.make_mpf((x_lower._mpi_[0], upper_excess._mpi_[1]))
        one_plus = one_local + excess
        return {
            "status": "positive_excess_certified",
            "x": _interval_record(context, x),
            "excess_interval": _interval_record(context, excess),
            "one_plus_positive_excess": _interval_record(context, one_plus),
            "theorem": "x <= exp(x)-1 <= x/(1-x) for 0<x<1",
            "delta_strict_positive_assumption": True,
            "delta_source": "saved delta interval",
            "exp_difference_not_evaluated": True,
            "direct_exp_x_minus_one_subtraction": False,
        }

    excess_certificate = positive_excess_certificate()
    excess_factor = excess_certificate.get("one_plus_positive_excess")
    stable_excess_interval = None
    if excess_factor is not None:
        stable_excess_interval = _interval_from_record(
            context, excess_factor, "one plus positive excess"
        )

    def stable_delta_factor(base_factor: Any, direct_factor: Any) -> Any:
        if stable_excess_interval is None:
            return direct_factor
        return base_factor * stable_excess_interval

    def identity_factor(
        factor: Any,
        identity: str,
        source_name: str | None = None,
    ) -> dict[str, Any]:
        certified = True if source_name is None else bool(source_nonzero.get(source_name))
        return {
            "factor": _interval_record(context, factor) if certified else None,
            "candidate_factor": _interval_record(context, factor),
            "identity": identity,
            "source_nonzero_certified": certified,
            "signed_uncertain": not certified if source_name is not None else False,
            "imposed_exact_scale_identity": certified,
            "status": "imposed_exact_scale_identity" if certified else "future_uncertified_nonzero",
        }

    return {
        "tau_pair": {"from": tau_a, "to": tau_b},
        "tau_ratio": _interval_record(context, tau_ratio),
        "radius": identity_factor(lam_ratio, "lambda_b/lambda_a"),
        "axialscale": identity_factor(kap_ratio, "kappa_b/kappa_a"),
        "aspect": {
            "factor": _interval_record(context, stable_excess_interval)
            if stable_excess_interval is not None
            else None,
            "direct_interval_not_used_for_claim": _interval_record(
                context, kap_ratio / lam_ratio
            ),
            "identity": "(kappa_b/kappa_a)/(lambda_b/lambda_a)=exp(x)",
            "imposed_exact_scale_identity": stable_excess_interval is not None,
            "status": "imposed_exact_scale_identity"
            if stable_excess_interval is not None
            else "future_uncertified_positive_excess",
            "positive_excess_certificate": excess_certificate,
        },
        "aspect_excess": excess_certificate,
        "ur_magnitude": identity_factor(
            tau_ratio ** (-one / context.mpf("2")),
            "(tau_b/tau_a)^(-1/2)",
            "Q",
        ),
        "uz_and_ut_magnitude_growth": identity_factor(
            stable_delta_factor(
                tau_ratio ** (-one / two),
                tau_ratio ** (-(one + delta) / two),
            ),
            "each of |uz| and |ut| scales as (tau_b/tau_a)^(-(1+delta)/2)",
            "U_and_Phi",
        ),
        "uz_over_ut": identity_factor(
            context.mpf("1"),
            "uz/ut is time-invariant because both carry lambda^(-1-delta)",
            "U_and_Phi",
        ),
        "omega_r": identity_factor(
            tau_ratio ** (-one),
            "(tau_b/tau_a)^(-1)",
            "omega_r",
        ),
        "omega_z": identity_factor(
            stable_delta_factor(
                tau_ratio ** (-one),
                tau_ratio ** (-(two + delta) / two),
            ),
            "(tau_b/tau_a)^(-(2+delta)/2)",
            "omega_z",
        ),
        "same_source_atoms": {
            "cancelled_without_interval_subtraction": [
                "F0",
                "Q",
                "U",
                "Phi",
                "A_Z+logF0_Z*A",
                "A+rho*A_rho",
            ],
            "nonzero_source_certification_required": True,
            "source_nonzero_status": _json_value(source_nonzero),
            "delta_speed_factor_ref": "positive_excess_certificate.one_plus_positive_excess",
            "per_quantity": {
                "ur_magnitude": ["Q"],
                "uz_and_ut_magnitude_growth": ["U", "Phi", "F0"],
                "uz_over_ut": ["U", "Phi", "F0"],
                "omega_r": ["A_Z+logF0_Z*A", "F0"],
                "omega_z": ["A+rho*A_rho", "F0"],
            },
            "not_an_independent_measurement": True,
        },
        "omega_theta": {
            "status": "direct_value_ratio_not_single_exponent",
            "fit": None,
            "future_selector": "packets.<rho>.<tau>.jets.omega_theta",
        },
    }


def _validate_core_shape(context: Any, payload: Mapping[str, Any]) -> None:
    for field in ("Phi", "Uz", "Mz_over_R"):
        for jet in JET_KEYS:
            for rho in RHO_POINTS:
                _core_cell(context, payload, field, jet, rho)
    p_field = "P_scaled" if _is_mapping(payload.get("ordinary_jets", {}).get("P_scaled")) else "P"
    for jet in JET_KEYS:
        for rho in RHO_POINTS:
            _core_cell(context, payload, p_field, jet, rho)
    for name in ("Q", "Q_rho", "Q_Z"):
        for rho in RHO_POINTS:
            _recovery_cell(context, payload, name, rho)
    mean = payload.get("mean_minus_axial")
    if not _is_mapping(mean) or not _is_mapping(mean.get("values")):
        raise BindingInputError("mean_minus_axial values are missing")
    for jet in JET_KEYS:
        for rho in RHO_POINTS:
            cell = mean["values"].get(jet, {}).get(rho)
            if not _is_mapping(cell) or "enclosure" not in cell:
                raise BindingInputError("mean_minus_axial enclosure is missing")


def derive(
    data: Mapping[str, Any],
    emit: Callable[[str, Mapping[str, Any]], Any] | None = None,
) -> dict[str, Any]:
    """Bind saved core enclosures to physical Z=0 geometry once.

    This routine intentionally does not construct Lambda, N, F0, a trajectory,
    or a full physical field.  A malformed capsule raises after preserving all
    callback events already emitted; counters are never reset.
    """

    if emit is not None and not callable(emit):
        raise BindingInputError("emit must be callable or None")
    events: list[dict[str, Any]] = []
    _emit(
        emit,
        events,
        "family_entered",
        {
            "entry_before_interval_math": True,
            "source_result_sha256": SOURCE_RESULT_SHA256,
            "source_only": True,
            "N_materialized": False,
        },
    )
    try:
        if not _is_mapping(data):
            raise BindingInputError("data must be a mapping")
        capsule, payload = _payload_from_data(data)
        if capsule.get("schema") != CAPSULE_SCHEMA:
            raise BindingInputError("unexpected physical binding capsule schema")
        source_binding = _validate_source(capsule, payload)
        capsule_guards = _validate_capsule_guards(capsule)
        capsule_grid = _validate_grid(capsule)

        _emit(
            emit,
            events,
            "physical_geometry_entered",
            {"rho_points": list(RHO_POINTS), "tau_points": list(TAU_POINTS), "Z": "0"},
        )
        _emit(
            emit,
            events,
            "physical_values_entered",
            {"F0_formal": True, "exp_materialized": False, "point_values_selected": False},
        )
        _emit(
            emit,
            events,
            "physical_jets_entered",
            {"cylindrical_first_order_jets": True, "vorticity": ["omega_r", "omega_theta", "omega_z"]},
        )
        _emit(
            emit,
            events,
            "common_source_growth_entered",
            {"tau_ratio_identities": True, "trajectory": False, "cumulative_winding": False},
        )

        context = _new_interval_context()
        _validate_core_shape(context, payload)
        epsilon = _parameter_interval(context, payload, "epsilon_core")
        delta = _parameter_interval(context, payload, "delta")
        logF0 = _interval_from_record(
            context,
            _path(payload, ("core_amplitude_Z0", "amplitude", "logF0_Z0"), "logF0_Z0"),
            "logF0_Z0",
        )
        gauge = _path(payload, ("core_amplitude_Z0", "ordinary_Z0_gauge"), "ordinary Z0 gauge")
        logF0_Z = _interval_from_record(
            context,
            _lookup(gauge, ("saved_ordinary_DZ_logF0", "logF0_Z1"), "ordinary DZ logF0"),
            "ordinary DZ logF0",
        )
        logF0_ZZ = _interval_from_record(
            context,
            _lookup(gauge, ("saved_ordinary_DZZ_logF0", "logF0_Z2"), "ordinary DZZ logF0"),
            "ordinary DZZ logF0",
        )

        geometries: dict[str, dict[str, Any]] = {}
        for rho in RHO_POINTS:
            geometries[rho] = {}
            for tau in TAU_POINTS:
                geometries[rho][tau] = _geometry_packet(
                    context, rho, tau, epsilon, delta
                )
        _emit(
            emit,
            events,
            "physical_geometry_complete",
            {"packet_count": len(RHO_POINTS) * len(TAU_POINTS), "coordinate": "rho=Lambda*R; Z=0"},
        )

        packets: dict[str, dict[str, Any]] = {}
        for rho in RHO_POINTS:
            source = _source_records(context, payload, rho)
            packets[rho] = {}
            for tau in TAU_POINTS:
                geometry = geometries[rho][tau]
                values = _physical_values(
                    context, source, geometry, epsilon, delta, logF0
                )
                packets[rho][tau] = {
                    "geometry": geometry,
                    "source": source,
                    "values": values,
                    "cartesian": _cartesian_values(context, values),
                }
        _emit(
            emit,
            events,
            "physical_values_complete",
            {"packet_count": len(RHO_POINTS) * len(TAU_POINTS), "F0_exp_materialized": False},
        )

        for rho in RHO_POINTS:
            source = packets[rho][TAU_POINTS[0]]["source"]
            for tau in TAU_POINTS:
                packets[rho][tau]["jets"] = _cylindrical_jets(
                    context,
                    source,
                    geometries[rho][tau],
                    epsilon,
                    delta,
                    logF0,
                    logF0_Z,
                )
                packets[rho][tau]["logF0_ZZ_saved"] = _interval_record(context, logF0_ZZ)
        _emit(
            emit,
            events,
            "physical_jets_complete",
            {"packet_count": len(RHO_POINTS) * len(TAU_POINTS), "divergence_is_diagnostic": True},
        )

        growth: dict[str, dict[str, Any]] = {}
        for rho in RHO_POINTS:
            growth[rho] = {}
            for index, tau_a in enumerate(TAU_POINTS):
                for tau_b in TAU_POINTS[index + 1:]:
                    key = tau_a + "->" + tau_b
                    base_packet = packets[rho][tau_a]
                    base_source_values = base_packet["source"]["values"]
                    base_jets = base_packet["jets"]
                    source_nonzero = {
                        "F0": True,
                        "Q": _sign_status(context, base_source_values["Q"])
                        in {"positive", "negative"},
                        "U_and_Phi": (
                            _sign_status(context, base_source_values["Uz"])
                            in {"positive", "negative"}
                            and _sign_status(context, base_source_values["Phi"])
                            in {"positive", "negative"}
                        ),
                        "omega_r": base_jets["omega_r"].get("coefficient_sign")
                        in {"positive", "negative"},
                        "omega_z": base_jets["omega_z"].get("coefficient_sign")
                        in {"positive", "negative"},
                    }
                    growth[rho][key] = _growth_factor(
                        context,
                        tau_a,
                        tau_b,
                        delta,
                        geometries[rho][tau_a],
                        geometries[rho][tau_b],
                        source_nonzero,
                    )
        _emit(
            emit,
            events,
            "common_source_growth_complete",
            {
                "rho_count": len(RHO_POINTS),
                "identity_kind": "imposed_exact_scale_identity_when_certified",
                "stable_positive_excess_certificate": True,
                "omega_theta_single_exponent_fit": False,
                "all_quantities_exact": False,
            },
        )

        # Publish shallow packet copies with complete endpoint records for the
        # private source intervals.  Keep ``packets`` private until all growth
        # calculations finish; no endpoint is selected as a value.
        published_packets: dict[str, dict[str, Any]] = {}
        for rho in RHO_POINTS:
            public_source = _source_public_record(
                context, packets[rho][TAU_POINTS[0]]["source"]
            )
            published_packets[rho] = {}
            for tau in TAU_POINTS:
                packet = dict(packets[rho][tau])
                packet["source"] = public_source
                published_packets[rho][tau] = packet

        proof = {
            "kind": "saved_core_physical_geometry_binding_v1",
            "source_binding": source_binding,
            "capsule_guards": capsule_guards,
            "capsule_grid": capsule_grid,
            "all_ordinary_values_from_complete_enclosures": True,
            "F0_formal_atom_only": True,
            "no_Lambda_or_N_reconstruction": True,
            "Lambda_scaling": "Lambda=1/epsilon_core source identity only; not materialized",
            "no_old_provider_or_checker": True,
            "coordinate": "rho=Lambda*R; Z=0",
            "Q_scope": "saved Z=0 local recovery only; do not extend Q=-M_Z to full Z",
            "theta_rotation": "u=ur*cos(theta)-ut*sin(theta), v=ur*sin(theta)+ut*cos(theta)",
            "vorticity": {
                "omega_theta": "ur_z-uz_r",
                "omega_r": "-ut_z",
                "omega_z": "ut_r+ut_over_r",
            },
            "growth": {
                "same_source_atoms_cancelled": [
                    "F0",
                    "Q",
                    "U",
                    "Phi",
                    "A_Z+logF0_Z*A",
                    "A+rho*A_rho",
                ],
                "uz_and_ut_magnitude_growth": "separate lambda^(-1-delta) magnitude factor",
                "uz_over_ut": "time-invariant factor 1 when U and Phi are certified nonzero",
                "per_quantity_atoms": {
                    "ur_magnitude": ["Q"],
                    "uz_and_ut_magnitude_growth": ["U", "Phi", "F0"],
                    "omega_r": ["A_Z+logF0_Z*A", "F0"],
                    "omega_z": ["A+rho*A_rho", "F0"],
                },
                "identities_are_not_independent_measurements": True,
                "omega_theta_single_exponent_fit": False,
            },
            "future": {
                "global_streamline_integration": True,
                "cumulative_winding": True,
                "local_streamline_slope": (
                    "defined per packet only when source U is certified nonzero; "
                    "no integrated curve"
                ),
            },
            "preservation": _json_value(payload.get("preservation")),
        }
        _emit(
            emit,
            events,
            "proof_bundle",
            {"proof_count": 1, "physical_field_materialized": False, "full_goal_complete": False},
        )
        _emit(
            emit,
            events,
            "family_completed",
            {"proof_count": 1, "source_only": True, "metadata_only": False},
        )

        payload_out = {
            "component": COMPONENT,
            "source_only": True,
            "numeric_evaluation": True,
            "metadata_only": False,
            "source_binding": source_binding,
            "capsule_guards": capsule_guards,
            "capsule_grid": capsule_grid,
            "source_identity": _json_value(payload.get("source_identity")),
            "parameters": {
                "epsilon_core": _interval_record(context, epsilon),
                "delta": _interval_record(context, delta),
                "logF0_Z0": _interval_record(context, logF0),
                "ordinary_DZ_logF0": _interval_record(context, logF0_Z),
                "ordinary_DZZ_logF0": _interval_record(context, logF0_ZZ),
                "F0_materialized": False,
            },
            "geometries": geometries,
            "packets": published_packets,
            "common_source_growth": growth,
            "coordinate": "rho=Lambda*R; Z=0",
        "Lambda_scaling": "Lambda=1/epsilon_core is a saved source identity; no Lambda value is constructed",
            "rho_points": list(RHO_POINTS),
            "tau_points": list(TAU_POINTS),
            "theta_cases": _json_value(THETA_CASES),
            "units": _json_value(payload.get("row_units")),
            "axis_packets": _json_value(payload.get("axis_packets")),
            "preservation": _json_value(payload.get("preservation")),
            "proof": proof,
            "events": events,
            "runtime_model_metadata": None,
            "no_physical_huge_field_materialized": True,
            "no_NS_claim": True,
        }
        result_out = {
            "schema": SCHEMA,
            "component": COMPONENT,
            "passed": True,
            "metadata_only": False,
            "source_only": True,
            "numeric_evaluation": True,
            "runtime_model_metadata": None,
            "physical_core_packets_evaluated": True,
            "physical_layer_completed": False,
            "full_physical_field_materialized": False,
            "full_goal_complete": False,
            "proof_count": 1,
            "actual_counts": dict(SUCCESS_COUNTS),
            "expected_counts": dict(EXPECTED_COUNTS),
            "new_flags": dict(SUCCESS_FLAGS),
            "flags": dict(SUCCESS_FLAGS),
            "expected_flags": dict(EXPECTED_FLAGS),
            "api": dict(API),
            "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
            "count_event_map": {key: list(value) for key, value in COUNT_EVENT_MAP.items()},
            "payload": payload_out,
        }
        # Metadata-only serialization validation: all exact interval endpoint
        # tuples and formal records must remain JSON-shaped before publication.
        try:
            json.dumps(result_out, ensure_ascii=False, allow_nan=False)
        except (TypeError, ValueError) as error:
            raise BindingInputError(
                "physical binding result is not JSON-serializable"
            ) from error
        return result_out
    except BindingInputError:
        raise


__all__ = [
    "API",
    "BindingInputError",
    "CALLBACK_EVENT_COUNTS",
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "EXPECTED_COUNTS",
    "EXPECTED_FLAGS",
    "INTERVAL_DPS",
    "RHO_POINTS",
    "SCHEMA",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "TAU_POINTS",
    "derive",
]
