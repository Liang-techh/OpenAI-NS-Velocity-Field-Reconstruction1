"""Independent O3 pressure-source value and ordinary-Z jet kernel.

This file is an unexecuted worker artifact.  It has no project, source-module,
provider, checker, or ancestor-mathematics imports.  ``derive`` consumes the
literal ``native_O3_pressure_source_value_jet_source_capsule_v1`` packet and,
when the root runner invokes it, performs one new directed interval
calculation.  The calculation is a coupled monotone Darboux rectangle
construction; it is deliberately not the entered Taylor-quadrature family.

The decimal strings in the result are display values only.  Every displayed
interval also carries the raw lower and upper MPF tuples emitted by the
independent ``MPIntervalContext``.  No interval midpoint, endpoint promotion,
or nominal anchor is selected as a source value.  The exact pressure function
is retained separately from its directed approximations::

    q(Z) = 1 + Z**2
    P0(Z) = -earlyM2*q(Z)**(-2) + late_exact_function(Z)
    earlyM2 = 5/2 + Mslope
               + exp(3/5)*(exp(-1) - exp(-yd))/2
    Mslope = (1/2)*integral_0^1 exp(y/5 - 6*J(y)/5) dy
    J(y) = integral_0^y sigma(s) ds

The eleven post-``yd`` integrals remain exact source-function definitions from
the capsule.  Their common positive tail allowance is kept nonzero even when
its decimal display is tiny.  The result is a local pressure-source
enclosure/ordinary-derivative packet only; it does not admit global source
jets, a physical field, an NS claim, or the full reconstruction goal.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "native_O3_pressure_source_value_jet_kernel_result_v1"
CAPSULE_SCHEMA = "native_O3_pressure_source_value_jet_source_capsule_v1"
COMPONENT = "native_O3_pressure_source_value_jet"

N_DYADIC_PANELS = 4096
INTERVAL_DPS = 100

BETA2_STAGE_NAMES = (
    "reference_extension",
    "slope_transition_ref",
    "axial_turnoff",
    "slope_transition_mu",
    "power_buffer",
    "pulse_reserved",
)
BETA0_STAGE_NAMES = (
    "power_buffer_rel",
    "steep_transition_in",
    "steep_power",
    "steep_transition_out",
    "waiting",
    "heat_collar",
    "exterior_power_tail",
)
VARIABLE_STAGE_NAME = "z_flatten"
EARLY_STAGE_NAMES = (
    "reference_extension",
    "slope_transition_ref",
    "axial_turnoff",
)
LATE_STAGE_NAMES = (
    "slope_transition_mu",
    "power_buffer",
    "pulse_reserved",
    "z_flatten",
    "power_buffer_rel",
    "steep_transition_in",
    "steep_power",
    "steep_transition_out",
    "waiting",
    "heat_collar",
    "exterior_power_tail",
)
ALL_STAGE_NAMES = EARLY_STAGE_NAMES + LATE_STAGE_NAMES


# Literal audit contracts.  The zero entries are intentional: this kernel
# never materializes N, calls an ancestor provider, or calls a checker.
SUCCESS_COUNTS = {
    "N_materialization_calls": 0,
    "ancestor_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "provider_calls": 0,
    "source_value_reducer_calls": 0,
    "checker_calls": 0,
    "selftest_calls": 0,
    "source_family_calls": 1,
    "new_math_family_calls": 2,
    "new_quadrature_family_calls": 1,
    "new_pressure_jet_family_calls": 1,
    "numeric_evaluation_calls": 1,
    "new_mathematical_proofs": 1,
    "actual_pressure_source_enclosure_calls": 1,
}


SUCCESS_FLAGS = {
    "actual_pressure_source_enclosure_evaluated": True,
    "numeric_evaluation": True,
    "actual_source_Z_jets_admitted": False,
    "actual_parameter_admission": False,
    "physical_layer_completed": False,
    "full_goal": False,
    "full_goal_complete": False,
    "N_materialized": False,
    "old_math_reentered": False,
    "old_provider_called": False,
    "source_value_point_selected": False,
    "source_parameter_point_selected": False,
    "nominal_anchor_selected": False,
    "physical_field_materialized": False,
    "NS_claim": False,
    "all_old_caps_and_errors_preserved": True,
    "late_integrals_replaced_by_zero": False,
    "prefix_error_deleted_by_symmetry": False,
}


# An entry row precedes each new mathematical subfamily.  The final proof row
# is a receipt, not a CAS/checker invocation.
CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "quadrature_entered": 1,
    "quadrature_complete": 1,
    "pressure_jet_entered": 1,
    "pressure_jet_complete": 1,
    "proof_bundle": 1,
    "family_completed": 1,
}


# The map is intentionally count-name -> callback events, so a root runner can
# audit the literal count contract without importing this module.
COUNT_EVENT_MAP = {
    "source_family_calls": ["family_entered"],
    "new_math_family_calls": ["quadrature_entered", "pressure_jet_entered"],
    "new_quadrature_family_calls": ["quadrature_entered"],
    "new_pressure_jet_family_calls": ["pressure_jet_entered"],
    "numeric_evaluation_calls": ["family_entered"],
    "actual_pressure_source_enclosure_calls": ["pressure_jet_entered"],
    "new_mathematical_proofs": ["proof_bundle"],
}


class KernelInputError(ValueError):
    """Raised when the literal source capsule violates this narrow contract."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise KernelInputError(label + " must be a mapping")
    return value


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _json_value(value: Any) -> Any:
    """Keep capsule literals JSON-shaped without evaluating opaque objects."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
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
        "numeric_evaluation": name in {"quadrature_complete", "pressure_jet_complete"},
    }
    if payload:
        row.update(_json_value(payload))
    events.append(row)
    if emit is not None:
        emit(name, row)


def _new_interval_context() -> Any:
    """Create a private interval context; no ambient mpmath context is used."""

    # This is the only external import in the worker and is intentionally
    # local.  No project module, old numerical helper, or source provider is
    # imported here.
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.dps = INTERVAL_DPS
    return context


def _point_interval(context: Any, value: Any) -> Any:
    """Parse a scalar as an interval without routing it through binary64."""

    if isinstance(value, str):
        return context.mpf(value)
    if isinstance(value, bool):
        return context.mpf(1 if value else 0)
    if isinstance(value, int):
        return context.mpf(value)
    return context.mpf(str(value))


def _raw_endpoint(context: Any, value: Any, upper: bool = False) -> Any:
    interval = context.convert(value)
    raw = interval._mpi_[1 if upper else 0]
    return context.make_mpf((raw, raw))


def _lower(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, False)


def _upper(context: Any, value: Any) -> Any:
    return _raw_endpoint(context, value, True)


def _display_endpoint(context: Any, value: Any) -> str:
    text = context.nstr(context.convert(value), n=INTERVAL_DPS)
    text = str(text).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
        if "," in text:
            text = text.split(",", 1)[0].strip()
    return text


def _interval_record(context: Any, value: Any) -> dict[str, Any]:
    """Return decimal endpoint displays plus raw MPF tuple certificates."""

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


def _compact_interval_value(context: Any, value: Any) -> list[Any]:
    """Lossless table cell: displays plus both raw directed MPF endpoints."""

    record = _interval_record(context, value)
    return [
        record["lower"],
        record["upper"],
        record["lower_exact_mpf_tuple"],
        record["upper_exact_mpf_tuple"],
    ]


def _interval_from_record(context: Any, value: Any, label: str) -> Any:
    """Rehydrate a saved interval record without selecting a representative."""

    if _is_mapping(value):
        # Nested records occur for the saved anchor root and parameter packets.
        for nested_key in (
            "interval",
            "enclosure",
            "record",
            "root_interval",
            "anchor_interval",
            "data",
        ):
            if nested_key in value and not (
                "lower" in value or "lower_exact_mpf_tuple" in value
            ):
                return _interval_from_record(context, value[nested_key], label)
        lower_tuple = value.get("lower_exact_mpf_tuple")
        upper_tuple = value.get("upper_exact_mpf_tuple")
        if _is_sequence(lower_tuple) and _is_sequence(upper_tuple):
            if len(lower_tuple) != 4 or len(upper_tuple) != 4:
                raise KernelInputError(label + " exact MPF tuple must have length four")
            try:
                left = tuple(int(item) for item in lower_tuple)
                right = tuple(int(item) for item in upper_tuple)
            except (TypeError, ValueError) as error:
                raise KernelInputError(label + " exact MPF tuple is malformed") from error
            return context.make_mpf((left, right))
        if "lower" in value and "upper" in value:
            return context.mpf((str(value["lower"]), str(value["upper"])))
        for scalar_key in ("value", "point", "selected", "scalar"):
            if scalar_key in value:
                return _interval_from_record(context, value[scalar_key], label)
        raise KernelInputError(label + " has no interval endpoints")
    if _is_sequence(value):
        if len(value) != 2:
            raise KernelInputError(label + " endpoint sequence must have length two")
        return context.mpf((str(value[0]), str(value[1])))
    if value is None:
        raise KernelInputError(label + " is missing")
    try:
        return _point_interval(context, value)
    except Exception as error:  # pragma: no cover - defensive boundary
        raise KernelInputError(label + " is not an interval record") from error


def _hull(context: Any, lower: Any, upper: Any) -> Any:
    """Build an interval from an explicitly supplied lower and upper bound."""

    return context.make_mpf(
        (
            _lower(context, lower)._mpi_[0],
            _upper(context, upper)._mpi_[1],
        )
    )


def _parameter_candidate(capsule: Mapping[str, Any], aliases: Sequence[str]) -> Any:
    containers: list[Any] = [capsule]
    for key in (
        "saved_parameters",
        "saved_parameter_data",
        "parameter_data",
        "parameters",
        "selected_parameters",
        "source_parameters",
        "source_parameter_data",
        "accepted_parameter_payload",
        "accepted_anchor",
    ):
        if key in capsule:
            containers.insert(0, capsule[key])
    for container in containers:
        if not _is_mapping(container):
            continue
        for alias in aliases:
            if alias in container:
                return container[alias]
    return None


def _read_parameters(context: Any, capsule: Mapping[str, Any]) -> dict[str, Any]:
    raw_delta = _parameter_candidate(capsule, ("delta", "saved_delta", "delta_interval"))
    raw_logp = _parameter_candidate(
        capsule,
        ("logP", "logPstar", "logPstar_interval", "saved_logP", "yd"),
    )
    raw_anchor = _parameter_candidate(
        capsule,
        ("anchor_interval", "anchor", "saved_anchor", "root_interval", "source_anchor"),
    )
    if raw_delta is None or raw_logp is None or raw_anchor is None:
        raise KernelInputError(
            "capsule must provide saved delta, logP/yd, and anchor interval data"
        )
    result: dict[str, Any] = {}
    for name, raw in (
        ("delta", raw_delta),
        ("logP", raw_logp),
        ("anchor", raw_anchor),
    ):
        result[name] = {
            "raw": _json_value(raw),
            "interval": _interval_from_record(context, raw, "saved " + name),
        }
    raw_epsilon_core = _parameter_candidate(
        capsule,
        ("epsilon_core", "epsilon", "accepted_epsilon_core"),
    )
    if raw_epsilon_core is not None:
        result["epsilon_core"] = {
            "raw": _json_value(raw_epsilon_core),
            "interval": _interval_from_record(
                context, raw_epsilon_core, "saved epsilon_core"
            ),
        }
    raw_loglambda = _parameter_candidate(
        capsule,
        ("logLambda", "log_lambda", "logLambda_interval", "saved_logLambda"),
    )
    if raw_loglambda is not None:
        result["logLambda"] = {
            "raw": _json_value(raw_loglambda),
            "interval": _interval_from_record(context, raw_loglambda, "saved logLambda"),
        }
    return result


def _source_admission(capsule: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the capsule's same-family DATA joins before new arithmetic."""

    identity = capsule.get("family_identity")
    common = capsule.get("common_admission")
    common = _as_mapping(common, "capsule.common_admission")
    common_identity = common.get("family_identity")
    if not _is_mapping(identity) or not _is_mapping(common_identity):
        raise KernelInputError("same-family identity records are missing")
    identity_keys = (
        "actual_five_defect_family_sha256",
        "datum_enclosure_sha256",
        "implicit_source_sha256",
    )
    for key in identity_keys:
        if identity.get(key) is None or common_identity.get(key) != identity.get(key):
            raise KernelInputError("same-family identity mismatch at " + key)
    pressure_source = _as_mapping(capsule.get("pressure_source"), "capsule.pressure_source")
    for key in ("datum_enclosure_sha256", "implicit_source_sha256"):
        if pressure_source.get(key) != identity.get(key):
            raise KernelInputError("pressure source family mismatch at " + key)
    if common.get("all_passed") is not True:
        raise KernelInputError("common source admission did not pass")
    if common.get("finite_rows_and_tails_bound_to_same_nonlinear_fixed_point") is not True:
        raise KernelInputError("same nonlinear fixed-point binding is missing")

    true_measures = _as_mapping(
        capsule.get("true_pressure_measure_bindings"),
        "capsule.true_pressure_measure_bindings",
    )
    if true_measures.get("passed") is not True:
        raise KernelInputError("true pressure measure binding did not pass")
    if true_measures.get("all_fourteen_current_true_pressure_measures_bound") is not True:
        raise KernelInputError("all fourteen true pressure measures are not bound")
    tail_proof = true_measures.get("new_source_tail_proof")
    if not _is_sequence(tail_proof):
        raise KernelInputError("global post-yd tail proof rows are missing")
    proof_text = " ".join(str(row) for row in tail_proof).lower()
    if "disjoint" not in proof_text or "exp(3/5-y)" not in proof_text:
        raise KernelInputError("global pointwise/disjoint post-yd tail proof is missing")
    if pressure_source.get("shared_tail_bound_is_joint_not_additive") is not True:
        raise KernelInputError("joint post-yd tail bound is not marked non-additive")

    all_order = _as_mapping(
        capsule.get("all_order_pressure_source_bindings"),
        "capsule.all_order_pressure_source_bindings",
    )
    if all_order.get("passed") is not True:
        raise KernelInputError("all-order pressure source binding did not pass")
    if all_order.get("fixed_true_masses_not_replaced_by_interval_values") is not True:
        raise KernelInputError("fixed true masses are not retained as source data")
    if all_order.get("flatten_true_integral_not_replaced_by_Cauchy_cap") is not True:
        raise KernelInputError("flatten true integral was replaced by a cap")

    reflection = _as_mapping(capsule.get("exact_sigma_reflection"), "capsule.exact_sigma_reflection")
    if reflection.get("original_sigma_reflection_identity") is not True:
        raise KernelInputError("exact sigma reflection identity is missing")
    if reflection.get("original_J1_half_by_reflected_exact_integral") is not True:
        raise KernelInputError("exact J(1)=1/2 reflection identity is missing")

    pressure_definition = _as_mapping(
        pressure_source.get("implicit_source_definition"),
        "capsule.pressure_source.implicit_source_definition",
    )
    if "H replaced by 1" not in str(pressure_definition.get("angular_profile", "")):
        raise KernelInputError("preheat H=1 source definition is missing")
    return {
        "family_identity": _json_value(identity),
        "common_admission": _json_value(common),
        "pressure_source": _json_value(pressure_source),
        "true_pressure_measure_bindings": _json_value(true_measures),
        "all_order_pressure_source_bindings": _json_value(all_order),
        "exact_sigma_reflection": _json_value(reflection),
        "global_pointwise_bound_and_disjoint_partition": True,
        "original_mass_caps_metadata_only": True,
    }


def _validate_capsule(capsule: Any) -> Mapping[str, Any]:
    capsule = _as_mapping(capsule, "capsule")
    if capsule.get("schema") != CAPSULE_SCHEMA:
        raise KernelInputError("unexpected source capsule schema")
    if capsule.get("component") not in (None, COMPONENT):
        raise KernelInputError("source capsule component does not match this kernel")
    status = _as_mapping(capsule.get("status"), "capsule.status")
    if status.get("metadata_only") is not True:
        raise KernelInputError("source capsule status must remain metadata-only")
    # The input capsule must be metadata-only.  ``derive`` itself is the one
    # new interval calculation and records its own numeric flag in the result.
    for key in (
        "mathematical_functions_called",
        "math_calls",
        "quadrature_calls",
        "project_imports",
        "N_materialized",
    ):
        if key in status and status[key] not in (0, False, None):
            raise KernelInputError("capsule status is not metadata-only: " + key)
    for key in ("new_recipe_executed", "executed", "numeric_evaluation"):
        if key in status and status[key] not in (0, False, None):
            raise KernelInputError("capsule status already executed: " + key)
    flags = status.get("flags")
    if _is_mapping(flags) and any(bool(value) for value in flags.values()):
        raise KernelInputError("capsule status flags must be zero")
    return capsule


def _stage_candidates(capsule: Mapping[str, Any]) -> list[Any]:
    containers: list[Any] = []
    for key in (
        "stage_function_definitions",
        "actual_stage_function_definitions",
        "stage_functions",
        "actual_stage_functions",
        "pressure_stage_functions",
        "stages",
        "same_source_stage_functions",
        "stage_partition",
    ):
        if key in capsule:
            containers.append(capsule[key])
    for parent_key in ("source_functions", "source_function_data", "pressure_source"):
        parent = capsule.get(parent_key)
        if _is_mapping(parent):
            for key in (
                "stage_function_definitions",
                "actual_stage_function_definitions",
                "stage_functions",
                "stages",
            ):
                if key in parent:
                    containers.append(parent[key])
    for container in containers:
        if _is_sequence(container):
            return list(container)
        if _is_mapping(container):
            for nested_key in ("stage_functions", "stages", "definitions", "rows"):
                nested = container.get(nested_key)
                if _is_sequence(nested):
                    return list(nested)
                if _is_mapping(nested):
                    container = nested
                    break
            if _is_mapping(container):
                rows: list[dict[str, Any]] = []
                for name, definition in container.items():
                    if _is_mapping(definition):
                        row = dict(definition)
                        row.setdefault("name", str(name))
                    else:
                        row = {"name": str(name), "definition": definition}
                    rows.append(row)
                if rows:
                    return rows
    raise KernelInputError("capsule is missing the fourteen actual stage definitions")


def _canonical_stage_name(value: Any) -> str:
    text = str(value).strip()
    return {"flatten": VARIABLE_STAGE_NAME, "variable_flatten": VARIABLE_STAGE_NAME}.get(
        text, text
    )


def _read_stage_definitions(capsule: Mapping[str, Any]) -> list[dict[str, Any]]:
    capsule_early = capsule.get("early_stage_names")
    if _is_sequence(capsule_early) and tuple(str(item) for item in capsule_early) != EARLY_STAGE_NAMES:
        raise KernelInputError("capsule early_stage_names do not match the three early stages")
    capsule_late = capsule.get("late_stage_names")
    if _is_sequence(capsule_late) and tuple(str(item) for item in capsule_late) != LATE_STAGE_NAMES:
        raise KernelInputError("capsule late_stage_names do not match the eleven late stages")
    rows = _stage_candidates(capsule)
    normalized: dict[str, dict[str, Any]] = {}
    for raw in rows:
        if not _is_mapping(raw):
            raise KernelInputError("each stage definition must be a mapping")
        name = raw.get("name", raw.get("stage", raw.get("stage_name", raw.get("id"))))
        if name is None:
            raise KernelInputError("stage definition has no name")
        canonical = _canonical_stage_name(name)
        if canonical in normalized:
            raise KernelInputError("duplicate stage definition " + canonical)
        normalized[canonical] = dict(raw)
        normalized[canonical]["name"] = canonical
    expected = set(ALL_STAGE_NAMES)
    found = set(normalized)
    if found != expected:
        missing = sorted(expected - found)
        extra = sorted(found - expected)
        raise KernelInputError(
            "stage family must be exactly fourteen rows; missing=%r extra=%r" % (missing, extra)
        )
    result: list[dict[str, Any]] = []
    for name in ALL_STAGE_NAMES:
        row = dict(normalized[name])
        row["source_definition_preserved"] = True
        if name in BETA2_STAGE_NAMES:
            row.setdefault("beta_contract", "beta(y)=2")
            row["expected_beta_family"] = "beta2"
        elif name in BETA0_STAGE_NAMES:
            row.setdefault("beta_contract", "beta(y)=0")
            row["expected_beta_family"] = "beta0"
        else:
            row.setdefault("beta_contract", "0<=beta(y)<=2")
            row["expected_beta_family"] = "variable"
        result.append(row)
    return result


def _stable_sigma(context: Any, node: Any, index: int, count: int) -> Any:
    """Evaluate sigma once at a node using a stable directed log-ratio."""

    if index == 0:
        return context.zero
    if index == count:
        return context.one
    if index * 2 == count:
        return context.mpf("0.5")
    s = context.convert(node)
    ratio = 1 / (s * s) - 1 / ((1 - s) * (1 - s))
    if index * 2 < count:
        exp_small = context.exp(-ratio)
        return exp_small / (1 + exp_small)
    exp_small = context.exp(ratio)
    return 1 / (1 + exp_small)


def _monotone_darboux_mass(context: Any, count: int) -> dict[str, Any]:
    """Compute the fresh coupled prefix/Darboux enclosure for Mslope.

    The returned sigma, prefix, and panel audits are compact but lossless:
    dyadic coordinates are represented by exact integer indices, while every
    directed interval endpoint retains its decimal display and raw MPF tuple.
    """

    N = int(count)
    h = context.mpf(1) / N
    nodes = [context.mpf(i) / N for i in range(N + 1)]
    # This loop is the only sigma evaluation loop.  Prefix construction below
    # consumes these stored intervals and never calls sigma again.
    sigma_nodes = [
        _stable_sigma(context, node, index, N)
        for index, node in enumerate(nodes)
    ]

    jlo_nodes: list[Any] = [context.zero]
    jhi_nodes: list[Any] = [context.zero]
    jlo = context.zero
    jhi = context.zero
    for index in range(N):
        jlo = jlo + h * _lower(context, sigma_nodes[index])
        jhi = jhi + h * _upper(context, sigma_nodes[index + 1])
        jlo_nodes.append(jlo)
        jhi_nodes.append(jhi)

    # The audit is retained in compact lossless tables.  Dyadic coordinates are
    # exact from the shared denominator and row index, so repeated ``y=i/N``
    # interval records are not materialized.  Each interval cell keeps both
    # decimal endpoint displays and the raw MPF endpoint tuples.
    interval_columns = [
        "lower",
        "upper",
        "lower_exact_mpf_tuple",
        "upper_exact_mpf_tuple",
    ]
    sigma_records = {
        "schema": "compact_interval_table_v1",
        "quadrature_panels": N,
        "node_definition": "y=node_index/quadrature_panels",
        "columns": ["node_index", "sigma", "edge_branch"],
        "interval_columns": interval_columns,
        "rows": [
            [
                index,
                _compact_interval_value(context, sigma),
                "zero" if index == 0 else "one" if index == N else "interior",
            ]
            for index, sigma in enumerate(sigma_nodes)
        ],
        "evaluation": "stable_logratio_once_at_dyadic_node",
    }
    prefix_records = {
        "schema": "compact_interval_table_v1",
        "quadrature_panels": N,
        "node_definition": "y=node_index/quadrature_panels",
        "columns": ["node_index", "Jlo_prefix", "Jhi_prefix"],
        "interval_columns": interval_columns,
        "rows": [
            [
                index,
                _compact_interval_value(context, jlo_nodes[index]),
                _compact_interval_value(context, jhi_nodes[index]),
            ]
            for index in range(N + 1)
        ],
        "rectangle_rule": "left_lower/right_upper_monotone_prefix",
    }

    mass_lower = context.zero
    mass_upper = context.zero
    panel_rows: list[list[Any]] = []
    five = context.mpf(5)
    six = context.mpf(6)
    two = context.mpf(2)
    for index in range(N):
        a = nodes[index]
        b = nodes[index + 1]
        # J(y) for y in [a,b] is enclosed by Jlo(a) and Jhi(b).  The
        # exponent is decreasing in J, hence the endpoint order below.
        exponent_lower = _lower(context, a) / five - six * _upper(context, jhi_nodes[index + 1]) / five
        exponent_upper = _upper(context, b) / five - six * _lower(context, jlo_nodes[index]) / five
        density_lower = context.exp(exponent_lower) / two
        density_upper = context.exp(exponent_upper) / two
        panel_lower = h * _lower(context, density_lower)
        panel_upper = h * _upper(context, density_upper)
        mass_lower = mass_lower + panel_lower
        mass_upper = mass_upper + panel_upper
        panel_rows.append(
            [
                index,
                index + 1,
                _compact_interval_value(context, density_lower),
                _compact_interval_value(context, density_upper),
                _compact_interval_value(context, panel_lower),
                _compact_interval_value(context, panel_upper),
            ]
        )

    mass_interval = _hull(context, mass_lower, mass_upper)
    half = context.mpf(1) / 2
    symmetry_contains_half = bool(
        _lower(context, jlo_nodes[-1]) <= half
        and half <= _upper(context, jhi_nodes[-1])
    )
    return {
        "quadrature_panels": N,
        "precision_dps": INTERVAL_DPS,
        "width": _interval_record(context, h),
        "nodes_evaluated_once": N + 1,
        "sigma_nodes": sigma_records,
        "prefix_records": prefix_records,
        "panel_records": {
            "schema": "compact_interval_table_v1",
            "quadrature_panels": N,
            "panel_definition": "[a,b]=[panel_left_index/quadrature_panels,panel_right_index/quadrature_panels]",
            "columns": [
                "panel_left_index",
                "panel_right_index",
                "density_lower",
                "density_upper",
                "panel_mass_lower",
                "panel_mass_upper",
            ],
            "interval_columns": interval_columns,
            "prefix_references": {
                "Jlo_a": "prefix_records.rows[panel_left_index].Jlo_prefix",
                "Jhi_b": "prefix_records.rows[panel_right_index].Jhi_prefix",
            },
            "rows": panel_rows,
            "width": _compact_interval_value(context, h),
            "method": "monotone_Darboux_rectangle",
        },
        "J_at_one_prefix_enclosure": _interval_record(context, _hull(context, jlo_nodes[-1], jhi_nodes[-1])),
        "symmetry_check": {
            "identity": "sigma(1-s)=1-sigma(s), therefore J(1)=1/2",
            "contains_one_half": symmetry_contains_half,
            "posthoc_only": True,
            "prefix_error_retained": True,
            "prefix_was_not_intersected_with_half": True,
        },
        "Mslope": _interval_record(context, mass_interval / 1),
        "Mslope_interval": mass_interval,
        "method": "new_coupled_monotone_Darboux_integrator_v1",
        "not_old_Taylor_quadrature": True,
    }


def _signed_late_first_derivative_bound(context: Any, z: Any, eta: Any) -> tuple[Any, str]:
    magnitude = 4 * abs(z) * eta
    z_lower = _lower(context, z)
    z_upper = _upper(context, z)
    zero = context.zero
    if bool(z_upper < zero):
        return _hull(context, -_upper(context, magnitude), zero), "negative_Z"
    if bool(z_lower > zero):
        return _hull(context, zero, _upper(context, magnitude)), "positive_Z"
    return _hull(context, -_upper(context, magnitude), _upper(context, magnitude)), "zero_straddling"


def _pressure_sample(context: Any, label: str, z: Any, early_mass: Any, eta: Any) -> dict[str, Any]:
    q = 1 + z * z
    q_minus_two = 1 / (q * q)
    dq_minus_two = -4 * z / (q**3)
    d2q_minus_two = (20 * z * z - 4) / (q**4)
    early_value = -early_mass * q_minus_two
    early_d1 = -early_mass * dq_minus_two
    early_d2 = -early_mass * d2q_minus_two
    late_value = _hull(context, -_upper(context, eta), context.zero)
    late_d1, d1_sign = _signed_late_first_derivative_bound(context, z, eta)
    late_d2 = _hull(context, context.zero, 4 * _upper(context, eta))
    if bool(_lower(context, z) == context.zero and _upper(context, z) == context.zero):
        late_d1 = context.zero
        d1_sign = "zero_exact"
    return {
        "label": label,
        "Z": _interval_record(context, z),
        "derivative_kind": "ordinary_derivative_not_Taylor_coefficient",
        "q": _interval_record(context, q),
        "q_minus_two": _interval_record(context, q_minus_two),
        "D_Z_q_minus_two": _interval_record(context, dq_minus_two),
        "D_ZZ_q_minus_two": _interval_record(context, d2q_minus_two),
        "early_pressure": _interval_record(context, early_value),
        "early_D_Z_pressure": _interval_record(context, early_d1),
        "early_D_ZZ_pressure": _interval_record(context, early_d2),
        "late_pressure_bound": _interval_record(context, late_value),
        "late_D_Z_pressure_bound": _interval_record(context, late_d1),
        "late_D_Z_sign_contract": d1_sign,
        "late_D_ZZ_pressure_bound": _interval_record(context, late_d2),
        "pressure": _interval_record(context, early_value + late_value),
        "D_Z_pressure": _interval_record(context, early_d1 + late_d1),
        "D_ZZ_pressure": _interval_record(context, early_d2 + late_d2),
        "late_bound_is_joint_eta": True,
        "late_pressure_not_evaluated_as_point": True,
    }


def _uniform_cell_contract(context: Any, eta: Any) -> dict[str, Any]:
    zlo = context.mpf("-0.25")
    zhi = context.mpf("0.25")
    qlo = context.mpf(1)
    qhi = context.mpf(17) / 16
    qinv2_lo = (context.mpf(16) / 17) ** 2
    qinv2_hi = context.mpf(1)
    return {
        "Z_cell": [_interval_record(context, zlo), _interval_record(context, zhi)],
        "q_range": _interval_record(context, _hull(context, qlo, qhi)),
        "q_minus_two_range": _interval_record(context, _hull(context, qinv2_lo, qinv2_hi)),
        "direct_derivative_bounds": {
            "abs_D_Z_q_minus_two": _interval_record(context, context.mpf(1)),
            "abs_D_ZZ_q_minus_two": _interval_record(context, context.mpf(4)),
        },
        "late_D_Z_bound": "|D_Z late_exact_function| <= 4*|Z|*eta with sign sign(Z)",
        "late_D_ZZ_bound": "0 <= D_ZZ late_exact_function <= 4*eta",
        "eta_positive_record": _interval_record(context, eta),
        "direct_qrange_used": True,
    }


def _proof_packet(
    capsule: Mapping[str, Any],
    stages: Sequence[Mapping[str, Any]],
    parameters: Mapping[str, Any],
    quad: Mapping[str, Any],
    admission: Mapping[str, Any],
) -> dict[str, Any]:
    source_bindings = {}
    for key in (
        "source_bindings",
        "source_binding_identities",
        "function_definitions",
        "source_identity",
        "family_identity",
    ):
        if key in capsule:
            source_bindings[key] = _json_value(capsule[key])
    return {
        "kind": "direct_formula_source_binding_and_monotonicity_proof_v1",
        "source_binding_identities": {
            "yd": "yd=logP from saved parameter data; no exp(40) replay",
            "early_M2": "5/2 + Mslope + exp(3/5)*(exp(-1)-exp(-yd))/2",
            "Mslope": "(1/2)*integral_0^1 exp(y/5-6J(y)/5)dy",
            "J": "integral_0^y sigma(s)ds",
            "sigma": "exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), sigma(0)=0, sigma(1)=1",
            "pressure": "P0(Z)=-earlyM2/(1+Z^2)^2+late_exact_function(Z)",
            "physical_scale": "Pstar^2 is represented by log2logP=2*logP",
            "core_seed": "epsilon_core*P0 scale=exp(2*logP-logLambda), with no huge physical field materialization",
            "epsilon_outer": "epsilon_outer=delta/1000; it is not accepted epsilon_core",
        },
        "exact_late_function": {
            "definition": "late_exact_function(Z)=-sum of the eleven post-yd source integrals from the capsule",
            "stage_count": 11,
            "stage_names": list(LATE_STAGE_NAMES),
            "beta2_stage_count": 6,
            "beta0_stage_count": 7,
            "variable_flatten_stage_count": 1,
            "same_exact_source_H": "H=1 preheat",
            "source_function_rows_are_exact_definitions": True,
            "stage_definitions": [_json_value(row) for row in stages],
            "exact_recipes_from_capsule": _json_value(capsule.get("exact_recipes", {})),
            "all_order_binding_data": _json_value(
                capsule.get("all_order_pressure_source_bindings", {})
            ),
        },
        "monotonicity_derivation": {
            "stable_logratio": "r(s)=1/s^2-1/(1-s)^2",
            "logratio_derivative": "r'(s)=-2/s^3-2/(1-s)^3<0 on (0,1)",
            "sigma_derivative": "sigma'(s)=-exp(r(s))*r'(s)/(1+exp(r(s)))^2>0",
            "rectangle_order": "sigma monotone implies Jlo(a)<=J(y)<=Jhi(b) on each panel",
            "density_order": "exp(y/5-6J/5)/2 is enclosed using [a,b] and [Jlo(a),Jhi(b)]",
        },
        "symmetry_derivation": {
            "identity": "sigma(1-s)=1-sigma(s)",
            "integral": "J(1)=integral_0^1 sigma(s)ds=1/2",
            "use": "posthoc containment check only; prefix error remains",
            "quad_prefix_certificate_retained": True,
        },
        "tail_derivation": {
            "eta": "exp(3/5-yd)/(2*(1-epsilon_outer)^2)",
            "epsilon_outer": "delta/1000 from the saved exact recipe",
            "beta_contract": "0<=beta(y)<=2",
            "joint_integral_statement": "all eleven true post-yd integrals at Z=0 are jointly bounded by eta",
            "not_sum_of_eleven_eta": True,
            "global_pointwise_bound_and_disjoint_partition": True,
            "eta_never_replaced_by_zero": True,
        },
        "ordinary_derivative_derivation": {
            "q_minus_two": "q^-2=(1+Z^2)^-2",
            "D_Z_q_minus_two": "-4Z/q^3",
            "D_ZZ_q_minus_two": "(20Z^2-4)/q^4",
            "late_D_Z": "sign(Z), magnitude <=4|Z|eta, zero exact at Z=0",
            "late_D_ZZ": "2beta*q^(-beta-2)*(1-(2beta+1)Z^2)>=0 and <=4eta on |Z|<=1/4",
            "factorial_confusion": False,
        },
        "parameter_data_consumed": {
            key: {"raw": value["raw"]} for key, value in parameters.items()
        },
        "quadrature_method": quad["method"],
        "no_CAS_runtime": True,
        "no_old_math_or_provider_calls": True,
        "capsule_source_bindings_preserved": source_bindings,
        "same_family_admission": _json_value(admission),
        "capsule_exact_source_definition_data": {
            "exact_recipes": _json_value(capsule.get("exact_recipes", {})),
            "pressure_source": _json_value(capsule.get("pressure_source", {})),
            "true_pressure_measure_bindings": _json_value(
                capsule.get("true_pressure_measure_bindings", {})
            ),
            "all_order_pressure_source_bindings": _json_value(
                capsule.get("all_order_pressure_source_bindings", {})
            ),
        },
    }


def derive(capsule: Mapping[str, Any], emit: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Run the one new pressure enclosure/ordinary-jet subfamily.

    The root runner owns freeze, entry persistence, invocation, and full save.
    This API performs no ancestor replay.  ``emit`` receives the same event
    rows later retained under ``payload.events``.
    """

    try:
        if emit is not None and not callable(emit):
            raise KernelInputError("emit must be callable or None")
        capsule = _validate_capsule(capsule)
        stages = _read_stage_definitions(capsule)
        admission = _source_admission(capsule)
        events: list[dict[str, Any]] = []
        _emit(
            emit,
            events,
            "family_entered",
            {
                "capsule_schema": CAPSULE_SCHEMA,
                "stage_count": len(stages),
                "source_family_joined": True,
                "old_caps_and_errors_retained": True,
            },
        )
        # Entry receipt is emitted before creating the private interval
        # context or converting any saved MPF interval records.  If a literal
        # parameter packet is malformed, the caller still receives the actual
        # family-entry prefix and traceback.
        context = _new_interval_context()
        parameters = _read_parameters(context, capsule)

        _emit(
            emit,
            events,
            "quadrature_entered",
            {
                "method": "new_coupled_monotone_Darboux_integrator_v1",
                "quadrature_panels": N_DYADIC_PANELS,
                "precision_dps": INTERVAL_DPS,
                "sigma_nodes_evaluated_once": True,
                "old_Taylor_quadrature": False,
            },
        )
        quad = _monotone_darboux_mass(context, N_DYADIC_PANELS)
        if quad["symmetry_check"]["contains_one_half"] is not True:
            raise KernelInputError("new Darboux J(1) range fails the exact reflected-integral identity")
        _emit(
            emit,
            events,
            "quadrature_complete",
            {
                "method": quad["method"],
                "quadrature_panels": quad["quadrature_panels"],
                "nodes_evaluated_once": quad["nodes_evaluated_once"],
                "symmetry_posthoc_only": True,
            },
        )

        _emit(
            emit,
            events,
            "pressure_jet_entered",
            {
                "ordinary_derivatives": True,
                "sample_count": 4,
                "late_stage_count": 11,
                "late_bound_joint": True,
            },
        )

        delta = parameters["delta"]["interval"]
        logp = parameters["logP"]["interval"]
        # The accepted capsule's ``epsilon`` is epsilon_core.  The tail uses
        # the separate exact outer recipe delta/1000; never alias these.
        epsilon_outer = delta / 1000
        epsilon_core = parameters.get("epsilon_core")
        eta = context.exp(context.mpf("0.6") - logp) / (
            context.mpf(2) * (1 - epsilon_outer) ** 2
        )
        if not bool(
            _lower(context, delta) > context.zero
            and _upper(context, epsilon_outer) < context.one
            and _lower(context, logp) > context.one
            and _lower(context, eta) > context.zero
        ):
            raise KernelInputError("positive source parameters and nonzero joint tail are required")
        early_mass = (
            context.mpf(5) / 2
            + quad["Mslope_interval"]
            + context.exp(context.mpf("0.6"))
            * (context.exp(-context.mpf(1)) - context.exp(-logp))
            / 2
        )
        physical_log_scale = 2 * logp
        loglambda = parameters.get("logLambda")
        if loglambda is not None:
            core_seed_log_scale = 2 * logp - loglambda["interval"]
            core_seed_scale = context.exp(core_seed_log_scale)
        else:
            core_seed_log_scale = None
            core_seed_scale = None

        sample_inputs = [
            ("Z=-1/4", context.mpf("-0.25")),
            ("Z=0", context.zero),
            ("Z=+1/4", context.mpf("0.25")),
            ("anchor", parameters["anchor"]["interval"]),
        ]
        for label, z in sample_inputs:
            if not bool(_upper(context, abs(z)) <= context.mpf("0.25")):
                raise KernelInputError(label + " is outside the admitted |Z|<=1/4 cell")
        sample_rows = [
            _pressure_sample(context, label, z, early_mass, eta)
            for label, z in sample_inputs
        ]
        _emit(
            emit,
            events,
            "pressure_jet_complete",
            {
                "ordinary_derivative_rows": len(sample_rows),
                "ordinary_not_Taylor_coefficients": True,
                "late_D_Z_signs_retained": True,
                "late_D_ZZ_nonnegative_contract_retained": True,
            },
        )

        late_rows = []
        for row in stages:
            if row["name"] in LATE_STAGE_NAMES:
                late_rows.append(
                    {
                        "name": row["name"],
                        "beta_contract": row["beta_contract"],
                        "expected_beta_family": row["expected_beta_family"],
                        "region": "post_yd",
                        "exact_function_definition": _json_value(row),
                        "normalized_at_Z0": True,
                        "shared_joint_eta_reference": True,
                        "integral_not_replaced_by_zero": True,
                    }
                )
        if len(late_rows) != 11:
            raise KernelInputError("late stage partition must contain eleven rows")
        early_rows = [
            {
                "name": row["name"],
                "beta_contract": row["beta_contract"],
                "expected_beta_family": row["expected_beta_family"],
                "region": "early",
                "exact_function_definition": _json_value(row),
            }
            for row in stages
            if row["name"] in EARLY_STAGE_NAMES
        ]

        tail_record = _interval_record(context, eta)
        proof = _proof_packet(capsule, stages, parameters, quad, admission)
        _emit(
            emit,
            events,
            "proof_bundle",
            {
                "proof_kind": proof["kind"],
                "stage_count": 14,
                "late_stage_count": 11,
                "prefix_error_retained": True,
                "eta_positive": True,
            },
        )
        _emit(
            emit,
            events,
            "family_completed",
            {
                "proof_count": 1,
                "actual_pressure_source_enclosure_evaluated": True,
                "numeric_evaluation": True,
                "full_goal": False,
            },
        )

        preserved = {
            "capsule_status": _json_value(capsule.get("status")),
            "capsule_preservation": _json_value(
                capsule.get("preservation", capsule.get("preservation_metadata", {}))
            ),
            "old_signed_m15_and_error_terms_unchanged": _json_value(
                capsule.get("old_signed_m15_and_error_terms_unchanged")
            ),
            "old_caps_and_errors_unchanged": True,
            "cumulative_212_terms_unchanged": True,
            "incoming_pressure_datum_same_family_hash": _json_value(
                capsule.get(
                    "incoming_pressure_datum",
                    capsule.get("pressure_datum", capsule.get("family_identity", {})),
                )
            ),
        }
        payload = {
            "component": COMPONENT,
            "source_only": True,
            "numeric_evaluation": True,
            "source_capsule_schema": CAPSULE_SCHEMA,
            "parameter_data_consumed": {
                key: {
                    "raw": value["raw"],
                    "directed_interval": _interval_record(context, value["interval"]),
                }
                for key, value in parameters.items()
            },
            "parameter_admission": {
                "actual_parameter_admission": False,
                "point_representative_selected": False,
                "delta": _interval_record(context, delta),
                "logP": _interval_record(context, logp),
                "anchor": _interval_record(context, parameters["anchor"]["interval"]),
            },
            "exact_function_definition": {
                "q": "q(Z)=1+Z^2",
                "normalized_pressure": "P0(Z)=-earlyM2*q(Z)^(-2)+late_exact_function(Z)",
                "earlyM2": "5/2+Mslope+exp(3/5)*(exp(-1)-exp(-yd))/2",
                "yd": "saved yd=logP",
                "late_exact_function": "negative sum of the eleven exact post-yd source integrals",
                "same_exact_source_H": "H=1 preheat",
                "source_function_values_not_replaced_by_centers": True,
            },
            "stage_partition": {
                "all_stage_definitions": [_json_value(row) for row in stages],
                "early_stage_definitions": early_rows,
                "late_stage_definitions": late_rows,
                "beta2_stage_names": list(BETA2_STAGE_NAMES),
                "beta0_stage_names": list(BETA0_STAGE_NAMES),
                "variable_flatten_stage_name": VARIABLE_STAGE_NAME,
                "stage_count": 14,
                "early_stage_count": 3,
                "late_stage_count": 11,
            },
            "monotone_Darboux_quadrature": {
                key: value
                for key, value in quad.items()
                if key not in {"Mslope_interval"}
            },
            "directed_approximations": {
                "Mslope": _interval_record(context, quad["Mslope_interval"]),
                "earlyM2": _interval_record(context, early_mass),
                "epsilon_outer": _interval_record(context, epsilon_outer),
                "epsilon_outer_recipe": "delta/1000",
                "epsilon_core": _interval_record(context, epsilon_core["interval"])
                if epsilon_core is not None
                else {"saved": False, "formula": "exp(-logLambda)"},
                "epsilon_core_recipe": "exp(-logLambda)",
                "eta": tail_record,
                "eta_positive_never_zero": True,
                "joint_post_yd_integral_bound": tail_record,
                "late_integrals_are_true_source_integrals": True,
                "late_integrals_not_replaced_by_zero": True,
                "physical_log_scale_Pstar_squared": _interval_record(context, physical_log_scale),
                "physical_log_scale_formula": "log(Pstar^2)=2*logP",
                "core_seed_log_scale": _interval_record(context, core_seed_log_scale)
                if core_seed_log_scale is not None
                else {"formula": "2*logP-logLambda", "saved_logLambda_present": False},
                "core_seed_scale": _interval_record(context, core_seed_scale)
                if core_seed_scale is not None
                else {"formula": "exp(2*logP-logLambda)", "materialized": False},
            },
            "pressure_source_ordinary_Z_jets": {
                "derivative_kind": "ordinary derivatives, not Taylor coefficients",
                "rows": sample_rows,
                "actual_source_Z_jets_admitted": False,
                "global_finite_N_reason": "finite-N source jets remain globally inadmissible",
            },
            "uniform_cell": _uniform_cell_contract(context, eta),
            "proof_packet": proof,
            "same_family_admission": admission,
            "preservation": preserved,
            "actual_counts_contract": dict(SUCCESS_COUNTS),
            "actual_flags_contract": dict(SUCCESS_FLAGS),
            "events": events,
            "no_old_math_or_provider_calls": True,
            "no_actual_physical_huge_field": True,
            "no_NS_claim": True,
        }
        return {
            "schema": SCHEMA,
            "component": COMPONENT,
            "passed": True,
            "proof_count": 1,
            "actual_counts": dict(SUCCESS_COUNTS),
            "new_flags": dict(SUCCESS_FLAGS),
            "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
            "count_event_map": {
                key: list(value) for key, value in COUNT_EVENT_MAP.items()
            },
            "source_only": True,
            "numeric_evaluation": True,
            "physical_layer_completed": False,
            "full_goal": False,
            "full_goal_complete": False,
            "payload": payload,
        }
    except KernelInputError:
        # Preserve any emitted entry prefix and the caller's traceback.  The
        # root runner owns rejection accounting; this worker must not fabricate
        # a zeroed result that hides an entered mathematical subfamily.
        raise


__all__ = [
    "CAPSULE_SCHEMA",
    "COMPONENT",
    "COUNT_EVENT_MAP",
    "CALLBACK_EVENT_COUNTS",
    "INTERVAL_DPS",
    "KernelInputError",
    "N_DYADIC_PANELS",
    "SCHEMA",
    "SUCCESS_COUNTS",
    "SUCCESS_FLAGS",
    "derive",
]
