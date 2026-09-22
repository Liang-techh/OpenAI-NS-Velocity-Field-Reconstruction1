"""Exact external-runtime rebind contract for the current Kokuno I4 composite.

Kokuno Agent 2 owns only the dependency/runtime identity boundary here.  The
current branch already contains the exact A2 #960 oscillatory differential
runtime and the A3 exact-backend bind implementation.  The exact A2 #1080
current-I4 composite and its exact A1 #1079 leading backend remain pinned
external repository snapshots; they are not copied, modified, or relabelled as
source-exact data.

This module accepts those already-constructed exact runtime objects and asks
``ExactCurrentI4NonlinearBackend.bind`` to authenticate them.  It does not run
RF30--RF39, apply a mean correction, form Cartesian delta-u, or assess an NS
residual.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
from typing import Any, Callable

from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT1_LEADING_CLASS,
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_MODULE,
    AGENT1_LEADING_SOURCE_BLOB,
    AGENT2_COMPOSITE_CLASS,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_MODULE,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_FUNCTION,
    AGENT2_DIFFERENTIAL_HEAD,
    AGENT2_DIFFERENTIAL_MODULE,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    ExactCurrentI4NonlinearBackend,
    _git_blob_sha1,
)

TASK = "K2-OSC-116"
SCHEMA = "kokuno-a2-current-i4-exact-external-runtime-rebind-v1"

PARENT_A2_PR = 1234
PARENT_A2_HEAD = "f90bdb9c232fb16df881dfbcac28be1cd5c29780"
A2_COMPOSITE_PR = 1080
A1_LEADING_PR = 1079
A2_DIFFERENTIAL_PR = 960

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class ExactRuntimeRebindError(RuntimeError):
    """Raised when the pinned runtime identity cannot be authenticated."""


def _sha256(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _source_blob(obj: Any, label: str) -> str:
    source = inspect.getsourcefile(obj)
    if source is None:
        raise ExactRuntimeRebindError(f"{label} source path is unavailable")
    return _git_blob_sha1(source)


@dataclass(frozen=True)
class ExactCurrentI4RuntimeRebindReceipt:
    composite_module: str
    composite_class: str
    composite_source_blob: str
    composite_semantic_sha256: str
    oscillatory_runtime_sha256: str
    leading_module: str
    leading_class: str
    leading_source_blob: str
    differential_module: str
    differential_function: str
    differential_source_blob: str
    differential_semantic_sha256: str
    exact_backend_rebind_executed: bool
    exact_external_a2_1080_runtime_authenticated: bool
    exact_external_a1_1079_runtime_authenticated: bool
    local_exact_a2_960_runtime_authenticated: bool
    rf30_rf39_correction_executed_here: bool
    rf44_postupdate_state_materialized: bool
    cartesian_delta_u_materialized: bool
    heldout_ns_residual_assessed: bool
    residual_reduction_claimed: bool
    paper_exact: bool
    pde_validated: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "schema": SCHEMA,
            "task": TASK,
            "parent_a2_pr": PARENT_A2_PR,
            "parent_a2_head": PARENT_A2_HEAD,
            "a2_composite_pr": A2_COMPOSITE_PR,
            "a2_composite_head": AGENT2_COMPOSITE_HEAD,
            "a1_leading_pr": A1_LEADING_PR,
            "a1_leading_head": AGENT1_LEADING_HEAD,
            "a2_differential_pr": A2_DIFFERENTIAL_PR,
            "a2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
            "truth_boundary": truth_boundary(),
        }


def bind_exact_current_i4_runtime(
    composite_field: object,
    differential_function: Callable[[Any, Any], Any],
) -> tuple[ExactCurrentI4NonlinearBackend, ExactCurrentI4RuntimeRebindReceipt]:
    """Authenticate the pinned #1080/#1079 runtime against local exact #960.

    No correction or residual routine is invoked by this function.
    """
    try:
        backend = ExactCurrentI4NonlinearBackend.bind(
            composite_field,
            differential_function,
        )
    except Exception as exc:
        raise ExactRuntimeRebindError(
            "exact A2 #1080 / A1 #1079 / A2 #960 runtime rebind failed"
        ) from exc

    field_type = type(composite_field)
    leading = getattr(composite_field, "leading_backend", None)
    if leading is None:
        raise ExactRuntimeRebindError("authenticated composite lost its leading backend")
    leading_type = type(leading)

    checks = (
        (field_type.__module__, AGENT2_COMPOSITE_MODULE, "composite module"),
        (field_type.__name__, AGENT2_COMPOSITE_CLASS, "composite class"),
        (leading_type.__module__, AGENT1_LEADING_MODULE, "leading module"),
        (leading_type.__name__, AGENT1_LEADING_CLASS, "leading class"),
        (
            getattr(differential_function, "__module__", None),
            AGENT2_DIFFERENTIAL_MODULE,
            "differential module",
        ),
        (
            getattr(differential_function, "__name__", None),
            AGENT2_DIFFERENTIAL_FUNCTION,
            "differential function",
        ),
    )
    for actual, expected, label in checks:
        if actual != expected:
            raise ExactRuntimeRebindError(f"{label} drifted")

    composite_blob = _source_blob(field_type, "A2 #1080 composite")
    leading_blob = _source_blob(leading_type, "A1 #1079 leading")
    differential_blob = _source_blob(differential_function, "A2 #960 differential")
    if composite_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ExactRuntimeRebindError("A2 #1080 composite source blob drifted")
    if leading_blob != AGENT1_LEADING_SOURCE_BLOB:
        raise ExactRuntimeRebindError("A1 #1079 leading source blob drifted")
    if differential_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ExactRuntimeRebindError("A2 #960 differential source blob drifted")

    backend_receipt = backend.to_receipt()
    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "composite_source_blob": composite_blob,
        "composite_semantic_sha256": backend.composite_semantic_sha256,
        "oscillatory_runtime_sha256": backend.oscillatory_runtime_sha256,
        "leading_source_blob": leading_blob,
        "differential_source_blob": differential_blob,
        "differential_semantic_sha256": backend.differential_semantic_sha256,
        "backend_receipt": backend_receipt,
        "exact_backend_rebind_executed": True,
        "rf30_rf39_correction_executed_here": False,
        "heldout_ns_residual_assessed": False,
    }
    receipt_sha = _sha256(payload)

    receipt = ExactCurrentI4RuntimeRebindReceipt(
        composite_module=field_type.__module__,
        composite_class=field_type.__name__,
        composite_source_blob=composite_blob,
        composite_semantic_sha256=backend.composite_semantic_sha256,
        oscillatory_runtime_sha256=backend.oscillatory_runtime_sha256,
        leading_module=leading_type.__module__,
        leading_class=leading_type.__name__,
        leading_source_blob=leading_blob,
        differential_module=str(differential_function.__module__),
        differential_function=str(differential_function.__name__),
        differential_source_blob=differential_blob,
        differential_semantic_sha256=backend.differential_semantic_sha256,
        exact_backend_rebind_executed=True,
        exact_external_a2_1080_runtime_authenticated=True,
        exact_external_a1_1079_runtime_authenticated=True,
        local_exact_a2_960_runtime_authenticated=True,
        rf30_rf39_correction_executed_here=False,
        rf44_postupdate_state_materialized=False,
        cartesian_delta_u_materialized=False,
        heldout_ns_residual_assessed=False,
        residual_reduction_claimed=False,
        paper_exact=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )
    enforce_receipt(receipt.to_dict())
    return backend, receipt


def enforce_receipt(receipt: dict[str, object]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task") != TASK:
        raise AssertionError("runtime rebind receipt identity drifted")
    expected_true = (
        "exact_backend_rebind_executed",
        "exact_external_a2_1080_runtime_authenticated",
        "exact_external_a1_1079_runtime_authenticated",
        "local_exact_a2_960_runtime_authenticated",
    )
    for key in expected_true:
        if receipt.get(key) is not True:
            raise AssertionError(f"required runtime authentication missing at {key}")
    expected_false = (
        "rf30_rf39_correction_executed_here",
        "rf44_postupdate_state_materialized",
        "cartesian_delta_u_materialized",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    )
    for key in expected_false:
        if receipt.get(key) is not False:
            raise AssertionError(f"truth boundary promoted unexpectedly at {key}")
    if receipt.get("composite_source_blob") != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise AssertionError("composite source blob mismatch")
    if receipt.get("leading_source_blob") != AGENT1_LEADING_SOURCE_BLOB:
        raise AssertionError("leading source blob mismatch")
    if receipt.get("differential_source_blob") != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise AssertionError("differential source blob mismatch")


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(bind_exact_current_i4_runtime)
    forbidden = {
        "residual", "defect", "forcing", "pressure", "gain", "optimizer",
        "threshold", "viscosity", "nu", "coefficient", "correction", "target",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "a2_1080_runtime_is_pinned_external_repository_snapshot": True,
        "a1_1079_runtime_is_pinned_external_repository_snapshot": True,
        "a2_960_differential_runtime_is_local_exact_blob": True,
        "exact_backend_bind_is_the_authenticator": True,
        "external_runtime_bytes_copied_or_modified_by_this_increment": False,
        "oscillatory_complete_curl_reimplemented": False,
        "leading_profile_reimplemented": False,
        "agent3_mean_correction_reimplemented": False,
        "rf30_rf39_correction_executed_here": False,
        "rf44_postupdate_state_materialized": False,
        "cartesian_delta_u_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
