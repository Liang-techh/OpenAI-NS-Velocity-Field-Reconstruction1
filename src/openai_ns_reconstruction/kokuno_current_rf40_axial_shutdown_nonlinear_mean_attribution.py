"""Current-lineage nonlinear mean attribution through RF40 axial shutdown.

This Agent-3 increment consumes the exact Agent-2 PR #999 field

    u_current_to_X2 = u_lead,current-through-X2 + u_osc,frozen-complete-curl

through the RF40 axial-shutdown endpoint X_2, together with the already-frozen
Agent-2 PR #960 oscillatory Jacobian.  It deliberately reuses the existing
Agent-3 nonlinear decomposition and rotating cylindrical m=0 projector rather
than reproducing Agent-2 curl/Jacobian machinery or Agent-1 RF40 construction.

This remains a scoped nonlinear-transport attribution.  It does not perform the
radial inverse, materialize pressure or restricted forcing, define a complete
Navier--Stokes defect, construct a correction velocity, run a finite correction
cycle, or assess held-out full-PDE residual reduction.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_current_partial_nonlinear_mean_attribution import (
    FIXED_LEADING_SPATIAL_STEP,
    _components_from_interfaces,
    _git_blob_sha1,
    _valid_sha256,
)
from .kokuno_oscillatory_nonlinear_mean_attribution import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    OscillatoryNonlinearMeanAttributionWitness,
    _materialize_from_provider,
)

TASK = "KOKUNO-A3-CURRENT-RF40-AXIAL-SHUTDOWN-NONLINEAR-MEAN-109"
SCHEMA = "kokuno-a3-current-rf40-axial-shutdown-nonlinear-mean-v1"

PARENT_AGENT3_PR = 1006
PARENT_AGENT3_HEAD = "2c10faab7c229b801cb868355e0a1de0a729445c"

AGENT2_COMPOSITE_PR = 999
AGENT2_COMPOSITE_HEAD = "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5"
AGENT2_COMPOSITE_MODULE = (
    "openai_ns_reconstruction.kokuno_current_rf40_axial_shutdown_leading_oscillatory_identity"
)
AGENT2_COMPOSITE_CLASS = "CurrentRF40AxialShutdownLeadingOscillatoryField"
AGENT2_COMPOSITE_SOURCE_BLOB = "255e65ce0c37652858509c33e7c9fad40b73ca97"

AGENT2_COMPOSITE_PARENT_PR = 992
AGENT2_COMPOSITE_PARENT_HEAD = "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377"
AGENT2_COMPOSITE_PARENT_SOURCE_BLOB = "c126dace161e5dae1745185c79152736817994b8"

AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_MODULE = (
    "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
)
AGENT2_DIFFERENTIAL_FUNCTION = "evaluate_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

AGENT1_LEADING_PR = 993
AGENT1_LEADING_HEAD = "2ac6460b483efb1c07f2fa65e7fed781a32f2718"
AGENT1_LEADING_SOURCE_BLOB = "057a0514c8a941c3e59922b8158f480434b4441e"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_ZENODO_RECORD = "22678406"


def _expected_source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_READER_HEAD,
        "path": SOURCE_PATH,
        "blob_sha1": SOURCE_BLOB_SHA1,
        "corrected_release_date": SOURCE_READER_DATE,
        "zenodo_record": SOURCE_ZENODO_RECORD,
        "license": None,
        "classification": "reimplement_math_only",
        "migration_scope": "provenance_only",
        "scope": "localized waves / complete curls plus RF40 axial-shutdown structure",
        "paper_exact_claim": False,
    }


def _validate_axial_shutdown_configuration(
    configuration: Mapping[str, Any],
    *,
    runtime_semantic_sha256: str,
    runtime_oscillatory_sha256: str,
) -> None:
    """Validate the exact #999 identity/truth contract without recomputing A2 work."""
    if not isinstance(configuration, Mapping):
        raise ValueError("A2 #999 configuration must be a mapping")

    parent = configuration.get("parent_agent2")
    if not isinstance(parent, Mapping):
        raise ValueError("A2 #999 parent identity is missing")
    if parent.get("pr") != AGENT2_COMPOSITE_PARENT_PR:
        raise ValueError("A2 #999 parent PR identity drifted")
    if parent.get("head") != AGENT2_COMPOSITE_PARENT_HEAD:
        raise ValueError("A2 #999 parent head identity drifted")
    if parent.get("source_blob_sha1") != AGENT2_COMPOSITE_PARENT_SOURCE_BLOB:
        raise ValueError("A2 #999 parent source identity drifted")

    a1 = configuration.get("agent1_leading")
    if not isinstance(a1, Mapping):
        raise ValueError("A2 #999 Agent-1 leading identity is missing")
    if a1.get("pr") != AGENT1_LEADING_PR or a1.get("head") != AGENT1_LEADING_HEAD:
        raise ValueError("A2 #999 Agent-1 leading lineage drifted")
    if a1.get("source_blob_sha1") != AGENT1_LEADING_SOURCE_BLOB:
        raise ValueError("A2 #999 Agent-1 leading source identity drifted")
    if not _valid_sha256(a1.get("semantic_sha256")):
        raise ValueError("A2 #999 Agent-1 semantic identity is malformed")
    if not _valid_sha256(a1.get("configuration_sha256")):
        raise ValueError("A2 #999 Agent-1 configuration identity is malformed")

    oscillatory = configuration.get("oscillatory_runtime")
    if not isinstance(oscillatory, Mapping):
        raise ValueError("A2 #999 oscillatory runtime identity is missing")
    payload_sha = oscillatory.get("payload_sha256")
    if not _valid_sha256(payload_sha):
        raise ValueError("A2 #999 oscillatory runtime SHA is malformed")
    if payload_sha != runtime_oscillatory_sha256:
        raise ValueError("A2 #999 oscillatory runtime SHA disagrees with runtime")

    if not _valid_sha256(runtime_semantic_sha256):
        raise ValueError("A2 #999 composite semantic identity is malformed")

    if configuration.get("source_provenance") != _expected_source_provenance():
        raise ValueError("A2 #999 corrected-source provenance drifted")

    truth = configuration.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("A2 #999 truth boundary is missing")
    required_true = (
        "current_leading_plus_oscillatory_velocity_through_RF40_axial_shutdown_materialized",
        "RF40_axial_shutdown_leading_plus_oscillatory_materialized",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_RF40_axial_shutdown_composite_save_load_available",
        "reload_recomputes_leading_and_oscillatory_identities",
        "a2_semantic_identity_binds_corrected_source_provenance_block",
    )
    required_false = (
        "agent1_internal_external_source_semantic_binding_repaired_by_this_increment",
        "velocity_after_RF40_axial_shutdown_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "RF40_lambda_turn_current_lineage_materialized",
        "RF40_power_law_current_lineage_materialized",
        "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "agent3_correction_velocity_composed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_velocity_pressure_forcing_api",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "residual_reduction_claimed",
        "velocity_export_ready",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"A2 #999 truth boundary drifted at {key}")
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"A2 #999 truth boundary drifted at {key}")


@dataclass(frozen=True)
class ExactCurrentRF40AxialShutdownNonlinearBackend:
    """Checksum-bound A2 #999/#960 bridge into the existing A3 mean projector."""

    composite_field: object
    differential_function: Callable[[Any, Any], Any]
    composite_semantic_sha256: str
    oscillatory_runtime_sha256: str
    differential_semantic_sha256: str
    composite_source_blob: str
    differential_source_blob: str

    @classmethod
    def bind(
        cls,
        composite_field: object,
        differential_function: Callable[[Any, Any], Any],
    ) -> "ExactCurrentRF40AxialShutdownNonlinearBackend":
        field_type = type(composite_field)
        if field_type.__module__ != AGENT2_COMPOSITE_MODULE:
            raise ValueError("unexpected A2 #999 composite module identity")
        if field_type.__name__ != AGENT2_COMPOSITE_CLASS:
            raise ValueError("unexpected A2 #999 composite class identity")
        if not callable(getattr(composite_field, "velocity", None)):
            raise TypeError("A2 #999 composite must expose velocity(x,y,z,t)")
        if not callable(getattr(composite_field, "configuration", None)):
            raise TypeError("A2 #999 composite must expose configuration()")

        leading = getattr(composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise TypeError("A2 #999 composite must retain its authenticated leading backend")

        field_source = inspect.getsourcefile(field_type)
        if field_source is None:
            raise ValueError("A2 #999 composite source is unavailable")
        field_blob = _git_blob_sha1(field_source)
        if field_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise ValueError("A2 #999 composite source blob drifted")

        composite_semantic = str(getattr(composite_field, "semantic_sha256", ""))
        oscillatory_semantic = str(
            getattr(composite_field, "oscillatory_runtime_sha256", "")
        )
        _validate_axial_shutdown_configuration(
            composite_field.configuration(),
            runtime_semantic_sha256=composite_semantic,
            runtime_oscillatory_sha256=oscillatory_semantic,
        )

        if not callable(differential_function):
            raise TypeError("A2 #960 differential function must be callable")
        if differential_function.__module__ != AGENT2_DIFFERENTIAL_MODULE:
            raise ValueError("unexpected A2 #960 differential module identity")
        if differential_function.__name__ != AGENT2_DIFFERENTIAL_FUNCTION:
            raise ValueError("unexpected A2 #960 differential function identity")
        diff_source = inspect.getsourcefile(differential_function)
        if diff_source is None:
            raise ValueError("A2 #960 differential source is unavailable")
        diff_blob = _git_blob_sha1(diff_source)
        if diff_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
            raise ValueError("A2 #960 differential source blob drifted")

        diff_module = inspect.getmodule(differential_function)
        diff_semantic_fn = getattr(diff_module, "batch_differentials_sha256", None)
        if not callable(diff_semantic_fn):
            raise ValueError("A2 #960 differential semantic identity is unavailable")
        diff_semantic = str(diff_semantic_fn())
        if not _valid_sha256(diff_semantic):
            raise ValueError("A2 #960 differential semantic identity is malformed")

        return cls(
            composite_field=composite_field,
            differential_function=differential_function,
            composite_semantic_sha256=composite_semantic,
            oscillatory_runtime_sha256=oscillatory_semantic,
            differential_semantic_sha256=diff_semantic,
            composite_source_blob=field_blob,
            differential_source_blob=diff_blob,
        )

    def components(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return _components_from_interfaces(
            self.composite_field.leading_backend.velocity,
            self.composite_field.velocity,
            self.differential_function,
            x,
            y,
            z,
            t,
        )

    def to_receipt(self) -> dict[str, object]:
        return {
            "agent2_composite_pr": AGENT2_COMPOSITE_PR,
            "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
            "agent2_composite_source_blob": self.composite_source_blob,
            "agent2_composite_semantic_sha256": self.composite_semantic_sha256,
            "agent2_oscillatory_runtime_sha256": self.oscillatory_runtime_sha256,
            "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
            "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
            "agent2_differential_source_blob": self.differential_source_blob,
            "agent2_differential_semantic_sha256": self.differential_semantic_sha256,
            "agent1_leading_pr": AGENT1_LEADING_PR,
            "agent1_leading_head": AGENT1_LEADING_HEAD,
            "fixed_leading_spatial_step": FIXED_LEADING_SPATIAL_STEP,
            "leading_derivative_realization": "centered_cartesian_fd2_repository_fixed",
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_blob": SOURCE_BLOB_SHA1,
        }


def materialize_current_rf40_axial_shutdown_nonlinear_mean_attribution(
    backend: ExactCurrentRF40AxialShutdownNonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project exact through-X2 nonlinear pieces; no surrogate defect is accepted."""
    if not isinstance(backend, ExactCurrentRF40AxialShutdownNonlinearBackend):
        raise TypeError("backend must be ExactCurrentRF40AxialShutdownNonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ValueError("A2 #999 composite provenance drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ValueError("A2 #960 differential provenance drifted")
    return _materialize_from_provider(
        backend.components,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact-current-a2-pr999-through-rf40-axial-shutdown-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_current_rf40_axial_shutdown_nonlinear_mean_attribution
    )
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "angular_order",
        "spatial_step", "derivative_step", "time_step", "viscosity", "nu",
        "normalized_score", "scientific_threshold", "correction", "stage_budget",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_blob": SOURCE_BLOB_SHA1,
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "current_leading_plus_oscillation_through_RF40_axial_shutdown_consumed": True,
        "RF40_axial_shutdown_velocity_materialized": True,
        "current_RF40_axial_shutdown_mixed_nonlinear_mean_materialized": True,
        "current_RF40_axial_shutdown_quadratic_nonlinear_mean_materialized": True,
        "current_RF40_axial_shutdown_aggregate_nonlinear_mean_materialized": True,
        "current_leading_spatial_derivative_is_fixed_repository_fd": True,
        "caller_tunable_leading_spatial_step": False,
        "agent2_oscillatory_jacobian_reimplemented_by_agent3": False,
        "velocity_after_RF40_axial_shutdown_materialized": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "RF40_lambda_turn_current_lineage_materialized": False,
        "RF40_power_law_current_lineage_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "radial_inverse_performed_in_this_increment": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "complete_ns_defect": False,
        "scoped_nonlinear_mean_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "caller_supplied_surrogate_defect_allowed": False,
        "public_inputs": list(signature.parameters),
        "forbidden_public_parameters_present": sorted(
            forbidden.intersection(signature.parameters)
        ),
        "forbidden_public_parameters_absent": not bool(
            forbidden.intersection(signature.parameters)
        ),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "final_normalized_momentum_target_met": False,
        "final_normalized_divergence_target_met": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
        "blowup_proved": False,
    }
