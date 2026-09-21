"""Current-lineage nonlinear mean attribution through the first public RF40 turn.

This Agent-3 increment consumes the exact Agent-2 PR #992 field

    u_current_to_X1 = u_lead,current-through-X1 + u_osc,frozen-complete-curl

with X_1 = e X_R, together with the already-frozen Agent-2 PR #960
oscillatory Jacobian.  It deliberately reuses the existing Agent-3 nonlinear
decomposition and rotating cylindrical m=0 projector rather than reproducing
Agent-2 curl/Jacobian machinery or Agent-1 RF40 leading construction.

This remains a scoped nonlinear-transport diagnostic.  It does not materialize
pressure, restricted forcing, a complete Navier--Stokes defect, a correction
velocity, held-out residual reduction, or a blow-up proof.
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

TASK = "KOKUNO-A3-CURRENT-RF40-FIRST-TURN-NONLINEAR-MEAN-106"
SCHEMA = "kokuno-a3-current-rf40-first-turn-nonlinear-mean-v1"

PARENT_AGENT3_PR = 988
PARENT_AGENT3_HEAD = "fa95dee3709a82326167af1a2cb42b9b5f9a88a5"

AGENT2_COMPOSITE_PR = 992
AGENT2_COMPOSITE_HEAD = "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377"
AGENT2_COMPOSITE_MODULE = (
    "openai_ns_reconstruction.kokuno_current_rf40_first_turn_leading_oscillatory_identity"
)
AGENT2_COMPOSITE_CLASS = "CurrentRF40FirstTurnLeadingOscillatoryField"
AGENT2_COMPOSITE_SOURCE_BLOB = "c126dace161e5dae1745185c79152736817994b8"

AGENT2_COMPOSITE_PARENT_PR = 987
AGENT2_COMPOSITE_PARENT_HEAD = "4be2c9ee898c44dd1ad2217e90601b161fe81964"

AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_MODULE = (
    "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
)
AGENT2_DIFFERENTIAL_FUNCTION = "evaluate_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

AGENT1_LEADING_PR = 986
AGENT1_LEADING_HEAD = "23c00b98526e187ff04d432745300a084ec859f2"
AGENT1_LEADING_SOURCE_BLOB = "bf71948cbdb827d791b86ab8f00496bd1b231961"

SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


def _validate_first_turn_configuration(
    configuration: Mapping[str, Any],
    *,
    runtime_semantic_sha256: str,
    runtime_oscillatory_sha256: str,
) -> None:
    """Validate the exact #992 identity/truth contract without recomputing A2 work."""
    if not isinstance(configuration, Mapping):
        raise ValueError("A2 #992 configuration must be a mapping")

    parent = configuration.get("parent_agent2")
    if not isinstance(parent, Mapping):
        raise ValueError("A2 #992 parent identity is missing")
    if parent.get("pr") != AGENT2_COMPOSITE_PARENT_PR:
        raise ValueError("A2 #992 parent PR identity drifted")
    if parent.get("head") != AGENT2_COMPOSITE_PARENT_HEAD:
        raise ValueError("A2 #992 parent head identity drifted")

    a1 = configuration.get("agent1_leading")
    if not isinstance(a1, Mapping):
        raise ValueError("A2 #992 Agent-1 leading identity is missing")
    if a1.get("pr") != AGENT1_LEADING_PR or a1.get("head") != AGENT1_LEADING_HEAD:
        raise ValueError("A2 #992 Agent-1 leading lineage drifted")
    if a1.get("source_blob_sha1") != AGENT1_LEADING_SOURCE_BLOB:
        raise ValueError("A2 #992 Agent-1 leading source identity drifted")
    if not _valid_sha256(a1.get("semantic_sha256")):
        raise ValueError("A2 #992 Agent-1 semantic identity is malformed")
    if not _valid_sha256(a1.get("configuration_sha256")):
        raise ValueError("A2 #992 Agent-1 configuration identity is malformed")

    oscillatory = configuration.get("oscillatory_runtime")
    if not isinstance(oscillatory, Mapping):
        raise ValueError("A2 #992 oscillatory runtime identity is missing")
    payload_sha = oscillatory.get("payload_sha256")
    if not _valid_sha256(payload_sha):
        raise ValueError("A2 #992 oscillatory runtime SHA is malformed")
    if payload_sha != runtime_oscillatory_sha256:
        raise ValueError("A2 #992 oscillatory runtime SHA disagrees with runtime")

    if not _valid_sha256(runtime_semantic_sha256):
        raise ValueError("A2 #992 composite semantic identity is malformed")

    truth = configuration.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("A2 #992 truth boundary is missing")
    required_true = (
        "current_leading_plus_oscillatory_velocity_through_RF40_first_turn_materialized",
        "post_XR_RF40_first_turn_leading_plus_oscillatory_materialized",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_RF40_first_turn_composite_save_load_available",
        "reload_recomputes_leading_and_oscillatory_identities",
    )
    required_false = (
        "velocity_after_RF40_first_turn_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "RF40_axial_shutdown_current_lineage_materialized",
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
            raise ValueError(f"A2 #992 truth boundary drifted at {key}")
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"A2 #992 truth boundary drifted at {key}")


@dataclass(frozen=True)
class ExactCurrentRF40FirstTurnNonlinearBackend:
    """Checksum-bound A2 #992/#960 bridge into the existing A3 mean projector."""

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
    ) -> "ExactCurrentRF40FirstTurnNonlinearBackend":
        field_type = type(composite_field)
        if field_type.__module__ != AGENT2_COMPOSITE_MODULE:
            raise ValueError("unexpected A2 #992 composite module identity")
        if field_type.__name__ != AGENT2_COMPOSITE_CLASS:
            raise ValueError("unexpected A2 #992 composite class identity")
        if not callable(getattr(composite_field, "velocity", None)):
            raise TypeError("A2 #992 composite must expose velocity(x,y,z,t)")
        if not callable(getattr(composite_field, "configuration", None)):
            raise TypeError("A2 #992 composite must expose configuration()")

        leading = getattr(composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise TypeError("A2 #992 composite must retain its authenticated leading backend")

        field_source = inspect.getsourcefile(field_type)
        if field_source is None:
            raise ValueError("A2 #992 composite source is unavailable")
        field_blob = _git_blob_sha1(field_source)
        if field_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise ValueError("A2 #992 composite source blob drifted")

        composite_semantic = str(getattr(composite_field, "semantic_sha256", ""))
        oscillatory_semantic = str(
            getattr(composite_field, "oscillatory_runtime_sha256", "")
        )
        configuration = composite_field.configuration()
        _validate_first_turn_configuration(
            configuration,
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
        }


def materialize_current_rf40_first_turn_nonlinear_mean_attribution(
    backend: ExactCurrentRF40FirstTurnNonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project the exact through-X1 nonlinear pieces; no surrogate defect is accepted."""
    if not isinstance(backend, ExactCurrentRF40FirstTurnNonlinearBackend):
        raise TypeError("backend must be ExactCurrentRF40FirstTurnNonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ValueError("A2 #992 composite provenance drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ValueError("A2 #960 differential provenance drifted")
    return _materialize_from_provider(
        backend.components,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact-current-a2-pr992-through-rf40-first-turn-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_current_rf40_first_turn_nonlinear_mean_attribution
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
        "source_reader_date": SOURCE_READER_DATE,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "current_leading_plus_oscillation_through_RF40_first_turn_consumed": True,
        "post_XR_RF40_first_turn_velocity_materialized": True,
        "current_RF40_first_turn_mixed_nonlinear_mean_materialized": True,
        "current_RF40_first_turn_quadratic_nonlinear_mean_materialized": True,
        "current_RF40_first_turn_aggregate_nonlinear_mean_materialized": True,
        "current_leading_spatial_derivative_is_fixed_repository_fd": True,
        "caller_tunable_leading_spatial_step": False,
        "agent2_oscillatory_jacobian_reimplemented_by_agent3": False,
        "velocity_after_RF40_first_turn_materialized": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "RF40_axial_shutdown_current_lineage_materialized": False,
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
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "paper_exact": False,
        "pde_validated": False,
        "blowup_proved": False,
    }
