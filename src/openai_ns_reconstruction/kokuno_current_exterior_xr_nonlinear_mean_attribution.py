"""Current-lineage nonlinear mean attribution on the A2 through-X_R composite.

This Agent-3 increment consumes the exact Agent-2 PR #987 field

    u_current_to_XR = u_lead,current-through-XR + u_osc,frozen-complete-curl

together with the already-frozen Agent-2 PR #960 oscillatory Jacobian.  It
deliberately reuses the existing Agent-3 nonlinear decomposition and rotating
cylindrical m=0 projector rather than reproducing Agent-2's curl machinery or
Agent-1's leading profile.

The only new responsibility is a provenance-locked adapter from the current
through-X_R composite into the already-reviewed Agent-3 mean machinery.  This is
still a scoped nonlinear-transport diagnostic: pressure, restricted forcing,
the complete Navier--Stokes defect, a correction velocity, and held-out residual
reduction are not materialized here.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_current_partial_nonlinear_mean_attribution import (
    COMPOSITION_CLOSURE_GATE,
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

TASK = "KOKUNO-A3-CURRENT-EXTERIOR-XR-NONLINEAR-MEAN-105"
SCHEMA = "kokuno-a3-current-exterior-xr-nonlinear-mean-v1"

PARENT_AGENT3_PR = 982
PARENT_AGENT3_HEAD = "dd11348e8eecf298656d7da58aee4526c6102306"

AGENT2_COMPOSITE_PR = 987
AGENT2_COMPOSITE_HEAD = "4be2c9ee898c44dd1ad2217e90601b161fe81964"
AGENT2_COMPOSITE_MODULE = (
    "openai_ns_reconstruction.kokuno_current_exterior_xr_leading_oscillatory_identity"
)
AGENT2_COMPOSITE_CLASS = "CurrentExteriorXRLeadingOscillatoryField"
AGENT2_COMPOSITE_SOURCE_BLOB = "93a3da5c3aed4be8310cc510edfab6af940b64dc"

AGENT2_COMPOSITE_PARENT_PR = 981
AGENT2_COMPOSITE_PARENT_HEAD = "9997fc55455d126f935643da36bf17eaa0491aa4"

AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_MODULE = (
    "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
)
AGENT2_DIFFERENTIAL_FUNCTION = "evaluate_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

AGENT1_LEADING_PR = 980
AGENT1_LEADING_HEAD = "d3c971f2c62e272333e124e532212d23cca4908d"
AGENT1_LEADING_SOURCE_BLOB = "2ba28636a56b252fa485719e7b3e8da754d7e588"

SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"


def _validate_through_xr_configuration(
    configuration: Mapping[str, Any],
    *,
    runtime_semantic_sha256: str,
    runtime_oscillatory_sha256: str,
) -> None:
    """Validate the exact #987 identity/truth contract without recomputing A2 work."""
    if not isinstance(configuration, Mapping):
        raise ValueError("A2 #987 configuration must be a mapping")

    parent = configuration.get("parent_agent2")
    if not isinstance(parent, Mapping):
        raise ValueError("A2 #987 parent identity is missing")
    if parent.get("pr") != AGENT2_COMPOSITE_PARENT_PR:
        raise ValueError("A2 #987 parent PR identity drifted")
    if parent.get("head") != AGENT2_COMPOSITE_PARENT_HEAD:
        raise ValueError("A2 #987 parent head identity drifted")

    a1 = configuration.get("agent1_leading")
    if not isinstance(a1, Mapping):
        raise ValueError("A2 #987 Agent-1 leading identity is missing")
    if a1.get("pr") != AGENT1_LEADING_PR or a1.get("head") != AGENT1_LEADING_HEAD:
        raise ValueError("A2 #987 Agent-1 leading lineage drifted")
    if a1.get("source_blob_sha1") != AGENT1_LEADING_SOURCE_BLOB:
        raise ValueError("A2 #987 Agent-1 leading source identity drifted")
    if not _valid_sha256(a1.get("semantic_sha256")):
        raise ValueError("A2 #987 Agent-1 semantic identity is malformed")
    if not _valid_sha256(a1.get("configuration_sha256")):
        raise ValueError("A2 #987 Agent-1 configuration identity is malformed")

    oscillatory = configuration.get("oscillatory_runtime")
    if not isinstance(oscillatory, Mapping):
        raise ValueError("A2 #987 oscillatory runtime identity is missing")
    payload_sha = oscillatory.get("payload_sha256")
    if not _valid_sha256(payload_sha):
        raise ValueError("A2 #987 oscillatory runtime SHA is malformed")
    if payload_sha != runtime_oscillatory_sha256:
        raise ValueError("A2 #987 oscillatory runtime SHA disagrees with runtime")

    if not _valid_sha256(runtime_semantic_sha256):
        raise ValueError("A2 #987 composite semantic identity is malformed")

    truth = configuration.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("A2 #987 truth boundary is missing")
    required_true = (
        "current_leading_plus_oscillatory_velocity_through_XR_materialized",
        "velocity_beyond_Xh_through_XR_materialized",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_through_XR_composite_save_load_available",
        "reload_recomputes_leading_and_oscillatory_identities",
    )
    required_false = (
        "post_XR_velocity_materialized",
        "post_XR_RF40_current_lineage_materialized",
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
        "paper_exact",
        "pde_validated",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"A2 #987 truth boundary drifted at {key}")
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"A2 #987 truth boundary drifted at {key}")


@dataclass(frozen=True)
class ExactCurrentExteriorXRNonlinearBackend:
    """Checksum-bound A2 #987/#960 bridge into the existing A3 mean projector."""

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
    ) -> "ExactCurrentExteriorXRNonlinearBackend":
        field_type = type(composite_field)
        if field_type.__module__ != AGENT2_COMPOSITE_MODULE:
            raise ValueError("unexpected A2 #987 composite module identity")
        if field_type.__name__ != AGENT2_COMPOSITE_CLASS:
            raise ValueError("unexpected A2 #987 composite class identity")
        if not callable(getattr(composite_field, "velocity", None)):
            raise TypeError("A2 #987 composite must expose velocity(x,y,z,t)")
        if not callable(getattr(composite_field, "configuration", None)):
            raise TypeError("A2 #987 composite must expose configuration()")

        leading = getattr(composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise TypeError("A2 #987 composite must retain its authenticated leading backend")

        field_source = inspect.getsourcefile(field_type)
        if field_source is None:
            raise ValueError("A2 #987 composite source is unavailable")
        field_blob = _git_blob_sha1(field_source)
        if field_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise ValueError("A2 #987 composite source blob drifted")

        composite_semantic = str(getattr(composite_field, "semantic_sha256", ""))
        oscillatory_semantic = str(
            getattr(composite_field, "oscillatory_runtime_sha256", "")
        )
        configuration = composite_field.configuration()
        _validate_through_xr_configuration(
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


def materialize_current_exterior_xr_nonlinear_mean_attribution(
    backend: ExactCurrentExteriorXRNonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project the exact through-X_R nonlinear pieces; no surrogate defect is accepted."""
    if not isinstance(backend, ExactCurrentExteriorXRNonlinearBackend):
        raise TypeError("backend must be ExactCurrentExteriorXRNonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ValueError("A2 #987 composite provenance drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ValueError("A2 #960 differential provenance drifted")
    return _materialize_from_provider(
        backend.components,
        radius,
        z,
        t,
        backend=backend,  # runtime protocol is to_receipt(); generic witness is reusable
        backend_kind="exact-current-a2-pr987-through-xr-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_current_exterior_xr_nonlinear_mean_attribution
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
        "current_leading_plus_oscillation_through_XR_consumed": True,
        "velocity_beyond_Xh_through_XR_materialized": True,
        "current_through_XR_mixed_nonlinear_mean_materialized": True,
        "current_through_XR_quadratic_nonlinear_mean_materialized": True,
        "current_through_XR_aggregate_nonlinear_mean_materialized": True,
        "current_leading_spatial_derivative_is_fixed_repository_fd": True,
        "caller_tunable_leading_spatial_step": False,
        "agent2_oscillatory_jacobian_reimplemented_by_agent3": False,
        "post_XR_velocity_materialized": False,
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
