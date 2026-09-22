"""Current-I4 nonlinear mean attribution for the exact A1/A2 candidate identity.

This Kokuno Agent-3 increment consumes Agent-2 PR #1080

    u_current_I4 = u_lead,A1#1079 + u_osc,frozen-complete-curl

and the already-frozen Agent-2 PR #960 oscillatory differential interface.  It
forms only the existing Agent-3 nonlinear decomposition

    A = (u_lead . grad) u_osc,
    B = (u_osc . grad) u_lead,
    Q = (u_osc . grad) u_osc,
    N = A + B + Q,

then reuses the repository cylindrical m=0 projector.  The corrected Kokuno
reader identifies I4 as the reserved mean-correction interval (RF40a), but the
source five-row mean correction (RF30--RF49) is *not* executed here.  This
increment first materializes the actual current-I4 nonlinear mean input on one
semantic candidate identity; it does not accept a caller-supplied surrogate
residual/defect.
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

TASK = "KOKUNO-A3-CURRENT-I4-NONLINEAR-MEAN-119"
SCHEMA = "kokuno-a3-current-i4-nonlinear-mean-v1"

PARENT_AGENT3_PR = 1090
PARENT_AGENT3_HEAD = "5e8512392df683375fabd132f35577e1b94773ff"

AGENT2_COMPOSITE_PR = 1080
AGENT2_COMPOSITE_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
AGENT2_COMPOSITE_MODULE = "openai_ns_reconstruction.kokuno_current_i4_leading_oscillatory_identity"
AGENT2_COMPOSITE_CLASS = "CurrentI4LeadingOscillatoryField"
AGENT2_COMPOSITE_SOURCE_BLOB = "2a0a5aa5966b02da856bdcf51940f3c186042802"

AGENT2_COMPOSITE_PARENT_PR = 1071
AGENT2_COMPOSITE_PARENT_HEAD = "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
AGENT2_COMPOSITE_PARENT_MODULE = "openai_ns_reconstruction.kokuno_current_i2_leading_oscillatory_identity"
AGENT2_COMPOSITE_PARENT_SOURCE_BLOB = "be68248aca97173a15474324f63efff2c9ffce56"

AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_MODULE = "openai_ns_reconstruction.kokuno_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_FUNCTION = "evaluate_oscillatory_batch_differentials"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

AGENT1_LEADING_PR = 1079
AGENT1_LEADING_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
AGENT1_LEADING_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_i4_leading_preservation"
AGENT1_LEADING_CLASS = "KokunoPA16CurrentCartesianI4LeadingPreservation"
AGENT1_LEADING_SOURCE_BLOB = "6f04ce0a856b44430402576dad88438da90d1ebb"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_MEAN_INTERVAL_FORMULA = "RF40a: I_mean = I4"
SOURCE_MEAN_DEFECT_FORMULAS = "RF30--RF49"


def _expected_source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_READER_HEAD,
        "path": SOURCE_PATH,
        "blob_sha1": SOURCE_BLOB_SHA1,
        "corrected_release_date": SOURCE_READER_DATE,
        "classification": "structural_math_provenance_only",
        "scope": "localized waves / complete curls plus I3/I4 source-role organization",
        "paper_exact_claim": False,
    }


def _validate_current_i4_configuration(
    configuration: Mapping[str, Any],
    *,
    runtime_semantic_sha256: str,
    runtime_oscillatory_sha256: str,
) -> None:
    """Fail closed unless the runtime is exact A2 #1080 / A1 #1079."""
    if not isinstance(configuration, Mapping):
        raise ValueError("A2 #1080 configuration must be a mapping")
    if configuration.get("schema") != "kokuno-a2-current-i4-leading-oscillatory-identity-v1":
        raise ValueError("A2 #1080 configuration schema drifted")
    if configuration.get("task") != "K2-OSC-097":
        raise ValueError("A2 #1080 task identity drifted")

    parent = configuration.get("parent_agent2")
    if not isinstance(parent, Mapping):
        raise ValueError("A2 #1080 parent identity is missing")
    expected_parent = {
        "pr": AGENT2_COMPOSITE_PARENT_PR,
        "head": AGENT2_COMPOSITE_PARENT_HEAD,
        "module": AGENT2_COMPOSITE_PARENT_MODULE,
        "source_blob_sha1": AGENT2_COMPOSITE_PARENT_SOURCE_BLOB,
    }
    for key, expected in expected_parent.items():
        if parent.get(key) != expected:
            raise ValueError(f"A2 #1080 parent identity drifted at {key}")

    leading = configuration.get("agent1_leading")
    if not isinstance(leading, Mapping):
        raise ValueError("A2 #1080 Agent-1 identity is missing")
    expected_leading = {
        "pr": AGENT1_LEADING_PR,
        "head": AGENT1_LEADING_HEAD,
        "module": AGENT1_LEADING_MODULE,
        "class": AGENT1_LEADING_CLASS,
        "source_blob_sha1": AGENT1_LEADING_SOURCE_BLOB,
    }
    for key, expected in expected_leading.items():
        if leading.get(key) != expected:
            raise ValueError(f"A2 #1080 Agent-1 identity drifted at {key}")
    if not _valid_sha256(leading.get("semantic_sha256")):
        raise ValueError("A2 #1080 Agent-1 semantic identity is malformed")
    if not _valid_sha256(leading.get("configuration_sha256")):
        raise ValueError("A2 #1080 Agent-1 configuration identity is malformed")

    oscillatory = configuration.get("oscillatory_runtime")
    if not isinstance(oscillatory, Mapping):
        raise ValueError("A2 #1080 oscillatory runtime identity is missing")
    payload_sha = oscillatory.get("payload_sha256")
    if not _valid_sha256(payload_sha):
        raise ValueError("A2 #1080 oscillatory runtime SHA is malformed")
    if payload_sha != runtime_oscillatory_sha256:
        raise ValueError("A2 #1080 oscillatory runtime SHA disagrees with runtime")

    if not _valid_sha256(runtime_semantic_sha256):
        raise ValueError("A2 #1080 composite semantic identity is malformed")
    if configuration.get("source_provenance") != _expected_source_provenance():
        raise ValueError("A2 #1080 corrected-source provenance drifted")
    if configuration.get("composition") != "u_current_I4 = u_lead_A1_1079 + u_osc_frozen_complete_curl":
        raise ValueError("A2 #1080 composition contract drifted")

    truth = configuration.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("A2 #1080 truth boundary is missing")
    required_true = (
        "current_leading_through_I4_consumed",
        "current_I4_leading_plus_frozen_complete_curl_oscillation_materialized",
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_current_I4_composite_save_load_available",
        "reload_recomputes_leading_and_oscillatory_identities",
        "a2_semantic_identity_binds_corrected_source_provenance_block",
    )
    required_false = (
        "field_summands_changed_by_this_increment",
        "source_positive_order_I3_profiles_materialized",
        "I3_positive_order_correction_materialized",
        "source_I4_mean_correction_materialized",
        "current_I4_mean_correction_materialized",
        "leading_I4_overlay_invented",
        "source_terminal_tail_schedule_bound_into_current_velocity",
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
        "source_forcing_provider_materialized",
        "source_background_provider_materialized",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"A2 #1080 truth boundary drifted at {key}")
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"A2 #1080 truth boundary drifted at {key}")


@dataclass(frozen=True)
class ExactCurrentI4NonlinearBackend:
    """Checksum-bound A2 #1080/#960 bridge into the existing A3 mean projector."""

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
    ) -> "ExactCurrentI4NonlinearBackend":
        field_type = type(composite_field)
        if field_type.__module__ != AGENT2_COMPOSITE_MODULE:
            raise ValueError("unexpected A2 #1080 composite module identity")
        if field_type.__name__ != AGENT2_COMPOSITE_CLASS:
            raise ValueError("unexpected A2 #1080 composite class identity")
        if not callable(getattr(composite_field, "velocity", None)):
            raise TypeError("A2 #1080 composite must expose velocity(x,y,z,t)")
        if not callable(getattr(composite_field, "configuration", None)):
            raise TypeError("A2 #1080 composite must expose configuration()")

        leading = getattr(composite_field, "leading_backend", None)
        if leading is None or not callable(getattr(leading, "velocity", None)):
            raise TypeError("A2 #1080 composite must retain its authenticated leading backend")

        field_source = inspect.getsourcefile(field_type)
        if field_source is None:
            raise ValueError("A2 #1080 composite source is unavailable")
        field_blob = _git_blob_sha1(field_source)
        if field_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
            raise ValueError("A2 #1080 composite source blob drifted")

        composite_semantic = str(getattr(composite_field, "semantic_sha256", ""))
        oscillatory_semantic = str(getattr(composite_field, "oscillatory_runtime_sha256", ""))
        _validate_current_i4_configuration(
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
            "agent1_leading_source_blob": AGENT1_LEADING_SOURCE_BLOB,
            "fixed_leading_spatial_step": FIXED_LEADING_SPATIAL_STEP,
            "leading_derivative_realization": "centered_cartesian_fd2_repository_fixed",
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_blob": SOURCE_BLOB_SHA1,
            "source_mean_interval_formula": SOURCE_MEAN_INTERVAL_FORMULA,
            "source_mean_defect_formulas": SOURCE_MEAN_DEFECT_FORMULAS,
        }


def materialize_current_i4_nonlinear_mean_attribution(
    backend: ExactCurrentI4NonlinearBackend,
    radius: Any,
    z: Any,
    t: Any,
) -> OscillatoryNonlinearMeanAttributionWitness:
    """Project exact current-I4 nonlinear pieces; no surrogate defect is accepted."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    if backend.composite_source_blob != AGENT2_COMPOSITE_SOURCE_BLOB:
        raise ValueError("A2 #1080 composite provenance drifted")
    if backend.differential_source_blob != AGENT2_DIFFERENTIAL_SOURCE_BLOB:
        raise ValueError("A2 #960 differential provenance drifted")
    return _materialize_from_provider(
        backend.components,
        radius,
        z,
        t,
        backend=backend,
        backend_kind="exact-current-a2-pr1080-through-i4-nonlinear-attribution",
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_nonlinear_mean_attribution)
    forbidden = {
        "residual", "defect", "mean", "source", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "angular_order",
        "spatial_step", "derivative_step", "time_step", "viscosity", "nu",
        "normalized_score", "scientific_threshold", "correction", "stage_budget",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_composite_source_blob": AGENT2_COMPOSITE_SOURCE_BLOB,
        "agent2_differential_pr": AGENT2_DIFFERENTIAL_PR,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent2_differential_source_blob": AGENT2_DIFFERENTIAL_SOURCE_BLOB,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "agent1_leading_source_blob": AGENT1_LEADING_SOURCE_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_blob": SOURCE_BLOB_SHA1,
        "source_reader_date": SOURCE_READER_DATE,
        "source_mean_interval_formula": SOURCE_MEAN_INTERVAL_FORMULA,
        "source_mean_defect_formulas": SOURCE_MEAN_DEFECT_FORMULAS,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(signature.parameters)),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "current_I4_leading_plus_oscillatory_identity_consumed": True,
        "current_I4_mixed_nonlinear_mean_materialized": True,
        "current_I4_quadratic_oscillatory_mean_materialized": True,
        "current_I4_aggregate_nonlinear_mean_materialized": True,
        "source_I4_reserved_mean_correction_interval_recorded": True,
        "agent2_oscillatory_curl_or_jacobian_reimplemented": False,
        "source_I3_positive_order_correction_materialized": False,
        "source_I4_five_row_mean_correction_materialized": False,
        "current_I4_correction_velocity_materialized": False,
        "post_I4_or_pulse_leading_identity_consumed": False,
        "outer_global_velocity_consumed": False,
        "radial_inverse_performed_in_this_increment": False,
        "compact_radial_stress_materialized_in_this_increment": False,
        "radial_force_materialized_in_this_increment": False,
        "scoped_nonlinear_mean_authorized_as_complete_ns_defect": False,
        "scoped_nonlinear_mean_authorized_as_correction_target": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_defect": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proof_claimed": False,
        "pde_validated": False,
    }
