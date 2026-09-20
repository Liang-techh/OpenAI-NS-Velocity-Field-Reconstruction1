"""Fail-closed CR002 audit for PR #875 first-cell extrapolation observability.

PR #875 improves the autonomous first-cell quadrature by integrating a quadratic
Lagrange extrapolant of the unweighted source.  It also records a degree-1 versus
degree-2 disagreement as a local sensitivity diagnostic.  This audit protects a
narrower semantic boundary: agreement of two extrapolants built only from
positive-radius nodes is not, without an independently justified regularity or
shape bound on the real source inside ``[0, r_min]``, a formal error certificate
for the unsampled first-cell integral.

Nothing here changes the #875 operator.  The manufactured hidden-cell function
is mechanics only and is deliberately not asserted to be a Kokuno/OpenAI source
profile.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import kokuno_strict_inner_transport_first_cell_polynomial_radial_stress as audited

SCHEMA = "cr002-kokuno-first-cell-extrapolation-observability-scope-v1"
TASK = "CR002-KOKUNO-FIRST-CELL-EXTRAPOLATION-OBSERVABILITY-SCOPE-119"
AUDITED_HEAD = "4b302412fe5733690394002a947eeedaa61a1835"
AUDITED_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"
CANONICAL_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
CONTRACT_REL = Path("configs/kokuno_first_cell_extrapolation_observability_scope.json")
AUDITED_REL = Path(
    "src/openai_ns_reconstruction/"
    "kokuno_strict_inner_transport_first_cell_polynomial_radial_stress.py"
)
CONSTRAINTS_REL = Path("configs/constraints.json")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    framed = b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload
    return hashlib.sha1(framed).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_contract(path: Path | None = None) -> dict[str, Any]:
    target = path or (_repo_root() / CONTRACT_REL)
    return json.loads(target.read_text())


def _formal_hidden_cell_integral(a: float, exponent: int, amplitude: float = 1.0) -> float:
    """Integral of amplitude*r^e*(a-r)^4 over [0,a]."""

    _require(exponent in (1, 2), "hidden-cell witness only registers e=1,2")
    return float(
        amplitude
        * a ** (exponent + 5)
        * math.factorial(exponent)
        * math.factorial(4)
        / math.factorial(exponent + 5)
    )


def hidden_cell_mechanics_witness(a: float = 0.2, amplitude: float = 1.0) -> dict[str, Any]:
    """Show a sample-only null direction for both #875 extrapolants.

    The hypothetical function is ``amplitude*(a-r)^4`` on ``0<=r<a`` and zero
    at/above ``a``.  Every positive node supplied to #875 is at least ``a``, so
    the degree-1 and degree-2 estimators both see only zeros.  The actual weighted
    first-cell integral is nevertheless positive.  This is an observability
    witness only; it is not a statement that the real source belongs to this
    piecewise-polynomial class.
    """

    radii = np.asarray([a, 1.5 * a, 2.0 * a, 2.5 * a], dtype=float)
    samples = np.zeros_like(radii)
    channels: dict[str, Any] = {}
    for exponent in (1, 2):
        affine = audited._first_cell_polynomial_integral(
            samples, radii, exponent=exponent, degree=1
        )
        quadratic = audited._first_cell_polynomial_integral(
            samples, radii, exponent=exponent, degree=2
        )
        formal = _formal_hidden_cell_integral(a, exponent, amplitude)
        channels[f"e{exponent}"] = {
            "degree_1_estimate": float(affine),
            "degree_2_estimate": float(quadratic),
            "degree_1_degree_2_disagreement": abs(float(affine - quadratic)),
            "formal_weighted_integral": formal,
            "formal_minus_degree_2": formal - float(quadratic),
        }
    return {
        "mechanics_only": True,
        "real_source_claim": False,
        "public_source_fact": False,
        "r_min": float(a),
        "amplitude": float(amplitude),
        "sample_radii": radii.tolist(),
        "sample_values": samples.tolist(),
        "channels": channels,
    }


def _audit_provenance(contract: Mapping[str, Any]) -> None:
    provenance = contract["provenance"]
    expected = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
    _require(set(provenance) == expected, "four provenance classes must remain exact and separate")
    flattened: list[str] = []
    for bucket in expected:
        values = provenance[bucket]
        _require(isinstance(values, list) and values, f"provenance bucket {bucket} must be nonempty")
        _require(all(isinstance(value, str) and value for value in values), f"invalid {bucket} entry")
        flattened.extend(values)
    _require(len(flattened) == len(set(flattened)), "provenance entries must not be duplicated across buckets")
    autonomous = "\n".join(provenance["autonomous_design"])
    public = "\n".join(provenance["public_source_fact"])
    pending = "\n".join(provenance["pending_unknown"])
    _require("degree-2 Lagrange" in autonomous, "#875 quadratic extrapolation must remain autonomous")
    _require("hidden-cell mechanics witness" in autonomous, "mechanics witness must remain autonomous")
    _require("hidden-cell" not in public.lower(), "mechanics witness cannot be laundered into public-source fact")
    _require("real-source first-cell convergence" in pending, "real-source convergence must remain pending")


def _audit_cr001(contract: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    c = contract["canonical_cr001"]
    _require(c["nu"] == constraints["nu"] == 0.01, "nu drift")
    _require(c["physical_domain"] == constraints["domain"]["physical"] == "R^3", "domain drift")
    _require(c["evaluation_box"] == constraints["domain"]["evaluation_box"], "evaluation box drift")
    _require(c["support"] == constraints["domain"]["support"], "support drift")
    _require(c["time_interval"] == constraints["domain"]["time_interval"] == [0.25, 0.75], "time interval drift")
    _require(c["forcing_mode"] == constraints["forcing"]["mode"] == "restricted_two_parameter_family", "forcing family drift")
    _require(c["forcing_parameter_bounds"] == constraints["forcing"]["parameters"], "forcing bounds drift")
    _require(c["reference_energy"] == constraints["nontriviality"]["reference_energy"] == 1.0, "energy target drift")
    _require(
        c["reference_energy_abs_tolerance"]
        == constraints["nontriviality"]["reference_energy_abs_tolerance"]
        == 0.001,
        "energy tolerance drift",
    )
    validation = constraints["validation"]
    _require(c["validation_seed"] == validation["seed"] == 914027, "validation seed drift")
    _require(c["held_out_points"] == validation["held_out_points"] == 4096, "held-out size drift")
    _require(c["derivative_steps"] == validation["derivative_steps"] == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(c["quadrature_orders_per_axis"] == validation["quadrature_orders_per_axis"] == [24, 48, 96], "quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key, expected in (
        ("divergence_max", 1e-5),
        ("divergence_L2", 1e-5),
        ("pde_residual_max", 1e-3),
        ("pde_residual_L2", 1e-3),
    ):
        _require(c[key] == thresholds[key] == expected, f"CR001 threshold drift: {key}")
    _require(c["residual_defined_pointwise_free_force_allowed"] is False, "free residual-defined force enabled")
    _require(c["amplitude_collapse_allowed"] is False, "amplitude collapse enabled")
    _require(c["post_hoc_threshold_relaxation_allowed"] is False, "post-hoc threshold relaxation enabled")


def audit_contract(contract: Mapping[str, Any] | None = None) -> dict[str, Any]:
    root = _repo_root()
    data = copy.deepcopy(dict(contract)) if contract is not None else load_contract()

    _require(data["schema"] == SCHEMA, "scope schema drift")
    _require(data["task"] == TASK, "scope task drift")
    _require(data["audited_pr"] == 875, "audited PR drift")
    _require(data["audited_head"] == AUDITED_HEAD, "audited exact-head drift")
    _require(data["audited_source_blob"] == AUDITED_SOURCE_BLOB, "audited source-blob declaration drift")
    _require(data["canonical_constraints_blob"] == CANONICAL_CONSTRAINTS_BLOB, "canonical constraints-blob declaration drift")
    _require(_git_blob_sha(root / AUDITED_REL) == AUDITED_SOURCE_BLOB, "#875 audited source bytes drifted")
    _require(_git_blob_sha(root / CONSTRAINTS_REL) == CANONICAL_CONSTRAINTS_BLOB, "canonical CR001 bytes drifted")

    _audit_provenance(data)

    realization = data["audited_realization"]
    _require(realization["production_first_cell_degree"] == 2, "production first-cell degree drift")
    _require(realization["sensitivity_first_cell_degree"] == 1, "sensitivity first-cell degree drift")
    _require(realization["both_extrapolants_use_only_positive_radius_samples"] is True, "positive-node observability changed")
    _require(realization["samples_inside_open_interval_0_rmin_used"] is False, "contract falsely claims sub-first-cell samples")
    _require(realization["linear_quadratic_disagreement_recorded"] is True, "#875 sensitivity diagnostic lost")
    for key in (
        "linear_quadratic_disagreement_is_formal_error_bound",
        "polynomial_exactness_through_degree_2_is_real_source_error_bound",
        "real_source_first_cell_convergence_verified",
        "formal_first_cell_accuracy_verified",
        "formal_axis_based_inverse_verified",
        "formal_full_domain_moment_verified",
        "axis_regularity_verified",
    ):
        _require(realization[key] is False, f"forbidden numerical/formal promotion: {key}")

    upstream_truth = audited.truth_boundary()
    for key in (
        "real_source_first_cell_convergence_verified",
        "formal_first_cell_accuracy_verified",
        "formal_axis_based_inverse_verified",
        "formal_full_domain_moment_verified",
        "axis_regularity_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(upstream_truth[key] is False, f"#875 upstream truth boundary drifted: {key}")

    witness = hidden_cell_mechanics_witness(
        float(data["observability_boundary"]["witness_r_min"]),
        float(data["observability_boundary"]["witness_amplitude"]),
    )
    registered = data["observability_boundary"]
    _require(registered["hidden_cell_witness_is_real_source_claim"] is False, "mechanics witness promoted to real-source claim")
    _require(registered["hidden_cell_witness_is_public_source_fact"] is False, "mechanics witness promoted to public-source fact")
    _require(registered["hidden_cell_witness_is_mechanics_only"] is True, "mechanics-only label lost")
    _require(witness["sample_values"] == [0.0, 0.0, 0.0, 0.0], "hidden-cell witness must be invisible at all positive nodes")
    for exponent in (1, 2):
        channel = witness["channels"][f"e{exponent}"]
        _require(channel["degree_1_estimate"] == 0.0, f"e={exponent} affine estimator should see zero")
        _require(channel["degree_2_estimate"] == 0.0, f"e={exponent} quadratic estimator should see zero")
        _require(channel["degree_1_degree_2_disagreement"] == 0.0, f"e={exponent} estimator disagreement should be zero")
        _require(channel["formal_weighted_integral"] > 0.0, f"e={exponent} hidden formal contribution must be positive")
    _require(
        math.isclose(witness["channels"]["e1"]["formal_weighted_integral"], registered["formal_weighted_integral_e1"], rel_tol=2e-15, abs_tol=0.0),
        "registered e=1 hidden-cell integral drift",
    )
    _require(
        math.isclose(witness["channels"]["e2"]["formal_weighted_integral"], registered["formal_weighted_integral_e2"], rel_tol=2e-15, abs_tol=0.0),
        "registered e=2 hidden-cell integral drift",
    )

    constraints = json.loads((root / CONSTRAINTS_REL).read_text())
    _audit_cr001(data, constraints)

    status = data["status_firewall"]
    for key in (
        "velocity_bytes_changed",
        "pressure_changed",
        "forcing_changed",
        "radial_stress_values_changed",
        "source_formula_changed",
        "candidate_coefficients_changed",
        "scientific_thresholds_changed",
        "velocity_export_ready_promoted_by_this_scope",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pde_pending_blocks_callable_velocity_delivery",
    ):
        _require(status[key] is False, f"forbidden status promotion/change: {key}")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "audited_head": AUDITED_HEAD,
        "audited_source_blob": AUDITED_SOURCE_BLOB,
        "canonical_constraints_blob": CANONICAL_CONSTRAINTS_BLOB,
        "hidden_cell_mechanics": witness,
        "scope_ok": True,
        "truth_boundary": {
            "affine_quadratic_agreement_is_error_certificate": False,
            "formal_first_cell_accuracy_verified": False,
            "formal_axis_based_inverse_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main() -> int:
    print(json.dumps(audit_contract(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
