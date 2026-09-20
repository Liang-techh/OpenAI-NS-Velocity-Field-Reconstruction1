"""Fail-closed CR002 audit for Kokuno inner divergence norm semantics.

PR #823 deliberately reports a finite held-out sampled RMS of divergence.  The
canonical CR001 contract instead names a volume-weighted spatial L2 norm at
each validation time.  Equal numerical gates (currently 1e-5) do not make the
metrics or validation protocols interchangeable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


TASK_ID = "CR002-KOKUNO-INNER-DIVERGENCE-NORM-SCOPE-118"
UPSTREAM_HEAD = "56f4cf63e79fd08fe005ea3dce061f9a2f8ce50e"
UPSTREAM_BLOB = "9e6a4c95ef0f52d952c3278d028a201d85efa578"
UPSTREAM_PATH = (
    "src/openai_ns_reconstruction/"
    "kokuno_a4_pa10_cartesian_center_spatial_derivatives_independent_audit.py"
)
CLASSIFICATION_KEYS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class ScopeAuditError(RuntimeError):
    """Raised when the frozen norm-scope contract drifts or overclaims."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ScopeAuditError(message)


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(payload, dict), f"{path} must contain a JSON object")
    return payload


def _git_blob_sha(text: str) -> str:
    data = text.encode("utf-8")
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def audit(config_path: Path | None = None) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    cfg_path = config_path or root / "configs" / "kokuno_inner_divergence_norm_scope.json"
    cfg = _load_json(cfg_path)
    canonical = _load_json(root / "configs" / "constraints.json")

    _require(cfg.get("task_id") == TASK_ID, "task_id drift")
    upstream = cfg.get("upstream", {})
    _require(upstream.get("pr") == 823, "upstream PR drift")
    _require(upstream.get("head_sha") == UPSTREAM_HEAD, "upstream exact-head drift")
    _require(upstream.get("source_path") == UPSTREAM_PATH, "upstream source-path drift")
    _require(upstream.get("source_blob_sha") == UPSTREAM_BLOB, "upstream blob binding drift")

    classification = cfg.get("classification", {})
    _require(set(classification) == CLASSIFICATION_KEYS, "four-way provenance classification drift")
    _require(classification["public_source_fact"] == [], "project norm semantics cannot be laundered into public-source facts")
    for key in CLASSIFICATION_KEYS:
        _require(isinstance(classification[key], list), f"classification.{key} must be a list")

    source_path = root / UPSTREAM_PATH
    source = source_path.read_text(encoding="utf-8")
    _require(_git_blob_sha(source) == UPSTREAM_BLOB, "upstream source bytes no longer match frozen #823 blob")
    for token in (
        "SEED = 9173481",
        "SAMPLE_COUNT = 384",
        "DIVERGENCE_MAX_GATE = 1.0e-5",
        "DIVERGENCE_L2_GATE = 1.0e-5",
        '"independent_divergence_sampled_l2_rms"',
        "np.sqrt(np.mean(np.square(fine_divergence)))",
        "The reported divergence L2 quantity is a held-out sampled RMS, not a whole-domain volume integral norm.",
    ):
        _require(token in source, f"#823 sampled-divergence semantics drift: missing {token!r}")

    metric = cfg.get("metric_scope", {})
    _require(metric.get("upstream_metric_name") == "independent_divergence_sampled_l2_rms", "upstream metric-name drift")
    _require(metric.get("upstream_formula") == "sqrt(mean(divergence**2)) on PR #823 finite held-out samples", "upstream metric formula drift")
    _require(metric.get("upstream_random_sample_count") == 384, "upstream sample-count drift")
    _require(metric.get("upstream_sample_seed") == 9173481, "upstream sample-seed drift")
    _require(metric.get("upstream_threshold") == 1.0e-5, "upstream divergence gate drift")
    _require(metric.get("canonical_metric_name") == "divergence_L2", "canonical metric-name drift")
    _require(metric.get("canonical_norm_definition") == "volume-weighted L2 spatial norm at each time", "canonical L2 definition drift")
    _require(metric.get("canonical_threshold") == 1.0e-5, "canonical divergence_L2 gate drift")
    _require(metric.get("same_numeric_threshold") is True, "same-threshold fact drift")
    _require(metric.get("metrics_equivalent") is False, "sampled RMS must not be declared equivalent to canonical volume-weighted L2")
    _require(metric.get("upstream_pass_can_discharge_cr001_divergence_L2") is False, "#823 PASS cannot discharge canonical CR001 divergence_L2")
    _require(metric.get("upstream_pass_can_discharge_complete_candidate_divergence_acceptance") is False, "inner derivative audit cannot discharge complete-candidate divergence acceptance")
    _require(metric.get("upstream_pass_can_discharge_pde_validation") is False, "inner derivative audit cannot discharge PDE validation")

    snapshot = cfg.get("canonical_cr001_snapshot", {})
    expected_snapshot = {
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "forcing_parameter_bounds": {"a": [0.0, 10.0], "c": [0.0, 10.0]},
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "validation_times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
        "pde_residual_max": 1.0e-3,
        "pde_residual_L2": 1.0e-3,
        "residual_defined_pointwise_free_force_allowed": False,
        "amplitude_collapse_allowed": False,
        "post_hoc_threshold_relaxation_allowed": False,
    }
    _require(snapshot == expected_snapshot, "frozen CR001 snapshot drift")

    _require(canonical.get("nu") == snapshot["nu"], "canonical nu drift")
    domain = canonical.get("domain", {})
    _require(domain.get("physical") == snapshot["physical_domain"], "canonical physical-domain drift")
    _require(domain.get("evaluation_box") == snapshot["evaluation_box"], "canonical evaluation-box drift")
    _require(domain.get("support") == snapshot["support"], "canonical support drift")
    _require(domain.get("time_interval") == snapshot["time_interval"], "canonical time-window drift")
    forcing = canonical.get("forcing", {})
    _require(forcing.get("mode") == snapshot["forcing_mode"], "canonical forcing-mode drift")
    _require(forcing.get("parameters") == snapshot["forcing_parameter_bounds"], "canonical forcing bounds drift")
    restriction = str(forcing.get("restriction", ""))
    _require("No residual-dependent basis or pointwise free force" in restriction, "canonical free-force prohibition drift")
    nontriviality = canonical.get("nontriviality", {})
    _require(nontriviality.get("reference_energy") == snapshot["reference_energy"], "canonical reference-energy drift")
    _require(nontriviality.get("reference_energy_abs_tolerance") == snapshot["reference_energy_abs_tolerance"], "canonical energy tolerance drift")
    _require("reject collapsed candidates" in str(nontriviality.get("enforcement", "")), "canonical anti-collapse rule drift")
    validation = canonical.get("validation", {})
    _require(validation.get("seed") == snapshot["validation_seed"], "canonical validation-seed drift")
    _require(validation.get("held_out_points") == snapshot["held_out_points"], "canonical held-out-count drift")
    _require(validation.get("times") == snapshot["validation_times"], "canonical validation-time drift")
    _require(validation.get("derivative_steps") == snapshot["derivative_steps"], "canonical derivative-ladder drift")
    _require(validation.get("quadrature_orders_per_axis") == snapshot["quadrature_orders_per_axis"], "canonical quadrature-ladder drift")
    _require("volume-weighted L2 spatial norm at each time" in validation.get("norms", []), "canonical volume-weighted L2 semantics missing")
    thresholds = validation.get("thresholds", {})
    _require(thresholds.get("divergence_max") == snapshot["divergence_max"], "canonical divergence_max drift")
    _require(thresholds.get("divergence_L2") == snapshot["divergence_L2"], "canonical divergence_L2 drift")
    _require(thresholds.get("pde_residual_max") == snapshot["pde_residual_max"], "canonical momentum-max drift")
    _require(thresholds.get("pde_residual_L2") == snapshot["pde_residual_L2"], "canonical momentum-L2 drift")
    _require("changing thresholds requires a new experiment version" in str(validation.get("failure_policy", "")), "canonical no-post-hoc-relaxation rule drift")

    status = cfg.get("status_separation", {})
    _require(status.get("inner_spatial_derivative_admission") == "pending", "#823 admission must remain pending until exact-head evidence completes")
    for key in (
        "canonical_cr001_divergence_validated",
        "complete_kokuno_candidate_validated",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "callable_velocity_delivery_blocked_by_pde_pending",
    ):
        _require(status.get(key) is False, f"forbidden status promotion: {key}")

    forbidden = cfg.get("forbidden_promotions", [])
    _require(len(forbidden) >= 5, "forbidden-promotion registry unexpectedly weakened")

    return {
        "ok": True,
        "task_id": TASK_ID,
        "upstream_head": UPSTREAM_HEAD,
        "upstream_metric": metric["upstream_metric_name"],
        "canonical_metric": metric["canonical_metric_name"],
        "same_numeric_threshold": True,
        "metrics_equivalent": False,
        "scientific_promotions": [],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    args = parser.parse_args(argv)
    print(json.dumps(audit(args.config), sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
