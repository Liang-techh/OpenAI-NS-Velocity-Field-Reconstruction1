"""Replay Agent-7's energy-neutral swirl redistribution on frozen Agent-9 paths.

CR-A9-050 reuses the exact frozen path engine from Agent-9 PR #435, reconstructs
ST048-S from PR #390, applies Agent-7 PR #447's first clean inner/mid swirl
redistribution crossing (gain=.025), restores E(.25)=1 with one positive common
scale, and measures the unchanged 48 material paths.  No pressure/forcing/PDE
receipt is transferred and no OpenAI-hidden numerical target is used.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

TASK_ID = "CR-A9-050"
SCHEMA = "st048s_energy_neutral_swirl_redistribution_material_path_v1"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"

ST048_SOURCE_HEAD = "97695a86f85ce68fb4ae70c41fc81c904d655183"
ST048S_RAW_SHA256 = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
AGENT9_PATH_SOURCE_PR = 435
AGENT9_PATH_SOURCE_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
AGENT7_SOURCE_PR = 447
AGENT7_SOURCE_HEAD = "527e20e7eaa506d144c725d2c90c468da6aa5d34"
AGENT7_TASK_ID = "CR003-ST048S-ENERGY-NEUTRAL-SWIRL-REDISTRIBUTION-065"

REDISTRIBUTION_GAIN = 0.025
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
ENERGY_ORDER = 64
AGENT7_REPORTED_ALPHA = 2.415185613025038
AGENT7_REPORTED_SCALE = 0.9985868857
SEED_RADII = (0.6, 0.9, 1.2)

ST006_REFERENCE = {
    "mean_radius_change": -0.08601065505093726,
    "mean_absolute_turns": 0.024485286231482325,
    "maximum_absolute_turns": 0.0436185396886357,
    "mean_pair_separation_change": 0.07130825328092755,
    "mean_pair_separation_ratio": 1.1188470888015458,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
}
ST048S_REFERENCE = {
    "mean_radius_change": -0.0853545661,
    "mean_absolute_turns": 0.0221757297,
    "maximum_absolute_turns": 0.0346689104,
    "mean_pair_separation_change": 0.0778727503,
    "mean_pair_separation_ratio": 1.1297879172,
    "pair_growth_count": 16,
    "pair_shrink_count": 8,
}

SCIPY_SOURCE = {
    "repository": "scipy/scipy",
    "commit": "eff78058ed0d4cb8e1f0b5c5585d62d712e948d2",
    "api": "scipy.integrate.solve_ivp / DOP853",
    "license": "BSD-3-Clause",
    "classification": "direct migration / public API only (through pinned Agent-9 path engine)",
    "copied_upstream_implementation": False,
}
INTERNAL_SOURCES = [
    {
        "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "pr": AGENT9_PATH_SOURCE_PR,
        "commit": AGENT9_PATH_SOURCE_HEAD,
        "classification": "direct internal evidence reuse",
        "migrated_scope": "frozen 48-path contract, ST048-S loader, temporal Piola transform, energy quadrature, DOP853 path measurement",
    },
    {
        "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "pr": AGENT7_SOURCE_PR,
        "commit": AGENT7_SOURCE_HEAD,
        "task_id": AGENT7_TASK_ID,
        "classification": "direct internal method reuse / independent reimplementation",
        "migrated_scope": "C-infinity inner/outer radial bumps, first-order energy-neutral balance, swirl redistribution multiplier",
        "difference": "Eulerian proxy screen is replaced by actual frozen material trajectories",
    },
]
TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "production_candidate_selected": False,
    "production_redistribution_gain_selected": False,
    "pressure_or_force_transferred_from_parent": False,
    "held_out_pde_residual_evaluated": False,
    "held_out_pde_residual_transferred_from_parent": False,
    "visual_acceptance_threshold_defined": False,
    "comparison_is_descriptive_not_acceptance": True,
    "hidden_openai_time_camera_seed_or_velocity_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


def _git_head(root: str | Path) -> str:
    done = subprocess.run(["git", "-C", str(Path(root).resolve()), "rev-parse", "HEAD"],
                          check=True, capture_output=True, text=True)
    return done.stdout.strip()


def _hash(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _load_module_from_file(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_agent9_path_engine(agent9_root: str | Path):
    root = Path(agent9_root).resolve()
    if _git_head(root) != AGENT9_PATH_SOURCE_HEAD:
        raise ValueError("Agent-9 path checkout is not pinned PR #435 head")
    path = root / "src/openai_ns_reconstruction/constrained_st048s_piola_swirl_material_path.py"
    engine = _load_module_from_file("agent9_048_pinned", path)
    expected = {
        "ST048_SOURCE_HEAD": ST048_SOURCE_HEAD,
        "SEED_RADII": SEED_RADII,
        "SEED_Z": (-0.3, 0.3),
        "SEED_ANGLES": 8,
        "OUTPUT_SAMPLES": 33,
        "SOLVER_METHOD": "DOP853",
        "SOLVER_RTOL": 1e-9,
        "SOLVER_ATOL": 1e-11,
        "SOLVER_MAX_STEP": 0.01,
        "REGISTERED_TIME_INTERVAL": (0.25, 0.75),
    }
    for key, value in expected.items():
        if getattr(engine, key) != value:
            raise ValueError(f"pinned path contract drift: {key}")
    return engine


def compact_interval_bump(radius, lo: float, hi: float) -> np.ndarray:
    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius < 0):
        raise ValueError("radius must be finite and nonnegative")
    lo, hi = float(lo), float(hi)
    if not 0 <= lo < hi <= 2:
        raise ValueError("window must lie in [0,2]")
    mid, half = 0.5*(lo+hi), 0.5*(hi-lo)
    x = (radius-mid)/half
    out = np.zeros_like(radius)
    mask = np.abs(x) < 1
    xm = x[mask]
    out[mask] = np.exp(1.0 - 1.0/(1.0-xm*xm))
    return out


def inner_profile(radius):
    return compact_interval_bump(radius, *INNER_WINDOW)


def outer_profile(radius):
    return compact_interval_bump(radius, *OUTER_WINDOW)


def redistribution_profile(radius, alpha: float):
    alpha = float(alpha)
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive")
    return inner_profile(radius) - alpha*outer_profile(radius)


class RedistributionVelocity:
    def __init__(self, parent: VelocityCallable, *, gain: float, alpha: float, scale: float = 1.0):
        self.parent = parent
        self.gain = float(gain)
        self.alpha = float(alpha)
        self.scale = float(scale)
        if not callable(parent) or not 0 <= self.gain <= 0.05 or self.alpha <= 0 or self.scale <= 0:
            raise ValueError("invalid redistribution wrapper")
        rr = np.linspace(0, 2, 4001)
        if np.min(1 + self.gain*redistribution_profile(rr, self.alpha)) <= 0:
            raise ValueError("redistribution multiplier lost positivity")

    def __call__(self, points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        values = np.asarray(self.parent(points, float(time)), dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or values.shape != points.shape or not np.isfinite(values).all():
            raise ValueError("invalid points/velocity")
        r = np.hypot(points[:, 0], points[:, 1])
        out = values.copy()
        mask = r > 1e-14
        if np.any(mask):
            rx, ry = points[mask, 0]/r[mask], points[mask, 1]/r[mask]
            ux, uy = values[mask, 0], values[mask, 1]
            ur = rx*ux + ry*uy
            ut = -ry*ux + rx*uy
            ut *= 1 + self.gain*redistribution_profile(r[mask], self.alpha)
            out[mask, 0] = rx*ur - ry*ut
            out[mask, 1] = ry*ur + rx*ut
        return self.scale*out


def _swirl_moment(velocity: VelocityCallable, profile, *, order: int = ENERGY_ORDER) -> float:
    node, weight = leggauss(order)
    r, z = node+1, 2*node
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    values = np.asarray(velocity(points, 0.25), dtype=float)
    ut2 = (values[:, 1]**2).reshape(rr.shape)
    return float(np.sum((weight[:,None]*(2*weight)[None,:])*(np.pi*rr*ut2*profile(rr))))


def balance_coefficient(temporal_velocity: VelocityCallable) -> dict:
    inner = _swirl_moment(temporal_velocity, inner_profile)
    outer = _swirl_moment(temporal_velocity, outer_profile)
    if inner <= 0 or outer <= 0:
        raise RuntimeError("energy moments must be positive")
    alpha = float(inner/outer)
    defect = float(inner-alpha*outer)
    return {"inner_swirl_energy_moment": inner, "outer_swirl_energy_moment": outer,
            "outer_balance_coefficient_alpha": alpha, "balanced_moment_defect": defect,
            "kinetic_energy_first_derivative_at_zero": 2*defect}


def _summarize_by_seed_radius(measurement: dict) -> dict:
    groups = {f"{r:.1f}": [] for r in SEED_RADII}
    for row in measurement["per_path"]:
        groups[f"{float(row['seed']['radius']):.1f}"].append(row)
    pairs = {f"{r:.1f}": [] for r in SEED_RADII}
    for row in measurement["pair_rows"]:
        pairs[f"{float(row['radius']):.1f}"].append(row)
    out = {}
    for key, rows in groups.items():
        pr = pairs[key]
        out[key] = {
            "path_count": len(rows),
            "mean_absolute_turns": float(np.mean([r["absolute_turns"] for r in rows])),
            "maximum_absolute_turns": float(np.max([r["absolute_turns"] for r in rows])),
            "mean_radius_change": float(np.mean([r["radius_change"] for r in rows])),
            "mean_pair_axial_separation_change": float(np.mean([r["axial_separation_change"] for r in pr])),
        }
    return out


def _comparison(summary: dict, ref: dict) -> dict:
    return {
        "absolute_radial_contraction_ratio": abs(summary["mean_radius_change"])/abs(ref["mean_radius_change"]),
        "mean_absolute_turns_ratio": summary["mean_absolute_turns"]/ref["mean_absolute_turns"],
        "maximum_absolute_turns_ratio": summary["maximum_absolute_turns"]/ref["maximum_absolute_turns"],
        "mean_pair_separation_change_ratio": summary["mean_pair_axial_separation_change"]/ref["mean_pair_separation_change"],
        "pair_growth_count_delta": summary["pair_axial_separation_growth_count"]-ref["pair_growth_count"],
        "pair_shrink_count_delta": summary["pair_axial_separation_shrink_count"]-ref["pair_shrink_count"],
    }


def crosscheck_agent7(agent7_root: str | Path) -> dict:
    root = Path(agent7_root).resolve()
    if _git_head(root) != AGENT7_SOURCE_HEAD:
        raise ValueError("Agent-7 checkout is not pinned PR #447 head")
    exp = root / "experiments/root_st048"
    old = list(sys.path)
    names = ("agent7_st048s_energy_neutral_swirl_redistribution", "agent7_st048s_piola_swirl_gain_screen",
             "agent7_st048s_piola_temporal_screen", "agent7_st048s_piola_warp_screen", "agent7_st048_morphology_transfer")
    try:
        sys.path.insert(0, str(exp))
        for name in names:
            sys.modules.pop(name, None)
        screen = importlib.import_module(names[0])
        if screen.TASK_ID != AGENT7_TASK_ID:
            raise ValueError("Agent-7 task drift")
        r = np.linspace(0, 2, 1001)
        profile_error = float(np.max(np.abs(redistribution_profile(r, AGENT7_REPORTED_ALPHA)-screen.redistribution_profile(r, AGENT7_REPORTED_ALPHA))))
    finally:
        sys.path[:] = old
    if profile_error > 2e-15:
        raise RuntimeError("Agent-7 radial formula mismatch")
    return {"source_pr": AGENT7_SOURCE_PR, "source_head_sha": AGENT7_SOURCE_HEAD,
            "task_id": AGENT7_TASK_ID, "profile_max_abs_difference": profile_error,
            "crosscheck_passed": True}


def build_report(candidate_root: str | Path, agent9_root: str | Path, agent7_root: str | Path) -> dict:
    engine = load_agent9_path_engine(agent9_root)
    parent, candidate = engine.load_st048s_velocity(candidate_root)
    if candidate["source_head_sha"] != ST048_SOURCE_HEAD or candidate["original_raw_candidate_sha256"] != ST048S_RAW_SHA256:
        raise ValueError("ST048-S provenance drift")
    temporal = engine.TemporalPiolaSwirlVelocity(parent, kappa=0.0, scale=1.0)
    balance = balance_coefficient(temporal)
    alpha = balance["outer_balance_coefficient_alpha"]
    raw_base = RedistributionVelocity(temporal, gain=0.0, alpha=alpha)
    raw_child = RedistributionVelocity(temporal, gain=REDISTRIBUTION_GAIN, alpha=alpha)
    base_scale = float(np.sqrt(1.0/engine.axisymmetric_energy(raw_base, order=ENERGY_ORDER)))
    child_scale = float(np.sqrt(1.0/engine.axisymmetric_energy(raw_child, order=ENERGY_ORDER)))
    base = RedistributionVelocity(temporal, gain=0.0, alpha=alpha, scale=base_scale)
    child = RedistributionVelocity(temporal, gain=REDISTRIBUTION_GAIN, alpha=alpha, scale=child_scale)
    base_m = engine.measure_material_paths(base)
    child_m = engine.measure_material_paths(child)
    base_m["per_seed_radius"] = _summarize_by_seed_radius(base_m)
    child_m["per_seed_radius"] = _summarize_by_seed_radius(child_m)
    same_ref = {
        "mean_radius_change": base_m["mean_radius_change"], "mean_absolute_turns": base_m["mean_absolute_turns"],
        "maximum_absolute_turns": base_m["maximum_absolute_turns"], "mean_pair_separation_change": base_m["mean_pair_axial_separation_change"],
        "pair_growth_count": base_m["pair_axial_separation_growth_count"], "pair_shrink_count": base_m["pair_axial_separation_shrink_count"],
    }
    report = {
        "schema": SCHEMA, "task_id": TASK_ID, "base_main_sha": BASE_MAIN_SHA,
        "candidate": candidate, "external_source": SCIPY_SOURCE, "internal_sources": INTERNAL_SOURCES,
        "frozen_path_source": {"pr": AGENT9_PATH_SOURCE_PR, "head": AGENT9_PATH_SOURCE_HEAD},
        "frozen_contract": {"time_interval": [0.25,0.75], "seed_radii": list(SEED_RADII), "seed_z": [-0.3,0.3],
                            "seed_angles": 8, "path_count": 48, "paired_material_line_count": 24,
                            "output_samples": 33, "solver_method": "DOP853", "solver_rtol": 1e-9,
                            "solver_atol": 1e-11, "solver_max_step": 0.01},
        "agent7_formula_crosscheck": crosscheck_agent7(agent7_root),
        "normalization": {"redistribution_gain": REDISTRIBUTION_GAIN, "balance": balance,
                          "beta_at_025": engine.beta_at_time(0.25), "beta_at_050": engine.beta_at_time(0.50), "beta_at_075": engine.beta_at_time(0.75),
                          "profile_at_seed_radii": {f"{r:.1f}": float(redistribution_profile(np.array([r]),alpha)[0]) for r in SEED_RADII},
                          "baseline_scale": base_scale, "child_scale": child_scale,
                          "baseline_energy": engine.axisymmetric_energy(base, order=ENERGY_ORDER),
                          "child_energy": engine.axisymmetric_energy(child, order=ENERGY_ORDER)},
        "measurements": {"redistribution_gain_0": base_m, "redistribution_gain_0025": child_m},
        "comparisons": {"gain_0025_vs_same_run_gain_0": _comparison(child_m, same_ref),
                        "gain_0025_vs_unwarped_st048s": _comparison(child_m, ST048S_REFERENCE),
                        "gain_0025_vs_st006": _comparison(child_m, ST006_REFERENCE)},
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    report["report_sha256"] = _hash(report)
    return report


def audit_truth_boundary(report: dict) -> None:
    if report.get("task_id") != TASK_ID or report.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("task/truth boundary drift")
    if report["candidate"].get("pde_validated") is not False or report["agent7_formula_crosscheck"].get("crosscheck_passed") is not True:
        raise ValueError("candidate/source truth drift")
    for key, value in TRUTH_BOUNDARY.items():
        if key != "comparison_is_descriptive_not_acceptance" and value is False and report["truth_boundary"][key] is not False:
            raise ValueError(f"forbidden promotion: {key}")
    if report["truth_boundary"]["comparison_is_descriptive_not_acceptance"] is not True:
        raise ValueError("comparison must remain descriptive")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-root", required=True); p.add_argument("--agent9-root", required=True)
    p.add_argument("--agent7-root", required=True); p.add_argument("--output", required=True)
    a = p.parse_args(argv)
    report = build_report(a.candidate_root, a.agent9_root, a.agent7_root)
    audit_truth_boundary(report)
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    b, c = report["measurements"]["redistribution_gain_0"], report["measurements"]["redistribution_gain_0025"]
    print(json.dumps({"task_id": TASK_ID, "report_sha256": report["report_sha256"],
                      "alpha": report["normalization"]["balance"]["outer_balance_coefficient_alpha"],
                      "child_scale": report["normalization"]["child_scale"],
                      "baseline_mean_turns": b["mean_absolute_turns"], "child_mean_turns": c["mean_absolute_turns"],
                      "baseline_mean_dr": b["mean_radius_change"], "child_mean_dr": c["mean_radius_change"],
                      "per_seed_radius": c["per_seed_radius"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
