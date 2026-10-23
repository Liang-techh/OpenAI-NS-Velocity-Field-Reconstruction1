"""Small, reproducible ansatz-capacity probe around Agent 1's compact candidate.

This lane does not replace or extend ``CompactCandidate``.  It asks a narrower
question: does a nested set of compact divergence-free higher-order corrections
span useful directions of the independently evaluated Navier--Stokes residual?
Pressure and the preregistered restricted force are frozen.  Coefficients are fit
only on one sample and the *actual nonlinear* residual is evaluated on disjoint
held-out samples.

The experiment is diagnostic only.  It is not CR003/CR005 acceptance and it does
not turn a linearized improvement into a validated Navier--Stokes solution.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_candidate import CompactCandidate
from openai_ns_reconstruction.constrained_force import RestrictedForce, compact_bump
from openai_ns_reconstruction.constrained_validation import residual

MODE_NAMES = (
    "poloidal_r4",
    "swirl_r4",
    "poloidal_r2z2",
    "swirl_r2z2",
    "poloidal_z4",
    "swirl_z4",
    "poloidal_r2z2_time",
    "swirl_r2z2_time",
)
# k=7 isolates the first time-dependent (poloidal) direction; k=8 then adds
# exactly one time-dependent swirl direction.
BASIS_SIZES = (0, 2, 4, 6, 7, 8)


@dataclass(frozen=True)
class ProbeConfig:
    points_per_time: int = 96
    train_times: tuple[float, ...] = (0.3125, 0.5, 0.6875)
    holdout_times: tuple[float, ...] = (0.375, 0.5625, 0.75)
    seed_pairs: tuple[tuple[int, int], ...] = (
        (7017, 17017),
        (7018, 17018),
        (7019, 17019),
    )
    box_halfwidth: float = 1.5
    derivative_step: float = 0.005
    jacobian_step: float = 2.0e-4
    nu: float = 0.01


def _load_candidate(path: str | Path) -> CompactCandidate:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or data.get("status") != "candidate":
        raise ValueError("expected a schema-v1 candidate artifact")
    return CompactCandidate(**data["parameters"])


def _load_force(path: str | Path) -> RestrictedForce:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    force = data.get("force")
    if not isinstance(force, dict) or set(force) != {"a", "c"}:
        raise ValueError("training artifact must contain exactly force coefficients a,c")
    return RestrictedForce(**force)


def _phi_and_derivatives(kind: int, X: np.ndarray, Z: np.ndarray, time: float):
    """Return phi, dphi/dX, dphi/dZ for the normalized higher-order modes."""
    x = X / 4.0
    z2 = Z * Z / 4.0
    if kind == 0:  # r^4
        return x * x, X / 8.0, np.zeros_like(X)
    if kind == 1:  # r^2 z^2
        return x * z2, z2 / 4.0, x * Z / 2.0
    if kind == 2:  # z^4
        return z2 * z2, np.zeros_like(X), Z**3 / 4.0
    if kind == 3:  # centered time-dependent r^2 z^2
        temporal = (time - 0.5) / 0.25
        return temporal * x * z2, temporal * z2 / 4.0, temporal * x * Z / 2.0
    raise ValueError("unknown basis kind")


def correction_mode(candidate: CompactCandidate, mode_index: int, points, time):
    """Evaluate one compact, axis-regular, divergence-free correction mode.

    Even modes are the curl of ``(-y H, x H, 0)`` with
    ``H=tau^-A Z b phi``.  Odd modes are smooth axisymmetric swirl
    ``(-y s, x s, 0)``.  No free forcing is introduced.
    """
    if isinstance(mode_index, bool) or not isinstance(mode_index, int):
        raise TypeError("mode_index must be an int")
    if not 0 <= mode_index < len(MODE_NAMES):
        raise ValueError("invalid mode index")
    p = np.asarray(points, dtype=float)
    if p.ndim < 1 or p.shape[-1] != 3 or not np.isfinite(p).all():
        raise ValueError("points must be finite with shape (...,3)")
    t = float(time)
    if not np.isfinite(t) or not 0.25 <= t <= 0.75:
        raise ValueError("time must lie in [0.25,0.75]")

    x, y, z = np.moveaxis(p, -1, 0)
    tau = 1.0 - t
    lr = candidate.radial_width * np.sqrt(tau)
    lz = candidate.axial_width * tau**0.495
    X = (x * x + y * y) / lr**2
    Z = z / lz
    br, dbr = compact_bump(X / 4.0)
    bz, dbz = compact_bump(Z * Z / 4.0)
    b = br * bz
    # X*bX_radial = r*d_r b; bZ = d_Z b.
    bX_radial = dbr * bz / 2.0
    bZ = br * dbz * Z / 2.0
    kind = mode_index // 2 if mode_index < 6 else 3
    phi, phiX, phiZ = _phi_and_derivatives(kind, X, Z, t)
    scale = tau**(-0.505)

    if mode_index % 2 == 0:
        radial = -scale / lz * (b * phi + Z * (bZ * phi + b * phiZ))
        axial = scale * Z * (
            2.0 * b * phi + X * bX_radial * phi + 2.0 * X * b * phiX
        )
        return np.stack((x * radial, y * radial, axial), axis=-1)

    swirl = scale / lr * b * phi
    return np.stack((-y * swirl, x * swirl, np.zeros_like(x)), axis=-1)


class CorrectedVelocity:
    def __init__(self, candidate: CompactCandidate, coefficients):
        self.candidate = candidate
        self.coefficients = np.asarray(coefficients, dtype=float)
        if self.coefficients.ndim != 1 or len(self.coefficients) > len(MODE_NAMES):
            raise ValueError("invalid correction coefficient vector")
        if not np.isfinite(self.coefficients).all():
            raise ValueError("correction coefficients must be finite")

    def __call__(self, points, time):
        value = self.candidate.velocity(points, time).copy()
        for j, coefficient in enumerate(self.coefficients):
            if coefficient != 0.0:
                value += coefficient * correction_mode(self.candidate, j, points, time)
        return value


def _sample(seed: int, n: int, halfwidth: float):
    if n <= 0 or halfwidth <= 0:
        raise ValueError("positive sample size and halfwidth required")
    return np.random.default_rng(seed).uniform(-halfwidth, halfwidth, size=(n, 3))


def _stack_momentum(candidate, force, points, times, cfg, coefficients=None):
    velocity = candidate.velocity if coefficients is None else CorrectedVelocity(candidate, coefficients)
    chunks = []
    for time in times:
        chunks.append(
            residual(
                velocity,
                candidate.pressure,
                force,
                points,
                time,
                nu=cfg.nu,
                step=cfg.derivative_step,
                time_bounds=(0.25, 0.75),
            )["momentum"]
        )
    return np.stack(chunks, axis=0)


def _cartesian_rms(momentum):
    array = np.asarray(momentum, dtype=float).reshape(-1, 3)
    return [float(np.sqrt(np.mean(array[:, j] ** 2))) for j in range(3)]


def _cylindrical_rms(momentum, points):
    """RMS of radial, azimuthal, axial momentum residual components."""
    points = np.asarray(points, dtype=float)
    x, y = points[:, 0], points[:, 1]
    radius = np.hypot(x, y)
    safe = np.where(radius > 1e-14, radius, 1.0)
    er = np.stack((x / safe, y / safe, np.zeros_like(x)), axis=-1)
    etheta = np.stack((-y / safe, x / safe, np.zeros_like(x)), axis=-1)
    radial = np.einsum("tnj,nj->tn", momentum, er)
    azimuthal = np.einsum("tnj,nj->tn", momentum, etheta)
    axial = momentum[..., 2]
    return {
        "radial": float(np.sqrt(np.mean(radial * radial))),
        "azimuthal": float(np.sqrt(np.mean(azimuthal * azimuthal))),
        "axial": float(np.sqrt(np.mean(axial * axial))),
    }


def _velocity_perturbation_ratio(candidate, coefficients, points, times):
    numerator = 0.0
    denominator = 0.0
    corrected = CorrectedVelocity(candidate, coefficients)
    for time in times:
        base = candidate.velocity(points, time)
        delta = corrected(points, time) - base
        numerator += float(np.sum(delta * delta))
        denominator += float(np.sum(base * base))
    return float(np.sqrt(numerator / denominator)) if denominator > 0 else None


def _fit_minimum_norm(J, r0):
    coefficients, *_ = np.linalg.lstsq(J, -r0, rcond=1e-10)
    return coefficients


def run_replication(candidate, force, cfg: ProbeConfig, train_seed: int, holdout_seed: int):
    train = _sample(train_seed, cfg.points_per_time, cfg.box_halfwidth)
    holdout = _sample(holdout_seed, cfg.points_per_time, cfg.box_halfwidth)
    train0 = _stack_momentum(candidate, force, train, cfg.train_times, cfg)
    holdout0 = _stack_momentum(candidate, force, holdout, cfg.holdout_times, cfg)
    r0 = train0.ravel()

    columns = []
    zero = np.zeros(len(MODE_NAMES))
    eps = cfg.jacobian_step
    for j in range(len(MODE_NAMES)):
        plus = zero.copy()
        minus = zero.copy()
        plus[j] = eps
        minus[j] = -eps
        rp = _stack_momentum(candidate, force, train, cfg.train_times, cfg, plus).ravel()
        rm = _stack_momentum(candidate, force, train, cfg.train_times, cfg, minus).ravel()
        columns.append((rp - rm) / (2.0 * eps))
    jacobian = np.column_stack(columns)

    baseline_holdout_rms = float(np.sqrt(np.mean(holdout0 * holdout0)))
    rows = []
    for k in BASIS_SIZES:
        coefficients = np.zeros(len(MODE_NAMES))
        if k:
            Jk = jacobian[:, :k]
            singular = np.linalg.svd(Jk, compute_uv=False)
            rank = int(np.linalg.matrix_rank(Jk, tol=singular[0] * 1e-10))
            condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else None
            coefficients[:k] = _fit_minimum_norm(Jk, r0)
        else:
            singular = np.array([])
            rank = 0
            condition = None

        train_residual = _stack_momentum(candidate, force, train, cfg.train_times, cfg, coefficients)
        heldout_residual = _stack_momentum(
            candidate, force, holdout, cfg.holdout_times, cfg, coefficients
        )
        train_rms = float(np.sqrt(np.mean(train_residual * train_residual)))
        heldout_rms = float(np.sqrt(np.mean(heldout_residual * heldout_residual)))
        rows.append(
            {
                "basis_dimension": k,
                "modes": list(MODE_NAMES[:k]),
                "coefficients": [float(v) for v in coefficients[:k]],
                "coefficient_l2": float(np.linalg.norm(coefficients)),
                "linearized_rank": rank,
                "jacobian_condition_number": condition,
                "singular_values": [float(v) for v in singular],
                "train_actual_nonlinear_rms": train_rms,
                "heldout_actual_nonlinear_rms": heldout_rms,
                "heldout_ratio_to_baseline": heldout_rms / baseline_holdout_rms,
                "heldout_cartesian_component_rms": _cartesian_rms(heldout_residual),
                "heldout_cylindrical_component_rms": _cylindrical_rms(heldout_residual, holdout),
                "heldout_velocity_perturbation_rms_ratio": _velocity_perturbation_ratio(
                    candidate, coefficients, holdout, cfg.holdout_times
                ),
            }
        )

    return {
        "train_seed": train_seed,
        "holdout_seed": holdout_seed,
        "baseline_train_actual_nonlinear_rms": float(np.sqrt(np.mean(train0 * train0))),
        "baseline_holdout_actual_nonlinear_rms": baseline_holdout_rms,
        "baseline_holdout_cylindrical_component_rms": _cylindrical_rms(holdout0, holdout),
        "rows": rows,
    }


def _summary(replications):
    by_dimension = {}
    for k in BASIS_SIZES:
        rows = [next(row for row in rep["rows"] if row["basis_dimension"] == k) for rep in replications]
        ratios = np.asarray([row["heldout_ratio_to_baseline"] for row in rows])
        conditions = [row["jacobian_condition_number"] for row in rows]
        by_dimension[str(k)] = {
            "heldout_ratio_median": float(np.median(ratios)),
            "heldout_ratio_min": float(np.min(ratios)),
            "heldout_ratio_max": float(np.max(ratios)),
            "full_column_rank_all_replications": all(row["linearized_rank"] == k for row in rows),
            "condition_number_range": None
            if k == 0
            else [float(min(conditions)), float(max(conditions))],
        }
    return by_dimension


def run_probe(candidate, force, cfg=ProbeConfig()):
    replications = [
        run_replication(candidate, force, cfg, train_seed, holdout_seed)
        for train_seed, holdout_seed in cfg.seed_pairs
    ]
    summary = _summary(replications)
    return {
        "schema_version": 1,
        "status": "capacity_probe_not_candidate_validation",
        "question": "Does nested higher-order compact divergence-free basis growth reduce independent residual, and which new direction matters?",
        "pressure_and_force": "frozen; no residual-defined/free forcing; pressure unchanged",
        "fit": "minimum-norm linearized least squares on training sample; actual nonlinear residual evaluated on disjoint held-out sample",
        "config": {
            "points_per_time": cfg.points_per_time,
            "train_times": list(cfg.train_times),
            "holdout_times": list(cfg.holdout_times),
            "seed_pairs": [list(pair) for pair in cfg.seed_pairs],
            "box_halfwidth": cfg.box_halfwidth,
            "derivative_step": cfg.derivative_step,
            "jacobian_step": cfg.jacobian_step,
            "nu": cfg.nu,
            "basis_sizes": list(BASIS_SIZES),
            "mode_order": list(MODE_NAMES),
        },
        "summary_by_basis_dimension": summary,
        "replications": replications,
        "interpretation_rule": "A lower held-out ratio is useful only as a local capacity signal; it is not a validation pass.",
        "limitations": [
            "optimized-v4 base candidate only",
            "local linearization for coefficient selection",
            "pressure and restricted force intentionally frozen",
            "finite random samples and finite-difference residual",
            "no full nonlinear reoptimization, energy renormalization, or project acceptance claim",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate", default="artifacts/constrained/optimized_v4/candidate.json"
    )
    parser.add_argument(
        "--training", default="artifacts/constrained/optimized_v4/training.json"
    )
    parser.add_argument(
        "--output", default="artifacts/constrained/basis_capacity/probe.json"
    )
    args = parser.parse_args()
    candidate = _load_candidate(args.candidate)
    force = _load_force(args.training)
    result = run_probe(candidate, force)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    compact = {
        k: {
            "median": row["heldout_ratio_median"],
            "range": [row["heldout_ratio_min"], row["heldout_ratio_max"]],
            "condition_range": row["condition_number_range"],
        }
        for k, row in result["summary_by_basis_dimension"].items()
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
