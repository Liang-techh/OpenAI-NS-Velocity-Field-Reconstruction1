"""Bounded low-degree ST006 continuation with joint momentum and moments.

This is a research experiment only.  The current checkout intentionally does
not contain the historical ``root_st030`` source, so the script extracts the
exact pinned source blobs into a temporary directory at runtime.  Nothing from
that directory is retained in the repository.  The final report is written
next to this file as ``st006_low_degree_joint.json``.  The optimization starts
from the complete 9 by 9 by 8 ST006 candidate and changes only an explicitly
listed low-degree raw coefficient subset; all other baseline coefficients are
frozen.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).with_suffix(".json")
SOURCE_COMMIT = "749999cb093906db32b859eb983cf9e8387da67f"
BASELINE_REL = "artifacts/research/ST006/candidate.json"
SOURCE_FILES = (
    "experiments/root_st030/spacetime.py",
    "experiments/root_st030/hybrid_basis.py",
    "experiments/root_st030/gauss_newton.py",
    "experiments/root_st030/harmonic_fit.py",
    "experiments/root_st030/validate.py",
)


def git_blob(relative: str) -> bytes:
    """Read a pinned historical blob without changing the checkout."""
    return subprocess.check_output(
        ["git", "show", f"{SOURCE_COMMIT}:{relative}"], cwd=ROOT
    )


def materialize_source(tmp: Path) -> tuple[Path, Path]:
    """Materialize source and baseline into a disposable import directory."""
    for relative in SOURCE_FILES:
        (tmp / Path(relative).name).write_bytes(git_blob(relative))
    baseline = tmp / "ST006_candidate.json"
    baseline.write_bytes(git_blob(BASELINE_REL))
    return tmp, baseline


def import_historical(tmp: Path):
    """Import the pinned modules under their original module names."""
    for name in ("spacetime", "hybrid_basis", "gauss_newton", "harmonic_fit", "validate"):
        sys.modules.pop(name, None)
    sys.path.insert(0, str(tmp))
    spacetime = importlib.import_module("spacetime")
    hybrid_basis = importlib.import_module("hybrid_basis")
    harmonic_fit = importlib.import_module("harmonic_fit")
    validate = importlib.import_module("validate")
    return spacetime, hybrid_basis, harmonic_fit, validate


def low_degree_warm_start(source_family, source_raw, target_family):
    """Project the leading spatial/time coefficient block to the target basis.

    ``Tp``, ``Tw`` and ``Tq`` map transformed coefficients to the physical
    compact basis.  Truncating in that physical modal representation preserves
    the source field's low-degree block, then a single velocity rescaling makes
    the target's exact initial-energy normalization factor equal to one.
    """
    aa, bb, qq, force_coeff, _ = source_family.coefficients(source_raw)
    projected = []
    for coeff, source_transform, target_transform in (
        (aa, source_family.Tp, target_family.Tp),
        (bb, source_family.Tw, target_family.Tw),
        (qq, source_family.Tq, target_family.Tq),
    ):
        physical = (source_transform @ coeff.reshape(source_family.ns, source_family.nt)).reshape(
            source_family.nr, source_family.nz, source_family.nt
        )
        cropped = np.zeros((target_family.nr, target_family.nz, target_family.nt))
        cropped[: target_family.nr, : target_family.nz, : target_family.nt] = physical[
            : target_family.nr, : target_family.nz, : target_family.nt
        ]
        transformed = np.linalg.solve(
            target_transform, cropped.reshape(target_family.ns, target_family.nt)
        )
        projected.append(transformed.ravel())
    velocity = np.concatenate(projected[:2])
    q0 = target_family.q0
    qmat = velocity.reshape(2 * target_family.ns, target_family.nt)
    norm = float(np.sum((qmat @ q0) ** 2))
    if not np.isfinite(norm) or norm <= 1e-20:
        raise ValueError("low-degree projection collapsed the initial velocity")
    velocity *= np.sqrt(2.0 / norm)
    return np.r_[velocity[: target_family.n], velocity[target_family.n :], projected[2], force_coeff]


class JointObjective:
    """Analytic PDE Jacobian plus finite-difference Jacobian for small extras.

    The historical optimizer used PyTorch for the extra residual block.  The
    execution image for this branch has no torch, so this class keeps the
    historical analytic PDE Jacobian and differentiates only the inexpensive
    energy/core/moment/torque block by central differences.  This leaves the
    full three-component momentum residual in the objective and keeps the
    numerical experiment reproducible with NumPy/SciPy alone.
    """

    def __init__(
        self,
        family,
        spacetime,
        harmonic_fit,
        seed=9172800,
        count=512,
        moment_weight=10.0,
        jac_indices=None,
    ):
        self.f = family
        self.st = spacetime
        self.moment_weight = float(moment_weight)
        self.jac_indices = None if jac_indices is None else np.asarray(jac_indices, dtype=int)
        self.history: list[dict] = []
        self.start = time.time()
        rng = np.random.default_rng(seed)
        n = int(count)
        s = rng.uniform(0.0, 4.0, n)
        z = rng.uniform(-2.0, 2.0, n)
        t = rng.uniform(0.25, 0.75, n)
        s[: n // 6] = rng.uniform(0.0, 0.5, n // 6) ** 2
        z[: n // 6] = rng.uniform(-0.5, 0.5, n // 6)
        t[n // 3 : n // 3 + n // 10] = 0.25
        t[n // 2 : n // 2 + n // 10] = 0.75
        self.D = family.bundle(s, z, t)
        self.s = s
        self.sqrt_s = np.sqrt(s)
        self.N = n
        points = np.column_stack((self.sqrt_s, np.zeros(n), z))
        self.FA = spacetime.force(points, t, 1, 0)
        self.FC = spacetime.force(points, t, 0, 1)

        # These are the same low-order admissibility diagnostics used by the
        # historical objective, plus all four harmonic moment degrees.
        tc = np.linspace(0.25, 0.75, 11)
        tau = 1.0 - tc
        sc = 0.01 * tau
        zc = 0.1 * tau**0.495
        core = family.bundle(sc, zc, tc)
        self.core = {k: core[k] for k in ("A", "B", "C", "s")}
        from numpy.polynomial import Chebyshev

        self.tc = tc
        self.CT = np.column_stack([Chebyshev.basis(k)(4 * tc - 2) for k in range(family.nt)])
        self.CTd = np.column_stack([4 * Chebyshev.basis(k).deriv()(4 * tc - 2) for k in range(family.nt)])
        self.scales = np.column_stack((np.sqrt(tau), tau**0.505, tau**0.505))
        self.q0 = family.q0
        self.GZ = family.Gz
        self.JB = family.Jb
        sq, zq, wq = spacetime.quad(96)
        cforce = float(
            (wq * sq)
            @ (spacetime.bump_derivatives(sq / 4.0)[0] * spacetime.bump_derivatives(zq * zq / 4.0)[0])
        )
        self.torque_t = cforce * spacetime.bump_derivatives((2 * tc - 1) ** 2)[0]
        seedvel = family.fields(
            family.initial(), np.column_stack((np.sqrt(sc), np.zeros(11), zc)), tc
        )[0] * self.scales
        self.floor = 0.5 * float(np.linalg.norm(seedvel[0]))
        self.Hp, self.Hb, self.harmonic_gram = harmonic_fit.harmonic_tensors(
            family, order=80, degrees=(2, 4, 6, 8)
        )
        self.cache_x: np.ndarray | None = None
        self.cache_r: np.ndarray | None = None
        self.cache_j: np.ndarray | None = None

    def extra(self, x: np.ndarray) -> np.ndarray:
        f = self.f
        n = f.n
        al = np.asarray(x[:n]).reshape(f.ns, f.nt)
        be = np.asarray(x[n : 2 * n]).reshape(f.ns, f.nt)
        initial = np.vstack((al, be)) @ self.q0
        den = float(initial @ initial)
        if den <= 1e-20:
            raise ValueError("collapsed initial velocity")
        scale = np.sqrt(2.0 / den)
        aa = scale * al.ravel()
        bb = scale * be.ravel()
        ap = scale * al @ self.CT.T
        bp = scale * be @ self.CT.T
        energy = 0.5 * (np.sum(ap * ap, axis=0) + np.sum(bp * bp, axis=0))
        cv = np.column_stack(
            (
                np.sqrt(self.core["s"]) * (self.core["A"] @ aa),
                np.sqrt(self.core["s"]) * (self.core["B"] @ bb),
                self.core["C"] @ aa,
            )
        ) * self.scales
        ref = cv[0]
        ref_norm = float(np.linalg.norm(ref))
        drift = np.sqrt(np.sum((cv - ref) ** 2, axis=1) + 1e-28) / (ref_norm + 1e-12)
        moments = np.einsum("at,kab,bt->kt", ap, self.Hp, ap) + np.einsum(
            "at,kab,bt->kt", bp, self.Hb, bp
        )
        D = self.moment_weight * moments.ravel() / np.sqrt(len(energy))
        torque = self.JB @ (scale * be @ self.CTd.T) + x[-1] * self.torque_t
        m = len(energy)
        return np.concatenate(
            (
                10.0 * np.maximum(0.1 - energy, 0.0) / np.sqrt(m),
                10.0 * np.maximum(energy - 10.0, 0.0) / np.sqrt(m),
                np.sqrt(10.0 / m) * np.maximum(drift - 0.045, 0.0),
                np.array([np.sqrt(10.0) * max(self.floor - ref_norm, 0.0) / self.floor]),
                np.sqrt(10.0 / m) * np.maximum(cv[:, 0] + 1e-6, 0.0),
                np.sqrt(10.0 / m) * np.maximum(-cv[:, 1] + 1e-6, 0.0),
                np.sqrt(10.0 / m) * np.maximum(-cv[:, 2] + 1e-6, 0.0),
                D,
                np.sqrt(10.0) * torque / np.sqrt(m * 32.0 * np.pi),
            )
        )

    def evaluate(self, x: np.ndarray) -> None:
        if self.cache_x is not None and np.array_equal(x, self.cache_x):
            return
        f = self.f
        D = self.D
        n = f.n
        N = self.N
        S = self.s
        r = self.sqrt_s
        aa, bb, qq, fc, amp = f.coefficients(x)
        vals = {
            k: D[k] @ (bb if k.startswith("B") else qq if k.startswith("Q") else aa)
            for k in D
            if k not in ("s", "z", "t", "Q")
        }
        A, B, C = vals["A"], vals["B"], vals["C"]
        As, Az = vals["As"], vals["Az"]
        Bs, Bz = vals["Bs"], vals["Bz"]
        Cs, Cz = vals["Cs"], vals["Cz"]
        Rr = r * (vals["At"] + A * A + 2 * S * A * As + C * Az - B * B + 2 * vals["Qs"] - self.st.NU * vals["AL"])
        Rr -= fc[0] * self.FA[:, 0] + fc[1] * self.FC[:, 0]
        Rt = r * (vals["Bt"] + 2 * A * B + 2 * S * A * Bs + C * Bz - self.st.NU * vals["BL"])
        Rt -= fc[0] * self.FA[:, 1] + fc[1] * self.FC[:, 1]
        Rz = vals["Ct"] + 2 * S * A * Cs + C * Cz + vals["Qz"] - self.st.NU * vals["CL"]
        Rz -= fc[0] * self.FA[:, 2] + fc[1] * self.FC[:, 2]
        J = np.zeros((3 * N, 3 * n + 2))
        J[:N, :n] = r[:, None] * (
            D["At"]
            + (2 * A + 2 * S * As)[:, None] * D["A"]
            + (2 * S * A)[:, None] * D["As"]
            + Az[:, None] * D["C"]
            + C[:, None] * D["Az"]
            - self.st.NU * D["AL"]
        )
        J[:N, n : 2 * n] = (-2 * r * B)[:, None] * D["B"]
        J[N : 2 * N, :n] = r[:, None] * ((2 * B + 2 * S * Bs)[:, None] * D["A"] + Bz[:, None] * D["C"])
        J[N : 2 * N, n : 2 * n] = r[:, None] * (
            D["Bt"] + (2 * A)[:, None] * D["B"] + (2 * S * A)[:, None] * D["Bs"] + C[:, None] * D["Bz"] - self.st.NU * D["BL"]
        )
        J[2 * N :, :n] = D["Ct"] + (2 * S * Cs)[:, None] * D["A"] + (2 * S * A)[:, None] * D["Cs"] + Cz[:, None] * D["C"] + C[:, None] * D["Cz"] - self.st.NU * D["CL"]
        J[2 * N :, n : 2 * n] = 0.0
        J[:N, 2 * n : 3 * n] = 2 * r[:, None] * D["Qs"]
        J[2 * N :, 2 * n : 3 * n] = D["Qz"]
        J[:N, -2] = -self.FA[:, 0]
        J[N : 2 * N, -2] = -self.FA[:, 1]
        J[2 * N :, -2] = -self.FA[:, 2]
        J[:N, -1] = -self.FC[:, 0]
        J[N : 2 * N, -1] = -self.FC[:, 1]
        J[2 * N :, -1] = -self.FC[:, 2]

        # Pull the PDE Jacobian back through the exact initial-energy map.
        q = x[: 2 * n]
        initial = q.reshape(2 * f.ns, f.nt) @ f.q0
        den = float(initial @ initial)
        v = (initial[:, None] * f.q0[None, :]).ravel()
        directional = J[:, : 2 * n] @ q
        J[:, : 2 * n] = amp * (J[:, : 2 * n] - directional[:, None] * (v / den)[None, :])
        R = np.concatenate((Rr, Rt, Rz)) / np.sqrt(N)
        J /= np.sqrt(N)
        extras = self.extra(x)

        # Extra residuals are only 100 or so rows; central differences avoid a
        # dependency on torch while retaining the exact PDE Jacobian above.
        diff_indices = np.arange(len(x), dtype=int) if self.jac_indices is None else self.jac_indices
        je = np.empty((len(extras), len(diff_indices)))
        for column, k in enumerate(diff_indices):
            step = 2e-6 * max(1.0, abs(float(x[k])))
            xp = x.copy(); xm = x.copy()
            xp[k] += step; xm[k] -= step
            je[:, column] = (self.extra(xp) - self.extra(xm)) / (2.0 * step)
        pde_jacobian = J if self.jac_indices is None else J[:, self.jac_indices]
        self.cache_x = x.copy()
        self.cache_r = np.r_[R, extras]
        self.cache_j = np.vstack((pde_jacobian, je))
        row = {
            "evaluation": len(self.history) + 1,
            "loss": float(self.cache_r @ self.cache_r),
            "pde_mse": float(R @ R),
            "extra_loss": float(extras @ extras),
            "elapsed_seconds": time.time() - self.start,
        }
        self.history.append(row)
        if len(self.history) == 1 or len(self.history) % 5 == 0:
            print(json.dumps(row), flush=True)

    def fun(self, x: np.ndarray) -> np.ndarray:
        self.evaluate(np.asarray(x, dtype=float))
        assert self.cache_r is not None
        return self.cache_r

    def jac(self, x: np.ndarray) -> np.ndarray:
        self.evaluate(np.asarray(x, dtype=float))
        assert self.cache_j is not None
        return self.cache_j


def low_degree_indices(family, spatial_degree=2, temporal_degree=2):
    """Select low-degree transformed coefficients in each field block."""
    selected = []
    per_field = []
    for block in range(3):
        field_spatial_degree = spatial_degree if block < 2 else 1
        offset = block * family.n
        field_indices = []
        for i in range(field_spatial_degree):
            for j in range(field_spatial_degree):
                for k in range(temporal_degree):
                    field_indices.append(offset + (i * family.nz + j) * family.nt + k)
        per_field.append(field_indices)
        selected.extend(field_indices)
    return np.asarray(selected, dtype=int), per_field


def bounded_fit(family, warm, spacetime, harmonic_fit, training_seed, training_points, max_nfev):
    """Fit only low-degree increments while retaining every ST006 coefficient."""
    baseline = np.asarray(warm, dtype=float).copy()
    selected, per_field = low_degree_indices(family)
    # Keep the perturbation small in raw coordinates.  The pressure cap is
    # larger because pressure enters linearly and is the least scaled block.
    caps = np.concatenate(
        (
            np.full(len(per_field[0]) + len(per_field[1]), 0.03),
            np.full(len(per_field[2]), 1.0),
        )
    )
    objective = JointObjective(
        family,
        spacetime,
        harmonic_fit,
        seed=training_seed,
        count=training_points,
        moment_weight=10.0,
        jac_indices=selected,
    )

    def expand(increment):
        full = baseline.copy()
        full[selected] += np.asarray(increment, dtype=float)
        return full

    def fun(increment):
        return objective.fun(expand(increment))

    def jac(increment):
        return objective.jac(expand(increment))

    zero = np.zeros(len(selected), dtype=float)
    initial_residual = fun(zero)
    start = time.time()
    result = least_squares(
        fun,
        zero,
        jac=jac,
        bounds=(-caps, caps),
        method="trf",
        tr_solver="lsmr",
        x_scale="jac",
        max_nfev=max_nfev,
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-8,
        tr_options={"maxiter": 80, "atol": 1e-5, "btol": 1e-5},
        verbose=0,
    )
    best = expand(result.x)
    return result, objective, baseline, best, selected, caps, initial_residual, time.time() - start


def compact_validation(validation: dict) -> dict:
    rows = validation["spatial_refinement"]
    fine = [row for row in rows if abs(row["space_step"] - 0.005) < 1e-15]
    return {
        "validation_seed": validation["validation_seed"],
        "uniform_box_points": validation["uniform_box_points"],
        "times": sorted({row["time"] for row in fine}),
        "fine_space_step": 0.005,
        "worst_momentum_max": max(row["momentum_max"] for row in fine),
        "worst_momentum_L2": max(row["momentum_L2"] for row in fine),
        "worst_divergence_max": max(row["divergence_max"] for row in fine),
        "worst_divergence_L2": max(row["divergence_L2"] for row in fine),
        "energy_at_order96": [row["energy"] for row in validation["energy_quadrature"] if row["order"] == 96],
        "energy_relative_change_48_to_96": validation["energy_relative_change_48_to_96"],
        "boundary_max": validation["boundary_max"],
        "derived_amplitude": validation["derived_amplitude"],
        "gates": validation["gates"],
        "all_numeric_gates_pass": validation["all_numeric_gates_pass"],
    }


def run(output: Path = OUT, max_nfev: int = 8, training_points: int = 256) -> dict:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="st006_low_degree_") as disposable:
        tmp = Path(disposable)
        materialize_source(tmp)
        spacetime, hybrid_basis, harmonic_fit, validate = import_historical(tmp)
        source_family, source_raw = spacetime.Family.load(tmp / "ST006_candidate.json")
        # The warm state is the complete published 9 by 9 by 8 candidate.  No
        # projection or modal truncation occurs here.
        target = source_family
        result, objective, baseline, best, selected, caps, initial_residual, fit_seconds = bounded_fit(
            target,
            source_raw,
            spacetime,
            harmonic_fit,
            training_seed=9172800,
            training_points=training_points,
            max_nfev=max_nfev,
        )
        best = np.asarray(best, dtype=float)
        baseline = np.asarray(baseline, dtype=float)
        # Replay the unchanged baseline and then validate the bounded child
        # through the same independent historical Cartesian FD routine.  Both
        # candidates remain disposable; only the report is retained.
        baseline_validation_path = tmp / "baseline_validation.json"
        baseline_validation = validate.validate(
            tmp / "ST006_candidate.json",
            baseline_validation_path,
            seed=9172801,
        )
        candidate_path = tmp / "best_candidate.json"
        target.save(
            best,
            candidate_path,
            metadata={
                "experiment": "ST006-low-degree-joint",
                "training_seed": 9172800,
                "training_points": training_points,
                "active_indices": selected.tolist(),
                "bounded_increment_max_abs": float(np.max(np.abs(best - baseline))),
                "pde_validated": False,
            },
        )
        validation_path = tmp / "independent_validation.json"
        validation = validate.validate(candidate_path, validation_path, seed=9172801)
        candidate_payload = json.loads(candidate_path.read_text())
        candidate_sha256 = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
        baseline_metrics = compact_validation(baseline_validation)
        candidate_metrics = compact_validation(validation)
        published_baseline = {
            "validation_seed": 9172801,
            "worst_momentum_max": 0.10822893050403214,
            "worst_momentum_L2": 0.10758432876171865,
        }
        baseline_replay_discrepancy = {
            key: float(baseline_metrics[key] - published_baseline[key])
            for key in ("worst_momentum_max", "worst_momentum_L2")
        }
        baseline_training_residual = objective.fun(baseline)
        best_training_residual = objective.fun(best)

        report = {
            "schema": "st006_low_degree_joint_experiment_v1",
            "experiment_id": "ST006-low-degree-joint",
            "status": "research_experiment_unaccepted",
            "source": {
                "baseline": BASELINE_REL,
                "baseline_commit": SOURCE_COMMIT,
                "baseline_basis": "hybrid_legendre7_inverse_even_v1",
                "historical_validation": "experiments/root_st030/validate.py",
            },
            "configuration": {
                "target_basis": {"nr": target.nr, "nz": target.nz, "nt": target.nt},
                "warm_start_basis_retained": True,
                "active_parameter_count": int(len(selected)),
                "active_indices": selected.tolist(),
                "training_seed": 9172800,
                "training_points": training_points,
                "validation_seed": 9172801,
                "validation_points": 4096,
                "validation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
                "times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
                "nu": 0.01,
                "objective": "all three cylindrical momentum components plus normalized harmonic weak moments degrees 2,4,6,8, energy/core/torque diagnostics",
                "harmonic_moment_weight": 10.0,
                "increment_caps": {"velocity_raw_abs": 0.03, "pressure_raw_abs": 1.0, "force_raw_abs": 0.0},
                "frozen_parameter_count": int(len(baseline) - len(selected)),
            },
            "optimizer": {
                "method": "bounded trust-region least_squares with analytic momentum Jacobian and central-difference diagnostic Jacobian",
                "success": bool(result.success),
                "message": str(result.message),
                "status_code": int(result.status),
                "nfev": int(result.nfev),
                "njev": int(result.njev) if result.njev is not None else None,
                "optimality": float(result.optimality),
                "fit_seconds": fit_seconds,
                "initial_joint_loss": float(initial_residual @ initial_residual),
                "best_joint_loss": float(2.0 * result.cost),
                "initial_pde_mse": float(baseline_training_residual[: 3 * training_points] @ baseline_training_residual[: 3 * training_points]),
                "best_pde_mse": float(best_training_residual[: 3 * training_points] @ best_training_residual[: 3 * training_points]),
                "max_abs_raw_increment": float(np.max(np.abs(best - baseline))),
                "max_abs_active_increment": float(np.max(np.abs(best[selected] - baseline[selected]))),
                "max_abs_velocity_increment": float(np.max(np.abs(best[selected[:16]] - baseline[selected[:16]]))),
                "max_abs_pressure_increment": float(np.max(np.abs(best[selected[16:]] - baseline[selected[16:]]))),
                "max_abs_force_increment": 0.0,
                "evaluations": objective.history,
            },
            "candidate": {
                "sha256": candidate_sha256,
                "schema": candidate_payload["schema"],
                "nr": candidate_payload["nr"],
                "nz": candidate_payload["nz"],
                "nt": candidate_payload["nt"],
                "basis_kind": candidate_payload.get("basis_kind"),
                "coefficients": candidate_payload["coefficients"],
                "metadata": candidate_payload["metadata"],
            },
            "baseline_replay": baseline_metrics,
            "published_baseline_headline": published_baseline,
            "baseline_replay_discrepancy": baseline_replay_discrepancy,
            "independent_validation": candidate_metrics,
            "claims": {
                "baseline_overwritten": False,
                "pde_validated": False,
                "accepted": False,
                "continuum_supremum_proved": False,
                "comparison_scope": "same six-time, 4096-point, seed9172801 Cartesian FD protocol for baseline replay and child; training used a bounded 256-point screening subset with 18 active coefficients and all other ST006 coefficients frozen",
            },
        }
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(output), "best_joint_loss": report["optimizer"]["best_joint_loss"], "validation": report["independent_validation"]}, indent=2), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--max-nfev", type=int, default=8)
    parser.add_argument("--training-points", type=int, default=256)
    args = parser.parse_args()
    run(args.out, args.max_nfev, args.training_points)
