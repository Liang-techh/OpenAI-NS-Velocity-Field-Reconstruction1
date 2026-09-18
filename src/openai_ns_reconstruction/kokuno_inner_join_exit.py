"""Executable later half of the corrected Kokuno PA.16 inner-to-outer join.

Pinned provenance is KokunoYumeto/yang-mills-interacting-workbench at commit
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.
The source fixes X_i=110, X_R=110(CP_*)^10, x=X/X_R and prescribes the
angular transition, axial restoration on -8<log x<-7, then the PA.16 two-U /
three-E five-moment repair on -6<log x<-5.  The upstream ell_i, G_i, incoming
moment discrepancy and admissible lower-bound choice of T_sh are outputs of
earlier Appendix-B existence choices, so this component takes them explicitly
and never invents hidden values.

At the repository's existence-style P_*~3e13 the raw PA.16 Jacobian rows span
many orders of magnitude.  Raw float64 matrix-rank testing therefore falsely
reports rank loss even though an invertible row scaling gives full rank.  This
module uses row scaling only for numerical conditioning; the exact nonlinear
moment map and targets are unchanged.  A final row-scaled Newton refinement is
used to close the original unscaled five moments rather than accepting an
optimizer stopping criterion in scaled coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from .kokuno_five_moment_repair import KokunoFiveMomentRepair, FiveMomentSolveResult
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_reference_continuation import source_smooth_step

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-inner-join-exit-v3"
X_I = 110.0
LOG_RESTORE_START = -8.0
LOG_RESTORE_END = -7.0
LOG_REPAIR_START = -6.0
LOG_JOIN_EXIT = -5.0

_SOURCE_FORMULAS = {
    "scales": "X_i=110; X_R=110*(C*P_*)^10; x=X/X_R; y_i=log(X/X_i)",
    "angular_transition": "log E=-log C+y_i/10+(1-sigma(y_i/T_sh))*ell_i+sigma(y_i/T_sh)*log f; U=G_i",
    "ideal_pair": "U_0=4 eta; E_0=P_* f x^(1/10); f=(1+eta^2)^(-1)",
    "axial_restoration": "restore G_i to 4 eta with a fixed flat step on -8<log x<-7",
    "PA16_repair": "two U bumps and three E bumps on -6<log x<-5 restore M,I,J,S,C_p",
    "join_exit": "X_h=X_R*exp(-5); matched fields+five moments preserve the exterior",
}
_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_later_inner_join_formula_executable": True,
    "source_PA16_inverse_reused": True,
    "scale_conditioned_PA16_solve": True,
    "stable_log_x_discrepancy_propagation": True,
    "caller_supplied_T_sh": True,
    "source_T_sh_lower_bound_verified": False,
    "caller_supplied_upstream_boundary_data": True,
    "actual_upstream_appendix_B_boundary_data_bound": False,
    "actual_inner_join_target_from_full_appendix_B": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _sigma_prime(s: np.ndarray) -> np.ndarray:
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    mask = (s > 0.0) & (s < 1.0)
    if np.any(mask):
        q = s[mask]
        sigma = source_smooth_step(q)
        out[mask] = sigma * (1.0 - sigma) * (2.0 / q**3 + 2.0 / (1.0 - q) ** 3)
    return out


@dataclass(frozen=True)
class InnerJoinSolveResult:
    eta: float
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    pre_repair_scaled_discrepancy: tuple[float, float, float, float, float]
    repair: FiveMomentSolveResult
    exit_scaled_discrepancy: tuple[float, float, float, float, float]
    max_abs_exit_discrepancy: float
    scaled_jacobian_condition: float


@dataclass(frozen=True)
class KokunoInnerJoinExit:
    T_sh: float
    outer_schedule: KokunoOuterReservedPatchSchedule = field(default_factory=KokunoOuterReservedPatchSchedule)
    repair: KokunoFiveMomentRepair = field(default_factory=KokunoFiveMomentRepair)
    quadrature_points: int = 96
    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        T = float(self.T_sh)
        if not math.isfinite(T) or T <= 0.0:
            raise ValueError("T_sh must be positive and finite")
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if not isinstance(self.repair, KokunoFiveMomentRepair):
            raise TypeError("repair must be a KokunoFiveMomentRepair")
        if isinstance(self.quadrature_points, bool) or not isinstance(self.quadrature_points, (int, np.integer)):
            raise TypeError("quadrature_points must be an integer")
        n = int(self.quadrature_points)
        if not 32 <= n <= 256:
            raise ValueError("quadrature_points must lie in [32,256]")
        object.__setattr__(self, "T_sh", T)
        object.__setattr__(self, "quadrature_points", n)
        if not self.log_x_sep < LOG_RESTORE_START:
            raise ValueError("source join requires x_sep<exp(-8); increase C/P_* or use a shorter caller-supplied T_sh")
        nodes, weights = leggauss(n)
        nodes = np.asarray(nodes, dtype=float); nodes.setflags(write=False)
        weights = np.asarray(weights, dtype=float); weights.setflags(write=False)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def log_x_i(self) -> float:
        return -10.0 * (math.log(float(self.outer_schedule.C)) + float(self.outer_schedule.log_P_star))

    @property
    def log_x_sep(self) -> float:
        return self.log_x_i + self.T_sh

    @property
    def log_X_i(self) -> float:
        return math.log(X_I)

    @property
    def log_X_sep(self) -> float:
        return self.log_X_i + self.T_sh

    @property
    def log_X_h(self) -> float:
        return float(self.outer_schedule.log_X_R) + LOG_JOIN_EXIT

    @property
    def P_star(self) -> float:
        return math.exp(float(self.outer_schedule.log_P_star))

    @staticmethod
    def source_f(eta: float) -> float:
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        return 1.0 / (1.0 + eta * eta)

    def geometry_report(self) -> dict[str, Any]:
        return {"X_i": X_I, "log_X_i": self.log_X_i, "T_sh": self.T_sh,
                "log_X_sep": self.log_X_sep, "log_x_i": self.log_x_i,
                "log_x_sep": self.log_x_sep,
                "axial_restore_log_x": [LOG_RESTORE_START, LOG_RESTORE_END],
                "PA16_repair_log_x": [LOG_REPAIR_START, LOG_JOIN_EXIT],
                "log_X_h": self.log_X_h, "source_T_sh_lower_bound_verified": False}

    def _ideal(self, log_x: Any, eta: float) -> tuple[np.ndarray, np.ndarray]:
        y = np.asarray(log_x, dtype=float)
        f = self.source_f(eta)
        E = np.exp(float(self.outer_schedule.log_P_star) + math.log(f) + 0.1 * y)
        return np.full_like(y, 4.0 * float(eta), dtype=float), E

    def profile_values_scaled(self, x: Any, *, eta: float, ell_i: float, G_i: float,
                              coefficients: Any | None = None) -> dict[str, np.ndarray]:
        raw = _finite(x, "x"); shape = raw.shape
        values = np.atleast_1d(raw).astype(float, copy=False)
        if np.any(values <= 0.0): raise ValueError("x must be positive")
        eta = float(eta); f = self.source_f(eta); ell_i = float(ell_i); G_i = float(G_i)
        if not math.isfinite(ell_i) or not math.isfinite(G_i): raise ValueError("ell_i and G_i must be finite")
        y = np.log(values)
        if np.any(y < self.log_x_i - 2e-14) or np.any(y > LOG_JOIN_EXIT + 2e-14):
            raise ValueError("scaled x is outside the executable inner-join interval")
        U0, E0 = self._ideal(y, eta)
        U = U0.copy(); E = E0.copy(); U_x = np.zeros_like(values); E_x = 0.1 * E0 / values
        mask = y <= self.log_x_sep
        if np.any(mask):
            s = (y[mask] - self.log_x_i) / self.T_sh
            sigma = source_smooth_step(s); sp = _sigma_prime(s)
            mismatch = ell_i - math.log(f)
            local_E = E0[mask] * np.exp((1.0 - sigma) * mismatch)
            U[mask] = G_i; E[mask] = local_E
            E_x[mask] = local_E * (0.1 - sp * mismatch / self.T_sh) / values[mask]
        mask = (y > self.log_x_sep) & (y <= LOG_RESTORE_START)
        if np.any(mask): U[mask] = G_i
        mask = (y > LOG_RESTORE_START) & (y < LOG_RESTORE_END)
        if np.any(mask):
            s = y[mask] - LOG_RESTORE_START; sigma = source_smooth_step(s); sp = _sigma_prime(s)
            U[mask] = G_i + sigma * (4.0 * eta - G_i)
            U_x[mask] = sp * (4.0 * eta - G_i) / values[mask]
        mask = (y > LOG_REPAIR_START) & (y < LOG_JOIN_EXIT)
        if np.any(mask):
            if coefficients is None: raise ValueError("PA.16 coefficients are required inside the repair interval")
            repaired = self.repair.corrected_scaled_profiles(values[mask], eta, self.P_star, f, coefficients)
            U[mask] = repaired["U"]; E[mask] = repaired["E"]
            U_x[mask] = repaired["U_x"]; E_x[mask] = repaired["E_x"]
        result = {"U": U, "E": E, "H": np.sqrt(2.0 * values) * E,
                  "F_scaled": E / np.sqrt(2.0 * values), "U_x": U_x, "E_x": E_x}
        return {k: np.asarray(v).reshape(shape) for k, v in result.items()}

    @staticmethod
    def _density(x: np.ndarray, U: np.ndarray, E: np.ndarray,
                 U0: np.ndarray, E0: np.ndarray) -> np.ndarray:
        H = np.sqrt(2.0*x)*E; H0 = np.sqrt(2.0*x)*E0
        return np.stack([U-U0, H-H0, U*H-U0*H0,
                         U*U-U0*U0-0.5*(E*E-E0*E0), (E*E-E0*E0)/(2.0*x)], axis=0)

    def _integrate(self, lo: float, hi: float, *, eta: float, ell_i: float, G_i: float) -> np.ndarray:
        if hi <= lo: return np.zeros(5)
        y = 0.5*(hi-lo)*self._nodes + 0.5*(hi+lo); x = np.exp(y)
        p = self.profile_values_scaled(x, eta=eta, ell_i=ell_i, G_i=G_i)
        U0, E0 = self._ideal(y, eta)
        d = self._density(x, np.asarray(p["U"]), np.asarray(p["E"]), U0, E0)
        return 0.5*(hi-lo) * ((d*x[np.newaxis,:]) @ self._weights)

    def pre_repair_scaled_discrepancy(self, *, eta: float, ell_i: float, G_i: float,
                                      incoming_scaled_discrepancy: Any) -> np.ndarray:
        incoming = _finite(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,): raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        eta = float(eta); self.source_f(eta); ell_i = float(ell_i); G_i = float(G_i)
        delta = self._integrate(self.log_x_i, self.log_x_sep, eta=eta, ell_i=ell_i, G_i=G_i)
        a, b = math.exp(self.log_x_sep), math.exp(LOG_RESTORE_START); dU = G_i - 4.0*eta
        if b > a and dU != 0.0:
            delta[0] += dU*(b-a)
            h = math.sqrt(2.0)*self.P_star*self.source_f(eta)
            delta[2] += dU*h*(b**1.6-a**1.6)/1.6
            delta[3] += (G_i*G_i-(4.0*eta)**2)*(b-a)
        delta += self._integrate(LOG_RESTORE_START, LOG_RESTORE_END, eta=eta, ell_i=ell_i, G_i=G_i)
        return incoming + delta

    @staticmethod
    def _row_matrix(eta: float) -> np.ndarray:
        T = np.zeros((5,5)); T[0,0]=1; T[1,1]=-4*eta; T[1,2]=1; T[2,1]=1; T[3,0]=-8*eta; T[3,3]=1; T[4,4]=1
        return T

    def _conditioned_pa16_solve(self, target: np.ndarray, *, eta: float,
                                absolute_tolerance: float) -> tuple[FiveMomentSolveResult, float]:
        target = np.asarray(target, dtype=float); f = self.source_f(eta); T = self._row_matrix(eta); zero = np.zeros(5)
        target_t = T @ target
        J0 = T @ self.repair.moment_jacobian(zero, eta=eta, P_star=self.P_star, f_eta=f)
        rs = np.linalg.norm(J0, axis=1)
        if np.any(~np.isfinite(rs)) or np.any(rs <= 0): raise RuntimeError("PA.16 Jacobian has zero/nonfinite row")
        Js = J0 / rs[:,None]; sv = np.linalg.svd(Js, compute_uv=False)
        if sv[-1] <= 1e-12*sv[0]: raise RuntimeError("source five-moment linearized map lost rank after row scaling")
        condition = float(sv[0]/sv[-1])
        coeff = np.linalg.solve(Js, target_t/rs)
        coeff = np.clip(coeff, -0.9*self.repair.coefficient_limit, 0.9*self.repair.coefficient_limit)
        fit_scale = np.maximum(np.abs(target_t), 1e-10)
        def residual(c: np.ndarray) -> np.ndarray:
            return (T @ (self.repair.moment_increment(c, eta=eta, P_star=self.P_star, f_eta=f)-target))/fit_scale
        def jac(c: np.ndarray) -> np.ndarray:
            return (T @ self.repair.moment_jacobian(c, eta=eta, P_star=self.P_star, f_eta=f))/fit_scale[:,None]
        fit = least_squares(residual, coeff, jac=jac,
                            bounds=(-self.repair.coefficient_limit, self.repair.coefficient_limit),
                            ftol=1e-14, xtol=1e-14, gtol=1e-14, max_nfev=400)
        coeff = np.asarray(fit.x, dtype=float)
        # Newton refinement uses only invertible row scaling of each current
        # Jacobian.  It therefore seeks zero of the original five equations.
        for _ in range(10):
            achieved = self.repair.moment_increment(coeff, eta=eta, P_star=self.P_star, f_eta=f)
            error = achieved-target
            if float(np.max(np.abs(error))) <= absolute_tolerance: break
            J = T @ self.repair.moment_jacobian(coeff, eta=eta, P_star=self.P_star, f_eta=f)
            e = T @ error; scale = np.linalg.norm(J, axis=1)
            step = np.linalg.solve(J/scale[:,None], -e/scale)
            trial = coeff + step
            if np.any(np.abs(trial) >= self.repair.coefficient_limit): break
            coeff = trial
        achieved = self.repair.moment_increment(coeff, eta=eta, P_star=self.P_star, f_eta=f)
        error = achieved-target; max_abs = float(np.max(np.abs(error)))
        sample = np.exp(np.linspace(LOG_REPAIR_START, LOG_JOIN_EXIT, 257))
        if np.any(np.asarray(self.repair.corrected_scaled_profiles(sample, eta, self.P_star, f, coeff)["E"]) <= 0):
            raise RuntimeError("PA.16 conditioned solve violated E positivity")
        if max_abs > absolute_tolerance:
            raise RuntimeError(f"conditioned PA.16 solve did not close: max residual {max_abs:.3e}")
        receipt = FiveMomentSolveResult(coefficients=tuple(map(float,coeff)),
            target_scaled_moments=tuple(map(float,target)), achieved_scaled_moments=tuple(map(float,achieved)),
            residual_scaled_moments=tuple(map(float,error)), max_abs_residual=max_abs,
            nfev=int(fit.nfev), success=True)
        return receipt, condition

    def solve_at_eta(self, *, eta: float, ell_i: float, G_i: float,
                     incoming_scaled_discrepancy: Any, absolute_tolerance: float = 2e-11) -> InnerJoinSolveResult:
        tol = float(absolute_tolerance)
        if not math.isfinite(tol) or not 0 < tol <= 1e-6: raise ValueError("absolute_tolerance must lie in (0,1e-6]")
        incoming = _finite(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,): raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        pre = self.pre_repair_scaled_discrepancy(eta=eta, ell_i=ell_i, G_i=G_i, incoming_scaled_discrepancy=incoming)
        repair, condition = self._conditioned_pa16_solve(-pre, eta=float(eta), absolute_tolerance=tol)
        exit_error = pre + np.asarray(repair.achieved_scaled_moments)
        return InnerJoinSolveResult(float(eta), tuple(map(float,incoming)), tuple(map(float,pre)), repair,
            tuple(map(float,exit_error)), float(np.max(np.abs(exit_error))), condition)

    def to_payload(self) -> dict[str, Any]:
        payload = {"schema": SCHEMA,
            "source": {"repository": SOURCE_REPOSITORY, "commit": SOURCE_COMMIT, "path": SOURCE_PATH,
                       "corrected_release": CORRECTED_RELEASE, "corrected_release_date": CORRECTED_RELEASE_DATE},
            "source_formulas": dict(_SOURCE_FORMULAS),
            "parameters": {"T_sh": self.T_sh, "quadrature_points": self.quadrature_points, "X_i": X_I,
                           "outer_schedule": self.outer_schedule.to_payload(), "repair": self.repair.to_payload()},
            "geometry": self.geometry_report(), "truth_boundary": dict(_TRUTH_BOUNDARY)}
        payload["sha256"] = hashlib.sha256(_json(payload).encode()).hexdigest(); return payload

    @property
    def sha256(self) -> str: return str(self.to_payload()["sha256"])

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoInnerJoinExit":
        if not isinstance(payload, dict): raise ValueError("payload must be a mapping")
        claimed = payload.get("sha256"); unsigned = {k:v for k,v in payload.items() if k!="sha256"}
        if claimed != hashlib.sha256(_json(unsigned).encode()).hexdigest(): raise ValueError("payload SHA mismatch")
        expected_source = {"repository": SOURCE_REPOSITORY, "commit": SOURCE_COMMIT, "path": SOURCE_PATH,
                           "corrected_release": CORRECTED_RELEASE, "corrected_release_date": CORRECTED_RELEASE_DATE}
        if payload.get("schema") != SCHEMA or payload.get("source") != expected_source: raise ValueError("source/schema metadata mismatch")
        if payload.get("source_formulas") != _SOURCE_FORMULAS: raise ValueError("source formula metadata mismatch")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY: raise ValueError("truth-boundary metadata mismatch")
        p = payload.get("parameters", {})
        if p.get("X_i") != X_I: raise ValueError("X_i metadata mismatch")
        obj = cls(T_sh=p.get("T_sh"), quadrature_points=p.get("quadrature_points"),
                  outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(p.get("outer_schedule")),
                  repair=KokunoFiveMomentRepair.from_payload(p.get("repair")))
        if obj.to_payload() != payload: raise ValueError("payload does not replay the exact inner-join exit")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target=Path(path); target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_payload(),indent=2,sort_keys=True)+"\n",encoding="utf-8"); return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoInnerJoinExit":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
