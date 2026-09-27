"""Bounded wave/momentum co-design with moment and cone-preserving inner solves.

The wave coefficients are the outer variables.  For each wave, all 180 tangent
controls are re-fit: the first 36 controls carry the four mean-moment and 81
cone constraints, while the remaining 144 controls are eliminated by a fixed
weighted least-squares projection.  The outer objective uses the envelope
gradient of this inner problem and a directional finite-difference check.

This is an assembled training diagnostic.  It does not establish a PDE
trajectory, finite-energy scale recursion, or acceptance of the construction.
"""

import hashlib
import json
import os
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import LinearConstraint, linprog, lsq_linear, minimize

from constrained_tangent_projection import ConstrainedTangent


ROOT = Path(__file__).resolve().parent


def decode(value):
    value = np.asarray(value, float)
    return value[..., 0] + 1j * value[..., 1]


def pack(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def real_matrix(matrix):
    matrix = np.asarray(matrix)
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class InnerProjection:
    """Moment/cone constrained tangent fit with harmonic controls eliminated."""

    def __init__(self, cache, moment, cones, seed_controls, rcond=1.0e-8,
                 safety_margin=1.0e-4):
        self.weights = np.asarray(cache["weights"], float)
        self.sqrtw = np.repeat(np.sqrt(self.weights), 3)
        self.design = np.asarray(cache["tangent_design"], float)
        self.weighted_design = self.design * self.sqrtw[:, None]
        self.D0 = self.weighted_design[:, :36]
        self.D1 = self.weighted_design[:, 36:]
        self.moment = moment
        self.cones = cones
        self.E = np.asarray(moment["moment_rows"], float)
        self.C = np.asarray(cones["cone_control_rows"], float)
        if self.E.shape != (4, 180) or self.C.shape != (81, 180):
            raise ValueError(f"Unexpected constraint rows {self.E.shape} {self.C.shape}")
        if max(float(np.max(np.abs(self.E[:, 36:]))),
               float(np.max(np.abs(self.C[:, 36:])))) > 1.0e-12:
            raise ValueError("Expected only first 36 tangent controls to affect moments/cones")
        self.E0 = self.E[:, :36]
        self.C0 = self.C[:, :36]
        self.lower = np.asarray(cones["cone_lower"], float)
        self.baseline = np.asarray(cones["cone_baseline"], float)
        self.wave_forms = np.asarray(cones["cone_wave_forms"], float)
        self.moment_baseline = np.asarray(moment["baseline_moments"], float)
        self.moment_forms = np.asarray(moment["wave_moment_forms"], float)
        self.safety_margin = float(safety_margin)

        # Project the 144 unconstrained controls out of the weighted residual.
        # Keeping only the numerical column space avoids allocating a 27k by
        # 27k complement basis.
        # Normalize harmonic columns before rank selection.  The pressure and
        # velocity tangent columns have different physical units; selecting a
        # rank from the unscaled matrix could discard a valid pressure block.
        self.harmonic_column_scales = np.maximum(
            np.linalg.norm(self.D1, axis=0), 1.0e-30)
        u, s, vh = np.linalg.svd(
            self.D1 / self.harmonic_column_scales[None, :], full_matrices=False)
        keep = s > float(rcond) * max(float(s[0]), 1.0e-300)
        self.U1 = u[:, keep]
        self.s1 = s[keep]
        self.Vh1 = vh[keep]
        self.D0p = self.D0 - self.U1 @ (self.U1.T @ self.D0)
        self.rank_harmonic = int(np.sum(keep))
        self.harmonic_rank_threshold = float(rcond) * max(float(s[0]), 1.0e-300)
        self.rcond = float(rcond)

        # A scaled nullspace parameter makes the mode-0 inner SLSQP well
        # conditioned while retaining exact physical control coefficients.
        tangent = ConstrainedTangent(self.D0, self.E0, rcond=rcond)
        self.tangent = tangent
        self.control_scales = tangent.scales
        self.particular_map = tangent.particular / self.control_scales[:, None]
        self.null_map = tangent.nullspace / self.control_scales[:, None]
        self.null_rank = self.null_map.shape[1]

        self.seed_controls = np.asarray(seed_controls, float)
        if self.seed_controls.shape != (180,):
            raise ValueError("Expected 180 saved tangent controls")
        self.last_z = None
        self.last_x = None
        self.last_result = None

    def wave_target(self, x_original):
        x_original = np.asarray(x_original, float)
        target = -self.moment_baseline - np.einsum(
            "i,kij,j->k", x_original, self.moment_forms, x_original)
        target_jac = -2.0 * np.einsum("kij,j->ki", self.moment_forms, x_original)
        return target, target_jac

    def wave_cone(self, x_original):
        x_original = np.asarray(x_original, float)
        return np.einsum("i,kij,j->k", x_original, self.wave_forms, x_original)

    def _maps(self, x_original):
        target, target_jac_original = self.wave_target(x_original)
        ypart = self.particular_map @ target
        cone_quad = self.wave_cone(x_original)
        rhs = (self.lower + self.safety_margin - self.baseline - cone_quad
               - self.C0 @ ypart)
        return target, target_jac_original, ypart, cone_quad, rhs

    def _seed_z(self, ypart):
        # nullspace columns are orthonormal in the scaled control coordinates.
        scaled = self.control_scales * (self.seed_controls[:36] - ypart)
        return self.tangent.nullspace.T @ scaled

    def _phase(self, A, rhs, z_guess):
        norms = np.linalg.norm(A, axis=1)
        responsive = norms > max(float(np.max(norms)) * 1.0e-10, 1.0e-14)
        if np.any(rhs[~responsive] > 1.0e-8):
            return None, dict(success=False, reason="fixed_cone_row_violation",
                               fixed_rows=np.flatnonzero(~responsive).tolist(),
                               fixed_rhs=rhs[~responsive].tolist(),
                               responsive_rows=np.flatnonzero(responsive).tolist())
        if np.all(A[responsive] @ z_guess >= rhs[responsive] - 2.0e-8):
            return z_guess.copy(), dict(success=True, method="warm_start", slack=0.0)
        An = A[responsive] / norms[responsive, None]
        bn = rhs[responsive] / norms[responsive]
        n = A.shape[1]
        phase = linprog(
            np.r_[np.zeros(n), 1.0],
            A_ub=np.column_stack((-An, -np.ones(len(bn)))),
            b_ub=-bn,
            bounds=[(None, None)] * n + [(0.0, None)], method="highs")
        if not phase.success:
            return None, dict(success=False, method="linprog", message=str(phase.message))
        return phase.x[:-1], dict(success=True, method="linprog",
                                   slack=float(phase.x[-1]))

    def evaluate(self, x_original, rw, jw_original, warm=True, maxiter=160):
        """Return the inner fit, envelope gradient, and feasibility evidence."""
        x_original = np.asarray(x_original, float)
        rw = np.asarray(rw, float)
        jw_original = np.asarray(jw_original, float)
        target, target_jac_original, ypart, cone_quad, rhs = self._maps(x_original)
        Acone = self.C0 @ self.null_map
        z_seed = self._seed_z(ypart)
        if warm and self.last_z is not None:
            z_start = self.last_z.copy()
        else:
            z_start = z_seed
        # _maps already subtracts C0 ypart, so rhs is the null-parameter
        # right hand side for C0 N z directly.
        bcone = rhs
        z_start, phase = self._phase(Acone, bcone, z_start)
        if z_start is None:
            return dict(success=False, phase=phase, objective=np.inf,
                        gradient=np.zeros(jw_original.shape[1]))

        p = rw - self.U1 @ (self.U1.T @ rw)
        Dp = self.D0p
        # Input derivative is with respect to original 54-vector coordinates.
        jp = jw_original - self.U1 @ (self.U1.T @ jw_original)
        yp = Dp @ ypart
        residual0 = p + yp
        # A scaled objective makes SLSQP termination independent of the huge
        # physical residual magnitude. The outer caller applies its own scale.
        obj_scale = max(float(np.linalg.norm(residual0)), 1.0)
        obj_scale = obj_scale * obj_scale

        # Physical null parameters are of order 10^6 because the moment rows
        # and tangent columns use different units. Optimize a unit-scale
        # variable and row-normalize the cone inequalities so neighboring
        # waves have a meaningful inner stationarity test.
        zscale = max(float(np.linalg.norm(z_start)), 1.0)
        Nscaled = self.null_map * zscale
        DpN = Dp @ Nscaled
        Aphysical = Acone * zscale
        row_norm = np.linalg.norm(Aphysical, axis=1)
        responsive = row_norm > max(float(np.max(row_norm)) * 1.0e-10, 1.0e-14)
        if np.any(~responsive) and np.any(bcone[~responsive] > 1.0e-8):
            return dict(success=False,
                        phase=dict(success=False,
                                   reason="fixed_cone_row_violation_after_scaling"),
                        objective=np.inf,
                        gradient=np.zeros(jw_original.shape[1]))
        Aconstraint = Aphysical[responsive] / row_norm[responsive, None]
        bconstraint = bcone[responsive] / row_norm[responsive]
        responsive_index = np.flatnonzero(responsive)
        w_start = z_start / zscale

        def objective(z):
            svec = residual0 + DpN @ z
            return 0.5 * float(svec @ svec) / obj_scale

        def objective_jac(z):
            svec = residual0 + DpN @ z
            return (DpN.T @ svec) / obj_scale

        constraints = [LinearConstraint(Aconstraint, bconstraint, np.inf)]
        fit = minimize(objective, w_start, jac=objective_jac, method="SLSQP",
                       constraints=constraints,
                       options=dict(maxiter=int(maxiter), ftol=1.0e-11, disp=False))
        w = np.asarray(fit.x if fit.x is not None else w_start, float)

        # Polish the convex inner QP with an active-set KKT solve.  This is
        # useful when SLSQP stops at a feasible point whose physical control
        # gradient is still large because the cone rows differ by units.
        H = (DpN.T @ DpN) / obj_scale
        g = (DpN.T @ residual0) / obj_scale
        # Select active rows in physical units. A normalized slack can be
        # misleading for nearly fixed rows with tiny control response.
        physical_slack = Aphysical @ w - bcone
        active_global = np.flatnonzero(physical_slack < 2.0e-5)
        active_qp = set(np.flatnonzero(np.isin(responsive_index, active_global)).tolist())
        qp_success = False
        qp_iterations = 0
        qp_mu = np.zeros(0)
        for qp_iterations in range(1, 101):
            active_list = np.array(sorted(active_qp), dtype=int)
            Aa = Aconstraint[active_list] if len(active_list) else np.zeros((0, len(w)))
            if len(active_list):
                K = np.block([[H, Aa.T], [Aa, np.zeros((len(active_list), len(active_list)))]])
                rhs_kkt = np.r_[-g, bconstraint[active_list]]
            else:
                K = H
                rhs_kkt = -g
            solution = np.linalg.lstsq(K, rhs_kkt, rcond=1.0e-12)[0]
            w_qp = solution[:len(w)]
            nu = solution[len(w):] if len(active_list) else np.zeros(0)
            slack = Aconstraint @ w_qp - bconstraint
            physical_slack = Aphysical @ w_qp - bcone
            violated_global = np.flatnonzero(physical_slack < -2.0e-8)
            violated = np.flatnonzero(
                np.isin(responsive_index, violated_global))
            if len(violated):
                # Add the most violated inactive row.
                candidates = [i for i in violated.tolist() if i not in active_qp]
                if not candidates:
                    # A currently active row can miss equality by roundoff;
                    # use the most negative row as an additional active row
                    # only when it is not already represented.
                    break
                active_qp.add(min(candidates,
                                  key=lambda i: physical_slack[responsive_index[i]]))
                continue
            if len(active_list):
                bad = np.flatnonzero(nu > 2.0e-8)
                if len(bad):
                    # For A w >= b, the KKT multiplier is -nu and must be
                    # nonnegative; remove the most negative physical multiplier.
                    active_qp.remove(int(active_list[bad[np.argmax(nu[bad])]]))
                    continue
            w = w_qp
            qp_mu = -nu
            qp_success = True
            break
        z_candidate = zscale * w
        candidate_y = ypart + self.null_map @ z_candidate
        candidate_margins = (self.baseline + self.C0 @ candidate_y + cone_quad
                             - self.lower - self.safety_margin)
        if np.min(candidate_margins) >= -2.0e-7 and np.isfinite(objective(w)):
            z = z_candidate
            if qp_success:
                fit = type("PolishedResult", (), {
                    "success": True, "status": 0,
                    "message": "SLSQP plus active-set KKT polish",
                    "fun": objective(w)})()
        else:
            z = z_start
            fit = type("FallbackResult", (), {
                "success": False,
                "status": int(getattr(fit, "status", 1)),
                "message": "SLSQP infeasible; retained phase-I feasible point",
                "fun": objective(w_start)})()
        y0 = ypart + self.null_map @ z
        # Reconstruct the eliminated 144 controls and the full weighted residual.
        shifted = rw + self.D0 @ y0
        projection = self.U1.T @ shifted
        y1 = -(self.Vh1.T @ (projection / self.s1)) / self.harmonic_column_scales
        controls = np.r_[y0, y1]
        full_residual = shifted + self.D1 @ y1
        margins = self.baseline + self.C @ controls + cone_quad - self.lower - self.safety_margin
        moment_error = self.E @ controls - target
        # The physical objective uses the unscaled full residual.  For the
        # inner gradient below, p and jp are both in original coordinates.
        physical_obj = 0.5 * float(full_residual @ full_residual)
        g_y = Dp.T @ full_residual
        active = np.flatnonzero(margins < 2.0e-5)
        if len(active):
            station = np.column_stack((self.E0.T, -self.C0[active].T))
            lower = np.r_[np.full(4, -np.inf), np.zeros(len(active))]
            upper = np.full(4 + len(active), np.inf)
            multipliers = lsq_linear(station, -g_y, bounds=(lower, upper),
                                     max_iter=200, lsmr_tol="auto")
            lam = multipliers.x[:4]
            mu = multipliers.x[4:]
            stationarity = float(np.linalg.norm(station @ multipliers.x + g_y))
        else:
            lam = np.linalg.lstsq(self.E0.T, -g_y, rcond=None)[0]
            mu = np.zeros(0)
            stationarity = float(np.linalg.norm(self.E0.T @ lam + g_y))
        # x derivatives are returned in original coefficient coordinates;
        # callers multiply by their normalized-to-original transform.
        # Use full-control KKT multipliers below, so the direct derivative is
        # taken at fixed physical y0.  The moment target derivative is then
        # supplied exactly once by -lambda.T @ target_jac_original.
        direct = jp
        cone_jac_original = 2.0 * np.einsum("kij,j->ki", self.wave_forms, x_original)
        grad_original = direct.T @ full_residual - target_jac_original.T @ lam
        if len(active):
            # d = lower + margin - baseline - cone_quad - C ypart, so the
            # envelope term is mu^T d_x = -mu^T cone_quad_x.
            grad_original -= cone_jac_original[active].T @ mu
        self.last_z = z.copy()
        self.last_x = x_original.copy()
        self.last_result = dict(
            success=bool(fit.success) and np.all(margins >= -2.0e-7),
            optimizer_success=bool(fit.success), optimizer_status=int(fit.status),
            optimizer_message=str(fit.message), phase=phase,
            objective=physical_obj, objective_scaled=float(fit.fun),
            gradient_original=grad_original,
            controls=controls, mode0_controls=y0, harmonic_controls=y1,
            residual=full_residual, moments=moment_error, margins=margins,
            active_cones=active, lambda_moments=lam, mu_cones=mu,
            kkt_stationarity=stationarity, cone_quad=cone_quad,
            target=target, target_jac_original=target_jac_original,
            ypart=ypart, null_parameter=z,
            harmonic_rank=self.rank_harmonic,
            harmonic_column_scale_min=float(np.min(self.harmonic_column_scales)),
            harmonic_column_scale_max=float(np.max(self.harmonic_column_scales)),
            active_set_qp_success=bool(qp_success),
            active_set_qp_iterations=int(qp_iterations),
            active_set_qp_count=int(len(active_qp)),
        )
        return self.last_result


class WaveProblem:
    """Outer source constraints and exact assembled momentum replay."""

    def __init__(self, source, seed, cache, moment, cones):
        self.source = source
        self.problem = source["problem"]
        self.cache = cache
        self.transform_source = decode(self.problem["whitening"])
        seed_white = decode(seed["selected"]["coefficients_whitened"])
        self.variable_scale = float(np.linalg.norm(seed_white))
        self.transform = real_matrix(self.transform_source) * self.variable_scale
        self.x0 = np.r_[seed_white.real, seed_white.imag] / self.variable_scale
        self.n = len(seed_white)
        self.x0_original = self.transform @ self.x0
        self.A = np.asarray(cache["A"], float)
        self.B = np.asarray(cache["B"], float)
        self.L = np.asarray(cache["L"], float)
        self.R0 = np.asarray(cache["R0"], float)
        self.sqrtw = np.repeat(np.sqrt(np.asarray(cache["weights"], float)), 3)
        self.flux_target = np.asarray(self.problem["target"], float).reshape(-1)
        self.flux_scales = np.maximum(np.abs(self.flux_target), 1.0)
        flux = decode(self.problem["flux_forms_whitened"]).reshape(-1, self.n, self.n)
        self.flux_forms = np.array([real_matrix(h) * self.variable_scale ** 2 for h in flux])
        self.mass = real_matrix(decode(self.problem["mass_whitened"])) * self.variable_scale ** 2
        self.growth = real_matrix(decode(self.problem["growth_whitened"])) * self.variable_scale ** 2
        self.floor = float(source["energy_threshold_lambda"])
        self.margin_matrix = self.growth - self.floor * self.mass
        self.margin_scale = max(abs(float(self.x0 @ self.growth @ self.x0)), 1.0)
        self.bound = float(self.problem["coefficient_bound"]) / self.variable_scale
        self.inner = InnerProjection(cache, moment, cones,
                                     np.asarray(seed["selected"]["tangent_coefficients"]),
                                     rcond=1.0e-8, safety_margin=1.0e-4)
        self.last = None
        self.objective_scale = 1.0

    def original(self, x):
        return self.transform @ np.asarray(x, float)

    def wave_residual(self, x):
        x_original = self.original(x)
        wave = self.A @ x_original
        gradient = self.B @ x_original
        residual = self.R0 + self.L @ x_original + np.einsum("nij,nj->ni", gradient, wave)
        derivative = self.L + np.einsum("nijq,nj->niq", self.B, wave)
        derivative += np.einsum("nij,njq->niq", gradient, self.A)
        rw = residual.reshape(-1) * self.sqrtw
        jw_original = derivative.reshape(-1, len(x_original)) * self.sqrtw[:, None]
        return x_original, rw, jw_original, residual, wave, gradient

    def source_equalities(self, x):
        x = np.asarray(x, float)
        return (np.einsum("i,kij,j->k", x, self.flux_forms, x) - self.flux_target) / self.flux_scales

    def source_equalities_jac(self, x):
        return 2.0 * np.einsum("kij,j->ki", self.flux_forms, x) / self.flux_scales[:, None]

    def source_inequalities(self, x):
        x = np.asarray(x, float)
        return np.array([x @ self.margin_matrix @ x / self.margin_scale,
                         1.0 - x @ x / self.bound ** 2])

    def source_inequalities_jac(self, x):
        x = np.asarray(x, float)
        return np.array([2.0 * self.margin_matrix @ x / self.margin_scale,
                         -2.0 * x / self.bound ** 2])

    def evaluate(self, x, warm=True, maxiter=160):
        x = np.asarray(x, float)
        if self.last is not None and np.array_equal(x, self.last["x"]):
            return self.last
        x_original, rw, jw_original, residual, wave, gradient = self.wave_residual(x)
        inner = self.inner.evaluate(x_original, rw, jw_original, warm=warm,
                                    maxiter=maxiter)
        if not inner.get("success", False):
            result = dict(x=x.copy(), x_original=x_original,
                          success=False, objective=np.inf, inner=inner,
                          flux=self.source_equalities(x),
                          inequalities=self.source_inequalities(x))
            self.last = result
            return result
        grad_original = np.asarray(inner["gradient_original"])
        grad = self.transform.T @ grad_original / self.objective_scale
        result = dict(x=x.copy(), x_original=x_original, success=True,
                      objective=float(inner["objective"] / self.objective_scale),
                      gradient=grad, inner=inner,
                      flux=self.source_equalities(x),
                      inequalities=self.source_inequalities(x),
                      residual=residual, wave=wave, wave_gradient=gradient)
        self.last = result
        return result


def _replay(problem, result, moment, cones):
    """Independent direct recomputation of all assembled constraints."""
    x = result["x"]
    x_original = result["x_original"]
    controls = result["inner"]["controls"]
    target = -np.asarray(moment["baseline_moments"]) - np.einsum(
        "i,kij,j->k", x_original, np.asarray(moment["wave_moment_forms"]), x_original)
    moment_error = np.asarray(moment["moment_rows"]) @ controls - target
    cone_margins = (np.asarray(cones["cone_baseline"])
                    + np.asarray(cones["cone_control_rows"]) @ controls
                    + np.einsum("i,kij,j->k", x_original,
                                np.asarray(cones["cone_wave_forms"]), x_original)
                    - np.asarray(cones["cone_lower"]))
    source_flux = problem.source_equalities(x)
    source_ineq = problem.source_inequalities(x)
    return dict(
        moment_max_abs=float(np.max(np.abs(moment_error))),
        moment_error=moment_error.tolist(),
        cone_min_margin=float(np.min(cone_margins)),
        cone_inequality_pass_count=int(np.sum(cone_margins >= 0.0)),
        cone_location_pass_count=int(np.sum(np.all(cone_margins.reshape(-1, 3) >= 0.0, axis=1))),
        source_flux_max_normalized=float(np.max(np.abs(source_flux))),
        source_flux_normalized=source_flux.tolist(),
        source_inequality_min=float(np.min(source_ineq)),
        source_inequalities=source_ineq.tolist(),
        coefficients_original=pack(decode(problem.problem["whitening"]) @
                                   (problem.variable_scale *
                                    (x[:problem.n] + 1j * x[problem.n:]))),
        coefficients_whitened=pack(problem.variable_scale *
                                   (x[:problem.n] + 1j * x[problem.n:])),
    )


def run(max_outer=10, output_path=None):
    started = time.perf_counter()
    paths = {name: ROOT / filename for name, filename in dict(
        source="wave_stress_growth_codesign.json",
        seed="wave_moment_cone_tangent.json",
        moment="wave_dynamics_mean_compatibility.json",
        cones="wave_mean_cone_projection.json").items()}
    raw = {name: path.read_bytes() for name, path in paths.items()}
    sources = {name: json.loads(value) for name, value in raw.items()}
    if any(report.get("status", "completed") not in ("completed", "selected")
           for report in sources.values()):
        raise ValueError("Source reports are not complete")
    with np.load(ROOT / "wave_momentum_projection.npz", allow_pickle=False) as loaded:
        cache = {key: loaded[key] for key in loaded.files}
    problem = WaveProblem(sources["source"], sources["seed"], cache,
                          sources["moment"]["reusable_moment_linearization"],
                          sources["cones"])
    output = Path(output_path) if output_path else ROOT / "wave_cone_codesign.json"
    report = dict(
        status="optimizing", accepted=False, pde_validated=False,
        constraints_maintained=False, scale_recursion_established=False,
        source=paths["source"].name,
        seed=paths["seed"].name, moment_source=paths["moment"].name,
        cone_source=paths["cones"].name,
        source_sha256=_sha(paths["source"]), seed_sha256=_sha(paths["seed"]),
        moment_sha256=_sha(paths["moment"]), cone_sha256=_sha(paths["cones"]),
        projection_cache_sha256=_sha(ROOT / "wave_momentum_projection.npz"),
        outer_dimension=int(len(problem.x0)), outer_maxiter=int(max_outer),
        variable_scale=problem.variable_scale,
        source_flux_count=int(len(problem.flux_target)), cone_count=81,
        moment_count=4, tangent_control_count=180,
        harmonic_elimination_rank=int(problem.inner.rank_harmonic),
        harmonic_elimination_columns=144,
        scope="Assembled wave/momentum co-design only; no PDE or recursive acceptance.",
        iterations=[])

    # Establish the normalization after the first actual inner solve.
    initial = problem.evaluate(problem.x0, warm=False, maxiter=220)
    if not initial.get("success", False):
        report.update(status="inner_initial_failure", reason=initial["inner"].get("phase"))
        report["elapsed_seconds"] = time.perf_counter() - started
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report
    problem.objective_scale = max(float(initial["inner"]["objective"]), 1.0)
    # Re-evaluate once so the stored gradient carries the final normalization.
    problem.last = None
    initial = problem.evaluate(problem.x0, warm=True, maxiter=220)

    # Directional derivative check at the fixed-wave seed. Two fresh inner
    # solves validate the envelope term including the moving moment/cone RHS.
    direction = np.linspace(1.0, 2.0, len(problem.x0))
    direction /= np.linalg.norm(direction)
    step = 2.0e-5
    plus = problem.evaluate(problem.x0 + step * direction, warm=False, maxiter=180)
    minus = problem.evaluate(problem.x0 - step * direction, warm=False, maxiter=180)
    fd = ((plus["objective"] - minus["objective"]) / (2.0 * step)
          if plus.get("success") and minus.get("success") else None)
    analytic = float(initial["gradient"] @ direction)
    report["objective_gradient_check"] = dict(
        method="central finite difference with independent inner re-solves",
        step=step, analytic_directional=analytic,
        finite_difference_directional=fd,
        absolute_error=(abs(analytic - fd) if fd is not None else None),
        relative_error=(abs(analytic - fd) / max(abs(fd), 1.0e-12)
                        if fd is not None else None),
        plus_inner_success=bool(plus.get("success", False)),
        minus_inner_success=bool(minus.get("success", False)),
        seed_kkt_stationarity=float(initial["inner"]["kkt_stationarity"]),
        seed_active_cones=np.asarray(initial["inner"]["active_cones"], int).tolist())
    problem.last = None

    best = problem.x0.copy()
    best_result = problem.evaluate(best, warm=False, maxiter=220)
    best_value = best_result["objective"]
    report["initial"] = dict(objective=best_value,
                              inner_optimizer_success=best_result["inner"]["optimizer_success"],
                              inner_kkt_stationarity=best_result["inner"]["kkt_stationarity"],
                              inner_active_cones=np.asarray(best_result["inner"]["active_cones"], int).tolist(),
                              assembled=_replay(problem, best_result, sources["moment"]["reusable_moment_linearization"], sources["cones"]))

    def objective(x):
        result = problem.evaluate(x, warm=True, maxiter=160)
        if not result.get("success", False):
            return 1.0e12
        return result["objective"]

    def objective_jac(x):
        result = problem.evaluate(x, warm=True, maxiter=160)
        if not result.get("success", False):
            return np.zeros_like(x)
        return result["gradient"]

    def callback(x):
        nonlocal best, best_result, best_value
        result = problem.evaluate(x, warm=True, maxiter=160)
        if result.get("success", False):
            assembled = _replay(problem, result, sources["moment"]["reusable_moment_linearization"], sources["cones"])
            feasible = (assembled["moment_max_abs"] <= 2.0e-6 and
                        assembled["cone_min_margin"] >= -2.0e-7 and
                        assembled["source_flux_max_normalized"] <= 2.0e-6 and
                        assembled["source_inequality_min"] >= -2.0e-7)
            if feasible and result["objective"] < best_value:
                best, best_result, best_value = np.asarray(x).copy(), result, result["objective"]
            report["iterations"].append(dict(
                iteration=len(report["iterations"]) + 1,
                objective=float(result["objective"]),
                feasible=bool(feasible), assembled=assembled,
                inner_optimizer_success=bool(result["inner"]["optimizer_success"]),
                inner_kkt_stationarity=float(result["inner"]["kkt_stationarity"]),
                inner_active_cones=np.asarray(result["inner"]["active_cones"], int).tolist()))
        else:
            report["iterations"].append(dict(iteration=len(report["iterations"]) + 1,
                                              objective=None, feasible=False,
                                              inner_failure=result["inner"].get("phase")))

    constraints = [
        {"type": "eq", "fun": problem.source_equalities,
         "jac": problem.source_equalities_jac},
        {"type": "ineq", "fun": problem.source_inequalities,
         "jac": problem.source_inequalities_jac},
    ]
    outer = minimize(objective, problem.x0, jac=objective_jac,
                     method="SLSQP", constraints=constraints,
                     callback=callback,
                     options=dict(maxiter=int(max_outer), ftol=1.0e-9, disp=False))
    candidate = np.asarray(outer.x if outer.x is not None else problem.x0, float)
    candidate_result = problem.evaluate(candidate, warm=False, maxiter=240)
    if candidate_result.get("success", False):
        candidate_assembled = _replay(problem, candidate_result,
                                       sources["moment"]["reusable_moment_linearization"],
                                       sources["cones"])
        candidate_feasible = (candidate_assembled["moment_max_abs"] <= 2.0e-6 and
                              candidate_assembled["cone_min_margin"] >= -2.0e-7 and
                              candidate_assembled["source_flux_max_normalized"] <= 2.0e-6 and
                              candidate_assembled["source_inequality_min"] >= -2.0e-7)
        if candidate_feasible and candidate_result["objective"] < best_value:
            best, best_result, best_value = candidate, candidate_result, candidate_result["objective"]
    best_result = problem.evaluate(best, warm=False, maxiter=240)
    report.update(status="completed", optimizer_success=bool(outer.success),
                  optimizer_status=int(outer.status), optimizer_message=str(outer.message),
                  selected_objective=float(best_result["objective"]),
                  selected=_replay(problem, best_result,
                                   sources["moment"]["reusable_moment_linearization"], sources["cones"]),
                  elapsed_seconds=time.perf_counter() - started,
                  outer_final_source_flux_max=float(np.max(np.abs(problem.source_equalities(candidate)))),
                  outer_final_source_inequality_min=float(np.min(problem.source_inequalities(candidate))))
    # This experiment is diagnostic by construction; retain explicit false
    # acceptance fields even when every sampled assembled constraint passes.
    report["accepted"] = False
    report["pde_validated"] = False
    report["constraints_maintained"] = False
    report["scale_recursion_established"] = False
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "outer_success": report["optimizer_success"],
                      "objective": report["selected_objective"],
                      "assembled": report["selected"],
                      "gradient_check": report["objective_gradient_check"]}), flush=True)
    return report


if __name__ == "__main__":
    run()
