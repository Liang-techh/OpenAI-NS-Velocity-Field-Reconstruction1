"""Fixed-slice positive E-stage capacity for the delayed ST073 mean.

Each row is an independent normalized eta slice.  The five E directions are
analytic axisymmetric swirl modes, so the corresponding three-dimensional
mode is exactly divergence free.  Their radial supports are disjoint and
start at X=1.04, which leaves the sampled inner cone through X=1.03 alone.

The optimization enforces the I and Cp equalities and a sampled relative E
floor.  Its S objective is an optimistic lower bound obtained by allowing an
arbitrary U correction on 1.04 <= X <= 2.98 while preserving target M and J.
The independent slice coefficients are deliberately not assembled into a
global physical field.
"""

from __future__ import annotations

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from azimuthal_capacity_optimize import INTERVALS as SOURCE_INTERVALS
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


TAU = 0.5 * 2 ** (-5.5)
TRAIN_ETAS = (0.2, 0.225, 0.25, 0.275, 0.3)
RELATIVE_E_FLOOR = 0.11
INNER_CONE_CUTOFF = 1.03
U_SUPPORT = (1.04, 2.98)

# The source first interval starts at 1.005.  Moving only that endpoint to
# 1.04 is required by the inner-cone preservation condition.  The remaining
# four source intervals are already outside the cutoff and remain unchanged.
RADIAL_INTERVALS = tuple(
    (max(INNER_CONE_CUTOFF + 0.01, float(lo)), float(hi))
    for lo, hi in SOURCE_INTERVALS
)

# This eta window is shared by the mode definitions only to make every mode
# nonzero on all training slices.  Coefficients below are still solved
# independently at each eta and are not a global lift.
MODE_ETA_INTERVAL = (0.15, 0.35)
QUADRATURE_ORDER = 96
START_COUNT = 6
MAX_ITERATIONS = 500


def quadrature(order=QUADRATURE_ORDER):
    """Piecewise Gauss rule resolving every support and U boundary."""
    edges = sorted(
        set(
            (
                0.0,
                3 / 32,
                1.0,
                INNER_CONE_CUTOFF,
                U_SUPPORT[0],
                3.0,
                3.5,
                *(edge for interval in RADIAL_INTERVALS for edge in interval),
            )
        )
    )
    nodes, weights = leggauss(order)
    X = np.concatenate(
        [(lo + hi) / 2 + (hi - lo) * nodes / 2
         for lo, hi in zip(edges[:-1], edges[1:])]
    )
    W = np.concatenate(
        [(hi - lo) * weights / 2
         for lo, hi in zip(edges[:-1], edges[1:])]
    )
    return X, W


def optimistic_continuum_slack(U, E, target, X, weights):
    """Return target S minus the unrestricted U minimum on U_SUPPORT.

    Outside U_SUPPORT, U is frozen at the current ST073 mean.  Inside it,
    the minimum of integral U^2 subject to M and J is the two-by-two Gram
    solve for the functions 1 and H=sqrt(2X) E.  This is an optimistic
    slice capacity calculation; it is evaluated by the displayed
    quadrature and imposes no smoothness or physical lifting conditions.
    """
    inside = (X >= U_SUPPORT[0]) & (X <= U_SUPPORT[1])
    outside = ~inside
    H = np.sqrt(2 * X) * E
    wi = weights[inside]
    hi = H[inside]
    gram = np.array(
        [[wi.sum(), wi @ hi],
         [wi @ hi, wi @ (hi * hi)]],
        dtype=float,
    )
    required = np.array(
        [
            target[0] - weights[outside] @ U[outside],
            target[2] - weights[outside] @ (U[outside] * H[outside]),
        ],
        dtype=float,
    )
    inside_minimum = float(required @ np.linalg.solve(gram, required))
    s_min = float(
        weights[outside] @ (U[outside] * U[outside])
        + inside_minimum
        - 0.5 * (weights @ (E * E))
    )
    return dict(
        S_min=s_min,
        S_slack=float(target[3] - s_min),
        target_S=float(target[3]),
        required_M=float(required[0]),
        required_J=float(required[1]),
        gram_condition=float(np.linalg.cond(gram)),
    )


def solve_eta(mean, target_field, base, X, weights, eta, rng):
    """Solve one independent eta slice and return a serializable row."""
    U, E = profile(mean, X, eta, TAU)
    target_U, target_E = profile(target_field, X, eta, TAU)
    target = moment_vector(target_U, target_E, X, weights)
    modes = [
        SimilaritySwirlMode(
            base,
            interval,
            eta_interval=MODE_ETA_INTERVAL,
        )
        for interval in RADIAL_INTERVALS
    ]
    mode_profiles = [mode_profile(mode, base, X, eta, TAU)
                     for mode in modes]
    # SimilaritySwirlMode has no axial component, hence the U columns are
    # exactly zero and only the E columns enter this stage.
    mode_U = np.column_stack([item[0] for item in mode_profiles])
    mode_E = np.column_stack([item[1] for item in mode_profiles])
    if float(np.max(np.abs(mode_U))) > 1e-12:
        raise RuntimeError("SimilaritySwirlMode unexpectedly changed U")

    def evaluate(coefficients):
        corrected_E = E + mode_E @ coefficients
        H = np.sqrt(2 * X) * corrected_E
        defects = np.array(
            [
                weights @ H - target[1],
                weights @ (corrected_E * corrected_E / (2 * X)) - target[4],
            ],
            dtype=float,
        )
        slack = optimistic_continuum_slack(
            U, corrected_E, target, X, weights
        )
        ratio = np.divide(
            corrected_E,
            target_E,
            out=np.full_like(corrected_E, np.inf),
            where=np.abs(target_E) > 1e-12,
        )
        return corrected_E, defects, slack, ratio

    def objective(coefficients):
        return -evaluate(coefficients)[2]["S_slack"]

    def equalities(coefficients):
        return evaluate(coefficients)[1]

    def positivity(coefficients):
        return evaluate(coefficients)[3] - RELATIVE_E_FLOOR

    constraints = [
        dict(type="eq", fun=equalities),
        dict(type="ineq", fun=positivity),
    ]
    starts = [np.zeros(len(RADIAL_INTERVALS))]
    starts.extend(
        rng.normal(0.0, 0.3, len(RADIAL_INTERVALS))
        for _ in range(START_COUNT - 1)
    )
    fits = []
    for start in starts:
        fits.append(
            minimize(
                objective,
                start,
                method="SLSQP",
                constraints=constraints,
                options=dict(
                    maxiter=MAX_ITERATIONS,
                    ftol=1e-11,
                    disp=False,
                ),
            )
        )

    candidates = []
    for index, fit in enumerate(fits):
        corrected_E, defects, slack, ratio = evaluate(fit.x)
        if (
            bool(fit.success)
            and float(np.max(np.abs(defects))) <= 1e-8
            and float(np.min(ratio)) >= RELATIVE_E_FLOOR - 1e-8
        ):
            candidates.append((slack["S_slack"], index, fit,
                               corrected_E, defects, slack, ratio))

    if candidates:
        _, selected_index, selected, corrected_E, defects, slack, ratio = max(
            candidates, key=lambda item: item[0]
        )
        row = dict(
            eta=float(eta),
            feasible=True,
            coefficients=selected.x.tolist(),
            selected_start=int(selected_index),
            starts_tried=len(starts),
            optimizer_success=bool(selected.success),
            optimizer_message=str(selected.message),
            I_Cp_defect=defects.tolist(),
            I_defect=float(defects[0]),
            Cp_defect=float(defects[1]),
            positivity=dict(
                floor=RELATIVE_E_FLOOR,
                min_relative_E=float(np.min(ratio)),
                min_E=float(np.min(corrected_E)),
                sampled_points=int(len(ratio)),
            ),
            optimistic_continuum_S_slack=slack,
        )
        return row

    # Keep a diagnostic of the closest attempted fit when the hard floor or
    # equalities are infeasible, but do not present it as a selected profile.
    best = min(
        fits,
        key=lambda fit: float(np.linalg.norm(evaluate(fit.x)[1]))
        + max(0.0, RELATIVE_E_FLOOR - float(np.min(evaluate(fit.x)[3])))
    )
    corrected_E, defects, slack, ratio = evaluate(best.x)
    return dict(
        eta=float(eta),
        feasible=False,
        coefficients=None,
        selected_start=None,
        starts_tried=len(starts),
        optimizer_success=bool(best.success),
        optimizer_message=str(best.message),
        I_Cp_defect=defects.tolist(),
        I_defect=float(defects[0]),
        Cp_defect=float(defects[1]),
        positivity=dict(
            floor=RELATIVE_E_FLOOR,
            min_relative_E=float(np.min(ratio)),
            min_E=float(np.min(corrected_E)),
            sampled_points=int(len(ratio)),
        ),
        optimistic_continuum_S_slack=slack,
        infeasibility_reason=(
            "No SLSQP start satisfied both exact I/Cp defects <= 1e-8 "
            "and the sampled E/E_target floor >= 0.11."
        ),
    )


def run():
    if any(lo <= INNER_CONE_CUTOFF for lo, _ in RADIAL_INTERVALS):
        raise ValueError("an E support reaches the unchanged inner cone")
    if any(RADIAL_INTERVALS[i][1] >= RADIAL_INTERVALS[i + 1][0]
           for i in range(len(RADIAL_INTERVALS) - 1)):
        raise ValueError("E radial supports are not ordered and disjoint")

    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    X, weights = quadrature()
    rng = np.random.default_rng(73073)
    rows = []
    for eta in TRAIN_ETAS:
        row = solve_eta(mean, target, base, X, weights, eta, rng)
        rows.append(row)
        print(
            json.dumps(
                dict(
                    eta=row["eta"],
                    feasible=row["feasible"],
                    I_Cp_defect=row["I_Cp_defect"],
                    min_relative_E=row["positivity"]["min_relative_E"],
                    S_slack=row["optimistic_continuum_S_slack"]["S_slack"],
                )
            ),
            flush=True,
        )

    report = dict(
        source_mean="delayed_remote_moment_repair.current_mean",
        source_target="high_frequency_shear_screen.make_field(16, 2.)",
        tau=TAU,
        moment_order=["M", "I", "J", "S", "Cp"],
        source_radial_intervals=[list(item) for item in SOURCE_INTERVALS],
        radial_intervals=[list(item) for item in RADIAL_INTERVALS],
        mode_family="SimilaritySwirlMode",
        mode_eta_interval=list(MODE_ETA_INTERVAL),
        inner_cone_cutoff=INNER_CONE_CUTOFF,
        u_correction_support=list(U_SUPPORT),
        training_etas=list(TRAIN_ETAS),
        relative_E_floor=RELATIVE_E_FLOOR,
        quadrature_order=QUADRATURE_ORDER,
        quadrature_points=int(len(X)),
        starts_per_eta=START_COUNT,
        rows=rows,
        scope=(
            "Fixed-slice E-stage capacity only. Each row uses independent "
            "coefficients on five disjoint radial SimilaritySwirlMode "
            "supports, all above X=1.03, and therefore leaves the sampled "
            "inner-cone velocity unchanged. I and Cp are solved on the "
            "displayed Gauss grid with a hard sampled E/E_target >= 0.11. "
            "The reported S slack is optimistic: it uses the quadrature "
            "minimum over an arbitrary U correction on X=1.04..2.98 "
            "preserving target M and J, without endpoint smoothness, "
            "eta/time compatibility, pressure, momentum, or PDE checks. "
            "Independent slice coefficients are not a global physical "
            "field or a solenoidal eta/time lift."
        ),
        accepted=False,
    )
    output = ROOT / "delayed_remote_positive_e.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            dict(
                output=str(output),
                feasible_rows=sum(bool(row["feasible"]) for row in rows),
                max_abs_I_Cp_defect=max(
                    max(abs(value) for value in row["I_Cp_defect"])
                    for row in rows
                ),
                min_relative_E=min(
                    row["positivity"]["min_relative_E"] for row in rows
                ),
                S_slacks=[
                    row["optimistic_continuum_S_slack"]["S_slack"]
                    for row in rows
                ],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    run()
