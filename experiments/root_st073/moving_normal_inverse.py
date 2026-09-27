"""Moving-normal principal amplitude inverse.

For a full wavevector ``n(t)`` (including any harmonic multiplier), this
module solves the physical-time principal equation

    a' = -P (K a + f) - n (n_dot . a) / |n|^2 - nu |n|^2 a,

where ``P = I - n n^T / |n|^2``.  The associated complex pressure amplitude is

    p = i (n . (K a + f) - n_dot . a) / |n|^2.

This is a local numerical inverse only.  It does not construct a spatial
field, a curl correction, a force, or a scale recursion.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline


PAPER_EQUATION = "7.13"
PAPER_URL = (
    "https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/"
    "navier-stokes.pdf"
)


def _finite(value, name):
    array = np.asarray(value)
    if np.iscomplexobj(array):
        valid = np.all(np.isfinite(array.real)) and np.all(np.isfinite(array.imag))
    else:
        valid = np.all(np.isfinite(array))
    if not valid:
        raise ValueError(f"{name} must contain only finite values")
    return array


def _real_vector(value, name):
    array = _finite(value, name)
    if array.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {array.shape}")
    if np.iscomplexobj(array) and np.max(np.abs(array.imag)) > 1e-14:
        raise ValueError(f"{name} must be real-valued")
    return np.asarray(array.real, dtype=float)


def _shape_value(value, shape, name, complex_allowed=True):
    array = _finite(value, name)
    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}")
    if not complex_allowed and np.iscomplexobj(array) and np.max(np.abs(array.imag)) > 1e-14:
        raise ValueError(f"{name} must be real-valued")
    dtype = complex if np.iscomplexobj(array) else float
    return np.asarray(array, dtype=dtype)


def _accessor(value, interval, shape, name, real=False):
    """Normalize a scalar callable or sampled array to a scalar accessor."""

    t0, t1 = interval
    if callable(value):
        def evaluate(t):
            result = value(float(t))
            if real:
                return _real_vector(result, name)
            return _shape_value(result, shape, name)

        return evaluate, "callable"

    array = _finite(value, name)
    if real:
        if array.shape == (3,) and (
            not np.iscomplexobj(array) or np.max(np.abs(array.imag)) <= 1e-14
        ):
            constant = np.asarray(array.real, dtype=float)
            return lambda _t: constant, "constant"
    elif array.shape == shape:
        constant = np.asarray(array, dtype=complex if np.iscomplexobj(array) else float)
        return lambda _t: constant, "constant"

    # A two-row or multi-row sampled prescription is interpolated smoothly in
    # physical time.  The spline is only queried after the public interval
    # check, and therefore cannot silently extrapolate.
    if array.ndim >= 1 and array.shape[1:] == shape:
        sample_count = array.shape[0]
        if sample_count < 2:
            raise ValueError(f"{name} sampled data needs at least two rows")
        sampled = np.asarray(array.real if real else array,
                             dtype=float if real else (
                                 complex if np.iscomplexobj(array) else float))
        spline = CubicSpline(
            np.array([t0, t1]) if sample_count == 2 else
            np.linspace(t0, t1, sample_count),
            sampled,
            axis=0,
            extrapolate=False,
        )

        def evaluate(t):
            result = spline(float(t))
            if real:
                return _real_vector(result, name)
            return _shape_value(result, shape, name)

        return evaluate, "sampled"

    raise ValueError(f"{name} must be a callable or have shape {shape}")


class MovingNormalSolution(Mapping):
    """Dense, non-extrapolating amplitude and pressure solution."""

    def __init__(self, *, interval, time_scale, dense, normal, normal_dot,
                 matrix, source, viscosity, metadata):
        self.interval = tuple(float(x) for x in interval)
        self.time_scale = float(time_scale)
        self._dense = dense
        self._normal = normal
        self._normal_dot = normal_dot
        self._matrix = matrix
        self._source = source
        self.nu = float(viscosity)
        self.metadata = metadata

    def _query(self, times):
        query = np.asarray(times, dtype=float)
        if not np.all(np.isfinite(query)):
            raise ValueError("query times must be finite")
        if np.any(query < self.interval[0]) or np.any(query > self.interval[1]):
            raise ValueError(
                f"query time outside closed interval {self.interval}"
            )
        scalar = query.ndim == 0
        shape = query.shape
        flat = query.reshape(-1)
        scaled = (flat - self.interval[0]) / self.time_scale
        return flat, scaled, shape, scalar

    def amplitude(self, times):
        """Evaluate the complex 3-vector amplitude at scalar or vector times."""

        flat, scaled, shape, scalar = self._query(times)
        values = np.asarray(self._dense(scaled), dtype=complex).T
        if scalar:
            return values[0]
        return values.reshape(shape + (3,))

    def pressure(self, times):
        """Evaluate the complex pressure amplitude at scalar or vector times."""

        flat, _scaled, shape, scalar = self._query(times)
        amplitudes = self.amplitude(flat).reshape((-1, 3))
        normals = np.asarray([self._normal(t) for t in flat], dtype=float)
        normal_dots = np.asarray([self._normal_dot(t) for t in flat], dtype=float)
        matrices = np.asarray([self._matrix(t) for t in flat])
        sources = np.asarray([self._source(t) for t in flat], dtype=complex)
        n_squared = np.einsum("ni,ni->n", normals, normals)
        if np.any(n_squared <= 0.0):
            raise ValueError("normal must be nonzero over the query interval")
        ka_f = sources + np.einsum("nij,nj->ni", matrices, amplitudes)
        numerator = np.einsum("ni,ni->n", normals, ka_f)
        numerator -= np.einsum("ni,ni->n", normal_dots, amplitudes)
        values = 1j * numerator / n_squared
        if scalar:
            return values[0]
        return values.reshape(shape)

    def __getitem__(self, key):
        if key == "amplitude":
            return self.amplitude
        if key == "pressure":
            return self.pressure
        if key == "metadata":
            return self.metadata
        raise KeyError(key)

    def __iter__(self):
        return iter(("amplitude", "pressure", "metadata"))

    def __len__(self):
        return 3


def solve_path(time_interval, normal, normal_dot, K, source, nu,
               initial_amplitude=None, *, rtol=1e-10, atol=1e-12,
               method="DOP853", max_step=np.inf):
    """Solve the moving-normal principal amplitude equation.

    ``time_interval`` is a strictly increasing pair of physical times.  Each
    prescription may be a callable ``fn(t)`` or sampled values with rows over
    that interval.  The returned :class:`MovingNormalSolution` evaluates only
    inside the closed interval and never extrapolates.
    """

    interval = np.asarray(time_interval, dtype=float)
    if interval.shape != (2,) or not np.all(np.isfinite(interval)):
        raise ValueError("time_interval must be two finite physical times")
    t0, t1 = map(float, interval)
    time_scale = t1 - t0
    if time_scale <= 0.0:
        raise ValueError("time_interval must be strictly increasing")
    viscosity = float(nu)
    if not np.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("nu must be finite and nonnegative")
    if rtol <= 0.0 or atol <= 0.0:
        raise ValueError("rtol and atol must be positive")

    normal_fn, normal_kind = _accessor(normal, (t0, t1), (3,), "normal", real=True)
    normal_dot_fn, normal_dot_kind = _accessor(
        normal_dot, (t0, t1), (3,), "normal_dot", real=True
    )
    matrix_fn, matrix_kind = _accessor(K, (t0, t1), (3, 3), "K")
    source_fn, source_kind = _accessor(source, (t0, t1), (3,), "source")

    def geometry(t):
        nv = normal_fn(t)
        ndv = normal_dot_fn(t)
        n_squared = float(np.dot(nv, nv))
        if not np.isfinite(n_squared) or n_squared <= 0.0:
            raise ValueError("normal must be nonzero throughout the interval")
        return nv, ndv, n_squared

    initial = np.zeros(3, dtype=complex) if initial_amplitude is None else np.asarray(
        initial_amplitude, dtype=complex
    )
    if initial.shape != (3,) or not np.all(np.isfinite(initial.real)) or not np.all(np.isfinite(initial.imag)):
        raise ValueError("initial_amplitude must be a finite complex vector of shape (3,)")
    initial_normal, _initial_dot, initial_n_squared = geometry(t0)
    initial_longitudinal = abs(np.dot(initial_normal, initial))
    initial_scale = max(float(np.linalg.norm(initial_normal) * np.linalg.norm(initial)), 1.0)
    if initial_longitudinal > 1e-9 * initial_scale:
        raise ValueError("initial_amplitude must be transverse to normal(t0)")

    def physical_rhs(t, amplitude):
        nv, ndv, n_squared = geometry(t)
        matrix = _shape_value(matrix_fn(t), (3, 3), "K")
        forcing = _shape_value(source_fn(t), (3,), "source")
        combined = matrix @ amplitude + forcing
        projection = np.eye(3) - np.outer(nv, nv) / n_squared
        return (
            -projection @ combined
            - nv * (np.dot(ndv, amplitude) / n_squared)
            - viscosity * n_squared * amplitude
        )

    def normalized_rhs(s, amplitude):
        # solve_ivp supplies normalized endpoints exactly, but floating-point
        # reconstruction of t0 + time_scale * s can land one ulp outside a
        # sampled CubicSpline's closed domain.  Clamp only this internal RHS
        # evaluation; public dense queries still reject extrapolation.
        physical_time = min(max(t0 + time_scale * float(s), t0), t1)
        return time_scale * physical_rhs(physical_time, amplitude)

    solved = solve_ivp(
        normalized_rhs,
        (0.0, 1.0),
        initial,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        dense_output=True,
    )
    if not solved.success:
        raise RuntimeError(solved.message)

    solver_times = (t0 + time_scale * np.asarray(solved.t)).tolist()
    metadata = dict(
        success=True,
        equation=PAPER_EQUATION,
        paper_url=PAPER_URL,
        time_interval=[t0, t1],
        normalized_interval=[0.0, 1.0],
        time_scale=time_scale,
        viscosity=viscosity,
        solver_method=method,
        solver_steps=int(len(solved.t)),
        solver_times=[float(t) for t in solver_times],
        normal_input=normal_kind,
        normal_dot_input=normal_dot_kind,
        matrix_input=matrix_kind,
        source_input=source_kind,
        initial_amplitude=[[float(x.real), float(x.imag)] for x in initial],
        initial_transversality=float(initial_longitudinal),
        extrapolation=False,
    )
    return MovingNormalSolution(
        interval=(t0, t1),
        time_scale=time_scale,
        dense=solved.sol,
        normal=normal_fn,
        normal_dot=normal_dot_fn,
        matrix=matrix_fn,
        source=source_fn,
        viscosity=viscosity,
        metadata=metadata,
    )


solve_moving_normal = solve_path


def _rotating_manufactured_case():
    t0, t1 = 0.0, 0.08
    viscosity = 0.07
    harmonic = 3.0
    omega = 1.7
    lam = 0.4 + 0.2j
    c = 0.7 - 0.2j
    d = -0.35 + 0.25j

    def frame(t):
        angle = omega * float(t)
        return np.array([np.cos(angle), np.sin(angle), 0.0])

    def transverse_frame(t):
        angle = omega * float(t)
        return np.array([-np.sin(angle), np.cos(angle), 0.0])

    def n(t):
        return harmonic * frame(t)

    def n_dot(t):
        return harmonic * omega * transverse_frame(t)

    def amplitude_exact(t):
        return np.exp(lam * float(t)) * (
            c * transverse_frame(t) + d * np.array([0.0, 0.0, 1.0])
        )

    def amplitude_dot_exact(t):
        value = float(t)
        return np.exp(lam * value) * (
            lam * (c * transverse_frame(value) + d * np.array([0.0, 0.0, 1.0]))
            - c * omega * frame(value)
        )

    def matrix(t):
        value = float(t)
        return np.array([
            [0.2 + 0.3 * value, -0.4, 0.1],
            [0.5, -0.1 + 0.2 * value, 0.25],
            [-0.15, 0.35, 0.3 - 0.1 * value],
        ])

    def gamma(t):
        return (0.35 + 0.15j) + (0.2 - 0.1j) * float(t)

    def source_fn(t):
        value = float(t)
        nv = n(value)
        ndv = n_dot(value)
        av = amplitude_exact(value)
        av_dot = amplitude_dot_exact(value)
        n_squared = float(np.dot(nv, nv))
        transverse_rhs = (
            -av_dot
            - nv * (np.dot(ndv, av) / n_squared)
            - viscosity * n_squared * av
        )
        g = transverse_rhs + nv * gamma(value)
        return g - matrix(value) @ av

    solution = solve_path(
        (t0, t1), n, n_dot, matrix, source_fn, viscosity,
        initial_amplitude=amplitude_exact(t0),
    )
    return solution, dict(
        interval=(t0, t1), n=n, n_dot=n_dot, matrix=matrix,
        source=source_fn, amplitude_exact=amplitude_exact,
        pressure_exact=lambda t: 1j * (
            gamma(t)
            - np.dot(n_dot(t), amplitude_exact(t))
            / np.dot(n(t), n(t))
        ),
    )


def _fixed_manufactured_case():
    t0, t1 = -0.03, 0.03
    viscosity = 0.11
    harmonic = 2.5
    lam = -0.3 + 0.35j
    c = -0.45 + 0.15j
    d = 0.25 + 0.4j
    nv = np.array([0.0, 0.0, harmonic])
    zero = np.zeros(3)

    def n(_t):
        return nv

    def n_dot(_t):
        return zero

    def amplitude_exact(t):
        value = np.exp(lam * float(t))
        return value * np.array([c, d, 0.0])

    def matrix(t):
        value = float(t)
        return np.array([
            [0.1, 0.3 + value, -0.2],
            [-0.5, 0.4, 0.15 - value],
            [0.2, -0.1, 0.25],
        ])

    def gamma(t):
        return -0.25 + 0.2j + (0.15 + 0.05j) * float(t)

    def source_fn(t):
        value = float(t)
        av = amplitude_exact(value)
        av_dot = lam * av
        n_squared = float(np.dot(nv, nv))
        transverse_rhs = -av_dot - viscosity * n_squared * av
        return transverse_rhs + nv * gamma(value) - matrix(value) @ av

    solution = solve_path(
        (t0, t1), n, n_dot, matrix, source_fn, viscosity,
        initial_amplitude=amplitude_exact(t0),
    )
    return solution, dict(
        interval=(t0, t1), n=n, n_dot=n_dot, matrix=matrix,
        source=source_fn, amplitude_exact=amplitude_exact,
        pressure_exact=lambda t: 1j * gamma(t),
    )


def _finite_difference_defect(solution, prescriptions, count=13):
    t0, t1 = prescriptions["interval"]
    scale = t1 - t0
    half_step = scale * 5e-4
    candidates = np.linspace(t0 + 0.02 * scale, t1 - 0.02 * scale, count)
    nodes = np.asarray(solution.metadata["solver_times"], dtype=float)
    candidates = np.asarray([
        t for t in candidates
        if np.min(np.abs(nodes - t)) > 2.0 * half_step
    ])
    rows = []
    for t in candidates:
        derivative = (
            solution.amplitude(t + half_step)
            - solution.amplitude(t - half_step)
        ) / (2.0 * half_step)
        av = solution.amplitude(t)
        nv = prescriptions["n"](t)
        ndv = prescriptions["n_dot"](t)
        kv = prescriptions["matrix"](t)
        fv = prescriptions["source"](t)
        n_squared = float(np.dot(nv, nv))
        projection = np.eye(3) - np.outer(nv, nv) / n_squared
        defect = (
            derivative
            + projection @ (kv @ av + fv)
            + nv * (np.dot(ndv, av) / n_squared)
            + solution.nu * n_squared * av
        )
        rows.append((t, derivative, defect))
    if not rows:
        raise RuntimeError("no finite-difference query times remained outside solve nodes")
    return np.asarray([row[2] for row in rows]), candidates, half_step


def run_checks():
    rotating, rotating_data = _rotating_manufactured_case()
    fixed, fixed_data = _fixed_manufactured_case()

    def errors(solution, data):
        interval = data["interval"]
        query = np.linspace(interval[0], interval[1], 23)
        amplitude_error = np.max(np.linalg.norm(
            solution.amplitude(query)
            - np.asarray([data["amplitude_exact"](t) for t in query]), axis=1
        ))
        pressure_error = np.max(np.abs(
            solution.pressure(query)
            - np.asarray([data["pressure_exact"](t) for t in query])
        ))
        source_longitudinal = np.max(np.abs([
            np.dot(data["n"](t), data["source"](t))
            / np.dot(data["n"](t), data["n"](t))
            for t in query
        ]))
        transverse_error = np.max(np.abs(
            np.einsum(
                "ni,ni->n",
                np.asarray([data["n"](t) for t in query]),
                solution.amplitude(query),
            )
        ))
        defects, check_times, half_step = _finite_difference_defect(solution, data)
        return dict(
            amplitude_error_max=float(amplitude_error),
            pressure_error_max=float(pressure_error),
            source_longitudinal_component_max=float(source_longitudinal),
            transverse_constraint_max=float(transverse_error),
            finite_difference_defect_max=float(np.max(np.linalg.norm(defects, axis=1))),
            independent_check_count=int(len(check_times)),
            finite_difference_half_step=float(half_step),
            excluded_from_solver_nodes=True,
        )

    rotating_errors = errors(rotating, rotating_data)
    fixed_errors = errors(fixed, fixed_data)

    outside_errors = {}
    for label, solution, interval in (
        ("rotating", rotating, rotating_data["interval"]),
        ("fixed", fixed, fixed_data["interval"]),
    ):
        try:
            solution.amplitude(interval[0] - 1e-12)
        except ValueError:
            outside_errors[label] = True
        else:
            outside_errors[label] = False

    report = dict(
        equation=PAPER_EQUATION,
        paper_url=PAPER_URL,
        method="normalized physical-time solve with dense non-extrapolating output",
        rotating_normal=dict(
            metadata=rotating.metadata,
            checks=rotating_errors,
            source_longitudinal_nonzero=True,
        ),
        fixed_normal=dict(metadata=fixed.metadata, checks=fixed_errors),
        extrapolation_rejected=outside_errors,
        scope=(
            "Manufactured complex amplitudes on two tiny physical intervals; "
            "moving-normal and fixed-normal pressure identities checked with "
            "finite differences at query times excluded from solve nodes. "
            "No spatial field, curl reconstruction, full residual, or scale recursion."
        ),
        accepted=False,
        pde_validated=False,
        scale_recursion_established=False,
    )
    return report


if __name__ == "__main__":
    result = run_checks()
    output = Path(__file__).with_name("moving_normal_inverse_check.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
