"""Read-only prototype of a globally localized current ST073 candidate.

The existing local candidate is wrapped before its compact Fourier wave is
added.  The mean meridional velocity is represented by an axisymmetric
Stokes streamfunction, so an axial cutoff changes ``psi`` and preserves
divergence analytically.  Pure swirl directions may be multiplied directly
by the same scalar cutoff.  The pressure is multiplied by the scalar cutoff
explicitly; the resulting cutoff pressure-gradient terms are diagnostic and
are not a Navier--Stokes closure.

This module is deliberately a bounded construction check.  It does not
modify the local candidate or claim a PDE, finite-time, or scale-recursion
acceptance result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import Legendre
from scipy.special import expit


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_constrained import load_saved_field  # noqa: E402
from broad_shear_dynamic_control import load_saved_field as load_dynamic  # noqa: E402
from enriched_shape_replay import Mode0TangentCorrection  # noqa: E402
from full_wave_tangent import LocalPotentialField, _unpack_full  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from joined_field import coordinates, independent_fd  # noqa: E402
from poloidal_bridge import coefficients  # noqa: E402
from separated_moment_modes import flat_bump, scale_partition  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402
from wave_higher_harmonic_tangent import (  # noqa: E402
    Mode2TangentCorrection,
    _decode_mode_block,
)


CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
SNAPSHOT_PATH = ROOT / "full_wave_frozen_cache.json"
OUTPUT_PATH = ROOT / "global_axial_extension.json"

# The frozen wave support reaches |eta| ~= 0.278 at the reference state.
# Keep the full support on the exact streamfunction plateau while retaining a
# nonempty smooth collar before the registered outer edge.
ETA_FLAT = 0.30
ETA_OUTER = 0.49
RADIAL_CUTOFF_RATIO = 2.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.isfinite(points).all():
        raise ValueError("points must be finite")
    return points


def _joined_streamfunction(joined, points, tau):
    """Return the existing JoinedField meridional streamfunction.

    This is the same polynomial construction used by
    ``axial_compact_join.AxiallyCompactField.streamfunction``.  It uses the
    inner coefficient primitive and ``poloidal_bridge.coefficients``; no
    radial quadrature is introduced.
    """

    points = _points(points)
    tau = float(tau)
    inner = joined.inner
    sn = math.sqrt(joined.nu)
    radius = np.hypot(points[:, 0], points[:, 1])
    coord = coordinates(radius / sn, points[:, 2] / sn, tau, inner.h)
    result = np.zeros(len(points), dtype=float)
    for index, (r, X, eta, q) in enumerate(
        zip(radius, coord["X"], coord["eta"], coord["q"])
    ):
        ri, ro, _, bridge = coefficients(
            inner, float(eta), tau, joined.join_X, joined.outer_ratio
        )
        if r <= ri:
            c = inner.coefficients(float(eta), float(q))[2, :, 0]
            primitive = np.r_[0.0, c / np.arange(1, len(c) + 1)]
            result[index] = joined.nu**1.5 * q * np.polynomial.polynomial.polyval(
                X, primitive
            )
        elif r < ro:
            y = (r - ri) / (ro - ri)
            result[index] = np.polynomial.polynomial.polyval(y, bridge)
    return result


def _wide_poloidal_increment(mode, points, tau):
    """Analytic psi increment for ``wide_modes.WidePoloidalModes``."""

    points = _points(points)
    tau = float(tau)
    radius = np.hypot(points[:, 0], points[:, 1])
    sn = math.sqrt(mode.nu)
    coord = coordinates(radius / sn, points[:, 2] / sn, tau, mode.inner.h)
    q = np.asarray(coord["q"], dtype=float)
    eta = np.asarray(coord["eta"], dtype=float)
    ri = np.sqrt(2.0 * mode.nu * q * mode.join_X)
    width = (mode.ratio - 1.0) * ri
    raw_y = (radius - ri) / width
    active = (raw_y > 0.0) & (raw_y < 1.0)
    y = np.clip(raw_y, 0.0, 1.0)
    psi_scale = mode.nu**1.5 * (2.0 * mode.join_X) * q ** (1.0 - mode.inner.A)
    result = np.zeros(len(points), dtype=float)
    for j in range(3):
        legendre = Legendre.basis(j).convert(kind=Poly)(Poly([-1.0, 2.0]))
        bubble = 256.0 * Poly([0.0, 1.0]) ** 4 * Poly([1.0, -1.0]) ** 4
        bubble = bubble * legendre
        radial_shape = np.polynomial.polynomial.polyval(y, bubble.coef)
        radial_shape = np.where(active, radial_shape, 0.0)
        for axial_index in range(2):
            axial = (eta / 0.3) ** axial_index
            result += (
                mode.a[0, j, axial_index]
                * psi_scale
                * radial_shape
                * axial
            )
    return result


def _separated_poloidal_increment(mode, points, tau):
    """Analytic psi increment for ``SeparatedMomentModes``."""

    points = _points(points)
    tau = float(tau)
    radius = np.hypot(points[:, 0], points[:, 1])
    sn = math.sqrt(mode.nu)
    coord = coordinates(radius / sn, points[:, 2] / sn, tau, mode.inner.h)
    q = np.asarray(coord["q"], dtype=float)
    eta = np.asarray(coord["eta"], dtype=float)
    ri = np.sqrt(2.0 * mode.nu * q * mode.join_X)
    width = (mode.ratio - 1.0) * ri
    y = (radius - ri) / width
    k = -math.log2(2.0 * tau)
    time_weights = scale_partition(k, mode.knots)
    psi_scale = mode.nu**1.5 * (2.0 * mode.join_X) * q ** (1.0 - mode.inner.A)
    result = np.zeros(len(points), dtype=float)
    for radial_index, (lo, hi) in enumerate(mode.windows):
        bump, _ = flat_bump(y, lo, hi)
        for axial_index, power in enumerate(mode.axial_powers):
            axial = (eta / 0.3) ** power
            coefficient = float(
                np.dot(time_weights, mode.a[:, 1, radial_index, axial_index])
            )
            result += coefficient * psi_scale * bump * axial
    return result


def _streamfunction(node, points, tau, trace=None):
    """Reconstruct psi through the current mean wrapper hierarchy.

    Known wrappers either preserve the lower meridional streamfunction or
    add one of the two explicit polynomial streamfunctions above.  Pure swirl
    and pressure wrappers contribute zero psi.  Unknown wrappers raise rather
    than silently dropping a meridional component.
    """

    points = _points(points)
    name = type(node).__name__
    if trace is not None:
        trace.add(name)

    if name in {"ZeroBackground", "CompactPressureDirection"}:
        return np.zeros(len(points), dtype=float)
    if name == "JoinedField":
        return _joined_streamfunction(node, points, tau)
    if name == "GroupedJoinedField":
        return _joined_streamfunction(node.joined, points, tau)
    if name == "CachedAdaptiveBridge":
        return _streamfunction(node.joined, points, tau, trace)
    if name == "WidePoloidalModes":
        return _streamfunction(node.base, points, tau, trace) + _wide_poloidal_increment(
            node, points, tau
        )
    if name == "WideJointModes":
        # WideJointModes.field is WideAngularModes(WidePoloidalModes(...));
        # the angular wrapper adds swirl only.
        return _streamfunction(node.base, points, tau, trace) + _wide_poloidal_increment(
            node.field.base, points, tau
        )
    if name == "WideAngularModes":
        return _streamfunction(node.base, points, tau, trace)
    if name == "SeparatedMomentModes":
        return _streamfunction(node.base, points, tau, trace) + _separated_poloidal_increment(
            node, points, tau
        )
    if name == "MeridionalSlope":
        k = -math.log2(2.0 * float(tau))
        return (k - node.k0) * _streamfunction(node.unit, points, tau, trace)
    if name == "CorrectedField":
        result = _streamfunction(node.baseline, points, tau, trace)
        for coefficient, mode in zip(node.coefficients, node.modes):
            if coefficient:
                result = result + float(coefficient) * _streamfunction(
                    mode, points, tau, trace
                )
        return result
    if name == "StatefulMean":
        k = -math.log2(2.0 * float(tau))
        amplitudes = node.state + (k - node.kref) * node.slope
        result = _streamfunction(node.base, points, tau, trace)
        for coefficient, unit in zip(amplitudes, node.values):
            if coefficient:
                result = result + float(coefficient) * _streamfunction(
                    unit, points, tau, trace
                )
        return result
    if name == "BroadCentrifugalPressure":
        return _streamfunction(node.modified, points, tau, trace)
    if name == "DynamicControlledMean":
        # The dynamic wrapper changes only broad/swirl amplitudes and scalar
        # pressure around its lower fixed mean.
        return _streamfunction(node.base, points, tau, trace)
    if name == "BroadShearSlope":
        # BroadShearSlope.fields multiplies its lower field by (k-k0).  Keep
        # the same factor in the streamfunction traversal so a nonzero
        # meridional lower base remains divergence-preserving as well.
        k = -math.log2(2.0 * float(tau))
        return (k - node.k0) * _streamfunction(node.base, points, tau, trace)
    if name in {
        "OuterPressure",
        "PressureBubbleField",
        "SwirlValue",
        "OuterSwirlSlope",
        "BroadAnnularShear",
    }:
        return _streamfunction(node.base, points, tau, trace)
    raise TypeError(f"No analytic streamfunction adapter for {name}")


class GlobalAxialExtension:
    """Divergence-preserving axial and optional outer-radial localization."""

    def __init__(
        self,
        base,
        *,
        eta_flat=ETA_FLAT,
        eta_outer=ETA_OUTER,
        radial_cutoff_ratio=RADIAL_CUTOFF_RATIO,
    ):
        self.base = base
        self.joined = _find_joined(base)
        self.inner = base.inner
        self.nu = float(base.nu)
        self.join_X = float(base.join_X)
        self.ratio = float(base.ratio)
        self.eta_flat = float(eta_flat)
        self.eta_outer = float(eta_outer)
        self.radial_cutoff_ratio = float(radial_cutoff_ratio)
        if not 0.0 < self.eta_flat < self.eta_outer < self.inner.p.eta_max:
            raise ValueError("Axial cutoff must lie inside the registered eta slab")
        if self.radial_cutoff_ratio <= 1.0:
            raise ValueError("Radial cutoff ratio must exceed one")
        self.tau_min = 0.5 * 2.0 ** (-float(self.inner.p.k_max))
        self.tau_max = 0.5

    def cutoff(self, eta):
        """Return ``chi(eta)`` and its derivative, matching AxiallyCompactField."""

        eta = np.asarray(eta, dtype=float)
        absolute = np.abs(eta)
        chi = np.zeros_like(eta)
        derivative = np.zeros_like(eta)
        chi[absolute <= self.eta_flat] = 1.0
        middle = (absolute > self.eta_flat) & (absolute < self.eta_outer)
        s = (absolute[middle] - self.eta_flat) / (
            self.eta_outer - self.eta_flat
        )
        value = expit(1.0 / s - 1.0 / (1.0 - s))
        chi[middle] = value
        derivative[middle] = (
            value
            * (1.0 - value)
            * (-1.0 / s**2 - 1.0 / (1.0 - s) ** 2)
            * np.sign(eta[middle])
            / (self.eta_outer - self.eta_flat)
        )
        return chi, derivative

    def radial_radii(self, tau):
        qmax = float(tau) / (1.0 - self.eta_outer**2)
        r1 = self.ratio * math.sqrt(2.0 * self.nu * qmax * self.join_X)
        return r1, self.radial_cutoff_ratio * r1

    def radial_cutoff(self, radius, tau):
        radius = np.asarray(radius, dtype=float)
        r1, r2 = self.radial_radii(tau)
        rho = np.ones_like(radius)
        derivative = np.zeros_like(radius)
        middle = (radius > r1) & (radius < r2)
        s = (radius[middle] - r1) / (r2 - r1)
        value = expit(1.0 / (1.0 - s) - 1.0 / s)
        rho[middle] = 1.0 - value
        value_derivative = value * (1.0 - value) * (
            1.0 / (1.0 - s) ** 2 + 1.0 / s**2
        )
        derivative[middle] = -value_derivative / (r2 - r1)
        rho[radius >= r2] = 0.0
        return rho, derivative

    def fields(self, points, tau):
        points = _points(points)
        times = np.broadcast_to(np.asarray(tau, dtype=float), (len(points),))
        if len(points) == 0:
            return np.empty((0, 3), dtype=float), np.empty(0, dtype=float)
        if not np.all(times == times[0]):
            raise ValueError("GlobalAxialExtension requires one common time")
        tau = float(times[0])
        if not self.tau_min <= tau <= self.tau_max:
            raise ValueError("Outside the registered finite time slab")

        sn = math.sqrt(self.nu)
        radius = np.hypot(points[:, 0], points[:, 1])
        coord = coordinates(radius / sn, points[:, 2] / sn, tau, self.inner.h)
        chi, chi_eta = self.cutoff(coord["eta"])
        rho, rho_r = self.radial_cutoff(radius, tau)
        active = (chi > 0.0) & (rho > 0.0)
        velocity = np.zeros_like(points)
        pressure = np.zeros(len(points), dtype=float)
        if not np.any(active):
            return velocity, pressure

        chosen = points[active]
        base_velocity, base_pressure = self.base.fields(chosen, tau)
        psi = _streamfunction(self.base, chosen, tau)
        chosen_radius = radius[active]
        safe = np.maximum(chosen_radius, 1.0e-300)
        radial_unit = np.column_stack((chosen[:, 0], chosen[:, 1])) / safe[:, None]
        azimuthal_unit = np.column_stack((-chosen[:, 1], chosen[:, 0])) / safe[:, None]
        radial_velocity = np.sum(base_velocity[:, :2] * radial_unit, axis=1)
        swirl_velocity = np.sum(base_velocity[:, :2] * azimuthal_unit, axis=1)
        axial_velocity = base_velocity[:, 2]
        chi_z = chi_eta[active] * np.asarray(coord["eta_z"])[active] / sn

        # This is the cylindrical curl of psi_tilde = chi(eta) rho(r) psi:
        # ur = -psi_tilde_z/r, uz = psi_tilde_r/r.  Swirl is divergence-free
        # by axisymmetry and can use the same scalar cutoff directly.
        radial = (
            chi[active] * rho[active] * radial_velocity
            - rho[active] * chi_z * psi / safe
        )
        axial = (
            chi[active] * rho[active] * axial_velocity
            + chi[active] * rho_r[active] * psi / safe
        )
        swirl = chi[active] * rho[active] * swirl_velocity
        velocity[active, :2] = radial[:, None] * radial_unit + swirl[:, None] * azimuthal_unit
        velocity[active, 2] = axial
        pressure[active] = chi[active] * rho[active] * np.asarray(base_pressure)
        return velocity, pressure


def _find_joined(node):
    """Find the actual JoinedField below the current wrappers."""

    current = node
    visited = set()
    while id(current) not in visited:
        visited.add(id(current))
        name = type(current).__name__
        if name == "JoinedField":
            return current
        if name == "GroupedJoinedField":
            current = current.joined
        elif name == "CachedAdaptiveBridge":
            current = current.joined
        elif name == "CorrectedField":
            current = current.baseline
        elif name == "BroadCentrifugalPressure":
            current = current.modified
        elif hasattr(current, "base"):
            current = current.base
        else:
            break
    raise TypeError("Could not locate JoinedField below current mean")


def _decode_candidate(candidate_path=CANDIDATE_PATH):
    candidate = json.loads(Path(candidate_path).read_text(encoding="utf-8"))
    snapshot = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    if candidate.get("status") != "completed" or not candidate.get("assembled_feasible"):
        raise ValueError("The source candidate is not a completed feasible report")
    control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    if control.shape != (264,):
        raise ValueError(f"Expected 264 tangent controls, got {control.shape}")
    packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    if packed.shape != (27, 2):
        raise ValueError(f"Expected 27 complex wave coefficients, got {packed.shape}")
    return candidate, snapshot, control, packed[:, 0] + 1j * packed[:, 1]


def build_candidate(candidate_path=CANDIDATE_PATH):
    candidate, snapshot, control, wave = _decode_candidate(candidate_path)
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    geometry = snapshot["inputs"]["wave"]
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    localized_mean = GlobalAxialExtension(mean)

    # The degree-2 LocalPotentialField retains the locked mode-1 wave and its
    # mode-1 tangent block.  Degree-3 mode-0 and mode-2 corrections are the
    # existing exact compact basis wrappers; they are added after localization
    # so they are not componentwise tapered.
    legacy = np.r_[np.zeros(36), control[64:136], np.zeros(72)]
    derivatives, pressures = _unpack_full(legacy, 9)
    field = LocalPotentialField(
        localized_mean,
        geometry["center"],
        geometry["widths"],
        2,
        {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier},
        (
            np.zeros(27, dtype=complex),
            wave,
            np.zeros(27, dtype=complex),
        ),
        derivatives,
        pressures,
        tau0,
    )
    mode0 = Mode0TangentCorrection(field, geometry, control[:64], tau0)
    mode2_velocity, mode2_pressure = _decode_mode_block(control[-128:], 3)
    full = Mode2TangentCorrection(
        mode0,
        geometry["center"],
        geometry["widths"],
        carrier,
        3,
        mode2_velocity,
        mode2_pressure,
        tau0,
    )
    full.nu = mean.nu
    return full, localized_mean, mean, candidate, snapshot, mean_report


def _divergence_fd(field, points, tau, h):
    points = _points(points)
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        step = np.zeros(3, dtype=float)
        step[axis] = h
        plus = field.fields(points + step, tau)[0]
        minus = field.fields(points - step, tau)[0]
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)
    return divergence


def _json_float(value):
    return float(value) if np.isfinite(value) else None


def run(candidate_path=CANDIDATE_PATH, output_path=OUTPUT_PATH):
    started = time.perf_counter()
    candidate_path = Path(candidate_path)
    output_path = Path(output_path)
    raw_candidate = candidate_path.read_bytes()
    raw_snapshot = SNAPSHOT_PATH.read_bytes()
    full, localized, mean, candidate, snapshot, mean_report = build_candidate(candidate_path)
    tau = float(snapshot["inputs"]["mean"]["tau"])
    k = float(snapshot["inputs"]["mean"]["k"])
    geometry = snapshot["inputs"]["wave"]
    center = np.asarray(geometry["center"], dtype=float)
    widths = np.asarray(geometry["widths"], dtype=float)
    r1, r2 = localized.radial_radii(tau)
    qmax = tau / (1.0 - localized.eta_outer**2)
    zcap = math.sqrt(localized.nu) * qmax ** (0.5 - localized.inner.h) * localized.eta_outer
    report = {
        "status": "initialized",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "global_field_ready": False,
        "scope": (
            "Bounded instantaneous global-localization prototype for the frozen "
            "enriched mean/wave candidate. The mean is localized through an "
            "axisymmetric streamfunction and the pressure is explicitly tapered. "
            "The compact Fourier wave corrections remain exact curls and are "
            "added after mean localization. No PDE, forcing, finite-time, or "
            "scale-recursion acceptance is claimed."
        ),
        "sources": {
            "candidate": {"path": candidate_path.name, "sha256": _sha256(candidate_path)},
            "frozen_cache": {"path": SNAPSHOT_PATH.name, "sha256": _sha256(SNAPSHOT_PATH)},
            "mean_source": {
                "path": mean_report.get("source_dynamic_report", "broad_shear_dynamic_control.json"),
                "coefficients": mean_report.get("coefficients"),
            },
        },
        "inputs": {
            "k": k,
            "tau": tau,
            "physical_time": "t=-tau",
            "nu": localized.nu,
            "eta_flat": localized.eta_flat,
            "eta_outer": localized.eta_outer,
            "registered_tau_interval": [localized.tau_min, localized.tau_max],
            "join_X": localized.join_X,
            "join_ratio": localized.ratio,
            "radial_cutoff_ratio": localized.radial_cutoff_ratio,
            "wave_center": center.tolist(),
            "wave_widths": widths.tolist(),
            "wave_mode": int(geometry["mode"]),
            "wave_degree": int(geometry["degree"]),
            "tangent_control_count": int(len(candidate["selected"]["tangent_coefficients"])),
        },
        "streamfunction": {
            "construction": (
                "JoinedField coefficient primitive + poloidal_bridge.septic, "
                "WidePoloidalModes polynomial bubble, SeparatedMomentModes "
                "flat_bump psi; MeridionalSlope contributes (k-k0) times its "
                "explicit unit psi. Pure swirl and pressure wrappers contribute zero."
            ),
            "radial_formula": "u_r=-(chi*rho*psi)_z/r; u_z=(chi*rho*psi)_r/r",
            "wave_handling": "LocalPotentialField and Mode2TangentCorrection are added after mean localization; supported_fourier_basis supplies exact compact curls.",
            "numerical_radial_integral_added": False,
        },
        "support": {
            "eta_interval": [-localized.eta_outer, localized.eta_outer],
            "z_halfwidth_at_reference": zcap,
            "radial_interval": [0.0, r2],
            "radial_plateau_end": r1,
            "radial_tail_energy_bound_beyond_r2": 0.0,
            "axial_tail_energy_bound_beyond_zcap": 0.0,
            "finite_spatial_support_at_each_registered_tau": True,
            "support_reason": "chi and rho are flat C-infinity cutoffs; exact-curl wave basis is compact in (r,z).",
        },
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    # Reconstruct the streamfunction on a small representative point set.  A
    # failure here is recorded as a construction gap instead of silently
    # dropping a meridional wrapper.
    trace = set()
    reference_points = np.array(
        [
            [center[0], 0.0, center[1]],
            [center[0] + 0.5 * widths[0], 0.0, center[1]],
            [center[0], 0.0, center[1] + 0.5 * widths[1]],
        ],
        dtype=float,
    )
    try:
        psi_reference = _streamfunction(mean, reference_points, tau, trace)
        report["streamfunction"].update(
            wrapper_trace=sorted(trace),
            unknown_wrapper_count=0,
            reference_values=psi_reference.tolist(),
        )
    except Exception as error:  # noqa: BLE001 - preserve explicit construction gap
        report["status"] = "streamfunction_construction_gap"
        report["streamfunction"].update(
            wrapper_trace=sorted(trace),
            unknown_wrapper_count=1,
            error=type(error).__name__ + ": " + str(error),
        )
        report["elapsed_seconds"] = time.perf_counter() - started
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report

    # The frozen wave support lies strictly inside the axial/radial plateau.
    # Compare both mean and full field against their unlocalized counterparts.
    raw_mean_wave = mean.fields(reference_points, tau)
    loc_mean_wave = localized.fields(reference_points, tau)
    # Build the original local field by replacing only the mean wrapper.
    original_mean, _ = load_saved_field()
    install_in_field(original_mean)
    raw_geometry = snapshot["inputs"]["wave"]
    raw_control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    raw_wave = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    raw_wave = raw_wave[:, 0] + 1j * raw_wave[:, 1]
    raw_legacy = np.r_[np.zeros(36), raw_control[64:136], np.zeros(72)]
    raw_derivatives, raw_pressures = _unpack_full(raw_legacy, 9)
    original = LocalPotentialField(
        original_mean,
        raw_geometry["center"],
        raw_geometry["widths"],
        2,
        {0: np.zeros(2), 1: np.asarray(raw_geometry["carrier"]), 2: 2.0 * np.asarray(raw_geometry["carrier"])},
        (np.zeros(27, complex), raw_wave, np.zeros(27, complex)),
        raw_derivatives,
        raw_pressures,
        tau,
    )
    original = Mode0TangentCorrection(original, raw_geometry, raw_control[:64], tau)
    raw_dv, raw_dp = _decode_mode_block(raw_control[-128:], 3)
    original = Mode2TangentCorrection(
        original, raw_geometry["center"], raw_geometry["widths"],
        np.asarray(raw_geometry["carrier"]), 3, raw_dv, raw_dp, tau,
    )
    original.nu = mean.nu
    # Check the *full* rectangular support of the frozen wave, including all
    # four (r,z) corners.  The wave's axial cutoff is defined in physical z,
    # while the localization plateau is defined in similarity eta.
    support_r = np.asarray(
        [center[0] - widths[0], center[0] + widths[0]], dtype=float
    )
    support_z = np.asarray(
        [center[1] - widths[1], center[1] + widths[1]], dtype=float
    )
    support_points = np.asarray(
        [[radius, 0.0, height] for radius in support_r for height in support_z],
        dtype=float,
    )
    support_coord = coordinates(
        np.hypot(support_points[:, 0], support_points[:, 1]) / math.sqrt(localized.nu),
        support_points[:, 2] / math.sqrt(localized.nu),
        tau,
        localized.inner.h,
    )
    support_eta = np.asarray(support_coord["eta"], dtype=float)
    wave_eta_min = float(np.min(support_eta))
    wave_eta_max = float(np.max(support_eta))
    wave_support_inside_eta_flat = bool(
        np.max(np.abs(support_eta)) <= localized.eta_flat + 1.0e-14
    )
    wave_support_inside_radial_plateau = bool(
        np.max(np.abs(support_r)) <= r1 + 1.0e-14
    )

    # Include support-edge probes in addition to interior probes.  Boundary
    # values of a compact basis vanish, but the support geometry still must be
    # contained in the plateau for arbitrary admissible coefficients.
    plateau_points = np.vstack(
        (
            reference_points,
            support_points,
            np.asarray(
                [[center[0] + widths[0], 0.0, center[1] - 0.999 * widths[1]],
                 [center[0] + widths[0], 0.0, center[1] + 0.999 * widths[1]],
                 [center[0] - widths[0], 0.0, center[1] - 0.999 * widths[1]],
                 [center[0] - widths[0], 0.0, center[1] + 0.999 * widths[1]]],
                dtype=float,
            ),
        )
    )
    old_mean_u_plateau, old_mean_p_plateau = mean.fields(plateau_points, tau)
    new_mean_u_plateau, new_mean_p_plateau = localized.fields(plateau_points, tau)
    old_u_plateau, old_p_plateau = original.fields(plateau_points, tau)
    new_u_plateau, new_p_plateau = full.fields(plateau_points, tau)
    report["plateau_check"] = {
        "point_count": int(len(plateau_points)),
        "mean_velocity_max_abs_difference": float(np.max(np.abs(old_mean_u_plateau - new_mean_u_plateau))),
        "mean_pressure_max_abs_difference": float(np.max(np.abs(old_mean_p_plateau - new_mean_p_plateau))),
        "full_velocity_max_abs_difference": float(np.max(np.abs(old_u_plateau - new_u_plateau))),
        "full_pressure_max_abs_difference": float(np.max(np.abs(old_p_plateau - new_p_plateau))),
        "wave_patch_inside_eta_flat": wave_support_inside_eta_flat,
        "wave_patch_inside_radial_plateau": wave_support_inside_radial_plateau,
        "full_support_probe_points": support_points.tolist(),
        "full_support_eta_range": [wave_eta_min, wave_eta_max],
        "full_support_radial_range": [float(np.min(support_r)), float(np.max(support_r))],
        "full_support_plateau_containment": bool(
            wave_support_inside_eta_flat and wave_support_inside_radial_plateau
        ),
    }
    report["support"].update(
        wave_full_support_eta_range=[wave_eta_min, wave_eta_max],
        wave_full_support_radial_range=[float(np.min(support_r)), float(np.max(support_r))],
        wave_full_support_inside_eta_flat=wave_support_inside_eta_flat,
        wave_full_support_inside_radial_plateau=wave_support_inside_radial_plateau,
    )

    # Explicitly test the zero side of both cutoffs, including the exact wave
    # wrapper.  The source mean is never called in this branch.
    outside_eta = mean.inner.from_similarity([0.54], [0.60], tau)
    outside_eta2 = mean.inner.from_similarity([0.54], [-0.60], tau)
    outside_radial = np.array([[1.10 * r2, 0.0, 0.0]])
    outside_points = np.vstack((outside_eta, outside_eta2, outside_radial))
    outside_u, outside_p = full.fields(outside_points, tau)
    report["off_support_check"] = {
        "point_count": int(len(outside_points)),
        "points": outside_points.tolist(),
        "velocity_max_abs": float(np.max(np.abs(outside_u))),
        "pressure_max_abs": float(np.max(np.abs(outside_p))),
        "mean_and_wave_zero_outside_cutoffs": bool(
            np.max(np.abs(outside_u)) < 1.0e-14
            and np.max(np.abs(outside_p)) < 1.0e-14
        ),
    }

    # Sample divergence at a plateau point and in both smooth cutoff collars.
    collar_points = np.vstack(
        (
            mean.inner.from_similarity([0.54], [0.35], tau),
            mean.inner.from_similarity([0.54], [-0.35], tau),
            np.array([[(r1 + r2) / 2.0, 0.0, 0.0]]),
            np.array([[center[0], 0.0, center[1]]]),
        )
    )
    h0 = 5.0e-4 * math.sqrt(localized.nu * tau)
    divergence_rows = []
    for factor in (4.0, 2.0, 1.0):
        h = factor * h0
        mean_div = _divergence_fd(localized, collar_points, tau, h)
        full_div = _divergence_fd(full, collar_points, tau, h)
        divergence_rows.append(
            {
                "space_step": h,
                "mean_max_abs": float(np.max(np.abs(mean_div))),
                "full_max_abs": float(np.max(np.abs(full_div))),
                "mean_values": mean_div.tolist(),
                "full_values": full_div.tolist(),
            }
        )
    report["divergence_check"] = {
        "points": collar_points.tolist(),
        "rows": divergence_rows,
        "analytic_mean_divergence_zero": True,
        "wave_basis_divergence_identity": "supported_fourier_basis exact curl",
    }

    # A small full-momentum sample records cutoff-induced residuals.  This is
    # intentionally diagnostic; no forcing is introduced and no acceptance
    # threshold is inferred from it.
    residual_points = np.vstack(
        (
            mean.inner.from_similarity([0.54], [0.30], tau),
            mean.inner.from_similarity([0.54], [0.42], tau),
            np.array([[(r1 + r2) / 2.0, 0.0, 0.0]]),
            np.array([[0.75 * r2, 0.0, 0.0]]),
        )
    )
    hspace = 5.0e-4 * math.sqrt(localized.nu * tau)
    htime = 1.0e-4 * tau
    residual, divergence = independent_fd(full, residual_points, tau, hspace, htime)
    report["cutoff_residual_samples"] = {
        "points": residual_points.tolist(),
        "space_step": hspace,
        "time_step": htime,
        "momentum_norms": np.linalg.norm(residual, axis=1).tolist(),
        "momentum_max": float(np.max(np.linalg.norm(residual, axis=1))),
        "divergence": divergence.tolist(),
        "cutoff_pressure_gradient_included": True,
        "forcing": "none",
    }
    report["limitations"] = {
        "cutoff_geometry": (
            "This bounded prototype uses a similarity eta cutoff and a tau-dependent "
            "radial cutoff. It is therefore not the paper's fixed spatial c(r^2,z) "
            "terminal localization through a critical-time neighborhood."
        ),
        "pressure": (
            "Pressure is extended as chi*rho*p_local and its cutoff-gradient terms "
            "are included in the diagnostic residual; no pressure Poisson or "
            "Navier--Stokes closure is asserted."
        ),
        "acceptance": (
            "Finite support, plateau equality, and finite-difference divergence "
            "convergence are construction evidence only; PDE, finite-time, and "
            "scale-recursion acceptance remain unverified."
        ),
    }
    report["global_field_ready"] = True
    report["status"] = "completed"
    report["elapsed_seconds"] = time.perf_counter() - started
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(output_path),
                "plateau_velocity_error": report["plateau_check"]["full_velocity_max_abs_difference"],
                "off_support_velocity": report["off_support_check"]["velocity_max_abs"],
                "divergence_finest": report["divergence_check"]["rows"][-1]["full_max_abs"],
                "cutoff_residual_max": report["cutoff_residual_samples"]["momentum_max"],
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=CANDIDATE_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.candidate, args.output)
