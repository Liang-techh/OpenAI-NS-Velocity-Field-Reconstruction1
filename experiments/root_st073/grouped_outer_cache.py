"""Grouped radial cache for the outer pressure stress-cone samples.

This is a drop-in replacement for :func:`outer_pressure_modes.outer_cache`.
Centers at the same physical axial coordinate share one radial quadrature.
The shared quadrature is cut at every target radius and at the moving support
cuts, so a target panel remains a prefix of its group's points.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from affine_momentum import jets
from joined_field import coordinates
from outer_pressure_modes import BREAKS, WINDOWS


def _default_locations():
    return [(eta, y) for eta in (-0.2, 0.0, 0.2) for y in (0.5, 0.75)]


def _group_indices_by_z(centers):
    """Return insertion-ordered groups of center indices with equal physical z."""

    groups = {}
    for index, center in enumerate(centers):
        # ``from_similarity`` is deterministic for a fixed eta.  Keeping the
        # float key exact avoids merging nearby eta-neighborhood samples.
        key = float(center[2])
        groups.setdefault(key, []).append(index)
    return tuple(groups.values())


def _group_quadrature(centers, group, inner, tau, order, points):
    """Append one group's radial nodes and return its per-center panels."""

    radii = np.asarray([float(centers[index][0]) for index in group], dtype=float)
    max_radius = float(np.max(radii))
    z = float(centers[group[0]][2])
    q = float(coordinates(0.0, z / np.sqrt(inner.nu), tau, inner.h)["q"])
    ri = float(np.sqrt(2.0 * inner.nu * q * inner.p.X_max))

    # Include every target radius, the inner join, and every moving support
    # cut.  Clipping is the same convention used by outer_cache for targets
    # that lie inside or outside a support cut.
    cuts = [0.0, max_radius, ri]
    cuts.extend(radii.tolist())
    cuts.extend(ri * (1.0 + 15.0 * float(b)) for b in BREAKS)
    edges = sorted(set(float(x) for x in np.clip(cuts, 0.0, max_radius)))

    g, w = leggauss(order)
    all_radii = []
    all_weights = []
    prefix_at_edge = {}
    prefix = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        rr = (lo + hi) / 2.0 + (hi - lo) / 2.0 * g
        ww = (hi - lo) / 2.0 * w
        all_radii.append(rr)
        all_weights.append(ww)
        prefix += len(rr)
        prefix_at_edge[hi] = prefix

    if all_radii:
        all_radii = np.concatenate(all_radii)
        all_weights = np.concatenate(all_weights)
    else:
        all_radii = np.empty(0, dtype=float)
        all_weights = np.empty(0, dtype=float)

    group_start = len(points)
    points.extend((float(radius), 0.0, z) for radius in all_radii)

    panels = {}
    for index, radius in zip(group, radii):
        count = prefix_at_edge.get(float(radius), 0)
        panels[index] = (
            slice(group_start, group_start + count),
            all_radii[:count],
            all_weights[:count],
            float(radius),
        )
    return panels


def outer_cache(field, units, k=11, order=48, locations=None):
    """Build a grouped replacement for ``outer_pressure_modes.outer_cache``.

    The returned ``baseline``, ``modes``, ``indices``, ``centers`` and
    ``panels`` entries have the same meaning as the legacy helper.  Additional
    ``points``, ``point_count``, ``group_count`` and ``quadrature_point_count``
    entries expose the benchmark geometry without requiring reconstruction.
    """

    inner = field.inner
    tau = 0.5 * 2.0 ** (-k)
    locations = _default_locations() if locations is None else list(locations)
    units = list(units)

    centers = [
        np.asarray(inner.from_similarity(
            [inner.p.X_max * (1.0 + 15.0 * y) ** 2], [eta], tau
        )[0], dtype=float)
        for eta, y in locations
    ]
    points = list(centers)
    groups = _group_indices_by_z(centers)
    panels_by_center = {}
    for group in groups:
        panels_by_center.update(
            _group_quadrature(centers, group, inner, tau, order, points)
        )
    panels = [panels_by_center[index] for index in range(len(centers))]

    points = np.asarray(points, dtype=float)
    args = (
        points,
        tau,
        0.0005 * np.sqrt(inner.nu * tau),
        0.0001 * tau,
    )
    baseline = jets(field, *args)
    unit_jets = [jets(unit, *args) for unit in units]
    if unit_jets:
        modes = tuple(
            np.stack([unit_jet[index] for unit_jet in unit_jets])
            for index in range(3)
        )
    else:
        # The legacy np.stack path rejects units=[]; preserve the natural
        # affine shape so callers can still use a baseline-only cache.
        modes = tuple(
            np.empty((0,) + np.asarray(base).shape, dtype=np.asarray(base).dtype)
            for base in baseline
        )

    point_count = int(len(points))
    return dict(
        baseline=baseline,
        modes=modes,
        panels=panels,
        indices=np.arange(len(units), dtype=int),
        centers=centers,
        points=points,
        tau=float(tau),
        point_count=point_count,
        group_count=int(len(groups)),
        quadrature_point_count=point_count - len(centers),
    )


# A descriptive alias lets benchmark callers avoid shadowing the legacy name.
grouped_outer_cache = outer_cache


def benchmark(order=24, locations=None):
    """Compare legacy and grouped cone replay for the staged k=11 field."""

    from adaptive_bridge_recursive_defect import build_fields
    from midplane_connected_cone_repair import ZeroBackground, cached_cones
    from radial_continuation import ROOT
    from separated_moment_modes import SeparatedMomentModes
    from outer_pressure_modes import OuterPressure, outer_cache as legacy_outer_cache

    inner, fields = build_fields()
    seed = json.loads((Path(ROOT) / "midplane_outer_pressure_staged_k11.json").read_text())
    velocity = SeparatedMomentModes(
        fields["two_sided_cone"],
        seed["amplitudes"],
        windows=WINDOWS,
        knots=(11.0, 15.0, 19.0),
        axial_powers=(0, 1, 2),
    )
    zero = ZeroBackground(fields["two_sided_cone"])
    units = [OuterPressure(zero, np.eye(6)[0])]
    if locations is None:
        locations = _default_locations()

    started = time.perf_counter()
    legacy = legacy_outer_cache(
        velocity, units, k=11, order=order, locations=locations
    )
    legacy_seconds = time.perf_counter() - started
    started = time.perf_counter()
    grouped = outer_cache(
        velocity, units, k=11, order=order, locations=locations
    )
    grouped_seconds = time.perf_counter() - started

    zero_delta = np.zeros(1)
    unit_delta = np.ones(1)
    old_zero = cached_cones(legacy, zero_delta)
    new_zero = cached_cones(grouped, zero_delta)
    old_unit = cached_cones(legacy, unit_delta)
    new_unit = cached_cones(grouped, unit_delta)
    return dict(
        order=int(order),
        location_count=len(locations),
        group_count=int(grouped["group_count"]),
        legacy_point_count=int(legacy["baseline"][0].shape[0]),
        grouped_point_count=int(grouped["point_count"]),
        legacy_quadrature_point_count=int(
            legacy["baseline"][0].shape[0] - len(legacy["centers"])
        ),
        grouped_quadrature_point_count=int(grouped["quadrature_point_count"]),
        legacy_seconds=float(legacy_seconds),
        grouped_seconds=float(grouped_seconds),
        speedup=float(legacy_seconds / grouped_seconds)
        if grouped_seconds > 0.0 else None,
        zero_cached_cones_max_abs_diff=float(np.max(np.abs(old_zero - new_zero))),
        unit_cached_cones_max_abs_diff=float(np.max(np.abs(old_unit - new_unit))),
        zero_cached_cones_old=old_zero.tolist(),
        zero_cached_cones_grouped=new_zero.tolist(),
        unit_cached_cones_old=old_unit.tolist(),
        unit_cached_cones_grouped=new_unit.tolist(),
    )


def benchmark_locations():
    """Return the 22-location original-plus-neighborhood benchmark set."""

    locations = _default_locations()
    locations.extend(
        (eta + de, 0.75 + dy)
        for eta in (-0.2, 0.2)
        for de in (-0.005, 0.0, 0.005)
        for dy in (-0.005, 0.0, 0.005)
    )
    result = []
    seen = set()
    for location in locations:
        key = tuple(float(x) for x in location)
        if key not in seen:
            seen.add(key)
            result.append(key)
    return result


if __name__ == "__main__":
    report = benchmark(order=24, locations=benchmark_locations())
    output = Path(__file__).with_name("grouped_outer_cache_check.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
