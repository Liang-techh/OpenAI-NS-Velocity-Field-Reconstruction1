"""Spatial-channel audit for the energy-capped bipolar candidate.

This is a diagnostic only.  It samples frozen public velocity candidates and
measures where radial, swirl and axial kinetic energy live; it does not tune
or promote a candidate.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .eq45_supported_delivery import Eq45SupportedDeliveryField


REGION_NAMES = ("core", "inner_plateau", "radial_collar", "axial_collar", "corner")


def _midpoint_grid(n: int, half_width: float = 2.0):
    dx = 2.0 * half_width / n
    a = np.linspace(-half_width + 0.5 * dx, half_width - 0.5 * dx, n)
    x, y, z = np.meshgrid(a, a, a, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    return points, dx


def _cylindrical(points: np.ndarray, velocity: np.ndarray):
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    safe = np.where(r > 0.0, r, 1.0)
    ur = (x * velocity[:, 0] + y * velocity[:, 1]) / safe
    ut = (-y * velocity[:, 0] + x * velocity[:, 1]) / safe
    ur = np.where(r > 0.0, ur, 0.0)
    ut = np.where(r > 0.0, ut, 0.0)
    return r, ur, ut, velocity[:, 2]


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    weights = np.asarray(weights, dtype=float)
    if not np.any(weights > 0.0):
        return 0.0
    order = np.argsort(values)
    v = np.asarray(values, dtype=float)[order]
    w = weights[order]
    c = np.cumsum(w)
    idx = int(np.searchsorted(c, q * c[-1], side="left"))
    return float(v[min(idx, len(v) - 1)])


def _regions(r: np.ndarray, z: np.ndarray):
    az = np.abs(z)
    core = (r < 0.8) & (az < 0.8)
    plateau = (r < 1.6) & (az < 1.6)
    return {
        "core": core,
        "inner_plateau": plateau & ~core,
        "radial_collar": (r >= 1.6) & (az < 1.6),
        "axial_collar": (r < 1.6) & (az >= 1.6),
        "corner": (r >= 1.6) & (az >= 1.6),
    }


def _metrics(field, points: np.ndarray, dx: float, time: float):
    vel = field.at_points(points, time)
    r, ur, ut, uz = _cylindrical(points, vel)
    az = np.abs(points[:, 2])
    dv = dx**3
    sq = {"radial": ur**2, "swirl": ut**2, "axial": uz**2}
    energy = {name: float(0.5 * np.sum(v) * dv) for name, v in sq.items()}
    energy["total"] = float(sum(energy.values()))
    masks = _regions(r, points[:, 2])
    region = {}
    for name, mask in masks.items():
        channel = {k: float(0.5 * np.sum(v[mask]) * dv) for k, v in sq.items()}
        channel["total"] = float(sum(channel.values()))
        region[name] = channel
    swirl_w = sq["swirl"]
    poloidal_w = sq["radial"] + sq["axial"]
    return {
        "energy": energy,
        "channel_fraction": {k: energy[k] / energy["total"] for k in ("radial", "swirl", "axial")},
        "regions": region,
        "region_fraction_of_total": {k: region[k]["total"] / energy["total"] for k in REGION_NAMES},
        "angular_momentum_rms": float(np.sqrt(np.mean((r * ut) ** 2))),
        "swirl_concentration": {
            "radial_q50": _weighted_quantile(r, swirl_w, 0.50),
            "radial_q90": _weighted_quantile(r, swirl_w, 0.90),
            "axial_q50": _weighted_quantile(az, swirl_w, 0.50),
            "axial_q90": _weighted_quantile(az, swirl_w, 0.90),
        },
        "poloidal_concentration": {
            "radial_q50": _weighted_quantile(r, poloidal_w, 0.50),
            "radial_q90": _weighted_quantile(r, poloidal_w, 0.90),
            "axial_q50": _weighted_quantile(az, poloidal_w, 0.50),
            "axial_q90": _weighted_quantile(az, poloidal_w, 0.90),
        },
    }


def _relative_change(new: float, old: float) -> float:
    return float((new - old) / old) if old != 0.0 else 0.0


def run(
    resolutions=(24, 32, 40),
    times=(0.3125, 0.5, 0.75),
    output: str | Path | None = None,
):
    odd = Eq45SupportedDeliveryField.load_candidate("artifacts/bipolar_joint_odd3/candidate.json")
    capped = Eq45SupportedDeliveryField.load_candidate("artifacts/bipolar_joint_capped/candidate.json")
    ob = odd.candidate.parent.profile_basis
    cb = capped.candidate.parent.profile_basis
    rows = []
    for n in resolutions:
        points, dx = _midpoint_grid(int(n))
        for t in times:
            om = _metrics(odd, points, dx, float(t))
            cm = _metrics(capped, points, dx, float(t))
            rows.append({"resolution": int(n), "time": float(t), "odd3": om, "capped": cm})

    fine_n = int(resolutions[-1])
    fine = [r for r in rows if r["resolution"] == fine_n]
    comparison = []
    for row in fine:
        o, c = row["odd3"], row["capped"]
        comparison.append({
            "time": row["time"],
            "total_energy_relative_change": _relative_change(c["energy"]["total"], o["energy"]["total"]),
            "radial_energy_relative_change": _relative_change(c["energy"]["radial"], o["energy"]["radial"]),
            "swirl_energy_relative_change": _relative_change(c["energy"]["swirl"], o["energy"]["swirl"]),
            "axial_energy_relative_change": _relative_change(c["energy"]["axial"], o["energy"]["axial"]),
            "radial_collar_total_relative_change": _relative_change(c["regions"]["radial_collar"]["total"], o["regions"]["radial_collar"]["total"]),
            "axial_collar_total_relative_change": _relative_change(c["regions"]["axial_collar"]["total"], o["regions"]["axial_collar"]["total"]),
            "corner_total_relative_change": _relative_change(c["regions"]["corner"]["total"], o["regions"]["corner"]["total"]),
            "angular_momentum_rms_relative_change": _relative_change(c["angular_momentum_rms"], o["angular_momentum_rms"]),
        })

    probe = np.array([[0.1 * np.sqrt(0.75), 0.0, 0.1 * 0.75**0.495]])
    center = []
    for t in times:
        ov = odd.at_points(probe, float(t))[0]
        cv = capped.at_points(probe, float(t))[0]
        _, our, out, ouz = _cylindrical(probe, ov[None, :])
        _, cur, cut, cuz = _cylindrical(probe, cv[None, :])
        center.append({
            "time": float(t),
            "odd3": {"radial": float(our[0]), "swirl": float(out[0]), "axial": float(ouz[0])},
            "capped": {"radial": float(cur[0]), "swirl": float(cut[0]), "axial": float(cuz[0])},
        })

    medium_n = int(resolutions[-2])
    stability = []
    for t in times:
        med = next(r for r in rows if r["resolution"] == medium_n and r["time"] == float(t))
        fin = next(r for r in rows if r["resolution"] == fine_n and r["time"] == float(t))
        entry = {"time": float(t)}
        for candidate in ("odd3", "capped"):
            entry[candidate] = {
                "total_energy_relative_change": abs(_relative_change(fin[candidate]["energy"]["total"], med[candidate]["energy"]["total"])),
                "angular_momentum_rms_relative_change": abs(_relative_change(fin[candidate]["angular_momentum_rms"], med[candidate]["angular_momentum_rms"])),
            }
        stability.append(entry)

    report = {
        "scope": "Public-[u,v,w] spatial-channel redistribution audit; diagnostic only, no coefficient selection or PDE/visual promotion.",
        "resolutions": [int(n) for n in resolutions],
        "times": [float(t) for t in times],
        "odd3_sha256": odd.sha256,
        "capped_sha256": capped.sha256,
        "phi_coefficient_max_abs_delta": float(np.max(np.abs(np.asarray(cb.phi_coefficients) - np.asarray(ob.phi_coefficients)))),
        "swirl_coefficient_max_abs_delta": float(np.max(np.abs(np.asarray(cb.swirl_coefficients) - np.asarray(ob.swirl_coefficients)))),
        "rows": rows,
        "finest_comparison": comparison,
        "central_signed_components": center,
        "medium_to_fine_stability": stability,
        "truth_boundary": {"pde_validated": False, "visualization_ready": False, "visual_correspondence_verified": False},
    }
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    import sys

    result = run(output=sys.argv[1] if len(sys.argv) > 1 else None)
    print(json.dumps(result, indent=2))
