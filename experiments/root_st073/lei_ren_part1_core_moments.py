"""Finite local five-moment receipt for the reusable PaperCoreReference.

This script wires the finite-interval Lei--Ren Part I moment interface to the
repository's scalar ``PaperCoreReference`` profile.  The reference profile is
vectorized here because its public ``F``, ``U`` and ``Pi`` methods intentionally
accept scalar coordinates.  Parameters and the finite near-axis interval are
explicit and autonomous; the result is a local numerical receipt, not the
nonlinear Part I field or an outer/heat-exterior construction.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = Path(__file__).resolve().parent
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
from lei_ren_part1_moments import (
    MOMENT_NAMES,
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
    cumulative_five_moments,
)


def _vectorize_profile(function):
    """Adapt a scalar profile method to NumPy arrays."""

    return np.vectorize(function, otypes=[float])


def _array_dict(values: dict[str, np.ndarray]) -> dict[str, list[float]]:
    """Convert result arrays to JSON-native lists while preserving all keys."""

    return {name: np.asarray(value, dtype=float).tolist() for name, value in values.items()}


def run() -> dict[str, object]:
    reference = PaperCoreReference(
        h=0.001,
        j=0.02,
        sigma=0.3,
        Lambda=10.0,
        C=2.0,
        pressure_scale=1.0,
    )
    rmax = 0.05
    z = np.linspace(-0.5, 0.5, 9)

    F = _vectorize_profile(reference.F)
    Uz = _vectorize_profile(reference.U)
    P0 = _vectorize_profile(lambda eta: reference.Pi(0.0, eta))(z)
    Pedge = _vectorize_profile(lambda eta: reference.Pi(rmax, eta))(z)

    results = {}
    for order in (32, 64):
        results[order] = cumulative_five_moments(
            F,
            Uz,
            rmax,
            z,
            n=order,
            P0=P0,
            P=Pedge,
        )

    refinement = {
        name: (
            np.asarray(results[64][name], dtype=float)
            - np.asarray(results[32][name], dtype=float)
        ).tolist()
        for name in (*MOMENT_NAMES, "pressure_defect")
    }
    pressure_defect_max = {
        str(order): float(np.max(np.abs(results[order]["pressure_defect"])))
        for order in (32, 64)
    }
    metadata = reference.metadata()
    metadata.update(
        {
            "autonomous_pressure": True,
            "pressure_comparison": (
                "P0(Z)=Pi(0,Z), P(Rmax,Z)=Pi(0.05,Z); both are supplied by "
                "the same finite reference profile"
            ),
            "nonlinear_part1_field": False,
            "LR1_04_status": (
                "partial: finite local profile moments are wired; outer/heat "
                "profile binding and global terminal matching remain open"
            ),
            "outer_global_certificate": False,
            "morphology_claim": False,
        }
    )
    report = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "profile": "openai_ns_reconstruction.paper_core_reference.PaperCoreReference",
        "parameters": {
            "h": 0.001,
            "j": 0.02,
            "sigma": 0.3,
            "Lambda": 10.0,
            "C": 2.0,
            "pressure_scale": 1.0,
        },
        "Rmax": rmax,
        "Z": z.tolist(),
        "quadrature_orders": [32, 64],
        "moments_n32": _array_dict(results[32]),
        "moments_n64": _array_dict(results[64]),
        "refinement_n64_minus_n32": refinement,
        "pressure_defect_max_abs": pressure_defect_max,
        "profile_metadata": metadata,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Finite local near-axis PaperCoreReference moment receipt on "
            "0<=R<=0.05 and |Z|<=0.5. The pressure datum is autonomous and "
            "unmatched to an exterior. This is not the nonlinear Part I field; "
            "no outer/global, heat-exterior, PDE, or morphology certificate "
            "is claimed."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()
