"""Strict local structure identities for the Eq. (4.5) visualization backbone.

These identities concern only the public Cartesian mixing and similarity
coordinates. They do not identify the unknown profiles, validate Navier--Stokes,
or imply paper-exact/OpenAI-field recovery.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as s


def build_eq45_structure_receipt() -> dict:
    x, y = s.symbols("x y", real=True)
    q = s.symbols("q", positive=True)
    h = s.symbols("h", real=True)
    eta = s.symbols("eta", real=True)
    v0, F, U = s.symbols("v0 F U", finite=True)

    swirl_scale = q ** (-1 - h)
    axial_scale = q ** (-s.Rational(1, 2) - h)
    u = x * v0 / (2 * q) - y * swirl_scale * F
    v = y * v0 / (2 * q) + x * swirl_scale * F
    w = axial_scale * U
    r2 = x * x + y * y

    radial_numerator = s.simplify(x * u + y * v - r2 * v0 / (2 * q))
    swirl_numerator = s.simplify(-y * u + x * v - r2 * swirl_scale * F)
    axis_value = tuple(s.simplify(expr.subs({x: 0, y: 0})) for expr in (u, v, w))

    angle = s.symbols("angle", real=True)
    c, sn = s.cos(angle), s.sin(angle)
    xr = c * x - sn * y
    yr = sn * x + c * y
    ur = s.simplify(xr * v0 / (2 * q) - yr * swirl_scale * F)
    vr = s.simplify(yr * v0 / (2 * q) + xr * swirl_scale * F)
    rotated_original = (
        s.simplify(c * u - sn * v),
        s.simplify(sn * u + c * v),
        w,
    )
    rotation_defect = tuple(
        s.trigsimp(s.simplify(a - b))
        for a, b in zip((ur, vr, w), rotated_original)
    )

    D = s.Rational(1, 2) - h
    z_from_eta = eta * q**D
    coordinate_defect = s.simplify(
        q - z_from_eta**2 * q ** (2 * h) - q * (1 - eta**2)
    )

    assert radial_numerator == 0
    assert swirl_numerator == 0
    assert axis_value[:2] == (0, 0)
    assert rotation_defect == (0, 0, 0)
    assert coordinate_defect == 0

    return {
        "schema": "eq45_structure_receipt_v1",
        "status": "four_strict_local_identities_verified",
        "engine": f"SymPy {s.__version__}",
        "claim_scope": "eq45_leading_backbone_structure_only",
        "identities": [
            {
                "id": "cylindrical_radial_projection",
                "result": str(radial_numerator),
                "statement": "x*u+y*v=(x^2+y^2)*v0/(2*q)",
                "assumptions": "q>0; finite profile values; Cartesian Eq. (4.5) mixing",
            },
            {
                "id": "cylindrical_swirl_projection",
                "result": str(swirl_numerator),
                "statement": "-y*u+x*v=(x^2+y^2)*q^(-1-h)*F",
                "assumptions": "q>0; finite profile values; Cartesian Eq. (4.5) mixing",
            },
            {
                "id": "axis_regularity_of_cartesian_mixing",
                "result": [str(value) for value in axis_value],
                "statement": "at x=y=0, (u,v,w)=(0,0,q^(-1/2-h)*U)",
                "assumptions": "q>0 and finite U; profile provider itself is finite on X=0",
            },
            {
                "id": "rotation_equivariance_and_coordinate_identity",
                "result": {
                    "rotation_defect": [str(value) for value in rotation_defect],
                    "coordinate_defect": str(coordinate_defect),
                },
                "statement": "the transverse velocity rotates with (x,y), and z=eta*q^(1/2-h) turns q-z^2*q^(2h)=tau into tau=q*(1-eta^2)",
                "assumptions": "profiles depend only on rotation-invariant X and eta; q>0",
            },
        ],
        "states": {
            "visualization_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
        "limitations": "Local algebra for the leading Eq. (4.5) backbone only. No claim about the unknown numerical profiles, physical-space taper, momentum equation, visual correspondence, singularity, or blow-up.",
    }


def write_eq45_structure_receipt(path: str | Path) -> dict:
    result = build_eq45_structure_receipt()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    write_eq45_structure_receipt("artifacts/constrained/eq45_structure_receipt.json")
