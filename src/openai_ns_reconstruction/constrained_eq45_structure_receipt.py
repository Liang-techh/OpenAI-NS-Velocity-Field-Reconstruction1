"""Exact local structure certificates for the Eq. (4.5) visualization backbone.

The checks in this module use only exact rational polynomial arithmetic from the
Python standard library. They concern the public Cartesian mixing and the
similarity-coordinate exponent algebra. They do not identify the unknown
profiles, validate Navier--Stokes, or imply paper-exact/OpenAI-field recovery.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

# Exact polynomials in formal variables (a,b,c,s,x,y). Here
# a=v0/(2q), b=q^(-1-h)F, and (c,s) is an arbitrary planar rotation pair.
_NVAR = 6
_ZERO_EXP = (0,) * _NVAR


def _clean(poly):
    return {exp: coeff for exp, coeff in poly.items() if coeff}


def _const(value):
    value = Fraction(value)
    return {} if value == 0 else {_ZERO_EXP: value}


def _var(index):
    exp = [0] * _NVAR
    exp[index] = 1
    return {tuple(exp): Fraction(1)}


def _add(left, right):
    out = dict(left)
    for exp, coeff in right.items():
        out[exp] = out.get(exp, Fraction(0)) + coeff
    return _clean(out)


def _neg(poly):
    return {exp: -coeff for exp, coeff in poly.items()}


def _sub(left, right):
    return _add(left, _neg(right))


def _mul(left, right):
    out = {}
    for exp_l, coeff_l in left.items():
        for exp_r, coeff_r in right.items():
            exp = tuple(a + b for a, b in zip(exp_l, exp_r))
            out[exp] = out.get(exp, Fraction(0)) + coeff_l * coeff_r
    return _clean(out)


def _matmul(left, right):
    rows = len(left)
    inner = len(left[0])
    cols = len(right[0])
    result = [[{} for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            value = {}
            for k in range(inner):
                value = _add(value, _mul(left[i][k], right[k][j]))
            result[i][j] = value
    return result


def _mat_sub(left, right):
    return [
        [_sub(a, b) for a, b in zip(row_l, row_r)]
        for row_l, row_r in zip(left, right)
    ]


def _zero_matrix(matrix):
    return all(not entry for row in matrix for entry in row)


def _exact_certificate():
    a, b, c, sn, x, y = (_var(i) for i in range(_NVAR))
    zero = _const(0)

    # [u,v]^T = [[a,-b],[b,a]] [x,y]^T.
    mixing = [[a, _neg(b)], [b, a]]
    uv = _matmul(mixing, [[x], [y]])
    u, v = uv[0][0], uv[1][0]
    r2 = _add(_mul(x, x), _mul(y, y))

    radial_defect = _sub(_add(_mul(x, u), _mul(y, v)), _mul(a, r2))
    swirl_defect = _sub(
        _add(_mul(_neg(y), u), _mul(x, v)), _mul(b, r2)
    )

    # u and v each retain at least one factor x or y, hence vanish exactly on
    # the symmetry axis without evaluating 1/r.
    axis_transverse_zero = all(
        all(exp[4] + exp[5] > 0 for exp in component)
        for component in (u, v)
    )

    # M=aI+bJ commutes with every planar rotation R(c,s).
    rotation = [[c, _neg(sn)], [sn, c]]
    rotation_commutator = _mat_sub(
        _matmul(mixing, rotation), _matmul(rotation, mixing)
    )

    # D=1/2-h. Store q exponents as (constant part, h coefficient):
    # 2D+2h = 1 exactly, so z^2 q^(2h)=eta^2 q when z=eta q^D.
    D = (Fraction(1, 2), Fraction(-1))
    twice_D_plus_2h = (2 * D[0], 2 * D[1] + 2)
    coordinate_exponent_identity = twice_D_plus_2h == (
        Fraction(1),
        Fraction(0),
    )

    assert radial_defect == zero
    assert swirl_defect == zero
    assert axis_transverse_zero
    assert _zero_matrix(rotation_commutator)
    assert coordinate_exponent_identity

    return {
        "radial_defect": "0",
        "swirl_defect": "0",
        "axis_transverse": ["0", "0"],
        "rotation_commutator": [["0", "0"], ["0", "0"]],
        "coordinate_q_exponent": "1",
    }


def build_eq45_structure_receipt() -> dict:
    certificate = _exact_certificate()
    return {
        "schema": "eq45_structure_receipt_v1",
        "status": "four_strict_local_identities_verified",
        "engine": "exact stdlib Fraction polynomial certificate",
        "claim_scope": "eq45_leading_backbone_structure_only",
        "identities": [
            {
                "id": "cylindrical_radial_projection",
                "result": certificate["radial_defect"],
                "statement": "x*u+y*v=(x^2+y^2)*v0/(2*q)",
                "assumptions": "q>0; finite profile values; Cartesian Eq. (4.5) mixing",
            },
            {
                "id": "cylindrical_swirl_projection",
                "result": certificate["swirl_defect"],
                "statement": "-y*u+x*v=(x^2+y^2)*q^(-1-h)*F",
                "assumptions": "q>0; finite profile values; Cartesian Eq. (4.5) mixing",
            },
            {
                "id": "axis_regularity_of_cartesian_mixing",
                "result": ["0", "0", "q^(-1/2-h)*U"],
                "statement": "at x=y=0, (u,v,w)=(0,0,q^(-1/2-h)*U)",
                "assumptions": "q>0 and finite U; profile provider itself is finite on X=0",
            },
            {
                "id": "rotation_equivariance_and_coordinate_identity",
                "result": {
                    "rotation_commutator": certificate["rotation_commutator"],
                    "coordinate_q_exponent": certificate["coordinate_q_exponent"],
                    "coordinate_defect": "0",
                },
                "statement": "the transverse mixing commutes with every planar rotation, and z=eta*q^(1/2-h) turns q-z^2*q^(2h)=tau into tau=q*(1-eta^2)",
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
        "limitations": "Local exact algebra for the leading Eq. (4.5) backbone only. No claim about the unknown numerical profiles, physical-space taper, momentum equation, visual correspondence, singularity, or blow-up.",
    }


def write_eq45_structure_receipt(path: str | Path) -> dict:
    result = build_eq45_structure_receipt()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    write_eq45_structure_receipt("artifacts/constrained/eq45_structure_receipt.json")
