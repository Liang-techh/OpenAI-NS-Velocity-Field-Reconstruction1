"""Source physical-evaluation differential algebra for Kokuno oscillatory fields.

The corrected 2026-09-09 reconstruction introduces an auxiliary two-torus
variable Y and only afterwards evaluates it on physical space-time,

    Y(r,t) = v_r r**d_r + v_t t  (mod Z^2).

For a smooth torus-dependent quantity F(r,theta,z,t,Y), differentiation after
this evaluation is therefore not the plain partial derivative. The source
operators are

    mathsf_r = partial_r + d_r r**(d_r-1) v_r . partial_Y,
    mathsf_t = partial_t + v_t . partial_Y,

and in the source wave chart

    D_r       = sqrt(Q) mathsf_r,
    D_z       = sqrt(Q) partial_z,
    mathsf_t* = Q**(1+h) mathsf_t.

The same corrected reader also displays the numerical torus/group constants
used by this pullback:

    J_g = [[3,1],[1,5]],  b_g = sqrt(2)-1,
    v_r = (1,-b_g),       v_t = (b_g,1),
    Lambda_g = 4-sqrt(2), T_g = 4+sqrt(2), kappa_s = 1e-5.

This module makes that differential-algebra seam executable and binds those
*displayed source constants* through :meth:`from_corrected_source_constants`.
It still does not reconstruct the positive-order background, pulse/mode labels,
partition centers, or a public source oscillatory xyz-t velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any

import numpy as np

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-physical-evaluation-differential-algebra-v2"

_SOURCE_SQRT2 = float(np.sqrt(2.0))
SOURCE_J_G = ((3.0, 1.0), (1.0, 5.0))
SOURCE_B_G = _SOURCE_SQRT2 - 1.0
SOURCE_V_R = (1.0, -SOURCE_B_G)
SOURCE_V_T = (SOURCE_B_G, 1.0)
SOURCE_LAMBDA_G = 4.0 - _SOURCE_SQRT2
SOURCE_T_G = 4.0 + _SOURCE_SQRT2
SOURCE_KAPPA_S = 1.0e-5

_SOURCE_IDENTITIES = {
    "phase_map": "Y=v_r*r^d_r+v_t*t mod Z^2",
    "mathsf_r": "partial_r+d_r*r^(d_r-1)*v_r dot partial_Y",
    "mathsf_t": "partial_t+v_t dot partial_Y",
    "D_r": "sqrt(Q)*mathsf_r",
    "D_z": "sqrt(Q)*partial_z",
    "mathsf_t_star": "Q^(1+h)*mathsf_t",
    "group_matrix": "J_g=[[3,1],[1,5]]",
    "torus_vectors": "b_g=sqrt(2)-1; v_r=(1,-b_g); v_t=(b_g,1)",
    "group_eigenvalues": "Lambda_g=4-sqrt(2); T_g=4+sqrt(2)",
    "group_relations": "J_g*v_r=Lambda_g*v_r; J_g*v_t=T_g*v_t; v_r dot v_t=0; det(J_g)=14",
    "radial_exponent": "d_r=2*((1+h)*rho_g-h*kappa_s), rho_g=log(Lambda_g)/log(T_g), kappa_s=1e-5",
}

_TRUTH_BOUNDARY = {
    "auxiliary_torus_physical_evaluation_map_executable": True,
    "source_chain_rule_operators_executable": True,
    "source_normalized_Dr_Dz_tstar_executable": True,
    "displayed_source_torus_vectors_bound": True,
    "displayed_source_group_constants_bound": True,
    "source_torus_vectors_recovered": True,
    "source_group_constants_recovered": True,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

_BINDINGS = {"caller_supplied", "corrected_source_displayed_constants"}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_real(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _two_vector(value: Any, name: str) -> np.ndarray:
    out = _finite_real(value, name)
    if out.shape != (2,):
        raise ValueError(f"{name} must be a finite two-vector")
    return out


def corrected_source_torus_constants() -> dict[str, Any]:
    """Return the auxiliary-torus constants displayed by the corrected reader.

    The returned numerical values are source-bound, not autonomous fit
    parameters. Diagnostics are recomputed from the displayed constants so a
    caller can fail closed if the binding is edited inconsistently.
    """
    jg = np.asarray(SOURCE_J_G, dtype=float)
    vr = np.asarray(SOURCE_V_R, dtype=float)
    vt = np.asarray(SOURCE_V_T, dtype=float)
    eigen_r = jg @ vr - SOURCE_LAMBDA_G * vr
    eigen_t = jg @ vt - SOURCE_T_G * vt
    return {
        "J_g": [list(row) for row in SOURCE_J_G],
        "b_g": SOURCE_B_G,
        "v_r": list(SOURCE_V_R),
        "v_t": list(SOURCE_V_T),
        "Lambda_g": SOURCE_LAMBDA_G,
        "T_g": SOURCE_T_G,
        "kappa_s": SOURCE_KAPPA_S,
        "identity_diagnostics": {
            "Jg_vr_minus_Lambda_vr_max_abs": float(np.max(np.abs(eigen_r))),
            "Jg_vt_minus_T_vt_max_abs": float(np.max(np.abs(eigen_t))),
            "vr_dot_vt": float(vr @ vt),
            "det_J_g": float(np.linalg.det(jg)),
        },
    }


@dataclass(frozen=True)
class KokunoPhysicalEvaluationMap:
    """Execute the source torus pullback and post-evaluation derivatives.

    Use :meth:`from_corrected_source_constants` to bind the numerical torus
    directions and group constants displayed by the corrected reader. Direct
    construction remains available for independent manufactured tests or other
    caller-supplied source-compatible realizations, with provenance recording
    that distinction explicitly.
    """

    v_r: tuple[float, float]
    v_t: tuple[float, float]
    d_r: float
    h: float
    binding: str = "caller_supplied"

    def __post_init__(self) -> None:
        vr = _two_vector(self.v_r, "v_r")
        vt = _two_vector(self.v_t, "v_t")
        d_r = float(self.d_r)
        h = float(self.h)
        binding = str(self.binding)
        if not np.isfinite(d_r) or d_r <= 0.0:
            raise ValueError("d_r must be finite and positive")
        if not np.isfinite(h) or not (0.0 < h < 0.01):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        if binding not in _BINDINGS:
            raise ValueError(f"binding must be one of {sorted(_BINDINGS)}")
        object.__setattr__(self, "v_r", (float(vr[0]), float(vr[1])))
        object.__setattr__(self, "v_t", (float(vt[0]), float(vt[1])))
        object.__setattr__(self, "d_r", d_r)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "binding", binding)

    @staticmethod
    def source_radial_exponent(
        lambda_g: float,
        t_g: float,
        h: float,
        *,
        kappa_s: float = SOURCE_KAPPA_S,
    ) -> dict[str, float]:
        """Evaluate the corrected-reader formula for ``rho_g`` and ``d_r``."""
        lam = float(lambda_g)
        tg = float(t_g)
        hh = float(h)
        ks = float(kappa_s)
        if not all(np.isfinite(x) for x in (lam, tg, hh, ks)):
            raise ValueError("source group constants must be finite")
        if lam <= 0.0 or tg <= 0.0 or tg == 1.0:
            raise ValueError("lambda_g and t_g must be positive and t_g must differ from 1")
        if not (0.0 < hh < 0.01):
            raise ValueError("h must satisfy 0<h<1/100")
        if ks < 0.0:
            raise ValueError("kappa_s must be nonnegative")
        rho_g = float(np.log(lam) / np.log(tg))
        d_r = float(2.0 * ((1.0 + hh) * rho_g - hh * ks))
        if not np.isfinite(d_r) or d_r <= 0.0:
            raise ValueError("source radial-exponent formula must produce d_r>0")
        return {"rho_g": rho_g, "d_r": d_r, "kappa_s": ks}

    @classmethod
    def from_source_group_constants(
        cls,
        *,
        v_r: tuple[float, float],
        v_t: tuple[float, float],
        lambda_g: float,
        t_g: float,
        h: float,
        kappa_s: float = SOURCE_KAPPA_S,
        binding: str = "caller_supplied",
    ) -> "KokunoPhysicalEvaluationMap":
        exponent = cls.source_radial_exponent(lambda_g, t_g, h, kappa_s=kappa_s)
        return cls(v_r=v_r, v_t=v_t, d_r=exponent["d_r"], h=h, binding=binding)

    @classmethod
    def from_corrected_source_constants(cls, *, h: float) -> "KokunoPhysicalEvaluationMap":
        """Bind the corrected reader's displayed torus/group constants."""
        return cls.from_source_group_constants(
            v_r=SOURCE_V_R,
            v_t=SOURCE_V_T,
            lambda_g=SOURCE_LAMBDA_G,
            t_g=SOURCE_T_G,
            h=h,
            kappa_s=SOURCE_KAPPA_S,
            binding="corrected_source_displayed_constants",
        )

    @property
    def v_r_array(self) -> np.ndarray:
        return np.asarray(self.v_r, dtype=float)

    @property
    def v_t_array(self) -> np.ndarray:
        return np.asarray(self.v_t, dtype=float)

    def torus_phase(self, r: Any, t: Any) -> np.ndarray:
        """Return ``Y(r,t)`` represented in the half-open unit square."""
        rr, tt = np.broadcast_arrays(_finite_real(r, "r"), _finite_real(t, "t"))
        if np.any(rr <= 0.0):
            raise ValueError(
                "source physical evaluation is defined for r>0; "
                "use the wave-shell guard before this map"
            )
        raw = rr[..., None] ** self.d_r * self.v_r_array + tt[..., None] * self.v_t_array
        return np.mod(raw, 1.0)

    def differentiate(
        self,
        *,
        r: Any,
        q: Any,
        partial_r: Any,
        partial_t: Any,
        partial_z: Any,
        grad_y: Any,
    ) -> dict[str, np.ndarray]:
        """Apply the exact post-evaluation chain rule and source scalings.

        All scalar derivative inputs broadcast to one sample shape. ``grad_y``
        must broadcast to that shape plus a final torus dimension of size two.
        """
        rr, qq, pr, pt, pz = np.broadcast_arrays(
            _finite_real(r, "r"),
            _finite_real(q, "q"),
            _finite_real(partial_r, "partial_r"),
            _finite_real(partial_t, "partial_t"),
            _finite_real(partial_z, "partial_z"),
        )
        if np.any(rr <= 0.0):
            raise ValueError("r must be strictly positive on the source oscillatory shell")
        if np.any(qq <= 0.0):
            raise ValueError("q must be strictly positive")
        gy = _finite_real(grad_y, "grad_y")
        try:
            gy = np.broadcast_to(gy, rr.shape + (2,))
        except ValueError:
            raise ValueError("grad_y must broadcast to sample_shape+(2,)") from None

        torus_r = np.einsum("...i,i->...", gy, self.v_r_array)
        torus_t = np.einsum("...i,i->...", gy, self.v_t_array)
        radial_weight = self.d_r * rr ** (self.d_r - 1.0)

        mathsf_r = pr + radial_weight * torus_r
        mathsf_t = pt + torus_t
        sqrt_q = np.sqrt(qq)
        d_r_normalized = sqrt_q * mathsf_r
        d_z_normalized = sqrt_q * pz
        t_star = qq ** (1.0 + self.h) * mathsf_t

        return {
            "mathsf_r": mathsf_r,
            "mathsf_t": mathsf_t,
            "partial_z_after_evaluation": pz,
            "D_r": d_r_normalized,
            "D_z": d_z_normalized,
            "mathsf_t_star": t_star,
            "radial_torus_weight": radial_weight,
            "v_r_dot_grad_y": torus_r,
            "v_t_dot_grad_y": torus_t,
        }

    def provenance(self) -> dict[str, Any]:
        source_constants = corrected_source_torus_constants()
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_identities": dict(_SOURCE_IDENTITIES),
            "source_displayed_torus_constants": source_constants,
            "numeric_binding": {
                "kind": self.binding,
                "v_r": list(self.v_r),
                "v_t": list(self.v_t),
                "d_r": self.d_r,
                "h": self.h,
                "uses_displayed_source_torus_constants": self.binding
                == "corrected_source_displayed_constants",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload

    def as_dict(self) -> dict[str, Any]:
        return self.provenance()
