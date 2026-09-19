"""Fail-closed PA.10 source choice-order guard for the Kokuno leading path.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

A fresh source audit exposes a dependency ordering that matters for the
executable leading profile.  The reader first fixes the amplitude-independent
bounds ``B_k`` and ``T_sh`` and only *afterward* enlarges ``C`` and defines
``X_R``.  In particular, with ``F=phi/C`` and ``E=sqrt(2X)F``, the PA.10
boundary datum satisfies the exact cancellation

    ell_i = log(C E_i)
          = log(phi_i) + 1/2 log(2 X_i),

so changing the later normalization constant C cannot be used, on the source
route, to improve the already-fixed ``B_0/T_sh`` datum.  The source also gives

    ||log phi(X_i)||_{C^k_eta}
      <= ||log phi_* + log Phi(4,.)||_{C^k_eta}
         + 1/2 (M_k + 0.8) L_0,

and obtains ``||ell_i||_{C^k_eta} <= B_k`` after adding
``1/2 log(2 X_i)``.  Then

    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)

is fixed before the later choice of C.

This module makes those public dependencies executable and prevents the
repository-autonomous coupled-C experiment from being promoted to a source
PA.10 certificate.  It does not reject that experiment as an autonomous
candidate; it only records that its C-dependent upstream repropagation is not
the source choice order and therefore cannot discharge the source B_0/T_sh
handoff to PA.16.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_coupled_c_appendix_b_normalization import (
    KokunoCoupledCAppendixBNormalization,
)
from .kokuno_rescaled_inner_join_tsh_certificate import LOG_F_SUP, SIGMA_PRIME_SUP


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-choice-order-guard-v1"
X_I = 110.0

_SOURCE_FORMULAS = {
    "normalization": "F=phi/C; E=sqrt(2X)F",
    "boundary_log": (
        "ell_i=log(C E_i)=log(phi_i)+0.5 log(2 X_i); "
        "the later C cancels exactly"
    ),
    "B_k_bound": (
        "||log phi(X_i)||_{C^k_eta} <= "
        "||log phi_*+log Phi(4,.)||_{C^k_eta} + 0.5(M_k+0.8)L_0; "
        "adding 0.5 log(2X_i) gives ||ell_i||_{C^k_eta}<=B_k"
    ),
    "PA10_T_sh": "T_sh>=20||sigma'||_inf(B_0+||log f||_inf)",
    "choice_order": (
        "(M_d,T_d,P_*,lambda,h) -> moment tolerance,j_0 -> "
        "delta_*,sigma_*,Lambda -> (B_k),T_sh -> C,X_R -> "
        "kappa_0,t_1,final widths"
    ),
}

_SOURCE_CHOICE_ORDER = (
    "base_outer_parameters",
    "moment_tolerance_and_j0",
    "delta_sigma_Lambda",
    "B_k_and_T_sh",
    "C_and_X_R",
    "kappa0_t1_and_final_widths",
)

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_choice_order_recorded": True,
    "source_B_k_fixed_before_C": True,
    "source_T_sh_fixed_before_C": True,
    "source_ell_i_C_cancellation_executable": True,
    "displayed_B0_bound_formula_executable": True,
    "autonomous_coupled_C_candidate_may_remain_executable": True,
    "C_oversize_may_discharge_source_B0_or_T_sh": False,
    "coupled_C_pointwise_pass_is_source_T_sh_certificate": False,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return out


@dataclass(frozen=True)
class KokunoPA10SourceChoiceOrderGuard:
    """Executable source-order guard for the PA.10 -> PA.16 handoff."""

    X_i: float = X_I
    sigma_prime_sup: float = SIGMA_PRIME_SUP
    log_f_sup: float = LOG_F_SUP

    def __post_init__(self) -> None:
        X_i = float(self.X_i)
        sigma_prime_sup = float(self.sigma_prime_sup)
        log_f_sup = float(self.log_f_sup)
        if not math.isfinite(X_i) or X_i <= 0.0:
            raise ValueError("X_i must be finite and positive")
        if not math.isfinite(sigma_prime_sup) or sigma_prime_sup <= 0.0:
            raise ValueError("sigma_prime_sup must be finite and positive")
        if not math.isfinite(log_f_sup) or log_f_sup < 0.0:
            raise ValueError("log_f_sup must be finite and nonnegative")
        object.__setattr__(self, "X_i", X_i)
        object.__setattr__(self, "sigma_prime_sup", sigma_prime_sup)
        object.__setattr__(self, "log_f_sup", log_f_sup)

    @property
    def half_log_2X_i(self) -> float:
        return 0.5 * math.log(2.0 * self.X_i)

    @property
    def pa10_prefactor(self) -> float:
        return 20.0 * self.sigma_prime_sup

    def source_log_E_i(self, log_phi_i: float, log_C: float) -> float:
        """Return the normalized source ``log E_i`` for a *later* fixed C."""
        log_phi_i = float(log_phi_i)
        log_C = float(log_C)
        if not math.isfinite(log_phi_i):
            raise ValueError("log_phi_i must be finite")
        if not math.isfinite(log_C) or log_C <= 0.0:
            raise ValueError("log_C must be finite and positive")
        return self.half_log_2X_i + log_phi_i - log_C

    def source_ell_i(self, log_phi_i: float) -> float:
        """Source PA.10 boundary datum after the exact C cancellation."""
        log_phi_i = float(log_phi_i)
        if not math.isfinite(log_phi_i):
            raise ValueError("log_phi_i must be finite")
        return self.half_log_2X_i + log_phi_i

    def verify_C_cancellation(
        self, log_phi_i: float, log_C: float, *, atol: float = 1.0e-12
    ) -> dict[str, float | bool]:
        """Numerically replay the exact algebraic cancellation for diagnostics."""
        if not math.isfinite(float(atol)) or float(atol) < 0.0:
            raise ValueError("atol must be finite and nonnegative")
        log_E_i = self.source_log_E_i(log_phi_i, log_C)
        via_normalized = float(log_C) + log_E_i
        direct = self.source_ell_i(log_phi_i)
        error = abs(via_normalized - direct)
        return {
            "log_E_i": log_E_i,
            "ell_i_via_log_C_plus_log_E_i": via_normalized,
            "ell_i_direct_C_independent": direct,
            "absolute_cancellation_error": error,
            "cancellation_passed": bool(error <= float(atol)),
        }

    def displayed_B0_upper_bound(
        self,
        combined_profile_C0_norm: float,
        M0: float,
        L0: float,
    ) -> float:
        """Evaluate the displayed *formula* for the k=0 bound.

        The three inputs are deliberately caller supplied.  Evaluating this
        formula does not set ``source_B0_dependencies_machine_bound`` or prove
        that the values are the source values.
        """
        profile_norm = _finite_nonnegative(
            combined_profile_C0_norm, "combined_profile_C0_norm"
        )
        M0 = _finite_nonnegative(M0, "M0")
        L0 = _finite_nonnegative(L0, "L0")
        return profile_norm + 0.5 * (M0 + 0.8) * L0 + self.half_log_2X_i

    def displayed_T_sh_lower_bound_from_B0(self, B0: float) -> float:
        B0 = _finite_nonnegative(B0, "B0")
        return self.pa10_prefactor * (B0 + self.log_f_sup)

    def classify_coupled_candidate(
        self, candidate: KokunoCoupledCAppendixBNormalization
    ) -> dict[str, Any]:
        """Classify #521-style autonomous coupled-C evidence fail closed."""
        if not isinstance(candidate, KokunoCoupledCAppendixBNormalization):
            raise TypeError(
                "candidate must be a KokunoCoupledCAppendixBNormalization"
            )
        point = candidate.pointwise_report()
        return {
            "candidate_schema": candidate.to_payload()["schema"],
            "candidate_sha256": candidate.sha256,
            "log_C_factor": float(candidate.log_C_factor),
            "autonomous_coupled_C_pointwise_screen_passed": bool(
                point["not_excluded_by_necessary_screen"]
            ),
            "autonomous_pointwise_strict_margin": float(point["strict_margin"]),
            "autonomous_pointwise_ell_i": float(point["ell_i"]),
            "source_choice_order_allows_C_to_reparameterize_upstream_phi": False,
            "source_choice_order_allows_C_to_discharge_B0_or_T_sh": False,
            "source_B0_dependencies_machine_bound": False,
            "source_B0_analytic_bound_proved": False,
            "source_T_sh_lower_bound_verified": False,
            "selected_pa16_handoff_allowed": False,
            "interpretation": (
                "candidate-only pointwise evidence; later source C is chosen after "
                "B_k/T_sh and cancels from ell_i"
            ),
        }

    def report(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "source_choice_order": list(_SOURCE_CHOICE_ORDER),
            "geometry": {
                "X_i": self.X_i,
                "half_log_2X_i": self.half_log_2X_i,
                "sigma_prime_sup": self.sigma_prime_sup,
                "log_f_sup": self.log_f_sup,
                "pa10_prefactor": self.pa10_prefactor,
            },
            "unresolved_source_B0_dependencies": [
                "||log phi_*+log Phi(4,.)||_{C^0_eta}",
                "M_0",
                "L_0",
            ],
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": self.report()["source"],
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "source_choice_order": list(_SOURCE_CHOICE_ORDER),
            "parameters": {
                "X_i": self.X_i,
                "sigma_prime_sup": self.sigma_prime_sup,
                "log_f_sup": self.log_f_sup,
            },
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10SourceChoiceOrderGuard":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected PA.10 source-choice-order schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("PA.10 source formulas changed")
        if payload.get("source_choice_order") != list(_SOURCE_CHOICE_ORDER):
            raise ValueError("PA.10 source choice order changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("PA.10 source-choice truth-boundary metadata changed")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing PA.10 source-choice-order parameters")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("PA.10 source-choice-order payload SHA-256 mismatch")
        obj = cls(**params)
        if obj.to_payload() != payload:
            raise ValueError("PA.10 source-choice-order replay changed payload")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPA10SourceChoiceOrderGuard":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit the source-ordered PA.10 guard for the coupled-C leading path"
    )
    parser.add_argument("--log-c-factor", type=float, default=32.0)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    guard = KokunoPA10SourceChoiceOrderGuard()
    candidate = KokunoCoupledCAppendixBNormalization(log_C_factor=args.log_c_factor)
    report = guard.report()
    report["guard_sha256"] = guard.sha256
    report["selected_coupled_candidate"] = guard.classify_coupled_candidate(candidate)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
