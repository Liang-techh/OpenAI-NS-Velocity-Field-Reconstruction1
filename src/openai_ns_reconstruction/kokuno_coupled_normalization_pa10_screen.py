"""Necessary PA.10 feasibility band for a coupled upstream normalization.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected 2026-09-09 reader
(Zenodo 22678406).

The corrected reconstruction uses one fixed ``C`` in both

    ell_i = log(C E_i)
    X_R = 110 (C P_*)^10,

and PA.10 requires

    T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf),
    ||ell_i||_inf <= B_0,
    X_sep / X_R < exp(-8).

For any one eta witness, these imply the *necessary* strict inequality

    20 ||sigma'||_inf (|ell_i| + ||log f||_inf)
        < 10 (log C + log P_*) - 8.

This module turns that inequality into an executable upstream target.  It
computes the open admissible band for ``ell_i`` and therefore the open band of
``log E_i`` that a newly propagated Appendix-B normalization must reach before
PA.16 can even be considered.  Passing this screen is deliberately weaker
than the source theorem: the unknown analytic ``B_0`` may be larger than a
point witness, so a passing point is only *not excluded*.  Failure, however,
safely excludes that realization.

No post-hoc C multiplier is applied.  A future realization must propagate its
chosen C through the upstream Appendix-B path and then rerun this screen.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_rescaled_inner_join_tsh_certificate import (
    LOG_F_SUP,
    SIGMA_PRIME_SUP,
    KokunoRescaledInnerJoinTshCertificate,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-coupled-normalization-pa10-necessary-screen-v1"
LOG_X_RESTORE_START = -8.0

_SOURCE_FORMULAS = {
    "boundary_log": "ell_i=log(C E_i)=log C + log E_i",
    "B0_bound": "||ell_i||_{C^0_eta} <= B_0",
    "PA10_T_sh": "T_sh >= 20 ||sigma'||_inf (B_0 + ||log f||_inf)",
    "shared_C_outer_scale": "X_R=110(C P_*)^10",
    "separation": "X_sep=110 exp(T_sh); X_sep/X_R<exp(-8)",
    "pointwise_necessary": (
        "20||sigma'||_inf(|ell_i|+||log f||_inf) "
        "< 10(log C+log P_*)-8"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "coupled_shared_C_dependency_recorded": True,
    "pointwise_PA10_necessary_screen_executable": True,
    "required_log_E_i_open_band_executable": True,
    "posthoc_C_multiplier_applied": False,
    "passing_pointwise_screen_is_source_T_sh_certificate": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
    "source_hidden_C_recovered": False,
    "source_hidden_numeric_choices_recovered": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _json_value(value: np.ndarray) -> float | list[float]:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return float(arr)
    return arr.tolist()


@dataclass(frozen=True)
class KokunoCoupledNormalizationPA10Screen:
    """Necessary pointwise PA.10 screen for one coupled ``(C,P_*)`` scale."""

    log_C: float
    log_P_star: float
    sigma_prime_sup: float = SIGMA_PRIME_SUP
    log_f_sup: float = LOG_F_SUP
    log_x_restore_start: float = LOG_X_RESTORE_START

    def __post_init__(self) -> None:
        values = (
            self.log_C,
            self.log_P_star,
            self.sigma_prime_sup,
            self.log_f_sup,
            self.log_x_restore_start,
        )
        if not all(math.isfinite(float(v)) for v in values):
            raise ValueError("PA.10 coupled-normalization parameters must be finite")
        if self.log_C <= 0.0:
            raise ValueError("log_C must be positive")
        if self.sigma_prime_sup <= 0.0:
            raise ValueError("sigma_prime_sup must be positive")
        if self.log_f_sup < 0.0:
            raise ValueError("log_f_sup must be nonnegative")
        if self.log_x_restore_start >= 0.0:
            raise ValueError("log_x_restore_start must be negative")

    @classmethod
    def from_selected_certificate(
        cls, certificate: KokunoRescaledInnerJoinTshCertificate
    ) -> "KokunoCoupledNormalizationPA10Screen":
        if not isinstance(certificate, KokunoRescaledInnerJoinTshCertificate):
            raise TypeError("certificate must be a KokunoRescaledInnerJoinTshCertificate")
        binding = certificate.incoming.binding
        return cls(
            log_C=float(binding.selected_log_C),
            log_P_star=float(binding.template_outer_schedule.log_P_star),
            sigma_prime_sup=SIGMA_PRIME_SUP,
            log_f_sup=LOG_F_SUP,
            log_x_restore_start=LOG_X_RESTORE_START,
        )

    @property
    def pa10_prefactor(self) -> float:
        return 20.0 * self.sigma_prime_sup

    @property
    def t_sh_strict_upper(self) -> float:
        """Strict upper bound implied by ``X_sep/X_R < exp(-8)``."""
        return 10.0 * (self.log_C + self.log_P_star) + self.log_x_restore_start

    @property
    def ell_abs_strict_cap(self) -> float:
        """Largest pointwise ``|ell_i|`` compatible with the necessary screen."""
        return self.t_sh_strict_upper / self.pa10_prefactor - self.log_f_sup

    @property
    def pointwise_feasible_band_nonempty(self) -> bool:
        return bool(self.ell_abs_strict_cap > 0.0)

    def required_log_E_i_open_interval(self) -> tuple[float, float]:
        """Open interval required of ``log E_i`` at every screened witness.

        This interval is only a necessary pointwise target.  It is not a proof
        of the source ``C^0`` bound on ``ell_i``.
        """
        cap = self.ell_abs_strict_cap
        if not cap > 0.0:
            raise ValueError("no pointwise PA.10-feasible ell_i band exists at this scale")
        center = -self.log_C
        return center - cap, center + cap

    def evaluate_log_E_i(self, log_E_i: Any) -> dict[str, Any]:
        values = np.asarray(log_E_i, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("log_E_i must be finite")
        ell = self.log_C + values
        lower = self.pa10_prefactor * (np.abs(ell) + self.log_f_sup)
        margin = self.t_sh_strict_upper - lower
        not_excluded = np.logical_and(self.ell_abs_strict_cap > 0.0, margin > 0.0)
        return {
            "log_E_i": values,
            "ell_i": ell,
            "pointwise_T_sh_lower_bound": lower,
            "t_sh_strict_upper": np.full(values.shape, self.t_sh_strict_upper),
            "strict_margin": margin,
            "not_excluded_by_necessary_screen": not_excluded,
        }

    def evaluation_report(self, log_E_i: Any) -> dict[str, Any]:
        evaluation = self.evaluate_log_E_i(log_E_i)
        return {
            key: (
                bool(np.asarray(value))
                if key == "not_excluded_by_necessary_screen" and np.asarray(value).ndim == 0
                else (
                    np.asarray(value, dtype=bool).tolist()
                    if key == "not_excluded_by_necessary_screen"
                    else _json_value(np.asarray(value))
                )
            )
            for key, value in evaluation.items()
        }

    def geometry_report(self) -> dict[str, Any]:
        interval = None
        if self.pointwise_feasible_band_nonempty:
            interval = list(self.required_log_E_i_open_interval())
        return {
            "log_C": self.log_C,
            "log_P_star": self.log_P_star,
            "sigma_prime_sup": self.sigma_prime_sup,
            "log_f_sup": self.log_f_sup,
            "log_x_restore_start": self.log_x_restore_start,
            "pa10_prefactor": self.pa10_prefactor,
            "t_sh_strict_upper": self.t_sh_strict_upper,
            "ell_abs_strict_cap": self.ell_abs_strict_cap,
            "pointwise_feasible_band_nonempty": self.pointwise_feasible_band_nonempty,
            "required_log_E_i_open_interval": interval,
            "interval_center": -self.log_C,
            "source_T_sh_lower_bound_verified": False,
        }

    def selected_realization_report(
        self, certificate: KokunoRescaledInnerJoinTshCertificate
    ) -> dict[str, Any]:
        selected = self.from_selected_certificate(certificate)
        if selected.to_payload() != self.to_payload():
            raise ValueError("screen scale does not match the selected certificate")
        ell_i = float(certificate.pointwise_report()["ell_i_at_witness"])
        log_E_i = ell_i - self.log_C
        evaluation = self.evaluation_report(log_E_i)
        cap = self.ell_abs_strict_cap
        ratio = math.inf if cap <= 0.0 else abs(ell_i) / cap
        interval = self.required_log_E_i_open_interval() if cap > 0.0 else (math.nan, math.nan)
        return {
            "witness_eta": float(certificate.witness_eta),
            "selected_ell_i": ell_i,
            "selected_log_E_i": log_E_i,
            "selected_abs_ell_over_strict_cap": ratio,
            "required_log_E_i_open_interval": list(interval),
            "evaluation": evaluation,
            "selected_realization_excluded": not bool(
                evaluation["not_excluded_by_necessary_screen"]
            ),
            "source_T_sh_lower_bound_verified": False,
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "parameters": {
                "log_C": self.log_C,
                "log_P_star": self.log_P_star,
                "sigma_prime_sup": self.sigma_prime_sup,
                "log_f_sup": self.log_f_sup,
                "log_x_restore_start": self.log_x_restore_start,
            },
            "geometry": self.geometry_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoCoupledNormalizationPA10Screen":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected coupled-normalization PA.10 screen schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("coupled-normalization PA.10 source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("coupled-normalization PA.10 truth-boundary metadata changed")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing coupled-normalization PA.10 parameters")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("coupled-normalization PA.10 payload SHA-256 mismatch")
        obj = cls(**params)
        if obj.to_payload() != payload:
            raise ValueError("coupled-normalization PA.10 replay changed payload")
        return obj

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoCoupledNormalizationPA10Screen":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit the coupled-normalization PA.10 necessary feasibility band"
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    certificate = KokunoRescaledInnerJoinTshCertificate()
    screen = KokunoCoupledNormalizationPA10Screen.from_selected_certificate(certificate)
    report = {
        "schema": SCHEMA,
        "screen_sha256": screen.sha256,
        "geometry": screen.geometry_report(),
        "selected_realization": screen.selected_realization_report(certificate),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
