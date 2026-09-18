"""Corrected Kokuno heat-replacement discrepancy after independent tail audit.

Agent-1 #280 made the public heat-replacement moment discrepancy executable.
Agent-4 #283 then evaluated the same source formula through an independent
small-z series plus adaptive log-space quadrature and isolated one local
implementation error: the first-order ``y>=3`` contribution to ``Delta C_p``
was missing one factor of two.

This module keeps the original v1 object immutable as historical evidence and
provides a fail-closed v2 revision for all future repair/global-assembly work.
No outer schedule scale is inferred here: ``X_star`` and ``e_star`` remain
caller supplied, and no PDE/velocity candidate is promoted by this correction.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy.special import roots_laguerre

from .kokuno_heat_replacement_discrepancy import (
    KokunoHeatReplacementDiscrepancy,
)


SCHEMA_V2 = "kokuno-heat-replacement-discrepancy-v2"
INDEPENDENT_AUDIT_PR = 283
INDEPENDENT_AUDIT_HEAD = "ecf1be298e880afe62b1f568a8f42bd27bc54d4a"
PARENT_IMPLEMENTATION_HEAD = "51076990abdd36bcb7dc598a7abd85d677dee53d"


@dataclass(frozen=True)
class KokunoHeatReplacementDiscrepancyV2(KokunoHeatReplacementDiscrepancy):
    """Heat discrepancy with the audited first-order ``Delta C_p`` tail fixed.

    For ``y>=3`` the source terminal factors satisfy ``chi_K=f_o=1``.  Writing
    ``H(z)-1=z R(z)``, ``z=2d/X`` and ``A=1/2+h``, the first-order contribution
    to ``Delta C_p`` is

    ``2 c_inf^2 d integral X^(-(2A+1)) R(2d/X) d(log X)``.

    The v1 implementation included only one copy of this term.  V2 adds the
    missing copy while retaining the already-correct quadratic tail, ``Delta S``
    and ``Delta I_sub`` channels.
    """

    def _missing_first_order_cp_tail_scalar(self, eta: float) -> float:
        eta_value = float(eta)
        if not math.isfinite(eta_value) or abs(eta_value) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        d = 1.0 - eta_value * eta_value
        if d == 0.0:
            return 0.0

        # This is exactly the one copy absent from v1.  The final v2
        # first-order coefficient is therefore 2*c_inf^2*d, matching the
        # source expansion and Agent-4's independent log-space audit.
        s, weights = roots_laguerre(self.tail_quadrature_order)
        rate = 2.0 * self.A + 1.0
        X0 = self.X_b
        z0 = 2.0 * d / X0
        ratio = self._heat_minus_one_over_z(z0 * np.exp(-s / rate))
        return float(
            self.c_inf**2
            * d
            * X0 ** (-rate)
            / rate
            * np.sum(weights * ratio)
        )

    def missing_first_order_cp_tail(self, eta: Any) -> np.ndarray:
        """Return the v1->v2 ``Delta C_p`` correction on a vectorized eta input."""

        values = np.asarray(eta, dtype=float)
        if not np.all(np.isfinite(values)) or np.any(np.abs(values) > 1.0):
            raise ValueError("eta must contain finite values in [-1,1]")
        flat = values.reshape(-1)
        corrected = np.asarray(
            [self._missing_first_order_cp_tail_scalar(float(value)) for value in flat],
            dtype=float,
        )
        return corrected.reshape(values.shape)

    def _scalar_discrepancy(self, eta: float) -> np.ndarray:
        result = np.asarray(super()._scalar_discrepancy(eta), dtype=float).copy()
        result[0] += self._missing_first_order_cp_tail_scalar(eta)
        return result

    def _unsigned_payload(self) -> dict[str, Any]:
        payload = super()._unsigned_payload()
        payload["schema"] = SCHEMA_V2
        payload["source_formulas"] = dict(payload["source_formulas"])
        payload["source_formulas"]["terminal_delta_cp_first_order"] = (
            "for y>=3 and H-1=zR, Delta C_p first-order tail = "
            "2*c_inf^2*d*int X^(-(2A+1))*R(2d/X)d(log X)"
        )
        payload["implementation_revision"] = {
            "parent_schema": "kokuno-heat-replacement-discrepancy-v1",
            "parent_implementation_head": PARENT_IMPLEMENTATION_HEAD,
            "independent_audit_pr": INDEPENDENT_AUDIT_PR,
            "independent_audit_head": INDEPENDENT_AUDIT_HEAD,
            "change": (
                "restore the missing second copy of the first-order y>=3 "
                "Delta C_p tail term"
            ),
            "delta_s_changed": False,
            "delta_i_sub_changed": False,
        }
        return payload

    @classmethod
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> "KokunoHeatReplacementDiscrepancyV2":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA_V2:
            raise ValueError("unexpected corrected heat-replacement-discrepancy schema")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("parameters are missing")
        obj = cls(**params)
        if obj.to_payload() != payload:
            raise ValueError(
                "corrected heat replacement discrepancy payload hash or content mismatch"
            )
        return obj
