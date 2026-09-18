"""Scale-aware selected realization of the Appendix-B stress-activation path.

The corrected public reconstruction leaves ``kappa_0`` as a sufficiently small
existence choice.  Reusing the earlier finite-scale diagnostic value ``0.01``
after the Appendix-A pressure datum is raised to source scale is not numerically
or structurally neutral: at the selected ``Lambda~2.5e27`` it drives ``log F``
out of the binary64 materialization range long before ``X_i=110``.

For the executable selected-pressure path we therefore make one explicit,
source-compatible autonomous choice

    kappa_0 = kappa0_lambda_multiplier / Lambda,

with default multiplier one.  This is a conditioning/existence realization,
not a recovered Kokuno/OpenAI hidden parameter.  The returned object is the
fully serializable :class:`KokunoSourceRescaledAppendixBBoundary`; its payload
records both the numerical ``kappa_0`` and the complete reference payload, so
the relation can be independently replayed.
"""

from __future__ import annotations

import math
from typing import Any

from .kokuno_rescaled_appendix_b_boundary import (
    KokunoSourceRescaledAppendixBBoundary,
)
from .kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)


def make_source_scale_aware_appendix_b_boundary(
    *,
    reference: KokunoSourceRescaledReferenceContinuation | None = None,
    kappa0_lambda_multiplier: float = 1.0,
    **kwargs: Any,
) -> KokunoSourceRescaledAppendixBBoundary:
    """Return the selected-pressure Appendix-B path with ``kappa_0=m/Lambda``.

    ``kappa0_lambda_multiplier`` is deliberately user-visible and guarded.  It
    is not sourced from a hidden reconstruction coefficient.  Callers may pass
    the ordinary numerical controls accepted by
    :class:`KokunoSourceRescaledAppendixBBoundary` via ``kwargs``; supplying a
    second raw ``kappa0`` is rejected so the scale-aware contract cannot be
    silently bypassed.
    """

    if "kappa0" in kwargs:
        raise ValueError(
            "raw kappa0 is incompatible with the scale-aware selected factory"
        )
    if reference is None:
        reference = KokunoSourceRescaledReferenceContinuation()
    if not isinstance(reference, KokunoSourceRescaledReferenceContinuation):
        raise TypeError(
            "reference must be a KokunoSourceRescaledReferenceContinuation"
        )
    multiplier = float(kappa0_lambda_multiplier)
    if not math.isfinite(multiplier) or not 0.0 < multiplier <= 1.0e6:
        raise ValueError("kappa0_lambda_multiplier must lie in (0,1e6]")
    kappa0 = multiplier / float(reference.rescaling_lambda)
    if not math.isfinite(kappa0) or not 0.0 < kappa0 < 0.5:
        raise ValueError("derived kappa0 must satisfy 0<kappa0<1/2")
    return KokunoSourceRescaledAppendixBBoundary(
        reference=reference,
        kappa0=kappa0,
        **kwargs,
    )


def selected_kappa0_relation(
    candidate: KokunoSourceRescaledAppendixBBoundary,
) -> dict[str, float | bool]:
    """Expose the replayable scale relation for routing/audit metadata."""

    if not isinstance(candidate, KokunoSourceRescaledAppendixBBoundary):
        raise TypeError("candidate must be a KokunoSourceRescaledAppendixBBoundary")
    Lambda = float(candidate.reference.rescaling_lambda)
    multiplier = float(candidate.kappa0 * Lambda)
    return {
        "rescaling_lambda": Lambda,
        "kappa0": float(candidate.kappa0),
        "kappa0_lambda_multiplier": multiplier,
        "source_hidden_numeric_choice_recovered": False,
    }
