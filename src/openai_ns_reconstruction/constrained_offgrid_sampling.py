"""Independent off-grid validation sampling for constrained reconstruction.

Truth boundary
--------------
This module only generates deterministic held-out sample locations.  It does
not evaluate a Navier--Stokes residual, set or relax an acceptance threshold,
or turn a sampled check into a PDE proof.  In particular, a green test of this
sampler is not evidence that any candidate satisfies the PDE.

External method provenance
--------------------------
We use SciPy's public ``scipy.stats.qmc.Sobol`` API and ``random_base2`` rather
than copying SciPy implementation code.  Sobol points are used here only as a
space-filling *independent diagnostic* to reduce structured-grid blind spots.
"""

from __future__ import annotations

from dataclasses import dataclass
import inspect

import numpy as np
from scipy.stats import qmc


@dataclass(frozen=True)
class SobolValidationSample:
    """Metadata-bound held-out Sobol sample in a rectangular box."""

    points: np.ndarray
    seed: int
    m: int
    scramble: bool

    @property
    def count(self) -> int:
        return int(self.points.shape[0])

    @property
    def dimension(self) -> int:
        return int(self.points.shape[1])


def _sobol_engine(dimension: int, *, scramble: bool, seed: int) -> qmc.Sobol:
    """Construct Sobol compatibly across the repository's SciPy >=1.10 range.

    New SciPy releases prefer ``rng=`` while older supported releases expose
    ``seed=``.  Signature inspection avoids emitting deprecation warnings in
    warning-as-error test runs while preserving deterministic scrambling.
    """
    parameters = inspect.signature(qmc.Sobol).parameters
    if "rng" in parameters:
        return qmc.Sobol(
            d=dimension,
            scramble=scramble,
            rng=np.random.default_rng(seed),
        )
    return qmc.Sobol(d=dimension, scramble=scramble, seed=seed)


def sobol_validation_sample(
    bounds: np.ndarray,
    *,
    m: int,
    seed: int,
    scramble: bool = True,
) -> SobolValidationSample:
    """Generate ``2**m`` deterministic low-discrepancy points in ``bounds``.

    ``bounds`` has shape ``(dimension, 2)``.  ``random_base2`` is deliberately
    used so the Sobol balance property is not silently discarded by an
    arbitrary sample count.  Callers must keep ``seed`` distinct from
    optimizer/training seeds.

    This function is a sampling primitive only; it does not evaluate or accept
    a PDE candidate.
    """
    box = np.asarray(bounds, dtype=float)
    if box.ndim != 2 or box.shape[1] != 2 or box.shape[0] < 1:
        raise ValueError("bounds must have shape (dimension, 2)")
    if not np.all(np.isfinite(box)):
        raise ValueError("bounds must be finite")
    if np.any(box[:, 1] <= box[:, 0]):
        raise ValueError("each upper bound must exceed its lower bound")
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or m < 0:
        raise ValueError("m must be a nonnegative integer")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")

    engine = _sobol_engine(box.shape[0], scramble=bool(scramble), seed=int(seed))
    unit = engine.random_base2(m=int(m))
    points = qmc.scale(unit, box[:, 0], box[:, 1])

    if points.shape != (2 ** int(m), box.shape[0]) or not np.all(np.isfinite(points)):
        raise RuntimeError("Sobol sampler returned malformed/nonfinite points")
    if np.any(points < box[:, 0]) or np.any(points > box[:, 1]):
        raise RuntimeError("Sobol sampler returned points outside declared bounds")

    points.setflags(write=False)
    return SobolValidationSample(
        points=points,
        seed=int(seed),
        m=int(m),
        scramble=bool(scramble),
    )


def assert_held_out_seed(*, training_seed: int, validation_seed: int) -> None:
    """Fail closed if training and validation reuse the same declared seed."""
    if training_seed == validation_seed:
        raise ValueError("validation_seed must differ from training_seed")
