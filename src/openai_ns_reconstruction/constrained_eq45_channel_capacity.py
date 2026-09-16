"""Pointwise channel-capacity diagnostics for the paper Eq. (4.5) velocity mixing.

This module does not implement or identify the unknown two-dimensional profiles.
For fixed positive q and finite h it studies only the linear map from profile
values ``(v0, F, U)`` to Cartesian velocity ``(u, v, w)``:

    u = x/(2q) v0 - y q^(-1-h) F
    v = y/(2q) v0 + x q^(-1-h) F
    w = q^(-1/2-h) U.

The result is a representation-capacity diagnostic, not a PDE validation and
not evidence that any chosen profiles reproduce the OpenAI field.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Eq45ChannelCapacity:
    radial_distance: np.ndarray
    singular_values: np.ndarray
    numerical_rank: np.ndarray
    raw_condition_number: np.ndarray
    channel_norms: np.ndarray
    max_normalized_cross_dot: np.ndarray


def _broadcast_inputs(x, y, q, h):
    arrays = [np.asarray(value, dtype=float) for value in (x, y, q, h)]
    if not all(np.all(np.isfinite(arr)) for arr in arrays):
        raise ValueError("x, y, q and h must be finite")
    x, y, q, h = np.broadcast_arrays(*arrays)
    if np.any(q <= 0.0):
        raise ValueError("q must be strictly positive")
    return x, y, q, h


def eq45_profile_mixing_matrix(x, y, q, h):
    """Return the pointwise map ``(v0,F,U) -> (u,v,w)``.

    The output has shape ``broadcast(x,y,q,h).shape + (3,3)``. Columns are the
    radial-profile, swirl-profile and axial-profile velocity directions.
    """

    x, y, q, h = _broadcast_inputs(x, y, q, h)
    shape = x.shape
    matrix = np.zeros(shape + (3, 3), dtype=float)
    radial_factor = 1.0 / (2.0 * q)
    swirl_factor = q ** (-1.0 - h)
    axial_factor = q ** (-0.5 - h)

    matrix[..., 0, 0] = x * radial_factor
    matrix[..., 1, 0] = y * radial_factor
    matrix[..., 0, 1] = -y * swirl_factor
    matrix[..., 1, 1] = x * swirl_factor
    matrix[..., 2, 2] = axial_factor
    return matrix


def diagnose_eq45_channel_capacity(x, y, q, h, *, rank_rtol=1e-12):
    """Diagnose pointwise rank/conditioning of the Eq. (4.5) profile channels.

    ``raw_condition_number`` is infinite where the pointwise map is rank
    deficient. Because the first two columns vanish exactly on the symmetry
    axis, rank one there is expected and should not be "fixed" by adding an
    arbitrary Cartesian basis direction.
    """

    if not np.isfinite(rank_rtol) or not 0.0 < rank_rtol < 1.0:
        raise ValueError("rank_rtol must be finite and in (0,1)")

    x, y, q, h = _broadcast_inputs(x, y, q, h)
    matrix = eq45_profile_mixing_matrix(x, y, q, h)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    largest = singular_values[..., 0]
    threshold = rank_rtol * largest
    numerical_rank = np.sum(singular_values > threshold[..., None], axis=-1)
    smallest = singular_values[..., -1]
    raw_condition_number = np.full_like(largest, np.inf)
    full_rank = numerical_rank == 3
    raw_condition_number[full_rank] = largest[full_rank] / smallest[full_rank]

    channel_norms = np.linalg.norm(matrix, axis=-2)
    gram = np.swapaxes(matrix, -1, -2) @ matrix
    denom = channel_norms[..., :, None] * channel_norms[..., None, :]
    normalized = np.zeros_like(gram)
    np.divide(gram, denom, out=normalized, where=denom > 0.0)
    diag = np.arange(3)
    normalized[..., diag, diag] = 0.0
    max_normalized_cross_dot = np.max(np.abs(normalized), axis=(-2, -1))

    return Eq45ChannelCapacity(
        radial_distance=np.hypot(x, y),
        singular_values=singular_values,
        numerical_rank=numerical_rank,
        raw_condition_number=raw_condition_number,
        channel_norms=channel_norms,
        max_normalized_cross_dot=max_normalized_cross_dot,
    )


def analytic_eq45_channel_norms(x, y, q, h):
    """Return exact column norms ``(||v0||, ||F||, ||U||)`` for regression."""

    x, y, q, h = _broadcast_inputs(x, y, q, h)
    r = np.hypot(x, y)
    return np.stack(
        (r / (2.0 * q), r * q ** (-1.0 - h), q ** (-0.5 - h)),
        axis=-1,
    )
