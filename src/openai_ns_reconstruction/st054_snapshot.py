"""Repository-local continuous-time Python evaluator for the published ST054 snapshot.

This module evaluates the same immutable finite-basis coefficients consumed by the
native MATLAB viewer merged in PR #689.  It is a delivery adapter, not a new fit,
not an interpolation of the sampled ST052 grid, and not PDE validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

ST054_SOURCE_COMMIT = "c77492a48e9c0f13d4d51244987c57c28519409b"
ST054_SCHEMA = "ns_matlab_spectral_v1"
ST054_TESTED_MAT_SHA256 = "59d43fdc40cd4df1675c2191f85810b9351802b2f4b1e0ac25a3bccfd9914c62"
ST054_REFERENCE_ATOL = 1.0e-8  # engineering parity only; not a scientific/PDE gate
_ST054_MODEL_IDS = ("ST054-Q2", "ST054-M3")
_DEFAULT_MAT_PATH = (
    Path(__file__).resolve().parents[2]
    / "visualization"
    / "matlab"
    / "data"
    / "st054_models.mat"
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _bump(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    d = np.ones_like(q)
    inside = q < 1.0
    d[inside] = 1.0 - q[inside]
    y = 1.0 / d
    b = np.zeros_like(q)
    b[inside] = np.exp(1.0 - y[inside])
    return np.column_stack(
        (
            b,
            -b * y**2,
            b * (y**4 - 2.0 * y**3),
            b * (-y**6 + 6.0 * y**5 - 6.0 * y**4),
        )
    )


def _rational(q: np.ndarray, k: int) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    d = np.ones_like(q)
    inside = q < 1.0
    d[inside] = 1.0 - q[inside]
    y = 1.0 / d
    b = np.zeros_like(q)
    b[inside] = np.exp(1.0 - y[inside])
    yc = np.minimum(y, 1000.0)
    c = np.zeros(k + 8, dtype=float)
    c[k] = 1.0
    out = np.zeros((q.size, 4), dtype=float)
    for derivative in range(4):
        value = np.zeros_like(q)
        for power, coefficient in enumerate(c):
            if coefficient != 0.0:
                value += coefficient * yc**power
        out[:, derivative] = b * value
        out[y > 1000.0, derivative] = 0.0
        nxt = np.zeros_like(c)
        for power in range(len(c) - 2):
            nxt[power + 1] += power * c[power]
            nxt[power + 2] -= c[power]
        c = nxt
    return out


def _zchain(b: np.ndarray, z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float).reshape(-1)
    return np.column_stack(
        (
            b[:, 0],
            b[:, 1] * z / 2.0,
            b[:, 2] * z**2 / 4.0 + b[:, 1] / 2.0,
            b[:, 3] * z**3 / 8.0 + 3.0 * b[:, 2] * z / 4.0,
        )
    )


def _legendre_jets(x: np.ndarray, n: int) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    jets = np.zeros((x.size, n + 1, 4), dtype=float)
    jets[:, 0, 0] = 1.0
    if n == 0:
        return jets
    jets[:, 1, 0] = x
    jets[:, 1, 1] = 1.0
    for degree in range(2, n + 1):
        for derivative in range(4):
            value = x * jets[:, degree - 1, derivative]
            if derivative > 0:
                value += derivative * jets[:, degree - 1, derivative - 1]
            jets[:, degree, derivative] = (
                (2 * degree - 1) * value
                - (degree - 1) * jets[:, degree - 2, derivative]
            ) / degree
    return jets


def _basis(v: np.ndarray, kind: str, odd: int, n: int) -> tuple[np.ndarray, ...]:
    """Port of visualization/matlab/ns_basis.m, including derivative orders 0:3."""
    v = np.asarray(v, dtype=float).reshape(-1)
    result: list[np.ndarray] = []
    if kind == "r":
        q = v / 4.0
        base = _bump(q)
        for derivative in range(4):
            base[:, derivative] /= 4.0**derivative
        legendre = _legendre_jets(v / 2.0 - 1.0, 6)
        for derivative in range(4):
            values = np.zeros((v.size, n), dtype=float)
            for column in range(min(n, 7)):
                for split in range(derivative + 1):
                    values[:, column] += (
                        math.comb(derivative, split)
                        * base[:, split]
                        * legendre[:, column, derivative - split]
                        / 2.0 ** (derivative - split)
                    )
            for column in range(7, n):
                rational = _rational(q, 2 * (column - 6))
                values[:, column] = rational[:, derivative] / 4.0**derivative
            result.append(values)
        return tuple(result)

    if kind != "z":
        raise ValueError(f"Unsupported ST054 basis kind: {kind!r}")
    if n != 12:
        raise ValueError("The published ST054 snapshot expects exactly 12 axial columns")
    base = _zchain(_bump(v**2 / 4.0), v)
    legendre = _legendre_jets(v / 2.0, 13)
    for derivative in range(4):
        values = np.zeros((v.size, 12), dtype=float)
        for column in range(7):
            degree = 2 * column + odd
            for split in range(derivative + 1):
                values[:, column] += (
                    math.comb(derivative, split)
                    * base[:, split]
                    * legendre[:, degree, derivative - split]
                    / 2.0 ** (derivative - split)
                )
        for column in range(7, 9):
            rational = _zchain(_rational(v**2 / 4.0, 2 * (column - 6)), v)
            if odd:
                values[:, column] = v / 2.0 * rational[:, derivative]
                if derivative > 0:
                    values[:, column] += derivative / 2.0 * rational[:, derivative - 1]
            else:
                values[:, column] = rational[:, derivative]
        for offset in range(3):
            degree = 2 * offset + (1 - odd)
            for split in range(derivative + 1):
                values[:, offset + 9] += (
                    math.comb(derivative, split)
                    * base[:, split]
                    * legendre[:, degree, derivative - split]
                    / 2.0 ** (derivative - split)
                )
        result.append(values)
    return tuple(result)


def _chebyshev_values(t: float, count: int) -> np.ndarray:
    x = 4.0 * float(t) - 2.0
    values = np.zeros(count, dtype=float)
    values[0] = 1.0
    if count > 1:
        values[1] = x
    for column in range(2, count):
        values[column] = 2.0 * x * values[column - 1] - values[column - 2]
    return values


def _as_model_list(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, dict):
        return [value]
    if isinstance(value, np.ndarray):
        items = value.reshape(-1).tolist()
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        items = [value]
    if not all(isinstance(item, dict) for item in items):
        raise ValueError("Malformed ST054 model table")
    return items


@dataclass(frozen=True)
class ST054Snapshot:
    """Immutable repository snapshot exposing continuous ``velocity(x,y,z,t)``."""

    model_id: str
    F: np.ndarray
    G: np.ndarray
    tmin: float
    tmax: float
    nu: float
    support_radius: float
    support_z: float
    source_commit: str
    mat_path: Path
    mat_sha256: str

    def velocity(self, x: Any, y: Any, z: Any, t: float) -> np.ndarray:
        """Evaluate the published field at broadcastable Cartesian coordinates.

        ``t`` is a scalar physical time in [0.25, 0.75].  The return shape is the
        broadcast coordinate shape with a final component axis of length three.
        """
        time = float(t)
        if not np.isfinite(time) or not (self.tmin <= time <= self.tmax):
            raise ValueError(f"t must lie in [{self.tmin}, {self.tmax}]")
        bx, by, bz = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        if not (np.isfinite(bx).all() and np.isfinite(by).all() and np.isfinite(bz).all()):
            raise ValueError("Coordinates must be finite")
        shape = bx.shape
        xf = bx.reshape(-1)
        yf = by.reshape(-1)
        zf = bz.reshape(-1)
        s = xf * xf + yf * yf

        nr, nz, nt = self.F.shape
        if self.G.shape != self.F.shape:
            raise ValueError("Malformed ST054 coefficient shapes")
        radial = _basis(s, "r", 0, nr)
        axial_f = _basis(zf, "z", 1, nz)
        axial_g = _basis(zf, "z", 0, nz)
        temporal = _chebyshev_values(time, nt)
        coeff_f = np.tensordot(self.F, temporal, axes=([2], [0]))
        coeff_g = np.tensordot(self.G, temporal, axes=([2], [0]))

        def evaluate(coeff: np.ndarray, axial: tuple[np.ndarray, ...], dr: int, dz: int) -> np.ndarray:
            return np.sum((radial[dr] @ coeff) * axial[dz], axis=1)

        f00 = evaluate(coeff_f, axial_f, 0, 0)
        f10 = evaluate(coeff_f, axial_f, 1, 0)
        f01 = evaluate(coeff_f, axial_f, 0, 1)
        swirl = evaluate(coeff_g, axial_g, 0, 0)
        a = -f01
        vertical = 2.0 * f00 + 2.0 * s * f10
        velocity = np.column_stack(
            (xf * a - yf * swirl, yf * a + xf * swirl, vertical)
        )
        return velocity.reshape(shape + (3,))


def available_st054_models() -> tuple[str, ...]:
    """Return the two immutable model identifiers bundled with main."""
    return _ST054_MODEL_IDS


@lru_cache(maxsize=4)
def load_st054(model_id: str = "ST054-Q2", path: str | Path | None = None) -> ST054Snapshot:
    """Load a checksum-bound ST054 snapshot from the tested MAT-v5 payload.

    Custom paths are allowed only when their bytes exactly match the native-MATLAB
    tested payload merged in PR #689.  This intentionally does not bless arbitrary
    re-exports or sampled-grid interpolants as the same continuous field identity.
    """
    if model_id not in _ST054_MODEL_IDS:
        raise ValueError(f"Unknown ST054 model {model_id!r}; choose one of {_ST054_MODEL_IDS}")
    mat_path = Path(path) if path is not None else _DEFAULT_MAT_PATH
    mat_path = mat_path.resolve()
    if not mat_path.is_file():
        raise FileNotFoundError(f"Published ST054 MAT payload not found: {mat_path}")
    digest = _sha256(mat_path)
    if digest != ST054_TESTED_MAT_SHA256:
        raise ValueError("ST054 MAT checksum mismatch; refusing unverified field bytes")

    payload = loadmat(mat_path, simplify_cells=True)
    if str(payload.get("schema")) != ST054_SCHEMA:
        raise ValueError("Unexpected ST054 MAT schema")
    if str(payload.get("source_commit")) != ST054_SOURCE_COMMIT:
        raise ValueError("Unexpected ST054 source commit")
    if int(np.asarray(payload.get("pde_validated")).reshape(-1)[0]) != 0:
        raise ValueError("Published ST054 truth boundary was unexpectedly promoted")

    models = _as_model_list(payload.get("models"))
    selected = next((model for model in models if str(model.get("id")) == model_id), None)
    if selected is None:
        raise ValueError(f"Model {model_id!r} is absent from the published ST054 payload")
    if str(selected.get("source_commit")) != ST054_SOURCE_COMMIT:
        raise ValueError("Model/source provenance mismatch")
    if int(np.asarray(selected.get("pde_validated")).reshape(-1)[0]) != 0:
        raise ValueError("Model PDE truth boundary was unexpectedly promoted")

    F = np.asarray(selected["F"], dtype=float)
    G = np.asarray(selected["G"], dtype=float)
    if F.ndim != 3 or G.shape != F.shape or not (np.isfinite(F).all() and np.isfinite(G).all()):
        raise ValueError("Malformed ST054 velocity coefficients")
    tmin = float(selected["tmin"])
    tmax = float(selected["tmax"])
    nu = float(selected["nu"])
    support_radius = float(selected["support_radius"])
    support_z = float(selected["support_z"])
    if (tmin, tmax, nu, support_radius, support_z) != (0.25, 0.75, 0.01, 2.0, 2.0):
        raise ValueError("ST054 physical metadata drift")

    F.setflags(write=False)
    G.setflags(write=False)
    return ST054Snapshot(
        model_id=model_id,
        F=F,
        G=G,
        tmin=tmin,
        tmax=tmax,
        nu=nu,
        support_radius=support_radius,
        support_z=support_z,
        source_commit=ST054_SOURCE_COMMIT,
        mat_path=mat_path,
        mat_sha256=digest,
    )
