"""Divergence-free poloidal profile coupling for the Eq. (4.5) backbone.

This module is a representation-capacity diagnostic. It does not solve the
implicit q relation and does not identify the unknown paper profiles.

For the corrected similarity coordinates

    q - z^2 q^(2h) = 1 - t,
    X = r^2/(2q),
    eta = z/q^(1/2-h),

write an axisymmetric streamfunction as

    psi(r,z,t) = q^(1/2-h) X Phi(X, eta).

Then the poloidal part of Eq. (4.5),

    u_r = r/(2q) v0(X,eta),
    u_z = q^(-1/2-h) U(X,eta),

is divergence-free by construction when v0 and U are coupled through Phi.
The swirl profile F remains an independent divergence-free axisymmetric channel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Eq45PoloidalJet:
    """Profile values and first derivatives needed by the divergence identity."""

    v0: np.ndarray
    U: np.ndarray
    v0_X: np.ndarray
    U_X: np.ndarray
    U_eta: np.ndarray


def _broadcast_finite(*values):
    arrays = [np.asarray(value, dtype=float) for value in values]
    if not all(np.all(np.isfinite(arr)) for arr in arrays):
        raise ValueError("all profile arguments must be finite")
    return np.broadcast_arrays(*arrays)


def eq45_poloidal_divergence_numerator(
    X,
    eta,
    h,
    v0,
    v0_X,
    U,
    U_X,
    U_eta,
):
    """Return the dimensionless numerator of the Eq. (4.5) poloidal divergence.

    With the corrected positive-exponent similarity relation, direct
    differentiation gives

        q div(u_poloidal) = N / (1 - 2 h eta^2),

    where

        N = (1 - 2 h eta^2)(v0 + X v0_X)
            + (1 - eta^2) U_eta
            - 2 eta X U_X
            - (1 + 2h) eta U.

    Hence ``N == 0`` is the profile-space incompressibility condition. This
    function diagnoses that condition only; it is not a Navier--Stokes residual.
    """

    X, eta, h, v0, v0_X, U, U_X, U_eta = _broadcast_finite(
        X, eta, h, v0, v0_X, U, U_X, U_eta
    )
    if np.any(X < 0.0):
        raise ValueError("X must be nonnegative")
    if np.any((h <= 0.0) | (h >= 0.5)):
        raise ValueError("h must lie in the source coordinate range (0, 1/2)")
    if np.any(np.abs(eta) >= 1.0):
        raise ValueError("this diagnostic uses the similarity interior |eta| < 1")

    S = 1.0 - 2.0 * h * eta * eta
    return (
        S * (v0 + X * v0_X)
        + (1.0 - eta * eta) * U_eta
        - 2.0 * eta * X * U_X
        - (1.0 + 2.0 * h) * eta * U
    )


def eq45_streamfunction_poloidal_jet(
    X,
    eta,
    h,
    phi,
    phi_X,
    phi_eta,
    phi_XX,
    phi_Xeta,
):
    """Map one scalar similarity streamfunction jet to ``(v0,U)`` and derivatives.

    ``Phi`` is defined by ``psi = q^(1/2-h) X Phi(X,eta)``. For finite Phi and
    derivatives, the map has no division by X, so it is regular at the axis
    ``X=0``. The returned profiles satisfy the Eq. (4.5) incompressibility
    identity to floating-point roundoff for any internally consistent smooth Phi.

    The required formulas are

        U = Phi + X Phi_X,

        v0 = [2 eta ((1/2+h) Phi + X Phi_X)
              - (1-eta^2) Phi_eta] / (1 - 2 h eta^2).

    First derivatives used by the profile-space divergence identity are returned
    as well. No q solver or physical-space velocity evaluator is implemented.
    """

    X, eta, h, phi, phi_X, phi_eta, phi_XX, phi_Xeta = _broadcast_finite(
        X, eta, h, phi, phi_X, phi_eta, phi_XX, phi_Xeta
    )
    if np.any(X < 0.0):
        raise ValueError("X must be nonnegative")
    if np.any((h <= 0.0) | (h >= 0.5)):
        raise ValueError("h must lie in the source coordinate range (0, 1/2)")
    if np.any(np.abs(eta) >= 1.0):
        raise ValueError("this diagnostic uses the similarity interior |eta| < 1")

    A = 0.5 + h
    S = 1.0 - 2.0 * h * eta * eta

    U = phi + X * phi_X
    U_X = 2.0 * phi_X + X * phi_XX
    U_eta = phi_eta + X * phi_Xeta

    numerator = (
        2.0 * eta * (A * phi + X * phi_X)
        - (1.0 - eta * eta) * phi_eta
    )
    v0 = numerator / S

    numerator_X = (
        2.0 * eta * ((A + 1.0) * phi_X + X * phi_XX)
        - (1.0 - eta * eta) * phi_Xeta
    )
    v0_X = numerator_X / S

    return Eq45PoloidalJet(v0=v0, U=U, v0_X=v0_X, U_X=U_X, U_eta=U_eta)
