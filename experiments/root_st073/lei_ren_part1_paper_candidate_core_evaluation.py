"""Evaluate candidate gauge rows in scaled R without losing the axis datum.

This module only evaluates the supplied finite coefficient enclosures. It does
not identify them with the analytic fixed point or add an infinite tail.
"""

import math


def scaled_mixed_derivative(ctx, rows, *, Lambda, scaled_radius,
                            radial_order=0, axial_order=0, component='Phi'):
    """Enclose a finite derivative at the Taylor center in Z.

    Physical gauge rows are A_n and U_n. With s=Lambda*R,
    Phi=sum A_n*s^n/Lambda^n and
    Psi=Lambda*(U-U0)=sum_{n>=1} U_n*s^n/Lambda^(n-1).
    The exact same stored axis U0 is removed algebraically, avoiding interval
    subtraction of two copies of a non-point datum.
    """
    i, k = radial_order, axial_order
    if not isinstance(i, int) or not isinstance(k, int) or min(i, k) < 0:
        raise ValueError('Nonnegative integer derivative orders required')
    if component not in ('Phi', 'Psi'):
        raise ValueError('Component must be Phi or Psi')
    lam = ctx.mpf(Lambda)
    if lam._mpi_[0][0] or lam._mpi_[0][1] == 0:
        raise ValueError('Strictly positive Lambda required')
    s = ctx.mpf(scaled_radius)
    first = max(i, 1 if component == 'Psi' else 0)
    value = ctx.mpf(0)
    # Horner evaluation of the differentiated polynomial, retaining gaps.
    for n in range(len(rows) - 1, first - 1, -1):
        if len(rows[n]) <= k:
            raise ValueError('Insufficient axial derivative depth')
        factor = math.factorial(k) * math.prod(range(n-i+1, n+1))
        exponent = n if component == 'Phi' else n-1
        coefficient = rows[n][k] * factor / lam**exponent
        value = value * s + coefficient
    return value * s**(first-i)


def self_check():
    """Independent physical/scaled polynomial checks including derivatives."""
    import mpmath as mp
    from mpmath.ctx_iv import MPIntervalContext
    ctx = MPIntervalContext()
    ctx.dps = 70
    with mp.workdps(100):
        lam = ctx.mpf(10)
        # Phi=1+2*s+3*s^2+4*zeta+5*s*zeta+6*s^2*zeta^2.
        a = [[ctx.mpf(v) for v in row] for row in
             ((1,4,0), (20,50,0), (300,0,600))]
        # Psi=2*s+3*s^2+5*s*zeta+6*s^2*zeta^2; U0 is uncertain.
        u = [[ctx.mpf([7,9]),ctx.mpf(0),ctx.mpf(0)],
             [ctx.mpf(2),ctx.mpf(5),ctx.mpf(0)],
             [ctx.mpf(30),ctx.mpf(0),ctx.mpf(60)]]
        references = {
            'Phi': {(0,0):mp.mpf('2.75'),(1,0):mp.mpf(5),
                    (2,0):mp.mpf(6),(0,1):mp.mpf('6.5'),
                    (1,1):mp.mpf(5),(0,2):mp.mpf(3),
                    (1,2):mp.mpf(12),(2,2):mp.mpf(24)},
            'Psi': {(0,0):mp.mpf('1.75'),(1,0):mp.mpf(5),
                    (2,0):mp.mpf(6),(0,1):mp.mpf('2.5'),
                    (1,1):mp.mpf(5),(0,2):mp.mpf(3),
                    (1,2):mp.mpf(12),(2,2):mp.mpf(24)}}
        count = 0
        for name, expected in references.items():
            for (i,k), target in expected.items():
                result = scaled_mixed_derivative(ctx, a if name=='Phi' else u,
                    Lambda=lam, scaled_radius='.5', radial_order=i,
                    axial_order=k, component=name)
                lo, hi = [mp.make_mpf(t) for t in result._mpi_]
                if not lo <= target <= hi:
                    raise AssertionError((name, i, k, result, target))
                count += 1
        print('Independent scaled-coordinate polynomial checks:', count, flush=True)


if __name__ == '__main__':
    self_check()
