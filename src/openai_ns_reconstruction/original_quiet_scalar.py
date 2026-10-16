"""Stable original quiet own-rate kernels, with separate nonzero-mu terms."""
import math

from .original_core_numeric import _bounds


class OriginalQuietScalar:
    def __init__(self, ctx, iv, mu, mui):
        self.ctx, self.iv, self.mu, self.mui = ctx, iv, mu, mui

    def moments(self, s, n, ctx, *, forward):
        """Positive K_n=int exp(-(s-u))*u^n du or B_n=int exp(-u)*u^n du.

        Short intervals use positive convergent series and a geometric tail.
        Long intervals use exact recurrences/finite incomplete-gamma sums.
        No exponential power series is extended over the long quiet length.
        """
        directed = ctx is self.iv
        high = _bounds(self.ctx, s)[1] if directed else s
        if high <= 4:
            term = ctx.mpf(1) / (n + 1)
            total = term
            for j in range(1, 128):
                term *= (s * (n + j) / (j * (n + j + 1)) if forward
                         else s / (n + j + 1))
                total += term
            next_term = term * (s * (n + 128) / (128 * (n + 129)) if forward
                                else s / (n + 129))
            tail = next_term / (1 - s / 129)
            if directed:
                total += ctx.mpf([0, _bounds(self.ctx, tail)[1]])
            return ctx.exp(-s) * s**(n + 1) * total
        if forward:
            value = -ctx.expm1(-s)
            for j in range(1, n + 1):
                value = s**j - j * value
            return value
        return ctx.mpf(math.factorial(n)) - ctx.exp(-s) * sum(
            (ctx.mpf(math.factorial(n)) * s**j / math.factorial(j) for j in range(n + 1)), ctx.mpf(0))

    def _exprel(self, z, ctx):
        directed = ctx is self.iv
        bound = max(map(abs, _bounds(self.ctx, z))) if directed else abs(z)
        if bound > ctx.mpf('.25'):
            return ctx.expm1(z) / z
        value = sum((z**j / math.factorial(j + 1) for j in range(12)), ctx.mpf(0))
        if directed:
            rem = ctx.exp(bound) * bound**12 / math.factorial(13)
            value += ctx.mpf([-_bounds(self.ctx, rem)[1], _bounds(self.ctx, rem)[1]])
        return value

    def evaluate(self, s, ctx):
        directed = ctx is self.iv
        mu = self.mui if directed else self.mu
        if (_bounds(self.ctx, s)[0] if directed else s) < 0:
            raise ValueError('Nonnegative original quiet log radius required')
        K = {n: self.moments(s, n, ctx, forward=True) for n in range(1, 5)}
        B = {n: self.moments(s, n, ctx, forward=False) for n in range(1, 5)}
        U = -s * self._exprel(-mu * s, ctx)
        # These are exact finite-mu divided differences up to the directed
        # Taylor-Lagrange remainder. The positive moments retain cancellation
        # at s=0 and at arbitrarily small nonzero s.
        Ch = -K[1] + mu * K[2] / 2 - mu**2 * K[3] / 6
        Ce = -s**2 + ctx.mpf(2) / 3 * mu * s**3 - mu**2 * s**4 / 3
        Cp = -2 * B[1] + 2 * mu * B[2] - ctx.mpf(4) / 3 * mu**2 * B[3]
        if directed:
            for name, rem in (('h', mu**3 * K[4] / 24),
                              ('e', ctx.mpf(2) / 15 * mu**3 * s**5),
                              ('p', ctx.mpf(2) / 3 * mu**3 * B[4])):
                positive = ctx.mpf([0, _bounds(self.ctx, rem)[1]])
                if name == 'h': Ch += positive
                elif name == 'e': Ce += positive
                else: Cp += positive
        return dict(decay1=ctx.exp(-s), decay_half=ctx.exp(-s/2), decay_three_half=ctx.exp(-3*s/2),
                    mass1=-ctx.expm1(-s), E_mu_divided=U, h_mu_divided=Ch,
                    e_mass_mu_divided=Ce, p_mass_mu_divided=Cp,
                    **{'K'+str(n): v for n, v in K.items()}, **{'B'+str(n): v for n, v in B.items()})
