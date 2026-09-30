# Seeded continuous angular moments through the preheat exterior

`ContinuousIncomingProfile.angular_moments_jet(logR, Z)` now returns the
actual angular primitive, full swirl-square integral, and both Z derivatives
from Rh through Rtail. It reads the actual corrected inner moments and raw
quadratic integrals at Rh. The inner raw `swirl` is half of the square
integral, so this API uses `2*raw['swirl']` and its corresponding derivative.
The actual inner offsets are added exactly once; no terminal target replaces
a measured primitive.

Define E0 = log U_base(R,0). The reference normalization is Z independent:

- X = Mtheta / (sqrt(2) R^(3/2) exp(E0)).
- Y = integral(Utheta^2 dR) / (R exp(2 E0)).

Before the compact angular correction, the forward equations are
X' = G - (3/2 + E0') X and Y' = G^2 - (1 + 2 E0') Y, where
G = U_base(R,Z)/U_base(R,0). Their Z equations differentiate the same G.
Constant stages use exact exponential transfers with `expm1` at small
rates; finite transition stages use declared Gauss nodes and continuous MP
schedule evaluations. Astronomically long constant stages are never sampled
on a uniform grid.

The two compact bumps are integrated separately using the installed owned
beta atoms. Their theta weight is exp((1-mu)t), and their energy weight is
exp(-2mu*t), where t=log(R/Rrel). Nonoverlap removes quadratic cross terms.
The returned `bump_normalized` and `inner_offsets` preserve separate signed
terms even when adding a tiny correction to the baseline loses it at the
working precision.

The executable checks actual Rh seed agreement, radial derivatives of the
four cumulative quantities against actual velocity integrands in flattening,
and the separately retained tiny bump primitives against point corrections.
It writes a reproducible receipt with samples through Rtail.

This is nominal numerical propagation. Transition quadrature and inherited
inner uncertainty are not enclosed. Heat continuation is explicitly rejected;
complete five moments, exact mean cancellation, finite energy, stress and
scale recursion remain open. A small derivative defect is not a proof of any
of those global properties.
