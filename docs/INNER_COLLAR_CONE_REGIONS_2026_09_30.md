# Conditional inner collar cone regions — 2026-09-30

The paper-defined strict-to-relaxed marker is now computed with directed
logarithmic inversion, rather than inferred from the coarse stress samples.
Run `python experiments/root_st073/lei_ren_part1_paper_collar_cone_regions.py`.
The exact interval receipt includes input provenance and explicit open gates.

Cached Lei–Ren text lines 5702–5747 defines K from global C3 core/frozen
norms, including reciprocal F and normalized inertial stresses. It sets
gamma=.01 and h_b=epsilon_b=c_star K^-100. Lines 5891–5913 and
16770–16790 define

s_an=sigma_inverse(gamma/(10 K^10)),

R_an=R_a exp(h_b s_an).

The numerical source runs currently inject log K=1e152 and log c_star=-100,
equivalently h_b=epsilon_b=exp(-100-100*1e152). These inputs have not been
derived from certified uniform C3 norms or certified theorem constants.
The conditional marker is s_an approximately 3.16227766016837933e-77.
Ordinary samples at s=.25,.5,.75,1,2 therefore miss the whole strict collar.

The receipt preserves log(log(R_an/R_a))=log(h_b)+log(s_an). It deliberately
does not compute R_an-R_a by subtracting almost equal radii. The inverse
bracket is directed and also checked against an independent resolvable
sigma inverse. The switch convention is exp(-1/s^2), not exp(-1/s).

For the retained conditional inlet at Z=.3, F and I_theta are positive.
The ratio H=(D^2+E^2)/D has log approximately 1e152. Thus log(epsilon_b H)
is approximately -9.9e153: the inlet direction is compatible with the
paper's small-kappa relaxed region after the switch. This is an inlet
diagnostic, not a uniform frozen-comparison bound or an actual-field cone
certificate. Pointwise necessary contributions log(1/F), log|E| and
log(1/R_a) are below the assumed log K; C3 sup norms remain unproved.

## Edge-factored margins

The new `factored_margins` helper assumes I/F=q+e, S/F=-chi q,
q=(D,E), and |e|<=rho(1-chi). Let Q=|q| and kappa=chi Q^2/D.
After dividing out positive vanishing-edge factors, it returns

negative dot lower bound: Q^2-rho Q;

strict angular lower bound: 2(Q^2-rho Q)^2-(kappa-2)rho^2 Q^2;

relaxed small-kappa lower bound: Q^2-2D-(1-chi)rho Q.

The angular bound is only supplied when the dot lower bound is positive
and kappa>2. At chi=0 or 1 these are limiting coefficients, not strict
physical inequalities. A rational D=4, E=3, chi=.9, rho=.1 fixture checks
all four returned quantities independently. The helper does not substitute
an invented rho for the missing actual error estimate.

## Next dependencies

1. Bound the global C3 core/frozen quantities entering K, including
   reciprocal F/D and the full axial domain. A single Z=.3 inlet is insufficient.
2. Establish K1 and the c_star restrictions involving epsilon0 and eta_tol.
3. Enclose the normalized comparison error with its essential (1-chi)
   factor; the paper gives rho=K1 K^20(h_b+epsilon_b). An absolute sampled
   stress error cannot replace this edge-factored bound.
4. Apply strict margins only on (R_a,R_an], and relaxed margins on the
   later connection. Retain source/pressure/width remainder contributions.
5. Continue through later angular/axial switches and all other regions,
   then flat remainder and true n-dependent time recursion.

No finite-width cone, global matching, infinite-order flatness, or temporal
recursion completion is claimed by this calculation.
