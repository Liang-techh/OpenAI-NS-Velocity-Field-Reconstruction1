# First-collar continuation of the finite first Volterra word, 2026-10-30

This checkpoint continues `n1_core_word1_z2_v2` from the same core exit into the actual first collar. It transports the saved, signed-sector enclosure for ordinary `Z`-derivative `W1_Z2` without resetting the core history. It also bounds the first-collar contribution separately. The latter is a symmetric analytic error radius, not a computed signed function. This is a finite first-word enclosure on this collar, not a completed `n1` solution, an m15 calculation, or global five-moment closure.

## Definition and no-reset transport

On the selected first collar, with the saved same-source scale `sc` and Z-independent `hb`, set `s = sc*q/4`, `X = 2 exp(hb*s/2)` and `rho = X^2`, for `0 <= q <= 1` and real `Z` in `[-1,1]`. These are the original collar coordinates. `hb` remains a formal same-source factor; it is not replaced by an error-bound value. At fixed phase `s`, the homogeneous part transports the full core-exit `W1_Z2` data `(F,U,K,P,Fx,Ux)` as

```
F(X)  = F(2) + (1-exp(-hb*s))*Fx(2)
U(X)  = U(2) + hb*s*Ux(2)
K(X)  = exp(-hb*s)*K(2) - (1-exp(-hb*s))*Ux(2)
P(X)  = P(2)
Fx(X) = exp(-3*hb*s/2)*Fx(2)
Ux(X) = exp(-hb*s/2)*Ux(2)
```

Every input sector keeps its original source metadata and exact directed coefficient interval. The result has 808, 1056, 1056, 24, 404 and 528 transported sectors for `F1`, `Uz1`, `K1`, `P1`, `d_x_F1` and `d_x_Uz1`, respectively. These are algebraic sectors, not independent samples or new moment terms. Bounds for the exponential multipliers preserve their common `hb` factors; at `q=0`, the transport contains the core inlet exactly. The pre-existing core errors remain in the transported intervals.

## Separate local feedback bound

The saved same-source tube bounds give `H = Mbase*(Cp + C1/loss)` for `h = C0 W0 + B1 W0_Z` on the relevant analytic tube. With the separately verified gap `eta/16`, Cauchy's estimate gives `H2 = 2*H/(eta/16)^2` for `h_ZZ` on the real interval. This bound does not identify the sign or value of `h_ZZ`.

For the local Green integral, `u=2 exp(hb*v/2)` gives `du=hb*u*dv/2`. The kernels include `u*(1-(u/X)^2)/2`, `u*log(X/u)`, their negative or the unit/decay kernels. The first two are bounded after the change of variable by `X-u`, and the latter by 1. Consequently local `F`, `U`, `K` have symmetric radii `H2*(Delta x)^2`, while `P`, `Fx`, `Ux` have `H2*Delta x`, where `Delta x = X-2 <= hb*s/(1-sc/16)`. The exported radii retain `hb^2` or `hb` explicitly. They must be added as errors, not treated as signed source values or used to cancel other errors. The all-word tail remains separate and has not been added here.

## Run and provenance

The uniquely entered runner `scripts/run_n1_word1_z2_first_collar_once.py` returned successfully after raw-hash and same-family/seam/gap checks. Its independent synthetic fixture passed 18 checks: six exact homogeneous propagator comparisons, six inlet identities and six physical `Delta x` bounds. The fixture is not an actual field reconstruction. All directed output sectors, six local radii, exact input bindings and the entry record are preserved in `experiments/root_st073/n1_word1_z2_first_collar_v1_receipt.json.gz` (SHA-256 `89e4266ceaa89aa69c8ac41ac680b3f8fc92c6e6ee13e7b10f616b825b15cce9`); its one-entry result SHA-256 is `fb8944f0be388e1abb1493c1bb5e7e54c29b9f438523fa7e1f6cb5f670c9f393`.

No m15 integral or five-moment terminal identity was computed. `P1` is a field component and is not m15. The local bound can be very wide and is not evidence of endpoint cancellation. No old moment term was added or removed.

## Next executable step

Construct actual same-source, signed first-collar `W0`, `C0` and `B1` Z-jets, and integrate the resulting local `W1_Z2` Green source while preserving the transported core history. Keep its independent analytic tail and prior errors; only then connect the relevant signed contributions to the independently defined m15 comparison. Do not infer a signed local value from the radius or claim terminal closure from interval overlap.
