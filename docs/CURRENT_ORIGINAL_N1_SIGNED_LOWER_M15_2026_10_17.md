# Actual signed lower m15 and the correlated q1 pressure seam

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The actual lower m15 increment on q in [1,2] is now enclosed with its ordinary Z derivative at the saved axial interval near .371. Its signed value is negative. The actual core/cutoff/lower prefix is joined by a direct pressure primitive over q in [0,2], so the shared q1 pressure endpoint cancels as a function before boxing. This does not identify the old finite numerical rows with exact functions or complete forward Rp/all-Z five-moment matching.

## Existing domain and actual density

Write a=sc/4 and rho(s)=4*exp(hb*s). The lower interval s in [a,2a] is q in [1,2], inside the already accepted native first window [0,sc/2]. This local task needs no later-window identity. The exact selected phase, common axial source, positive hb<1/2 margin, rho(2a)<4/(1-1/2048)<4.1, and Z-independent endpoints are checked against the pinned native/lower source packets.

The own n1 cutoff is eta=1-sigma(4q-2). Its argument is in [2,6] on this lower interval, so eta and all radial derivatives vanish by the exact constant-one branch of sigma. The leading N1p is not cut. The own mean and pressure retain the common core inlet, giving the physical density

```
DeltaM15_lower = -epsilon^2 * integral_rho(a)^rho(2a) rho*N1p_leading drho
               = -epsilon^2*hb * integral_a^2a rho(s)^2*N1p_leading ds.
N1p_leading = -S/2 + 2*b_rho + rho*b_rhorho
Q = Qaxis + epsilon*b, b=Qpsi
S = Saxis + epsilon*S1 + epsilon^2*S2.
```

These are physical N1p units. A normalized N1p/Pstar^2 would require exactly one Pstar^2 restoration. The core radial epsilon is exp(-logLambda), distinct from the pressure-datum epsilon. Current lower function handles are [36592,36593]; the actual C_inner handles are [53161,53162].

## Same-inlet paired primitive

With pA=Psi_core,rho(4,Z), L=1-delta*Z^2 and the original axis Saxis, define

```
K = ((3+delta)*Z*pA-(1-Z^2)*pA_Z)/(2*L)
V = 2*Z*pA/L
Cw_delta = (rho(2a)^2-rho(a)^2)/(2*hb)
         = 8*exp(2*hb*a)*(exp(2*hb*a)-1)/hb.
```

The signed leading approximation is

```
-2*hb*epsilon^2*Cw_delta*K
+hb*epsilon^2*Cw_delta*Saxis/2
+epsilon^2*(1-hb)*V*(rho(2a)^2*sigma(2a)-rho(a)^2*sigma(a)).
```

The exact geometry remains a source function. The directed cover 16a<=Cw_delta<=16a/(1-a)^2 follows from exp(x),exprel(x)<=1/(1-x) and 0<2hb*a<a. It is not selected as the geometry's value. The same actual pA and ordinary derivatives are used at both ends; independently boxed Psi/mean/Qpsi endpoint values are not subtracted. The six-factor basis is hb,F0_anchor,Pstar_squared,Lambda,sigma(a),sigma(2a). Both nonzero prescribed-shear flat factors remain separate from the own velocity cutoff sigma(4q-2).

The standalone lower error includes the sum of the two complete correlated endpoint-model radii E(a)+E(2a), then the lower S1/S2 inertial volume radius with weight hb*Cw_delta/2 once. Whole-error ordinary Z1 uses the same eta/16 Cauchy gap. Uniform source bounds cover real phase<=1/2 and the complete common Z tube, so the q1 endpoint row is not extrapolated as a value at q2.

## Actual shared-endpoint cancellation in the joined prefix

The new joined source is evaluated as

```
core[0,4] + cut_axial[0,a] + cut_swirl[0,a] + leading_pressure[0,2a]
  = C_inner + lower.
```

The release gate independently pins the accepted velocity/core assembly source ASTs and verifies both ordinary orders: complete_first_collar equals the identical q1 pressure signed/error rows plus the two cut-velocity rows, and C_inner equals core plus that collar. The canonical lower source must retain its direct differential definition and identical Z-independent geometry. Both flat factors must have zero ordinary Z derivative.

The two accepted cut-velocity integrals and all their raw/local/all-word errors are used once. The core and its errors are retained once. The pressure part is freshly evaluated directly from 0 to 2a, using only E(2a) and the S1/S2 volume over that full interval. The shared q1 pressure endpoint, including its same-source model error, telescopes before boxing. No interval error is deleted based on overlap, a midpoint or a numerical cap. The old standalone endpoint errors and ledger caps remain untouched in inherited packets.

## Numerical result

The following are rounded displays of directed enclosures. Exact endpoints are stored as MPF tuples.

| ordinary Z order | lower / (hb*epsilon^2*Pstar^2*a) | joined prefix / (epsilon^2*Pstar^2) | lower relative width upper |
|---|---:|---:|---:|
| 0 | approximately -26.755789772317654675 | approximately -13.377894886158827338 | 7.788e-29 |
| 1 | approximately -537.97662348022233944 | approximately -268.98831174011116972 | 2.093e-29 |

The joined prefix retains the same displayed core-dominated digits with widths below 9.531e-29 and 2.548e-29. This does not discard the lower increment: it remains explicit in its own physical source factors.

In the common reference*A*Rv^2 unit, the negative joined prefix has a one-sided cover [-2^-4096,0]. Replacing only the current inner+lower source inside the fresh reserved-target expression leaves the original preRp enclosure and independent negative complete-pulse target unchanged. The new correlated difference bounds are below 1.523e-434 for Z0 and 1.016e-434 for Z1. The pulse still cancels as a function before boxing and is consumed once. These are bounds for this source comparison, not a zero identity or the full forward/Cartesian NS residual.

Production exits 0 in 2.849 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_signed_lower_m15.json` is 13,022,554 bytes, SHA256 `08a5c9cde3fbc33272de00f97deb88bd92832f7df3b55f51001c1632fc6e3e06`, with 1990 input and 336 implementation pins. Focused compilation, six new symbolic identities, directed scalar source arithmetic, saved-only source recomposition and scoped raw index checks pass. Fixed Luna/max workers verified the domain/source ports and independently reviewed the lower factors, paired errors and direct joined primitive. Runtime worker model metadata is not exposed by the tools and is not claimed.

## Preservation and next source interface

No accepted owner/fixture/recurrence/norm/old quadrature was replayed. Every inherited packet is preserved at JSON-value level except stage/hash/contract/check/time metadata. The original six signed rows, Rc178/quiet179/raw213/reduced183, slots 19/20/41/42, all finite-N errors and target/quiet one-use contracts remain active. The old finite prefix is not promoted to the exact raw function. Actual histories remain separate from the direction/comparison histories; equality between those different histories is not a missing condition to impose.

Next bind this actual q2 native prefix to the next first_micro source inlet and its accumulated preRp history, keeping the old six-row approximation/error certificate separate. The already accepted native D/drive and raw/cut identities must not be repeated because a historical inherited packet still has false flags. Certify any remaining finite-approximation difference against the true raw function, without trying to prove a finite polynomial equals its infinite function. Continue actual signed chart integrals and common-N axial coverage before full forward Rp/terminal closure.

Full forward Rp matching, all-five/common-N functional terminal closure, exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2 recursion/summation, mean/oscillatory corrections and independent corrected Cartesian NS remain incomplete. Goal ACTIVE.
