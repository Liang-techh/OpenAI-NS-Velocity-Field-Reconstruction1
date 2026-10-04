# Whole original steep-exit admissible two-vector cone — 2026-10-04

The original steep-exit phase now has a whole-domain cone proof against its same complete waiting/collar/Gamma future moments. The accepted similarity and physical right joins compose this region with the already accepted waiting and collar cones. The unchanged physical field has a receipt-checked cone wrapper. More inward regions, global flat-remainder/energy bounds and genuine coefficient recursion remain unfinished.

## Source domain and varying shear strength

The original phase is `t in [0,1]`, the source axial variable is `Z in [-1,1]`, and the exact radius relation is

```
q = log(R/Rtail) = -wait-1+t
K = (1-epsilon)*exp(k*f(t))
f(t) = integral_t^1 (1-original_sigma(v))dv
delta=2*a, k=1-a, bh=1/2+a
```

The source shape is independent of Z. The actual flat sigmoid satisfies `0<=sigma<=1`, so

```
Kq/K = k*(sigma-1)
kappa-2 = delta+2*k*(1-sigma)
delta <= kappa-2 <= 2
```

The upper 2 follows from the exact source identity `delta+2*k=2`; it is not obtained by adding independently rounded parameter endpoints. The new source bridge guards the actual parameter/rate assignments and stress helper calls, and consumes the accepted full-history/pressure/meridional-zero bridges.

The angular shear is

```
Stheta = Qtheta*2*S*exp(-q)*[Kq-(1+a)*K]
Kq-(1+a)*K = -[1+a+k*(1-sigma)]*K < 0
```

Here exact `S=1/Rtail>0` and exact K are original sources. The directed shear enclosure includes zero because its S enclosure has zero as its lower endpoint. That enclosure is not used to infer strict negativity. `Sz=0` follows from the same actual pure-swirl/zero-meridional source route.

## Whole-domain directed stress bounds

`CompliantSteepExitStressC3.steep_out([-1,1],[0,1])` supplies one directed enclosure covering the complete phase/Z rectangle, with the same full angular, energy and absolute-pressure moments. The cone companion takes the certified lower bound of its zeroth normalized theta row and the certified upper bound of the absolute axial row:

```
Ctheta >= ctheta_lower > 0
abs(Cz) <= cz_upper < infinity
```

For the current source family, the theta lower bound is approximately `2.15988653151524e18`, and the axial absolute upper bound is approximately `1.22501134604699`. These are bounds for normalized source coefficients. They are neither selected physical stresses nor exact extrema. The full physical stress retains its original positive amplitude/radius/nu/lambda factors.

## Actual amplitude and cone inequality

The common source factors satisfy

```
Ttheta=Qtheta*Ctheta
Tz=Qz*Cz
Qz=Qtheta*B
B=Ev0*theta_base*exp(-bh*q)
Bmax=Ev0*thetaT*exp(bh)/(1-epsilon)
```

B decreases across the entire original exit. The original log-source function is replayed symbolically at `q=-wait-1`. Its `waiting_and_current` component is exactly `+bh`; this correlated cancellation is performed before interval enclosure. No independent huge waiting-log boxes are subtracted, and no cap is materialized as a replacement field amplitude.

Because Ctheta is positive and the exact shear is negative, `T dot S<0`. The original two-vector directional cone reduces to

```
2*Ctheta^2 - (kappa-2)*B^2*Cz^2 > 0
```

A sufficient whole-domain inequality, using the exact upper 2, is

```
log(2) + 2*log(Bmax) + 2*log(cz_upper/ctheta_lower) < log(2)
```

The directed source log bound satisfies this strictly. All eleven positive source/stress margins pass. This proof covers every original phase and Z; finite-grid signs and sampled overlaps are not used as proof.

## Physical interface and composed tail

`CertifiedSteepExitPhysical.steep_out(Z,phase,log_tau,theta,viscosity)` validates the accepted checker receipt, all input hashes and the common family/source. It delegates to the unchanged `CompliantSteepExitPhysicalC2` and adds regional cone metadata. It preserves the completed physical tensor and the generally nonzero `Etheta=-nu*partial_zz(utheta)`.

The accepted source stress mixed3, pressure mixed4 and physical stress/remainder right joins meet the waiting region at phase1. Together with the accepted waiting/collar cone, the nonzero stress cone now covers

```
Rtail*exp(-wait-1) <= R < Rtail*exp(3)
```

Gamma exterior stress is exactly zero separately. The strict cone excludes that zero-stress endpoint/exterior. Physical scope remains fixed `nu>0`, `tau>0`, `r>0`, finite `|Z|<1`; source Z endpoints are physical infinity limits. The cone concerns the two original vectors, not the additional diagonal of the completed tensor.

## Reproduction and evidence

With current prerequisite receipts:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepexitcone
```

The stage orders the new cone producer then checker; `all` places them after steep-exit physical recovery. The focused producer/checker passed with 338 current checker input hashes. The checker replays exact cone and source-log identities, recomputes the entire actual phase/Z enclosure and all margins, and consumes accepted similarity/physical joins. A bounded read-only Luna/max review found no material normalization or source-transfer gap. A physical-wrapper smoke call on the actual interior source at `nu=.7` passed and retained the regional nonzero-remainder flag. The existing native Cartesian physical-decomposition fixture is inherited; it is not repeated or promoted into a global NS certificate.

## Next dependencies

Transport the SAME complete angular/energy/pressure moments backward through original Ts-long steep power, then steep entry. Bind phase/radius and amplitude units, retain all meridional primitive source identities and prove both functional joins. Continue the same construction through angular repair, following power and the original flatten. Each region needs its own physical stress/remainder and cone proof.

Upstream finite-width Ra-to-R100 feedback, independent global flat/volume estimates, required-domain physical energy, common leading inputs, the coupled n=1 solve, n>=2 recursion, oscillatory correction and the full corrected Cartesian residual remain open. No global/full-tensor/recursion completion flag is admitted by this increment.
