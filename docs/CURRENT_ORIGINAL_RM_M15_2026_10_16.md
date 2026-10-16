# Actual original Rm-to-Rh cumulative pressure moment - 2026-10-16

`OriginalRmM15` continues the complete156-term pressure state atRm through the actual active repair patch toRh. It appends12 signed/error sectors, producing168 cumulative terms. The implementation uses genuine third-source/implicit-control jets, actual partial bump integrals, exact F endpoints, positive native G kernels and both signed/finite-N RQy endpoints. The pressure moment is now integrated across this patch, rather than merely available through a field interface.

The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Selected-N signed correction centers/useful relative cancellation are unresolved. Outer/compact cumulative pressure toRp, its independent reserved-moment value/Z1 seam, global all-five closure/common N, exact heat/finite energy/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open.

## Same-branch mass constraint and four radial shapes

The original implicit mass row has L00=L01=1 and its nonlinear part is identically zero. The admitted unique full-weight branch therefore gives d0+h0+h1=0 as an axial functional identity. This identity comes from the same receipt-backed implicit branch, not from an interval residual merely containing zero.

Factor the original fixed axial normalization t=patch.units[0] from the two physical axial controls, writing h_i=t H_i. This is an unchanged, Z-independent source normalization, not epsilon and not a temporal coefficient. Let I_i(x) be the normalized bump cumulative integral and T_i=1-I_i, for the bumps at5/4 and7/4. Then

```
U = 4Z + t*(H0*beta0 + H1*beta2),
M = 4Z - t*(H0*T0 + H1*T2)/x.
```

The pressure recovery is linear in U,M,M_Z. Consequently its perturbation uses four radial shapes:

```
beta0, beta2, T0/x, T2/x.
```

The exact original half-S1 F/G densities are polynomials of degree0,1,2 in this fixed t. Base3 scalar codes1,3,9,27 distinguish the four shapes and their quadratic products without collisions. All control/source enclosures remain in the coefficients; t powers and physical radius atoms stay separate. The normalized mass diagnostic uses `OriginalRmControlZ3.defects_truth`, which is already normalized. Physical source defects remain under `source.truth.defects`.

At each endpoint the actual normalized pressure is recovered through Q_ZZ from M_ZZZ. Qy follows the original mean ODE My=U-M and the actual x*beta_x derivative. All perturbations vanish identically past71/40, so the upper leading RQy endpoint is exactly zero there; inherited finite-N histories do not vanish.

## Positive native scalar kernels

The physical measure is R*dR=Rm^2*x*dx. `OriginalRmMomentKernels` retains the exact original bump normalization/radius/support. For each bump it integrates positive A=integral x*beta dx, B=integral x*beta^2 dx, J=integral T dx, L=integral T dx/x and S=integral T^2 dx/x.

Before support, J=x-1 and L=S=log(x). On a completed support, normalization and symmetry give A=center, J=center-1 and B=center*integral beta^2 dx. L/S use direct positive cells on the support and then saturate. Their tail enclosures use the monotonic original partial-CDF bounds at both cell endpoints. The reported run uses32 cells per support and128 cells per partial CDF; midpoint quadrature values are proxies with full directed source/quadrature errors.

The other shape products have exact prefix reductions: integral beta_i*T_i dx=(1-T_i(x)^2)/2; beta0*T2 integrates to I0; beta2*T0 and beta0*beta2 vanish by disjoint support. Since T2=1 wherever T0 is nonzero, the tail cross-product reduces to L0. These are prefixes from1; a later implementation of arbitrary subintervals must keep both lower and upper primitives.

The constant U=M=4Z background is integrated as (2F0+G0)*(x^2-1)/2, avoiding subtraction of its nearly equal endpoint values at an empty interval. The complete patch increment remains

```
delta m15 = [R^2 F] + integral R^2 G dy - [R Qy],
dy=dx/x.
```

The existing F/G helpers already include the half-S1 normalization. There is no additional1/2 on G. Only the two nonconstant t degrees need signed F endpoints; their G integrals, the two signed linear-t RQy endpoints and three finite-N error sectors give12 new terms in total.

## Finite-N errors over the whole volume prefix

Both actual lower and upper R^2F/RQy error budgets are retained. For the volume, the final decayed history alone is insufficient: an earlier history can be larger. The new prefix envelope keeps the incoming correction magnitude and every visited source's positive Duhamel mass without later decay. Omitting that decay is conservative and provides a supremum over the entire prefix. Local B bounds cover every visited source cell. Slow By/ByZ and fast Bphi/BZphi bounds remain separate with the original N/Pstar units.

The complete leading-prefix U/Q enclosure is combined with these prefix finite-N errors through the inherited quadratic product/C1 error formulas. The exact physical positive mass is Rm^2*(x^2-1)/2. All incoming1/N factors, new driver1/N factors, Q_ZZ/mean_ZZZ errors and quadratic cross-errors remain. No extra epsilon, hb, x, N denominator or selected phase is introduced. Signed correction centers are still unresolved; magnitude bounds are only errors.

## Observed run and report

The run exited0 in150.394 seconds; saved-report inspection exited0. It evaluates Z=-.3,0,.3 at x=1,1.25,1.5,1.75,1.775,Rh.

- 18 cumulative packets and18 exact156-term inlet replays, each producing168 terms.
- 432 serialized increment/error rows, including72 exact empty rows atx1.
- 234 pressure-source seam comparisons and72 normalized mass-row diagnostic rows.
- 300 full-prefix finite-N conditions,60 lower endpoint conditions and270 positive scalar-kernel conditions.
- 12 exact terminal leading RQy rows. Thirty32/128 versus64/256 kernel overlaps pass; five tail-kernel widths strictly improve.
- An independent finite family with normalized polynomial bumps, Rm=3 and t=.07 passes12 direct physical C0/C1 pressure integrals. Its S1 operators differentiate Qy directly and its physical integrals use independent Gauss quadrature. Three invalid queries are rejected.

The32,847,661-byte report SHA256 is `cfcb87fd7af7451dc936409beb7d6515e1eabd52f8f85373d28ff667c8dc3452`. It pins1775 inputs and83 implementations. Inspect with `python scripts/run_original_Rm_m15.py --inspect-report`; accepted ancestors do not need rerunning.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without a new agent. It independently verified the survival-tail reductions, fixed t scaling, exact native volume mass, conservative prefix envelope and C1 derivative requirements. Its initial raw-versus-normalized diagnostic objection was corrected against the actual `OriginalRmControlZ3` return values. No concrete mathematical implementation issue remained in that review.

## Next implementation queue

1. Build pressure rows for the original Rh_reference chart, native offset[-5,0], from this same four-row Rh history inlet. The actual source in `original_outer_reference_numeric.py` has leading U=M=4Z, Uy=Qy=0; use its exact analytic mean third derivative and rational Q derivatives, not an old Z2 adapter padded into missing finite-N rows.
2. Transport all four C3 finite-N rows using the original Rh_reference source cover. Preserve local B/N and the correct receipt-bound slow/fast derivative chart. Merge coincident axial covers and carry each native offset width once; keep leading functional zero identities separate from nonzero errors.
3. Append Rh..Rref to the complete168-term inlet with the exact constant-background primitive and full endpoint/prefix error budgets. Geometry is R=110*C^10*Pstar^10*exp(offset), with offset=-5 atRh. Keep that unit offset separate from huge log-radius atoms.
4. Continue original outer_slope and outer_axial, lifting genuine pressure/mean derivatives from their source/scalar/memory equations where needed. Carry their changing-source means and all cumulative pressure terms through outer buffer/o3 and compact support toRp.
5. Independently compare source-correlated cumulative pressure value/Z1 with the reserved Rp moment; establish all-five/global joins/common N, flatten/preheat/exact heat and finite energy. Then admissible stress/flat remainder, genuine n-dependent temporal recursion/smooth summation, realizable oscillatory stress cancellation and independent corrected Cartesian NS remain.
