# Genuine original Rm source and implicit controls through Z3 - 2026-10-16

`OriginalRmControlZ3` restores the missing third axial Taylor rows of the leading Rm defects and implicit controls. It evaluates the actual private six-history transport before the public Z2 truncation. The admitted axial5 archive is an independent comparison enclosure; its interval endpoints or midpoints never supply the new source/controls.

The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. The cumulative pressure inlet remains the published156-term state atRm. This milestone prepares genuine higher source data; active-patch pressure rows, finite-N C3 transport and the new Rm-to-Rh pressure integral are still unimplemented. Global all-five closure/common N, exact heat/energy/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open.

## Private source path

The existing R110 implementation already has order3 internal jets for H,Mpsi,Kpsi,Apsi,B,C. Its public packets truncate these to2. The previous long/reference public APIs then truncate the histories again. The new source path uses the genuine internal order3 jets directly:

```
R110 private six3
  -> original full finite long variable-source backward kernels3
  -> reference positive relaxation3
  -> actual restoration mean/mixed/square kernels3
  -> postrestore decay3 atRm.
```

Targets/rates are unchanged: H/Phi tends to1.25 and Kpsi/Phi to1.25*Psi at rate1.6; B/Phi^2 tends to5/6 at1.2; C/Phi^2 tends to5 at.2; Mpsi tends toPsi and Apsi toPsi^2 at1. All finite incoming memories remain. The long kernel uses the original variable sigma source and Bell-polynomial derivative enclosures; its constant-rate point kernel is only a proxy with complete directed error.

The original changing j baseline is carried by gamma1,gamma2,Kalpha/Phi,Dpsi, with the same restoration/postrestore factors as the existing source. Scalar gamma histories are promoted to constant Taylor jets because they depend only on the radial coordinate, not Z. This analytic constant-jet promotion does not fill a missing spatial derivative.

For centered angular H/B/C memories, subtract the exact equilibrium at the genuine long endpoint and combine its coefficient magnitude with exp(-rate*(Lref+2)) in logarithms before any cap. The original `_small` condition proves the exp(-64logP) error bound. A generic1e-1000 cap would become far too broad after division by the small angular normalization unit; it is not used for those centered defects.

All five leading defects are then evaluated through Taylor order3 in the unchanged original units. P0 is excluded from the defect map. Finite-N perturbations stay separate and never determine these leading controls.

## Same original nonlinear inverse

The diagonal units are the fixed source axial t for the first two rows and exp(-2logP)*10^10*t^2 for the last three. Normalized alpha remains exp(6/5)*(1+Z^2)^2/10^10. The strict C0/C1 branch and original Z2 inverse are reused on the genuine lower rows.

For the third coefficient, temporarily set only h3 to zero in the full original quadratic Taylor map. Let b be the negative preconditioned d3+Qknown3, J the same full implicit Jacobian, R its fixed preconditioner and K=I-RJ. With the checked norm q<1,

```
h3 = b + K*h3,
||h3|| <= ||b||/(1-q).
```

Twelve directed intersections enclose this unique third row. The full order3 residual is diagnostic; existence and identification with the original control derivative come from the same branch theorem and the strict inverse. The point control uses the original root solve plus the same-Jacobian coefficient recurrence with genuine source derivatives. No finite differences or archive midpoint define the derivative.

Taylor coefficients are D_Z^k/k!; ordinary third derivatives require multiplication by6 exactly once. New output contains normalized defects/controls, physical controls, all six Rm histories and the physical mean through this order, plus exact source tuples and full errors.

## Acceptance receipt

The original-source run exited0 in75.846 seconds; saved-report inspection exited0. It evaluates genuine third-source/control data at Z=-.3,0,.3.

- 180 serialized coefficient/error conditions;80 independent axial5 defect overlaps and80 control overlaps, including both neighboring source cells at Z0.
- 45 retained lower-order defect overlap conditions;72 genuine six-history Taylor rows throughZ3.
- All15 third defect and15 third control enclosures are retained without zero-filling. This does not assert every exact third derivative is nonzero.
- 12 physical mean-to-pressure source seam comparisons; three strict third inverse conditions;36 centered angular coefficient-before-decay proof conditions.
- An independent known-cubic-control fixture executes the original nonlinear source and encloses60 known Taylor coefficients, with factorial conversion once. Two invalid axial queries are rejected.

The1,008,651-byte report SHA256 is `e893d24ef78b65f7b53b80b17d530d8ad7353eeb74eb5646262e1f3fcf94e616`. It pins1771 inputs and75 implementations. Inspect with `python scripts/run_original_Rm_control_Z3.py --inspect-report`; do not rerun accepted ancestors.

An early synthetic inverse fixture was outside the original common unit-box perturbative regime; its synthetic nonlinear weights were rescaled without changing the production source or parameters. The first original-source attempt then exposed an unsupported mpmath interval-scalar/Taylor addition. Scalar gamma histories are now promoted to analytically constant Z3 jets, matching the established centered-frame convention. Both terminal attempts are retained in the external recovery checkpoint.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without spawning another agent. It found the internal R110 order3 source path and confirmed the lifted source rates, baseline kernels, angular centering and same-Jacobian third recurrence; no concrete underbound was found.

## Next actual pressure work

Implementation unit note: `coefficients()['point'/'truth']` are physical controls; `normalized_point/normalized_truth` are normalized controls. `defects_point/defects_truth` are normalized defect jets, while `source['point'/'truth']['defects']` are physical defect jets. The old patch field reconstruction expects physical defects; use the latter when lifting it.

1. Build the Rm patch field/pressure adapter from `OriginalRmControlZ3.coefficients(Z)`, not the old Z2 public controls. Lift the existing exact partial-weight reconstruction to order3: keep U/m/own moments throughZ3, recover Q throughZ2 from the true mean derivative, and Uy/Qy throughZ1 from the actual bump derivatives.
2. Keep the same original x=R/Rm, beta normalization/support, five implicit controls, exact partial/full bump weights and dx/x measure. Exact terminal identities follow from the same unique implicit map; zero containment alone is not closure.
3. Extend the finite-N history transport through all four genuine C3 source rows. Preserve actual slow By/ByZ and fast Bphi/BZphi bounds for Uy/UyZ; fast phase terms have no extra N denominator or hb factor. Source ranges are errors, not signed field values.
4. Append only the active patch interval to the complete156-term moment inlet using [R^2F]+integral R^2G dy-[RQy]. Preserve all signed and finite-N endpoints, source/quadrature errors, phase/native measure and full upstream history.
5. Continue original outer/compact pressure toRp and independently compare source-correlated value/Z1 with the reserved pressure moment. Complete all-five/global joins/common N, original flatten/preheat/exact heat and finite-energy tail.
6. Establish admissible stress/flat remainder; implement genuine n-dependent temporal recursion and smooth summation, mean/oscillatory stress cancellation and independent corrected Cartesian NS validation.
