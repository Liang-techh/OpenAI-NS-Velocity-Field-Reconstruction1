# Original active Rm patch pressure with C3 finite-N histories - 2026-10-16

`OriginalRmPressureRows` now evaluates the actual partial-bump Rm-to-Rh patch from the genuine physical source and implicit control jets. It supplies U/Q through axial order2, the physical mean through order3 and Uy/Qy through order1, with separate complete finite-N absolute error bounds. This removes the missing pressure-derivative interface needed for the next active-patch moment integral.

The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. The accepted cumulative m15 inlet remains156 terms atRm. The active Rm-to-Rh signed m15 integral is still open. These spatial derivatives are not genuine temporal n>=2 recursion. Global all-five closure/common N, exact heat/energy/stress/flat remainder, oscillatory correction and full corrected NS remain open.

## Actual physical pressure fields

The adapter shares the genuine `OriginalRmControlZ3` source and contexts. It uses physical controls under `point/truth` and physical defects under `source.point/truth.defects`; normalized defect jets do not replace the latter. The original normalized compact bump has radius1/40 and centers5/4,3/2,7/4. Full and actual partial weight integrals, the original nonlinear five-control map and all three bump derivatives remain unchanged.

With x=R/Rm, y=log(x), g the two axial bumps and d0_partial the mass defect after partial integration,

```
U = 4Z + g,
M = 4Z + d0_partial/x,
Uy = x*g_x,
My = U-M.
```

Writing A=1-delta*Z^2 and D=1-Z^2, the original pressure recovery is

```
Q  = [2ZU -(1-delta)ZM -D*M_Z]/A,
Qy = [2ZUy -(1-delta)Z(U-M) -D*(U_Z-M_Z)]/A.
```

The same third source/control rows provide M_ZZZ and Q_ZZ. Taylor coefficients use D_Z^k/k! internally; serialized ordinary derivatives multiply by k! once. `R_Q_y_endpoint` keeps signed leading Qy, its complete finite-N coefficient error and the physical log-radius factor. Atomic log-radius terms retain the separate native offset; the physical radius is never materialized.

At exact x>=71/40 every compact bump and its derivative vanish. The unique full-weight implicit branch gives all five leading residuals zero, U=M=4Z, Uy=Qy=0. These exact source identities are applied before interval subtraction loses their correlation, including the exact support endpoint1.775. The unrefined full-map intervals remain diagnostics; zero containment alone does not establish closure. Finite-N incoming and newly generated histories remain nonzero at the terminal endpoint.

## Four genuine finite-N rows

The inlet reads the raw `owner._correction_logs` result on the genuine postrestore cover atRm, before the old adapter truncates its public serialization. All four ordinary Z rows for m,h,k,e,p already include1/N. They are not padded from the old two-row interface and are not divided by N a second time.

For every original native cell, the new local source is `actual_density_ranges.N_scaled` divided by N once. Coincident radial cells from overlapping axial covers are merged. A partial query clips its final cell exactly, and each width is log(right/left), so dy=dx/x occurs once. Rates remain1,1.5,1.5,1,0; the pressure history uses positive mass=width, while other histories use -expm1(-rate*width)/rate. Positive sums are enclosed directly in logarithms without materializing the enormous original N or exponentially small errors.

Physical pressure errors use deltaU=Pstar*B/N and deltaM=Pstar*history_m. The actual stored slow derivatives By,ByZ carry1/N; the fast Bphi,BZphi terms in Uy do not. No extra x, hb, epsilon or N factor is inserted. The pressure coefficient recovery uses the required U_Z/m_Z rows and retains the full mean third-derivative error for Q_ZZ. Magnitude bounds never define a signed finite-N field or a phase sample.

`uniform_background` provides a conservative leading-field enclosure over the whole native prefix. Positive partial weights lie between0 and their full integral; the explicitly negative weight row uses the corresponding signed hull. Beta and beta_x use their original sup bounds. This is only the leading enclosure. Subsequent total-field volume budgets must add the finite-N history/local/slow-fast errors separately.

## Observed execution

The original-source run exited0 in127.338 seconds; saved-report inspection exited0. It evaluates Z=-.3,0,.3 at x=1,1.25,1.26,1.4,1.5,1.75,1.775,Rh.

- 24 original packets;336 serialized pressure/mean rows with complete approximation errors and336 finite-N error rows.
- 42 original postrestore inlet seam comparisons,90 independently saved lower-pressure overlaps and three exact four-row inlet-history replays.
- 480 genuine finite-N history rows,630 positive native weights and24 merged-width conditions, including both axial covers at Z0.
- 432 conservative leading-prefix checks; six terminal branch identities,18 retained nonzero terminal Q error rows and120 unrefined full-map residual diagnostics.
- An independent finite fixture executes the original `patch_fields` source and compares78 pressure/field Taylor coefficients. Symbolic native mean and Qy derivative identities and five independently integrated Duhamel kernels pass. Five invalid queries are rejected.

The3,142,819-byte report SHA256 is `7d4ba34214321f03521cb04858542e6afa70436db458f93685dad1a0d42111dd`. It pins1773 inputs and77 implementations. Inspect with `python scripts/run_original_Rm_pressure_rows.py --inspect-report`; accepted ancestor computations do not need rerunning.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused; no new agent. It confirmed the raw four-row inlet equivalence, native derivative units, rate-zero/partial transport and the leading-only scope of the uniform prefix enclosure. The external recovery checkpoint retains terminal attempts and the exact final run handle.

## Next work: append only the active interval

1. Build an active-patch moment adapter from these actual pressure rows and the full156-term `OriginalRestoreM15` endpoint atRm. Copy that inlet exactly once; do not reset or replace its accumulated signed/error sectors.
2. Evaluate the original identity delta_m15=[R^2 F]+integral R^2 G dy-[R Qy]. The existing F/G helpers already include the half-S1 normalization; do not add another factor1/2. Integrate dy=dx/x on1..x and keep all source/quadrature errors.
3. Keep fixed Rm radius atoms separate from native x powers. Retain both signed leading endpoints and both finite-N endpoint errors, including RQy. Use the conservative prefix field only for leading magnitude bounds; add the complete finite-N C3 and slow/fast derivative budgets before claiming total C1 errors.
4. Continue original outer/compact pressure and cumulative m15 toRp. Compare a source-correlated independent reserved moment in both value and Z1. Selected-N signed centers, useful relative cancellation and all-five/global closure are still unresolved.
5. Complete common-global-N joins, flatten/preheat/exact heat and finite energy, admissible stress/flat remainder, genuine n-dependent temporal recursion/smooth summation, realizable oscillatory stress correction and independent corrected Cartesian NS.
