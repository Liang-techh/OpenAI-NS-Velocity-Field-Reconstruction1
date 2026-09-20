# ST057 — initial redistribution and full coupled velocity/pressure

Task #847. Base: ST056 PR754 at `4df0221b428441b4d92afe2f994b024ee66b2146`. This bounded numerical round is completed locally. The original full NS 1e-3 goal remains **unmet**. No default, publication, existing field, other-agent task or schedule was changed.

## Actual implementation

All1728existing velocity and864pressure coefficients are available in the same2594stored-parameter compact axisymmetric family; the two original force parameters stay fixed. Initial velocity is no longer frozen. Energy E(.25)=1 is still normalized by the original representation and differentiated in the analytic adjoint. A rank35 numerical restriction holds raw velocities only at the original single core probe at25trainingtimes, not throughout the initial spatial field. Global normalization can scale that probe.

A separable tensor evaluator and adjoint implement all three nonlinear momentum components. Initial-Jacobian whitening includes velocity-pressure cross blocks. Degree2-8 harmonic weak moments enter at13Gauss times plus endpoints. Exact polynomial Gram matrices normalize them, but velocity integrals remain floating quadrature. Full residual peak/L2, pressure directions, core velocity, shear/bias, scaled profile and enstrophy-distribution safeguards are retained at finite training samples. These added choices are autonomous, not exact source-paper data. Exact constrained pressure elimination is not implemented.

Original viscosity.01,time[.25,.75],smooth compact u AND p in r<2/|z|<2,prescribed restricted force,rawbounds and all original scientific gates are unchanged.

## Frozen independent results

KS and HS frozen at **2026-09-20T13:40:44.718761+00:00**, before both fresh seeds. No later coefficient tuning. Each4096Cartesianpoints,sixoriginaltimes,FULL separate original spatial/time/energy-quadrature ladders. Fullvector max and spatialvolumeL2=sqrt(64mean|R|²),worsttime at h=.005/timeh=.0025.

|Seed|Field|Full max|Volume L2|
|---|---|---:|---:|
|9205791|ST056-J2|.0383197678607861|.047454189610801806|
|9205791|ST057-KS|.03543187902240814|.04469103275662134|
|9205791|ST057-HS|.03555208742524539|.04522906674160093|
|9205792|ST056-J2|.035374015033719994|.04585654069858168|
|9205792|ST057-KS|.037138421000019135|.04203863600252839|
|9205792|ST057-HS|.03777423473769322|.042845451439170415|

KS improves L2 by5.82%/8.33%; its first max improves7.54%,second max **worsens4.99%**. HS has the same tradeoff. All six original scientific commands exited1 for momentum_max and momentum_L2; other original sampled gates passed. No uniform improvement or continuum certificate.

New59x97axis-inclusivegrid/15times: maxima J2.03886235,KS.03527385,HS.03518402. IndependentCartesianFD checks every gridpeak; largest vector discrepancy2.02e-7. Randomseed2 finds larger child peaks; fixedpeak refinement to h=.00125 gives KS.03715482,so the regression persists. Original reports are not overwritten.

## Initial compatibility and structure

InitialKS relative velocityL2change4.963%; independent96-order E0=.9999999999993818. Best uniform poloidal/swirl scaling still leaves nonuniformL2change.06809817 out of total.07018678: actual spatial redistribution,not simply global scaling.

Degree2-8 normalized harmonic lower-bound estimates at t=.25/.5/.75:
J2 .0107617664/.0090159338/.0093481796;
KS .0057014858/.0044778818/.0055422360;
HS .0071428088/.0055389034/.0065920267.
KS reductions47.0%/50.3%/40.7%. Separate48/72/96quadrature and nineNEWinteriortimes+endpoints support the change. Estimates remain above.001 and are NOT interval-certified lower bounds. The identity applies under this compact solenoidal u/f,compactp contract. It does not explain the full current residual.

Fresh800corepoints/17times: all five flow/pressure direction fractions1. Fresh1025shearprobes: minimumsignedretentionKS.9984313/HS.9983782,above.995. Upwardbias remainspositive but its minimumdecreases. Finalscaledprofiledrift is essentially unchanged J2.19215846→KS.19210694; HS slightlyworsens. No whole-domain/source-image correspondence or blow-up proof.

## Actual execution, failures and delivery boundary

C1200iterations failed auxiliary selection and fell back to the parent. H/K each reached2200iterations,not convergence; their unpromoted endpoints violate auxiliarystructureconditions. Separately registered halfsteps plus312-variable smooth nonlinear-feasibility repairs produce KS/HS; one successful repairiteration each,about4.11/4.14s. These are training-feasible checkpoints,not converged constrained optima. Unrestricted and linear-only repair failures are retained. Two initial gradient checks failed near a hinge boundary before fitting; refined and off-boundary checks passed the same tolerance.

Final focusedlocaltests11passed7.49s,warnings-as-errors;compileall0. Six full independentreplays and all supplementary audits completed. No cloud/nativeMATLAB/Lean/fullhistoricaltestclaim.

Original rawSHA256:
J2 `31961e6c7ff58358954967308fdaa3dcda3d52dfe05c796e35d898203bf0375d`;
KS `af06787a530f5152df01c7a361033bc939d5764498e2c2d04c9cb0fcb6a10237`;
HS `f2f5fd60c14d17d652d44e928c16cf9229686130446dab8ca92d80f1e65381cf`.

**This branch currently publishes the task/result documentation only.** An attempted numerical-source file upload was blocked by the tool. Complete executed code, immutable raw arrays, failed controls, registrations/checkpoints/histories, six full validation reports, audits and a Git patch are delivered in the conversation archive `NS_ST057_three_route_progress.zip`. Do not treat this README as evidence that those code/array files are already on GitHub. No PR or merge is claimed. The original raw files are preserved, not merely hashes of missing data.

`pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `blowup_proved=false`.
