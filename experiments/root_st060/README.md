# ST060 — quadratic residual-budget continuation

Task #892. Actual bounded computation from the complete ST059R-V checkpoint. The original full NS 1e-3 target remains UNMET. No default, visualization, publication, other-agent task or schedule was changed.

## Eight completed original independent comparisons

All three children frozen at **2026-09-20T19:35:44.867251+00:00**, BEFORE both new seeds. No later coefficient tuning or selection using holdouts. Each4096 Cartesian points, six original times, FULL original separate spatial/time/energy-quadrature ladders. Norms: full-vector maximum and fixed-time spatial volume L2=sqrt(64*mean(|R|^2)); worst time at h=.005/timeh=.0025.

| Seed | Field | Full sampled maximum | Spatial volume L2 |
|---|---|---:|---:|
|9206191|ST059R-V|.027988668119964047|.03964057938757012|
|9206191|ST060-C|.02800815224991933|.03830084384062102|
|9206191|ST060-T|.02771167137063021|.038486068368988845|
|9206191|ST060-Q|.027317018472915143|.03341285134684025|
|9206192|ST059R-V|.02782408269176644|.040307775419152904|
|9206192|ST060-C|.02789759156846098|.03920615334334708|
|9206192|ST060-T|.02763194954829432|.03915244214529982|
|9206192|ST060-Q|.027967382158514038|.03521816752894042|

Q paired L2 improves15.71%/12.63%; first maximum improves2.40%, second maximum **worsens0.515%**. C lowers L2 but slightly worsens both maxima. T has smaller L2 gains and mixed ordering against C. All8 actual scientific commands exit1 for momentum_max/momentum_L2; other original gates pass on these samples. No uniform or continuous-domain improvement claim.

## Actual implementation

Same2594stored-parameter compact axisymmetric family, viscosity.01,time[.25,.75],smooth compact u AND p,r<2/|z|<2,original fixed prescribed force,raw bounds,E0=1 and scientific thresholds. Initial spatial velocity can change with normalization differentiated.

C uses578 coupled directions; T uses737 after higher temporal modes within the existing8-term representation, numerical single-probe nullspace and whitening. Both use422800training points (80x151spatial,35times), mapped five-region quadrature and inherited finite core/pressure/shear/bias/profile/harmonic/morphology guards. New effective energy/enstrophy participation-volume ratio screen[.99,1.02] at training quadrature/times. This is not a changed physical support or weakened PDE threshold.

C/T first-order energy-budget constraints allowed only.01-.03 fractions of proposed steps after full nonlinear screening. Seven steps actually committed in each before redesign (the initial stop request named six; seventh states and amendment retained). C189.80s/T241.18s,not converged. Temporal enrichment alone did not beat C's training L2 at this budget.

Q warm-starts from C's seventh state with the enriched737 map. New supporting cuts account for the positive quadratic term in

`||R+Jd||^2 = ||R||^2 + 2R^T Jd + ||Jd||^2`.

The old constraint approximation omitted this last term, although its final full nonlinear screen did not. Q uses up to5supporting-cut exchanges for convex norms of the LINEARIZED residual, followed by unchanged nonlinear checks/backtracking. Accepted fractions.7,.3,.1,.03,.1,1.0 over6steps/209.08s. This is a finite local outer approximation, not exact conic optimization, exact pressure elimination, continuum minimax or convergence proof. Warm start,map and work budget also differ,so not a pure one-factor attribution.

## Independent structure and actual effective volume: parent/Q only

Fresh800corepoints/17times: all5direction fractions1 for parent/Q. RangeR[.035,.215],|Z|[.025,.215],not the full support domain. Fresh1025midplane probes: minimum signed shear retention.99722571 (above.995); positive bias minimum decreases.000194722 to.000171919. Final expanded drift.19082753 to.18955407,small improvement,not exact similarity.

Independent64x96quadrature at t=.25/.5/.75: energy participation-volume changes-.901%/-.526%/-.095%; enstrophy participation-volume changes+1.473%/+.165%/-.897%. These three checks remain within[.99,1.02],not an all-time guarantee. Initialenergy.9999999999994551; relative initialvelocityL2change3.645%,mostly nonuniform spatial redistribution. Support cylinder16pi unchanged. No amplitude/volume-collapse shortcut; concentration is not literally unchanged.

Fine five-slab independent quadrature: at t=.75 squared residual-energy falls29.14%in bulk|z|<=1.75,32.01%in adjacentcollars,6.08%in|z|>1.9. All25signed-slab/time pairs decrease on this finite comparison. .25/.5/.75 are shared times; .4173/.6327 are additional times. Not pointwise/all-time monotonicity.

Degree2-8 harmonic lower-bound ESTIMATES at t=.25/.5/.75: parent.00537529/.00392092/.00506504 toQ.00466022/.00305253/.00496685. Exact weak/Gram identities with floating velocity integrals,NOT interval-certified bounds; still above.001. Independent original-operator96-order L2 .03903020/.03173569/.03896470 to .03414959/.02926540/.03312392;48/72/96levels retained.

## Peak counterevidence and independent differentiation

New85x159spacegrid/20times:max.02939971 toQ.02850748. Every gridpeak independently differentiated; largest vector discrepancy1.412e-7. Not nested with old grids or a continuum bound. Second random peak Qrefines.02796738/.02797709/.02797771 versus parent.02782408/.02783289/.02783346; its small regression survives spatial refinement.

Additional post-freeze local maximization from13same starting locations found parent.03062200163 andQ.02924590045,about4.49%lower but ABOVE their finite-grid values. Qpeak atinitialt=.25,r=.28058,z=1.93583. IndependentFD at h=.000625/timeh=.000625 differs by<9.77e-9. Nine of13 local solves perfield report convergence; remaining stop records retained. No field coefficients changed, and no global maximum certificate is claimed.

## Recovery, testing and publication boundary

Every committed step has actual raw candidate, delta,map,parent,sourcebinding,damping and history. T actually stopped atstep3 and resumed in a newprocess; intermediate ZIPs were exported. C/T preservegenerations0-7,Q0-6. Frozen fit entry points now reject further fitting. New work needs a new registration/holdouts.

Main-directory40focusedtests passed69.99s,warnings-as-errors: full residual/normalization derivatives, effective-volume gradient/scaling, convex supporting-cut inequality, raw/source/force identity, recovery/corruption and frozen-phase tests. All8originalscience validations and both6-type supplementary audit sets completed. No cloudCI, nativeMATLAB, Lean or full historical suite run.

Failures retained: initial outer calibration timeout; test-module name collision; pre-validation driver import failure; duplicate read-only audit scheduling. Only own duplicate workers were stopped, their partial reports retained, completed report bytes compared/reused and verification resumed successfully. No candidate,seed,operator or threshold changed.

Actual optimizer sources and replay entry are committed: sourcehead `a5b902b2cc81aabcac0020ed9f3e409604dd2a18`. Returned Gitblob identities for both optimizers match executed local bytes exactly. Complete rawarrays, dependencies, all maps/generations and full reports are supplied in the offline conversation archive, not claimed already uploaded here. No PR/merge/defaultpromotion is claimed.

From the COMPLETE delivery:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
python experiments/root_st060/session_status.py
python experiments/root_st060/replay_st060.py --id ST060-Q --out outputs/Q --seed 9206191 --validate --structure
```
Expected scientificexit1; add--resume to verify/skip completed reports. No optimizer runs during replay.

QrawSHA256 `a865f5a7cdc1682b065769b01f822ffc6865fd3b74ce5942cb4e154c17141c3f`.
ParentrawSHA256 `9d76d3f6236e3bef78dfdbcbed30f69bb84d6190497e1b33092f96ddf27e412c`.

`pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `blowup_proved=false`.
