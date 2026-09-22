# ST070 — shear feasibility before wave construction

Task #1252. Manual user-authorized continuation from the actual complete ST069-M bundle. Paused schedules remain paused. No main/default/publication/visualization changes, and no new accepted NS candidate.

## Main result: two fixed-data incompatibilities

### A necessary monotone pressure-moment bound

For nonincreasing positive F on [a,b], let F(a)=M>m=F(b)>=0, I=integral(2 X F), and C=integral(F^2). Define

```
x_*^2 = a^2 + [I-m(b^2-a^2)]/(M-m).
C <= M^2(x_*-a) + m^2(b-x_*).
```

Proof: y^2-2 lambda X y is convex on [m,M]. For lambda=(M+m)/(2 x_*), its maximizing endpoint is M before x_* and m afterward. That two-level comparator has the prescribed I. Integrate the pointwise inequality. The comparator is ONLY a bound, never a velocity field.

At eta=0 the original ST069 data give:

| Quantity | Value |
|---|---:|
| Inner endpoint F | 0.8107143772822949 |
| Heat endpoint F | 0.04380938387260532 |
| Annular I | 2.061353199597424 |
| Required pressure moment C | 0.7867295538762102 |
| Monotone endpoint upper bound | 0.7322968118904365 |
| Excess | 0.05443274198577375 |

On 101 recorded eta points, 71 violate this necessary inequality. The functional bound is analytic, but the frozen profile/heat inputs are floating evaluations, not interval-certified enclosures. This diagnoses these fixed data, NOT impossibility of general NS reconstruction. Crossing the midpoint necessary-amplitude threshold (~0.277928) does not ensure the other moments or cone.

### The frozen inner seam also fails a strict cone prerequisite

At X=.25, eta=0, the original core has a=.46116232, b=.05932963, v=a+b^2/a=.46879522<2. Keeping these seam jets prevents v>2 on a neighboring collar of the proposed closed wave annulus. Moving the wave region requires separately resolving the uncovered collar, not ignoring it.

## Actual new control: ST070-S

A separately registered trial changes the autonomous heat amplitude from .25 to .35 (+40%), retaining the whole original inner core and axis pressure. Same h=.005, nu=.01, Xc=.25, Xb=4, |eta|<=.5 and local unforced physical operator. This is a DIFFERENT outer-data comparison, not a same-problem optimizer speedup; the old global fixed-force benchmark is unchanged.

F is built by integrating a nonpositive derivative with flat core/heat collars, nonnegative compact density mixtures and a small envelope density. Endpoint/angular/pressure moments determine the mixture. U uses eight compact bumps: M/J impose linear constraints and S fixes a quadratic norm in a predeclared nullspace direction. No full NS optimizer runs.

Seventeen eta interpolation nodes were the first registered level passing calibration (maximum 7.69265e-9; about 7.56 seconds). The model was frozen at **2026-09-22T21:04:20.564552+00:00**, before new independent seeds 9227091/9227092. No later parameter tuning.

On 41 new eta points and 24/40/64-point-per-segment quadrature, highest-order maximum scaled moment errors [M,I,J,S,Cp] are:

```
[7.53113691e-9, 1.47417050e-12, 5.35395588e-9,
 2.11603501e-9, 4.44844162e-12].
```

Outer leading-stress traces are actually computed, not zeroed: maxima about 1.24e-7 and 1.67e-7. These are integral/seam checks, not full NS errors or exact smooth global matching certificates.

## Shear improves, but the full cone and residual do not

Same independent 79x55 grid, X in [.27,3.98], eta in [-.48,.48]:

| Finite sample check | ST069-M | ST070-S |
|---|---:|---:|
| Positive shear a | 51.90% | 100% |
| Relaxed finite cone | 7.62% | 6.65% |
| Strict finite cone | 4.99% | 3.18% |

The new minimum analytic shear sample is only 1.46e-15. Thus 100% positive samples is NOT a robust uniform positive margin. Requiring a>1e-3 leaves 81.77% passing. An inherited large-stress directional ratio produced a divide-by-zero warning when finite differences rounded a weak derivative to zero; its log is retained. The actual finite cone uses a separately evaluated derivative and its stable algebraic form is compared to the source radical expression. Weak-margin samples remain fragile. Both mesh levels show widespread failure. No realizable waves are claimed.

The original physical Cartesian operator was evaluated on the same 12 fresh ring points, with three spatial steps per scale and fixed time step .001*tau. Finest-step sampled maxima:

| k | ST069-M | ST070-S |
|---|---:|---:|
| .4 | 147.67164654 | 421.93562966 |
| 3 | 2214.37034571 | 6306.25249721 |
| 5.5 | 29922.61982944 | 84951.03584038 |

These are local finite-strip point maxima, not whole-space suprema or volume L2. The new control is worse by roughly a factor 2.8 and is NOT promoted. Original 0.001 physical gates remain failed. This comparison does not claim a new independent time-step ladder or a global fixed-force validation.

## Actual local axis-data scan and a coarse-grid false positive

Registered controls changed only U0=m*eta+.02, for m=4,6,8,12,16,24, at two radial truncations and two candidate seam radii. A coarse 101-point seam grid suggested m=16 might satisfy v>2. Independent eta resolution and radial refinements were then checked before freezing a new LOCAL diagnostic core.

For m=16, orders24/32/40/48 and eta degrees128/128/160/192 were retained. Higher order eventually exhibits cancellation; there is no arbitrary-precision convergence claim. The frozen order48 core was tested on 4096 new points (seed9227094): leading F defect 2.54e-11, leading U defect 1.19e-7. These are not complete momentum residuals.

Location-only searches of that frozen core find U_X zeros at eta=-.3369709005, -.0012053458 and .3349164994. At these points b=0 and v=a, respectively .69409412, 1.79973986 and .69324248, all below2. Merely increasing axial slope does not remove the seam obstruction. This scan did not replace the ST070-S inner field or generate a matched global candidate.

The inherited core builder stores ST068-I as an internal label in these NEW diagnostic NPZ files. Their distinct paths, actual axis parameters and SHA256 identities are explicit in data/registry.json; they must not be confused with the original frozen ST068-I. Original bytes are preserved.

## Execution and delivery

Main focused tests: 13 passed in8.75s, warnings treated as errors. Clean complete-package tests: 13 passed in9.41s (13.43s process wall time); compileall succeeded. Clean-package highest-order 41-point moment replay actually ran in13.13s and matched all reference errors exactly. The finest-step 12-point k5.5 complete physical residual actually reran in29.06s, returned scientific exit1 and reproduced every vector exactly. A subsequent --resume verified/skipped the result in4.80s; no fit ran. Other cone/axis/truncation audits ran in the main work directory and are copied unchanged, not claimed rerun in full by that clean smoke test.

Failures retained: initial root bracket selection, pre-freeze list-abs construction/report errors, independent report list-abs error, inherited weak-shear ratio warning, coarse-grid axis false positive and high-order cancellation. No held-out retuning, changed thresholds or manufactured zero force. No native MATLAB, cloud CI, Lean or full historical test suite ran.

This branch publishes the actual moment-bound module and construction entry plus this result record. The complete monotone solver, original dependencies, raw ST070-S and diagnostic core arrays, solved eta nodes, all registrations/freezes, complete numerical reports and failures are in the conversation archive **NS_ST070_Shear_Feasibility.zip**. Do NOT claim all arrays/dependencies are already on GitHub. No PR or merge is claimed.

Original input archive SHA256: cd22bc537d6500300b53bff282e2bf72939da6cd25ee1af9b617951b2ca90746.
ST070-S raw SHA256: 4120333fba57752959d246b7d99e42efa3058a2ac6c0ac5d77e52437594867e0.

From the complete bundle:

```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --part matching --out outputs/matching.json
python replay.py --part full --k 5.5 --out outputs/k55.json
```

Moment replay currently exits0; physical momentum replay exits1. Add --resume only for identity-checked completed reports.

Next work must jointly select compatible axis swirl/axial/pressure data and exterior parameters, test the necessary moment bound and U_X-zero seam neighborhoods, then re-close five moments and finite stress-cone constraints. It still needs higher-order core dynamics and actual non-axisymmetric waves. All PDE/global/wave/source-correspondence/blow-up flags remain false; scale recursion is not newly demonstrated.
