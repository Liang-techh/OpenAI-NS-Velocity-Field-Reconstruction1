# ST053-Q — weak-moment progress with a peak tradeoff

Task #534, stacked on ST052/#508 at b3b8bfdbe1077f9ec967d158602951997d81e17d.
Original NS1e-3 acceptance is NOT achieved. The candidate reduces a global necessary
imbalance and slightly reduces held-out L2, but both sampled maxima and a new edge-grid
maximum worsen. Retain ST052-M as a lower-peak control; no default promotion.

## Fixed field identity and original contract

Actual parent is mathematical ST052-M. All10published velocity/pressure references
were reproduced with maximum difference0.0 from its checked518modifier recipe and
ST051-B. New parent raw SHA482471204a2c2217a95dde23017a8a09ff35f4ddd8c5b31b5ae48141edaa37ae
has different metadata from original e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da.
Generated metadata is not claimed byte-identical to the original training file.

Retain nu=.01,t[.25,.75],physical R3/box[-2,2]^3,smooth compact u AND p inside
r<2/|z|<2,original two-parameter divergence-free curl force a,c in[0,10],E(.25)=1,
and all original coefficient/energy/support/divergence/core/full-momentum gates.
Both full-vector Euclidean max and fixed-time spatial volume L2 must be <=.001.
L2=sqrt(64*mean(|R|^2)), not training MSE or space-time RMS.
Same2594stored field coefficients and518bounded modifiers; no new spatial basis,
unconstrained force or physical amplitude collapse.

## New necessary-condition objective, not an alternative acceptance test

For the inherited compact u,p and divergence-free compact forcing, let
phi=(x/2,y/2,-z), E=.5 integral|u|^2 and alpha=integral uz^2/integral|u|^2.
The inherited integration-by-parts identity is

```
<phi,R> = D = integral(uz^2-(ux^2+uy^2)/2) = E*(3*alpha-1)
||phi||_L2(cylinder)^2 = 88*pi/3
||R||_L2 >= |D|/sqrt(88*pi/3)
```

Compact support removes boundary terms, phi is the gradient of a harmonic quadratic,
and its pairing with the restricted divergence-free force vanishes. Exact R=0 requires
alpha=1/3, but satisfying that fraction does NOT suffice for NS. This is an inherited
fixed-velocity necessary condition, not a new whole-family impossibility theorem.
Floating quadrature values below are not interval-certified lower bounds.

`moment_step.py` adds .2*mean(D^2)/(88*pi/3) to the inherited full-momentum objective,
and bounds |D(t)| by .9 of its parent value at15training times. Moment matrices use
48x72spatialGauss; exact reduced derivatives include initial-energy normalization.
A new deterministic24586point training pool (seed9175350) supplies563active peak
constraints, capped by the parent's measured pool maximum .03770400148555011.
The rest of the pool is checked after fitting/restoration. This is not continuous minimax.

Original finite corepressure/flow,axisB/|A|>=1.75,bias/shear(.999parent/.995fixedS),
morphology and per-training-time L2cap1.0 remain. New anchor is .08relativeM;
profile and nonpressure-axial ratios1.0. These are autonomous auxiliary restrictions.
Training uses32x48Gauss space,13Gauss times+endpoints,17old/19newconstraint times,
3423inheritededge probes. No source-calibrated shape/amplitude data are claimed.

## Execution, restoration and explicit metadata corrections

The actual prefit ST053 registration is reproduced in `registrations.json`, originally
written locally at2026-09-19T02:05:37.574472+00:00. The unchanged inherited wrapper's
registration/summary contains old profile.995 and old seed9175191/2 labels. The actual
subclass overrides profile1.0, and the explicit ST053 prefit file and actual validator
use9175391/2. Preserve old records and read the clarification rather than interpreting
those inherited labels as the executed protocol. Original scientific gates are unaffected.

Main fit reached its420second budget:68iterations/135solvercalls, minconstraint
-8.5905e-5, NOT converged or feasible. Before any new holdout, a locally registered
conditional restoration permitted<=8minimum-norm steps at unchanged1e-7tolerance.
Actual restoration:ONE correction,2innerQPiterations,53.59s,minconstraint-2.4045e-8.
Five marginal peak violators were selected, but all were already in563active points:
zero UNIQUE points added, despite inherited summary key `added_points=5`.
The final whole training pool max .03770400153196702 is within declared tolerance.
This is a feasible finite-training checkpoint, not a converged optimum or PDE acceptance.
All infeasible intermediate checkpoints and the original records are retained.

Frozen2026-09-19T02:14:50.669273+00:00 BEFORE either held-out sample. No later fitting
or coefficient selection. An initial freeze-helper import failed before any freeze
or validation; retry loaded the bootstrap import, without changing numerical formulas.

## Four actual complete independent validations

Each seed4096freshCartesianpoints,sixoriginaltimes,FULL separately varied original
space/time/energyquadrature ladders. Worst-time norms at spaceh=.005,timeh=.0025:

|Seed|Field|Full max|Spatial volume L2|
|---|---|---:|---:|
|9175391|ST052-M|.03546763313766915|.050664143969171425|
|9175391|ST053-Q|.03699893531869021|.050288249184028584|
|9175392|ST052-M|.036094676526043326|.052079420700878655|
|9175392|ST053-Q|.03801579542111926|.05121256999038296|

L2 improves ONLY.74%/1.66%; maxima WORSEN4.32%/5.32%. All4original acceptanceCLI
calls actually exit1 for momentum_max and momentum_L2; other original numerical
gates pass on these samples, not continuously. No tradeoff is relabeled a joint success.

New42x138 axial-edge grid at13times: maximum .03825950059 -> .03922228865 (+2.52%).
Peak shifts from late upper edge to t=.25,r=.28293,z=-1.90493. Every time's gridpeak
was independently CartesianFD checked, maximum vector difference6.60e-8 across both
fields. This supplement has no radial-edge scan and is not a continuum upper bound.

## What improved structurally and what did not

The inherited global weak moment is substantially reduced on independent48x72,
64x96,80x120quadrature. At80x120:
- initial axial kinetic fraction .28676918 -> .31988339 (exactR=0 target1/3);
- final fraction .29159383 -> .32527862;
- initial weak-L2 lower estimate .01455182 -> .00420326;
- final lower estimate .01395654 -> .00270061.

The necessary-component estimate drops71-81% over six times, but full L2 only .74-1.66%.
Do NOT call .0027-.0042 the achieved residual. Even this necessary condition is not
at the .001 level. Direct residual pairing matches D within1.74e-9 at finest quadrature;
that is numerical checking, not interval certification or an all-family no-go result.

Comparable242point/six-time and NEW800point/seventeen-time core grids,seed9175393,
retain inwardradial,positiveswirl,bipolaroutflow,inwardradialpressure and axialpressure
towardmidplane at ALLcheckedpoints. Finite domain R[.035,.215],|Z|[.025,.215].
NEW1025midplane probes retain signed shear .99900809relativeM and1.00020631relativeS,
no reversals. Maximum shear t=.5 increases .0030884252 -> .0030972791 (~.29%).
Positive bias remains, but pressure/bias magnitudes do not all strengthen.

Comparable final profile drift slightly WORSENS17.92011% ->17.92411%; fresh offgrid
improves19.13685% ->19.09776%. Fine axial vorticity moments generally increase, but
not an image similarity score; one initial axial second moment slightly decreases.
The largest independent fine radial second-moment ratio is1.01002514 vs autonomous
trainingcap1.01: this small quadrature/generalization miss is preserved.
No new particle integration this round. No source annular pulses, matched heat exterior,
complete source identity or blow-up proof. Not all structural metrics improve.

## Reproduction and executed checks

From repository root, no fitting:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st053/test_moment.py
python experiments/root_st053/replay_st053.py --out outputs/ST053-Q --seed 9175391 --validate --structure
python experiments/root_st053/weak_audit.py outputs/ST053-Q/candidate.json --out outputs/ST053-Q/weak.json
```

The scientific replay currently exits1. Use an empty output directory. Raw-parent
absence on GitHub falls back through the checked ST052/ancestor recipes. Raw parent
is included in the self-contained user archive. Modifier SHA,bounds,flags and10u/p
references at1e-8 are checked; local reconstruction matches EVERY stored coefficient
exactly0.0. Regenerated metadata can change raw JSON SHA.

Initial6tests passed11.01s; final7passed12.97s,warnings-as-errors. Includes normalized
moment derivative,pressure/force independence,refined whole-support pairing,energy/
support,three symbolic identities,deterministic pool and frozen/claim/hash/parent tests.
Full inherited local suite and Lean were not run. Actual cloud results are recorded
in the PR only after execution; software CI success is not scientific acceptance.

Source/recipe/results/registrations and read-onlyCI are tracked. Full raw parent/child,
registered fit and restoration histories,4fullvalidationreports,2structure/weak/edge
sets,freeze/test/execution records and hashes are delivered in the user archive.
Do not treat all older historical runs as newly recovered. No old source/default or
agent schedule altered, no external numerical solver executed.

Original local child SHA256:
cb7a85c95d13fffcd8265c441fe156c26c3b53722807ea5f4374adc566b2eb29
Modifier SHA256:
074802381a4b3afb52dd0baacbb72405670c5ffa3fb536fa3086b0123acab01c

pde_validated=false; source_correspondence_verified=false; paper_exact=false;
blowup_proved=false. Review and merge are separate from implementation delivery.
