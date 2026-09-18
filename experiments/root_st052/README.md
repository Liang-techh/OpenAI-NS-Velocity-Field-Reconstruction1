# ST052 — explicit sampled minimax with physical structure safeguards

Task #502; stacked on #460 at `4b784f1b8457af2ead49295631d834d4e882000b`.
This is an actual new frozen candidate, NOT an accepted NS solution. No existing
field, default, physical configuration, threshold or other-agent schedule changed.

## Starting point and original contract

Actual ST051-B raw parent SHA256:
`0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d`.
All182payload files of the provided parent archive passed its stored SHAmanifest.

Keep nu=.01,t[.25,.75],physicalR3/evaluationbox[-2,2]^3,smooth compact u AND p in
r<2/|z|<2,original restricted two-parameter curl force a,c in[0,10],E(.25)=1,
and all original energy/support/bounds/divergence/core/max-L2 gates. Full-vector
max and fixed-time volumeL2=sqrt(64*mean(|R|^2)) must both be <.001. No arbitrary
residual-defined force, physical amplitude collapse or modified validator.

## Actual construction and bounded fits

Use the same2594stored-coefficient actual field basis and518bounded modifiers.
Add ONE nonphysical scalar q, with explicit constraints

```
q >= |R(x,t)|^2 / parent_training_pool_max^2
```

at active TRAINING points. Minimize q plus .10 times normalized inherited loss.
The deterministic23079point training pool includes39r x57z x9times and3072random
cylinder-space-time points,seed9175250. Select48highest residuals per gridtime plus
96random-pool peaks,then rescan the entire pool and cumulatively add new peaks.
Active set sizes528 and833. A finite sampled epigraph is NOT a continuum bound.

Inherited autonomous training safeguards: initialcoreanchor.08relativeB,expanded
profile ratio1(no-worsening),signed nonpressure axial ratio1,axis swirl/strain1.75,
positive corepressure/directions,shear.999relativeB plus.995fixedST048-S,normalized
axial vorticity moments>=.9999parent,radialmoment<=1.01parent,and per-training-time
L2cap1.0. Same32x48Gauss space,13Gauss times+endpoints,17/19constrainttimes.
These are finite autonomous choices,not recovered source coefficients.

Both minimax stages hit their540second wall budgets,at67 and57iterations. Neither
converged or produced a feasible final physical/structure vector: minima were
-2.8762e-4 and-2.3013e-5. Their active epigraph bounds are not accepted outcomes.
The first stage's fullpoolpeak.0407496 notably exceeds its activebound.0281024.
The second fullpoolpeak.0282885 nearly matches its activebound,but still fails
physical/structure feasibility. Actual records and allfailediterates are retained.
Wall-stop nfev counters192/167 include internalcallback/calibration/objective calls;
they are not returned SciPy nfev or proof of optimization convergence.

### Pre-holdout repair, not relaxed acceptance

`selection_amendment.json` was written at22:49:24UTC, BEFORE fitting finished or
any new holdout was inspected. It permits at most8minimum-norm linearized feasibility
projections with unchanged1e-7tolerance,original bounds and constraints.
Actual repair used ONE correction,whose small QP completed2iterations. Constraint
minimum became-3.2051e-9. The selected fullTRAININGpool maximum is.02829035887
versus parent.05285740325. The old infeasible summaries remain separately saved.
The epigraph is recomputed from this actual field; changing q alone does not alter u/p.
This is a feasible finite-training checkpoint,not a local/global optimum.

Both fields were frozen at2026-09-18T22:58:00.668368+00:00 BEFORE both independent
seeds9175291/9175292. No tuning or coefficientselection followed those results.

## Four complete independent reports

Each seed4096fresh uniformCartesianpoints,sixoriginaltimes,and FULL original separate
space/time/energyquadrature ladders. Worsttime at h=.005,timeh=.0025:

|Seed|Field|Full momentum max|Spatial volumeL2|
|---|---|---:|---:|
|9175291|ST051-B|.058161679441807|.05134763166762801|
|9175291|ST052-M|.03681687198184801|.05179943624413842|
|9175292|ST051-B|.054758531349444965|.05213559608216369|
|9175292|ST052-M|.03697350305423294|.05246325441518328|

Max improves36.70%/32.48%,but L2 WORSENS.88%/.63%. A per-training-time no-worsening
constraint is not a holdout guarantee. All4original scientific acceptance CLIs
actually exit1 for momentum_max/L2; remaining original gates pass on these samples.
No uniform pointwise improvement,accepted1e-3solution or automatic defaultpromotion.

Independent81x121axis-inclusivegrid/six times:max .05356755 ->.03729748 (~30.37%lower).
Separate finer axial41x102 and radial31x31edgegrids at13times:max .05816755 ->.03825784
(~34.23%lower),still at late upper axialedge. Each gridpeak is independently checked
with CartesianFD. Random/grid maxima are NOT continuumupperbounds. FullTRAININGpeak
.02829 is lower than the independentedgepeak.03826: exchange did not certify minimax.

Fixed randompeak refinement to h=.00125 yields child maxima.03682303/.03698382,
so the remaining residual does not vanish with differential error. A separate fresh
4096point random space-AND-time check,seed9175294,max .03210344 ->.02598123; its RMS
is a separate space-time metric,NOT fixed-timeL2 or an acceptance replacement.

## Structure and real particle paths

On comparable242point/six-time and NEW800point/seventeen-time scaledcoregrids,
all tested inwardradial/positiveswirl/bipolaroutflow/inwardradialpressure/axialpressure-
towardplane directions remain correct. RegionR[.035,.215],|Z|[.025,.215]; no global
statement near a biased midplane or across the entire support.

NEW1025shearprobes:min signedretention.99900194relativeB and.99713380relativefixedS,
no reversals. Maxshear at t=.5 rises.0030389243 ->.0030884252 (~1.63%). Minimum sampled
upward midplane speed remains positive. But expandedprofiledrift slightly WORSENS:
comparable17.87595% ->17.92011%,newgrid19.47162% ->19.50141%.

Fine64x96vorticityquadrature: finalnormalized axialsecondmoment grows2.57%,fourth3.77%,
axialRMSextent1.28%,and radialsecondmoment also grows~1%. This is redistribution,not
comparison to an actual OpenAI targetfield or proof alltailfeatures are preserved.

Actual x'=u(x,t) trajectories are integrated on36newinitialseeds with two DOP853
accuracy/step settings. All sampled paths retain radialcontraction,positive turn,
axialoutwarddisplacement offplane and upward displacement onplane. Mrotation is only
.0498-.0649radians; contraction/rotation slightlyweaken versusB. Integratoragreement
<8.3e-15 is not an independent-integrator test or sourcecorrespondence certificate.

## Remaining nonlocal constraint

The inherited weak identity is remeasured,not newly claimed as a proof:
<(x/2,y/2,-z),R>=integral(uz^2-(ux^2+uy^2)/2)=D.
On the compact support, ||R||L2 >= |D|/sqrt(88pi/3). For M the floating96point
quadrature estimates range .01372-.01455 acrossvalidationtimes; not intervalcertified.
At t=.25 axial kineticfraction is.286769 versus exactR=0 target1/3. This fixedvelocity
still has a global imbalance that pressure-only fitting cannot fix. It does NOT
prove the whole admissiblefamily infeasible or justify changing the forcing/threshold.

## Tests, persistence and replay

Seven final localtests passed16.85s,warnings-as-errors; compileall succeeded.
Includes epigraph and completeCartesian residual derivatives,energy/support/rotation,
deterministiccumulativeexchange,three symbolic identities and frozenrecipe mutations.
Full inheritedsuite andLean notrun. Actual cloud outcome is recorded in the PR,
not assumed from workflow existence.

From repository root, without fitting:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st052/test_minimax.py experiments/root_st052/test_frozen.py
python experiments/root_st052/replay_st052.py --out outputs/ST052-M --seed 9175291 --validate --structure
```

The last command currently returns1 for scientificfailure. It refuses a nonempty
outputdirectory. The legacy structure JSON key `minimum_signed_ratio_to_P` denotes
this round's ST051-B comparator,NOT ST050R-P; parenthash disambiguates it.

Recipe verifies immutableparent,modifierSHA,bounds,scientificfalseflags and10u/p
references at1e-8. On GitHub it reconstructs parent via ancestors; the self-contained
userarchive has actualrawparents. GeneratedJSONmetadata differs,so rawSHAidentity
is not asserted for regeneratedfiles. Refit scripts require the published actual
rawparent and saved evidencepaths,while frozenreplay needs no training.

GitHub has newreadablesource,exactfrozenmodifiers,summaryandread-onlyCI. The userarchive
additionally has actualrawfield,allregistrations/checkpoints/failedstagehistories,
repairamendment/results,freeze,fourfullvalidationreports,structure/dense/fineedge/
fixedpeak/spacetime/particledata and per-filehashes. No old historicalarray recovery
or externalnumericalsolver execution is claimed.

Raw child SHA256:
`e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da`.
Modifier SHA256:
`4d8b6e92f47151a7eed112985779a58fdf1188bd30c35c3a59263910629f1be9`.

pde_validated=false;source_correspondence_verified=false;paper_exact=false;
blowup_proved=false. Original1e-3 target remains UNMET.
