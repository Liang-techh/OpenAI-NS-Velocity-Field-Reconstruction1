# ST048: boundary peaks and signed shear preservation

**The original full NS 1e-3 target is NOT achieved. No default is promoted.**

Issue #384; stacked on #379 at `2c51cd20e036ab29954cf77a70814bc918d5c6a5`.
This round starts from the actual frozen ST047-E field, not ST006.

## Fixed physical and scientific contract

nu=.01, t in [.25,.75], physical R3, evaluation box [-2,2]^3. Smooth compact
velocity AND pressure in r<2 and |z|<2. Original two-parameter curl force with
a,c in [0,10], E(.25)=1, original bounds, energy, support, divergence, core signs
and scaled single-probe drift constraints. Both full-vector momentum max and
spatial volume L2 must be <=.001. No arbitrary force or amplitude collapse.
L2 is sqrt(64*mean(|R|^2)) separately at each original time, not training MSE.

## Actual implementation and fits

`boundary_shear.py` reuses the ST047 coupled objective and solver-coordinate
conditioning. The underlying field remains the same 2594-coefficient basis.
S uses398 reduced variables (5x5 spatial directions); B uses530 (6x6). Both
use4 correction time modes and start independently from ST047-E.

Add3423 axial/radial support-edge training probes including the axis and time
endpoints. A smooth logsumexp squared-residual peak term (weight.04,
temperature.0001) and .01mean edge residual squared complement the existing
full momentum objective. Training still uses32x48 spatial Gauss nodes and
13Gauss time nodes plus endpoints.

New actual signed radial-axial shear constraints at7physical radii x17times:
`sign(parent_shear)*child_shear >= .995*abs(parent_shear)` and an upper limit
`1.25*abs(parent_shear)+1e-5`. They use partial_r uz=2r*C_s and include the
normalization derivative. They do not constrain pressure. These finite probes,
strengths and ratios are autonomous choices, not recovered source data.

Original constraints and previous profile/pressure/nonpressure-Mz safeguards
remain. A best feasible TRAINING checkpoint is selected at declared1e-7
feasibility tolerance, with modifier/raw-parameter bounds separately checked.
Checkpoints are saved every5iterations.

S:196iterations/222function evaluations, optimizer success, minimum training
constraint -3.69e-10. B:300iterations/339evaluations, **iteration limit reached**,
NOT optimizer convergence; selected feasible checkpoint minimum -1.06e-10.
No local/global-optimality or feasibility-of-all-possible-fields claim.

Both were frozen at2026-09-18T11:05:01.762940Z, BEFORE all independent results
below. No fitting or coefficient selection followed the held-out checks.

## Paired full independent validation

Each seed:4096new Cartesian points, six original times, full original separate
spatial/time/quadrature ladders. Worst time at space step.005/time step.0025:

|seed|field|max|spatial volume L2|
|---|---|---:|---:|
|9174801|ST047-E|.0539967628896|.0598695387759|
|9174801|ST048-S|.0487906430674|.0546078095035|
|9174801|ST048-B|.0520067670456|.0505474714611|
|9174802|ST047-E|.0479178772475|.0640575108889|
|9174802|ST048-S|.0445511980112|.0581197199193|
|9174802|ST048-B|.0430654628439|.0536919394715|

B improves paired L2 by15.57%/16.18% and max by3.69%/10.13%. S gives a better
max on seed1; B is not uniformly superior. All SIX original scientific gate
commands actually exit1, failing momentum_max and momentum_L2. Other original
gates pass on these samples only. None is a continuous space-time certificate.

Axis-inclusive81x121/six-time grid maxima:parent.05304220,S.04970658,
B.04468835. But randomseed1 finds B.05200677 near z=-1.9363 at the initial
time, which that grid misses. Never present the smaller grid maximum as a
continuum bound or claim that the boundary peak is eliminated.
Separate fixed-peak spatial refinement yields Bseed1 .05200677/.05206546/
.05206929 at h=.005/.0025/.00125 while divergence falls to1.64e-8. The original
reports are preserved. This is fixed-point checking, not a new global search.

## Actual structure and new checks

On the same comparable242point core grid at t=.75, profile drift falls from
18.25043% to about18.1435% in BOTH children; only a small improvement. Adverse
axial-pressure RMS falls from.04940861 to.04831697(S)/.04872070(B), but axial
pressure still points OUTWARD at every checked core point/time. Inward radial
flow, positive swirl, bipolar outflow and inward radial pressure remain.

At t=.5, max|partial_r uz| rises .00218244 ->.00219328(S)/.00218925(B).
NEW651off-grid midplane probes (seed9174881) find minimum signed shear ratios
.99500735(S)/.99500228(B), no sign reversals and no sample below.995. This is
finite evidence, not a continuum bound. The older722point structure grid is
reused; its inherited label fresh_offgrid does NOT make it newly fresh here.

Actual nonautonomous x'=u(x,t) trajectories on36seeds retain inward contraction,
positive rotation, bipolar displacement and upward midplane displacement.
Rotation and contraction are weaker than the parent, so not all structural
measures improve. B rotations are only.06838-.09982 radians, not multiple turns.
No source annular pulses, matched exterior, source identity or blow-up proof.

`extra_audit.py` additionally tests4096new interior space-time points, with
random time per point (seed9174882) and independent vector-time CartesianFD4.
Parent/B sampled maxima .04032820/.03510318 and volume-scaled SPACE-TIME RMS
.04339766/.03717360. This RMS is NOT the original fixed-time spatial L2 gate.
Analytic/FD discrepancy stays below1.8e-6. No result is used for further fitting.

## Tests and replay

Final local tests:17passed in10.03s with warnings as errors. Three symbolic
identities verify divergence, shear chain rule and positive-energy normalization.
Frozen reconstruction reproduces both stored coefficient arrays exactly locally.
Full inherited suite and Lean are not run. Cloud results are recorded in the PR
only after actual execution; a green replay means faithful reproduction of the
FAILED scientific result, not acceptance.

From the repository/package root:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st048/test_boundary_shear.py
python experiments/root_st048/replay_st048.py --id ST048-B --out outputs/ST048-B --seed 9174801 --validate --structure
```

The last command currently exits1. Use ST048-S for the other frozen field;
nonempty output directories are refused. No refitting is needed. Complete fit
commands and histories are in the accompanying archive's registrations/logs.

The source checkout reconstructs ST047-E through the already-checked ancestor
recipes. The full archive contains its exact raw JSON instead. Frozen modifiers
are SHA256-checked .npy data, decoded with allow_pickle=False; numerical u/p
reference values are checked at1e-8. Regenerated metadata differs, so original
raw-file byte identity is NOT claimed for regenerated JSON.

Raw original SHA256:
- parent:dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0
- S:6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09
- B:b599856cf191be10d554bdfa260328c089c96875142b41466b2c601a44ba3262

GitHub tracks readable code, frozen recipes and a result summary. Complete raw
fields, registrations, all checkpoints/histories, six primary reports and all
supplementary checks are in the user archive. No claim that all historical
research archives are tracked. Existing source/defaults/configs are untouched.
