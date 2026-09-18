# ST051 — pressure-aligned continuation, not full NS acceptance

Task #453, stacked on ST050R/#448 at dc1d5476ea9979566b774293add76196083288c5.
The original 1e-3 target is NOT achieved. No default, old field, physical configuration,
threshold, other agent task or schedule is changed.

## What is actually new

Start from the real ST050R-P mathematical field, not its lower-L2 but pressure-wrong
C counterpart. P's raw SHA256 is f5b5a2869991387e503a522b599af86ea814ba02bdae2f8b56fd16674f7dab19.

A retains398 reduced modifiers; B adds120 high spatial support-edge directions for518.
Both use the SAME2594 stored-coefficient asymmetric axisymmetric field family. The
extra poloidal directions are odd in z, giving zero direct unnormalized midplane axial
velocity/shear. The energy normalization still couples their coefficients to the total
field; the code includes that derivative. No ideal unrepresented field is evaluated.

Both fits retain original nu=.01, t[.25,.75], R3/evaluation box[-2,2]^3, smooth compact
u AND p within r<2 and |z|<2, original two-parameter divergence-free curl force with
a,c in[0,10], E(.25)=1, original energy range, coefficient limits, support, divergence,
single-probe core signs/drift and full-vector momentum max/volume-L2 thresholds .001.
The L2 estimator is sqrt(64*mean(|R|^2)) separately at each time, NOT training MSE.

## Registered autonomous optimization choices

The core anchor is .12 relative to P. This is not the old .08 relative to S or the .60
used when constructing P. A per-training-time residual cap1.0 relative to P is restored.
Poisson weight .003, axis swirl/strain lower bound1.75, signed axial pressure target0,
expanded profile ratio.995, nonpressure axial ratio1, axial vorticity moments>=.9999
parent and radial second moment<=1.01parent, inherited signed-shear ratio.999relativeP.

A new950point training core uses R in[.035,.215], |Z| in[.025,.215],19times: signed
axial pressure margin1e-5, radial p_s margin1e-7 and scaled velocity direction margins
1e-7. A new247point physical midplane training grid requires signed shear>=.995 of
immutable ST048-S. These grids and margins are autonomous, not recovered source data.

Use32x48Gauss space,13Gauss time nodes plus endpoints and17inherited constraint times;
3423edge probes, doubled edge objective weight.08, smooth peak weight.03. The same
bounded SLSQP/whitening machinery is reused. This is not a controlled single-factor
ablation of only the added edge modes because their budgets also differ.

A:140iterations/185calls, optimizer success, minconstraint -3.58e-11,340.17seconds.
B:138iterations/225calls, registered wall-budget checkpoint at651.98seconds,
NOT optimizer convergence; selected feasible training minimum -1.50e-8 under declared
1e-7 tolerance. Both save actual raw candidate and modifiers every5iterations.
All training finished before the freeze receipt2026-09-18T18:16:45.276154+00:00.
No coefficient retuning or selection after either held-out sample.

## Eight actual full independent reports

Each fresh seed has4096Cartesian points, six original times and the FULL original
separate space/time/quadrature ladders. Worst-time values at h=.005, timeh=.0025:

|Seed|Field|Full vector max|Spatial volume L2|
|---|---|---:|---:|
|9175191|ST050R-P|.0533163081|.0610288898|
|9175191|ST050R-C|.0470099142|.0501130512|
|9175191|ST051-A|.0560898350|.0584628154|
|9175191|ST051-B|.0521649879|.0507727319|
|9175192|ST050R-P|.0437438555|.0608276747|
|9175192|ST050R-C|.0425168800|.0515825730|
|9175192|ST051-A|.0447198112|.0576431307|
|9175192|ST051-B|.0563478753|.0511819394|

A lowers L2 by4.20%/5.24% versus P but worsens both maxima. B lowers L2 by16.81%/15.86%
but its second maximum REGRESSES28.81% (first improves2.16%). B's L2 is close to C,
while its sampled maxima remain worse than C. Do not combine the best measures from
different fields, seeds or times into one fictional success. All8original scientific
acceptance commands actually return1 for momentum_max and momentum_L2. Other original
numeric gates pass on THESE samples, not a continuum certificate.

Independent81x121axis-inclusive grid/six times: maxP .0500781370, C .0498248921,
A .0482586448, B .0535675542. The largest B grid peak moves to t=.75,z=1.9285.
Independent CartesianFD agrees at every grid peak within5.60e-8. Grid maxima are not
upper bounds; random samples find larger maxima. Fixed-random-peak differentiation
refinement confirms B's second peak .05634788/.05635777/.05635840 as h halves from.005
to.00125. The regression does not disappear with numerical differentiation error.

A separate fresh4096point space-AND-time sample, seed9175194, yields maxP .04186243,
C .03619570,A .04138116,B .03841111. Its space-time RMS is a distinct auxiliary
quantity, not fixed-time volume L2. All four peak vectors were independently checked.

## Actual structural progress and counterevidence

On comparable242point/six-time and NEW800point/seventeen-time scaled core grids,
A/B retain inward radial flow, positive swirl, bipolar outflow and inward radial and
axial pressure at every checked point. On the SAME newgrid, parentP has worst-time
pressure-correct fraction98.25% and bipolar fraction97.875%; these misses are absent
in both children. Minimum signed pressure gradients on the newgrid are positive:
A1.5388e-4, B1.2628e-4. This is finite-grid direction evidence, not full-domain alignment.

NEW1025physical midplane shear probes pass both thresholds for A/B: minimum signed
retention relative P is1.00197/1.00201, relative S .995497/.995533. No sign reversal.
Maximum shear at t=.5 increases from P.00273030 to A.00301649/B.00303892 (~10-11%).
However minimum upward midplane speed decreases .000593763 ->.000225717, still positive.
No source-calibrated amplitude target is available; not all features strengthen.

Comparable final profile drift P17.78626%, A17.86942%, B17.87595%: slightly WORSE.
New offgrid drift P19.36467%, A19.12982%, B19.04861%: improves. Do not claim uniform
profile convergence or compare absolute percentages from different grids.

Fine64x96vorticity quadrature at t=.75: B's normalized axial second/fourth moments
grow12.23%/22.43%, axial RMS extent grows5.94% versus P. This is actual vorticity
redistribution, not target-image similarity. Small autonomous radial-cap misses remain:
max fine second-moment ratio A1.01000117057/B1.01000439494 versus1.01; original sampled
energy quadrature gate is a separate test. No interval certification.

New instantaneous zero-sheet diagnostic: at13scaled radii and3times, both children
have p_z=0 very near Z=0 (within7.09e-5) and u_z=0 at negativeZ between about-.025 and
-.0133. The parent u_z-zero extends to-.0403. These are instantaneous sign zeros,
NOT material/Lagrangian dividing surfaces. Positive midplane velocity implies nearby
negative-z points can also move upward; the finite bipolar checks exclude |Z|<.025.
No source annular pulses, independently matched heat exterior, full correspondence,
continuum NS certificate or blow-up proof is implemented.

## Executed tests, identity and replay

Final8focused tests passed15.14s with warnings as errors. Full momentum and constraint
Jacobians, independent Cartesian shear, energy/support/rotation, both recipes and
claim/hash/parent mutations are checked. Four local symbolic algebra identities and
the direct edge-parity check are not a complete PDE proof. Full inherited local suite
and Lean not run. CI execution is recorded separately in the PR, not assumed here.

From the repository root, without fitting:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st051/test_aligned.py
python experiments/root_st051/replay_st051.py --id ST051-B --out outputs/ST051-B --seed 9175191 --validate --structure
```

The last command currently returns1 for the failed original momentum gates. A nonempty
output directory is rejected. Use ST051-A for the other frozen field.

The self-contained user archive contains raw P/C/S parents and A/B children, actual
registrations, both histories/checkpoints, eight original validation reports and all
supplementary audits, freeze/exit receipts and per-file hashes. GitHub tracks readable
new source, both exact modifier arrays/reference values and the compact result summary;
it does NOT silently claim every old historical run is uploaded. No old missing ST050
arrays were used. Ancestor recipes reconstruct the parent on GitHub, with8u/p reference
checks at1e-8. Reconstructed JSON metadata differs from local raw files; field equivalence,
not raw-byte identity, is claimed for generated recipes.

Original local raw child SHA256:
A f34f021228f5f769c9a14c6030f1b29e50ec76bde4e443de2139df46c036236e
B 0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d

Source and physics boundaries remain explicit. pde_validated=false,
source_correspondence_verified=false, paper_exact=false, blowup_proved=false.
