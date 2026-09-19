# ST054 — full compact-pressure completion, fixed velocity and force

Task #626; stacked on #536 at fe63ec7fa197faaa4ba4e51f24b1b439411e9014.
Original NS 1e-3 target remains UNMET. Neither candidate is promoted to a default.

## Actual new result

All 864 existing pressure coefficients are optimized, instead of only the pressure
directions exposed by the previous reduced velocity/pressure model. The 1728 velocity
coefficients and two restricted-force coefficients are bit-for-bit unchanged relative
to each reconstructed parent. The actual field still has 2594 stored coefficients.
Pressure is compact in the same basis; no new velocity or arbitrary force is invented.

Fresh independent full comparisons, after all fields were frozen at
2026-09-19T11:46:37.438480+00:00, with no retuning or selection on either holdout:

| Seed | Field | Full momentum max | Spatial volume L2 |
|---|---|---:|---:|
|9175491|ST052-M|.03640749110968585|.05281790944039198|
|9175491|ST054-M3|.036407650088919226|.04832320460536207|
|9175491|ST053-Q|.03735078185041553|.05215427233401148|
|9175491|ST054-Q2|.037302189295728934|.04761902101078046|
|9175492|ST052-M|.03705927465615278|.051864030548944624|
|9175492|ST054-M3|.03712485619372925|.04789220830933563|
|9175492|ST053-Q|.038473126037781864|.050921843787282615|
|9175492|ST054-Q2|.038455392348817515|.0477042271226719|

M3 L2 improves8.51%/7.66%, but maxima worsen0.00044%/0.177%. Q2 L2 improves8.70%/6.32%,
while maxima improve only0.130%/0.046%. No uniform pointwise improvement or peak breakthrough.
Each seed uses4096 new Cartesian box points, six original times and FULL original,
separately varied spatial/time/energy-quadrature ladders. The full-vector Euclidean
maximum and fixed-time L2=sqrt(64*mean(|R|^2)) are worst over six times at h=.005,
time step=.0025. All eight actual scientific acceptance CLIs exit1 for both momentum
gates. Other original sampled gates pass, not a continuum certificate.

## Fixed original contract

nu=.01, t[.25,.75], physical R3/evaluation box[-2,2]^3, smooth compact u AND p in
r<2 and |z|<2, original two-parameter divergence-free curl force a,c in[0,10], E0=1,
all original bounds/energy/core/support/divergence and momentum max/L2<=.001 gates.
Pressure coefficient bounds are checked on saved candidates. No amplitude collapse.

## Method and completed fits

40x60 whole-cylinder Gauss quadrature,17Gauss time nodes plus both endpoints (each
endpoint weight.05), exact pressure-gradient Gram and Cholesky solver coordinates.
Unconstrained least-squares diagnostics first quantify pressure-only capacity; they are
TRAINING-only and are not accepted pressure-direction candidates.

Final constrained objectives are quadratic and the sampled feasible sets are convex:
linear radial/signed-axial pressure-gradient restrictions and Euclidean residual caps.
3078core training points: R[.03,.22], |Z|[.02,.22],19times; positive pressure-gradient
margin1e-6. A distinct27228point training pool includes a deterministic spatial/time
mesh and3072random cylinder-space-time points, seed9175450. Two cumulative active-set
exchanges start from384worst and64spread pool indices, then add the next384worst.

M3 cap=max(.97*parent poolmax,1.005*azimuthal poolmax)=.03128900829;
Q2 cap=max(.90*parent poolmax,1.005*azimuthal poolmax)=.03615944229.
M3 completed5+4=9iterations; Q2 completed9+1=10. Both stages report numerical success;
this is not a continuous minimax or global certificate. Each exchange selects its best
feasible TRAINING iterate (1e-9 normalized active feasibility), then checks the entire
pool. The old generic `selection` label says final constrained iterate; this paragraph
specifies the actual selection rule. Both final entire-pool caps are met within floating
tolerance. Training MSE falls by13.64%/16.23%, distinct from independent L2.

## Structural invariance and independent pressure audit

Identical velocity coefficients in the same basis imply identical u, vorticity, shear,
energy, scaled velocity profiles and exact solutions of x'=u(x,t). We did NOT reintegrate
particles or claim new structural amplification. Existing structural weaknesses remain.
The force is also identical, while p is a new function.

Old242point/six-time and NEW800point/seventeen-time core checks, seed9175493, retain all
sampled inward radial/positive swirl/bipolar outflow/inward radial pressure/axial pressure
toward midplane signs. The fresh grid uses R[.035,.215], |Z|[.025,.215],15random interior
times plus endpoints. Not full-support/continuous-time certification.1025fresh physical
midplane probes give velocity and shear differences exactly0. No source amplitude target
or complete source-field correspondence is asserted.

## Independent grid and pressure-only limitation

A new54x115axis-inclusive full/edge-refined grid at13times gives:
M max.03776085243 -> M3.03777948640 (+0.049%);
Q max.03879113343 -> Q2.03852793410 (-0.679%).
All child grid peaks were independently Cartesian-FD checked, vector differences below
4.40e-8. The training pool misses higher independent peaks; no continuum upper bound.

For fixed axisymmetric u/f and every axisymmetric pressure, R_theta is invariant since
(grad p)_theta=0. On this independent grid its max is.02811182640 for M and.03669764936
for Q; these floating samples show that PRESSURE-ONLY repair of these frozen velocities
cannot plausibly attain1e-3. This is not an interval certificate, not a new theorem about
all candidates, and not a reason to stop joint velocity-pressure construction. All weak
residual moments paired with divergence-free test fields are likewise pressure-invariant
under compact support; ST053's velocity imbalance is not repaired by this experiment.

## All failed/diagnostic attempts retained

The first cap1.0 M/Q solves and unconstrained outputs are training diagnostics only.
Q completed5+1steps and saved its field, then its summary writer failed on a duplicate
`iterations` keyword. Metadata was repaired from unchanged saved field/stages; traceback
retained. Its surviving callback history has5records, not6. Later source fixes the writer.

A tighter M2 .90cap generated new pool peaks and worse training loss, then was explicitly
stopped after65saved callbacks, without claiming convergence or a final accepted field.
Checkpoint/stage/logs are retained. M3 intermediate .97cap and final M3/Q2 selection were
registered BEFORE any independent validation; amendments and freeze receipt are archived.

An audit-only bug was caught after validation: multidimensional `jets` output was indexed
along time instead of its final component axis in two midplane summary fields. Flattening
before indexing fixed it. Old audits are retained; same frozen fields were re-audited and
a regression added. No candidate, fitting, full residual or core pressure-grid value changed.

## Replay and checks

From repository root, without fitting:

```bash
python -m pip install numpy scipy sympy pytest
python -m pytest -q -W error experiments/root_st054/test_pressure.py
python experiments/root_st054/replay_st054.py --id ST054-M3 --out outputs/ST054-M3 --seed 9175491 --validate --structure
```

Use ST054-Q2 for the other field. Scientific CLI currently exits1 for both momentum gates;
it refuses a nonempty output directory. Eight final local tests actually passed13.66s,
warnings-as-errors; compileall succeeded. Includes full Cartesian derivatives, pressure
Gram/quadratic identity, azimuthal invariance, energy/support, both frozen recipes and
mutations, four local symbolic identities and the multidimensional-audit regression.
Full inherited suite and Lean were NOT run. Actual cloud execution is recorded in the PR.

GitHub stores readable source, exact binary pressure increments, numerical references,
result summary and read-only CI. Parents are rebuilt through ancestor recipes if raw
files are absent. The user archive contains SHA-pinned raw parents and children, actual
registrations, all attempted fits including failures, eight full validation reports,
audits/grids, freeze and test receipts. Not every older historical run is re-uploaded.
Local raw parent metadata differs from original M/Q bytes; reconstruction is numerical
field equivalence, not original raw-file identity. The ten reference u/p values and all
pressure increment bytes are checked; NumPy loads with allow_pickle=False.

Raw child SHA256:
M3 2b2b571986297b52bc884a2ee4808a7ade1ade86f38f938604f2a64118882df3
Q2 d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5

No prior source, default, physical configuration, threshold or agent schedule changed.
No external numerical solver, source annular pulses, matched heat exterior, continuum
residual certificate or blow-up proof. pde_validated/source_correspondence_verified false.
