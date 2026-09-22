# ST065 — fitted scale-native profiles, not accepted NS recursion

Task #1091. Bounded numerical continuation of the FAILED ST064 imposed-scale scaffold. Actual profile fitting and independent audits are complete. No accepted/default candidate, publication, other-agent task or schedule is changed. Full physical 0.001 momentum gates remain failed.

## Main result: preservation-control ST065-C

In 800 fresh scale-aligned core probes, mean actual radial inflow divided by the prescribed shrinking-window speed increases from **0.05702819 to 0.96951650**; the child range is approximately 0.81095–1.04300. This changes the physical velocity, not just the displayed window.

The SAME twelve particles have mean signed cumulative winding **0.06364325 -> 0.32933506 turns** by k=6. The child peaks at 0.35377009 turns at k=3, then has small reverse increments. All twelve have left the shrinking observation core by k=4. A separate, special three-particle z=0 control makes 1.17645007 mean turns and remains inside; it is NOT representative of the whole core. Scale-aligned re-seeding uses different material particles and is reported separately. A shrinking Eulerian core need not retain every parcel; exit counts limit this winding claim rather than proving all Eulerian recursion impossible.

Core aspect at k=0/3/6 is 2.11625014 / 2.13836801 / 2.16071705, versus ST064 1.73146643 / 1.74956276 / 1.76784822. Cross-scale aspect growth is STILL only 2.10%, imposed by h=.005. Improved initial geometry is a fitted profile change, not a recovered recursion exponent. Characteristic diffusion-time turns 1.20019 -> 1.22541 are not actual particle turns completed in one interval.

## Physical residuals remain huge; peaks worsen

Same imposed time scales, viscosity and prescribed force. Full-domain integrals include R=-f outside shrinking velocity/pressure support. Values below use independent 96-order quadrature and a 91x151 reference-space peak grid.

| k | Physical time | ST064 L2 | ST065-C L2 | ST064 sampled max | ST065-C sampled max |
|---|---:|---:|---:|---:|---:|
|0|.5|6.34112235119503|4.012768217277348|2.1264706465994547|2.1753725524408822|
|3|.9375|29.738426015629027|19.263710182169003|45.87670674756071|51.15135848706049|
|6|.9921875|143.1487507998976|92.84198336347596|1039.991318628079|1168.9129712012757|

L2 falls 36.72% / 35.22% / 35.14%, but peaks INCREASE 2.30% / 11.50% / 12.40%. No sampled or continuous NS acceptance. This new ansatz/horizon cannot be compared with the old ST063 0.03 residual as though only coefficients changed. Normalized training residuals never replace physical acceptance values.

## Actual construction and phases

Static F/G/P have 108 coefficients each, 324 total. The old 9x12 bump-polynomial SPATIAL family and original ST064 physical operator are reused; the eight-term old ST063 temporal polynomial is not extrapolated. T=1, tau0=.5, k in[0,6], h=.005, nu=.01. ar=q^.5, az=q^(.5-h), q=2^-k. All scale-chain derivatives are included. Velocity and pressure are smoothly compact, with the unchanged original force in its unshrunk support. Energy is normalized only once at the NEW initial time t=.5. No arbitrary force=R or per-frame normalization.

Phase A fits R/E/H with autonomous core target weights0/1/10, A=-(1+h)/(2tau0), B=4, C=-2Az. Use 40x40 spatial quadrature, seven Gauss log-times, 96-order kinetic/pressure-gradient variable conditioning and analytic full-momentum/normalization adjoints. Residuals are multiplied by q^1.5 only inside the optimization objective. L-BFGS stopping counts R438/E357/H447 (2.63/2.71/3.37seconds), not global optimality. All frozen before fresh log/point seeds9226591/9226592.

Free E reduces L2 more, but initial energy participation volume falls about70% and enstrophy participation volume about91%; axial pressure direction reverses and bias nearly disappears. It is retained as a FAILED preservation control, not a preferred field.

After those failures, a NEW phase-B registration starts V again from ST064 with volume, pressure-direction and small-bias penalties. It reaches the2400iteration cap (23.88seconds), not convergence. Then a registered108-variable minimum-pressure-gradient repair fixes finite linear sign guards while preserving V's complete F/G arrays exactly; SLSQP reports success in2iterations (.0165seconds). Result C is frozen BEFORE new log/point seeds9226595/9226596. No old-holdout retuning is represented as new independent evidence. C means preservation-control, not NS-compatible.

## Independent structure and volume counterevidence

All five core flow/pressure fractions are1 for C on the fresh finite window: r=.025..15 and z=-.25...25 in reference coordinates; |z|<=.03 excluded for bipolar/axial-pressure tests. Positive midplane bias is retained. These do not certify the full support or exact source correspondence.

Initial energy=.9999999999783935; at k6 it is.1321759185387313. Same-time energy participation-volume ratio at k0/3/6 is.99020872/.99041123/.99061696. Enstrophy ratio is **1.02062404**/1.00573677/.99100626. The initial independent +2.0624% slightly EXCEEDS the declared +2% goal. No rounding into a pass; the soft training guard did not produce an exact independent all-time guarantee. No amplitude collapse, but concentration is not literally unchanged.

## Audits actually run

Phase A: four models(ST064/R/E/H) times12scales. Phase B: two models(ST064/C) times12freshscales. Total72model-scale reports, each with40/64/96full-domain quadrature, separate Cartesian spatial and temporal FD ladders, peak grid and fresh core data. Seven integer scales and five independent interior log-times per phase. These are extended-window physical audits, NOT72old-window4096-point ST063 validations. All physical momentum gates fail. No global peak/interval certificate.

True particle trajectories use DOP853, two tolerances, explicit similarity-coordinate drift, a twelve-seed ensemble and separately labeled midplane/re-seeded controls. A further independent physical-x/y/z/time integration of three particles through k2 agrees with transformed-coordinate ODEs within1.63e-13 position and2.35e-12 angle. These are numerical operator/trajectory consistency checks, not small NS residuals.

## Algebraic diagnostic and research direction

SymPy verifies zero residual for a local affine field u=(-ax-Omega*y,Omega*x-ay,2az), a=alpha/tau, Omega=omega0*tau^(-2alpha), with its matching quadratic pressure and ZEROforce. Its axial pressure curvature is negative for positivealpha, and the velocity has infinite energy/noncompact support. This is NOT a project candidate. It only demonstrates why matching radial strain alone need not give the desired axial-pressure balance; it does not exclude general nonlinear/asymmetric profiles.

Primary context read: Duraiswami arXiv:2609.17642 (2-D similarity profiles, porous walls, axial throughflow and resolution limits) and Maekawa–Miura–Prange arXiv:1807.10341 (Burgers vortices with imposed linear strain). Their domains/backgrounds differ; no ancillary solver was run or source success imported. Next work requires coupled nonlinear axial throughflow, pressure and transition/exterior matching, not simply higher swirl or more displayed lines.

## Persistence and actual tests

Three executed optimizer files are committed at **b05037ae9794f291cbe8cb8eb3be09d627cb1609**. Returned Gitblob IDs match local executed bytes exactly. Complete dependencies, all profile arrays, parameter snapshots/history, six-model records, all72audits, particle trajectories, figures, MATLAB data/UI and test evidence are supplied in the conversation archive; they are NOT all claimed already on this branch. No PR/merge/default promotion.

Main20focusedtests passed3.18s; final packaging rerun20passed2.87s; independent clean copy20passed3.32s. Clean C k6 replay ACTUALLY recomputed all40/64/96physical integrals, reference-core metrics and91x151gridpeak exactly, plus physical FD checks, returning scientificexit1. Resume verified/skipped the report, still exit1, no optimizer. The complete72-scale/particle audit was run in the primary directory and not fully rerun in that smoke test.

Periodic L-BFGS array snapshots are preserved but NOT the full optimizer state for bitwise recovery. Initial gradient calibration found a missing derivative term and it was fixed before any fit; full-operator parity tests passed before fitting. A packaging-only loader metadata revision labels the actual phase-B seeds/guards; its prior executed snapshot is preserved. No profile or physical operator changed after freezing. Plots use Python, not native MATLAB screenshots.

The new MATLAB viewer is included but **NOT natively executed**. Python MAT coefficient/reference roundtrip passed. No cloudCI, Lean or full historical suite run. Default MATLAB command is `start_scale` for C; `start_scale('ST064-S')` opens the old scaffold. A logarithmic slider and three coordinate views are included; instantaneous streamlines and actual material paths are separate modes.

Use COMPLETE offline delivery:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --candidate ST065-C --k 6 --out outputs/C_k6.json
```
Final scientific exit1 is expected; --resume verifies completed reports without fitting.

C profile SHA256: **8b4bf31eae7fb48455ab0dc196688b670521f5fa959be7f31d52dd4b27cba54c**.

All PDE, dynamically generated recursion, source-correspondence and blow-up achievement flags remain false.
