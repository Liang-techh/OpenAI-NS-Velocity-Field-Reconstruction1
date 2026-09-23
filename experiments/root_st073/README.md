# ST073 — full-momentum local recurrence, not global NS admission

Task #1255. Manual bounded continuation; paused schedules remain paused. This round implements a local radial recurrence for ALL three momentum equations, not only the leading similarity equations. The reported finite local samples pass 0.001, but the original full-space fixed-force/energy/support benchmark and scale-recursion capability remain UNMET. No default, main, viewer or schedule was changed.

## What actually changed

In source-viscosity-one coordinates, write s=x^2+y^2 and u=(x A-y B,y A+x B,C). Incompressibility gives A_n=-partial_z C_n/[2(n+1)]. The full angular/axial equations recursively generate B_(n+1),C_(n+1), then the radial equation generates p_(n+1). Time derivatives, ALL axial viscosity, radial inertia, convection and pressure enter before generating the next order. Pressure is not held to the leading centrifugal relation. A two-variable Lagrange expansion of the implicit q equation supplies consistent z/time jets. Physical nu=.01 scaling is explicit.

This is an autonomous local Taylor construction, not an implementation or proof of the source paper's global correction scheme. Finite radial and Lagrange truncations have no uniform infinite-series or interval remainder certificate.

## Two separately frozen axis-data controls

Local domain: X<=1/64, |eta|<=.5, h=.005, tau in [.5/64,.5], local f=0. No whole-space force is redefined and no global/per-frame energy normalization is applied.

F uses the smooth ST068 axis traces unchanged: swirl factor1, axial profile4eta+.02, pressure -1+eta^2/2. It is NOT a refit of sharp ST072-N. Calibration order2/4/6/8/10 yields complete maxima25.1044/.0317288/.000218249/5.44585e-7/3.77920e-10 on the same declared points. The first order meeting the preregistered1e-5 calibration screen is8, confirmed at10. Order12 was interrupted after two scales and excluded. F's radial pressure force points outward in the tested core despite small momentum residual.

After that failure, V was separately registered with swirl trials2,3,4. Amplitude4 and radialorder8 were the first passing fixed calibration choices, confirmed atorder10. Axis slope,bias,pressure remain unchanged. This amplitude change is explicit, not energy preservation or an amplitude-collapse trick. V is compared with a newly generated leading-only field with the SAME V axis traces. F and V use different fresh holdouts.

F frozen2026-09-23T19:44:57.514363+00:00; V frozen2026-09-23T19:48:32.992875+00:00. No post-holdout coefficient tuning.

## Independent full physical local results

|Field|Seed|Full-vector sampled maximum|Divergence maximum|
|---|---:|---:|---:|
|ST073-F|9237391|4.4891615248e-7|5.55e-17|
|ST073-F|9237392|5.0685733767e-7|5.55e-17|
|ST073-V|9237395|4.8692457014e-7|5.55e-17|
|ST073-V|9237396|4.7230327801e-7|5.55e-17|

Each row comprises4096 CLUSTERED space-time points, not IID whole-space samples or4096pointsper time. V uses32 independent axial positions x16radialpoints x8logtimes (0,3,6 plus5newinterior times); F uses64x8x8. V's five finite flow/pressure-direction fractions are all1; bipolar and axial-pressure checks exclude |eta|<=.02. These are sampled properties, not continuum or global claims.

Matched-domain and matched-V-axis quadrature,18radial x32axial nodes:

|k|Leading-only control full local L2|ST073-V full local L2|V quadrature-node maximum|
|---|---:|---:|---:|
|0|.00953151461199|1.4974424946e-12|1.3121516468e-9|
|3|.04511377600175|6.2605930210e-12|2.5915428936e-8|
|6|.21353240931223|2.6189318725e-11|5.1208435632e-7|

Radial and axial quadrature resolutions varied separately. Tiny L2 partly reflects the tiny integration volume: at k6 volume~8.90e-8 and localenergy~1.01e-7. Always report maxima as well. These values must NOT be compared to old whole-space norms or represented as E0=1. Same-axis local effective volumes change by only about-.25% (energy) and-.1% (enstrophy), but there is no global effective-volume acceptance.

Independent Cartesian finite differences recompute all terms at six physical points over k=.4,2.7,5.5, separately refining space and time in three levels. V's independent residual is at most~1.02e-6 over these levels; finest vector discrepancies<=1.11e-7. Roundoff plateaus are retained. This does not independently certify1e-12PDEerror or a continuous supremum. Float64 versus extended arithmetic,120vs160axis sums and larger jets were separately compared.

## Remaining interface and scientific failures

V's old leading-background midplane shear diagnostic v at X=1/64 is only~.0093-.0099, below the old wave-cone prerequisite2. That leading formula is not a universal NS condition and cannot be transplanted unchanged onto a full corrected background, but the old interface is plainly NOT preserved. New pressure, radial velocities and boundary jets require a fresh annular/exterior matching calculation. Old ST069 moment arrays must not be reused.

No new five-moment restoration, heat-exterior match, admissible stress cone, real non-axisymmetric waves, global axial localization, finite total energy, original fixed-force validation, material winding or dynamically generated time-scale recurrence. No fake full-field animation by zero-extending this local expression. All global/PDE/sourceidentity/blowup flags remainfalse.

Eight additional ST072 smoothing/pressure controls were actually run and retained. Wider axial transitions reduce sampled viscous peaks, but the tested old-radius monotone-pressure bounds still fail. No such control was promoted.

## Actual tests, replay and remote availability

Primary18focusedtests passed6.40s; independent copied bundle18passed6.39s,warnings-as-errors. Compile and166-file pre-receipt inventory verification passed. Packaged V k6/18x18 report actually reran in7.09s and matched the ENTIRE primary report dictionary exactly. --resume verified/skipped it in1.33s. Extra holdouts, derivative ladders and expanded quadratures ran in the primary directory and were copied identically, not all rerun in the clean smoke. No MATLAB,cloudCI,Lean or fullhistoricalsuite execution.

The actual evaluator, coordinate dependency and BOTH generative model JSONs are committed. Their returned Git blob identities match the executed local bytes. These four files also ran in a separate minimal directory. Complete audit dependencies, coefficient snapshots, all reports, negative controls and logs are in the conversation archive NS_ST073_Full_Local_Recurrence.zip, not claimed all imported to GitHub. No PR or merge.

Minimal evaluator from this directory:
```python
from full_radial import FullRadialField
f = FullRadialField.load('data/ST073-V.json')
tau = 0.5 / 64  # remaining time, NOT physical t
xyz = f.from_similarity([.005, .01], [-.2, .3], tau)
print(f.evaluate(xyz, tau)['residual'])
```

From the COMPLETE offline bundle:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --candidate ST073-V --k 6 --order 18 --out outputs/V_k6.json
```
The last command returns0 ONLY for its finite LOCAL check. Its report keeps pde_validated=false, global_field_ready=false and original global benchmark=NOT CONSTRUCTED. Add--resume to verify/skip completed data without fitting.

V JSON SHA2569066b03217fdfa1a83679afdded4a60d55bf1775867b8bc888360af9cc47b574.
F JSON SHA256653eaf03d9ce29d530eb26ec13c817ccbfce7a3f973e03de9444aaed58cce176.
Full evaluator SHA256d45518365431b03f42161ef0d202be5dcc6803a856de21b54a24d2d602db8643.

Next: preserve small full local momentum while adjusting actual axis/boundary data and rebuilding the exterior feedback. Do not revert to trading the governing equations against a soft seam penalty or claim global success from this small local domain.
