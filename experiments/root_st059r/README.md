# ST059R — tested interruption recovery and bulk-volume continuation

User task: first repair interruption recovery, then continue structure/residual/volume research. This is a NEW refit from the complete ST058-IP, not a recovered copy of historical ST059-B/E. The prior ST059 completion note exists, but its ZIP, patch and raw arrays were not found in the current runtime or Library. No numerical comparison to unavailable E is claimed.

## Recovery actually tested

The committed checkpoint.py and bulk_step.py implement complete per-step generations containing full candidate JSON, coefficient increments, optimizer damping/history, actual reduced-direction matrix and immutable parent/source bindings. Generation payloads and hashes are synced before a directory rename commits them. Incomplete temporary directories and stale latest pointers are handled; corrupt committed data or changed immutable inputs are rejected, not silently reused.

An actual two-step optimization was run uninterrupted and compared with one step followed by a new process resuming the checkpoint: maximum delta difference0.0, full candidate bytes identical and damping identical. The surrounding tool driver timed out during that test; the saved generation was subsequently recovered. Unit tests separately inject pre-commit interruption and corruption. Production was also actually stopped after step4, exported as a downloadable checkpoint ZIP, then resumed in a new process and completed18steps. This is tested local recovery, not a guarantee against host termination or temporary-storage loss.

Validation tasks have per-report completion receipts and input/hash checks. A phase-aware entry refuses further fitting once holdouts are frozen. The portable full ZIP, not a hash-only progress note, contains the actual arrays and state.

## Actual numerical continuation

Same2594stored-parameter compact axisymmetric family,578reduced coupled velocity-pressure directions. Initial spatial distribution remains free, E(.25)=1 normalization and its derivative remain exact in the representation; original prescribed force stays fixed. Originalnu=.01,time[.25,.75],compact smoothu ANDp in r<2/|z|<2,rawbounds and both0.001momentum gates unchanged. No residual-defined arbitrary force, new fit amplitude collapse or default promotion.

New mapped radial32Gauss quadrature and axial14/12/36/12/14orders on five signed slabs,13Gauss times+endpoints,31-time denser peak pool. Worst-time reweighted integrated full residual, harmonic compatibility, signed-slab residual-energy screens, active peaks and finite core/shear/morphology safeguards feed local quadratic steps with actual nonlinear screening. Streamed128spatial-point Jacobians limit memory. The slab screen explicitly allows ratio1.0005, not strict monotonicity. Eighteen bounded steps,176.96recordedseconds; no optimizer-convergence/global-minimax claim.

All coefficients frozen2026-09-20T18:10:59.864792+00:00 BEFORE new seeds9206091/9206092. No retuning after validation.

## Six completed original independent validations

Each4096fresh Cartesianpoints,sixoriginaltimes,FULL separately varied original spatial/time/energy-quadrature ladders. Fullvector Euclideanmaximum and spatialvolumeL2=sqrt(64mean|R|²),worsttime at h=.005/timeh=.0025.

|Seed|Field|Full sampled max|Spatial volume L2|
|---|---|---:|---:|
|9206091|ST058-IP|.02848094025174967|.042904470814721055|
|9206091|ST058-S|.0343863058246871|.04256038617220756|
|9206091|ST059R-V|.027624409988372185|.039932038980469244|
|9206092|ST058-IP|.028501961957705654|.04256673334946904|
|9206092|ST058-S|.03382246732007187|.04231146933233183|
|9206092|ST059R-V|.02750508849086926|.039596274927955255|

V paired L2 improves6.928%/6.978%, sampledmaximum3.007%/3.498% versusIP. All6original scientific CLIs exited1 for momentum_max and momentum_L2; other original gates pass on these samples. These are not continuous-domain bounds.

Independent81x153grid/21times: maxIP.0300820390 toV.0300691745, only~.043% improvement. Worst point remains initialt=.25,r=.2838,z=1.93835. Original CartesianFD crosschecks at all42gridpeaks agree within1.503e-7. No peak breakthrough. Further randompeak spatial refinement gives V .02763451/.02751134 at h=.00125; the small sampled advantage persists.

## Structure and volume

Fresh800corepoints/17times: all five velocity/pressure direction fractions1. Fresh1025midplane shear probes: minimum signedretention.9994624,above.995. Upwardbias remains positive but minimumweakens. Same-grid final scaledprofiledrift19.0340% to18.8377%,small improvement,not a fixed similarity profile. These probes cover a declared small scaledcore,not the whole support.

Initial relativevelocityL2change1.55194%; independent64x96energy E0=.9999999999994953. Energy participationvolume decreases~.64%-.91%, enstrophy participationvolume~2.13%-2.57% across3times: modest concentration change,not exactvolume invariance or observedcollapse. Original declaredsupport cylinder16pi remains unchanged.

At t=.75, independently mapped squaredresidualintegrals fall14.54% in bulk|z|<=1.75,15.48% in adjacentcollars,9.13% in outercollars. One signed outer slab at t=.5 increases0.00365%,retained; no all-region/all-time monotonicity. Combineddegree2-8 necessaryL2lower-bound ESTIMATES at .25/.5/.75 improve(.00549608,.00419555,.00533516) to(.00537529,.00392092,.00506504),still>.001. Exact identities plus floatingvelocityquadrature are not interval-certified bounds.

## Checks, corrections and publication boundary

30focused new/inherited tests actually passed14.71s,warnings-as-errors; clean deliverycopy30passed16.94s. Independent full operator/Jacobian/normalization, selected inherited symbolic identities, checkpoint corruption, frozen-phase and final parameter checks included. Full6validation reports and both complete independent audit sets finished locally. No cloudCI/nativeMATLAB/Lean/fullhistorical suite run.

A supplementary diagnostic import failed before calculation and was fixed. One regional metadata sentence incorrectly called all diagnostic times unused by training; .25/.5/.75 are shared and .4173/.6327 additional. Corrected metadata and full audits rerun with ALL numerical values unchanged. Original failedsource/logs and oldreports retained. No field/threshold changes.

Committed optimizer blob1e2f27a89a3adc530b41cdb2980279ee34fc3f53 matches the executed local file. Full complete dependency stack, rawparent/S/V, exactmap, generations0-18, allcandidatevectors/histories, reports and actual recovery receipts are provided in the conversation ZIP. Their presence is not inferred merely from this README. Fullraw arrays and checkpointmap are NOT yet uploaded to this branch. No PR/merge/default/viewer/publication or other-agent/schedule changes.

RawV SHA256:9d76d3f6236e3bef78dfdbcbed30f69bb84d6190497e1b33092f96ddf27e412c.
RawparentIP SHA256:5ea5d117d7330f34b0b32f7dd67828c72efd0e03ac007bf9ca3115ce9f21123f.

From the complete delivery root:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
python experiments/root_st059r/resume_session.py status
python experiments/root_st059r/replay_st059r.py --id ST059R-V --out outputs/V --seed 9206091 --validate --structure
```
The last command currently exits1 for failedmomentumgates. Replay neveroptimizes. Add --resume to verify/skip completed outputs. Further optimization after the freeze requires a NEW experiment and NEW holdouts.

OriginalfullNS1e-3 target UNMET. pde_validated/source_correspondence_verified/paper_exact/blowup_proved=false.
