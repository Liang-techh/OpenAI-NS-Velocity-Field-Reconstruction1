# CR-A9-037 — OpenAI public visual target contract

## Source screened

- Publisher: OpenAI
- Public source: `https://openai.com/index/navier-stokes-solution/`
- Publication date: 2026-09-08
- Retrieval/audit date: 2026-09-18
- Source type: public article text, displayed figure, and figure caption
- Repository / commit: not applicable; this increment consumes no external source code
- License treatment: no redistribution license is asserted and no OpenAI image/source implementation is copied. The repository records only short factual/qualitative observations needed to govern candidate review.
- Classification: **direct public-observable extraction / governance only**

## Public observations retained

The contract records only statements exposed by the public article/caption: the field is described as a spinning vortex/swirl; trajectories spiral inward and show axial stretching; the vortex becomes increasingly elongated; the central region shrinks while the flow speeds up; orange denotes faster angular rotation while teal denotes slower rotation; and circulating speed depends on radius.

These are qualitative observations. The source does **not** publish numerical velocity samples, a physical coordinate calibration for the displayed figure, a numerical frame time, camera pose/projection, streamline seed locations, hidden construction parameters, or quantitative visual-acceptance thresholds. Those remain explicitly unobserved in the contract.

## Local implementation

`constrained_public_visual_target.py` freezes those observations behind stable IDs and provides `audit_visual_evidence_receipt(...)`. Candidate diagnostics may mark an observable as measured, inconclusive, or not measured and may carry their own candidate-side numeric measurements. They may not inject numerical target values/thresholds, hidden frame time, camera registration, recovered field/parameters, or promote visualization/PDE/OpenAI-identity truth states.

The audit returns coverage only. It deliberately defines no visual pass threshold and cannot set `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved` true.

## Direct contribution to final `[u,v,w]`

The repository already contains or has open work for streamlines, pathlines, vorticity/Q/`lambda_ci`, helicity, core traces, temporal trends, and multiple export formats. Without one frozen public-observable target, those diagnostics can accumulate without a clear truth-bounded criterion for what each is meant to support. This contract gives future candidate review a single source-bound checklist: candidate-side evidence can be attached to the public observations without pretending that the public image reveals hidden numerical data.

This does not modify any velocity, pressure, forcing, support, energy normalization, optimizer, validation sample, residual operator, or scientific threshold.
