# CR-A9-023 — truth-bounded shared-frame multi-candidate visual comparison

## Goal

Close one visualization-candidate selection seam in the active Eq45 path: compare two or more already-frozen `velocity(points,t)->[u,v,w]` candidates on **exactly the same physical meridional frame, times, sampling grid, and speed normalization**. This is needed now that quartic, compact-C2, and bounded blend candidates are being screened concurrently. Per-panel autoscaling can otherwise make a weaker or decaying candidate appear equally strong.

The increment is candidate-agnostic and changes no velocity value.

## External screening / migration

Selected source:

- repository: `matplotlib/matplotlib`
- screened commit: `cfa473e42fa41572d9acbd8d6b7be6f1489dee46`
- source file: `lib/matplotlib/colors.py`
- public API: `matplotlib.colors.Normalize`
- license: Matplotlib License (permissive Python-style license)
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; Matplotlib is already the repository's `dev` visualization dependency

`Normalize(vmin, vmax)` publicly maps one declared scalar interval to `[0,1]`. This module uses exactly one instance with `vmin=0` and `vmax=max(|u|)` over **all candidates and all declared times**. It does not call `autoscale` separately for each panel.

The core sampler itself requires only NumPy. Matplotlib is imported lazily only by the optional PNG renderer.

## Repository-local behavior

`sample_shared_frame_candidates(...)`:

- requires at least two uniquely named frozen velocity callables plus non-empty provenance;
- builds one caller-declared `y=constant` meridional x-z frame and one strictly increasing time list;
- calls every candidate on the same point coordinates and same times;
- retains full Cartesian `[u,v,w]`, total speed, projected poloidal speed `sqrt(u^2+w^2)`, and `|v|` swirl magnitude;
- computes one global speed maximum shared across every candidate/time panel;
- fail-closes if any panel is numerically inactive, malformed, or nonfinite;
- performs no camera fit, registration, alignment, coefficient fit, threshold search, hidden-time alignment, or candidate selection.

`render_shared_frame_candidate_comparison(...)` optionally writes a side-by-side PNG with one shared Matplotlib `Normalize`. The white line overlay uses projected meridional `(u,w)` only and is explicitly labelled **not a true 3-D streamline** when swirl is present.

## Truth boundary

This is a fair-rendering/comparison primitive only. It does not:

- compare to an OpenAI image or infer hidden image coordinates;
- segment public imagery or choose a visual threshold;
- recover hidden velocity, pressure, forcing, coefficients, or frame times;
- declare any candidate visually accepted;
- change domain, `nu`, support, forcing family, energy normalization, optimization/validation split, residual definition, or acceptance thresholds;
- convert visual similarity into PDE validity.

The anti-zero activity floor is only a fail-closed rendering guard. It is not a scientific acceptance threshold and cannot make `u -> 0` succeed.

## Direct contribution to final velocity delivery

The active checkpoint explicitly calls for quartic, compact-C2, and bounded blend candidates to be compared through the same public `[u,v,w]` interface before choosing a visualization candidate. This module supplies the missing fair-view layer: one identical physical frame, one time list, and one global amplitude scale across candidates.

It complements, rather than replaces, whole-domain vorticity fingerprints, true-3D streamline diagnostics, public-observable mask metrics, and independent PDE validation. A visually preferable candidate can be identified for further review without being relabelled PDE-valid or as the hidden OpenAI field.
