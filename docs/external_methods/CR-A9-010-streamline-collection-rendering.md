# CR-A9-010 — collection-based 3D streamline rendering

## Current blocker

The final deliverable is a callable/savable 3D time-dependent velocity field, but a useful visual inspection commonly needs roughly 200 thin streamlines. Rendering every tiny streamline segment as an independent graphics object creates tens of thousands of artists and can make MATLAB/Python figures slow to open or rotate. Existing constrained work owns seed planning, streamline integration, morphology diagnostics, and arclength resampling; this increment owns only the final rendering-data packing layer.

## External result screened

- source repository: `matplotlib/matplotlib`
- screened commit: `cfa473e42fa41572d9acbd8d6b7be6f1489dee46`
- relevant public API: `mpl_toolkits.mplot3d.art3d.Line3DCollection`
- relevant source file: `lib/mpl_toolkits/mplot3d/art3d.py`
- license: Matplotlib License (permissive Python-style license)
- classification: **directly reusable public API**

At the screened commit, `Line3DCollection` is a collection type derived from `LineCollection` and accepts a collection of 3D lines. This repository does not copy Matplotlib rendering/projection code. It only packs already-computed streamlines into the segment/scalar arrays consumed by the public collection API, then optionally instantiates one `Line3DCollection`.

## Minimal migration scope

`constrained_streamline_collection.py`:

1. accepts only already-computed Cartesian streamline polylines plus optional point scalars;
2. converts all line segments into one `(n_segments,2,3)` array and midpoint scalar array for a single Matplotlib collection artist;
3. produces NaN-separated duplicated `(2,N)` x/y/z/CData arrays that MATLAB can draw with one `surface(...,'FaceColor','none','EdgeColor','interp')` object;
4. uses height as the default scalar, while caller-provided scalars can carry frozen-time speed or another visualization-only quantity;
5. changes no streamline geometry, velocity, pressure, forcing, optimizer, residual, threshold, or evidence state.

The MATLAB packing convention is a small autonomous interoperability design, not copied from Matplotlib and not a claim about OpenAI's visualization implementation.

## Truth boundary

Collection packing is visualization performance plumbing only. Fewer graphics objects, smooth interaction, attractive color gradients, or resemblance to a public image do not establish physical support, divergence, Navier–Stokes residual acceptance, OpenAI hidden-field identity, paper exactness, singularity, or blow-up. The pack therefore keeps all scientific/readiness claims false.
