# CR-A9-021 external method record — fixed mask connectivity

## Purpose

Provide one deterministic, truth-bounded descriptor for whether an explicitly
presegmented public/candidate projection is one continuous visible structure or
is split into detached layers/blobs.  This is a visualization diagnostic only;
it does not alter or validate `velocity(x,y,z,t)->[u,v,w]`.

## External source screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.ndimage.label`
- source file: `scipy/ndimage/_measurements.py`
- license: BSD-3-Clause terms in SciPy `LICENSE.txt`
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; this repository already requires `scipy>=1.10,<2`

SciPy's public `label` API supplies connected-component labeling.  Repository
code fixes connectivity to a four-neighbor cross and owns all input checking,
component census, axial-gap/run descriptors, physical pixel scaling,
provenance, and truth-state metadata.

## Migrated scope and differences

Only the public connected-component labeling call is reused.  The repository
wrapper deliberately does **not** expose connectivity as a fit parameter and
does not perform segmentation, threshold selection, component cleanup, gap
bridging, hole filling, translation/rotation/reflection/scaling, camera fitting,
registration, hidden-time inference, or hidden-velocity recovery.

The wrapper additionally reports:

- component count and sorted component sizes;
- largest-component and detached-foreground fractions;
- number of contiguous axial row segments;
- total and largest blank axial gaps;
- fraction of active rows with multiple transverse foreground runs;
- maximum transverse run count;
- image-frame contact.

Diagonal-only contact is intentionally disconnected under the frozen
four-neighbor convention.  This policy prevents topology from being changed
post hoc to improve a visual score.

## Truth boundary

A small component count or low detached fraction is only a public-observable
morphology measurement.  It is not evidence that OpenAI's hidden velocity was
recovered, does not make the field `visualization_ready`, and cannot promote
PDE validity, paper exactness, field identity, singularity, or blow-up claims.
No visual pass threshold is defined here.
