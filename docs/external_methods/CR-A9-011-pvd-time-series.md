# CR-A9-011 — ParaView PVD time-series index

## Current blocker

The constrained repository already has candidate-agnostic VTK snapshot export work, but a folder of separate snapshots does not preserve one tool-readable physical-time sequence by itself. That makes 3D time-evolution inspection more manual than necessary even after a frozen `velocity(points,time)->[u,v,w]` candidate exists.

## External result screened

- source repository: `pyvista/pyvista`
- screened commit: `5cb68204922cac06b18dd3587b92db05163a72be`
- relevant public behavior: `PVDReader` exposes collection time values and lets callers activate a time value; its PVD parser records `DataSet` attributes `timestep`, `part`, and `file`
- license: MIT
- classification: **suitable to reimplement minimally**

No PyVista/VTK source code is copied and no dependency is added. This repository writes only the small XML `Collection/DataSet` index needed to associate already-existing snapshot files with declared physical times. The local implementation is intentionally narrower than PyVista: it is a deterministic single-part writer, requires portable relative file references, and adds a separate SHA256/truth-boundary sidecar. It is not a general VTK reader/writer.

## Migration scope

`write_pvd_time_collection(...)` consumes existing visualization snapshot paths plus strictly increasing times and writes:

1. one `.pvd` collection with one `DataSet` entry per time; and
2. one `.pvd.meta.json` sidecar binding the PVD and every referenced snapshot by SHA256.

The sidecar keeps `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved` false. PVD generation changes no velocity values and does not inspect PDE residuals.

## Direct contribution to final `[u,v,w]`

Once a frozen candidate has been exported as a sequence of VTK-family snapshots, the same files can be opened as one physical-time series in a PVD-aware visualization tool instead of being selected manually frame by frame. This makes the time evolution of streamlines, vorticity and velocity geometry easier to inspect while retaining the original snapshot bytes unchanged.

## Remaining limitation

This increment starts from `main`, where the active Eq. 4.5 candidate and VTK exporter remain in open/integration work. The regression therefore indexes synthetic placeholder snapshot files and verifies the collection schema/hashes; it does not claim that the production Eq. 4.5 snapshots have been exported, visually matched, physically supported, or PDE validated.
