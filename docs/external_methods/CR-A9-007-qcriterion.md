# CR-A9-007 external method screen: VTK Q-criterion

## Classification

**Suitable to reimplement minimally.**

## Public source

- Repository: `Kitware/VTK`
- Screened commit: `723b76dba6baf6ccc7d16fe5b9b10d4971b4193d`
- Relevant file: `Filters/General/vtkGradientFilter.cxx`
- Public behavior: `vtkGradientFilter` can compute a Q-criterion scalar from a 3-component velocity gradient. Its implementation uses the incompressible form `Q = -0.5 * trace((grad u)^2)`.
- License: VTK's main source is BSD-style / BSD-3-Clause compatible; the screened source file carries the VTK project copyright/license regime.

## Migration scope and difference

No VTK source code or derivative machinery is copied. The repository only reimplements the mathematical Q definition in NumPy on a uniform Cartesian visualization grid, using centered second-order finite differences on interior points. This avoids adding the large VTK Python dependency while retaining the one vortex-identification quantity needed by the constrained visualization workflow.

The local diagnostic additionally reports Q-positive radial/axial quantiles, positive-Q volume fraction, and grid-resolution sensitivity. These are autonomous visualization diagnostics, not VTK behavior and not public facts about the OpenAI field.

## Why it lowers the current blocker

Vorticity magnitude alone cannot distinguish rotation-dominated vortex cores from strong shear/strain. The final callable `[u,v,w]` is intended to reproduce the public vortex morphology, so Q>0 gives an independent rotation-versus-strain core diagnostic that can be used alongside the existing vorticity and streamline fingerprints when tuning the Eq. (4.5) `Phi/F` profiles.

## Truth boundary

A stable Q-positive core is visualization morphology evidence only. It does not validate the Navier--Stokes momentum equation, divergence threshold, physical support, OpenAI hidden numerical data, paper exactness, singularity, or blow-up. No universal visual acceptance threshold is introduced in this increment.
