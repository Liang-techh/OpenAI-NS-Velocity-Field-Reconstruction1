# Shape direction of the cone-compatible wave tangent

The cone-compatible wave candidate does not yet point in all requested
geometric directions. On a fixed physical cylinder its affine tangent
expands the enstrophy-weighted radial RMS and decreases the axial/radial
ratio. Enstrophy-weighted angular speed increases, while sampled maximum
swirl decreases. This is a reason not to advance this tangent as a
successful recursive step, despite its sampled compatibility pass.

`wave_tangent_observables.py` reconstructs the saved candidate and evaluates
the existing vorticity/enstrophy observables at k0 plus or minus 1e-6,
with k0=11.0003 and tau=0.5*2^(-k). Increasing k advances physical time
t=-tau. The cylinder is fixed at k0 and uses 7,776 points (12 angles).
Thus the sampling domain itself does not prescribe contraction.

| Observable | At k0 | Central derivative per k | Desired direction |
|---|---:|---:|---|
| Enstrophy radial RMS | 0.0011358189 | +0.00379472 | Decrease: fails |
| Axial/radial RMS ratio | 0.20859563 | -2.469372 | Increase: fails |
| Enstrophy-weighted angular speed | 728003.90 | +7.687531e9 | Magnitude increase: passes |

The maximum sampled swirl falls from 6578.47 to 6560.36 between the two
offset samples. Weighted spin and peak swirl are different diagnostics;
one positive spin statistic must not be presented as universal winding
amplification.

These are finite-grid, finite-difference observations of an affine tangent
extension, not a resolved NS trajectory or a globally identified vortex
core. The small offsets are used for local direction, not as accepted time
steps. Spatial/time-step convergence and persistence of the trends are
unproven. Complete momentum remains far above the acceptance threshold.

Reproduction from repository root:

```powershell
python experiments/root_st073/wave_tangent_observables.py
```

Next incorporate instantaneous shape-rate conditions into the constrained
tangent solve. At frozen velocity and vorticity, the derivatives of these
integral observables are affine in velocity time-derivative coefficients;
pressure coefficients contribute no velocity. Any new feasible candidate
still needs actual spatial momentum and evolving-geometry assessment.
