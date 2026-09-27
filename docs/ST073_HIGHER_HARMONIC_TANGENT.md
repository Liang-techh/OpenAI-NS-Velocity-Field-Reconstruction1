# Spatial enrichment of the second harmonic

Replacing only the mode-2 time/pressure tangent basis with spatial degree 3
reduces the independent-grid momentum L2 from 2,847,158.20 to 2,632,913.46
(about 7.5%). Degree 4 lowers training error further but slightly worsens
the independent grid. Degree 3 is therefore the useful enrichment to carry
into a subsequent jointly constrained construction.

| Mode-2 degree | Training volume L2 | Independent volume L2 |
|---|---:|---:|
| 2 | 2,868,781.90 | 2,847,158.20 |
| 3 | 2,642,709.87 | 2,632,913.46 |
| 4 | 2,586,547.30 | 2,635,302.48 |

The wave velocity, mode-0 and mode-1 controls remain fixed at the
cone-compatible candidate. The correction changes only the mode-2 velocity
time derivative and pressure. Its zero instantaneous velocity and zero
angular-mean time/pressure response preserve the instantaneous mean
compatibility structure. This does not preserve an evolving trajectory or
guarantee the geometric shape-rate conditions.

The independent check uses 18,720 points. It combines the stored original
wave residual, a directly evaluated five-point finite-difference change
for the locked wave, and exact linear tangent response columns. It is not
a fresh full-field finite-difference replay of each enriched candidate.
The report records those source distinctions and coefficient snapshots.
The degree-2 reconstruction reproduces the earlier independent result.

`Mode2TangentCorrection` in the script provides the physical affine-time
wrapper for later combined-field replay. Neither this residual reduction
nor the wrapper's reference-time check establishes PDE or scale-recursion
acceptance. All residuals remain far above 1e-3.

```powershell
python experiments/root_st073/wave_higher_harmonic_tangent.py
```
