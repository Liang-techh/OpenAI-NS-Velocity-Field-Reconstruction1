# Shape-constrained degree-3 second harmonic

The enriched tangent combines the useful mode-2 spatial extension with
moment, cone and instantaneous shape-rate constraints. There are 236
controls: the original 108 mode-0/mode-1 controls and 128 real degree-3
mode-2 derivative/pressure controls. The wave and instantaneous velocity
remain unchanged. The older 180-control replay must not be used directly.

The feasible QP reduces training momentum L2 from 2,869,199.2442 to
2,643,162.9080 (about 7.88%). Its retained objective rank is 227, four
moments are eliminated exactly, and phase-I slack is zero. The moment
maximum is 4.89e-9 and all 27 locations / 81 cone inequalities pass, with
minimum margin 9.9999988e-5. The fractional shape-rate floor is 10 per k;
radial and aspect constraints are active to floating-point precision.

The degree-2 shape rows and reference observables recomputed during
enrichment match the saved rows exactly. This is a consistency check of
the sampled affine operator, not a new continuum or trajectory result.

`enriched_shape_replay.py` builds the actual 236-control field by combining
the first 108 controls with a separate `Mode2TangentCorrection` of degree 3.
The replay reads `enriched_shape_candidate_snapshot.json`, whose SHA-256
matches the frozen optimizer report. Independent actual-field momentum
replay now gives L2 2,633,294.4983 and maximum 1.6832019e11 on 18,720 points.
Compared with the larger-margin degree-2 candidate, these decrease about
7.52% and 19.95%, respectively. Mean harmonic squared-residual share is
80.33%, mode 2 is 18.92%, and mode 1 is 0.74%; the mean is the principal
remaining spatial correction target.

Actual fixed-cylinder forward geometry from k0 to k0+1e-6 also retains
all three desired signs:

| Observable | Reference | Forward endpoint |
|---|---:|---:|
| Radial RMS | 0.0011358189224 | 0.0011358083478 |
| Axial/radial RMS ratio | 0.2085956309 | 0.2085972545 |
| Weighted angular speed | 728003.8958 | 737094.7434 |

Peak swirl still decreases (6569.41 to 6559.18), and the physical increment
is about 1.692e-10. Reference-time matching for this enriched candidate is
assembled, not a new independent compatibility replay. No endpoint
matching, accepted NS step, continuous monotonicity or scale recursion is
claimed. Full momentum remains far above the required absolute threshold.

```powershell
python experiments/root_st073/enriched_shape_tangent.py
python experiments/root_st073/enriched_shape_replay.py --source experiments/root_st073/enriched_shape_candidate_snapshot.json --mode momentum
python experiments/root_st073/enriched_shape_replay.py --source experiments/root_st073/enriched_shape_candidate_snapshot.json --mode shape
```

The candidate snapshot preserves the replayed coefficients. To reproduce
a new optimizer run directly, use `--source experiments/root_st073/enriched_shape_tangent.json`
and separate output paths so earlier evidence is not overwritten.
