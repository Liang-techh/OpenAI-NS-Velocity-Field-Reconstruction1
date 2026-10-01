# F21: fresh Md40 degree-110 local core and controlled exit

Date: 2026-10-01. This supersedes F20's finite-core-pending status. Source and datum are unchanged.

## Concrete result

All 110 radial orders have been generated directly in r=Lambda R, from 114 ordinary axial Taylor coefficients over real centers Z in [0.49,0.51]. Final rows retain four axial coefficients. Cumulative recurrence compute time was 693.621 seconds; the resumable exact-endpoint state is about 12 MB. No old finite tensor is reused.

The seed uses actual directed Md40 delta, the same fourteen-atom pressure envelope and the implicit positive axis amplitude of the admitted analytic core. An independent seed acceptance companion checks the pole-free tube, scalar Gbar, unique anchor root, pressure normalization, all 114 fixed jets and positive-upper swirl bounds. The 110-order run completed successfully and its state hash is bound to the receipt.

Finite rows are combined with F20's same-source degree-110 mixed C3 tail. Twenty mixed derivative enclosures at scaled radii 4 and 4.1 are available. The full normalized swirl shape Phi is positive at both exits:

| scaled radius | full Phi enclosure | full Uz enclosure |
| --- | --- | --- |
| 4 | [0.04690343900766381, 0.5177119157013341] | approximately [1.96000000000001, 2.04000000000001] |
| 4.1 | [0.02570662042941482, 0.5153695327916064] | approximately [1.96000000000001, 2.04000000000001] |

These displayed endpoints are rounded summaries; exact directed endpoints are in the receipt. The width mainly describes the local axial family and conservative interval propagation. The finite truncation error is a separate analytic bound, not the displayed enclosure width.

## Amplitude and pressure fidelity

The physical F0 value is NOT materialized or numerically selected. It remains the strictly positive analytic function exp(-logC-Lambda G), with G anchored at the unique H root in [-j,0]. On the certified complex tube, logC=Lambda Gbar+2logLambda+1000 gives |F0| <= exp(-2logLambda-1000).

The scaled swirl source is Sscaled=epsilon^2 F0^2. Its complex supremum is bounded by exp(-6logLambda-2000). Ordinary axial Taylor coefficients are enclosed by Cauchy on eta/2 disks. The zeroth coefficient is [0, positive upper] and higher coefficients are symmetric nonzero-upper intervals. A zero lower bound is an enclosure limitation, not a choice F0=0. No small swirl term is discarded.

The scaled pressure seed is epsilon P0=(P0/Pstar^2) exp(2logPstar-logLambda), so enormous Pstar^2 need not be generated. The callable exit restores pressure from the same swirl primitive:

    epsilon P = epsilon P0 + epsilon^2 F0^2 integral_0^r Phi(s,Z)^2 ds.

No independent pressure fit or added tail is introduced. Six normalized radial integral quantities needed to recover the five matching moments are available through axial order two, including separate Uz-squared and weighted Phi-squared pieces for the combined fifth moment. Uniform analytic tail bounds are integrated and product errors retained. These are core moments; their terminal target identities have NOT yet been repaired.

## Artifacts and reproduction

All under `experiments/root_st073/`, prefix `lei_ren_part1_paper_`:

- logarithmic_interval_core.py: source-guarded resumable factory.
- logarithmic_interval_core_Z049_Z051_state.json: exact fixed jets and all finite rows.
- logarithmic_interval_core_Z049_Z051.json: completed receipt and state SHA.
- logarithmic_core_seed_check.py/.json: seed acceptance and root/tube/units guards.
- logarithmic_core_exit.py/.json: callable mixed derivatives, controlled core moments and restored scaled pressure.

Run in that directory, or use repository-relative script paths:

    python lei_ren_part1_paper_logarithmic_interval_core.py --seconds 30
    python lei_ren_part1_paper_logarithmic_core_seed_check.py
    python lei_ren_part1_paper_logarithmic_core_exit.py

The first command resumes the existing matching state, and exits without further updates if degree 110 is already complete. State and source identities must match. Downstream work should consume the F21 core/exit receipts, not promote F20's historical finite_core_generated=false field.

## Next tasks

- [x] Direct scaled recurrence with physical recurrence comparison at representable scales.
- [x] Same-source 114-jet seed, nonzero swirl envelope, root/tube/pressure guards.
- [x] Fresh 110-order finite local core with exact resumable endpoints.
- [x] Combine finite Phi/Uz with same-source mixed C3 analytic tails.
- [x] Provide controlled core moments and pressure primitive through axial order two.
- [ ] Adapt the comparison ODE to these normalized rows, delta and implicit amplitude; never load old Lambda/logC or require a materialized F0 point value.
- [ ] Integrate the switched comparison with a directed integral enclosure. Preserve the small-radius core tail and verify positive Phi throughout the matching segment.
- [ ] Build fresh transition/annulus data and normalized physical inlet transfers. Keep all symbolic amplitude/Lambda factors until algebraic cancellation is justified.
- [ ] Rebuild the five-bump functional inverse and terminal identities on the same axial family; validate functional errors and source hashes, not isolated point samples.
- [ ] Bind the repaired reference inlet to F19 O.2 and extend radial/axial cone scope. F19 remains conditional until this step is complete.
- [ ] Restore the angular correction, analytic heat-pressure compatibility and exact heat exterior with logarithmic/arbitrary-precision scales.
- [ ] Construct the admissible stress lift and flat remainder, then true temporal n-dependent recursion, mean/oscillatory corrections and independent Cartesian full residual checks.

This remains a local analytic coefficient enclosure. Physical-R derivatives, whole-axis validity, a fully selected physical velocity sampler, terminal moment repair, assembled background, admissible stress lift, temporal recursion and oscillatory correction are not certified here. The long-term objective remains unchanged.
