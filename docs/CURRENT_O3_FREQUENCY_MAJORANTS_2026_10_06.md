# Source-bound finite-frequency derivative and shear-error envelopes

Successor: [CURRENT_PATCH_SUPPORT_AND_MODULATION_BOUNDS_2026_10_06.md](CURRENT_PATCH_SUPPORT_AND_MODULATION_BOUNDS_2026_10_06.md), commit [cf89c259](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/cf89c2591a9e785a9f8063ca202c7d1474fef86b), closes the six patch support source-function gap and supplies actual whole-support modulation source norms. Global physical interface exports, common N/cones and recursion remain open.

Implementation: commit [2dd8180d](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2dd8180d2774757e8a0c2ee88d765c2db3e5c902).

The supported O2/O3 modulation now has callable analytic envelopes with explicit N dependence. They cover ordinary logR derivatives0..4 and axial derivatives0..5 of theta/Pstar, the swirl increment and axial/Pstar, and the errors of both actual local shear components against the same periodic loop. This is BOUND1a1. The source norms supplied to the formulas remain conditional; they have not been certified on the actual whole domain.

## Actual definitions and units

The source is current_O3_finite_frequency_profiles.profile. Its actual translated phase is N*t, with t=buffer_offset-11 on O2 and t=offset on O3. The cutoff is chi(t)=sigma(t+2)*(1-sigma(4*t-1)), supported on[-2,1/2], with plateau[-1,1/4]. N is a positive integer and mu is the same positive slow source parameter.

The actual exponent is A=-mu*chi^2*sin(4*pi*N*t)/(8*pi), with G=exp(A/N). The modified fields are E0*G and -sqrt(mu)*E0*chi*cos(2*pi*N*t)/(2*pi*N). Here E0 means the original theta/Pstar source. Ordinary derivative rows are converted to Taylor coefficients once for the exponential and converted back once; no second factorial or chart-width factor is introduced.

profile_majorants(mu,N,cutoff,theta_source,log_theta_y) takes five cutoff norms C_j, a5-by6 ordinary original-source norm grid F_jk and L1 bounding |d_logR log(E0)|, all on one common domain and independent of N. It rejects negative/infinite inputs, noninteger/nonpositive N and explicitly N-dependent symbolic slow norms.

Let a=mu*C0^2/(8*pi). The exact conditional inequalities include:

- |actual theta shear - periodic theta shear| <= mu*C0*C1/(2*pi*N).
- |actual axial shear - periodic axial shear| <= exp(a/N)/N * (2*sqrt(mu)*C0*a + sqrt(mu)*(C0*L1+C1)/pi).

They hold at every actual phase. The exponential difference uses its defining real integral and |exp(h)-1|<=|h|*exp(|h|). Small exponential increments are not obtained by subtracting rounded numbers near1.

All higher derivative envelopes use the actual ordinary Leibniz coefficients and a positive exponential Bell recurrence. After factoring exp(a/N), the largest generic N power in both swirl-increment and axial derivative envelopes is:

| logR derivative order | Largest N power |
|---:|---:|
|0|-1|
|1|0|
|2|1|
|3|2|
|4|3|

Thus increasing N reduces the local shear error while increasing some high derivative bounds. These formulas alone do not select a common admissible construction frequency; actual source parameters, implicit repair, pressure/tensor completion and cone margins must be combined.

## Execution and evidence

Files under experiments/root_st073 use prefix lei_ren_part1_paper_compliant_current_O3_frequency_majorants: .py, _check.py, .json, _check.json. Run the producer and focused checker directly. No accepted ancestor constructors are needed. This conditional module does not promote the reconstruction controller's whole-field or cone gates.

The receipt records61 source/derivative identities,276 independent direct-differentiation/shear comparisons and6 invalid-input rejections. Checks bind the actual N phase, harmonic rates, cutoff widths, Taylor conversions and shear programs. Arbitrary local exponential jets are compared with an independently truncated exponential series. Synthetic smooth nonconstant cutoff cases exercise the conditional formulas; they do not certify the current graph's actual whole-domain source norms. Working/index hashes matched all7 dependencies.

## Detailed next work

- [x] **BOUND1a1:** actual source program bindings, conditional spatial/axial derivative envelopes, exact finite-N shear deviations and explicit N-growth degrees.
- [x] **BOUND1a2:** certify actual cutoff C_j on the whole modulation support. Consume the squared-exponent original sigma definition and the analytic tail bounds from flat_pulse_derivatives/check; retain its factor4 on the right cutoff. Do not use a different sigmoid definition.
- [x] **BOUND1a3:** derive original theta/Pstar F_jk and logarithmic derivative L1 bounds on the same O2/O3 source and all Z[-1,1], with its actual amplitude and original nonnegative J. Retain the common amplitude at Rd and do not replace it with a sampled coefficient.
- [ ] **BOUND1a4:** substitute those certified source bounds into this envelope, preserve positive mu and one N parameter, and publish actual full-domain profile bounds. The source and implicit controls must match the checked modified graph.
- [ ] **BOUND1b:** add independent repair bump/control/partial-history bounds on every directed support piece. Local flatness does not erase cumulative moments or pressure.
- [ ] **BOUND2a/BOUND3:** apply the actual signed tensor and pressure operators and preserve phase correlations. Compare the completed stress errors with the admitted loop margins; generic absolute derivative caps alone cannot prove cone signs.
- [ ] **COMMONN:** combine repair contraction, supported shear, taper, pressure, radial/diagonal completion and joins into one sufficient finite integer. Keep N=10^12 labeled repair-only until this is done.

Whole modified cones, total physical energy, true n-dependent coefficient recursion, oscillatory/mean corrections and full corrected NS remain open. These changes do not establish a lower full NS residual or a blow-up theorem.

BOUND1a4a analytic substitution is complete in the successor. Actual runtime parameter/physical factor binding (BOUND1a4b), edge correlations and independent repair/completed stress/common N remain open.
