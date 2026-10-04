# Original pulse-end full meridional similarity stress — 2026-10-04

The pulse-end companion now computes the original similarity stress on the complete end chart `s ∈ [-4,0]`, `Z ∈ [-1,1]`, including both original beta supports centered at `-3` and `-1` with width `0.15`. It preserves the current selected source coefficients, original velocities, five partial moments, analytic absolute pressure datum, and signed incoming angular history. This extends the stress reconstruction into a region where radial and axial velocities are generally nonzero.

## Implementation and evidence

`experiments/root_st073/lei_ren_part1_paper_compliant_pulse_end_stress_C3.py` consumes the current native pulse mixed-C4 data and accepted flatten stress data. It reconstructs all terms of the original full meridional stress formulas (3.16)–(3.18), including both shear terms, axial transport, nonlinear meridional transport, backward linear moments, and selected quadratic energy loss. The result stores four theta and six axial stress sectors with mixed derivatives of total order at most three.

The physical factors `B`, `D`, `H`, and the original radius remain exact logarithmic factors. Their very small nonzero contributions are retained; runtime caps are only enclosures and never defining field values. Complete native future energy and its half-normalization are preserved. Pressure uses the accepted analytic absolute datum transported inward by its original differential identity.

The amplitude bridge replays the original entry, angular, power, and flatten source recipes. It establishes `K0(Z)=C0/(1+Z²)` and `C0*theta_base*exp(-bh*q0)=1`, then consumes the current production radius and absolute-pressure unit proofs. It explicitly replays `Ev0=Pstar*U*exp(-13/(2*mu)-13)` and the production `logEv2_parts`, proving that pulse and flatten use the same physical amplitude and pressure-square scale. This is a source-unit connection, not the full functional interface proof.

The focused checker recomputes the current source report and verifies 500 finite signed sector rows. An independent direct-integral fixture checks 30 full-stress and radial-source mixed derivatives at points inside each support and between them. It retains nonzero meridional histories and signed incoming angular memory; its tolerance is `1e-55`. Exact symbolic comparison also binds the full formulas and the ten sector derivative rates to the original evaluator.

Run only this new stage with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulseendstress
```

## Next dependencies and limits

1. Prove actual `s=0` / flatten `t=0` functional equality of K4, all five moments, stress mixed3, and pressure mixed4. Use original radius, selected moment equations, empty future beta supports, and the complete native future energy. Interval overlap and scalar unit equality are insufficient.
2. Prove the required derivative joins at all four original beta support endpoints. Keep the complete end chart and ordinary `s=logR` derivatives.
3. Build a physical adapter with all meridional terms, original viscosity units, completed tensor and diagonal, and actual three-component remainder. The downstream pure-swirl adapter alone is insufficient. Close the completed physical interface at flatten `t=0`.
4. Establish continuous admissibility through both supports with signed correlated source histories and exact scales. Extend the same companions through the pulse gap, main pulse, and entrance.
5. Resolve upstream finite-width bridge feedback, compose completed full-tensor global admissibility, and independently bound flat remainder, physical-volume norms, and required-domain energy.
6. Implement the actual coupled `n=1` and `n>=2` coefficient recovery, independent moment repair, finite-order remainder, and smooth summation. Then add mean/oscillatory stress correction and the independent full corrected Cartesian residual and measured dynamics.

The full pulse–flatten stress interface, pulse-end physical decomposition, pulse-end cone, full pulse coverage, global tensor admissibility, required-domain energy, temporal coefficient recursion, and oscillatory correction remain open. Existing downstream flatten and tail certificates retain their original regional scope.
