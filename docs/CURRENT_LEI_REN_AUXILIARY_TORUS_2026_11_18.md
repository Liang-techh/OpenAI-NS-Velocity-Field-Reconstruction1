# Lei–Ren Part II auxiliary correction operators — 2026-11-18

ST073 remains the reconstruction target. The canonical reserved `m15` source proof in the [previous ST073 handoff](CURRENT_RESERVED_ACTUAL_M15_C1_2026_11_03.md) remains the basis for the forward-chain comparison. The new code implements a part of the subsequent oscillatory-correction framework.

## Paper and implementation

The supplied *Finite-Time Blowup for Navier–Stokes with Smooth Forcing, Part II. Residual Correction via Oscillatory Pulses* is [arXiv:2611.06869v1](https://arxiv.org/abs/2611.06869v1), 105 pages, SHA-256 `979e65cc4fd0b236bab04cef0284fa5808000962fa8ffa83923f3345c3519e5a`. Section 1.1 requires a complete Part I background with an admissible stress and flat remainder. The project’s `m15` radial moment, `P1`, and Part II’s nonzero angular Fourier modes are distinct sources. Part II §4.3 corrects defects of a later, axisymmetric mean velocity and preserves its linear moments.

[`lei_ren_auxiliary_torus.py`](../src/openai_ns_reconstruction/lei_ren_auxiliary_torus.py) implements three operations on caller-supplied data at fixed slow coordinates:

- Eq. (2.45) and Lemma 2.5 (pp. 25–26). Pull a finite Fourier spectrum from level `i` to a common level `i0` by mapping the integer mode `n` to `(E^(i-i0))^T n`. The zero mode remains the normalized torus average.
- Eqs. (4.30)–(4.31) (p. 64). Invert the physical fast-direction derivative on zero-mean data with the eigenvalue factor `Λ^(-i)`. The implementation rationalizes near-cancelling quadratic irrational divisors using their exact integer norms. To generate the tangential correction, the caller supplies the negative torus oscillation of the angular-mean force.
- Eqs. (4.33)–(4.35) (pp. 65–66). Compute the translated radial primitive and its physical radial derivative. The derivative retains the term `-χ_r` times the full translated integral. For the axial correction, the input is `h = r B_aux` and the potential is `I₀h/r`. The numerical quadrature error estimates are diagnostics. The source and cutoff must have the support and regularity specified by the interface.

## Checked result

The new fixture uses two adjacent bands (`k = 64,65`, covering levels `47,48`) at the paper-admissible manufactured value `δ = 1/256`. It compares the pulled-back spectra, both directional inverses, a near-resonant Pell mode, and the transported primitive against an independent Gauss–Legendre rule along the physical radius path. It also compares the physical radial derivative to an independent centered difference and checks the zero extension and rejection of a nonzero source mean.

| Quantity | Measured error |
| --- | ---: |
| Cross-band time and radial inverse, relative | ≤ 2.3e-16 |
| Cross-band primitive, divided by radial scale | 2.08e-17 |
| Independent radial quadrature, divided by radial scale | 4.51e-17 |
| Cross-band physical radial derivative | 2.64e-16 |
| Physical radial derivative vs centered difference | 5.96e-11 |

[Frozen run receipt](../experiments/root_st073/lei_ren_auxiliary_torus_receipt.json). One new run completed. Historical calculations were not repeated. The run pinned the Windows input bytes; publication normalizes CRLF to LF. Only line endings differ from the frozen Python inputs. The local continuation retains both byte hashes and the invoked files.

## Remaining work

Provide actual common-torus Fourier data for the tangential angular-mean force before physical evaluation. Then assemble both components of the auxiliary velocity from the corrected azimuthal velocity and the curl of `(I₀(r B_aux)/r) e_θ` with the axial derivative retained. At present the fixture exercises an operator with manufactured input. The existing uncommitted C1 companion work stays separate from ST073.

The critical ST073 task remains binding the six inner/lower contributions and active/compact histories of the forward 212 ledger to the canonical C1 functions in `reference*A*Rv^2` units. Preserve the signed sources and all inherited finite-`N` errors, add the reserved own-swirl once, and keep the independently defined negative full pulse as the comparison target. Independent terminal matching and the full Navier–Stokes residual remain open.
