# Current O3 variable direction and finite-frequency shear profiles

Current successor: [CURRENT_O3_MODULATED_HISTORIES_2026_10_06.md](CURRENT_O3_MODULATED_HISTORIES_2026_10_06.md), commit [e6e06575](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e6e0657567e776fb609f0f1d42e5296566d8b821). Own cumulative five histories, frequency-uniform axial bounds, common-axis pressure and radial recovery are now constructed. The remaining repair/common-N/modified-cone tasks are detailed there; the original profile-only statements below describe the earlier milestone.

Implementation and current source/profile receipts: commit [43320b19](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/43320b19d5681b115c7c552f8ce1755caa6145d8).

## Constructed result

The original whole O3 variable transition now has a continuous positive-theta and direction proof on offset [0,1], Z [-1,1]. All eleven signed stress sectors, full energy, ordinary variable log-amplitude derivatives, same absolute functional pressure and radial remainder are retained. Its conservative direction expression is below 5.076e-435.

The original strict shear fails at offset 0: `vs-2=2*mu*sigma(offset)=0`, while stress is nonzero. This is a real original-source obstruction, not a zero exterior. On offset (0,1] the original source is strict. The complete closed transition remains unadmitted.

A mean-preserving periodic shear loop, its exact primitives, and an actual seam-crossing finite-frequency local profile candidate are now callable. The latter changes both O2 buffer and O3 source velocities with the same translated log-radius phase. It supplies ordinary logR derivatives through order 4 with axial Taylor rows and all five moment-increment densities. It does **not** export modified pressure, integrated moment histories, radial velocity, completed tensor or an admitted common frequency. These require the next construction.

The original checked graph still has **15 strict nonzero whole regions plus a separate exact zero exterior; 17 regions remain without whole-region cone admission**. These are regional counts, not project percentages. Inventory remains 33 regions / 32 adjacent / 14 internal tensor traces; primitive atlas 14 / 8. The original background receipts remain proofs for the original source, not automatic proofs for the modified candidate.

## Actual variable direction

Let C=1/(1+Z^2), L=1-delta*Z^2, and define

`A=(1-delta/2)*C+(1-delta)*Z^2*C^2`,

`Nshape=((2*delta*Z^2-1)*C+2*(1-Z^2)*Z^2*C^2)/L`.

At the actual O2 buffer11 source use Ua, Ma, EZa, EQa and the source deficit

`da=(1-Xs)*exp(-(logPstar-1))>0`.

The same current original transition kernels give

`U=Ua*exp(-v/2-mu*J(v))`, `D=da*exp(-v+mu*J(v))`,

`Q(v)=integral_0^v sigma(s)*exp(s-mu*J(s))ds>=0`,

`X=1-D+mu*Q(v)*exp(-v+mu*J(v))`,

`K/U=M*exp(mu*J(v))-4D`, `M=Ma*exp(-v)`.

The exact theta operator is

`Theta=(A*X-C)/L+C*M+Nshape*K/U-(2+2mu*sigma(v))*C/R`.

The same M/K/X correlation supplies a strictly positive deficit term. With cX=(1-delta)/4 and cD=(15/8-14delta)/8, the source-bound drift, delta and inverse-R terms are smaller than the deficit reserve. Thus

`Theta >= cX*Z^2+cD*D(v)+positive_reserve > 0`

throughout the closed transition. Independent interval history boxes are not subtracted to infer the deficit.

Full energy has coefficients

`AZ=(EZa/Ua^2)*exp(2mu*J)`,

`AQ=(EQa/Ua^2-KE(v)/2)*exp(2mu*J)`.

The same absolute pressure is composed backwards from the actual power phase0 and common signed Rv datum:

`Pabs=-U(v)^2*C^2*H(v)/2+U1^2*(Pv+C^2/(2p))*exp(-p*(13/mu+Tw))`,

where p=1+2mu, U1=Ua*exp(-1/2-mu/2), and

`H(v)=[f(1)^2/p+integral_v^1 f(s)^2 ds]/f(v)^2`, `1/p<=H<=1`.

The functional Pv derivative is retained in the axial operator. The energy/pressure part has an exact odd-Z factor. AM-GM uses the exact cancellation `U(v)^2/D(v)=Ua^2*exp(-3mu*J(v))/da`; v cancels. The jet helper's logistic sigma and the history-integral positive fraction are proved to define the same reflected cutoff.

## Periodic target and actual finite-frequency candidate

For the uncut target loop, a0=2+2mu*sigma(v), b0=0 and

`aL=a0+mu*cos(4*pi*phi)`, `bL=2*sqrt(mu)*sin(2*pi*phi)`.

The zero-mean primitives are

`A=-mu*sin(4*pi*phi)/(8*pi)`,

`B=-E0*sqrt(mu)*cos(2*pi*phi)/(2*pi)`.

Here E0 is the actual current swirl, with its exact Pstar log factor. It is not a scalar cap. The same original inviscid vector gives a **changed** stress target:

`TthetaL/F=t0-(aL-a0)`, `TzL/F=Tz0/F+bL`.

The bound uses the changed ratio, not the old stress ratio. The original t0 reciprocal is at most exp(-1000); mu<1/6 is required. The target alignment is positive, its shear excess exceeds mu/2, and its conservative direction expression is below 4.031e-217. This is the target loop proof, not finite-N profile admission.

The callable local construction uses the common coordinate

`t=O2_buffer_offset-11=log(R/Rd)`, or `t=O3_offset`.

Its cutoff is

`chi(t)=sigma(t+2)*(1-sigma(4*t-1))`.

It is supported in (-2,1/2), equals 1 on [-1,1/4], and is flat at O2 buffer9 and O3 offset1/2. The old O2/O3 seam is inside the modulation. An O3-only cutoff that is zero there would leave the original obstruction unchanged.

Use delta-a=mu*chi^2*cos(4*pi*phi), delta-b=2sqrt(mu)*chi*sin(2*pi*phi). The corresponding primitives are

`A=-mu*chi^2*sin(4*pi*phi)/(8*pi)`,

`B=-E0*sqrt(mu)*chi*cos(2*pi*phi)/(2*pi)`.

For any finite integer N>=1, evaluate them at phi=N*t:

`E_N=E0*exp(A/N)`, `U_N=B/N`.

This phase is exactly N*logR translated by the constant -N*logRd. It is not an independent phase sample. The actual original axial velocity vanishes on both source charts. The profiles have exact shears

`a_N=aL-2*D_X^slow(A)/N`,

`b_N=exp(-A/N)*(bL+2*D_X^slow(B)/(N*E0))`.

The slow derivative of B retains both chi' and E0'/E0. All higher ordinary radial derivatives include the fast N factors. No exp(logPstar) or absolute enormous radius is materialized; both modified velocity coefficients are divided by the common Pstar.

On the whole original O3 interval, the supported target loop has a positive uniform shear floor: for v<=1/4, chi=1 gives at least mu/2; for v>=1/4, source monotonicity gives at least 2mu*sigma(1/4). The latter is the stored floor. The O2 taper has no inherited strict collar and remains a separate open proof. This does not establish the actual pN or actual finite-N cone.

Eight worked profile evaluations use N=10^12. This integer exercises the callable formula and is **not a certified sufficient frequency**.

## Five moment sources

Write e=E_N-E0 and u=U_N-U0; U0=0 in the modulation support. The exact per-dX increments are

| Moment | Increment density | Factored source |
| --- | --- | --- |
| M | u | Pstar |
| I | sqrt(2X)*e | Pstar*sqrt(2R) |
| J | sqrt(2X)*u*(E0+e) | Pstar^2*sqrt(2R) |
| S | u^2-E0*e-e^2/2 | Pstar^2 |
| Cp | (E0*e+e^2/2)/X | Pstar^2/R |

The code stores signed coefficient functions and the powers separately. The swirl increment uses the exact signed integral for exp(A/N)-1; subtracting nearly equal exponentials would erase it. Densities are not cumulative integrals. Until they have been integrated and independently repaired, neither the original pressure nor original radial velocity can be used for the modified source.

## Entry points and evidence

Prefix: `lei_ren_part1_paper_compliant_`.

- `current_O3_transition_direction_operator.py`: exact variable source, sigma binding, full signed stress/pressure and target loop algebra.
- `current_O3_transition_direction.py` and `_check.py`: producer, checked native API and complete source views. Focused controller stage: `currento3transitiondirection`.
- `current_O3_finite_frequency_profiles.py` and `_check.py`: actual local finite-N E/U, derivatives and density source construction. This is a separate candidate API.
- Corresponding JSON receipts: 39 exact variable/loop identities and 23 positive directed bounds; 17 exact finite-N/coordinate/moment identities.
- Producer/checker, focused controller, original source endpoint/interior APIs, inherited power and exact zero exterior, candidate chain-rule consistency, common modified source seam, flat profile edges and compilation passed.
- Working/index source audit: 823 dependency files matched.
- Read-only formula and implementation reviewer: GPT-5.6 Luna / max. No Astra child and no worker file edits.

With the existing warm checked graph:

```python
direction = CurrentO3TransitionDirection(powercone=checked_current_O3_power)
candidate = CurrentO3FiniteFrequencyProfiles(direction)
point = candidate.profile('O3_slope_mu', '.7', '.1337', 10**12)
left = candidate.profile('O2_buffer', '.7', 11, 10**12)
```

Do not cold-run `--stage all` to continue this task. Warm session 64009 retains the checked graph. Source hashes must remain consistent. Candidate construction checks must never turn the global profile/cone gates true.

## Next executable tasks, in dependency order

- [x] O3-power: complete actual correlated whole-domain cone and full energy/pressure.
- [x] O3-variable-direction: complete original closed-domain direction, positive stress and exact nonzero zero-shear endpoint diagnosis.
- [x] O3-periodic-target: actual-source mean-preserving loop and primitives; same inviscid vector with changed stress.
- [x] O3-local-finite-profiles: actual supported E_N/U_N source coefficients, logR/Z rows, common seam phase, flat profile edges and all five increment densities.
- [x] **MOD1 — Integrate the five signed increments.** Construct delta-m_N(X,Z) as functions, not fitted point tables. Retain axial derivatives at least through the orders needed by original stress and pressure. Use t=log(R/Rd), the exact dX=Rdt factor and source log sectors. Do not reset incoming original histories.
- [x] **MOD2 — Bound the high-frequency integrals.** Obtain source-dependent O(1/N) constants for the five increments and Z derivatives. Naive fixed-cell quadrature cannot resolve N=10^12; use analytic phase integration/averaging or justified integration by parts with all cutoff and amplitude derivative terms.
- [x] **MOD3 — Restore the common-axis pressure.** Integrate the changed Cp with the original axis datum and retain its Z derivatives. Recompute the same absolute pressure and radial recovery. Do not copy original pressure or add a pressure tail to lower a residual.
- [ ] **MOD4 — Fix a separate repair interval in the quiet O3 power.** A concrete candidate is q=Tw*power_phase in [1,2], after all modulation ends; normalize x=exp(q-1) in [1,e]. There U0=0 and E0=K(Z)*X^(-1/2-mu). Check the actual finite Tw exceeds 2 before installing it. Fix all widths and bump centers before choosing N.
- [ ] **MOD5 — Rebuild two axial and three swirl bump blocks.** For the current source, form the matrix in row order (M,J;I,S,Cp). The U weights are x^0 and x^(-mu); the E weights are x^(1/2), x^(-1/2-mu), x^(-3/2-mu). Keep the strictly positive mu: the U block becomes ill-conditioned as mu tends to zero. Prove its determinant/inverse using correlated exponent differences, rather than coarse independent boxes or an old fixed matrix.
- [ ] **MOD6 — Build the exact linear-plus-quadratic repair map.** Reuse bump integral and interval inverse algebra where applicable, but rebuild the weights, defect family and receipts for this interval and source. Solve B(Z)c+Q_Z(c,c)=-delta-m_N(Y1,Z), with C2 or stronger axial bounds. Existing SharedFiveMomentRepair is tied to a different source patch and cannot be copied unchanged.
- [ ] **MOD7 — Choose one finite uniform N.** Derive profile/shear constants, pN-p0 constants, cone Lipschitz bounds and repair inverse/quadratic bounds from the fixed source data. Choose N to satisfy all profile positivity, target margin, moment repair and corrected-cone inequalities together. A working sample integer is insufficient.
- [x] **MOD8 — Recover the modified radial velocity from its own moments.** Use the original meridional recovery and common pressure. Preserve exact divergence through the recovery structure. Export physical source rows only after modified M and its Z derivatives exist.
- [ ] **MOD9 — Prove affected joins as function identities.** At O2 buffer9 the profiles are unchanged and increments start at zero; at O2/O3 the modified profiles and phase agree; at O3 offset1/2 only the profiles are already unchanged while moment increments remain; at the repair exit both profiles and all five cumulative moments must equal the original source. Recompute tensors/remainders from those identities.
- [ ] **MOD10 — Resolve the O2 taper separately.** Its original shear is zero-excess and there is no inherited positive collar at the cutoff edge. Establish actual finite-N shear/direction there with the modified histories; do not assert that compact support alone preserves a strict cone.
- [ ] **MOD11 — Admit the modified closed O3 cone.** Use actual pN/pf and full stress, pressure, energy and derivatives, after repair and common-N selection. Promote only gates proved by this completed modified source graph. The unmodified regional cone count must not be used as the modified count.
- [ ] **GLOBAL1 — Finish the remaining inner/O2/current cones and smooth tensor/lift.** Rh, slope/turnoff/buffer, patch/restore/reshape/switch/bridge/core/axis still need their respective current whole-function arguments. Preserve 33-region source identities and full remainder.
- [ ] **RECUR1 — Implement actual n-dependent coefficient recovery.** Distinguish n=1 from n>=2, use the common inner domain and a separate five-moment repair at each order; establish finite-order remainder bounds and smooth summation. Coordinate rescaling of the leading profile does not implement recursion.
- [ ] **WAVE1 — Build both actual oscillatory velocity families and means.** Use positive amplitudes and flat supports, recover covariance and finite-error estimates, and then demonstrate averaged quadratic stress cancellation.
- [ ] **FULL1 — Remove the nonflat leading origin remainder and verify the corrected field.** Resolve physical u/v/w, independent Cartesian NS residual, finite total energy, core scale exponents and material winding. Keep this full original long-term goal active.

The leading-origin flatness problem, actual coefficient recursion, two wave families, corrected NS/energy and quantitative dynamical measurements remain open. Completion of local E/U formulas does not narrow the original project goal.
