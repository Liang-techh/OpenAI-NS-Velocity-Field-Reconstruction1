# Actual pulse higher radial sources and pressure C3 joins

The actual main/gap/end sources now expose U_yyy/U_yyyy and m_yyy/m_yyyy through ordinary Z5, Q_yyy/Q_yyyy through Z4, and P1_yy/P1_yyy through Z3. Both internal seams and the Rp/Rv exterior seams match these higher pressure sources. This completes a missing dependency of the original n1 source assembly. The full goal remains **ACTIVE / INCOMPLETE**; genuine n>=2 coefficient recursion and corrected Navier–Stokes are not yet installed.

## Construction and conventions

The source uses the accepted current Rp native frame, selected coefficients, original gp integral and normalized end bumps. The original flat sigma and beta third/fourth derivatives are actual functions, with interior bodies and explicit inactive branches. Inactive quotients are not evaluated at support boundaries. The checker proves those interior bodies against repeated differentiation of the exact scalar definitions, then consumes their derivative relations in the actual Leibniz and velocity binomial formulas.

Main/gap physical-y differentiation is mu*d/dxi; end differentiation is d/ds. The velocity binomial expansion uses bp=1/2+mu. The saved original scalar Taylor APIs divide derivatives by n!; the new source bodies are already ordinary derivatives. The factors6/24 are recorded and not applied twice. End order3/4 bumps have ell^4/ell^5 denominators with the original normalization and end scale.

Higher mass uses the same actual ODE m_y=U-m:

```text
m_yyy  = U_yy  - m_yy
m_yyyy = U_yyy - m_yyy
Q_j = [2Z U_j - (1-delta)Z m_j - (1-Z^2) d_Z m_j] / (1-delta Z^2)
```

This is leading 1-delta recovery. The accepted own positive-order velocity recovery with 1+delta remains separate. No additional Q1_Z3 row is fabricated.

The same retained pressure FTC gives P1_y=-R*S/2+Q_y+Q_yy. Since physical R_y=R, its next derivatives are:

```text
P1_yy  = -R*(S+S_y)/2          + Q_yy  + Q_yyy
P1_yyy = -R*(S+2S_y+S_yy)/2   + Q_yyy + Q_yyyy
N1p_y  = (P1_yy-P1_y)/R
N1p_yy = (P1_yyy-2P1_yy+P1_y)/R
m15_yy  = -R*(P1_y+P1_yy)
m15_yyy = -R*(P1_y+2P1_yy+P1_yyy)
```

All terms and radius derivatives are retained. No extra mu multiplies R_y. Absolute pressure values are retained as the earlier actual m13/m15 primitives and finite preceding histories; derivative matching does not reset a pressure value.

## Quiet source and seam composition

The saved quiet modified flow and leading flow are distinct source packets. Their equality is now proved explicitly: actual compact profile corrections are flat at x=2, the pinned relative terminal certificates zero only the correction histories, and the actual modified/leading terminal input recipes agree. Applying those proved input identities to both saved power flows gives identical V/m rows through y4 and Z5.

The original linear-Z native source gives normalized terminal mean kappa*Z and K=R2*Pstar*kappa. Applying Pstar once yields m_j=(-1)^j*KZ/R and Q_j=(-1)^(j+1)*K/R, with every U_j=0. The checker evaluates the actual saved input/flow/recovery recipes before binding any recovered Q input in the pressure map. The complete incoming pressure is checked query independent and retained as an arbitrary nonzero function in the full closed-pressure formula.

Main/gap is xi=11; gap/end is xi=13-4mu versus s=-4. Accepted actual low-order fields and radius identities compose with the newly proved flat high velocity and mass/Q source identities. Q may depend on Z at these internal seams; its Z derivative is not set to zero. Six actual pressure/density functions match through Z3 at each seam. Rp uses the actual quiet endpoint and current main inlet. Rv keeps the same saved terminal pressure value and proves that its radial derivatives vanish from query independence.

This is exact source-function composition. Direct whole-function comparisons, direct ordinary-row comparisons and higher rows obtained by pinned ordinary derivative commutation are reported separately. It is not an original-scale numerical enclosure or a claim that every ancestor was expanded again.

## Scoped evidence

| Evidence | Scope |
| --- | --- |
| 8 scalar source bodies | Actual sigma/beta orders3/4, with exact support guards |
| 12 direct whole axial identities | B3/B4 and U3/U4 on three pulse charts |
| 60 higher axial rows | Pinned ordinary derivative commutation |
| 16 flat support-boundary identities | Main0/main11/end-4/end0, B3/B4 and U3/U4 |
| 36 actual high mass ODE rows | m3/m4 through ordinary Z5 |
| 15 whole Q recovery identities | Q0..4 on three charts; 60 higher rows by commutation |
| 6 actual radius identities | R_y=R and R_yy=R on three charts |
| 24 whole higher pressure identities | Eight pressure/density derivatives per chart; 72 higher rows by commutation |
| 48 compact terminal input rows | Modified versus leading E/V/m/h/k/e/p/P0 through Z5 |
| 12 flat compact correction rows | Actual delta_E/delta_V endpoint recipes |
| 60 quiet flow comparisons | Saved modified/leading V/m rows y0..4 through Z5 |
| 85 quiet physical source rows | U/m/Q y0..4 with their actual saved ordinary row lengths |
| 32 quiet higher pressure rows | Eight derivative functions through ordinary Z3 |
| 8 retained full quiet pressure rows | Incoming pressure plus the signed closed increments |
| 48 internal higher pressure rows | Six functions through Z3 on main/gap and gap/end |
| 48 Rp/Rv higher pressure rows | Same six functions through Z3, both actual sides |

Six proof groups pass in108.422 seconds. The archive has289183 exact source nodes. Exact modified/leading quiet row lengths and pinned ordinary derivative algebra are checked before seam jet abstraction. The public `--restore-source` command passes in a fresh process with checked=true and source_only_restored=true. Numerical oracle availability remains false.

## Recovery

```powershell
python scripts/run_original_pulse_higher_radial_pressure_C3_functions.py
python scripts/run_original_pulse_higher_radial_pressure_C3_functions.py --restore-source
```

Use `CurrentOriginalPulseHigherRadialPressureC3Functions()`. The source, checker, exact archive and receipt share `lei_ren_part1_paper_compliant_current_original_pulse_higher_radial_pressure_C3_functions` under `experiments/root_st073`. `charts`, `quiet`, `endpoints` and `outer_pressure_traces` retain the actual source functions. These are symbolic source accessors; they do not provide an admitted numerical oracle.

## Next tasks for agents

Each task needs its own source/receipt evidence. Mark a task complete only when its full stated scope is supported; record the commit and command. Consume accepted ancestors instead of replaying historical broad suites.

1. [x] Construct and identify actual sigma/beta third/fourth scalar derivatives and support guards.
2. [x] Construct actual pulse U/m/Q third/fourth radial source functions with ordinary derivative conventions.
3. [x] Construct actual P1_yy/P1_yyy, N1p_y/N1p_yy and weighted m15 derivatives through Z3.
4. [x] Prove quiet modified-to-leading input/flow equivalence, original K recovery and retained full pressure.
5. [x] Prove both internal and Rp/Rv higher pressure source joins; keep nonzero means and pressure values.
6. [x] Construct actual compact/reserved physical-y velocity and mean sources through y4, preserving the original controls, normalization and complete histories. Completed in [CURRENT_ORIGINAL_COMPACT_HIGHER_RADIAL_PRESSURE_C3_2026_10_14.md](CURRENT_ORIGINAL_COMPACT_HIGHER_RADIAL_PRESSURE_C3_2026_10_14.md).
7. [x] Recover actual compact leading Q through y4/C4 and reserved own Q through y4/C2 from genuine mean rows, keeping1-delta and1+delta separate. See the successor source/receipt.
8. [x] Construct actual total own P1_yy/P1_yyy and weighted m15 derivatives at2Rc, including all R and angular repair product terms. Sixteen total own pressure rows and24 leading higher pressure rows match the saved quiet side. Prior nonzero values/first derivatives remain.
9. [x] Extend actual own quiet/main/gap/end fields to y4 C3/QV C2, proving saved M1 query independence first and preserving the own1+delta recovery. See the successor source/receipt. The broader combined outer mixed contract remains task10 below.
10. [ ] Write an explicit finite derivative contract for the combined outer provider: each field, each physical-y order, each ordinary Z order, both sides of every seam, and the source theorem proving it. Separate finite function identities from original-scale numerical bounds and qualitative all-order smoothness.
11. [ ] Complete the required own outer mixed rows under that contract, including any missing primitive/density rows. Existing leading pressure matching alone is not complete own velocity smoothness.
12. [ ] Publish the combined outer seam receipt/accessor only after all contracted2Rc/Rp/Rv rows and the pressure conventions are proved. The historical normal outer receipt is still absent; do not create it from the204 ordinary rows alone.
13. [ ] Finish actual own Rin/rminus source joins with the original inner Volterra/cutoff/lower collar. Check both physical radius recipes, nonzero histories, leading versus own signs, pressure normalization and required mixed derivatives.
14. [ ] Finish the actual own Rc join between intermediate/compact sources, retaining the actual selected coefficients, common phase/N, all signed histories and every finite pressure boundary term.
15. [ ] Assemble a checked global n1 source provider over all original charts and own supports. Document query ownership, axis parity, finite-energy radial tails and the genuine derivative contract. Local chart acceptance is insufficient.
16. [ ] Restore the numerical oracle for the original selected branch and transfer source enclosures for all required derivatives. Keep exact source values distinct from bounds, coordinate scale factors and comparison covers.
17. [ ] Establish one common global N and the original quantitative containment/contraction requirements. Do not choose independent per-chart frequencies or use a scale-reduced numerical fixture as the original admission.
18. [ ] Complete axis regularity, physical-volume energy and radial-tail diagnostics, analytic preheat pressure data, heat collar and exact heat exterior from the same global source family.
19. [ ] Construct and check the admissible stress cone and separate flat remainder with norms/margins on every relevant chart. Background residual magnitude is not the final corrected NS gate.
20. [ ] Implement true n>=2 n-dependent recovery, independent moment repair and finite-order remainders on a common inner interval; then justify the smooth sum. Coordinate rescaling of n1 is not coefficient recursion.
21. [ ] Implement mean/oscillatory two-family correction from admissible stress and check averaged quadratic cancellation with independent error estimates.
22. [ ] Export corrected Cartesian [u,v,w](x,y,z,t) with the fixed forcing, and verify independent NS residual Linf and physical-volume L2 below1e-3. Report core radial/axial exponents, relative elongation and material winding across time scales.

Paper reference: supplied Lei–Ren Part I §13.4 and §15.1–15.2 for physical-y recovery, signed pressure histories and the qualitative smoothness assertion. The new source functions fill finite derivative dependencies; they do not turn that qualitative assertion into quantitative bounds or prove blowup numerically.
