# Actual own Rin cutoff/lower source join, 2026-10-14

Status: finite source accepted. Five scoped groups pass in46.922 seconds, with138 contracted join rows and296133 exact graph nodes. Independent public restoration passes with checked=true/source_only_restored=true and exit code0. The full research goal remains **ACTIVE / INCOMPLETE**.

This construction joins the original cutoff endpoint to the canonical lower own chart, retaining its nonzero core, mean, pressure and all five moment histories. It consumes the existing Rc and outer source joins. The finite contract uses the available actual inner radial sources; the missing higher radial pressure source is listed explicitly below.

## Common pressure and coordinates

The original inner first leading packet is exactly the saved lower leading packet, including its histories and radial derivatives. Put s0=sc/4, rho=4 exp(hb s), R=epsilon rho and Lambda epsilon=1. The raw inner source and the lower source are

```text
N_inner = -S/2 + Lambda (2 Q_rho + rho Q_rhorho)
N_lower = -S/2 + (Q_yy + Q_y)/R
Q_y = rho Q_rho
Q_yy = rho Q_rho + rho^2 Q_rhorho
```

They therefore agree as whole functions. The ordinary Z derivatives through order3 use the admitted same-source derivative algebra and Z-independent coordinate conversion. The proof first compares the actual N values without binding them to an imposed common output.

At Rin the original cutoff argument is2, so sigma=1 and every cutoff derivative is flat. F1 and Uz1 vanish there, while the mean, pressure and accumulated moments generally remain nonzero. The true cutoff pressure derivative becomes epsilon N. Since d rho/ds=hb rho, the signed cutoff densities push forward to

```text
dm13/ds = hb R N
dm15/ds = -hb R^2 N
P1_y = R N
m15_y = -R^2 N
```

The actual lower partial integrals are empty at s0. Their addition preserves the full nonzero cutoff m13/m15 endpoints. No lower base-pressure value replaces the own pressure history. The own recovery retains1+delta; the distinct leading pressure source retains1-delta.

## Finite source contract

| Own field or history | Physical y orders | Ordinary Z orders |
| --- | --- | --- |
| Uz1, Utheta1, F1, M1, mean1, K1 | 0 through2 | 0 through3 |
| Q1, V1 | 0 through2 | 0 through2 |
| Total own P1 | 0 through1 | 0 through3 |
| All five cumulative m11–m15 | 0 through1 | 0 through3 |

There are138 actual left/right rows. The conversion uses D_y=rho d_rho, M1=R mean1 and the genuine square-root factor in Utheta1. Both state sides retain the original node references. The checker consumes the original flat sigma source, defining mean/pressure FTC bodies, signed measures and actual empty lower integrals. Binding a retained source history in symbolic algebra does not set it to zero or choose its numerical value.

This is an exact source join in the stated finite scope. It is not the full state-y4/pressure-y3 contract, an original-scale numerical velocity oracle, quantitative smoothness, a global n1 field, genuine n>=2 recursion or a corrected Navier–Stokes solution.

## Recovery

The source/archive/checker/receipt share prefix `lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join` under `experiments/root_st073`.

```powershell
python scripts/run_original_n1_Rin_source_join.py
python scripts/run_original_n1_Rin_source_join.py --restore-source
```

Use `CurrentOriginalN1RinSourceJoin()`. `left` and `right` retain the chart source rows; `join` contains their actual evaluations at Rin. `rc` retains the accepted Rc/combined outer providers. Ordinary recovery uses saved sources and receipts and does not replay numerical ancestors.

## Next construction tasks

1. [x] Admit the actual Rin pressure-source identity, original FTC and complete138-row join, then independently restore and publish the finite source and task queue. All five scoped groups pass, and the fresh public process restores the same296133-node source with exit code0.
2. [x] Construct the missing genuine leading Q_rhorhorho and Q_rhorhorhorho from the same analytic core and lower ratio/Psi/mean transport. Completed by [the actual lower radial4 / Rin246 source](CURRENT_ORIGINAL_LOWER_HIGHER_RADIAL_RIN_JOIN_2026_10_14.md), preserving the accepted radial0..2 prefix, original sigma/hb/rho and histories.
3. [x] Derive the actual pressure source derivatives, including the derivative of1/R. The new source independently differentiates the original pressure operator. In the regular rho form,

   ```text
   N_rho = -S_rho/2 + Lambda (3 Q_rhorho + rho Q_rhorhorho)
   N_rhorho = -S_rhorho/2 + Lambda (4 Q_rhorhorho + rho Q_rhorhorhorho)
   ```

   Prove the same result by differentiating the physical-R pressure form. Keep all leading U, Q and Q_Z terms in S.
4. [x] Extend actual own state to y3/y4 on the flat cutoff neighborhood by the original support theorem and mean/history FTC. The whole-flat-source proof checks K/Q1+delta/V/M and the independent mean flow while retaining the nonzero leading pressure.
5. [x] Build P1_yy=R(N+rho N_rho), P1_yyy=R(N+3rho N_rho+rho^2 N_rhorho), and the corresponding signed m15 rows. All246 finite Rin rows now match in the new source above.
6. [ ] Finish own rminus from the lower right endpoint to first-micro left, with actual empty first-micro increments and the same leading/own pressure gauge. Reuse the accepted leading source theorem without replacing own recovery by leading Q.
7. [ ] Enumerate and join the remaining own intermediate charts; compose a checked global n1 source and Cartesian coefficient export. Regional source availability is insufficient for a global field.
8. [ ] Restore actual original-scale numerical evaluation and derivative enclosures for the analytic, Volterra, implicit and integral sources. Keep evaluated values, exact source handles and comparison bounds distinct.
9. [ ] Establish common global N, quantitative mixed bounds, axis regularity, physical-volume energy and radial-tail diagnostics.
10. [ ] Complete the analytic preheat pressure datum, heat collar and exact heat exterior with pressure-compatible matching.
11. [ ] Construct admissible stress and flat remainder separately, with regional signs, cone margins, norms and scale dependence.
12. [ ] Implement independent n-dependent coefficient recovery and moment repair for n>=2 on a common inner interval; then finite-order remainder control and smooth summation.
13. [ ] Construct mean corrections and two-family oscillatory pulses and their averaged quadratic stress cancellation.
14. [ ] Independently validate the corrected Cartesian field with fixed forcing and both Linf and physical-volume L2 residuals below1e-3. Measure critical-time scale laws, anisotropic core contraction, vorticity growth and material winding separately.

Tasks2..5 are complete in [the new actual lower/high-pressure/Rin246 source](CURRENT_ORIGINAL_LOWER_HIGHER_RADIAL_RIN_JOIN_2026_10_14.md). The next construction is own rminus (task6), followed by other joins/global n1. The implemented recurrence below remains the mathematical source route; numerical and quantitative bounds remain open.

## Implemented exact radial recurrence for tasks2..5

Implement a new derived leading-only radial provider. In `current_original_n1_analytic_core_lower_pressure_C3_functions.py`, `d.extend(...,5)` extends ordinary Z rows, not radial order. Its accepted radial arrays still stop at2. `inner_radial_functions.CoreRhoDerivative` handles the separate n1 even system and axis and does not supply these leading collar rows.

The original `current_original_n1_core_fixed_point_functions.py` component recipe denotes the same unique analytic twenty-term core solution at every radial order. Preserve its full operator/source/Green-inverse definitions and source family when constructing radial3/4 nodes in the derived provider. Use a guarded derivative adapter admitting those same component definitions for ordinary Z0..5; do not mutate a saved ancestor owner or its admission set. Numerical enclosures and their centers cannot define a component.

Let cPhi and cPsi be the original core functions evaluated at rho(s), and let a_ratio be the retained collar ratio (distinct from physical radius R). Reuse its original integral and write

```text
chi = 1 + (hb-1) sigma(s)
A = cPhi_rho / cPhi
g = (chi-1) A
(log a_ratio)_rho = g
(a_ratio)_rho = a_ratio g
B = a_ratio chi
```

Generate chi's first three rho derivatives using D_rho=(hb rho)^-1 D_s and the actual `flat_source.sigma_jets`; multiply Taylor coefficients by the appropriate factorial. Obtain the a_ratio and B derivative arrays by ordinary product rules. For n=1..4,

```text
Phi_n = sum_{k=0}^{n-1} binomial(n-1,k) B_k cPhi_{n-k}
Psi_n = sum_{k=0}^{n-1} binomial(n-1,k) B_k cPsi_{n-k}
meanPsi_n = (Psi_{n-1} - n meanPsi_{n-1}) / rho
```

Retain the original Phi/Psi/meanPsi order0 primitives and the nonzero rho4 mean inlet. Prove agreement with every saved radial0..2 row before accepting orders3/4. Recover U_n=epsilon Psi_n and m_n=epsilon meanPsi_n for n>=1; their order0 rows retain U0. The leading map

```text
Q_n = (2 Z U_n - (1-delta) Z m_n - (1-Z^2) (m_n)_Z) / (1-delta Z^2)
```

then supplies genuine Q radial3/4, with ordinary Z0..4 from U/m Z0..5. It uses the leading1-delta weight. The own1+delta map remains separate.

For S radial0..2 put q_n=Q_n, u_n=U_n, y_n=rho q_{n+1}+n q_n and L=1-delta Z^2, d=1-Z^2, h=1-delta. Use

```text
S_n = (q_n + h Z (q_n)_Z/2 + y_n)/L
    + sum_{i=0}^n binomial(n,i) q_i q_{n-i}/2
    + sum_{i=0}^n binomial(n,i) q_i y_{n-i}
    + sum_{i=0}^n binomial(n,i) u_i
        (-2 Z q_{n-i} + d (q_{n-i})_Z - 2 Z y_{n-i})/L
```

Combine these actual S_n with the N0/N1/N2 formulas in task3. Derive radial pressure source rows from those bodies and the signed physical measures; a serialized C3 integral output alone does not supply a new radial derivative. After constructing these source rows, use the flat cutoff theorem for high own state derivatives and finish the full Rin246 contract. Keep the nonzero pressure primitive and incoming m15 datum throughout.
