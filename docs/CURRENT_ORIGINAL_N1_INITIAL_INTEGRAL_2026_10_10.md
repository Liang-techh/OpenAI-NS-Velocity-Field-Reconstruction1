# Genuine n=1: actual inner operator and computed initial regular integrals

The actual first-collar operator inputs and the first regular source-integral field are now computed. This advances beyond the preceding analytic majorants: source-owned interval integration returns F1, Uz1, K1, P1, d_x F1 and d_x Uz1 at the canonical core exit and the actual inner collar endpoint. The long-term reconstruction goal remains **ACTIVE / INCOMPLETE**. These are initial iterates, not the complete coupled n=1 solution.

Read this checkpoint before the older [domain/tail checkpoint](CURRENT_ORIGINAL_N1_ACTUAL_COLLAR_DOMAIN_2026_10_10.md). That analytic source domain, pressure binding and tail majorant remain accepted. Do not restart their ancestor validation.

## New source artifacts

Both quartets are in `experiments/root_st073/`:

| Stem | Accepted gate | Scope |
| --- | --- | --- |
| `lei_ren_part1_paper_compliant_current_original_n1_inner_operator` | `current_original_n1_actual_inner_operator_input_enclosures_installed` | Actual leading mixed-rho inputs, n=1 forcing, B0/B1 and g |
| `lei_ren_part1_paper_compliant_current_original_n1_initial_integral` | `current_original_n1_initial_regular_source_integrals_computed` | Computed Wbase=(I+J)Gg interval enclosures |

Each has `.py`, `.json.gz`, `_check.py` and `_check.json`. The JSON contains enclosures and source log sectors; it cannot replace the live accepted owners.

The operator uses the same canonical original core on rho in [0,4]. For s in [0,1/2], rho=4 exp(hb*s), it uses `OriginalWholeZBridgeSource.actual.packet`: the actual prescribed-shear Phi/Uz and its own six cumulative moments. The analytic core extension is used only for the known comparison direction in the bridge ODE. It does not replace the actual velocity or actual mean.

The M input is already Mz/R. Its radial rows follow the exact FTC identities

```text
mean_rho = (Uz - mean)/rho
mean_rhorho = (Uz_rho - 2 mean_rho)/rho.
```

The pressure retains Pstar^2 P0 plus (rho/Lambda) F0^2 C with the actual own normalized pressure integral C. Ordinary amplitude derivatives are included before exporting Z jets. Source bases are `(hb, F0_anchor, Pstar_squared, Lambda)`; none of their production logarithms is exponentiated to select a numerical parameter.

Matrix entries carry ordinary Z derivatives through order two; known forcing and g through order one. The available leading rho rows and their Z orders are exported explicitly. Missing higher derivatives raise errors; they are never zero-padded. The actual point moment history remains an enclosure, with `actual_point_moment_history_recovered=False`.

## Computed initial field

For D=(0,0,2,0,3,1), define

```text
(Gg)_i(x,Z) = integral_0^x (s/x)^D_i g_i(s,Z) ds.
J = G Akin, with Akin_15=1, Akin_26=1, Akin_36=-1.
J^2=0, so Wbase=(I+J)Gg solves the source plus kinematic rows.
```

The double integrals are collected before interval evaluation:

| Component | Integrated source kernel |
| --- | --- |
| F1 | s/2 * [1-(s/x)^2] * g5 |
| Uz1 | s log(x/s) * g6, continuously zero at s=0 |
| K1 | -s/2 * [1-(s/x)^2] * g6 |
| P1 | g4 |
| d_x F1 | (s/x)^3 * g5 |
| d_x Uz1 | (s/x) * g6 |

Core integration uses 16 entire x cells. Actual collar integration uses four entire phase cells and

```text
dx = sqrt(4 exp(hb*s))/2 * hb * ds.
```

The hb factor is multiplied into the source polynomial before enclosure, cancelling every inverse-width term. Remaining positive width factors are retained. Every cell encloses the entire integrand; no midpoint quadrature or finite-axis polynomial supplies field values. Exact zero phase uses the formal core join rho=4, without rounding a positive collar radius to 4.

Published samples use Z=371/1000: rho=4 and s=1/2. There are 46 combined source sectors across the twelve field components (ten core sectors and 36 core-plus-collar sectors), each with Z value and first derivative. Production integration on the already prepared live source took about 1.3 and 2.1 seconds respectively. These timings exclude construction of the inherited source graph and are not a whole-project estimate.

The result contains **combined source and cell-enclosure uncertainty**. It does not yet separate a sharp quadrature error from source uncertainty, certify uniform complex-domain approximation error, or bound the full coupled field error. B0 minus Akin, B1 d_Z W and pressure backreaction through the unknown F1 remain to be iterated. The previous 64-term tail estimate is not the error of this initial field.

## Validation and live use

Focused independent checks cover 228 direct paper forcing/matrix value/derivative rows, three exact bridge/mean derivative identities, actual mean/Q consistency, three independently integrated kinematic kernels, 24 exact core/collar polynomial integration rows, and 36 logarithmic-kernel bounds. All six axis traces vanish exactly. The source-integral core exit and exact first inlet coincide. Owner, chart, derivative, width, cached-amplitude-frame and live packet mutation guards pass.

The interval consistency checks are not a proof that the exact point moment history was numerically reconstructed. One reused read-only worker, **GPT-5.6 Luna / max**, found no definite algebraic error in the bridge rho rows, mean recovery or pressure normalization. No new subagent or Astra child was spawned.

```python
from lei_ren_part1_paper_compliant_current_original_n1_inner_operator import CurrentOriginalN1InnerOperator
from lei_ren_part1_paper_compliant_current_original_n1_initial_integral import CurrentOriginalN1InitialIntegral

inputs = CurrentOriginalN1InnerOperator(accepted_actual_domain_owner)
initial = CurrentOriginalN1InitialIntegral(inputs)
packet = initial.evaluate('.371', '.5', 'first')
view = initial.report(packet)
assert view['current_original_n1_initial_regular_source_integrals_computed']
assert not view['n1_full_inner_interval_solution_certified']
```

## Next executable tasks

- [x] **CURRENT-RP-N1-ACTUAL-COLLAR-OPERATOR-INPUTS** Bind actual Phi/Uz, own mean, independent pressure, mixed rho rows, genuine N1theta/N1z/N1p and B0/B1/g. Preserve ordinary amplitude derivatives and formal source sectors. Reject comparison-core continuation as actual post-core input.
- [x] **CURRENT-RP-N1-INITIAL-REGULAR-INTEGRAL-ENCLOSURES** Compute the six Wbase components and their first Z derivatives with entire-cell integration, exact nilpotent kinematic feedback, zero-axis traces and core/first-inlet identity. This completes the initial integral adapter only.
- [ ] **CURRENT-RP-N1-COMPLETED-LEADING-SUPPORT-GUARD** Bind the actual completed-leading dispatch, Rm and r_minus to this exact source family. Prove Rin below every later modification support using strict formal radii. Keep its current false gate until the support relation is established.
- [ ] **CURRENT-RP-N1-SOURCE-VS-INTEGRATION-ERROR** Retain an explicit ledger separating inherited source uncertainty from the integration enclosure over each cell. Add controlled refinement or analytic quadrature with radial derivative/error bounds. Demonstrate useful normalized widths at core/first-collar endpoints without changing hb, Lambda, P0 or Cstar.
- [ ] **CURRENT-RP-N1-HIGHER-Z-INITIAL-FUNCTION** Extend the mixed source jet budgets or use the accepted holomorphic representation/Cauchy margins so successive B1 applications have genuine derivatives. Provide a callable initial function on the common Z domain, not just the two published samples. Missing higher jets may be bounded explicitly; they must not be replaced by zero.
- [ ] **CURRENT-RP-N1-CONTROLLED-ITERATION-PREFIX** Compute the coupled integral terms from actual operator inputs. Preserve B0 feedback, B1 d_Z W, and the original 2F0F1 pressure term. The exact kinematic resummation may precondition the solve, but it needs its own compatible error majorant. Keep all retained terms and their error budgets. Never call a term zero because its scale underflows.
- [ ] **CURRENT-RP-N1-FULL-FIELD-ERROR** Combine a genuinely computed prefix, source/integration uncertainty, derivative loss and a compatible remaining tail. The preceding 64-term unpreconditioned tail cannot be attached to Wbase as though the retained 64 terms had been computed. Validate the original n=1 equations and recovery on one actual common interval before promoting coefficient-solved/full-inner gates.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Compute the solved order-one field's five whole-Z terminal defects; apply its own original bump controls and pressure-compatible repair. Truncate streamfunction/vector potential before curl. Deliver derivatives, joins and finite-energy tail checks. Leading-order repair is insufficient.
- [ ] **CURRENT-RP-N2-AFTER-REPAIRED-N1** Build actual temporal n=2 products and n-dependent recovery from solved/repaired lower orders on the same interval, then solve and repair it independently. Radial Taylor degree two is not temporal n=2.

Regional stress cone and flatness, measured core width and material winding, finite energy, global/axis physical coverage and oscillatory stress cancellation remain required long-term tasks.
