# Direct original q squared jets before continuous phase transport

Successor: [CURRENT_NATIVE_PERIODIC_C1_SUPPORT_2026_10_07.md](CURRENT_NATIVE_PERIODIC_C1_SUPPORT_2026_10_07.md) derives original whole-period Z support before nonlinear density, tightening80 source derivative ranges and all five terminal Z target upper ranges. Actual controls/global N/closure and recursion remain open.

Checked source: [20afd757](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/20afd757cc0f233ccad160d7178efc5a73bde1d2). Same original source family, full Z[-1,1], candidate N2048 and24 exact continuous radial cells as the preceding [primitive-support stage](CURRENT_NATIVE_PRIMITIVE_SUPPORT_2026_10_07.md). No strict final target upper reduction was established in this run. Rows without a strict final upper reduction: M, I, S, Cp, D=(J-M)/mu, M_Z, I_Z, S_Z, Cp_Z, D=(J-M)/mu_Z.

## Original formulas

Write R=q², Delta=kappa-2 and theta=Delta/eta. Eta has no slow derivatives. On Delta<=0, sigma=1 and R=(2eta-Delta)/(2a); every nonzero numerator derivative is -Delta_alpha. On Delta>=eta, all six R jets are exactly zero.

On the transition branch, let s=sigma(1-theta), H(theta)=s²(2-theta), N=eta H and R=N/(2a). H_i denotes the ordinary ith derivative with respect to theta. The original numerator rows are:

```text
N_0    = eta*H
N_y    = H1*Delta_y
N_Z    = H1*Delta_Z
N_yy   = H1*Delta_yy + H2*Delta_y²/eta
N_yZ   = H1*Delta_yZ + H2*Delta_y*Delta_Z/eta
N_yyZ  = H1*Delta_yyZ
         + H2*(Delta_yy*Delta_Z+2*Delta_y*Delta_yZ)/eta
         + H3*Delta_y²*Delta_Z/eta²
```

The ordinary positive-denominator quotient recurrence supplies all six R rows. It uses only the original positive a lower proof, with no positive lower bound for R at its vanishing endpoint. Sigma Taylor coefficients are multiplied by n! before constructing H_i. In particular,

```text
R_Z = -(2*s*sigma_prime*(2-theta)+s²)*Delta_Z/(2a) - R*a_Z/a
q*q_Z = R_Z/2
```

This collects the exact first-order eta/eta cancellation before directed log arithmetic. The second/third inverse eta factors in the displayed formulas are genuine and remain. A huge-log synthetic fixture checks that R_Z stays bounded independently of an artificial exp(-eta_log), while tiny positive R does not become exact zero. Synthetic fixtures are not original field data.

## Implementation and source identity

Files use `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_collected_q2_transport.py`: `conditional_q2_jet`, `NativeCollectedQ2Cover`, original five-site first-jet adapter, source-local density/oracle and complete transport.
- `current_native_collected_q2_transport.json`: actual original24-cell route, branch and adapter provenance, all ten terminal target range comparisons.
- `current_native_collected_q2_transport_check.py` / `_check.json`: independent original mixed derivatives and phase references, focused actual route/ledger/memory acceptance.

`NativeCollectedQ2Cover` leaves original q and its six derivative objects intact. It attaches separate q2_rows to each original cutoff branch, with the same source roots/context/basis/ledger. Conditional C0 coordinates are outer ranges, not selected source values. Nonlinear phase/density is evaluated before any cutoff or signed-u hull.

`compile_collected_first_jets` checks exact AST replacement counts [1,1,1,5]: the same signed-u constructor, the two accepted collected h factors and the five `q.square()` sites. Only the last five sites are newly adapted. Linear q/q_Z, q-derived u/h, the monotone inverse and implicit derivative rules are unchanged. C0 A/B support is still applied before the nonlinear density, keeping any tighter original formal range. E/E_Z, V/V_Z, P0/P0_Z and all source cross terms remain.

The new companion consumes accepted source receipts; earlier accepted modules and data are unchanged. Candidate N2048 is not a globally compatible N. These are source and integral range coordinates; neither bound endpoints nor midpoint samples define the field or a control.

## Actual target comparison

Strict reductions against the checked C0-support baseline:

| Rc target row | Strict absolute upper reduction |
|---|---|
| M | no |
| I | no |
| S | no |
| Cp | no |
| D=(J-M)/mu | no |
| M_Z | no |
| I_Z | no |
| S_Z | no |
| Cp_Z | no |
| D=(J-M)/mu_Z | no |

The first-order eta cancellation is established locally even where a final target upper bound remains dominated by another unchanged term. Final bounds remain insufficient for a useful five-control repair contraction. Do not convert the number of completed source checks into a percentage of completed scale recursion.

Focused acceptance:42 independent original q² derivative comparisons across negative/zero/transition/eta-endpoint/flat cases;54 original scalar C0/y/Z/phi comparisons with unchanged GenericShearLoop references and numerical allowances;210 actual branch mixed rows;69 five-site original first-jet adapters;240 live integral C0/Z rows. All arithmetic retains the original context/ledger and all quiet pressure memory. Producer 94.937s / checker 5.594s. Read-only review: **GPT-5.6 Luna / max**, scoped PASS.

## Executable next-agent tasks

- [x] **SOURCE3c1-q²-Z:** direct six mixed q² rows, original five-site phase consumer and complete24-cell continuous transport. Source precision is enforced for independent helper checks as well as the producer/checker.
- [ ] **SOURCE3c1-C1-support:** consume the accepted `serial.whole_period_C1`/`periodic_parameter_theorem` original identities. Rewrite nu_Z and the derivative of q² terms using R_Z; derive the corresponding source-function outer A_Z/B_Z ranges. Use branch-local original a/Delta derivatives and same family/hash/context. Retain an existing tighter derivative object. Compare target ranges after actual transport.
- [ ] **SOURCE3c1-paired-q/u:** use `(q*P)_Z=q_Z*(P+u*P_u)+q²*c_Z*P_u` and `(q²*H)_Z=q*q_Z*(2H+u*H_u)+q³*c_Z*H_u`, where u=c*q, c=p2/dstar, P=h^-1*W1 and H=h^-2*W2. Derive original signed/small Poisson factor bounds before any hull. Preserve genuinely linear q_Z terms and report them separately.
- [ ] **SOURCE3c1-phase-C0-conditioning:** if direct q² C0 is used in the inverse normalizations, adapt the original normalizations and phase numerator together, prove original monotonicity/source equivalence and retain the original field. A smaller interval used in only a subset of an identity is not a new source definition.
- [ ] **SOURCE3c2-joint-source:** construct full `Dk*Mk*f_k-Ac*Dm*Mm*f_m` and its first-Z derivative before separate key hulls. m has rate1 and k rate3/2; preserve distinct Duhamel masses, Ac_Z, C_Z and division by mu*Ac². One shared mass is invalid.
- [ ] **SOURCE3c3-refinement:** identify the currently dominant cell from actual weighted target contributions; refine only that source region, preserving exact selected sc/native endpoints, original global radius phase, microscopic widths and incoming memory. Continuous covers, not point samples, define source/integral range bounds.
- [ ] **SOURCE3d-compatible-N:** retain original N^-1/N^-2 coefficient functions before evaluation. Apply same-source support bounds through an all-N route; prove finite log N satisfies original repair inequalities. A fixed-N2048 report cannot be rescaled into new coefficients or choose controls.
- [ ] **CONTROL1a/1b/2:** exact original bump/Gram/quadratic/inverse/Picard function map and first-Z chain, actual convergent C1 controls, rebuilt velocity and whole-Z terminal moment identities. Preserve `B_exact*h+N*r+Q_exact/N=0`; source caps remain separate from defining functions.
- [ ] **REC/HIGH/OUTER/WAVE/PHYS:** genuine n=1 and n>=2 n-dependent recovery, common core interval and independent per-order repair, higher jets/joins and full cones/pressure/heat/energy, smooth summation/flat remainder, two pulse families/averaged stress cancellation, corrected Cartesian NS and measured material winding/scale recursion.

Gate: `current_original_collected_q_squared_ordinary_jets_and_fixed_N_transport_executed`. Actual control, terminal closure, compatible global N, coefficient recursion, oscillatory correction and corrected NS gates remain false. Keep the long-term goal active.
