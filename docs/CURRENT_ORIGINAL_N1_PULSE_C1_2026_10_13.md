# Original n1 outer-pulse pressure C1 source — 2026-10-13

This extends the accepted Rp-to-Rv pressure source with genuine first-Z derivative enclosures. A separate 4096-cell derivative grid is checked on every real Z in [-1,1], with maximum normalized errors 0.0005117131891907913 for m13 and 0.0001210322320544479 for m15, both below .001. The exact signed linear jets are reused; nonlinear derivative errors are derived from the original C5 coefficient/incoming-mean rows through order three. Full pressure C1 accuracy, actual n1 repaired fields, terminal five-moment closure and temporal recursion remain open.

Read the [physical n1 linear operator checkpoint](CURRENT_ORIGINAL_N1_LINEAR_OPERATOR_2026_10_13.md) for the repair matrices and actual angular amplitude. The outer-pulse derivative data supplies one of that operator's missing input segments; it does not provide the missing pre-Rp derivative integrals or choose the actual defect functions.

## Actual source and differentiated pressure

The producer restores one checked pulse owner using the constructor cache in `run_original_n1_five_defect_values.py`; it does not execute ancestor producers/checkers. The four original whole-Z coefficient packets supply B_j and Tinitial_j for j=0,...,3. Actual ordinary derivatives include their factorial conversion exactly once. The source rows retain ap, both end coefficients, the incoming mean, finite mu and the same original pressure units.

For L=1-delta*Z^2, d=1-Z^2 and a=1-delta, meridional recovery is q=(2ZU-aZm-d*m_Z)/L. The first and second derivatives are taken from this defining expression, not from a magnitude bound. The seven q-Z2 coefficients have uniform majorants [56,24,4,28,38,16,2] on |Z|<=1, 0<=delta<=1/2. From E(U_j)<=epsilon and E(m_j)<=10epsilon, seven-term Cauchy gives E(q_ZZ)<=200256epsilon. The positive epsilon=exp(-1000) is retained, including its directed endpoints.

Differentiating reduced S0/S1 also differentiates every Z-dependent coefficient: 1/L, aZ/(2L), d/L, Z^2/L^2 and Zd/L^2. The original product derivatives include

```text
D(U*q_Z) = U_Z*q_Z + U*q_ZZ
D(U*U_Z) = U_Z^2 + U*U_ZZ
D(q^2) = 2*q*q_Z
```

The energy bound for the second line is E(U_Z)+sqrt(E(U)E(U_ZZ))<=2epsilon. Multiplying scalar energy caps would incorrectly give epsilon squared. The independent checker differentiates the actual integrated pressure formulas and checks finite coefficient envelopes rather than testing a synthetic repeated identity.

The finite weighted viscous endpoint derivative remains signed: -qinit_Z*invRv/reference. The unweighted endpoint cancellation remains exact. Parent C0 value models are retained; first-derivative errors separately identify source, center-slope ambiguity and whole-cell second-derivative Taylor remainder.

## Current numerical accuracy and scope

The inherited 256-cell C0 grid has maximum first-derivative errors approximately .0085353753 for m13 and .0014673022 for m15. Their source errors are only about 5.49e-867 and 1.86e-867. The corresponding maximum center-slope ambiguity contributions are about 6.10e-6 and 6.05e-5. The principal limitation in this segment is therefore the Taylor remainder on the inherited grid, not an unresolved nonlinear pulse source term.

The finer derivative grid and its local 1e-3 accuracy gate require a first-derivative error estimate: for a constant center-slope approximation, the whole-cell Taylor contribution is 2*abs(pressure_Taylor_row_2)*half_width. The value estimate abs(row_2)*half_width^2 cannot be used for this derivative. An independent nonconstant-polynomial fixture checks that distinction. All thresholds use the original positive pressure-moment units. Local accuracy does not imply accuracy of the full Rin-to-infinity derivative: the pre-Rp q-Z2 source remains too wide and finite lower endpoint U_yZ/U_ZZ recovery still needs an adapter.

The installed gates are `current_original_n1_pulse_pressure_C1_source_installed` and `current_original_n1_pulse_pressure_C1_local_grid_accuracy_installed`. The 4096 lean derivative records keep the actual center-slope and whole-subcell row-two enclosures, physical source error, center-slope ambiguity and correct derivative Taylor budget. The independent p(Z)=Z^2 fixture checks the 2h derivative variation at center 1/2 and rejects the h^2 estimate. Root rechecked the corrected budget arithmetic on all 8192 moment rows and the 1307 exact-byte input pins. The earlier invalid derivative estimate is replaced and was never published as a completed gate.

## Reproduction and agent tasks

Run `python scripts/run_original_n1_pulse_pressure_C1.py`. The quartet is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1` with .py, .json.gz, _check.py and _check.json. The original 256 C0 models are preserved and genuine source samples remain separate from any numerical affine representation. No sample/range endpoint is selected as a field or coefficient function.

- [x] Extend actual U and initial-mean derivative source caps through order three.
- [x] Correct the q-Z2 seven-term energy factor to 200256epsilon.
- [x] Propagate differentiated reduced-S coefficient and product bounds in the original units.
- [x] Retain the finite weighted viscous derivative with its original minus sign.
- [x] Check actual C5 source rows, seven samples and inherited 256 models.
- [x] Install a finer separate 4096-cell derivative grid with local normalized absolute C1 error below 1e-3 on all 8192 moment rows.
- [ ] Add the inner cutoff's genuine signed row-one sectors and unchanged lower-collar derivatives.
- [ ] Replace the unusably wide pre-Rp q-Z2 magnitude cap with a correlated signed derivative integral or sharper genuine source estimate.
- [ ] Recover genuine lower endpoint U_yZ and U_ZZ rows; preserve the actual finite endpoints.
- [ ] Compose the complete pressure derivative in fixed original positive units without differentiating C0 affine errors.
- [ ] Bind complete genuine signed defect functions to the actual n1 linear operator.
- [ ] Select actual repair coefficients, recover repaired own M/V/P and prove all five terminal identities.

Full original goal **ACTIVE / INCOMPLETE**. A checked local derivative source is a construction dependency, not a completed velocity reconstruction or blow-up theorem.
