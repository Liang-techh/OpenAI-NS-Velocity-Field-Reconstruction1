# Actual micro physical U/Q increment bounds

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The
physical finite-N increments delta_U through ordinary Z3 and delta_Q through
ordinary Z2 now have directed error caps on both original micro charts for
Z in [37/100,93/250]. These caps are below 2^-1024. The total-y Q increment
and its R-weighted endpoint expression have C1 norm caps, but the current
majorants do not certify their smallness. All such bounds are retained.
The threshold applies only to these local finite-N increments, not to the
leading-background NS residual.

## Same-source recovery before taking norms

Write B0 for the C5 source B_over_Pstar. The saved mean packet already
contains physical Pstar*B0 caps and the full physical delta_m, including
both C/Cy histories and the second-chart incoming memory. This stage uses

```
delta_U = (Pstar*B0)/N
delta_m = Pstar/N^2*(C-JC_joint-JCy_joint)
L = 1-delta*Z^2 > 0
Q[U,m] = (2Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z)/L.
```

The Q operator is exactly linear. Thus delta_Q=Q[delta_U,delta_m] as a
function before taking absolute caps; no subtraction of two enormous
interval boxes is used. Pstar is applied once. The actual delta and positive
L enclosures remain source-bound; no zero-delta substitution is made.

The reconstructed mean and the original C5 density history solve the same
physical equation delta_m_y+delta_m=Pstar*B0/N, with the same source-owned
zero left inlet and nonreset second history. The C5 mean density rate is
exactly 1, checked from its AST. Linear ODE uniqueness identifies those mean
functions through ordinary Z3. Applying the same U/m operator therefore
identifies recovered Q with the original Q4 function. This is a source
replacement, not an additional ledger contribution or interval-overlap proof.

For h=1-delta and d=1-Z^2, the ordinary derivative recurrence is

```
n0 = 2Z*delta_U - hZ*delta_m - d*delta_m_Z
n1 = 2delta_U + 2Z*delta_U_Z - h*delta_m
     + (1+delta)Z*delta_m_Z - d*delta_m_ZZ
n2 = 4delta_U_Z + 2Z*delta_U_ZZ + 2delta*delta_m_Z
     + (3+delta)Z*delta_m_ZZ - d*delta_m_ZZZ
q0 = n0/L
q1 = (n1-L_Z*q0)/L
q2 = (n2-2L_Z*q1-L_ZZ*q0)/L.
```

Derivative inputs have exactly the required order: U0..2 and m0..3 for
Q0..2. No missing top row is padded by zero and no factorial normalization
is inserted. This is ordinary axial differentiation, not temporal
n-dependent coefficient recursion.

## Total-y chain and the retained outer bound

The actual source defines m_y=U-m and includes the fast phase in U_y:

```
delta_m_y = delta_U-delta_m
delta_U_y = Pstar*(slow_B_y/N+B0_phi)
delta_Q_y = Q[delta_U_y,delta_m_y]
delta(RQ_y) = R*delta_Q_y.
```

B0_phi has no further N divisor or multiplier: its N chain cancels the
1/N in delta_U. The fixed-phase slow B_y/N term remains separately.
The saved complete total-y mean caps are consumed without recomputing them.
Q_y through Z1 needs U_y0..1 and m_y0..2. The original physical R is applied
once using its actual shared-radius upper bound, with no chosen fast phase.
Pinned core recipes and the identical LOOP/Q4 graph prefix establish
R=Ra*exp(hb*s)=epsilon*rho, Ra=4*epsilon, epsilon*Lambda=1, and Z-independent
R on both charts. Q4 RQ_y source nodes are checked to multiply R once.

## Certified error caps

Write Csel for the selected positive logC. The new strict envelopes are
log||delta_U_Zk||<-600*Csel for k=0..3 and
log||delta_Q_Zk||<-2600*Csel for k=0..2. Each also lies below log(2^-1024).
The table gives rounded displays of the directed log-cap upper endpoint
divided by Csel; exact tuples in the report are authoritative.

| ordinary Z order | first delta_Q log-cap/Csel | second delta_Q log-cap/Csel |
|---|---:|---:|
| 0 | -5877.62334199342517903509 | -5977.62334199342517903509 |
| 1 | -4669.62334199342517903509 | -4968.62334199342517903509 |
| 2 | -2649.62334199342517903509 | -3147.62334199342517903509 |

The R-weighted Q_y bounds currently retain the following majorants:

| ordinary Z order | first Rdelta_Qy log-cap/Csel | second Rdelta_Qy log-cap/Csel |
|---|---:|---:|
| 0 | 354.0 | 254.0 |
| 1 | 2374.0 | 2075.0 |

These positive logarithmic majorants fail the local small-error check.
They are not lower bounds on actual error and do not prove that an actual
field or NS residual is large. No cap, midpoint or interval overlap is used
as a signed function value. No outer term is deleted, and the current
stage does not alter the signed ledger or its existing endpoint contractions.

All old 212 forward terms, six signed inner/lower sources and error packets
remain unchanged. Complete signed Q values, the six I1 C1 defect bounds,
outer-endpoint closure, terminal reserved-pulse functional matching and
global five-moment closure remain open. Heat/energy, stress/flat remainder,
true temporal recursion/summation, realizable stress corrections and
independent corrected Cartesian NS validation also remain incomplete.

## Execution and next work

New production exits 0 in 3.565 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_Q_increment_bounds.json` is
17,890,398 bytes, SHA256 `9c41dd5f3bf0dccb3e65582bc2452ede22bfe50ce223dadc57c6eccdf40e823c`, with 2004 input and
373 implementation pins. Saved-only inspection,
same-source AST/ODE identification, nonunit quotient and linearity controls,
directed bounds and independent bounded Luna/max math review pass. Separate
runtime model metadata is unavailable. A pre-acceptance writer gate rejected
a malformed singleton AST constraint; it was corrected to a tuple and
independently reviewed before any report was accepted. Numerical formulas
and enclosures were unchanged. No accepted producer, owner, norm,
recurrence, quadrature or seam checker was replayed.

Implementation: [U/Q error provider](../src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py),
[runner](../scripts/run_original_micro_Q_increment_bounds.py).

```
python scripts/run_original_micro_Q_increment_bounds.py --inspect-report
```

NEXT bind actual background U/Q norms and apply these source-linear errors
to the six I1 C1 defects. Keep the external RQy term until a same-source
endpoint identity or a certified correlated endpoint estimate resolves it.
Then combine the signed micro integral with the reserved-pulse terminal
functional comparison. Goal ACTIVE.
