# Actual micro fixed-phase inverse and loop mixed bounds

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. This
milestone supplies same-source fixed-phase inverse and loop A/B logarithmic
norm caps on both original micro charts, with physical slow-y order 0..1 and
ordinary Z order 0..3. It also supplies free-angle jets through total order
four, positive D/lambda conditioning, and the normalized implicit inverse
recurrence. It does not select an unknown phase, signed source value, finite-N
phase, or full C5/global closure.

## Source formulas and fixed-phase scope

For the original free angle psi, with D=1-2*r*cos(psi)+r^2, the source
recipe is retained as

```
r = u/sqrt(1+u^2),  alpha = 2*q/sqrt(1+u^2)
t = t0 + alpha*(cos(psi)-r)/D
lambda = K*(1+t^2) > 0
Phi(psi) = K*(psi + integral_0^psi t(theta)^2 dtheta)
```

The first line is the actual source's scaled Poisson parametrization; the
provider carries the exact source AST and positive directed D/lambda floors.
The free-angle phase is normalized by Phi(0)=0 and Phi(2*pi)=1, and strict
monotonicity supplies one inverse angle for every target phase. No target
phase is chosen for the norm provider.

The inverse recurrence uses ordinary derivatives. Its normalized increment is
a Taylor coefficient: phase terms are divided by parameter and angle
factorials, coefficient products carry no binomial factors, and the returned
ordinary rows restore the parameter factorials. The 24 free-angle rows per
chart are all (physical slow-y order i, ordinary Z order k, free-angle order
j) with i<=1, k<=3 and i+k+j <= 4. Each inverse, A, B, A_phi and B_phi mixed field has the
eight ordinary rows (i,k) with i in {0,1} and k in {0,1,2,3}.

The fixed slow-y derivative is not the total physical-y derivative. The
actual chain is

```
d_total_y = d_slow_y + N*d_phase,    phase_Z = 0.
```

The actual E product rule is retained. E is a positive point-frozen
reference in the cap arithmetic, not a selected value and not a division of
a function by itself. The source conventions are

```
A = a/2*(phase - psi/(2*pi))
B/Pstar = E/2*(-a*T1hat/(2*pi) - b*phase)
A_phi = a/2*(1 - 1/(2*pi*lambda_hat))
(B/Pstar)_phi = E/2*(-a*t_hat/(2*pi*lambda_hat) - b).
```

Here B is physical only after multiplying the source B/Pstar by Pstar. The
provider retains all E cross terms in y/Z products; cap signs are not signed
source values.

## Endpoint Fourier identities and flat branch

With H=(cos(psi)-r)/D,

```
H = sum(n>=1, r^(n-1)*cos(n*psi)),
mean(H) = 0,
mean(H^2) = 1/(2*(1-r^2)),
mean(t) = t0,
mean(t^2) = t0^2 + 2*q^2.
```

These are parameter-function Fourier identities, with accepted 1-r^2>0;
they are not sampled-angle endpoint substitutions. They yield
T1(2*pi)=2*pi*t0, T2(2*pi)=2*pi*(t0^2+2*q^2), Phi(0)=0 and Phi(2*pi)=1.
Together with t0=-b/a they prove A and B vanish at both phase endpoints.

The source-owned flat branch is applied before active inverse/sqrt formulas.
There psi=2*pi*phase, Phi=psi/(2*pi), lambda=1/(2*pi), and the slow inverse
and Z derivatives are exact zeros through the stated order. The proof is
one-sided at the first-chart boundary where required; it does not evaluate a
singular active formula on the flat body.

## Directed caps and limits

Write C for the selected positive logC; the saved parameter records also
prove logN > 6131*C. The strict logarithmic norm-cap envelopes are < 9000*C for all eight mixed loop/inverse fields, including yZ3, and < 5600*C for A and the B/Pstar derivative rows divided by frozen E0 at ordinary Z0..3. These are norm caps; they do not prove smallness of the full physical E or any signed component.

The following entries are read from the report's serialized `log_cap_ratios`
records, which already divide each logarithmic cap by selected positive logC. They are rounded displays only;
the directed records in the report remain authoritative. The y1/Z3 columns
show the highest requested slow-y/ordinary-Z row, while the Z3 columns also
show the y0 ordinary-Z3 row.

| field | first y0/Z0 | first y0/Z3 | first y1/Z3 | second y0/Z0 | second y0/Z3 | second y1/Z3 |
|---|---:|---:|---:|---:|---:|---:|
| inverse | 1.28469806572e-408906090034569674 | 5248.0 | 7467.0 | 1.28469806572e-408906090034569674 | 4651.0 | 6472.0 |
| A | 7.67943138321e-408906090034569675 | 5248.0 | 7467.0 | -100.0 | 4551.0 | 6372.0 |
| B/Pstar rows / E0 | 254.0 | 5502.0 | 7721.0 | 154.0 | 4805.0 | 6626.0 |
| A_phi | 100.0 | 6160.0 | 8379.0 | 1.0428728639e-79 | 5463.0 | 7284.0 |
| B_phi/Pstar rows / E0 | 354.0 | 6414.0 | 8633.0 | 254.0 | 5717.0 | 7538.0 |

The actual source/signed error packets, all old 212 forward terms, six signed
inner/lower sources, endpoint identities and inherited error caps remain
unchanged. This stage does not provide finite-N derivative errors, C_y Z3,
mean Z3, QZZ, outer RQy, terminal matching or global five-moment closure.
NEXT use these saved B and slow By caps in the source-owned C/Cy primitives
and joint mean correction, with Pstar/N^2 and the second-chart incoming
history retained. Then derive finite-N errors with the total-y chain, QZZ,
outer RQy and independent terminal reserved-pulse matching. Separate runtime worker-model
metadata is unavailable and is not claimed.

## Execution receipt

New production exits 0 in 3.994 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_phase_inverse_bounds.json`
is 17,617,354 bytes with SHA256 `0551d3b8acbcc4578cf97591a299e679f3e589efd4ce6d8e8e66b765e2cbabdc`, 2002 input pins and
369 implementation pins. Saved-only inspection and source
provenance checks are retained as receipts; no accepted producer, owner, norm
computation, recurrence, quadrature or seam checker was replayed.

Implementation: [phase-inverse provider](../src/openai_ns_reconstruction/original_micro_phase_inverse_bounds.py), [phase-jet arithmetic](../src/openai_ns_reconstruction/original_log_phase_jets.py), and [runner](../scripts/run_original_micro_phase_inverse_bounds.py).

Inspect the saved receipt without replay:

```
python scripts/run_original_micro_phase_inverse_bounds.py --inspect-report
```

Goal ACTIVE / INCOMPLETE.
