# Actual micro q mixed bounds and scaled Poisson conditioning

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The actual
q and scaled Poisson inputs u, hinv, r, alpha and K now have ordinary Z0..3
and physical slow-y0..1 logarithmic norm caps on both original micro charts,
for Z in [37/100,93/250]. The original first chart is s in [sc/2,1] and the
second is s in [1,2]; the source-owned flat left extension is zero. Positive
free-angle denominator and phase-derivative log floors are now quantitative.
The mixed derivatives of the phase/inverse, actual loop primitives and
finite-N derivative errors are still open.

## Source branch and mixed arithmetic

The current C5 source is retained exactly:

```
kappa=a+b^2/a; Delta=kappa-2
gamma=2eta-Delta; xcut=1-Delta/eta
q=sigma_flat(xcut)*sqrt(gamma/(2a)) on the active branch
q=0 on Delta>=eta or the original first-micro collar.
```

On the active branch Delta<eta, so gamma>eta>0, kappa<2+eta<3,
a<=kappa<3 and |b|<=kappa/2<2. These inequalities tighten only the value
rows of the source jets. All higher rows come from the already accepted
actual E/a and b/t0/p2 derivative providers, not derivatives of value caps.
The source-owned flat body is never evaluated through a singular square root.

Mixed jets have eight ordinary rows (y order 0 or1, Z order0..3), including
the total-order-four yZ3 row. Directed log-sum bounds, binomial products,
reciprocal and square-root recurrences, and scalar Taylor composition
propagate norm caps without exponentiating the underlying extreme functions.
Missing rows fail closed; None is reserved for a proved exact-zero derivative.
Exact nonunit SymPy controls validate the ordinary factorial convention.

For reciprocal a, the reference a0 is frozen at the point being bounded.
The local function a(x)/a0 has value1 there and derivatives a_derivative/a0.
The normalized inverse jet is then multiplied by the global positive 1/a0
bound. This is not a global a/a function or a choice of the actual a value.
The first-chart chi_y/hb effects remain in the source mixed input caps.

The original scalar sigma has global derivative caps through order 4:

| derivative order | absolute cap |
|---|---:|
| 0 | 1 |
| 1 | 1296 |
| 2 | 38880 |
| 3 | 19290960 |
| 4 | 2251925712 |

These follow analytically from the logistic derivative polynomials,
sigma(1-x)=1-sigma(x), and the exp(-1/x^2) flat-tail bound. No sampled
supremum is substituted. The scalar derivatives match the constant branches
through order 4; this supplies the shear-flat mixed gluing.

## Natural flat neighborhood of the first collar

The original forced collar boundary is 3sc/4, and accepted sc<2^-10 places
it strictly inside [0,1/512]. On that neighborhood alpha_comparison=1.
FTC identifies the comparison fields with the same core, while complete
core-integral inlets and the same linear history ODEs identify the histories.
The accepted functional first-interface lemma is consumed as a receipt;
its checker, moment producers and core computation are not replayed.

The same signed core bound gives Dbar_core>3.069, the actual rho argument
lies in [4,4.1], and

```
sigma_flat(s)<=exp(4-1/s^2)<1/1000 for 0<s<=1/512
chi_R>=999/1000; a=chi_R*Dbar_core>3
Delta=kappa-2>1>eta.
```

Thus the underlying shear recipe is already exactly flat in a neighborhood
of the forced collar boundary. q and all mixed derivatives through total
order4 are zero there. At s=0 the statement is one-sided relative to the
first-chart source domain; no lower-rho extension is claimed. The second
comparison is never identified with the core. Exact C6 core/inlet prefixes
and current q handles bind this argument to the same original source.

## Actual parameter and derivative bounds

Write C for the selected positive logCstar. The saved exact parameter logs
prove -245C<logeta<-243C, -101C<loghb<-99C and logN>6131C. Ratios near -244,
-100 and 6131.62334 are displays, not exact equalities. Eta, hb, N and the
fast phase remain unmaterialized.

Every q/u/hinv/r/alpha/K mixed log norm is strictly below 2048C. The following
are rounded displays of the directed q log-cap upper endpoint divided by C;
exact tuples in the report are authoritative. These are norm caps, not
actual signed functions or signed correction values.

| y order | Z order | first q log-cap/C display | second q log-cap/C display |
|---|---:|---:|---:|
| 0 | 0 | 50.0 | 50.0 |
| 0 | 1 | 395.0 | 196.0 |
| 1 | 0 | 594.0 | 196.0 |
| 0 | 2 | 768.0 | 416.0 |
| 1 | 1 | 967.0 | 416.0 |
| 0 | 3 | 1213.0 | 660.0 |
| 1 | 2 | 1412.0 | 660.0 |
| 1 | 3 | 1857.0 | 904.0 |

## Strict scaled Poisson and phase conditioning

The actual source uses

```
u=p2*q/dstar; hinv=1/sqrt(1+u^2)
r=u*hinv; alpha=2q*hinv
K=1/[2pi*(1+t0^2+2q^2)].
```

Both scalar functions hinv(u) and r(u) have global ordinary derivative caps
[1,1,3,15,105] through order 4. Differentiating gives a polynomial divided by
(1+u^2)^(k+1/2); each monomial quotient is at most 1, so its coefficient sum
bounds the derivative. This avoids divisions by a rounded near-unit r.

For every real free angle psi in [0,2pi], the exact recipe proves

```
1-r^2=1/(1+u^2)>0
D=1-2r*cos(psi)+r^2 >= (1-|r|)^2 >= (1-r^2)^2/4
lambda=K*(1+t^2)>=K>0.
```

The new strict logarithmic floors on both charts are
log(1-r^2)>-105C, logD>-211C and loglambda>-101C. The report retains the
actual directed floors separately. No actual q,u,r or hinv value is selected,
and r is never replaced by 1. Phase normalization from the accepted C0 source
is retained; these floors do not supply implicit-inverse derivative bounds.

## Execution and next work

New production exits 0 in 6.742 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_q_mixed_bounds.json`
is 17,020,832 bytes, SHA256 `a28ae027174ceefa92c80c7f7437689ae964c83ca5c116c4c088a2ec156d2c39`, with 2001 input and
366 implementation pins. Saved-only inspection,
source/AST checks, ordinary mixed and scalar controls, natural-collar proof,
and directed cap/conditioning checks pass. Fixed Luna/max workers reviewed
bounded math/source tasks; separate runtime model metadata is unavailable.

All old 212 forward terms, six signed inner/lower sources, endpoint identities
and derivative-error packets remain unchanged. Existing dirty files and the
unaccepted scaled-Poisson draft are preserved. No accepted producer, owner,
norm computation, recurrence, quadrature or seam checker was replayed.

Implementation: [q provider](../src/openai_ns_reconstruction/original_micro_q_mixed_bounds.py),
[log mixed arithmetic](../src/openai_ns_reconstruction/original_log_mixed_jets.py),
[flat sigma derivative proof](../src/openai_ns_reconstruction/original_flat_sigma_derivative_bounds.py),
[runner](../scripts/run_original_micro_q_mixed_bounds.py).

Inspect the saved receipt without replay:

```
python scripts/run_original_micro_q_mixed_bounds.py --inspect-report
```

NEXT bound the actual mixed free-angle t/Phi/lambda functions, then propagate
the implicit inverse and compose A/B through physical-y1/Z3. Only afterward
certify finite-N derivative errors, actual Cy_Z3, joint mean-Z3/QZZ and outer
RQy. Independent terminal reserved-pulse matching and global five-moment
closure remain open, along with heat/energy, stress/flat remainder, true
n>=2 recursion/summation, realizable corrections and independent corrected
Cartesian NS validation. Goal ACTIVE.
