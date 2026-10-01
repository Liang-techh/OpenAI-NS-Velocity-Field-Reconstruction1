# F28: shared fixed coefficients and whole-axis analytic exit

The same increased-Cstar analytic family selected in F27 now has a fixed,
input-independent Section 9.4 coefficient ledger and an analytic relaxed-cone
certificate on Z in [-1,1], Ra<R<=110. The inner collar is admissible; the
remaining exit is relaxed. This is an exact implicit construction with
analytic bounds, not a whole-axis finite velocity evaluator.

## Reproduce

From the repository root, with Python 3.11 and mpmath:

```powershell
python experiments/root_st073/lei_ren_part1_paper_shared_fixed_step_bound.py
python experiments/root_st073/lei_ren_part1_paper_shared_K1_ledger.py
python experiments/root_st073/lei_ren_part1_paper_shared_global_exit_certificate.py
```

Each producer saves its JSON receipt and binds the complete inherited input
hashes. The fixed-coefficient ledger is constructed before input profiles are
read. Existing finite core/exit producers are unchanged.

## Fixed step and coefficient ledger

The existing flat step exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)) satisfies
0<=sigma'<=8, sharply at s=1/2. Let t=|1-2s| and v=t^2. Then

```
sigma'/8 = (1+3v)/(1-v)^3 / cosh^2(8t/(1-v)^2)
cosh^2(x) >= 1+x^2
(1-v)^4+64v-(1+3v)(1-v) = v(58+9v-4v^2+v^3) >=54v
```

This proves the global bound without changing the cutoff or using a sampled
maximum. It prepares the paper's T=400A long reshape; no long reshape is
claimed in this checkpoint.

Norms are unweighted sums of separate derivative suprema. C1 products have
factor 1, C2 factor 2, and C3 factor 4. A Z derivative consumes an order.
The componentwise F derivative bounds K,K^2,2K^3,5K^4 give FV C3<=23K^5
and F^2 C3<=31K^5. Squaring a full C3 norm would incorrectly lose the
required K^7 stress bound. The comparison uses Aop:C3->C2 coefficient 14
and Pop:C3->C2 coefficient 20, not an Aop:C1->C1 assertion.

| Fixed coefficient | Value |
|---|---:|
| Cq, comparison normalized stress C2 <= Cq K^7 | 748136224 |
| Summed moment/pressure C1 difference coefficient on R<=110 | 140000 |
| Absolute moment/pressure C1 coefficient | 250000 |
| Cstress, normalized stress difference <= Cstress K^6 X | 402922648 |
| K1, one plus the sum of the displayed fixed coefficients | 192922259423111113597 |

For X equal to the summed C1 differences of F and V, the six individual
moment/pressure difference coefficients are 12100,110,72600,48840,880,880.
Their sum is 135410, rounded up to 140000. The radius domain includes both
short switches; the earlier proposed 120000 coefficient is not used.

The adopted eta_tol=min(epsilon0/4,e_star/100) is read directly from the
fixed-bump receipt. Its selected minimum branch and exact interval identity
8j=eta_tol are verified. The choice cstar=eta_tol/(24K1) satisfies both
Section 9.4 restrictions. Both h_b and epsilon_b are the identical positive
scalar cstar K^-100. Its logarithmic enclosure is recorded; the extremely
small physical width is not rounded to zero or materialized.

All 13 decreasing-K conditions pass at K>=10^6: comparison and actual log
changes, axial change, comparison stress/H errors, signed gradient, weighted
actual stress error, weak/strong cone errors, axial bridge, initial switch a,
small kappa, and continuation width.

## Actual exit and moments

The comparison freezes logF,V by (9.23), retaining all five moments. The
actual field is integrated by (9.26), followed by the two flat switches and
the constant-power segment with a=4/5,b=0. All actual moments integrate
continuously from the same core data, with P=P0+Mp. No stress or moment is
reset.

The actual error keeps omega=1-chi:

```
|I/F-q| <= K1 K^20(h_b+epsilon_b) omega.
```

This controls the vanishing stress at the flat inner endpoint. The physical
K sum charges Df and Ef, so |q_f|<=K. The comparison gate
K|qbar-q_f|<1/2 implies |qbar|<=K+1/(2K)<=2K, which justifies the weak
cone coefficient 4K^2. For the first short switch, the receipt explicitly
checks e2^2<1/16<3.5<=Dbar rather than assuming that squared-error gate.

Both axial quantities in (9.28) are checked separately. Core coefficient
averaging for Mz/R divides the nth radial coefficient by n+1<=1, so the
positive Xh embedding bounds the inherited mean. Subsequently,

```
Mz(R)/R = (Ra/R) Mz(Ra)/Ra + (1/R) integral_Ra^R V(s) ds.
```

Its weights are positive, Z independent, and sum to one. The same bridge
budget includes both short-switch velocity increments and hence their
exact moment integrals. No artificial reset or separate unaccounted switch
moment is used. Each velocity/mean C2 error is <epsilon0; their sum is
<2epsilon0. These are the inputs to the W estimate in (9.31)-(9.32).

After the second switch a=4/5 exactly. The angular source bound is about
1.7604999996667>7/5 and the derivative at a hypothetical D=3 crossing is
at least 173.04999996667>0. Therefore D>3 persists to R110. A nonempty
admissible collar is Ra<R<=Ra exp(h_b sigma^-1(gamma/(10K^10))).

A read-only Luna/max mathematical review accepted the updated coefficient,
q, squared-error, eta, and convex-mean arguments. The later amplitude budget
is still an analytic bound recorded in this exit receipt; it should be
recomputed when installing the long reshape.

## Scope and next dependency

The whole-axis exit relaxed cone and nonempty inner admissible collar are
certified analytically. Full Section 9 admission remains false: a corrected
outer/heat profile at the selected Rref is not assembled. Functional five
terminal moment identities, Section 11 admissible stress lifting, flat
remainder, and genuine temporal coefficient recursion remain unfinished.
The completed 144 radial core orders are spatial series, not time orders.

Next: recompute the R110 amplitude budget, install the same-family (9.30)
long reshape with derivative bound 8 and T from the admitted A bound, then
the (9.38) axial restoration. Use exact log-radius offsets -8,-7,-6,-5 from
logRref; finite offsets must not be lost by subtracting them from enormous
floating-point log radii. Preserve the inherited five moments and pressure.
After that, derive the actual normalized defect test and perform Section 10
functional repairs before attaching the corrected exterior.
