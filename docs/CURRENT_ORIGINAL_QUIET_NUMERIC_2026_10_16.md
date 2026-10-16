# Original quiet continuation from2Rc toRp - 2026-10-16

`OriginalQuietVelocity.evaluate(Z,q,log_tau="-2",azimuth=".25")` now continues the original physical velocity and pressure through the whole quiet segment toRp. q is the exact decimal fraction of the original quiet log length. It consumes the actual repaired2Rc inlet and the original positive mu, preserves all nonzero leading histories and P0, and retains the complete incoming n1 pressure plus its fixed quiet increment. The n1 velocity and own mass vanish on this tail by the checked original source identity. The full project remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This is a regional numerical value/derivative enclosure, not genuine temporal n>=2 recursion or a full corrected NS field.

## Original segment and unchanged inlet

```
s = q*Lq, 0<=q<=1
Lq = -60*log(mu)-2-log2
R = 2Rc*exp(s)
Rp = Rc*exp(-60*log(mu)-2).
```

`OriginalQuietNumeric` loads the accepted original compact velocity/pressure report at x=2. It verifies the source family, implementation and scientific input hashes, then consumes the checked original quiet flow, true quiet C5/pressure-C3/K source bindings, signed quiet-pressure source and repaired own-tail identity. Previously checked source digests are reused within this constructor rather than repeatedly rereading the same scientific ancestor bytes. A new Z query computes its actual2Rc inlet with the existing compact implementation.

The exact original selected N is retained. Only the relative finite-N corrections vanish at2Rc by the original compact theorem. Full leading m/h/k/e/p, angular velocity, radial velocity and the independent pressure datum remain. The quiet length and radius use absolute Rc and the original mu; no diagnostic scale or alternative profile is substituted.

## Exact own-rate flow and nonzero mu

Let A and H2=(m2,h2,k2,e2,p2) be the actual leading fields at2Rc. The consumed source equations are

```
E = A*exp(-(1/2+mu)*s), V=0
m = m2*exp(-s), k=k2*exp(-3s/2)
h = h2*exp(-3s/2)+A*exp(-s/2)*(exp(-mu*s)-exp(-s))/(1-mu)
e = e2*exp(-s)-A^2*exp(-s)/2 * integral_0^s exp(-2mu*u)du
p = p2+A^2/2 * integral_0^s exp(-(1+2mu)*u)du.
```

Directly adding mu to1/2 at ordinary precision would lose its contribution. The implementation instead exports each mu correction separately from the reference frame, with the original log(mu) as a separate physical atom. The reference frame is an algebraic decomposition of the same original field, not a selected mu=0 source. The actual inlet already contains its full previous nonzero-mu effects and errors.

Define positive own-rate moments

```
K_n(s) = integral_0^s exp(-(s-u))*u^n du
B_n(s) = integral_0^s exp(-u)*u^n du.
```

The normalized quiet corrections, with one outside factor mu, are

```
U = expm1(-mu*s)/mu = -s*exprel(-mu*s)
Ch = -K1+mu*K2/2-mu^2*K3/6 + [0,mu^3*K4/24]
Ce = -s^2+(2/3)*mu*s^3-(1/3)*mu^2*s^4 + [0,(2/15)*mu^3*s^5]
Cp = -2*B1+2*mu*B2-(4/3)*mu^2*B3 + [0,(2/3)*mu^3*B4].
```

Each interval remainder follows Taylor-Lagrange for exp(-mu*u) or exp(-2mu*u), integrated against a positive kernel; the omitted exponential factor lies in[0,1]. The source derivatives are not derivatives of a cap. Combining them with the complete original inlet gives

```
E = A*exp(-s/2) + mu*A*exp(-s/2)*U
h = h2*exp(-3s/2)+A*exp(-s/2)*(1-exp(-s)) + mu*A*exp(-s/2)*Ch
e = e2*exp(-s)-A^2*exp(-s)*s/2 - mu*A^2*exp(-s)*Ce/2
p = p2+A^2*(1-exp(-s))/2 + mu*A^2*Cp/2.
```

Ordinary Z Taylor rows through2 are transported from the accepted inlet, and products use the same Taylor algebra. The physical swirl and leading pressure preserve their separate mu terms even when addition to the base would round away the correction.

## Positive short kernels and stable long kernels

For0<=s<=4, the positive convergent series are

```
K_n = exp(-s)*s^(n+1)*sum_j s^j/(j!*(n+j+1))
B_n = exp(-s)*s^(n+1)*sum_j n!*s^j/(n+j+1)!.
```

They are different series. The B_n formula follows expansion of exp(s-u) inside exp(-s)*integral exp(s-u)*u^n du and a beta integral. The reader sums128 terms and includes a directed geometric upper bound for the remaining positive tail. Tiny s therefore retains its positive s^(n+1) scale without subtracting nearly equal terms.

For s>4, exact finite recurrences give K0=1-exp(-s), K_n=s^n-n*K_(n-1), and B_n=n!-exp(-s)*sum_(j=0)^n(n!/j!)*s^j. These do not extend an exponential power series over the long quiet length. The full original Lq is about5.65e19; no physical hyper-exponential radius or amplitude is materialized.

## Own pressure survives the zero-velocity tail

The same original mass identity gives q0=-K/R, with K exactly Z independent and k=K/Rc from the checked original turnoff kernel. On the quiet segment

```
N1p = K^2/(4R^2)
P1(s) = P1(2Rc)+K^2/4*(1/(2Rc)-1/(2Rc*exp(s)))
      = P1(2Rc)+Rc*k^2/8*(1-exp(-s)).
```

The increment is exactly zero at s=0, positive for s>0, and has exactly zero ordinary Z derivative. The factor1/8 follows the lower radius2Rc. The inherited core/lower/active/compact/swirl pressure terms are copied without reset or renormalization. Their current incoming function binding is explicitly an evaluation at the compact exit; the quiet pressure recipe and scope are recorded separately.

The n1 velocity/mass tail is zero from2Rc by the same accepted own-field source theorem. Own pressure is certified zero only beyondRv. It is not set to zero at2Rc orRp. The physical P1 weight remains lambda^-2, with no extra Lambda^-2 or sqrt(2R).

## Physical field and pulse input

The physical packet exports cylindrical and Cartesian velocity, streamfunction, radius/time coordinates and pressure through temporal orders0/1. The leading radial field and cumulative mass remain nonzero away from the axial symmetry zero; the local axial velocity is the original exact quiet zero. One shared meridional streamfunction supplies the dependent recovery. Finite differences of independently summarized centers are not certified physical-gradient or divergence measurements.

The native pulse input frame is recovered from these same current fields: u=E, m1=m/E, m2=k/E^2, X=h/E, energy=e/E^2, Mp=p, pressure=P0+p. There is no extra Pstar division. Exporting these enclosures prepares the next layer; the selected original pulse numerical join remains open.

## Focused evidence

The focused run exits0 in74.501s and saved-report inspection exits0. It exports21 original physical packets at3 axial positions over q=0,1e-160,1e-20,1/4,1/2,1, including source-scaled time atRp. Independent115-digit defining-integral comparisons cover128 positive short/long moments and56 nonzero-mu divided integrals against85-digit directed bounds. Other checks include117 actual2Rc numerical seam conditions,78 nonzero-mu sector conditions,144 positive scalar conditions,372 retained-history/same-N conditions,90 incoming/integrated-P1 conditions,144 physical-order/unit/source-zero-tail conditions,45 native pulse-frame unit conditions,3 originalRp radius-recipe conditions,30 scaled-time conditions and5 rejected requests.

The first independent check used a tiny unnormalized integral of order1e-240, below the adaptive integrator's absolute target. Normalizing its independent quadrature to[0,1] resolved the verifier issue before the actual source run. A first complete run passed in74.922s; the final run clarified pressure binding/scope metadata and checked theRp frame/radius before publication. Checked ancestors were not rerun.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_quiet_numeric.json`; 10,954,566 bytes; SHA256 `8aec66cb2a84a386fd52eb472028479443ed4d5c2811033d9a219e91a48ceb4b`; 1839 scientific input pins and76 implementation pins.

The existing read-only worker is `/root/n1_atlas_map`, GPT-5.6 Luna/max; no new agent was spawned. The root retained the original quiet endpoint, pressure coefficient and positive-kernel identities, with independent defining-integral comparisons. The positive B_n series was disputed during review and retained after an explicit beta-integral derivation; the reviewer retracted that concern. The implemented nonzero-mu expansions, physical units and source-tail ownership have no remaining material issue from that review. No accepted scientific ancestor or implementation bytes were modified.

## Next original construction work

1. Connect the independently selected original pulse atRp using this current eight-function input frame, selected U0 and exact physical scales. Check the same-source numerical inlet before extending towardRv.
2. Continue the original pulse pressure and other own histories toRv with its native implicit functions and directed errors. Do not use a free velocity interpolation, zero incoming pressure or fitted pressure tail.
3. Improve signed coefficient evaluation and pressure cancellation with a backward representation only after consuming the complete terminal theorem and all remaining forcing. Current factored P1 terms have directed errors but no useful combined relative-cancellation claim.
4. Close global high-order/mixed joins and common-global-N admission, axis/tail energy, compatible analytic preheat/exact heat exterior, admissible stress and independently measured flat remainder.
5. Build actual independent temporal n>=2 coefficient recursion, smooth sum and two-family oscillatory correction, followed by independent corrected Cartesian NS residual validation. These full-goal gates remain open.
