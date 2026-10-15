# Full original n1 core axial Z2 and recovered q_Z — 2026-10-15

## Evaluated interface

`OriginalN1CoreAxial2.evaluate(Z,rho)` returns complete-error numerical
core values through native axial Taylor order two for the six natural
state fields f,u,k,p,fx,ux and their own mean. The same mean recovers q
through order one, including q_Z. The original domain is `|Z|<1, 0<=rho<=4`.

```python
from openai_ns_reconstruction.original_n1_core_axial2 import OriginalN1CoreAxial2
field = OriginalN1CoreAxial2()
packet = field.evaluate('.3', '4', degree=20, kernel_degree=28)
f_jet = packet['normalized_original_n1']['f']
q_value, q_Z = packet['normalized_original_n1']['q']
```

All returned arrays are Taylor coefficients: entry k is `D_Z^k/k!`.
Multiply entry 2 and its absolute error by **2** for the ordinary second
derivative. The six state and mean arrays have three entries; q has two.
This is a separate API; the existing accepted Z1 APIs and receipts remain
unchanged. In particular, the older mean_u scalar view is not replaced.

The natural frame remains
`f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon, p=P1/(epsilon*Pstar^2)`,
with `fx=D_x f`, `ux=D_x u`, `x=sqrt(rho)`.

## Principal center and full feedback enclosure

The original forcing polynomial already has genuine axial coefficients
0..2, obtained from the coupled leading recurrence through Z4. The new
path keeps those orders during principal Volterra integration. The finite
source polynomial is integrated analytically with the Kchi inverse; it
is not treated as a radial Taylor truncation of the unknown solution.

The old nominal feedback step consumes one axial order. No Z2 feedback
iterate is supplied without the missing Z3 principal data. Instead the
new numerical center is the principal response, and **all** feedback words
are bounded from that same principal origin by the accepted index-zero
regular-series receipt. Feedback is not set to zero in the original field.

The error ledger separates:

1. true-leading-source versus finite forcing-polynomial tail;
2. finite versus infinite principal Kchi kernel tail;
3. infinite feedback-word tail relative to the exact principal response.

The accepted ordinary Z2 feedback remainder is divided by 2 for Taylor
units, then multiplied by rho for f,u,k,p and sqrt(rho) for fx,ux. Its
same-tube Cauchy proof is inherited and hash-pinned; it is not reused as
a radial Taylor remainder or as a principal-kernel remainder.

## Genuine Z2 source-error propagation

The mixed source-difference map takes leading radial Taylor order two
and axial order four. Each weighted T operator consumes one radial and
one axial order. Q is kept internally through Z3 so Q_Z is present
through Z2. The original nonzero implicit S and all pressure atoms remain
inside the directed coefficient/source errors. Missing mixed rows are
rejected; no unknown coefficient is padded with zero.

For the Bessel model tail, actual real chi jets through order four are
used. The exact real identity `chi=H^2/(H^2+(j/500)^2)` proves chi0 in
[0,1]; intersecting its interval with this range preserves the same
source. The differentiated factorial-tail bound retains small actual
chi derivatives rather than applying an inverse tube radius to every
model coefficient. The original nonlinear Xh correction still uses the
accepted differentiated tail_factor.

The model axial tail is already in Taylor units. Nonlinear ordinary
axial bounds are divided by k!, and both are divided by i! for the
radial Taylor coefficient i. The axial field tail additionally retains
its original epsilon factor. The own mean is taken from the same axial
polynomial with radial coefficients divided by n+1.

## Principal kernel Z2 tail

Let K be the retained kernel degree, j=K+1, q=rho/2, and

`t_n=q^n/(n!*(n+1)!)`, `T_m=sum_{n>=j} n^m*t_n`.

For m=0,1,2,3, the successive ratio

`r_m(n)=q*(n+1)^(m-1)/(n^m*(n+2))`

decreases with n. Thus `T_m <= j^m*t_j/(1-r_m(j))` when the denominator
is positive. This matches the existing m=0,1,2 tails and supplies m=3.
No extra power of 2 is introduced: q=rho/2 already includes it.

For actual chi Taylor coefficients c1,c2, the kernel coefficient bounds are

```
f:  x/2 * [T0, |c1|*T1, |c2|*T1 + |c1|^2*T2/2]
fx:       [T1, |c1|*T2, |c2|*T2 + |c1|^2*T3/2]
```

Positive convolution against the finite source-size jet includes the
external path factor x exactly once. The independent source-error
response uses `exp(q)` times the multipliers

`[1, q|c1|, q|c2|+(q|c1|)^2/2]`,

with original prefactors rho/8 for f and x/4 for fx. In particular the
Z2 cross term involving the source's Z1 coefficient is retained. Kernel,
source and feedback errors are not interchanged.

## Own mean and q_Z

With m=u+k, the original n1 recovery is

`q=[2Z*u-(1+delta)Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2)`.

The state/mean Z2 data now support q_Z. They do not support q_ZZ, which
would consume m_ZZZ. This is higher axial differentiation of temporal
n=1, not the temporal n=2 coefficient or genuine higher-order recursion.

## Executed evidence

The accepted run exits 0 in **71.283 seconds**. Two levels use
precision/source-degree/kernel-degree 55/14/18 and 85/20/28, with pressure
quadrature sizes 96 and 192. Each computes 12 original core points:
Z=-0.3,0,0.3 and rho=0,1,2,4. It passes:

- 18 independent raw source Z0/Z1/Z2 base-and-difference identities;
- 36 independent positive finite-kernel tail checks;
- 276 complete-error comparisons across precision/source/kernel levels;
- 42 comparisons with the separately invoked accepted Z1 interface;
- the exact axis with all returned Z coefficients zero and five rejected inputs.

At Z=0.3,rho=4, the refined f Taylor coefficient of order2 is approximately
**-8.76374656071860947480**, with error below **4.952e-33**. Thus ordinary
f_ZZ is approximately -17.52749312143721894960, with error below 9.904e-33.
The recovered q_Z is approximately **94.46400000000000000000086518**, with
absolute error below 2.989e-85. These are normalized sample values, not
physical derivatives or uniform relative-error claims.

The largest sampled Z2 coefficient relative errors are below 5.651e-34
for f, 3.135e-33 for fx and 2.720e-27 for p. The report stores the other
state fields separately. The 942309-byte report pins 1478 source inputs
and 11 implementation files; fresh saved-report inspection exits 0.

The existing read-only Luna/max worker independently checks the actual
source-order/factorial propagation and the kernel moment-tail proof. It
also confirms that the feedback tail uses Taylor units and q stops at Z1.

Driver: `scripts/run_original_n1_core_axial2.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_axial2.json`.

## Remaining construction

Tasks 14/15 remain IN PROGRESS; the full goal remains ACTIVE / INCOMPLETE.

1. Transport these genuine Z2 inlet values and their errors through the
   actual admitted first-half collar. Its source/kernel/regular-operator
   bounds must retain the extra required axial orders.
2. Recover collar q_Z, native and physical mixed derivatives, then complete
   quantitative physical joins and vorticity/strain diagnostics.
3. Admit/evaluate remaining collars, apply the source-owned cutoff and own
   five-moment repair, and complete global physical/chart assembly.
4. Complete five-moment numerical closure, heat exterior, finite energy,
   admissible stress/flat remainder, genuine temporal n>=2 coefficients
   and their own repair/summation, oscillatory correction and independent
   corrected Cartesian NS validation.

The new interface is core-only. It does not supply collar Z2 transport,
physical mixed derivatives, a complete Cartesian C2 join, a global velocity
field, or a blow-up theorem.
