# Actual angular repair: axial jets through order four

The two actual implicit angular correction functions now have directed axial
Taylor jets through order four on every real `Z in [-1,1]`. The source is the
same compliant `.001` family whose five absolute terminal moments are already
closed. The C0 branch is inherited from its proved contraction; higher jets
come from differentiating its exact moment equations.

This completes the angular-coefficient part of the higher-derivative work.
It does not certify the full velocity field's C4 interfaces, stress cone,
physical energy or temporal coefficient recursion.

## Files and reproduction

Under `experiments/root_st073/`, with prefix `lei_ren_part1_paper_`:

- `compliant_angular_high_jets.py/.json`: same-source higher jets, true Gamma
  defect derivative bounds and uniform branch derivative bounds.
- `compliant_angular_high_jets_check.py/.json`: independent equation, branch,
  Gamma tail and flatten derivative checks.
- `compliant_reconstruction.py --stage angularjets`: ordered reproduction.

The pipeline has 45 modules. All previous pressure/core/outer receipts remain
unchanged. The new layer verifies their transitive hashes, the independently
checked angular branch and the absolute five-moment closure prerequisite.

`CompliantAngularHighJets().coefficients(Z)` returns the scaled and physical
jets of both d_j functions, the correlated preheat difference, all three
scaled Gamma future-defect jets, the exact RHS enclosures and Jacobian.
Coefficients are ordinary Taylor coefficients, **derivative divided by n!**.
An interval Z input encloses the jet at every real center in that interval;
it is not a power series about an interval treated as a single center.

## Differentiating the actual branch

The rescaled branch satisfies

```text
p*x+y=b1(Z)
x+q*y+K*scale*(x^2+q*y^2)=b2(Z).
```

The exact Jacobian at the already selected C0 branch is

```text
J = [[p, 1], [1+2*K*scale*x0, q*(1+2*K*scale*y0)]].
```

For ordinary coefficient n, every higher jet uses this same Jacobian:

```text
J*(x_n,y_n) = (
    b1_n,
    b2_n-K*scale*sum_i=1..n-1(x_i*x_(n-i)+q*y_i*y_(n-i))
).
```

No binomial factors appear with factorial-normalized coefficients. No
contraction iterate is differentiated, and no midpoint branch is selected.
The first derivative is intersected with the independently admitted C1
enclosure of the same source. At Z=0, odd coefficients1 and3 remain exactly
zero by the evenness of the actual functions. Their C0 heat-driven values
remain nonzero; symmetry does not erase the heat correction.

Uniform physical derivative-sum bounds for orders0 through4 are below
`mu/100`. Their natural logarithms are approximately `-2.8246232020442e19`
for the admitted actual parameters. These are enclosures of positive source
scales, not replacements of the functions by zero. This smallness statement
covers axial derivatives of d_j; radial/mixed derivatives of the complete
corrected velocity and their stress-cone contributions remain separate.

## True Gamma derivatives and infinite tails

Write `a=delta/2`, `d=1-Z^2`, `S=1/Rtail>0`, and

```text
H(xi) = E_Gamma(1+a)[(1+xi*v)^(-a)]
xi = 2*d*S*exp(-t)
Dhat = (1-H)/(a*S).
```

The leading term is `2*d*(1+a)*exp(-t)`. Its value error is bounded by
`2*B*d^2*S*exp(-2t)`, and its first d-derivative error by
`4*B*d*S*exp(-2t)`, where `B=(1+a)^2*(2+a)`. For m>=2,

```text
|partial_d^m Dhat / m!|
 <= (1+a)_(m-1)*(1+a)_m*2^m*S^(m-1)*exp(-m*t)/m!.
```

These follow by differentiating the exact positive Gamma expectation and
bounding `(1+xi*v)^(-a-m)` by one. The calculation uses only finite derivative
bounds, and never treats an infinite expansion in S as convergent. Composition
with `delta_d=-2Z*h-h^2` produces finite exponential-polynomial error envelopes
for all requested axial coefficients.

The flat epsilon collar is integrated with directed cells. The pure exterior
is integrated analytically to infinity, coefficient by coefficient:

```text
Theta_hat_tail = integral_3^infinity a*exp((1-a)t)*Dhat dt
Pressure_hat_tail = integral_3^infinity exp(-(1+2a)t)*(Dhat-a*S*Dhat^2/2) dt
Energy_hat_tail = integral_3^infinity exp(-2a*t)*(2*Dhat-a*S*Dhat^2) dt.
```

The leading angular factor uses the exact cancellation `a/(1-(1-a))=1`.
All remaining exponential poles are positive. Ordinary coefficient
convolution bounds the quadratic term. Thus no exterior interval is omitted,
even when delta is too small for an ordinary finite radial cutoff to represent
the true tail. The same exact positive S is kept in function definitions; the
strong cap only encloses its effects.

The preheat difference keeps the previously correlated value and first
derivative. Higher coefficients differentiate the actual flatten integral
through the Q(Z) power. The finite logarithm jet retains every order using
`(log Q)'=Q'/Q`, rather than the older C1-only log helper.

## Evidence and remaining scope

All 15 symbolic coefficient/Gamma/composition identities pass. Independent
fixtures add20 closed-form implicit branch checks,20 true Gamma deficit and
angular/pressure/energy tail derivative checks, and10 direct flatten/log
checks. The flatten reference differentiates the real power analytically
before integration, avoiding expensive nested quadrature differentiation.
The pressure/energy reference tail has an explicit nonzero error bound.

A separate read-only GPT-5.6 Luna/max review confirms the derivative factors,
ordinary-coefficient recurrence, tail poles and energy factor two; it also
checked five Gamma deficit coefficients against independent direct quadrature.

The provider supplies C4 derivatives but **no fifth-derivative Taylor
remainder**. Do not use the fourth-degree jet as a finite-cell approximation
without an appropriate remainder bound. Full outer C4 acceptance, mixed
derivatives/Ur_Z in the pulse and stress-cone acceptance remain false.

Next recover C4 complete future energy, the actual incoming moment/energy
functions, and the selected ap/c1/c2 functions. Preserve their fixed row
normalization factors at every derivative order, signed bump energy changes,
q^2 energy units, separate epsilon atoms and formal nonzero end-energy factors.
Then propagate the jets through the full velocity field and prove interfaces
before accepting the global cone or starting genuine temporal coefficient
recursion.
