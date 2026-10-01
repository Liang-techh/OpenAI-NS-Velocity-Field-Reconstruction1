# Actual functional angular and pressure repair (F38-C)

The compliant source now has a unique small smooth two-bump angular
correction on the entire axial interval Z in [-1,1]. It restores both the
actual heat angular target and the prescribed analytic preheat pressure
datum. The coefficients, their first axial derivatives and the partial
angular/pressure/energy changes are available as directed enclosures.
Actual axial amplitude selection and full outer five-moment assembly remain
unfinished.

## Exact correlated right sides

The source and five-moment family remain those of F37/F38. All pre-flatten
swirl and cumulative angular moments carry the same (1+Z^2)^-1 factor.
Consequently Xv=Mtheta/(sqrt(2)*Rv^(3/2)*Utheta) is independent of Z,
including its retained pulse memory. At the flatten end, eta=Z^2 gives

```text
Xf(0)-Xf(Z)
 = 2*Xv*exp(-100*rate)*(1-1/(1+eta))
   + integral_0^100 exp(-rate*(100-s))*2^k(s)
       *[1-(1+eta)^(-k(s))] ds,
rate=1-mu, k(s)=1-sigma(s/100).
```

The bracket is evaluated by its exact positive divided difference
`[1-(1+eta)^(-k)]/eta = integral_0^1 k*(1+v*eta)^(-1-k) dv`.
This retains the exact zero at Z=0, evenness and the first derivative;
subtracting independently enlarged X boxes would lose this correlation.
The subsequent power interval multiplies this difference by
`scale=mu^(30*(1-mu))`, which is retained as a nonzero directed number.

The preheat waiting identity at Z=0 makes this correlated difference the
complete preheat angular mismatch. No additional independently fitted
target or pressure datum is introduced. The final right side also includes
the true Gamma heat mismatch, so its value at Z=0 is generally nonzero.

Let L=log(Rtail/Rrel), Et=log(Utail/Erel), a=delta/2, S=1/Rtail and
ell1=log(1-epsilon). The exact heat outputs imply

```text
rH = exp(3L/2+Et-ell1)*S*Theta_hat,
sH = exp(2Et-2ell1+log(a))*S*Pressure_hat.
```

The angular exponent is evaluated in its cancellation-free form
`(rate+1-a)/2+(1-a)*tau-ell1`. A stronger finite log cap for S is proved
before dividing by the coefficient scale. Its upper bound is not used as
the actual S. The formal selected inverse-radius terms and the exact heat
integrals remain part of the coefficient definitions.

## Recovery and small branch

The fixed normalized bump has logarithmic radius .15. Its two translated
supports are centered at -3 and -1 relative to Rrel and are disjoint.
Directed integration computes A, B, D (and energy weights E, F). The actual
Section 7.21 equations are

```text
A*(exp(-3rate)*d1+exp(-rate)*d2) = r,
B*(exp(3prate)*d1+exp(prate)*d2)
  + D/2*(exp(3prate)*d1^2+exp(prate)*d2^2) = sH,
prate=1+2mu.
```

There is no mixed quadratic term because the supports do not overlap.
With p=exp(-2rate), q=exp(-2prate), k=D/(2B), the rescaled linear matrix is
`[[p,1],[1,q]]`. Its determinant is approximately -0.981684361111266 and
its inverse row norm is bounded by 1.15651764274967. The full nonlinear
fixed-point map preserves its whole-Z box and is a contraction.

Current numeric bounds:

- `log(scale)` is approximately -2.8246232020442e19.
- Scaled whole-Z right-side sup: 1.2224746325318e-16.
- Scaled coefficient ball radius: 1.2224746325318e-15.
- Contraction Lipschitz upper is positive and far below .05.
- The relative swirl perturbation plus its radial derivative is below mu/100.

The very small quantities are stored as arbitrary-exponent directed mpf
data; they are not float underflows or zero masks. The coefficients are
implicitly defined actual functions, with enclosures rather than midpoint
substitutions. Their C1 bounds come from differentiating the actual
quadratic system and inverting its uniformly nonsingular Jacobian.

## Callable changes and checks

`CompliantAngularRepair.partial_corrections(Z,offset)` provides corrected
swirl and all five cumulative changes in Rrel/Erel units. Uz, Delta Mz and
Delta Mtheta_z are zero in these supports; Delta Mztheta is minus one half
of the swirl energy change. Angular and pressure primitives include their
partial bump integrals. The fixed bump vanishes to every order at its
support endpoints. These changes still require composition with the actual
axial-pulse field to produce a full outer velocity evaluator.

The independent checker passes twenty-three identities and source invariants,
five independent limiting bump quadratures, two correlated
flatten value/derivative fixtures, the exact axis-even derivatives and the
whole-Z contraction/derivative gates. The quadratures check formulas;
actual-source existence follows from directed bounds, not from fixtures.
The q^-1 factor is checked directly against the u/h source expressions in
all five upstream transport methods. Independent symbolic propagation
across the complete pulse proves Xv is scalar for arbitrary Z; the angular
moment equation contains no Uz and is unchanged by any axial correction.

Reproduce with:

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage angular
```

The standalone modules are `compliant_outer_angular_repair.py` and
`compliant_outer_angular_repair_check.py`, with separate source-bound JSON
receipts and the usual prefix `lei_ren_part1_paper_`.

## Remaining work

The complete future swirl energy must now be assembled through flatten,
power buffer, these angular bumps, steep transitions, waiting, collar and
infinite Gamma exterior. It feeds the positive ap root and the actual axial
end corrections. Then compose all five outer primitives and verify global
terminal closure and the whole-outer cone. C4 repair/cone bounds, the global
stress lift, flat remainder and true temporal recursion are not claimed by
this C1 repair layer.
