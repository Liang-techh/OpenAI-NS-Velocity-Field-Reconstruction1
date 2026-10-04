# Original preceding-power two-vector cone - 2026-10-04

## Accepted increment

The original `outer_power` chart, phase in [0,1], Z in [-1,1], and length Lrel-4, now has a uniform regional two-vector cone using the same native forward history, complete moments, absolute pressure, velocities and physical units as the accepted stress/physical companions.

The focused producer and current checker pass. The receipt binds 386 current input hashes, 12 exact symbolic identities, 18 source-bridge facts and 24 positive source margins. The independent correlated-flatten fixture has three direct native-integral identity checks and three positive-lemma comparisons at 1e-60 tolerance. Its moderate sample parameters validate formulas; the actual cone proof is a continuous source inequality over the whole original domain.

## Native cumulative history and the actual flatten source

Set r=1-mu, k=1-a, b=(1-2a)/2, u=(Lrel-4)*phase. The native power source gives

```text
X = 1/r + (Xf-1/r)*exp(-r*u)
W = k*X - b*Z*X_Z - 1
  = (mu-a)/r + Hf(Z)*exp(-r*u)
Hf = k*(Xf-1/r) - b*Z*Xf_Z.
```

Independent interval boxes for Xf and Xf_Z cannot establish Hf positivity. The proof preserves their correlation from the original 100-unit flatten construction.

Let f=(1+Z^2)/2, rho=log(f), sigma=sigma(v/100), omega=1-sigma, and j=2Z^2/(1+Z^2). The actual source uses F(v)=exp(rho*sigma), F(100)=f, and an unnormalized Xint. Source AST replay binds q, rho, F, Fc, the weighted Xint cell sum and Xf. It verifies explicitly

```text
Fc/F(100) = f^(-omega)
Xint = f*I_norm
Xf = exp(-100*r)*Xv/f + I_norm
I_norm = integral_0^100 exp(-r*(100-v))*f^(-omega) dv.
```

The original source Xint is therefore not silently identified with I_norm. The completed bridge replays the original F/Fc expressions and endpoint normalization before applying the continuous integral bound.

Differentiating this same integral yields

```text
Hf = exp(-100*r)*[(k+b*j)*Xv/f-k/r]
   + integral_0^100 exp(-r*(100-v))
       *[k*(f^(-omega)-1)+b*j*omega*f^(-omega)] dv.
```

Xv is positive; it need not exceed 1/r. Exponential/log tangent inequalities and the nonnegative spatial remainder

```text
k*(1-Z^2)/2 + b*2Z^2/(1+Z^2) - b
 = a*(1-Z^2)/2 + b*Z^2*(1-Z^2)/(1+Z^2)
```

give an integrand lower bound b*omega. The original sigma odds imply sigma(v/100)<=1/2 on v in [0,50]. Thus

```text
Hf >= exp(-100*r)*[b*expm1(50*r)-2*k]/(2*r) > 0.
W >= (mu-a)/(1-mu) > 0
```

on the entire original phase/Z domain. The actual strict parameter condition is checked. No incoming signed term is discarded before this bound is established.

## Theta, source shear and original inlet amplitude

With y=(Lrel-4)*(phase-1)<=0, K=Kright*exp((a-mu)*y) and L=1-2a*Z^2,

```text
Ctheta = K*[W/L - 2*S*exp(-q)*(1+mu)]
kappa-2 = 2*mu
qleft = -wait-Ts-2-Lrel.
```

The actual positive inverse-radius S establishes strict negative shear; its enclosure upper endpoint only bounds magnitude. The continuous W bound dominates the whole-domain subtraction. K_Z=0 on this chart does not remove the full A/E/P axial derivatives, and the unchanged full-moment Cz bound is retained.

Actual source B log parts at qleft are replayed and correlated before enclosure:

```text
steep = -k*Ts
flatten_power = -100*(1/2+mu)+(a-mu)*Lrel
waiting_and_current = 2*(1/2+a).
```

All other source parts and the full source sum remain identical. The resulting uniform sufficient directional inequality is

```text
log(2*mu) + 2*log(Bmax) + 2*log(Cz_abs_upper/Ctheta_lower) < log(2).
```

No phase grid, raw-theta interval positivity, shortened Lrel, or cap-selected field is used.

## Composition, interface and scope

Current same-function power/angular K4, complete moments, stress mixed3, pressure mixed4 and completed physical right joins compose the accepted nonzero two-vector tail:

```text
Rtail*exp(-wait-Ts-2-Lrel) <= R < Rtail*exp(3).
```

Gamma stress is exactly zero separately. `CertifiedOuterPowerPhysical` checks the current cone receipt/family/hashes and adds regional metadata around the unchanged physical evaluation. The completed diagonal and generally nonzero axial-viscosity remainder remain.

The controller exposes `--stage outerpowercone`. Compilation and stage listing pass.

The following remain false: original flatten stress/physical left join; completed full-tensor/global cone; independently bounded global flat remainder; required-domain physical energy; actual n-dependent recursion; oscillatory correction; final full corrected Cartesian residual. The use of flatten history above proves a boundary-history inequality, not the cone or physical decomposition throughout flatten itself.

A read-only GPT-5.6 Luna / max review identified the Xint normalization gap; the source bridge was corrected and the worker confirmed the connection. The producer and checker reported here were rerun after that correction.

## Next implementation

1. Restore the same full A/E/P and analytic absolute-pressure datum through the original 100-unit flatten, retaining its original sigmoid and variable K_Z/K_ZZ.
2. Close both complete-moment/stress/pressure function joins, then the physical velocity/completed stress/remainder mixed joins with the same radius and units.
3. Establish the whole original flatten theta, strict shear and directional bounds from actual variable rates and continuous source histories.
4. Finish finite-width bridge feedback and the global admissibility/flatness/volume/required-domain energy prerequisites.
5. Solve the actual n-dependent recovery and independent moment repair, then finite-order remainder/smooth summation, oscillatory quadratic stress cancellation and final corrected residual.
