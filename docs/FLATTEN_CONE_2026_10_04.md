# Whole original flatten two-vector cone — 2026-10-04

The whole original flatten region `t in [0,100], Z in [-1,1]` now passes the source-bound two-vector stress cone. It joins the accepted preceding-power, angular, entry, steep-power, exit, waiting and collar tail through the existing full-moment/pressure and completed physical right interface. Original fields, signed histories, sigmoid, lengths, radius, pressure datum and viscosity units are retained.

Nonzero cone coverage now starts at `Rtail*exp(-wait-Ts-102-Lrel)` and ends at `Rtail*exp(3)`; Gamma stress is exactly zero beyond. This is regional two-vector admissibility. The completed full-tensor cone, left pulse connection and global remainder/energy/recursion are still open.

## Native history comparison

Let `r=1-mu, k=1-a, b=(1-2a)/2, g=(mu-a)/r`,
`f=(1+Z^2)/2, rho=log(f), j=2Z^2/(1+Z^2)`,
`F=f^sigma(t/100), N=F*X` and
`K=Kright*exp((a-mu)*(t-100))*F/f`.
The actual inertial numerator divided by K is

`W=((k+b*j)*N-b*Z*N_Z)/F-1`.

Thus both K's axial dependence and the native history's axial derivative remain. The source equation `N_t=F-r*N` gives

`W_t+(r+rho*sigma_t)*W=mu-a+b*j*(1-sigma)-rho*sigma_t`.

For `V=W-g`, the shifted forcing is

`V_t+(r+rho*sigma_t)*V=b*j*(1-sigma)-rho*sigma_t*(1+g)>=0`.

The actual pulse source supplies the scalar
`Xp=buffer.power('0',1).Mtheta/buffer.power('0',1).Utheta` and
`Xv=1/r+(Xp-1/r)*exp(-13*r/mu)`.
The exact initial value is
`V(0)=b*j*Xv+k*(Xp-1/r)*exp(-13*r/mu)`.
Actual Xp and Xv positivity are checked before dropping the nonnegative first term. The signed second term is retained. Since F>=1/2,

`W>=g-2*k*abs(Xp-1/r)*exp(-13*r/mu)`.

The impossible-to-materialize exponent is controlled in logarithms. A certified memory cap g_lower*exp(-1000) is used only after bounding the actual memory below it. This cap is a bound, never a field value. The serialized near-one Xv interval is not chosen as the defining function.

## Whole-region margins and original units

Eight continuous source intervals cover the original sigmoid x in [0,1]. They bound actual sigma_t=sigma'_x/100 by approximately .266813290066476. Original logistic/reflection derivatives prove nonnegative slope independently. The flat-tail majorant constant is not used as a global derivative bound.

The actual variable shear satisfies
`0<kappa-2=2mu-2rho*sigma_t<=.369881759491001<2`.
The source shear is strictly negative from exact S>0.
K>=Kright, the W lower bound and the whole source axial stress enclosure prove theta positivity and the original directional inequality.

AST replay preserves the actual inlet B log at
`q=-wait-Ts-102-Lrel`. Its correlated pieces are
`-k*Ts+(a-mu)*(100+Lrel)+2*bh`, with the remaining source log terms unchanged. The inequality is tested in original two-vector normalization and physical viscosity units.

The generally nonzero physical remainder, including axial viscosity, remains. The completed diagonal `Ttheta_theta=r*partial_z(Tz)` is retained by the physical adapter; proving its full-tensor cone is a separate task.

## Reproduction and evidence

Run `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage flattencone`.

Focused producer/checker PASS: 398 current hashes, 10 symbolic identities, 14 source bridge facts, 26 strictly positive necessary source margins, and the continuous eight-interval slope cover. A moderate independent original-sigmoid fixture has actual Xp and Xv below equilibrium and nonzero negative incoming memory. Four direct-history/integrating-factor comparisons pass at tolerance 1e-60; largest error is below 4.1e-86. These samples check the formula; the actual cone proof uses the continuous source comparison and interval bounds. Existing unchanged Cartesian physical/operator evidence is consumed without rerun.

`CertifiedFlattenPhysical` wraps the unchanged current physical evaluator with receipt/family/hash-checked regional cone metadata.

## Next work

Recover original selected pulse-end stress on s in [-4,0], including both beta supports centered at -3 and -1, full five moments, nonzero meridional terms, and compatible pressure. Close its s=0/flatten t=0 source and physical interfaces. Existing pulse C4 velocity continuity alone is insufficient.

Then extend through original pulse/gap/main/entrance regions, resolve upstream finite-width bridge feedback and remaining implicit data, and establish completed full-tensor/global admissibility, independently bounded global flat remainder, required-domain energy and physical-volume norms. Actual coupled n=1, n-dependent higher orders, smooth summation, mean/oscillatory stress cancellation, full corrected Cartesian residual and measured multitime dynamics remain uncompleted.
