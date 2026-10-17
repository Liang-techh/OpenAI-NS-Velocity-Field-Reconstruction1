# Original same-generator swirl and all-sector source C1

Tasks 14/15 remain **IN PROGRESS**; the long-term goal remains **ACTIVE / INCOMPLETE**. This stage recovers signed factored swirl source jets through ordinary Z2, encloses all pressure/hydro/swirl components at three Z queries, and computes nineteen all-sector raw46 main-term C0/C1 queries. The complete weighted primitive, raw47 signed remainder and global m15/five-moment closure remain open.

## Same-source amplitude and derivatives

The known comparison is still seeded from the analytic core at rho=4, not from actual first-collar exit increments. Saved phase-zero Phi/B/C inlet Taylor rows and the dressed C rows through order three are converted by k!. B is its own square moment and C is its own pressure moment; neither is replaced by a target moment.

The new consumer also hydrates saved anchored G(Z), logF0(Z), selected-logCstar enclosures and gradient Taylor rows at Z=-.3,0,.3. These are the existing certified pole/anchored-primitive data; no root solver, primitive integration or producer is rerun. The exact source remains

```
A(Zbase)=F0(Zbase)^2=Cstar^-2*exp(-2Lambda*G(Zbase)) > 0.
```

Pointwise logA=2logF0 is contained in the saved global amplitude log enclosure. Gbar is retained only in that global bound. The atom has no selected scalar value and retains the actual anchored G source packet and digest. Eleven specific AST guards bind the swirl direction/mode split, relative normalization, saved C dressing, exact comparison histories, anchored amplitude/gradient/parameters and the original core gradient/sigma. The two H formulas agree algebraically:

```
Z*(1-delta)/2+(1-Z^2)*(4Z+j)
 = -4Z^3-jZ^2+((9-delta)/2)Z+j,    sigma=j/500.
```

The relative F0 squared jet is normalized at its base point. Its approximation comes from the same-source dressed-C/raw-C product identity, with r0=1 imposed before interval division. Its truth enclosure is independently built from saved ordinary Gprime/Gsecond/Gthird rows and Lambda=1/epsilon:

```
l1=-2Lambda*Gprime, l2=-2Lambda*Gsecond, l3=-2Lambda*Gthird,
r=[1,l1,l2+l1^2,l3+3l1*l2+l1^3].
```

The independent moment-quotient enclosure remains serialized. Ordinary R_n=(D_n-sum(j=0..n-1) binom(n,j)R_j C_(n-j))/C_0 is used, with a positive C_0 enclosure. No interval-overlap argument proves source equality. The external A(Zbase) is frozen during these local derivatives; its derivatives are already carried by r and are not applied twice.

## Finite-width comparison and radial modes

A new order-three weighted-jet proof controls the logarithmically smoothed comparison Phi. With M=4*XhPhi, Htail=M*sum((weight/Zradius)^k,k=1..3), the inverse-jet margin Phi_floor-Htail is positive, approximately 0.0722090076827. Cauchy bounds and the positive controls give

```
K=4*exp(2hcap)*M*sum((weight/Zradius)^k,k=0..3)/(rhoRadius*(Phi_floor-Htail)),
PhiCover=norm(Phi0)*exp(2hcap*K),
norm(DeltaPhi)/hb <= 2K*PhiCover,
norm(DeltaB)/hb <= 2*(2norm(B0)+PhiCover^2),
norm(DeltaC)/hb <= 2*(norm(C0)+PhiCover^2).
```

The last two formulas use the exact own ODEs B_y+2B=Phi_bar^2 and C_y+C=Phi_bar^2. All 36 new log-width proofs pass after conversion to ordinary k!/weight^k errors and combination with the original log(hb). Width caps only widen truth. The inherited native2 radius and same-Phi quotient proofs are hydrated unchanged.

Dress the common target S=Phi^2, dB=B-Phi^2/2 and dC=C-Phi^2 with r. The three swirl numerator modes use target (S/2,S), j1=(0,dC), and j2=(dB,0). They share denominator 2L and have delta-degree zero/one pieces

```
N0=(1-Z^2)*(b_Z-c_Z)+2Z*c,
Ndelta=2Z*(c-b),    L=1-delta*Z^2.
```

All six mode/delta pieces remain separate. For mode j=0,1,2 and delta degree k=0,1, their source scale before A is

```
100^(2-j)*rho2^j*epsilon^(j+k-1)*exp(970*k),
R0=epsilon*rho2, rho2=4exp(2hb), delta=epsilon*exp(970).
```

This preserves the original R0 origin and all three radial modes. The common positive 2L denominator is evaluated without a formal quotient truncation. The inherited macro exp(ell)-1 norm remainder is included only in truth, with ordinary derivative coupling.

## Signed coefficients and all-sector queries

| Z | Swirl relative B Z0 approximation | Relative B Z1 approximation | Relative B Z2 approximation |
| --- | --- | --- | --- |
| -.3 | 1.78068425e+817812180069139138 | 2.23402565e+1226718270103708706 | 2.80278246e+1635624360138278274 |
| 0 | -1.39049948e+817812180069139160 | 1.2396283e+1226718270103708750 | -1.10512686e+1635624360138278340 |
| .3 | -1.78068425e+817812180069139138 | 2.23402565e+1226718270103708706 | -2.80278246e+1635624360138278274 |

These signed coefficients multiply the exact positive A atom. They are not physical velocities or exact function values. The code combines pointwise logA with each coefficient upper norm before admitting an absolute 1e-1000 error box. All nine amplitude log-product proofs pass. The numerical approximation zero for this extremely small component is controlled by those proofs; the actual amplitude stays positive and the signed source definition remains factored; no exp(logA) underflow or cap is used to define the field.

All three source components now have ordinary B/B_Z/B_ZZ enclosures at the three queries. Nineteen C(B)*M(s) rows retain both signed factored swirl C0/C1 and the controlled absolute swirl approximation, then add them to the saved pressure-plus-hydro main terms. Saved masses and R100*epsilon*hbS^2 are unchanged. This covers the raw46 main component only; raw45 actual-radius variation and raw47 weighted-primitive remainder remain active.

Production exits 0 in 2.271 s. The report is 945,099 bytes, SHA256 `b870f9f8957295e9c22a96d077b144a564ffe97278ccec4fe1f13c26dbab14aa`, with 1961 input pins and 231 implementation pins. Owner-free serialized inspection passes.

All 30 new pure-kernel controls and 139 new transport/amplitude controls pass. Independent oracles differentiate nonconstant functions, explicit exponential Phi/B/C histories, relative-amplitude Bell jets and a finite nonzero absolute amplitude. They check normalized C1 without double amplitude differentiation and reject missing source proofs, Z3 rows, positive denominators or inverse-jet margins. These checks do not validate the original NS residual.

Pressure/hydro packets, saved mass, raw213/reduced183, all old errors, six signed inner/lower terms, the independent negative full-pulse m15 target with reference*A*Rv^2 and reserved n1 own-swirl once retain exact parent witnesses. Accepted ancestor calculations are not replayed. Existing unrelated modified/untracked files remain intact.

## Next executable task

Recover the signed weighted source primitive on first-switch phase0..s using these same native2 modal sources. The inherited `original_first_switch_m15.source_norms` bounds radial factors with rate=power-j+radius_weight and stores only a positive `weighted_primitive_remainder_per_hb`; that bound is not a signed function. Implement new signed weighted pulse moments/radial kernels and bind the actual angular multiplier, preserving ordinary Z1 and distinct hbB/hbS, before altering raw46/raw47. Keep raw45 actual R(s)-100 variation, raw23/24/43..47 and all earlier errors until a same-source identity justifies a change. No accepted mass quadrature or old producer should run.

Total signed m15 and functional all-five/common-N closure, exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2 recursion/summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS 1e-3 gate is not applied to this leading-source stage.
