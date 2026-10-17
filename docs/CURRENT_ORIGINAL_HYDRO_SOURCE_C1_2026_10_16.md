# Original same-generator hydro source C1

Tasks 14/15 remain **IN PROGRESS**; the long-term goal remains **ACTIVE / INCOMPLETE**. This stage evaluates the original macro hydro component through ordinary Z2, and nineteen hydro and pressure-plus-hydro raw46 main-component C0/C1 rows. Swirl and the complete weighted primitive remain open. No m15 term or inherited error is replaced or removed.

## Source and finite-width comparison

The known macro direction consumes the comparison at native phase 2. The pinned `prepare -> comparison.inlet -> comparison.micro` source edge rebuilds the analytic core at rho=4; it does not consume the actual prescribed-shear first-collar exit. The saved phase-0 first-collar inlet rows for Psi, Mpsi and Apsi are therefore the correct zero-width source leaves. Their Taylor rows through Z3 are multiplied by k!. Actual first-collar increments are not added to this separate comparison. Eleven specific AST guards bind source formulas, moment targets, radial weights, exact comparison kernels, normalization, the core rebuild edge and Cauchy bounds.

The exact comparison V is a positive integration-by-parts average of core V along rho=4exp(y). Its alpha and negative alpha derivative weights sum to one and are Z-independent. On the same core continuation tube, the saved Cauchy majorant gives Cpsi=4*XhPsi*sum((weight/Zradius)^k,k=0..3). The domain upper majorant dominates the separately saved same-core Psi norm. No lower norm endpoint is used as an upper bound.

With u0=4Z+j, the exact normalized histories are

```
V=u0+epsilon*psi_bar,
M=u0+epsilon*m_bar,
A=u0^2+2*epsilon*u0*m_bar+epsilon^2*a_bar,
m_bar_y=psi_bar-m_bar, a_bar_y=psi_bar^2-a_bar.
```

The independent Apsi leaf remains separate from Mpsi squared. Positive exponential kernels imply the before-width order-three jet bounds

```
Delta psi / hb: 8*exp(2*hcap)*Cpsi/rhoRadius,
Delta m_bar / hb: 2*(Cpsi+norm(Mpsi0)),
Delta a_bar / hb: 2*(Cpsi^2+norm(Apsi0)),
(rho2-4)/hb: 8*exp(2*hcap).
```

The original log(hb) is combined with each ordinary k! times coefficient/weight^k before admitting a representable 1e-1000 cap. All 39 log-width proofs pass. These positive bounds widen source truth; they never supply signed values, select hb or choose rho2. The original smoothing interval stays within rho<4.1.

## Untruncated hydro algebra

The source direction is

```
d=1-Z^2, L=1-delta*Z^2, W=1-(1-delta)*Z*M-d*M_Z,
Dhydro=[-W*V+(1-delta)*(M-Z*M_Z)/2+2*delta*Z*A-d*A_Z]/(2L).
```

The helper retains a formal epsilon polynomial, including delta=epsilon*exp(970), and keeps the common denominator 2L separate without a formal quotient truncation. The exact j=1/5721363714816000000000 is hydrated from both saved core levels and checked against the pinned fixed-integer definition.

At the second-micro exit, target moments are M=V and A=V^2. The j=1 mode is formed directly as a linear moment difference; its normalized inputs are Mpsi-Psi, 2*u0*(Mpsi-Psi), and Apsi-Psi^2. Its epsilon-zero sector vanishes by this shared-source identity before interval subtraction. Hydro j=2 is exactly zero by the pinned source formula. At R100 the weights are [100,R0,R0^2/100], with R0=epsilon*rho2. The implementation evaluates

```
Bhydro=-(Phi0/barPhi2)/epsilon*(100*Dtarget+R0*Ddifference)*exp(ell).
```

All numerator sectors 0..4 and moment-difference sectors 1..3 remain serialized; each numerator sector n carries epsilon^(n-1) after normalization. Ordinary derivative factorial coupling is preserved. The saved same-Phi quotient proof is hydrated unchanged. The point uses q=1, ell=0, rho2=4 and saved numerical core approximations; macro exp(ell)-1 derivative coupling is added only to the truth enclosure with the inherited weighted norm bound. Final scalar sums remain numerical approximations; the separate factored sectors retain contributions that can be lost in finite-precision sums.

## Numerical results and scope

| Z | Hydro B Z0 approximation | Hydro B Z1 approximation | Hydro B Z2 approximation |
| --- | --- | --- | --- |
| -.3 | -2.00072711e+408906090034569570 | 4.42528425e+408906090034569570 | 2.2438061e+408906090034569571 |
| 0 | 3.06390854e+408906090034569548 | 7.7909934e+408906090034569570 | -1.63408456e+408906090034569549 |
| .3 | 2.00072711e+408906090034569570 | 4.42528425e+408906090034569570 | -2.2438061e+408906090034569571 |

These are before-width source coefficients, not physical velocities, exact source values or stress-cone signs. Directed enclosures and approximation errors remain serialized. Nineteen C(B)*M(s) hydro and pressure-plus-hydro main-term queries use unchanged saved masses and physical R100*epsilon*hbS^2. The accepted pressure source and C1 packets retain exact digests. Raw47 is retained; the weighted primitive is incomplete.

Production exits 0 in 2.48 s. The report has 818,067 bytes, SHA256 `660f377e75998bc631298d431fb2388fc89101f7b9cf58e5dd566bf6b2262033`, 1960 input pins and 226 implementation pins. Owner-free serialized inspection passes.

All 33 new hydro algebra controls and 78 new transport controls pass. Independent oracles use mp.diff of nonconstant exp/sin/polynomial fields, finite epsilon evaluations, the correlated delta identity, and explicit exponential integrals of a y-dependent comparison field and its square. The transport controls test the original hb-times-coefficient bounds as well as the coarser representable boxes. These controls do not validate the original NS residual.

Raw213/reduced183, all prior error budgets, six signed inner/lower terms, the independent negative full-pulse m15 target with reference*A*Rv^2, and reserved n1 own-swirl once remain inherited by exact parent witnesses. Saved mass quadrature, pressure production and all accepted ancestor constructors/producers/checkers/fixtures/calculations are not replayed. Existing unrelated modified/untracked files are preserved.

## Next executable task

Recover the same-generator swirl B/B_Z/B_ZZ using the original comparison Phi/B/C histories and dressed relative F0^2 derivatives, retaining the formal exact positive absolute F0^2 scale and its source definition. Avoid differentiating the external base amplitude twice or selecting an underflow zero as the field. Combine the signed sectors, recover the full weighted primitive/remainder, then compare coupled contributions with the independent negative complete-pulse m15 target.

Total signed m15 and functional all-five/common-N closure, exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2 recursion/summation, mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS 1e-3 threshold is not applied to this leading source interface.
