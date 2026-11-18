# ST073 signed finite-pressure difference in the fourth m15 contribution — 2026-11-18

The fourth contribution now has a signed difference formula between its actual regular pressure and cutoff sources and the nominal finite producer on the same \(q=2\) seam. The earlier [whole-axis calculation](CURRENT_ST073_REGULAR_COLLAR_WHOLE_Z_2026_11_18.md) supplies three independent diagnostics against the saved finite coefficients. This step also resolves one signed pressure-datum sector over the full real axial interval.

## Exact reduction of the nominal finite recipe

The archived finite recipe uses a degree-24 radius-four point inlet. Its point seed has zero \(\epsilon^2 F_0^2\) Taylor coefficients. The nominal cutoff-swirl product therefore vanishes algebraically for every \(Z\). Its point pressure density is constant in the radial quadrature variable. At \(q=2\) the four mesh pieces have total weight 2. Converting the coefficient from the \(s_c/4\) unit to the \(S=s_c/2\) unit gives exactly \(-16 N_{1p,f}(4,Z)/P_\star^2\) in the nominal real-arithmetic recipe. The true cutoff-swirl product is retained with its own support. The archived floating-point rounding is separate from this algebraic reduction.

Write \(L=1-\delta Z^2\), \(d=1-Z^2\), and \(p_f=\Psi_{f,\rho}(4,Z)/P_\star^2\), with \(p_{2f}\) the corresponding second radial derivative. The finite mean-derivative recipe gives the viscous part of the normalized nominal term as
\[
16\bigl(d\,\partial_Zp_f-(3+\delta)Zp_f-8Zp_{2f}\bigr)/L.
\]
The mean cancels as a shared algebraic expression before any enclosure. The inlet curvature and its ordinary axial derivative both enter.

## Signed transport on the actual collar

In this section write \(T,K\) in \(P_\star^2\) units. Let \(T_0=16p\), \(K_0=32p+64p_2\), using the actual inlet, and define \(T_f,K_f\) in the same way at the finite inlet. With \(B(T,K)=((1-\delta)ZT+d T_Z-2Z K)/L\), the viscous difference density is
\[
\chi B(T-T_0,K-K_0)+(\chi-1)B(T_0,K_0)+B(T_0-T_f,K_0-K_f).
\]
Its average over \(0\le s\le S\) contributes directly in \(\epsilon^2 h_b P_\star^2 S\) units. Inertia uses
\[
\big[(\rho^2-16)S_f+\rho^2(S_{\mathrm{flow}}-S_f)\big]/(2P_\star^2).
\]
The implementation expands the flow difference into first differences and products with common factors and accepts the factored \(\rho^2-16\) as an input. The true cutoff product contributes the third signed sector. Both orders use the ordinary axial derivative and the same endpoints. Ten exact algebra checks cover the mean cancellation, the three density differences at both orders, and the differentiated pressure-mass coefficient.

## Measured pressure mass sector

A saved finite pressure proxy at \(Z=0\) in the separate first-collar jets archive supplies the reference \(\overline M_2\). Keep it separate from the finite producer's exact internal scalar. The actual normalized fixed-\(\beta=2\) mass minus this reference lies in the directed interval \([2.5481317473355501,\,2.5481317527821548]\,10^{-20}\). The analytic coefficient of this difference follows from the exact \(Z(2d/(1+Z^2)^3+(1+\delta)/(1+Z^2)^2)/L\) inlet pressure slope. Its contribution has a strict value sign on 124 of 128 closed cells covering \([-1,1]\), and a strict derivative sign on 110 cells.

Subtracting this measured sector from the saved complete point-gap enclosures yields the following outward-rounded magnitude bounds in the common unit.

| \(Z\) | Value of remaining point gap | Ordinary derivative of remaining point gap |
| --- | --- | --- |
| \(-0.3\) | \(\lvert\cdot\rvert <1.79\times10^{-27}\) | \(\lvert\cdot\rvert <1.24\times10^{-26}\) |
| \(0\) | \(\lvert\cdot\rvert <2.62\times10^{-27}\) | unchanged axis enclosure, since this mass derivative is exactly zero |
| \(0.3\) | \(\lvert\cdot\rvert <1.79\times10^{-27}\) | \(\lvert\cdot\rvert <1.24\times10^{-26}\) |

This point comparison leaves the remaining finite inlet discrepancy and floating-point rounding uncertified over the whole axis. No finite approximation, uncertainty radius or interval intersection is equated with a physical source. The original six contributions, the 212-row ledger, all inherited errors, the reserved contribution and the independent full-pulse target are preserved. The terminal match remains open.

The frozen entry exited zero. Its [receipt](../experiments/root_st073/original_collar_finite_defect_receipt.json.gz) contains the ten source-recipe guards, ten exact algebra checks, all closed cells and the point differences. SHA-256 b425948f2301a86cd481b43cd9a771f89b625431a8862f47ef24eeb2a9a64e74. No ancestor recurrence or quadrature was replayed.

Next bind the remaining finite inlet terms or independently bound them using the zero-swirl seed and the normalized radial recurrence, with archival rounding tracked separately. The cut axial response in row three is another distinct signed source. It carries the actual all-word \(U_z^{(1)}\) and survives on the same cutoff support.
