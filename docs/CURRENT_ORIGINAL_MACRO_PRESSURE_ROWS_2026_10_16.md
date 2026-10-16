# Original macro pressure source rows through R=100 - 2026-10-16

`OriginalMacroPressureRows.evaluate(Z,q)` now provides the missing pressure-integration rows on the actual macro chart, q in [0,1]. It exports background U through ordinary Z2, mean M through Z3, Q through Z2, and U_y/Q_y through Z1, each with separate original selected-N logarithmic error bounds. The background mean uses genuine upstream jets before the older velocity packet truncates them. No missing row is padded with zero.

The full goal is **ACTIVE / INCOMPLETE**; tasks 14/15 remain IN PROGRESS. This implements the macro pressure source and finite endpoints. It has **not yet integrated cumulative macro m15**. The checked seventeen-term second-micro moment remains the next integral's inlet. Genuine temporal n>=2 recursion and global corrected NS remain unfinished.

## Actual macro geometry and incoming mean

The source is the original `OriginalMacroNumeric` family and exact selected repair integer. Geometry is

```
Y = logLambda + log(25)
L = Y - 2hb
R(q) = 4epsilon*exp(2hb + L*q)
logR = log(100) - Y*(1-q) + 2hb*(1-q)
dy = L*dq; R*dR = R^2*L*dq
```

Physical epsilon and the very large normalized rho are combined in logR before exponentiation. At q=1 the defining physical radius is exactly 100. No microscopic rho approximately4 or extrapolated core radial derivative is used in the macro source. The original hb is kept as a source log atom; its cap only bounds errors.

The actual leading mean uses the nonzero second-micro exit Mpsi and exp(-L*q). Its core-source part is weighted by 1-exp(-L*q), and the positive Volterra kernel bounds the full changing-Psi source by the uniform original modal increment times that mass. The signed numerical proxy uses the corresponding zero-width inlet and exp(-Y*q). It changes on the actual large-Y relaxation scale, including at q=1e-18. The full first and second field histories are retained once.

## Genuine macro radial source

The original source equation is

```
U_y = -hb*(Phi/barPhi)*(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl).
```

`source_Uy` uses the actual frozen-comparison macro modes R0^j*R^(power-j), the original quotient, and all three source parts. The accepted six-coefficient weighted ell jet has norm at most6*cap, giving the explicit nonlinear coefficient error

```
||source|| * 6*cap*exp(6*cap).
```

The signed numerical coefficient and full directed error are exported before the one hb factor, as `leading_U_y_factored`. The raw physical derivative uses a zero numerical proxy and its complete width error; that proxy does not declare the factored source or actual field zero. No derivative of an integral majorant substitutes for the original drive. Neither this proxy nor a comparison source is a resolved selected-N field center.

## Pressure and finite endpoints

Recovery uses the exact identities

```
Mpsi_y = Psi - Mpsi
Q = [2Z*U - (1-delta)Z*M - (1-Z^2)M_Z] / (1-delta*Z^2)
Q_y = [2Z*U_y - (1-delta)Z*(U-M) - (1-Z^2)(U_Z-M_Z)] / (1-delta*Z^2).
```

Q_ZZ retains the genuine third mean derivative and its finite-N logarithmic error. The available conservative M_ZZZ error is not numerically small, so a small-error corrected M_ZZZ/Q_ZZ value is not claimed.

Normalized Psi, Mpsi, Qpsi and their radial rows are exported separately from the affine physical background. This preserves epsilon-weighted sectors that could otherwise disappear when added to an O(1) field at finite precision. In particular,

```
background R*Q_y = exp(logR-logLambda) * Qpsi_y
finite-N endpoint error = exp(logR) * error(deltaQ_y).
```

These endpoint factors differ by epsilon. Ordinary Z1 is retained for each. The macro fast-phase errors use the original bridge_macro source and physical B/Pstar units:

```
deltaU_y  = Pstar*(B_y/N + B_phi)
deltaU_yZ = Pstar*(B_yZ/N + B_Z_phi).
```

The fast terms have neither another N denominator nor an extra hb. The mean history integrates all original first/second/macro source cells, with their original measures; Y only upper-bounds the macro length in error caps. Distinct axial frames are preserved.

## Focused evidence and reproduction

The focused run exited0 in77.408 seconds; saved-report inspection exited0. It produced19 actual macro packets,266 ordinary physical rows,570 serialized-error conditions,38 fast-phase conditions,38 finite-endpoint unit conditions,39 genuine meanZ3/QZZ conditions,57 modal-source conditions, six mean-inlet/relaxation conditions, six exact-R100 geometry conditions and six rejected invalid requests. The original macro_cell source assignments were independently executed symbolically for value/Z1; Q+Q_y and fast chain identities passed. Twelve independent Volterra quadratures cover mean ordinary Z0..3 with evolving sources and nonzero incoming means.

The1,401,681-byte report SHA256 is `28fac21299ae929447b45894331e33538d51687918323137cbecc4ccf90aad13`; it pins1751 inputs and37 implementations. Inspect using `python scripts/run_original_macro_pressure_rows.py --inspect-report`. Do not rebuild checked ancestors. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. Its review found no material source, geometry or endpoint-scale underbound.

## Next actual construction

1. Load the checked complete second_micro s=2 seventeen-term moment unchanged, with its exact source and implementation pins. Add only the macro pressure increment, not its incoming moments again.
2. Implement the original reduced identity m15_increment=S1/2-[R*Q_y]. Retain both macro endpoints, all six weighted volume products and ordinary Z1. Use actual R*dR and L; there is no extra hb on this macroscopic measure.
3. Integrate the signed exponential mean proxy with stable kernels. Keep the affine U0, epsilon-linear and epsilon-squared sectors distinct where needed; retain nonzero incoming-mean exponential modes. Include directed leading-source, geometry and numerical errors. A repeated uniform box is not a quadrature-convergence test.
4. Add full selected-N product errors, including quadratic cross products, Q_ZZ/M_ZZZ logs, full history and both nonflat endpoint errors. Preserve the existing correction signs as unknown until evaluated from correlated source functions.
5. Continue the original switch/remaining active/compact pressure and source-correlated independent Rp value/Z1 comparison. An exact terminal theorem or interval overlap does not measure that cancellation.

All-five global joins, common global N, analytic preheat/exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2/smooth sum, oscillatory correction and independent corrected Cartesian NS remain open.
