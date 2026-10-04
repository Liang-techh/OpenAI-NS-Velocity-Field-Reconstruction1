> Update2026-10-04: R100..110 now has accepted signed source-integral enclosures and corrected width coefficients; see SWITCH_SIGNED_INTEGRALS_2026_10_04.md. Actual higher bridge orders and point recovery remain open.

# Leading signed core-to-R100 bridge coefficient — 2026-10-02

The comparison histories now feed the leading formal-width coefficient from $R_a$ through the original microswitch charts to $R=100$. This includes the first actual core chart and the signed macro integral from the two-chart endpoint $R_2=R_a e^{2h_b}$ to $R=100$. At order $h_b$, the actual cutoff weight on the first chart is $1-\sigma(s)$, whose exact integral is $1/2$ by the pulse reflection symmetry. The second micro chart has $\chi=h_b$, so its leading $h_b^1$ coefficient is zero. On the macro, the frozen comparison fields are their $h_b=0$ inlet values and all six comparison moments retain their exact rates and targets. Their directions from the original (9.13) are finite exponential modes in $\Delta=y-2h_b$; the axial drive uses the original hydrostatic, analytic-pressure, and swirl split.

The generator integrates those modes using $Y=\log(100/R_a)$ and simplifies $R_a e^Y=100$ before interval arithmetic. Here $R_a\approx5.13\times10^{-408906090034569568}$ and $Y\approx9.41541067348080989\times10^{17}$. Directly materializing $e^Y$ or $F_0^2$ creates unusable intermediate scales. Pressure and swirl amplitudes remain exact logarithmic factors paired with directed normalized axial-jet intervals; their sum is an explicit factored expression, not a floating-point representative.

Reproduce with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals.py
python experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals_check.py
```

Both runs pass. The independent checker rebuilds modal directions from the accepted inlet and terminal comparison packets, verifies all six $R=100$ moment identities and each signed drive mode at $Z=0$, $Z=0.5$, and the exact shared root, quadratures the first-chart pulse weight, and independently quadratures the full leading $J_F$ coefficient after regularizing the long logarithmic radius interval. The exact source, family, pressure datum, and pulse source hashes are carried through.

The accepted packet gives the leading bridge contributions

$$
\log(F_{100}/F_a)=h_bJ_F+O(h_b^2),\qquad
V_{100}-V_a=h_b\sum_j e^{\ell_j}J_{V,j}+O(h_b^2),
$$

where $\ell_j$ are retained pressure/swirl scale logs. This resolves only the leading $h_b^1$ coefficient of the prescribed-shear bridge through $R=100$. The second microscopic chart's order-$h_b^2$ term, the original $100\to110$ first switch, higher formal-width coefficients and their remainder propagation, actual moment feedback beyond leading order, and implicit five-bump recovery remain unresolved. Accordingly, `actual_signed_bridge_completed` remains false. Do not present the leading coefficient as a full finite-width $I_{bridge}$ or a solved physical axial velocity.

Required-domain energy, admissible stress, independent flat remainder, global dynamics, and true scale-indexed coefficient recursion/corrections also remain open. This advances an upstream dependency of scale recursion; it does not implement recursion.
