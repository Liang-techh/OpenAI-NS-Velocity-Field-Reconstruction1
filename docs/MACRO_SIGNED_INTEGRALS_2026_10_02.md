# Leading signed macro integrals — 2026-10-02

The comparison histories now feed the leading formal-width coefficient of the actual signed macro contribution from the exact two-microswitch endpoint (R_2=R_a e^{2h_b}) to (R=100). For the (h_b^1) coefficient, the frozen comparison fields are their (h_b=0) inlet values and the six comparison moments retain their exact rates and targets. Their directions from the original (9.13) are finite exponential modes in (Delta=y-2h_b); the axial drive uses the original hydrostatic, analytic-pressure, and swirl split.

The generator analytically integrates those modes using (Y=log(100/R_a)) and simplifies (R_a e^Y=100) before interval arithmetic. This matters here: (R_a) is about (5.13\times10^{-408906090034569568}), so (Y\approx9.41541067348080989\times10^{17}). Directly materializing (e^Y) or (F_0^2) creates unusable intermediate scales. The pressure and swirl amplitudes therefore remain exact logarithmic factors, paired with directed normalized axial-jet intervals; their sum is an explicit factored expression, not an evaluated floating-point representative.

Reproduce with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals.py
python experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals_check.py
```

Both runs pass. The independent checker rebuilds the modal directions from the accepted inlet and terminal comparison packets, verifies all six (R=100) moment identities and each signed drive mode at (Z=0), (Z=0.5), and the exact shared root, and independently quadratures the (J_F) integral after regularizing the long logarithmic radius interval. The exact source/family/pressure datum hashes are carried through.

This is only the macro contribution (R_2\to100) at leading formal order:

\[
\log(F_{100}/F_2)=h_bJ_F+O(h_b^2),\qquad
V_{100}-V_2=h_b\sum_j e^{\ell_j}J_{V,j}+O(h_b^2),
\]

where the (\ell_j) are the retained pressure/swirl scale logs. The two microscopic (R_a\to R_2) contributions, the original (100\to110) first switch, higher formal-width coefficients and their remainder propagation, actual moment feedback beyond leading order, and the implicit five-bump recovery are still unresolved. Accordingly `actual_signed_bridge_completed` remains false. Do not present this partial macro term as the full (I_{bridge}) or a solved physical axial velocity.

The broader required-domain energy, admissible stress, independent flat remainder, global dynamics, and true scale-indexed coefficient recursion/corrections also remain open. This milestone advances an upstream dependency of scale recursion; it does not implement recursion.
