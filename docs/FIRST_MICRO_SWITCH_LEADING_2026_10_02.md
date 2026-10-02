# First micro-switch leading correction — 2026-10-02

The first switch chart in the accepted construction is the short interval $R=100e^{h_bs}$, $0\le s\le1$. The implementation now evaluates its leading $h_b^2$ signed field corrections at the exact $R=100$ comparison-history endpoint for $Z=0$, $Z=0.5$, and the shared root. It uses the original switch equations: $\partial_s\log F=-\frac12h_b^2(1-\sigma(s))\bar D+O(h_b^3)$ and $\partial_sV=-h_b^2(1-\sigma(s))(\phi_{actual}/\bar\phi)G+O(h_b^3)$.

At this order, $R=100$, the actual/comparison quotient is $1$, and the six moment histories equal their zero-width comparison endpoint. Pulse reflection symmetry gives $\int_0^1(1-\sigma(s))ds=1/2$; an independent high-precision quadrature verifies this. Pressure and swirl prefactors remain logarithmic and factored. The checker reconstructs the directions and verifies the signed axial-jet coefficients and source logs at all three points.

Reproduce with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_first_switch_leading.py
python experiments/root_st073/lei_ren_part1_paper_compliant_first_switch_leading_check.py
```

This is only the leading $h_b^2$ term on the first micro chart. It does not cover the second chart $100e^{h_b}\to100e^{2h_b}$, the subsequent power segment to $R=110$, or controlled higher-order remainders. Full switch completion and `actual_signed_bridge_completed` remain false. The leading full bridge coefficient through $R=100$ remains documented separately in `MACRO_SIGNED_INTEGRALS_2026_10_02.md`.
