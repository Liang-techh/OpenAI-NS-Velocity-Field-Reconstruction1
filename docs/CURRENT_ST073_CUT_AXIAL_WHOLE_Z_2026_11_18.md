# ST073 whole-axis third m15 cutoff axial contribution â€” 2026-11-18

The actual third forward contribution now has directed numerical enclosures for its value and ordinary axial derivative on 256 closed cells covering \(Z\in[-1,1]\). Both signs are resolved on 253 of the 256 cells for each order. In the archived unit \(\epsilon^2 h_b(s_c/4)\), the earlier saved point \(Z=0.371\) tightens to approximately \(107.88570255232\) in value and \(450.86816768\) in ordinary derivative. Both new enclosures lie inside the earlier independent saved enclosures.

## Actual source and signed response

Let \(a=s_c/4\), \(q=s/a\), \(\rho=4\exp(h_b s)\), \(x=\sqrt{\rho}\), \(\eta=1-\sigma(4q-2)\) and \(u_0=4Z+j\). The actual axial source from the raw six-row system is
\[
g_6/\epsilon
 = -2T_{-2}T_{-1-\delta}(U_{\rm lead})
   -4\epsilon Z\rho N_{1p}/L .
\]
In the axial kernel, \(T_k f=(kZ f+(1-Z^2)f_Z-2Z\rho f_\rho)/L\) with \(L=1-\delta Z^2\). The radius-independent signed source \(C=-2T_{-2}T_{-1-\delta}(u_0)\) has the exact raw zeroth-term response \(U^{(1)}_{z,\mathrm{axis}}/\epsilon=C\rho/4\). This follows by integrating the row-six Green kernel \(t\log(x/t)\) across the uninterrupted core/collar path, whose integral is \(x^2/4\). The calculation evaluates \(C\) and \(C_Z\) as rational functions of the original \(Z,\delta,j\).

The actual \(g_6/\epsilon-C\) source is controlled by the accepted nonaxis \(L^1_x\) bound, including pressure and the radial-source terms. The old natural-scale receipt supplies this source bound only. The verified raw matrix and preconditioned simplex majorant supply the first feedback word and the full tail beginning at word two separately. Whole-complex-tube bounds and a common axial Cauchy gap give ordinary derivative radii. The leading \(U_{\mathrm{lead}}-u_0\) contribution uses its own core norm and the true prescribed-shear collar equation. The complete enclosure adds each radius to the signed seed while preserving the actual all-word source definition.

The product is integrated against the positive exact support of \(\eta\). The signed seed product is \(u_0 C \int \rho^2\eta\,dq/2\). This weight is \(5\) at the zero-width limit and is bracketed with the actual nonzero width. The whole-cell outputs retain the signed model and separate radii for leading nonaxis velocity, zeroth-word nonaxis forcing, the first feedback word and the remaining tail.

## Same-unit comparison and correction

The conversion from the raw physical response to the legacy variable is
\[
u_{\rm cut}= \eta U_z^{(1)}/\epsilon,\qquad
R=\epsilon\rho,\qquad ds=(s_c/4)\,dq.
\]
The earlier cutoff binding receipt's description string had this conversion reversed. Its physical integrand \(2\epsilon h_b\rho\,U_{\rm lead}\eta U_z^{(1)}\), physical unit and selected numeric values are valid. The new analytic computation and the separate finite-mesh check both use the division by \(\epsilon\).

At the three legacy inputs, the actual column below displays rounded values in the common archived unit. Each of the six new intervals lies inside its corresponding directed archive interval. The two signed difference columns compare to the stored finite decimal and retain its provenance.

| \(Z\) | Actual value | Actual ordinary \(Z1\) | Actual value minus stored finite value | Actual \(Z1\) minus stored finite \(Z1\) |
| --- | ---: | ---: | ---: | ---: |
| \(-0.3\) | \(76.032\) | \(-437.76\) | \(-1.572\times10^{-78}\) | \(+6.415\times10^{-78}\) |
| \(0\) | \(3.05492734\times10^{-43}\) | \(4.89393812\times10^{-20}\) | \(+5.009\times10^{-123}\) | \(+4.607\times10^{-100}\) |
| \(0.3\) | \(76.032\) | \(437.76\) | \(+1.360\times10^{-78}\) | \(-5.012\times10^{-80}\) |

The finite midpoint producer freezes both axial factors at radius four. In exact real arithmetic the 16 midpoint nodes inside the transition step pair by reflection, so their cutoff mass is \(1/8\). With the retained inner half it totals \(5/8\). The normalized nominal third row consequently reduces to \(5 U_{\rm finite}(4,Z)u_{\rm finite}(4,Z)\). The two factors remain the archived finite recipes. The new formula isolates \(C-u_{\rm finite}\) and \(u_0-U_{\rm finite}\) as separate signed sources for further comparison. A whole-axis binary roundoff bound and the independent terminal match remain open. All 212 forward terms and inherited bounds are preserved.

## Reproduction evidence

The [whole-axis v2 receipt](../experiments/root_st073/original_collar_cut_axial_whole_z_v2_receipt.json.gz) exited zero, SHA-256 c65abf05dbe4210c79fda560e210f0b4da127f75c091dbde44499122442307cf. The [first entry](../experiments/root_st073/original_collar_cut_axial_whole_z_receipt.json.gz) rejected an equivalent symbolic representation of the Green primitive. The separate passing v2 checks the derivative by subtracting the two expressions and tests both endpoints. The [finite-mesh and unit receipt](../experiments/root_st073/original_collar_cut_axial_finite_mesh_receipt.json.gz) exited zero, SHA-256 b01d4da6a535f8ced432948deb03d97e68fe670e51249e972719cd4dc19e6d5d. The runs used archived inputs. No ancestor radial recurrence, Volterra iteration or quadrature was replayed.

Next connect the first two core moment contributions to the whole-axis third and fourth sources, retaining the actual core-pressure endpoint pairing and all active/compact error sectors. The independent negative full-pulse target and reserved own-swirl contribution retain their separate roles.
