# ST073 whole-axis finite inlet remainder for the fourth m15 contribution — 2026-11-18

The signed fourth-row difference now separates the entire \(\beta=2\) pressure-mass coefficient from a uniformly enclosed remainder. This uses the finite zero-swirl source isolated in the [preceding calculation](CURRENT_ST073_FINITE_COLLAR_DEFECT_2026_11_18.md), interpreted in exact real arithmetic with an explicit nominal mass \(M_f\). The supported mass range is \(\lvert M_f\rvert\le2^{40}\). The saved reference mass is inside this range. The archived producer's internal scalar and floating-point rounding still need independent binding.

## Finite axial coefficients without replaying the recurrence

The degree-24 finite seed has \(U_0=4Z+j\) and \(P_0=-\epsilon P_\star^2 M_f/(1+Z^2)^2\). Its swirl seed is exactly zero, so the recurrence gives \(P_n=0\) for \(n\ge1\) and the angular coefficients drop out of the axial recurrence. In the complex capsule of radius \(1/4\) about the real interval \([-1,1]\), use \(\lvert Z\rvert\le2\), \(\lvert d\rvert\le5\), \(\lvert L\rvert\ge1/2\), and \(\lvert1+Z^2\rvert\ge1/2\). The first axial coefficient obeys
\[
\lvert U_1\rvert\le D,\qquad
D=2^{10}(1+2^{40})\,\epsilon P_\star^2<1.
\]

Reserve an axial loss \(1/200\) per radial step. In the \(n\)-th step, write \(M_0=10\) for the bound on \(U_0\), and \(M_i=D\) for positive indices. Cauchy gives \(\lvert U_i'\rvert\le200M_i\) on the next capsule. The original \(W_i,H_i\) weights then obey \(\lvert W_i\rvert\le1003M_i\) and \(\lvert H_i\rvert\le6M_i\). Every product for \(n\ge1\) has one positive index, so \(M_i M_{n-i}\le10D\). The normalized recurrence denominator has modulus at least \((n+1)^2\), and its numerator after dividing by \(\epsilon D\) is at most \(10(n+1)(1003n+1204)+1\). It follows that
\[
\lvert U_n\rvert\le 2^{14}\epsilon D\quad (2\le n\le24).
\]
The calculation uses this induction and the finite derivative weights \(n(n-1)\cdots(n-i+1)4^{n-i}\). No individual radial coefficients are regenerated. The final capsule has radius \(26/200\), allowing a \(1/8\) Cauchy radius for ordinary \(Z\) derivatives of the inlet jets. Ten source checks pin the recurrence and pressure seed used in this argument.

## Actual minus finite with the common inlet terms paired

In units \(\epsilon^2 h_b P_\star^2 S\), with \(S=s_c/2\), the exact-real nominal comparison has the form
\[
\mathcal M_{4,\mathrm{actual}}-\mathcal M_{4,\mathrm{nominal}}(M_f)
   =(M_2-M_f)\,\mathcal V_2(Z)+E(Z).
\]
Here \(\mathcal V_2\) is the fixed-\(\beta=2\) pressure-mass response. The two inlet first-coefficient expansions have the same hydrodynamic term, which cancels algebraically before enclosure. The radial-constant axis inertia is also shared, leaving its contribution weighted by the actual \(\rho^2-16\) difference. The finite non-axis inertia is controlled by the axial and own-mean polynomial bounds. The actual late pressure masses, nonlinear inlet, viscous radial change, nonlinear inertia and true cutoff swirl all retain separate terms.

In the stated unit, the whole-axis magnitudes of \(E,E_Z\) are bounded respectively by \(5.70\times10^{-102226522508642279}\) and \(9.11\times10^{-102226522508642278}\). The normalized high radial coefficient \(U_n/(\epsilon P_\star^2)\) contributes at most \(4.39\times10^{-408906090034569549}\) before the radial and ordinary axial derivative weights. The dominant remainder bound comes from the late pressure masses. It is a bound on a signed residual function, with every inherited error packet kept separately.

Using the earlier explicit pressure reference as \(M_f\) gives a complete comparison with that exact-real nominal model on 128 closed cells. The difference has strict value sign on 124 cells and strict ordinary derivative sign on 110. This result keeps the distinct finite reference, the actual mass interval and the binary archived coefficient separate. It does not certify a binary-rounding bound or identify the producer's internal mass with the explicit reference.

The [frozen receipt](../experiments/root_st073/original_finite_collar_inlet_bound_receipt.json.gz) exited zero, SHA-256 23f6477433c9898c1b529494600492809d9d1b58b76b5cdf2bf29acf4c1e1f5e. The 212 forward terms, earlier error budgets and independent terminal comparison remain in their existing roles. No radial recurrence or quadrature was replayed.


## Third contribution and the actual cutoff axial response

The actual lower-collar graph is the exact prefix of the cutoff graph used here. The raw all-word \(U_z^{(1)}\) value and ordinary first derivative are nodes [43328, 43329]. The corresponding cut response nodes [52641, 52642] are those same functions multiplied by the axial-variable-independent cutoff \(\eta\). The leading velocity uses its physical \(U,U_Z\) nodes [36205, 36206]. Hence the third contribution and its derivative are defined by the exact radial density
\[
2\epsilon h_b\,\rho(s)\,U_{\mathrm{lead}}(s,Z)\,\eta(s)\,U_z^{(1)}(s,Z),
\qquad 0\le s\le 3s_c/16.
\]
The derivative uses both terms of the product rule. The archival comparison unit is \(\epsilon^2 h_b(s_c/4)\), using \(u_{\rm cut}=\epsilon\eta U_z^{(1)}\) and \(ds=(s_c/4)\,dq\).

At the saved \(Z=371/1000\) query the previously integrated signed inlet model and all its error sectors were copied into the new, separate third-row receipt. An outward rescaling into the archival unit gives a value between 97.15 and 117.00 and an ordinary first derivative between 406.02 and 488.96. Both are strictly positive. The separate directed radii are retained, including the all-word tail. The width is primarily inherited from the signed coefficient enclosures. The new function definition is global in the real axial interval. The numeric bounds above are for the saved axial input only.

The [cutoff axial receipt](../experiments/root_st073/original_collar_cut_axial_binding_receipt.json.gz) exited zero, SHA-256 37373941750517804af4e099ea70d49ecd596492289d39ee7bf1d76e255f59a1. The next source tasks are to establish the archival fourth-row scalar and rounding relation, and to enclose the signed third-row raw response over all real \(Z\) for the forward-source comparison. The independent negative complete-pulse target and the separate reserved own-swirl contribution are each carried once.
