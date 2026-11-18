# ST073 fourth forward m15 source over the whole real axis — 2026-11-18

The actual fourth source now has directed enclosures for its value and ordinary first derivative on 128 closed axial cells covering \(Z\in[-1,1]\). It is the sum of the regular leading pressure and the actual cutoff swirl from the [canonical source receipt](CURRENT_ST073_REGULAR_COLLAR_PRESSURE_2026_11_18.md). All results below use \(U_4=\epsilon^2 h_b P_\star^2 S\), with \(S=s_c/2\).

## Inlet jet and signed terms

Write \(p=\Psi_{\mathrm{core},\rho}(4,Z)/P_\star^2\), \(d=1-Z^2\), and \(L=1-\delta Z^2\). The inlet regular viscous term in units \(U_4\) is
\[
V_0=16\bigl(d\,p_Z-(3+\delta)Zp\bigr)/L .
\]
The derivative uses \(p_{ZZ}\) and differentiates \(d\) and \(L\). Set \(x=1+Z^2\) and write \(I_k=\int(\beta)_k x^{-\beta-k}\,d\mu\) for the positive flatten measure. Its exact derivative is \(I_k'=-2ZI_{k+1}\). For \(k\leq3\), use \(0\le I_k\le (k+1)!\,\mu([0,2])\) and retain both rational pressure masses. The hydrodynamic inlet formula and the \(X_h\) fixed-point remainder each carry ordinary derivatives through order two. The latter uses the coefficient embedding at \(\rho=4\), with ordinary orders \((1,k)\).

For the radial viscous difference, \(T_s=h_b K\) gives \(|T-T_0|/P_\star^2\le h_b S\sup|K|/P_\star^2\) on the common complex tube. The difference uses the saved bound for \(K-2T_0\) and the positive cap for \(|1-\chi|\). One Cauchy loss controls the axial derivative of \(T-T_0\). A second controls the first derivative of the entire viscous difference. The signed radial-constant inertial term is evaluated from \(Q_{\rm axis}, U_{\rm axis}\) and multiplied by the positive average of \(\rho^2/(2P_\star^2)\). The remaining signed inertial contribution is bounded with the saved whole-collar \(S_1,S_2\) polynomials and their physical \(\epsilon,\epsilon^2\) factors.

For the cutoff product, the all-word raw \(F_1\) handle is the first component of the canonical regular system. The common analytic norm sums the degree 0–64 majorants and the certified omitted-sum bound. The same complex neighborhood contains the leading \(F_0\Phi\). Apply Cauchy to their product, using the true radial weight and support length \(3S/8\). This bounds the cutoff contribution and its derivative in \(U_4\) units. The receipt keeps separate directed radii for the viscous difference, nonlinear inertia and cutoff swirl. These radii bound unknown signed terms.

## Measurements and archived comparison

The value has strict sign on 124 cells (84 negative and 40 positive). The ordinary derivative has strict sign on 114 cells (57 of each sign). At \(Z=0.371\), the value and derivative round to \(-26.75578977231765\) and \(-537.9766234802223\) and the new enclosures intersect both earlier directed point enclosures.

The three legacy queries reuse the original q=2 six-row prefix of each unchanged 212-term \(R_p\) ledger. The old fourth coefficient uses \(\epsilon^2 h_b P_\star^2(s_c/4)\), so divide it by two before comparison. The four paired log enclosures intersect at each point. The old directed coefficient enclosure contains the new canonical enclosure for both orders at each point. The directed difference between the canonical function and the *stored finite approximation* gives the values below, rounded for display.

| \(Z\) | Order zero difference | Ordinary order one difference |
| --- | ---: | ---: |
| \(-0.3\) | \(+1.2174004\times10^{-19}\) | \(+4.9990168\times10^{-18}\) |
| \(0\) | \(+1.2231032\times10^{-18}\) | interval contains zero |
| \(0.3\) | \(+1.2174004\times10^{-19}\) | \(-4.9990168\times10^{-18}\) |

The differences measure the fourth row against its saved finite approximation at these inputs. Full function correspondence of the six forward contributions and independent terminal matching remain open. The original six rows, all 212 cumulative terms and their error packets, the once-only reserved own-swirl contribution and the independently defined negative full-pulse target are carried forward.

## Execution and next step

[Passing v2 receipt](../experiments/root_st073/original_collar_regular_whole_z_v2_receipt.json.gz), SHA-256 85dae8766190ffd39c4728d35c25d4d406143e1d794ed789feecc763e04b9861. [The initial attempt](../experiments/root_st073/original_collar_regular_whole_z_receipt.json.gz) stopped when converting an inertial log containing a positive power of \(h_b\) to an unbounded integer. The distinct v2 entry bounds extremely negative binary logarithms above by \(-10^{20}\) and retains the resulting positive dyadic radius. The accepted run exited zero after verifying the frozen input digests. No accepted calculation was replayed.

Use the shared regular viscous source and the actual cutoff product to compare the signed fourth-row difference **as a function** against the finite producer on the same geometry. Bind the remaining inner/lower and active/compact signed source histories before applying the terminal comparison.
