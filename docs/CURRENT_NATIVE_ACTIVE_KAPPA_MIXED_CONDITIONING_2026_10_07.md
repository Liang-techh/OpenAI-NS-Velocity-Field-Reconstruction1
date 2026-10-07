# Active kappa correlations, weighted periodic kernels and source-term attribution

> Successor: [CURRENT_NATIVE_CENTERED_PHASE_CONDITIONING_2026_10_07.md](CURRENT_NATIVE_CENTERED_PHASE_CONDITIONING_2026_10_07.md) ([4eb8b731](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4eb8b73150a190e05603ff006076541e0c7b1d30)) installs centered K/phase and original a*nu=v mixed A/B bounds, including full M inverse sensitivity. All five averaged caps improve again. Next priority: exact N-dependent source transport/value/target interface and actual controls; terminal closure remains open.

Checked implementation: [e6b266c9](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e6b266c9e698fce5972eb65c210c46ab8364d7f9). `NativeActiveKappaMixedConditioning.route(Z)` rebuilds the same original inlet-to-Rc route on **24 true cells /17 charts**, for Z=[-1,1] and Z=[.49,.51], with every integer N>=160. Support-only a/b correlations and joint q/kernel derivative bounds reduce all five averaged C1 target certificates in both domains. The maximum natural-log cap decreases from approximately `1.548972048e+408906090034569677` to `1.106919531e+408906090034569677`.

These are conservative logarithmic upper bounds. They are not actual field/error values, measured small residual or completion percentages. The new averaged certificate still does not beat the accepted direct certificate at N=160 for any target; retain the minimum of valid same-source certificates. Actual five controls, terminal closure and genuine scale recursion remain incomplete.

## Exact original active-source correlations

On q support, a+b^2/a<2+eta<=5/2, so a<=3 and b^2<=3a. The same original identities give

`Delta_i=(1-t0^2)*a_i-2*t0*b_i`,

`Delta_yZ=(1-t0^2)*a_yZ-2*t0*b_yZ+2*(b_y+t0*a_y)*(b_Z+t0*a_Z)/a`, with t0=-b/a.

The mixed cross term also equals -2*t0_Z*(b_y+t0*a_y) or -2*t0_y*(b_Z+t0*a_Z). Intersect these equivalent derivative covers with the native roots. This removes unrelated unrestricted b extremes from the q cutoff jets. The body, transition and flat branches retain their original meanings; the flat branch has exact zero modulation jets while the underlying velocity remains nonzero where originally defined. No derivative or source coefficient is selected from a cap.

## Joint weighted kernel derivative bounds

Let h=sqrt(1+u^2), r=u/h and theta be the original Mobius angle. Exact theta_u=2*sin(theta)/h and theta_uu=4*sin(theta)*cos(theta)/h^2-2*u*sin(theta)/h^3 give decay in the outer |r|>=1/2 sector. On that sector, P=(theta-psi)/(2u) gives weighted P derivatives. H=K/(4r^2), K=(2-3/h^2)*theta+psi/h^2+2r*sin(theta), and 0<=H<=pi give weighted H derivatives. On |r|<=1/2, h<=2/sqrt(3) converts the accepted Fourier bounds. The resulting conservative **all-signed-u, all-partial-psi** bounds are

`|P|<=4pi/h; |P_u|<=(80pi/3)/h^2; |P_uu|<=160pi/h^3; 0<=H<=pi; |H_u|<=60pi/h; |H_uu|<=800pi/h^2`.

Differentiate the joint functions F=q*P(chi*q,psi) and Q=q^2*H(chi*q,psi), chi=p2/dstar, before taking absolute bounds. Here Q is a kernel product, distinct from the zero-mean density antiderivative G_j. For example, F_q=P+u*P_u, F_qq=chi*(2*P_u+u*P_uu), Q_q=q*(2*H+u*H_u) and Q_qq=2*H+4*u*H_u+u^2*H_uu. Complete first/mixed q/chi chains use |u|/h^2<=1/2 and |u|/h<=1. They require no positive p2 lower and no division by q or chi.

The original T1=t0*psi+2F and T2=t0^2*psi+4t0*F+4Q, implicit inverse, A/B primitives and exact signed density graph remain the same. Every native derivative of chi is retained. The source route keeps original widths, inlet, phase, endpoints, incoming pressure memory, P0/P0_Z and Rc A/A_Z/mu. No actual control field is installed.

## Source attribution and remaining dominant term

All56 signed f1_y/f1_yZ product terms are exposed per active cell. The first bridge exports their separate contributions to each normalized terminal target, using true mass once and actual downstream transport. They are algebraic absolute budget contributions, not independent repaired history functions. Both queries have the following natural-log target caps at N=160; full directed records remain in the JSON artifacts.

| Z query | Target | Previous averaged log cap | New averaged log cap | Dominant slow term |
| --- | --- | --- | --- | --- |
| whole_Z | M | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:m:0 |
| whole_Z | D=(J-M)/mu | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:m:0 |
| whole_Z | I | 1.528943778e+408906090034569677 | 1.086891261e+408906090034569677 | Z:h:0 |
| whole_Z | S | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:e:0 |
| whole_Z | Cp | 1.548972048e+408906090034569677 | 1.106919531e+408906090034569677 | Z:p:0 |
| Z_interval | M | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:m:0 |
| Z_interval | D=(J-M)/mu | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:m:0 |
| Z_interval | I | 1.528943778e+408906090034569677 | 1.086891261e+408906090034569677 | Z:h:0 |
| Z_interval | S | 1.536096732e+408906090034569677 | 1.094044215e+408906090034569677 | Z:e:0 |
| Z_interval | Cp | 1.548972048e+408906090034569677 | 1.106919531e+408906090034569677 | Z:p:0 |

The dominant terms are B_yZ for M and (J-M)/mu, E*A_yZ for I, 2V*B_yZ for S, and E^2*A_yZ for Cp. Their common dominant inverse component is **normalized_curvature*L_y*L_Z**, with first-bridge log cap `1.252839779e+408906090034569677`. The fixed-angle mixed component is `1.100124225e+408906090034569677` and each slow cross component is about `1.071512412e+408906090034569677`. These are source-bound component magnitudes before the older certificate intersection. The next production task is a proved centered same-function L residual, followed by centered B or integrated slow-variation bounds.

## Focused evidence

- 960 C0/Z transport rows, 480 inherited rows, 40 exact quiet rows and 8 quiet pressure memory rows.
- 2576 ordinary signed leading derivative source terms, 112 weighted first-bridge terms and 40 normalized target order rows.
- 84 independent direct-kernel weighted comparisons and 28 signed angle derivative identities. Complete source product rules and joint F/Q parameter derivatives are checked symbolically. Uniform bounds rely on analytic sector proofs, not the finite reference set.
- 80 independent original scalar-loop C0/y/Z/yZ references cover signed-r, r crossing, transition and flat branches with nonzero underlying fields.
- Producer 67.766s; checker 101.032s on the warm original runtime; 1098 working/index dependency hashes. Read-only reviewer: **GPT-5.6 Luna / max**.

## Agent tasks and acceptance criteria

- [x] **COND6a signed source-term attribution:** expose all56 ordinary f1_y/f1_yZ product terms per active cell, their exact signs and source derivative indices. Weight the first-bridge terms with the true mass once, actual downstream decay and same Rc A/A_Z/mu. Report C0/Z separately and retain the divided-row joint numerator. These are algebraic absolute budget contributions, not independently installed history/control functions.
- [x] **COND6b1 active kappa correlation:** use Delta=a+b^2/a-2 and t0=-b/a to bound Delta_i=(1-t0^2)a_i-2t0*b_i and Delta_yZ=(1-t0^2)a_yZ-2t0*b_yZ+2*(b_y+t0*a_y)*(b_Z+t0*a_Z)/a. Intersect with original source caps only within q support; use b^2<=3a and original denominator positivity. Preserve the original flat cutoff branch and all background velocity/source terms.
- [x] **COND6b2 joint weighted kernels:** prove uniform signed-u partial-angle decay of P_u/P_uu/H_u/H_uu; differentiate qP(chi*q) and q^2H(chi*q) before bounding, with chi=p2/dstar. Use full ordinary first/mixed source chain and |u|/h^2<=1/2, |u|/h<=1. No q/chi/p2 division or positive p2 lower is needed. Intersect with accepted same-function bounds and rebuild the whole original route.
- [x] **COND6b3 correlated phase residual (function-domain bounds; see successor):** replace only the dominant first-bridge bound descendants with a proved same-source residual adapter. Use F=1+t^2, L_i=2pi*phi*nu_i-T2_i(psi)=phi*T2_i(2pi)-T2_i(psi), L_yZ analogously, and psi_i=L_i/F. Keep the complete implicit mixed recurrence, its F denominators, signed u and every cross term. Prove a whole-domain cap for the centered difference; writing the identity alone or taking midpoint samples is insufficient. Compare the resulting curvature*L_y*L_Z budget against the accepted cap.
- [x] **COND6b4 centered B primitive (function-domain bounds; see successor):** use R=2pi*t0*phi-hat(T1), B=E*a*R/(4pi). Recover ordinary R_y/R_Z/R_yZ and the full E*a*R product rule from the original functions. Preserve the exact flat q=0 cancellation and use no q division. Intersect with accepted original B caps; do not replace B by a chosen zero or drop the original b/p2 source.
- [ ] **COND6c weighted slow variation:** prove a smaller integral or total-variation bound for G_y/G_yZ on the original first-bridge [3sc/4,1] coordinate interval. Preserve the true microscopic width, ordinary derivatives, one fractional phase, endpoint terms and actual incoming histories. A sampled quadrature result is not a uniform function-domain bound.
- [ ] **COND4 actual partitions:** if the correlated residual needs domains, derive a no-gap body/transition/flat cover from the original signed packet on both required Z domains. Bind positive denominators and include endpoints. Apply each coordinate Jacobian once; sum actual weighted budgets. Keep the unsplit certificate as a competing valid bound.
- [ ] **COND5 exact seam cancellation:** establish equality of G/G_Z on all16 joins using the same physical source functions, inverse and global phase. Cancel only after an exact source-function identity. Include the initial flat collar; shared radii or neighboring boxes alone do not certify equality.
- [ ] **COND6d beat the best retained certificate:** compare the direct, preceding periodic and active-kappa/joint-kernel certificates at identical N over both complete Z domains. The new averaging cap still does not beat the direct cap at N=160 for any target. Tighten the dominant residual/curvature or weighted integral enough to improve the best retained bound, then obtain a useful five-row C1 repair budget. Keep logarithmic cap reduction separate from actual error and physical-field values.
- [ ] **CONTROL1 actual function controls:** solve B_exact(mu)*h+N*r+Q_exact(mu,h)/N=0 with the original five-bump/Gram functions, mu and Ac/S. Keep controls h distinct from -r and preserve the actual Rc inlet. No interval-midpoint controls, selected-frequency proxy or downstream zero reset.
- [ ] **CONTROL2 whole-Z terminal closure:** certify five terminal equalities and first derivatives for the actual repaired field throughout the required Z domain; carry it through original Rc->2Rc geometry, matching and right collar. Post-repair2Rc->Rb source/cone admission remains separate.
- [ ] **HIGH/OUTER/ENERGY:** original higher slow/fast/physical jets and C4 seams, one common finite N and whole-route cone margins, compatible analytic pressure, exact heat exterior, physical tail energy and flat remainder.
- [ ] **REC/WAVE/PHYS:** distinct n=1/n>=2 recovery equations and independent moment repairs, finite-order remainder/smooth sum, two actual stress-cancelling oscillatory pulse families and mean lift, corrected Cartesian NS and measured scale recursion/material winding.


Scoped gate: `current_original_native_active_kappa_mixed_source_conditioning_and_weighted_attribution_executed`. Active original kappa and joint-kernel mixed bounds, source-term attribution and new five-target certificates are implemented. Actual controls, terminal closure, point inverse/higher jets, global finite N/cone, heat/energy, coefficient recursion, oscillatory stress cancellation and full corrected NS remain open.
