# Complete original n1 pressure-moment C1 accuracy — 2026-10-13

The complete original m13/m15 first-Z derivatives now include the inner cutoff sectors and every external radial segment from Rin to infinity. The checked 4096-cell grid covers every real Z in [-1,1], with maximum normalized absolute derivative errors approximately 0.0005117131891907913 and 0.0001210322320544479. Both are below .001 in the same fixed positive original pressure-moment units. These are pressure-moment derivative enclosures; actual signed defect function binding, selection of n1 alpha/beta functions, repaired own M/V/P, five terminal zeros and genuine temporal recursion remain open.

## Source improvement consumed by pressure

The original loop is kept unchanged. Write `chi_u(psi)=2h sum(r^n cos((n+1)psi))`, with `h=(1+u^2)^(-1/2)` and `r=u*h`. For `u=sinh(alpha)`, its exact Fourier covariance is `integral chi(u)chi(v)/(4pi)=sech(alpha-beta)`. Mixed derivatives give squared hyperbolic L2 norms 1, 1, 5 and 61. The genuine u derivative norms are bounded by `2sqrt(pi)*(1,1,4,17)`. Geometric Fourier sums supply the required angle/mixed derivative bounds.

Use the unnormalized inverse `F=psi+integral(t^2)-2pi*phi*nu=0`, `t=t0+q*chi`, `nu=1+t0^2+2q^2`, with `F_psi=1+t^2>=1`. This retains the normalization cancellation before estimating inverse derivatives. The original q-root/sigma derivatives and all source function handles remain bound. No magnitude cap is differentiated, no new N is chosen and no function is replaced by an interval midpoint.

The new A/B enclosures enter the actual five-density algebra, original own-rate transports and normalized five targets. The same implicit control branch then obeys `J*h_ZZZ=-d_ZZZ-3*D2Q[h_Z,h_ZZ]/N`, with its already accepted J inverse and C1/C2 bounds. The updated bound is consumed by compact-band V and mean derivatives and by the differentiated compact pressure identity. Updating the active core alone would have left the compact-band third derivative bound unusable.

The density bounds retain the accepted uniform `1/N0` envelope, a conservative bound for the much larger selected N. Actual physical mean/profile corrections are divided by the unchanged selected N once; these two bounds are not interchangeable physical parameters.

For the active pre-Rc mean, every preceding chart is retained. An actual whole-chart rate-zero dlogR width bounds each contained native subcell; the positive rate-one history kernel is bounded by one. This gives a safe bound for every free-query prefix, including a query inside a native cell. The genuine baseline m_ZZZ is added with one Pstar and the same selected N. The original physical quotient then yields q_ZZ. All 4 axial cells × 17 charts and 228 native cells are covered. The physical mean and implicit control archives have identical selected-N logarithms; the compact profile uses that same exact selected integer.

The differentiated reduced-S integrals keep coefficient derivatives and the products `U_Z*q_Z+U*q_ZZ`, `U_Z^2+U*U_ZZ`, and `2q*q_Z`. The accepted lower/upper finite endpoint derivatives remain signed. Exact shared Rc cancellation precedes pressure-unit arithmetic. Quiet-band Z independence remains exact; bounded active/compact contributions receive nonzero allowances.

## Complete moment composition and checks

The accepted unchanged collar provides its genuine full-domain first-derivative ranges in fixed Z-independent source units. Its positive scale ratio is the existing correlated physical ratio. Inner cutoff terms use `ordinary_Z_derivatives[1]` from the actual signed sectors, with their actual signs and nonzero physical source allowances. C0 affine errors are not differentiated. The checked outer pulse's signed center slopes and correct linear-in-halfwidth derivative Taylor budgets are reused.

Focused checks cover exact covariance derivative norms, a general symbolic third-order chain rule, direct third derivatives from the original analytic loop at signed u=-3,0,3, a nonunit finite-radius pressure fixture with delta=.2, actual archived-source replay, all 8192 pressure-moment rows, source/grid seams, signed lower contributions, exact source pins and accepted public queries at precision 15. The full Navier–Stokes residual and stress/recursion gates remain false.

The three new quartets in `experiments/root_st073` are:

- `lei_ren_part1_paper_compliant_current_original_n1_preRp_C1_inputs`: genuine baseline, correlated units, compact inputs and quiet source.
- `lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1`: Fourier bounds consumed through actual targets/controls/fields and upper pressure derivatives.
- `lei_ren_part1_paper_compliant_current_original_n1_complete_pressure_C1`: full two pressure-moment derivative composition.

Each has `.py`, `.json.gz`, `_check.py` and `_check.json`. Reproduction commands are `python scripts/run_original_n1_preRp_C1_inputs.py`, `python scripts/run_original_n1_preRp_pressure_C1.py`, and `python scripts/run_original_n1_complete_pressure_C1.py`. The first restores checked owners without executing ancestor producers/checkers. The latter two reuse checked archives and do not rebuild the original field graph. Only regenerate a quartet and affected descendants after changing its actual source.

The complete local gate is `current_original_n1_complete_pressure_moment_C1_accuracy_installed`. It certifies numerical first-Z derivative accuracy of the complete original m13/m15 moments, not terminal closure, repaired fields or a blow-up theorem.

## Next agent tasks

- [x] Recover lower/upper finite endpoint first derivatives.
- [x] Replace the excessive original primitive C3 bounds with the same-function Fourier/implicit estimate.
- [x] Propagate through five targets and the same implicit repair branch into compact C3 rows.
- [x] Retain every active mean predecessor and consume genuine m_ZZZ in q_ZZ.
- [x] Differentiate active/compact pressure and retain finite signed endpoints and quiet-band identities.
- [x] Compose genuine inner-cut and lower-collar derivatives with the outer pulse on all 4096 full-domain cells.
- [ ] Bind genuine signed C1 functions for all five full n1 defects to the installed reserved linear operator. The numerical enclosures and center slopes cannot serve as selected coefficient functions.
- [ ] Preserve the original inner cutoff, collar, active, compact, quiet and pulse signed radial integrals in those definitions, with one physical measure and the original m15 sign.
- [ ] Evaluate/select the actual alpha/beta functions from the original 2+3 linear moment repair, including finite mu and amplitude quotient derivatives.
- [ ] Construct the repaired family's own M/V/P and retained P0 datum; do not reuse the pre-repair pressure/history arrays.
- [ ] Establish all five terminal moment identities as functions across [-1,1], including zeros and joins.
- [ ] Complete source/axis/heat-exterior and stress constraints needed before genuine n2/temporal recursion.
- [ ] Implement order-dependent recursion, oscillatory stress cancellation and flat remainder, then independently validate the corrected Cartesian NS residual with forcing fixed/restricted in advance.

Full original goal **ACTIVE / INCOMPLETE**. The next work is actual function binding and repair, rather than another refinement of this accepted derivative grid.
