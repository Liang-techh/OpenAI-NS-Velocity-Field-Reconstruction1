# Whole-Z analytic core and physical axis field — 2026-10-02

The same admitted `.001*delta` source now supplies enclosures of the nonlinear analytic core through mixed order five, Cartesian velocity/pressure spatial derivatives through order four including the exact axis, the first fixed-physical-position time derivative, and local core kinetic-energy bounds. These are enclosures of the admitted analytic fixed point, not a recomputed coefficient table or a substitution of the explicit linear model for the nonlinear solution.

Core-to-inner-annulus interfaces, the full assembled field, admissible stress and independent flat remainder, actual higher-order coefficient recursion, and oscillatory stress correction remain incomplete. The previous finding that the original unlocalized source has infinite whole-space kinetic energy remains unchanged. Local core bounds do not establish terminal-time integrability or repair that global source-domain issue.

## Reproduce and provenance

Run `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage corephysical`. This stage runs the new core producer and checker; the full ordered pipeline now contains 80 modules.

The new pair is `lei_ren_part1_paper_compliant_core_physical_field.py` and `lei_ren_part1_paper_compliant_core_physical_field_check.py`, each with a JSON receipt. The checker binds transitive input hashes, the admitted Cartesian outer and physical-energy receipts, the analytic core transfer/uniform bounds/tube, and the physical norm family's selected `Cstar`.

- Actual five-defect family: `3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894`.
- Actual implicit source: `5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`.
- Analytic preheat datum: `dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010`.

Direct equality checks bind the physical norm family's analytic core hash to the transfer, uniform bounds and shared analytic tube. No old finite coefficient state is used or extrapolated. The original pressure datum and its derivatives through six supply the core forcing; pressure is not fitted afterward.

## Actual core profiles and error enclosures

Use paper Section 8 notation

\[
\rho=\Lambda R,\quad\epsilon=\Lambda^{-1},\quad
d=1-Z^2,\quad L=1-\delta Z^2,\quad U_0=4Z+j,
\]
\[
F=F_0\Phi,\qquad U_z=U_0+\epsilon\Psi,\qquad
F_0=\exp(-\log C_* -\Lambda G).
\]

Here `Cstar` is the selected physical norm family constant. With

\[
H_0=(1-\delta)Z/2+dU_0,\quad\sigma=j/500,
\quad G'=LH_0/(H_0^2+\sigma^2),
\quad\chi=H_0^2/(H_0^2+\sigma^2),
\]

`G` is anchored at the unique root of `H0` in `[-j,0]` and is nonnegative on real `[-1,1]`. The positive amplitude is retained in formal logarithmic form; no huge exponential is evaluated or replaced by zero.

The leading fixed-point pair is

\[
\Phi^{(0)}=B(\chi\rho/2),\qquad
B(q)=\sum_{n\ge0}\frac{(-q)^n}{n!(n+1)!},\qquad
\Psi^{(0)}=-\frac{\rho g}{2L},
\]
\[
g=-\left[(1+\delta)(1-2ZU_0)U_0/2+4H_0
+dP_{0,Z}-2(1+\delta)ZP_0\right].
\]

The degree-80 evaluation of the infinite Bessel model includes a factorial tail for every mixed derivative through five. It does not define a finite-polynomial field. A decreasing positive tail ratio bounds all omitted terms. The separate alternating-series positivity bound intersects the resulting interval enclosure.

Let `D` be the conservative admitted scaled map size divided by `1-q`, where `q<1` is its Lipschitz bound. The analytic embedding supplies

\[
E_{ik}=\frac{(i+k)!}{20^i h^k(k+1)^2
(1-\rho_{\max}/20)^{i+k+1}}.
\]

The nonlinear correction is bounded by `D*E_ik`; its zero-axis condition gives the sharper zeroth-radial bound `D*rho_max*E_1k`. Thus the interval represents the unique admitted nonlinear fixed point rather than just its leading model. The domain includes all real `Z in [-1,1]` and `rho in [0,4.1]`, with a short analytic extension beyond the inner exit at `rho=4`.

The radial mean `M(Psi)=integral_0^1 Psi(s*rho) ds` retains its correct averaging factors. Exact axis values and slopes are imposed from the source:

\[
\Phi(0)=1,\quad\Psi(0)=0,\quad
\Phi_\rho(0)=-(\chi+\epsilon\beta)/4,\quad
\Psi_\rho(0)=-g/(2L),
\quad\beta=(3-\delta/2+jZ)/L.
\]

## Radial recovery and original pressure

Let `M=U0+epsilon*M(Psi)`. The exact mean identity is `Uz=M+rho*M_rho`. Recover

\[
U_r=\sqrt{R/2}\,Q,\qquad
Q=\frac{2ZU_z-(1-\delta)ZM-dM_Z}{L}.
\]

Both mean terms have minus signs. At the axis,

\[
Q_0=\frac{(1+\delta)ZU_0-4d}{L}.
\]

Retain the same analytic preheat pressure and forward increment:

\[
P=P_*^2P_{0,\mathrm{norm}}
+\epsilon F_0^2\int_0^\rho\Phi(s,Z)^2\,ds.
\]

Bell derivative ratios for `F0` and `F0^2` are included. Differentiating only `Phi` would omit actual axial amplitude derivatives. The pressure increment is exactly zero at the axis.

## Physical Cartesian field without axis division

Use `X=sqrt(Lambda)*x/lambda`, `Y=sqrt(Lambda)*y/lambda`, so `rho=(X^2+Y^2)/2`. Then

\[
u_x=\sqrt\epsilon[\lambda^{-1}XQ/2
-\lambda^{-1-\delta}YF_0\Phi],
\]
\[
u_y=\sqrt\epsilon[\lambda^{-1}YQ/2
+\lambda^{-1-\delta}XF_0\Phi],\qquad
u_z=\lambda^{-1-\delta}U_z.
\]

This representation has no inverse-radius division, including at `X=Y=0`. Symbolic derivative templates retain polynomial basis derivatives, profile derivatives and distinct lambda powers for the radial and swirl contributions. All 35 spatial multiindices through total order four are enclosed for velocity and pressure, with first time derivatives evaluated at fixed physical `(x,y,z)`.

The exact radial mean recovery proves the symbolic core divergence identity. Independently enclosed component intervals need not numerically cancel; no whole assembled-field divergence or momentum-residual certificate is inferred from these core identities.

## Local energy scope

The physical measure and energy functional from the accepted outer stage are used. For bounds `q=sup|Q|`, `phi=sup|Phi|`, `v=sup|Uz|`, the core radial integrals obey

\[
I_r\le\epsilon^2\rho_{\max}^2q^2/4,\quad
I_\theta\le\epsilon^2\rho_{\max}^2\phi^2F_{0,\mathrm{upper}}^2,
\quad I_z\le\epsilon\rho_{\max}v^2.
\]

Original `Lambda/F0` logarithmic factors remain in the ledger. The three components have bounds on closed sectors `|Z|<=.5,.9,.99` at `log(tau)=-1,-10,-100`, on fixed physical `|z|<=1` at positive time, and on that strip over `exp(-100)<=tau<=exp(-1)`. The strip estimate behaves as `tau^-1` and proves neither a uniform terminal-time energy bound nor integrability to `tau=0`. Missing annuli prevent a full-background local-energy claim.

## Validation and next work

The accepted checker reports all passed: 1,050 actual mixed profile bounds, 50 exact axis conditions, 294 nonlinear correction comparisons, 490 actual Cartesian contribution bounds and 14 physical-time brackets. Independent moderate-parameter diagnostics cover 126 infinite-model mixed derivatives, 280 Cartesian spatial derivatives including the exact axis, eight time derivatives and three symbolic recovery identities. Those fixtures explicitly have `actual_source_admission=false`; actual source admission comes from the analytic receipts and bound comparisons.

Completed in this limited source-bound scope: whole-Z core profile enclosures through five, core/axis Cartesian spatial-four and time-one enclosures, core local positive-time energy bounds.

Next requirements:

1. Implement the original inner frozen/reference/join/matching/shear annuli with the same core, pressure and moment data; restore mixed derivatives sufficient for radial recovery and physical order four.
2. Certify core exit and every annular interface as functional identities or controlled source-bound joins, then connect to the accepted O3 field. Build an explicit regional dispatcher with exact axis treatment.
3. Extend physical divergence, local energy and terminal-time domain analysis to the entire assembled background.
4. Construct admissible divergence-form stress and independent flat remainder on every region, with cone margin and physical norms.
5. Implement actual n-dependent recursion and independent moment repair before oscillatory correction and full corrected residual validation.

Whole-field, stress-cone, full-energy and temporal-recursion gates remain false.
