# Smoothed comparison and actual inner bridge axial enclosures — 2026-10-02

The same admitted source now has callable **source-bound interval enclosures** of the smoothed comparison in (9.23), its own five cumulative moments, and the actual prescribed-shear bridge in (9.26) over `Ra=4/Lambda <=R<=100`, for every real `Z in[-1,1]`. The comparison has axial order six and its factored inertial direction has axial order five. The actual bridge has axial order five for regular velocity/moment/pressure inputs; radial velocity recovery has axial order four.

These are enclosures of the original analytic integrations. They are not recomputed point coefficients or a precision point evaluator. Radial/phase mixed derivatives, functional high-order core interfaces, the switches at `100..110`, and the remaining matching chain are incomplete. Full Cartesian field, physical energy/stress lift, actual temporal coefficient recursion and oscillatory correction remain incomplete.

## Reproduce and source binding

Run `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage bridgeprofiles`. The complete ordered pipeline has 84 modules. The new pair is `lei_ren_part1_paper_compliant_inner_bridge_profiles.py` and its `_check.py`, each with its own JSON receipt.

The producer consumes the accepted nonlinear core/frozen comparison, `compliant_K1_ledger` and `compliant_global_exit_certificate`. Transitive hashes and direct equality checks bind the selected `Cstar`, base analytic core, implicit source and exact inner parameter family. The family and implicit source remain those of the accepted `.001*delta` construction. Original `P0` and core moments are unchanged.

## The width is formal, positive and unchanged

The exact source is

\[
h_b=\epsilon_b=c_*K^{-100},\qquad
\log h_b=\log c_*-100\log K,
\]

where `K` is the actual physical norm sum (9.16), not an independently chosen convenient width. The admitted ledger supplies a directed log enclosure. Its upper bound proves `hb<1e-180`; the numerical interval `[0,1e-180]` is only an enclosure. No source formula is defined with that cap, and neither `hb`, `F0` nor `exp(log Cstar)` is materialized.

For each width-weighted source factor `A`, first combine

\[
\log h_b+\log|A|+\text{retained amplitude log}.
\]

Only after proving this upper bound smaller than `log(1e-180)` is its contribution enclosed by a positive cap. In particular, pressure contributions use `log hb+2log Pstar`; swirl contributions use `log hb-2log Cstar` with the original relative amplitude derivatives. Multiplying an arbitrary numerical width cap by a huge pressure would lose this essential cancellation. The receipt stores the actual log operands and each proof for independent recomputation.

## Smoothed comparison and its own moments

Let `y=log(R/Ra)`, `s=y/hb` and `alpha=1-sigma(s-1)`. Retain exactly

\[
\partial_y\log\bar F=\alpha\partial_y\log F_c,
\qquad\partial_y\bar V=\alpha\partial_yV_c,
\quad (\bar F,\bar V)(0)=(f,v).
\]

The comparison equals the exact core for `y<=hb`, smoothly freezes on `[hb,2hb]`, and is constant afterward. It is not the unsmoothed frozen field from the previous stage. Integration by parts gives

\[
\bar g(y)=\alpha(y)g(y)+\int_0^y[-\alpha'(t)]g(t)\,dt,
\qquad g=\log F_c\text{ or }V_c.
\]

The weights are nonnegative, independent of `Z`, and sum to one. The admitted core continuation `rho in[4,4.1]` therefore encloses the comparison and its axial derivatives without evaluating the microscopic width. The radial-constant `log F0` factors out exactly; exponentiating the convex enclosure of `log Phi` restores `Fbar/F0`. The sixth core axial derivatives use the admitted all-order `Xh` embedding, rather than extrapolating old local coefficient tables.

All five comparison moments continue from the same core exit. The producer uses the positive radial integration weights to enclose their actual histories. With `theta=Ra/R`, normalized shapes have weights `theta` and `1-theta` for axial/pressure integrals, and `theta^2` and `1-theta^2` for angular integrals. Integrands are the varying comparison fields; using an interval covering their range **does not replace them by a constant frozen profile**. `Mztheta` retains separate axial and `F0^2` swirl parts.

The original (9.13) recovers the comparison direction using its own moments. Bell ratios retain all axial `F0/F0^2` derivatives. The axial drive is factored as

\[
\sqrt{R/2}\,\bar F\bar E
=R\,\mathrm{hydro}+RP_*^2\,\mathrm{pressure}
+R^2F_0^2\,\mathrm{swirl}.
\]

## Actual prescribed-shear bridge

Use the exact source multiplier

\[
\chi=1-(1-h_b)\sigma(y/h_b).
\]

The actual field remains

\[
F=f\exp\left[-\tfrac12\int_0^y\chi\bar D\,dt\right],
\]
\[
V=v-\int_0^y\chi\frac{\phi_{\rm actual}}{\bar\phi}
\left(R\,\mathrm{hydro}+RP_*^2\,\mathrm{pressure}
+R^2F_0^2\,\mathrm{swirl}\right)dt.
\]

Here `F0` cancels in `Factual/Fbar=phi_actual/barphi`. The source comparison direction is used, not the earlier unsmoothed frozen direction. The exact inequality

\[
\int_0^y\chi(t)dt\le h_b(1+y)
\]

follows from the initial interval of length at most `hb` and the constant value `chi=hb` afterward. Together with `R<=100`, it gives derivative enclosures through axial order five, with width/pressure/gradient factors combined in logs before numerical enclosure. Comparison positivity comes from the inherited same-source analytic certificate; the actual swirl remains strictly positive.

The **actual** five moment histories are integrated from the core using these actual field ranges. Recover

\[
U_r=\sqrt{R/2}\,Q,\qquad
Q=\frac{2ZV-(1-\delta)Z(M_z/R)-d\partial_Z(M_z/R)}{L},
\]

and retain `P=P0+Mp`. At `R=Ra`, the log correction is exactly zero and all original moments are inherited. Structural divergence follows from the exact cumulative identity `d_R Mz=V`; it is not inferred by cancellation of independently enclosed intervals. Radial mixed derivatives and a numerical whole-field residual certificate are not claimed.

## Evidence and remaining work

The checker independently recomputes 142 actual source width-product log inequalities and checks 245 complete finite bridge profile enclosures, positive swirl and the unchanged exit. Separate moderate-parameter fixtures check 168 smoothed-comparison axial derivatives, six integrated-chi inequalities, 252 directly integrated variable-profile moment derivatives and 72 inertial-direction derivatives. Two symbolic identities verify cumulative-mean divergence and axial amplitude cancellation. Fixture receipts explicitly have `actual_source_admission=false`.

Next restore phase/radial derivatives with formal inverse-`hb` factors and the original flat sigma endpoint logic; then establish the full core join, integrate both original short switches at `R=100`, and obtain the actual `R=110` inlet. Continue that inlet through reshape/reference/axial restoration/moment repair and the missing preceding O3 region before whole physical-field assembly. Admissible stress/independent flat remainder and real `n`-dependent recursion remain separate subsequent tasks. The original unlocalized global-energy obstruction remains unchanged.
