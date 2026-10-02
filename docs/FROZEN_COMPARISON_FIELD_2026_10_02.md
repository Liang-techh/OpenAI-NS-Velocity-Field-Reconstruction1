# Actual core exit primitives and frozen comparison — 2026-10-02

The admitted nonlinear analytic core now supplies all five core-exit radial primitives with axial derivatives through five. Their original positive `F0` amplitude factors remain formal. The explicit frozen comparison in paper (9.12) is callable across `Ra=4/Lambda <= R <=110`, with mixed profile derivatives through total order four and axial derivatives through four of the normalized inertial direction.

This is the **unsmoothed auxiliary frozen comparison**, not the actual connecting velocity. Its radial derivatives generally differ from those of the core at the exit. No smooth core interface, actual prescribed-shear bridge, stress lift/cone or coefficient recursion is promoted by this stage.

## Reproduce

Run `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage frozenfield`. The complete ordered pipeline contains 82 modules. The stage executes the producer `lei_ren_part1_paper_compliant_frozen_comparison_field.py` and its independent checker, consuming the accepted core source and its checker receipt.

The family remains `3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894`, implicit source `5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`. Transitive hashes bind the original preheat pressure, selected `Cstar` and unique nonlinear core. No old coefficient table or conditional legacy bridge seeds the field.

## Five actual exit primitives

Write `epsilon=1/Lambda`, `r=4*epsilon`, `V=Uz`, `F=F0*Phi`, and define

\[
H=\tfrac18\int_0^4\rho\Phi\,d\rho,\quad
K=\tfrac18\int_0^4\rho\Phi V\,d\rho,\quad
A=\tfrac14\int_0^4 V^2\,d\rho,
\]
\[
B=\tfrac1{16}\int_0^4\rho\Phi^2\,d\rho,\qquad
C=\tfrac14\int_0^4\Phi^2\,d\rho.
\]

At the core exit, these give

\[
m_\theta=F_0r^2H,\quad m_z=r\,\overline V_c,\quad
m_{\theta z}=F_0r^2K,\quad
m_{z\theta}=rA-r^2F_0^2B,\quad m_p=rF_0^2C.
\]

The actual core's radial mean `Vbar_c` is inherited directly. Each integral and all its axial derivatives through five are enclosed using the admitted full-radial core bounds and the exact nonnegative integration masses `int_0^4 rho=8`, `int_0^4 1=4`. These are conservative integral enclosures, not recomputed point coefficients or a new definition using constant integrands. Squared zeroth coefficients use real interval squares to retain nonnegativity. Higher derivatives use product rules.

The receipt separates derivatives of normalized integral shapes from true moment derivatives divided by their **base-point** amplitude. Bell ratios of `F0` and `F0^2` are applied to the latter. This distinction matters when restoring the axial stress.

## Exact frozen moment transport

Let `phi=Phi(4,Z)`, `v=V(4,Z)`, and `theta=r/R`. With `F_f=F0*phi`, `V_f=v`, the five shapes are

\[
\frac{M_{\theta,f}}{F_0R^2}=\theta^2H+(1-\theta^2)\phi,
\qquad\frac{M_{z,f}}R=\theta\overline V_c+(1-\theta)v,
\]
\[
\frac{M_{\theta z,f}}{F_0R^2}=\theta^2K+(1-\theta^2)\phi v,
\]
\[
\frac{M_{z\theta,f}}R=
\underbrace{\theta A+(1-\theta)v^2}_{a}
-RF_0^2\underbrace{[\theta^2B+(1-\theta^2)\phi^2/2]}_{b},
\]
\[
\frac{M_{p,f}}{RF_0^2}=\theta C+(1-\theta)\phi^2=:p.
\]

These are exactly (9.12); core histories are retained at `theta=1`. Pressure remains `Pstar^2*P0_normalized+RF0^2*p`, using the original prescribed axis datum.

Set `d=1-Z^2`, `L=1-delta*Z^2`, `m=Mz/R`. Recover

\[
Q=\frac{2Zv-(1-\delta)Zm-dm_Z}{L},\qquad
U_r=\sqrt{R/2}\,Q.
\]

The exact mean identity `m+m_y=v` implies structural physical divergence zero for this auxiliary frozen field. It does not establish smooth matching to the core or divergence of the still missing full assembly.

## Inertial direction without inverse-amplitude evaluation

The original inertial formulas (9.13) yield `Df/R` from the inherited moments. The axial direction `Ef=Iz/F` contains an inverse `F0`, which cannot be materialized at the selected source scales. The producer instead supplies the exact drive

\[
\sqrt{R/2}\,F_f E_f
=R\,D_{\rm hydro}+RP_*^2\,D_{\rm pressure}
+R^2F_0^2\,D_{\rm swirl}.
\]

Its three coefficients retain all pressure/moment/amplitude derivatives. This quantity is what enters the axial bridge after multiplication by `F_actual/F_comparison`; it avoids replacing a positive amplitude by numerical zero. The true bridge uses the **smoothed** comparison's direction, so the unsmoothed frozen drive must not be passed directly into (9.26) as a shortcut.

For mixed `y=log(R/r)` derivatives, keep the base-point prefactors separate: `sqrt(R/2)` for `Ur`, `sqrt(2R)*F0` for `Utheta`, one for `Uz`, `Pstar^2` for axis pressure and `RF0^2` for its increment. The radial part of `Q` has a constant and a `theta` contribution, giving factors `(1/2)^i` and `(-1/2)^i`; all `i+k<=4` derivatives are retained.

## Independent evidence and next dependency

The checker independently derives all five primitive radial equations and their exit values, the mean recovery and structural divergence identities, and the two positive integration masses. On a separate moderate-parameter polynomial core it integrates all five core primitives and reconstructs inertial stress directly from (9.13). It compares 150 moment derivatives, 60 direction derivatives and 450 mixed profile derivatives. It also checks 525 actual source mixed profile enclosures for finiteness and completeness. Independent fixtures explicitly have `actual_source_admission=false`; source admission remains bound to the accepted nonlinear analytic core and interval integration argument.

Next implement (9.23), including the comparison's own five accumulated moments: it agrees with the core through `y<=hb`, smoothly freezes on `[hb,2hb]`, and is constant afterward. Then integrate the actual bridge (9.26), retaining the weighted flat error `1-chi_b`, and carry its **actual** moments through the short switches at `100`, reshape from `110`, axial restoration and five-moment patch. Only that construction can close the missing core-to-outer interfaces. Physical energy/stress/flat remainder and genuine recursion remain open.
