# Physical Cartesian outer field and energy domains — 2026-10-02

The accepted SAME `.001` source now has physical Cartesian vector spatial derivatives through total order four and the first physical-time derivative on the local O3 power neighborhood, all six O4 pulse charts, and the complete postpulse chain through the unbounded Gamma exterior. The physical volume measure and kinetic-energy functional are restored. Complete postpulse energy is bounded on closed similarity sectors, at fixed positive time on a fixed physical axial strip, and on that strip over a time interval separated from the singular time.

The original, unlocalized paper field has **infinite whole-space physical kinetic energy**. This is a source/domain fact, not a numerical instability or a failed radial-tail estimate. No axial cutoff has been inserted to conceal it. Whole-field/core/axis smoothness, full-background local energy, terminal-time energy, admissible stress, flat remainder, genuine coefficient recursion and oscillatory corrections remain open.

## 1. Reproduce

Run the existing Python entry point with `--stage physicalfield`:

```text
experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage physicalfield
```

This runs the Cartesian producer/checker and physical-energy producer/checker in order, consuming accepted source-bound receipts. The complete pipeline contains 78 modules. It does not rerun unrelated experiments or change the old `.01` family.

New producer/checker pairs:

- `lei_ren_part1_paper_compliant_cartesian_field.py` / `_check.py`
- `lei_ren_part1_paper_compliant_physical_energy.py` / `_check.py`

The family is `3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894`; the implicit source is `5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`. Input hashes bind the actual pressure/core/outer/heat family and accepted interfaces.

## 2. Cartesian vector derivatives

Let `y_profile=log R`, `c=cos(theta)`, `s=sin(theta)`. Use the paper's map

\[
R=\frac{r^2}{2\lambda^2},\quad Z=\frac{z}{\lambda^{1-\delta}},\quad
\lambda^2(1-Z^2)=\tau=1-t_{\rm phys}.
\]

Physical component exponents are `beta=-1` for `ur`, `beta=-1-delta` for `utheta/uz`, and `beta=-2-2delta` for pressure. They must remain separate even in the same Cartesian output component.

The vector transform is

\[
u_x=c u_r-s u_\theta,\qquad u_y=s u_r+c u_\theta,\qquad u_z=u_z.
\]

The exact differential operators are

\[
D_x=cD_r-\frac{s}{r}D_\theta,\qquad
D_y=sD_r+\frac{c}{r}D_\theta,\qquad
D_\theta=-s\partial_c+c\partial_s.
\]

Templates retain every basis derivative and inverse-radius factor. For a term `C(c,s) r^-q dr^a dz^b(lambda^beta G)`, with `N=a+q=i+j`, the exact common scale is

\[
\lambda^{\beta-N+b(\delta-1)}R^{-N/2}2^{(a-q)/2}H_{ab}.
\]

The factor of two stays inside the bracket. Contributions from `ur` and `utheta` have different lambda exponents and cannot be combined under one prefactor. All 35 spatial multiindices `i+j+b<=4`, all four outputs, and all theta are covered in each of 16 admitted charts. Enclosing `c,s` independently in `[-1,1]` is conservative; the symbolic divergence check separately uses `c^2+s^2=1`.

The Cartesian/cylindrical divergence identity is preserved:

\[
\partial_xu_x+\partial_yu_y+\partial_zu_z
=\partial_ru_r+u_r/r+\partial_zu_z.
\]

This identity does not itself certify the missing core/axis assembler. No division by an interval containing `r=0` is allowed.

### Actual source scales

The O4 radial grid already differentiates the `sqrt(R)` dependence. Its common physical profile normalization is

\[
U_r=\sqrt{R_p/2}\,P_*e^{-\mu\log(R/R_p)}\,G_r.
\]

Multiplying it by `sqrt(R)` again would be wrong. Angular/axial pulse grids retain `Pstar*exp(-(.5+mu)*log(R/Rp))`; pressure retains `Pstar^2`. The local O3 chart covers only `log(R/Rp) in[-1,0]`, not the entire preceding construction.

After the pulse,

\[
R_v=R_pe^{13/\mu},\qquad E_{v0}=P_*Ue^{-13/(2\mu)-13}.
\]

The `-26` in the squared amplitude is **not** an extra `+26` in `log(Rv/Rp)`. Postpulse velocity grids use `U/Ev0`; their radial and axial histories are exactly zero. Pressure remains `P/Pstar^2` with the original analytic datum and forward history.

Minimum postpulse log-radius offsets relative to `Rv` are: flatten `0`, following power `100`, angular `100+Lrel-4`, steep entry `100+Lrel`, steep power `101+Lrel`, steep exit `101+Lrel+Ts`, waiting `102+Lrel+Ts`, collar `102+Lrel+Ts+wait`, Gamma `105+Lrel+Ts+wait`. The Gamma upper radius is unbounded.

Huge absolute scales stay as exact finite logarithmic parts. Numerical positive caps never define a replacement field. The three physical-time sectors use `log(tau)=-1,-10,-100` and `lambda>=sqrt(tau)` with proved negative derivative exponents.

## 3. Physical time derivative

For the stationary leading profile, at fixed physical `(x,y,z)`, paper Lemma 2.1 gives

\[
\partial_{t_{\rm phys}}(\lambda^\beta G)
=\frac{\lambda^{\beta-2}}{1-\delta Z^2}
\left[-\frac\beta2G+\frac{1-\delta}{2}ZG_Z+G_{y_{\rm profile}}\right].
\]

This is not a derivative of a radial stage parameter. The Cartesian basis is constant when time varies at fixed physical position. Future explicitly time-dependent profiles require their additional fixed-profile-coordinate time derivative; this leading map does not implement coefficient recursion.

## 4. Physical volume and kinetic energy

At fixed positive `tau`, set `d=1-Z^2`, `L=1-delta*Z^2`. The exact measure is

\[
dV=2\pi\lambda^{3-\delta}\frac Ld\,dR\,dZ.
\]

Thus, including the usual factor `1/2`,

\[
E_{\rm kin}=\pi\int_{-1}^1\int_0^\infty\frac Ld
\left[\lambda^{1-\delta}U_r^2+
\lambda^{1-3\delta}(U_\theta^2+U_z^2)\right]dR\,dZ.
\]

Writing `Ir=int Ur^2 dR` and `Iperp=int(Utheta^2+Uz^2)dR`, their weights are

\[
\tau^{(1-\delta)/2}L d^{-(3-\delta)/2},\qquad
\tau^{(1-3\delta)/2}L d^{-(3-3\delta)/2}.
\]

The signed profile moment `Mztheta=int(Uz^2-Utheta^2/2)dR` and centrifugal pressure moment `Mp=int Utheta^2/(2R)dR` are not this energy. A positive normalized remaining swirl integral is also insufficient to control physical volume energy.

### Complete postpulse local energy

The accepted flatten inlet contains the entire infinite future radial swirl integral `J(Z)` in `Rv*Ev0^2` units. `Ur=Uz=0` on this whole postpulse chain. The physical radial mass scale is

\[
R_vE_{v0}^2=R_pP_*^2U^2e^{-26}.
\]

The inverse-mu terms cancel analytically **before** enclosure. The new receipt uses the original whole-Z positive lower/upper bound on `J`; it does not subtract independently bounded long past integrals.

For each `Zstar=.5,.9,.99`, fixed-time energy bounds are recorded on `R>=Rv, |Z|<=Zstar`, including the actual physical axial halfwidth. These domains shrink in physical z as `tau` shrinks.

For the fixed physical strip `|z|<=1`, `0<tau<=1`, and original `delta<1/200`, the implicit map implies `lambda^2<3` and `d>=tau/3`. If `alpha=(3-3delta)/2`, the complete postpulse energy obeys the conservative bound

\[
E_{\rm postpulse,strip}\le
2\pi 3^\alpha\tau^{-1}(R_vE_{v0}^2)\sup_Z J(Z).
\]

It is finite at every fixed positive time. Integrating over `exp(-100)<=tau<=exp(-1)` gives the same constant times `99`. This bound does **not** prove uniform energy or time integrability through `tau=0`, and it does not include the missing core and preceding annuli.

### Whole-space obstruction from the original Gamma tail

The exact exterior is

\[
U_\theta=c_\infty R^{-(1+\delta)/2}
H_\delta(2d/R),\qquad R\ge R_b,
\]

where `Hdelta` is the full positive Gamma expectation, `Hdelta(0)=1`, and `Rb=Rtail*exp(3)`. Its radial energy lies between

\[
\frac{c_\infty^2R_b^{-\delta}}{\delta}
H_\delta(2/R_b)^2
\quad\text{and}\quad
\frac{c_\infty^2R_b^{-\delta}}{\delta}.
\]

The source-bound lower factor is strictly positive uniformly in Z. The physical endpoint weight is `d^-alpha` with `alpha=(3-3delta)/2>1`. Consequently, the full Z integral diverges at `Z=+-1`, which represent physical axial infinity, at every fixed positive tau. The original whole-space kinetic energy is infinite already from swirl.

In physical coordinates the heat swirl has the useful exact form

\[
u_\theta=c_\infty(\sqrt2/r)^{1+\delta}
H_\delta(4\tau/r^2)
\]

on the mapped exterior region. Its amplitude and heat argument no longer depend explicitly on Z. The region boundary still depends on the physical map.

Theorem 1.1 makes no global kinetic-energy claim. Equations (1.2) impose radial profile moments. The scale/radial cutoffs in (13.1)-(13.2) and (16.7) retain the leading heat tail. Endpoint comparisons (17.28)-(17.29) are restricted to interior Z sectors. An axially localized finite-global-energy variant would require new divergence-preserving potential/curl cutoffs and corresponding pressure, moment, stress and residual repairs; no such variant is claimed here.

## 5. Evidence and remaining work

The Cartesian receipt checks 3,360 actual contribution bounds and 64 first physical-time brackets. Independent checks cover 140 Cartesian spatial derivatives and four physical-time derivatives by differentiating a nonconstant three-component field after an independent positive implicit lambda root. Moderate-parameter checks are fixtures, distinct from the actual-source interval admission. There are 46 independent basis/commutator/divergence identities and 240 original normalization/scale comparisons. A read-only GPT-5.6 Luna/max audit found no material map, source-normalization, radius-origin or scope error; it independently checked the pulse radial source's hidden square-root cancellation. That source-specific identity is established algebraically rather than by the generic physical fixture alone.

The energy receipt checks 24 actual local energy/source inequalities, five symbolic Jacobian/mass/heat identities, an independent nonzero three-component kinetic-energy integral in physical coordinates, three volume points and 15 positive-Gamma floor points. Actual original Gamma and complete future energy remain source-bound. A read-only GPT-5.6 Luna/max derivation identified the source's local/global domain distinction.

Next work must restore the common-family core, near-axis limits, inner/matching annuli and remaining preceding O3 ranges in physical variables, then include their radial/axial/swirl contributions in local energy. Specify terminal-time and spacetime domains explicitly. Continue to the actual divergence-form stress cone and independent flat remainder, followed by genuine `n=1` and `n>=2` recovery and oscillatory correction. None of those stages is completed by this derivative/energy map.
