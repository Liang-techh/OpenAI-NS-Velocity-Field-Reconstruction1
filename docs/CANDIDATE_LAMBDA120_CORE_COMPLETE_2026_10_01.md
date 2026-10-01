# Lambda120 core complete and functional pressure input

The separate Lambda120 core is complete through radial order124, with128
initial axial Taylor coefficients and four coefficients in the final row.
Its physical pressure datum and anchored amplitude are unchanged during
recurrence. All20 normalized mixed C3 budgets pass on `s in [0,4.1]` at
`Z=.3`. The worst combined midpoint errors are:

- Phi: `1.84246060226628271e-23`.
- Psi: `7.51073122007974972e-13`.

These combine directed finite coefficient uncertainty and the infinite radial
tail relative to the accepted analytic datum. They are not the full physical
Navier–Stokes residual and do not certify a finite field over the whole axis.

The pipeline also generated separate Lambda120 receipts for the angular
identity, shared five core moments, recovered physical pressure and inlet
stress jets. The angular trace `Ttheta/F` enclosure is approximately
`+/-1.75611258335390762e-135`; the normalized axial stress enclosure is
`+/-1.50196320705033601e-16`. Containment of the analytic core's zero stress
is a consistency check, not annular or terminal moment matching.

## Pressure as an axial function

`candidate_pressure_function.py` now evaluates directed ordinary pressure
Taylor coefficients at any real center, or interval of centers, in `[-1,1]`.
It preserves all fourteen accepted radial mass contributions and validates
their transitive receipt hashes. No new pressure fit or tail term is added.

Fixed beta2 pressure is `M_beta2*(1+Z^2)^(-2)`; fixed beta0 pressure is constant.
For the flatten contribution, beta lies in `[0,2]`, so alpha=beta/2 lies in
`[0,1]`. Factoring `1+alpha*z^2` at its two imaginary roots gives, on a disk of
radius rho about any real center c,

`|1+alpha*z^2| >= (sqrt(1+alpha*c^2)-sqrt(alpha)*rho)^2 >= (1-rho)^2`.

With rho=1/4, the positive accepted flatten mass M supplies ordinary Taylor
coefficient bounds `M*(1-rho)^(-4)*rho^(-k)`. At order0 the contribution is
enclosed by `[0,M]`. The source profile is not replaced by the worst-case bound;
this is its analytic error enclosure. The functional jet routine reproduces
overlapping enclosures for all128 existing center-pressure coefficients and
passes17 exact coefficients of `(1+Z^2)^(-2)` at zero.

Run:

```powershell
python experiments/root_st073/lei_ren_part1_paper_candidate_pressure_function.py
```

## Next actual matching dependency

Lei–Ren (2.22) uses the five cumulative moments

`Mtheta=2 int rho F`, `Mz=int Uz`, `Mthetaz=2 int rho F Uz`,
`Mztheta=int (Uz^2-rho F^2)`, and `Mp=int F^2`.

Pressure restoration is the prescribed identity `P=P0+Mp` (4.26).
The terminal goals (1.2) are `Mz=Mthetaz=Mztheta=0` at infinity,
`Mp(infinity,Z)=-P0(Z)`, and the specified power-law asymptotic of Mtheta.
The computed core inlet supplies initial data for these functional conditions;
it does not supply their terminal values.

The next implementation must accumulate the actual transition moments up to
Rh and form the five reference defects (9.39). Section10 then supplies three
angular and two axial bumps with the nonlinear map (10.8),
`A*h+Q_Z(h,h)=-d(Z)`. Once all five moments agree with the reference at Rh,
the common outer continuation preserves their equality for every R>=Rh.
The pressure datum, waiting/heat schedule and outer profile must remain the
same construction data throughout this operation.

The paper fixes P0 before selecting Lambda and Cstar (Section8.4 and the
parameter ordering before Section9). Changing Lambda alone therefore does not
require fitting a new P0. Changing the outer profile generally does require
rebuilding P0 and its dependent choices. The present receipts prove the
accepted stored-datum scope, not original parameter or derivation uncertainty.

Still open: whole-axis finite core/inlet generation; actual transition and
five-bump repair; final outer five-moment closure; coherent heat collar;
admissible cone; n-dependent temporal recursion; oscillatory corrections;
independent full corrected Cartesian residual.

## Authoritative reference endpoint parameters

Use `coherent_pressure_source_alignment.json`'s `accepted_schedule`:
`logPstar=14`, `delta=1e-200`, and
`logRref=5e152+144.7004803657924162280799350326494935074...`.
The older `reference_moments.json` values `delta=1e-32`, `logRref=10` are
obsolete diagnostics. The accepted delta agrees with the Lambda120 core.

The physical reference endpoint is `Rh=exp(logRref-5)`; the paper and
`outer.py` use the relative log radius `y_h=-5`. It is not the core's
scaled radius `s=4`, which corresponds to `R=4e-120` for Lambda120.

At Rh, (9.3) supplies the analytic target functions:

`uh=exp(13.5)/(1+Z^2)`;
`Mtheta=(5/8)*Rh*sqrt(2Rh)*uh`;
`Mz=4Z*Rh`;
`Mthetaz=4Z*Mtheta`;
`Mztheta=16Z^2*Rh-(5/12)*Rh*uh^2`;
`Mp=(5/2)*uh^2`.

Evaluate these using factored or signed-log quantities. The existing joined
outer receipt only stores complete five moments at Rh and Z=.3; subsequent
rows do not provide the full five-moment profile. Thus the analytic Rh target
can be implemented immediately, while full transition/terminal outer moment
data must still be generated. The functional pressure receipt cannot substitute
for that missing coherent velocity/moment profile.
