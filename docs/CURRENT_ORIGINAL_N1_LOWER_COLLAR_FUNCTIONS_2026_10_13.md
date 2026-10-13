# Actual original n1 lower-collar functions

The Rin-to-r_minus collar now continues the admitted nonlinear core with
actual prescribed-shear integral bodies. Its own mean and swirl pressure
retain the complete history from rho=4. Both signed n1 pressure moments,
m13 and m15, have genuine value and first-Z functions, using the original
physical radius and both viscous radial derivatives.

This 2026-10-13 milestone is a function interface. The full original goal
remains **ACTIVE / INCOMPLETE**. Numerical point accuracy for this extremely
small original collar is not claimed. The solenoidal inner-cutoff sectors,
complete five defects and actual n1 repair still precede global matching,
heat/stress, true temporal recursion and oscillatory cancellation.

## Restore and consume

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gate is
`current_original_n1_Rin_to_rminus_true_lower_collar_C1_functions_installed`.

```powershell
python scripts/run_original_n1_lower_collar_functions.py
```

`CurrentOriginalN1LowerCollarFunctions()` restores the accepted true-core
and complete-pulse owners, verifies exact source hashes and constructs
the same saved graph. `.fields` supplies Utheta, U, m, Q, Ur, base
pressure and radial derivative rows. `.functions` supplies m13/m15 C1.
`.pressure_bindings` exposes the actual N1p/S and signed integrands.
`.parameters` and `.geometry_certificate` retain the original scales
and support proof. `.assert_graph()` protects accepted prefixes.

The common graph has 36,594 nodes. The actual core component nodes and
their ordinary mixed derivatives come from the checked analytic fixed
point. Older model-plus-error rows and comparison/cap centers are not
selected as these collar functions. Saved ancestors need no rebuild.

## Exact geometry

```
rho(s) = 4*exp(hb*s)
R(s) = epsilon_core*rho(s), epsilon_core=1/Lambda
Ra = 4/Lambda
Rin = Ra*exp(hb*sc/4)
r_minus = Ra*exp(hb*sc/2)
sc/4 <= s <= sc/2
dR = epsilon_core*hb*rho(s) ds
R*dR = epsilon_core^2*hb*rho(s)^2 ds
```

The phase endpoints are independent of Z. They are not the physical
radii. The source's selected sc is restored from its exact binary
mantissa/exponent tuple, not from the center of a numerical interval.
Its enormous negative exponent remains an exact `exp(exponent*log(2))`
expression, so no astronomically large power-of-two denominator is
allocated and ordinary floating underflow cannot silently select zero.

The bridge width is the same hbB node used by the accepted outer source,
with the original `hb=cstar*K^-100` definition. The original support
selector is executed at sc/2 and confirms that later modification
increments vanish throughout this collar. That zero support concerns
later increments; it does not erase the leading velocity or its history.

## True continuation and cumulative history

Let Phi_c, Psi_c and meanPsi_c be the actual normalized core. Set
`A_c=Phi_c,rho/Phi_c` and `chi(s)=1+(hb-1)*sigma(s)` using the genuine
original flat step. The installed bodies are

```
log(ratio)(s) = integral_0^s (chi(v)-1)*A_c(rho(v),Z)*hb*rho(v) dv
Phi = ratio*Phi_c
Psi = Psi_c(4,Z)+integral_0^s ratio(t)*chi(t)*Psi_c,rho(rho(t),Z)*hb*rho(t) dt
meanPsi = [4*meanPsi_c(4,Z)+integral_0^s Psi(w,Z)*hb*rho(w) dw]/rho(s)
U = 4Z+j+epsilon_core*Psi
m = 4Z+j+epsilon_core*meanPsi
Q = [2Z*U-(1-delta)Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2)
```

The rho derivatives retain
`chi_rho=(hb-1)*sigma_prime(s)/(hb*rho)` and the complete ratio terms.
Mean derivatives use `rho*meanPsi_rho=Psi-meanPsi` and
`rho*meanPsi_rhorho=Psi_rho-2meanPsi_rho`.
All finite histories begin at the true core exit, not at Rin, a source
cell or a pressure chart boundary. Independent dummy variables preserve
the nested integral scopes. Ordinary Z derivatives through order three
commute through these fixed-endpoint histories.

The absolute base pressure is
`P0+epsilon_core*F0^2*(integral_0^4 Phi_c^2 drho + integral_0^s Phi^2*hb*rho ds)`.
It retains the same native P0, anchored F0 and core swirl primitive.
This base pressure is distinct from the n1 pressure forcing moments.
U, m and Q receive no extra Pstar factor.

## Signed n1 pressure

For y=log R, the true radial rows recover
`Q_y=rho*Q_rho` and `Q_yy=rho*Q_rho+rho^2*Q_rhorho`.
The direct original operator is

```
S = [Q+(1-delta)*Z*Q_Z/2+Q_y]/L
    +Q*(Q/2+Q_y)+U*(-2Z*Q+(1-Z^2)*Q_Z-2Z*Q_y)/L
N1p = -S/2+(Q_yy+Q_y)/R
m13 = integral N1p dR
m15 = -integral R*N1p dR
```

The m15 minus sign belongs to the exported function already. Sum it
directly with the other signed m15 pieces. The direct differential form
retains both finite viscous endpoint contributions without subtracting
independently boxed endpoints. Because R and the phase endpoints are
Z-independent, the first-Z integrals require no endpoint-motion term.

## Focused evidence and limits

An independent exponential-core fixture with nonconstant phase sigma
integrates the logarithmic ODE and checks Phi through two rho derivatives
and three ordinary Z derivatives, including nonzero chi_rho. A separate
nonunit epsilon/hb/delta fixture checks the exact own-mean history,
ordinary C3 fields, Q through C2 and the core-plus-collar base pressure.
Its two pressure C1 functions must equal independently integrated
physical-R differential forcing, including both viscous derivatives.

Actual-source checks require the accepted unique analytic core nodes,
original exact binary sc, the positive physical measures, fixed axial
endpoints and executed zero-increment support selector. Source hashes
and checked CLI restoration complete the handoff. These checks install
actual function bodies; they do not select a numerical point field,
prove uniform original collar accuracy or close five moments.

## Next executable tasks

The true even radial jets, original flat cutoff and all five signed inner
endpoint functions are now installed in
[CURRENT_ORIGINAL_N1_INNER_CUTOFF_FUNCTIONS_2026_10_13.md](CURRENT_ORIGINAL_N1_INNER_CUTOFF_FUNCTIONS_2026_10_13.md).
Original numerical n1 evaluation, full exterior pressure composition and
five-moment repair remain open.

- [x] Bind the true core logarithmic/Psi collar ODEs and complete histories.
- [x] Recover ordinary C3 velocity/mean/base pressure and C2 Q/radial rows.
- [x] Preserve exact nonzero hb/sc and original physical measures/support.
- [x] Export signed m13/m15 C1 and compare independent physical forcing.
- [x] `N1-TRUE-INNER-VOLTERRA-FUNCTIONS`: bind the actual regular six-row
  solution from `current_original_n1_regular_series_function.py` and
  `current_original_n1_coupled_first_step.py`. Preserve the Green
  integral, complete Volterra series, true leading forcing and feedback
  matrices. Existing k=0,1 words plus symmetric tail boxes are enclosures,
  not selected F1/U1/K1/P1 functions. Recover its own M1 and Q1 with the
  correct physical rho/Lambda factors before constructing the cutoff.
- [x] `N1-INNER-CUTOFF-FUNCTIONS`: trace
  `current_original_n1_radial_cutoff_extension.py` and
  `current_original_n1_five_defect_values.py` to construct the actual
  cutoff endpoint sectors at first-fraction q=1 (phase sc/4). The flat
  cutoff is zero after q=3/4, but own M/Q and cumulative pressure are
  retained. Use the same true core/collar histories; do not select
  values from `original_signed_Z1_sectors` error enclosures.
- [x] Bind the three signed velocity-support functions m11/m12/m14
  and all five Rin cumulative functions from the actual cutoff source.
  Exterior m13/m15 pressure contributions still require composition.
- [ ] Transfer existing same-source pressure/C5 receipts as bounds to
  the complete function bodies, with all nonzero source errors.
- [ ] Compose full five defects across inner cutoff, this collar,
  active/compact/quiet, pulse and exact exterior once in the common graph.
- [ ] Select actual alpha/beta functions with the reserved linear
  operator; recover repaired own M/V/P and prove terminal five zeros
  as axial functions.
- [ ] Complete matching/finite-energy/heat and admissible stress/flat
  remainder before genuine temporal coefficient recursion.
- [ ] Implement oscillatory quadratic cancellation and independently
  validate corrected Cartesian NS with fixed/restricted forcing.

The full original objective remains in [PROJECT_GOAL.md](PROJECT_GOAL.md).
