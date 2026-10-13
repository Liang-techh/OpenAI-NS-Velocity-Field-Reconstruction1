# Exact original pulse kernels and C1 selection, 2026-10-13

The original finite-mu pulse kernels are now real integral expressions in the
same graph as the checked pre-Rp pressure functions and n1 repair operator.
This removes the need to select fixed kernel interval centers as coefficients.
The actual native inlet u/m1/m2/energy C1 handles are also imported, retaining their original own histories and common normalization. These are prerequisites for the missing full-pulse pressure function adapter.
The full original goal remains **ACTIVE / INCOMPLETE**.

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local receipt gate is
`current_original_n1_exact_pulse_kernel_functions_installed`.

## Source definitions retained

- The old pulse end supports remain centered at -3 and -1 with ell=.15.
  These are distinct from the five fresh positive-order n1 repair supports.
- The exact original normalization is the raw beta integral over [-1,1].
- A and D are their genuine finite-mu integral functions. D uses
  `(-s)*exprel(-mu*s)` instead of subtracting near-equal row enclosures.
- The negative determinant is
  `-2 exp(-2+5mu) J(.5-mu) J(.5-2mu) exprel(2mu)`.
  Its positive absolute value supplies the correct signed inverse action.
- `gp(xi)=integral_0^xi sigma(50v)dv * sigma(11-xi)`, including both startup
  and exit transitions. Its energy is the original integral
  `Kpulse=integral_0^11 exp(-2xi) gp(xi)^2 dxi`.
- Both scaled full main rows retain the actual finite-mu exponential body
  and explicit original common prefactor, not a saddle interval center.
- Gram, Kj, incoming exponential factors, the original log_end_scale and
  `nu_j=mu*Kj*exp(2log_end_scale)>0` are retained. No end energy is zeroed.

The C1 selector accepts two genuine incoming moment functions, their energy,
and the complete future energy function. It reproduces the original divided
inverse, completed-square positive root and signed end coefficients with
true product/quotient/square-root first-Z derivatives. The supplied original
future must still prove the discriminant and chosen positive branch conditions;
this algebraic API alone does not certify them for arbitrary inputs.

## Validation and use

Direct mathematical fixtures at nonzero mu=.07 compare A/D, Gram and Kj to
independent physical bump integrals and prove the determinant identity for
every positive mu. A separate nonconstant positive-amplitude fixture compares
the selected amplitude and two signed coefficients, including first-Z rows,
against an independently solved matrix/energy system. Numeric caps are rejected.
Fixtures are labelled mathematical tests; they do not substitute mu=.07 or a
test future into the actual original construction.

Run from the repository root:

```powershell
python scripts/run_original_n1_pulse_kernel_functions.py
```

Import `CurrentOriginalN1PulseKernelFunctions()` for checked handles. Its `.g`
and `.preRp.operator.graph` are identical. `.functions` contains the original
kernel references. The `.incoming` dictionary contains the actual native u/m1/m2/energy C1 functions. Use `.original_coefficient_functions(future)` to select from these original incoming functions. Call `.coefficient_functions(incoming,energy,future)` only
with genuine C1 functions in this graph and the admitted original positive
branch. The constructor rebuilds this small function adapter from saved
checked ancestors; it does not run ancestor producers or pressure grids.

## Current work queue

- [x] `N1-PULSE-KERNEL-FUNCTIONS`: bind original A/D, determinant, Gram/Kj,
  gp, main rows, Kpulse, incoming factors and positive nonzero nu in the common graph.
- [x] Preserve the original finite mu and all normalization/radius factors.
- [x] Install the original conditional positive-root C1 selector and focused
  direct-integral/first-Z checks.
- [x] `N1-PULSE-INCOMING-FUNCTIONS`: import all actual native-frame inlet
  m1/m2/energy rows into the common graph. Reuse the existing exact native-frame
  function definitions, not the copied selected runtime's interval constants.
- [ ] `N1-PULSE-FUTURE-FUNCTION`: bind the actual complete postpulse swirl
  energy through flattening and heat exterior as a C1 function, including
  its real body, source parameters, original P0 and an enclosure oracle.
- [ ] Prove the admitted original positive discriminant and (.9,1.2) branch
  from the same true inputs, then persist actual ap/C1/C2 function handles.
- [ ] Build main/gap/end U and their own mean from the same ap/end coefficients.
  Recover q, retain the actual lower mean and the empty-future terminal identity.
- [ ] Integrate both full pulse pressure functions with actual radial measures;
  retain nonlinear energy and viscous endpoint terms and their Z derivatives.
- [ ] Join pulse pressure with pre-Rp, lower collar and inner-cut moments once,
  select full positive-order alpha/beta and recover repaired own M/V/P.
- [ ] Verify the complete five terminal functional identities before n>=2.

The parallel read-only source review also confirmed a core dependency for
`Rin -> r_minus`: the Bessel Phi and linear Psi are models with nonzero
fixed-point corrections. They cannot replace the actual core functions.
Use the true admitted core fixed-point law and its numerical enclosure oracle;
then build the collar ratio/FTC histories from the exact equations below.

```
rho=4 exp(hb*s), R=rho/Lambda, s in [sc/4,sc/2]
dR=R*hb ds, R dR=R^2*hb ds
chi=1+(hb-1)*sigma(s)
log(Phi_actual/Phi_core)'=(chi-1)*Phi_core'/Phi_core
Psi_actual'=(Phi_actual/Phi_core)*chi*Psi_core'
rho*meanPsi_actual=4*meanPsi_core(4)+integral_4^rho Psi_actual drho
U=Pstar*(4Z+j+epsilon*Psi), m=Pstar*(4Z+j+epsilon*meanPsi)
q=[2ZU-(1-delta)Zm-(1-Z^2)m_Z]/(1-delta*Z^2)
```

Here prime denotes the appropriate rho derivative in the ratio/Psi equations.
The radius map and its endpoints are independent of Z. Keep all original
corrections and the common P0; no constant-q collar primitive is admitted.

Full-pulse C1 pressure, full five-defect binding, actual repaired n1 fields,
terminal closure, global heat/stress, temporal recursion, flat remainder,
oscillatory cancellation and corrected Cartesian NS validation remain open.
