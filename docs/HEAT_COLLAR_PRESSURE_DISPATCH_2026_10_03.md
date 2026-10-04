# Full collar absolute pressure and heat companion dispatch - 2026-10-03

The full original heat collar now has a source-bound absolute pressure
companion. Its domain is Z in[-1,1], t=log(R/Rtail) in[0,3]. Here t is the
radial chart offset, not physical time. The companion supplies axial Taylor
order5 and ordinary logR/Z mixed derivatives through total order4.

`CompliantCollarPressureC4().collar(Z,t)` uses the original remaining
integral, including the entire unbounded Gamma exterior:

```
K(v,Z) = (1-sigma(v))*(1-epsilon)
         + sigma(v)*H(2*(1-Z^2)*S*exp(-v))*(1-epsilon*phi(v))
C = Ev2*theta_base^2
Pi(t,Z) = integral_t^infinity exp(-(1+delta)*v)*K(v,Z)^2/2 dv
P/Pstar^2 = -C*Pi(t,Z)
```

H is the full positive Gamma expectation. The original flat sigma and phi
are retained. S is the exact positive inverse Rtail; interval caps only
enclose it. No exterior cutoff, finite S expansion, new pressure constant,
velocity adjustment or repair coefficient fit defines this pressure.

The bridge inherits the accepted actual Ptail/P3 history, common K-squared
density, canonical infinite Gamma tail and original absolute closure.
Splitting that same integral at arbitrary t gives
`P_forward(t)+C*Pi(t)=original_infinity_constant=0`. Production assignments,
remaining-integral/output routes and the companion's derivative helper are
bound explicitly. Interval overlap is a diagnostic, not this proof.

The zeroth pressure row uses reference Rtail units directly. Positive-order
rows follow from the exact FTC/Leibniz expression
`d_t P=C*exp(-(1+delta)*t)*K(t,Z)^2/2`. They use the full collar bracket,
including sigma and phi before t=3. This avoids introducing cancelling
interval exponential factors into the zeroth row. The original flat K jets
and common complete future integral establish the waiting/collar and
collar/exterior pressure joins through mixed4.

The focused producer/checker passes. It covers120 actual finite pressure
mixed bounds,120 retained-forward overlap diagnostics and30 endpoint
overlap diagnostics. An independent moderate-parameter full Gamma/flat
fixture checks3 complete future integrals and21 mixed derivatives. This
fixture exercises the integral evaluator; it is not project NS residual or
physical blowup validation data.

The main source routes now select:

| Chart | Provider | Domain | Acceptance receipt |
| --- | --- | --- | --- |
| heat_collar | CompliantCollarPressureC4.collar | t in[0,3] | collar_pressure_C4_check |
| heat_exterior | CompliantHeatStressC4.exterior | t>=3 | heat_stress_C4_check |

The physical source assembler unwraps the companions' original `heat`
object for radius offsets, preserving the original flatten/angular/steep/
waiting history. Pressure and velocity rows remain in their original
normalizations. The reconstruction workflow has an ordered `heatcompanions`
stage before source dispatch, so receipt rebuilding has an explicit route.
The affected dispatcher and physical source assembly receipts are rebuilt.
All33 dispatched charts pass, and the assembly checker passes its140
independent Cartesian derivative,4 fixed-x time-derivative and120
micro/macro factor checks. These checks validate source-bound operators and
enclosures; they do not certify the physical stress equation or full NS
residual.
The unchanged exterior stress receipt was produced before adoption and still
records `main_dispatcher_and_collar_adoption_complete=False`. Its theorem
acceptance remains current; the fresh dispatcher receipt records the later
adoption. That older status field is not a stress acceptance gate.

## Scope and next work

Absolute pressure and original similarity stress in the Gamma exterior are
available through the main source API. This does not yet transfer the
exterior stress identity to the complete physical momentum equation.
All-region admissible stress/cone, independently bounded flat remainder,
physical energy, finite-width bridge/second switch, nonlinear global point
field, genuine n-dependent coefficient recursion and oscillatory correction
remain unfinished. Global completion flags remain false.

Next: apply the actual physical coordinate map and original stress
prefactors to the admitted exterior moment identity, with an independent
check of the mapped heat equations. Then continue remaining leading-field
stress regions, cone margins and flat remainder. Keep the parallel
finite-width annular and coefficient-recursion admission blockers visible.

## Reproduction

From the repository root, run the new collar producer and checker, then the
affected dispatcher and physical source assembly receipts:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_collar_pressure_C4.py
python experiments/root_st073/lei_ren_part1_paper_compliant_collar_pressure_C4_check.py
python experiments/root_st073/lei_ren_part1_paper_compliant_source_dispatcher.py
python experiments/root_st073/lei_ren_part1_paper_compliant_global_physical_assembly.py
python experiments/root_st073/lei_ren_part1_paper_compliant_global_physical_assembly_check.py
```

For rebuilding all admitted heat companions after relevant source changes,
use `lei_ren_part1_paper_compliant_reconstruction.py --stage heatcompanions`.
Do not rerun unchanged accepted checks solely to increase check counts.
