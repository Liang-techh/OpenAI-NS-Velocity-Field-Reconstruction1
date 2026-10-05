# Current external Rh source join — 2026-10-05

The current nonlinear five-moment patch now joins the actual outer reference owner at Rh=e*Rm=e^-5*Rref. A separate CurrentRhReferenceDispatcher exposes the previous eight current charts plus Rh_reference. The accepted eight-chart dispatcher, legacy ROUTES and historical pre-pulse artifacts remain unchanged.

Run the focused stage on codex/st073-transition-next:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualrhjoin
```

For a checked source call:

```python
from lei_ren_part1_paper_compliant_actual_Rh_source_join import CurrentRhReferenceDispatcher
field = CurrentRhReferenceDispatcher()
packet = field.evaluate("Rh_reference", Z=".3", coordinate=-5)
```

Rh_reference has offset=log(R/Rref) in[-5,0]; its derivatives are ordinary logR/Z derivatives. The explicit-chart API returns directed source-field enclosures and original formal scales. It does not select pressure or velocity point representatives.

## Analytic pressure function and projection

The proof binds both constructor/import paths to exactly CompliantPressureDatum('40',160), with the repair alias using precision=160. The outer reference reads initial.repair.datum through SharedOuterBuffer; the current patch retains its core datum through the actual bridge/switch/reshape/restoration chain. Both use the compliant epsilon=.001*delta source. The legacy .01 constructor is excluded.

Runtime callable/constructor checks, source AST bindings and common source definitions establish one analytic P0 function. The fourteen stage sources, beta routing, early mass refinement, m2/m0, rho and nonzero flatten enclosure remain common. Equality of source hashes or interval boxes alone is not the defining-function proof.

In P/Pstar² units the common defining function is:

```text
P0(Z)/Pstar² = -(M2*(1+Z²)^-2 + M0 + F_flat(Z))
q_n(Z) = partial_Z^n (1+Z²)^-2 / n!
P0_n = -(M2*q_n + [n=0]*M0 + partial_Z^n F_flat/n!)
```

F_flat is one common implicit analytic source function. Its positive enclosure and derivative bounds remain enclosures; the six coefficients are not independently selected values. The same inherited normalized_jets callable uses an order-independent prefix recurrence. The exact recurrence is checked symbolically through order6. Order5 returns six coefficients0..5; order6 returns seven and IntervalTaylor.truncate(5) selects the same first six. Three common-center diagnostics actually execute that truncation, checking18 enclosure coefficients after the source-function proof.

The current actual_R100_trace_bindings also consumes the original inner_switch_profiles source/P0 projection AST, so the proof applies to the actual retained pressure call chain.

## Nine physical source functions

The new proof consumes the accepted Rm source-coordinate relation y=log(R/Rm)=offset+6 and x=exp(y). Thus offset=-5 is y=1,x=e. All original bump supports end at71/40<e. On the outgoing neighborhood, the same current unique implicit map supplies full-weight closure of all five defects as axial functions.

With am=exp(-.6)/(1+Z²), the canonical terminal formulas agree with pre_pulse.reference:

```text
u/Pstar = am*x^.1 = exp(offset/10)/(1+Z²)
V = 4Z
h = 5u/8, k = 4Z*h
e = 16Z²/Pstar² - 5u²/12
p = 5u²/2
```

Both physical wrappers are AST-bound, including P0+Mp, the Pstar factors, all five primitive units and physical radial prefactors. Nine exact physical function identities imply135 mixed derivative rows through total order4. The functions share the canonical analytic extension on the corresponding sides of Rh; no finite chart outside its original domain is substituted.

The right reference grid is converted into the patch Rm units only after physical differentiation. Fixed Rh factors are sqrt(e) for Ur, e for Mz/Mztheta and exp(3/2) for angular primitives. They are frozen at the basepoint. The radial operator is D_y=x D_x, with the original Stirling combinations for higher powers. Interval overlap is not the join proof.

## Acceptance ownership and next work

Focused acceptance passed:341 current dependency hashes,6 exact pressure source coefficients,6 q-recurrence identities,18 prefix enclosure diagnostic coefficients,9 physical function identities/135 implied mixed rows and135 whole-domain reference rows. A fresh checked lazy call loaded the new acceptance and exposed9 chart owners.

The producer records source_join_proved=True and external_neighbor_join_certified=False. Only a runtime that consumes actual_Rh_source_join_check.json exposes certification=True. The native eight-chart Rh flag remains false, preserving its original scope.

Next extend current ownership through the remaining original pre-pulse charts and their Rp input/history transfer; then inject the current providers into the existing physical/Cartesian coordinate operators. Shared full leading inputs/nonlinear remainders, production points, completed global tensor/flatness/physical-volume/required-domain energy, genuine n-dependent coefficient recursion, oscillatory corrections and corrected dynamics remain open. Their flags remain false.
