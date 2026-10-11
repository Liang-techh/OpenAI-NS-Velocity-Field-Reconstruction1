# Original real core functions and genuine n=1 forcing — 2026-10-10

The canonical original leading-core source is now directly callable over the real rectangle rho in [0,4.1], Z in [-1,1]. It delivers ordinary mixed derivative enclosures of F_(0), Uz_(0), regular V_(0)/R and both pressure pieces, with genuine hierarchy-n=1 known forcing. The larger endpoint corresponds to R_in=4.1/Lambda > R_a=4/Lambda. These are enclosures of the admitted nonlinear analytic source, with no selection of a point solution or parameter representative. Full goal **ACTIVE / INCOMPLETE**: positive-order solutions, moment repair, a common complex-domain Assumption 14.1 certificate, regional cone, time flatness and oscillatory correction remain open.

Source commit: [10b721fe](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/10b721fe9cedd55acde936129eb58c7a8aeb33af). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_hierarchy_source` with `.py`, `.json.gz`, `_check.py`, `_check.json`.

## Source and parameter ownership

The adapter consumes the live accepted `CurrentOriginalRpNativeConstants` to identify the current outer source, and `OriginalWholeZBridgeSource` to supply the canonical actual inner source. The former is an outer constant/function frame; it does not define the inner grid. The new `OriginalCoreContext` shares the canonical interval context and independent pressure datum, and has its own cache. It re-evaluates the unchanged `CompliantCorePhysicalField.axis_inputs`, `normalized_jets` and `profiles` callbacks using that canonical tuple.

- logPstar=exp(40)+11; logdelta=-4logPstar-30; delta>0 is retained.
- logLambda=4logPstar+1000; Lambda=exp(logLambda); epsilon_core=exp(-logLambda); Lambda*epsilon_core=1. Live bounds and the original AST/accepted exact identity are both checked.
- The pressure-datum epsilon is delta/1000. It is distinct from epsilon_core and is never substituted for it.
- Selected logCstar, j, sigma=j/500, nonlinear correction bounds and anchored G(Z) come from the same admitted analytic core family.
- The current Rp and inner P0 have the accepted exact analytic function, fourteen continuous atom/flatten sources, normalized-jets callable and hydration AST identity. The atom enclosures remain enclosures; no mass endpoint is selected.

The accepted `CurrentCoreFirstInterface` already closes the pressure 4C and original first-interface identities. Its newer receipt supersedes the older common-fixed-point receipt's remaining-dependency note. Do not repeat that completed proof as a new blocker. No new whole-chart physical overlay or global field is claimed by this direct real core source API.

## Actual leading derivative interface

Each source packet supplies all ordinary mixed derivatives with radial order i plus Z order k <= 4. The grid coordinate is rho=Lambda R. Their scale/coefficient representation is:

| Field or piece | Factored source | Ordinary R derivative factor |
| --- | --- | --- |
| F_(0) | exp(logF0) * normalized F derivative | Lambda^i |
| Uz_(0) | actual nonlinear axial derivative enclosure | Lambda^i |
| V_(0)/R | actual regular Q derivative enclosure | Lambda^i |
| P_(0) axis datum | exp(2logPstar) * normalized PD derivative | Lambda^i; radial derivatives vanish for i>0 |
| P_(0) centrifugal increment | epsilon_core exp(2logF0) * PI derivative | Lambda^i |

Thus the centrifugal log scale is 2logF0+(i-1)logLambda. The axis datum remains independent. `logF0=-selected_logCstar-Lambda*G(Z)` is evaluated through the actual anchored G enclosure; its extreme exponential is not materialized.

Regular Q is [2Z Uz0-(1-delta)Z(Mz0/R)-(1-Z^2)partial_Z(Mz0/R)]/(1-delta Z^2). Mz0/R alone is not Q. No interval division by R is used at the axis. Internally factorial-normalized Taylor data are converted by the original callback to ordinary derivative grids before this API exports them.

## Genuine hierarchy-n=1 known forcing

With L=1-delta Z^2 and d=1-Z^2, the paper operators are

```text
T_a G = [-a G/2 + (1-delta)Z G_Z/2 + R G_R]/L
Z_a G = [a Z G + d G_Z - 2 Z R G_R]/L
Z_a^[2] G = Z_(a-1+delta)(Z_a G)
```

The adapter uses e1=2delta, b0=-2-delta, c0=-1-delta and p0=-2-2delta, then emits

```text
N1theta = -Z_b0^[2] F_(0)
N1z     = -Z_c0^[2] Uz_(0)
N1p     = -Omega_(0)/(2R)
```

For V_(0)=R Q, Omega_(0)/R is evaluated regularly as

```text
[Q+(1-delta)Z Q_Z/2+R Q_R]/L
+ Q*(Q/2+R Q_R)
+ Uz0*(-2Z Q+d Q_Z-2Z R Q_R)/L
- 4 Q_R - 2R Q_RR
```

The R derivatives use Lambda and Lambda^2 conversions from the actual rho jets. At n=0 there is no negative-index axial-viscosity term. The pressure source is only the known part: the eventual equation is P1_R=2F0 F1+N1p. The generated radial Taylor row indexed 1 in the existing leading-core rebuild is not used as hierarchy order 1.

## Focused acceptance and scope

Five actual deliveries cover the entire real rectangle, the entire axis, axis Z=.371, inner rho=2/Z=.371, and rho=4.1/Z=.371 beyond R_a. A separate 800-digit canonical source context replays 225 leading derivative rows, 150 pressure pieces, 15 genuine N1 forcing rows and 10 axis rows. Independently collected operators and 9 exact polynomial/axis fixtures check the weighted composition and regular Omega0/R. Copied packets, changed original delta, changed forcing and out-of-domain requests are rejected. The source/data/receipt inputs are hash-bound.

This reuses the inherited analytic fixed-point admission and its model-tail/nonlinear correction bounds. It does not numerically solve the nonlinear fixed point anew. The common holomorphic extension certificate is not promoted merely because the real callback or an inherited Cauchy majorant exists. None of F1, Uz1, P1 or their repaired moments has been solved by this layer.

One existing read-only worker is reused: **GPT-5.6 Luna / max**. Its scoped review found no material defect. No new worker or Astra child was spawned.

## Executable use

In the accepted live source chain, `constants` is the `CurrentOriginalRpNativeConstants` owner. Append:

```python
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source import OriginalWholeZBridgeSource
from lei_ren_part1_paper_compliant_current_original_core_hierarchy_source import CurrentOriginalCoreHierarchySource

core_source = OriginalWholeZBridgeSource(500)
hierarchy_source = CurrentOriginalCoreHierarchySource(constants, core_source)
packet = hierarchy_source.evaluate('2', '.371')
report = hierarchy_source.report(packet)
assert set(report['genuine_order_one_known_forcing']) == {'N1theta', 'N1z', 'N1p'}
assert report['current_original_core_real_functions_and_genuine_n1_known_forcing_installed']
assert not report['coefficient_F1_Uz1_P1_solved']
```

For real-domain coverage use `evaluate(['0','4.1'], ['-1','1'])`. Report JSON cannot hydrate live typed source owners. Keep coordinates as exact source literals when replaying at another precision; a rounded endpoint enclosure must not be passed as a new uncertain coordinate outside the intended exact domain.

## Completed and next tasks

- [x] **CURRENT-ORIGINAL-CORE-REAL-FUNCTION-SOURCE** Bind canonical original parameters/P0, re-evaluate leading mixed derivative enclosures on a common real interval beyond R_a, expose Q regularly at R=0 and preserve pressure pieces.
- [x] **CURRENT-RP-N1-FORCING-PACKET** Emit the actual original three known forcing enclosures with ordered weighted axial operators, original delta and regular pressure source. Source-only completion; no positive-order solution flag.
- [ ] **CURRENT-RP-CORE-COMMON-HOLOMORPHIC-DOMAIN** Bind the existing analytic fixed-point/tube/Cauchy majorants to every required original leading radial derivative over this same interval. State the common complex neighborhood and original parameter dependence. Implement a complex source evaluator if the planned solve requires it. A real rectangle enclosure alone does not close Assumption 14.1.
- [ ] **CURRENT-RP-N1-REGULAR-AXIS-INITIALIZATION** Use positive-order F1(0,Z)=Uz1(0,Z)=P1(0,Z)=0 to derive and expose the original first radial slope source data, and n=1 regular V1 recovery. Preserve the actual F0 factor/log scales and source functions over Z; do not insert arbitrary axis jets or confuse them with the leading radial Taylor rows.
- [ ] **CURRENT-RP-N1-COUPLED-SOLVE** Use this live hierarchy source to solve the paper's regular system for F1, Uz1, K1, P1 and xi derivatives, xi=sqrt(R), K1=Mz1/R-Uz1, on the same common interval. Preserve P1_R=2F0 F1+N1p and all n-dependent weights. Bound discretization/truncation and validate all three equations plus recovered V1. The known forcing alone is not a solution.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Derive actual order-one defects as Z functions, bind bump controls/conditioned inverse, cut streamfunction or vector potential before curl, verify all five repaired moments, joins and pressure compatibility. Only after repair build n=2 forcing/recovery on the same R_in.

Continue the remaining regional cone, all-chart, finite-width/norms, measured-core/winding, finite-order and flat remainder, and oscillatory tasks in [the multitime task list](CURRENT_ORIGINAL_RP_MULTITIME_REMAINDER_2026_10_10.md). The full objective stays active.
