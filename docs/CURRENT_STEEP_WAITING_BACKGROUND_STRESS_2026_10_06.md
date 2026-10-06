# Current O7 power/exit/waiting tensors and three joins (2026-10-06)

The actual current steep-power, steep-exit and waiting regions now have full-moment stress mixed3, symmetric completed physical tensors, physical divergence/remainder mixed2 and the regional identity `residual=-div(T_B)+E_B`. Their entry-power, power-exit and exit-waiting completed-tensor joins use source-function identities and common physical contribution bounds. Together with the admitted angular and entry tensors, this gives five actual tensor regions and four adjacent completed-tensor joins. The velocity/absolute-pressure trace inventory remains separate: 14 adjacent and 8 internal.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_steep_waiting_background_stress.py`, matching `_check.py`, producer `.json` and scoped `_check.json`. Focused controller stage: `currentsteepstress`. Reuse the checked [current entry tensor owner](CURRENT_STEEP_ENTRY_BACKGROUND_STRESS_2026_10_06.md); it retains the same checked angular/physical/selected/future/pressure graph. The full reconstruction goal remains active.

## Callable construction and domain

```python
from lei_ren_part1_paper_compliant_current_steep_waiting_background_stress import CurrentSteepWaitingBackgroundStress
stress = CurrentSteepWaitingBackgroundStress(entry=checked_current_entry_tensor)
view = stress.chart("steep_power", Z=".537", x=".417", log_tau="-2.6", theta=".41", viscosity=".8")
join = stress.interface("steep_entry_power", Z=("-.8", ".8"), log_tau=("-5", "-2"), viscosity=".2")
```

Charts are `steep_power`, `steep_exit`, and `waiting`. Their source coordinate x is respectively power phase, exit log-radius offset t, or waiting phase, each in [0,1]. Ordinary radial offsets are y=Ts*phase and u=wait*phase. All returned radial derivative rows are ordinary logR derivatives; they are not phase derivatives. Physical time is supplied separately by log_tau. The finite physical domain is R>0, |Z|<1, tau>0 and constant viscosity nu>0. Source Z=+/-1 are infinity limits, not the axis; theta=None covers all angles. Compact positive-time estimates do not extend uniformly to tau=0.

Outputs are directed source enclosures with exact positive scale logs. They do not select interval centers or turn amplitude/inverse-radius caps into exact physical field values.

## Complete moments in the same KR units

Let a=delta/2, k=1-a, r=1-mu and p=1+delta. Use the exact angular-entry normalization KR already admitted upstream. Define

`LS = a-mu-r/2`, `LQ = LS-k*Ts`, `LT = LQ-k/2`.

The three exact normalized shapes are

- Power: `K/KR = exp(LS-k*y)`, y=Ts*phase.
- Exit: `K/KR = exp(LQ-k*t+k*J(t))`, with the original sigmoid primitive J.
- Waiting: `K/KR = exp(LT)`, independent of logR and Z.

Native theta and the new normalized K expressions are executed as symbolic source ASTs on the same arbitrary parameters. The endpoint shapes agree, their ordinary log rates are -k, k*(sigma-1), and 0, and the flat sigmoid endpoint jets match. The exact source logKR satisfies `LT+logKR=logone`, so the actual waiting shape is 1-epsilon; the normalized shape is not confused with the local heat-reference shape.

Full angular and energy moments are `A/KR=(K/KR)*native_X` and `E/KR^2=2*(K/KR)^2*native_half_energy`. Native X, complete Gamma future, both epsilon atoms and the nonzero cumulative histories remain included. The original power half-energy floor 1/4 is retained. Full ordinary rows obey A'=K-kA, E'=delta*E-K^2, P'=p*P-K^2/2. Ninety identities replay the original power/exit signed-defect recurrence ASTs and remove their unit baselines before KR scaling. The checked entry's arbitrary log-rate A=KX/correlated-stress theorem applies to all three Z-independent shapes. No legacy terminal constructor or replacement history is introduced.

## Stable original absolute pressure

P here is the positive remaining normalized moment; the absolute-pressure output restores the minus sign. Its derivative is pP-K^2/2. This follows by differentiating exp(p*v) times the remaining source integral. The current analytic Cp=0 identity fixes the original PR exactly once; no new pressure gauge or constant is selected.

Let Ih be the complete original heat remaining-pressure function at collar offset zero. Define

`Jwaiting = Ip(wait)/2 + exp(-p*wait-2*logone)*Ih`,

`JafterPower = Iout + exp(-3+k)*Jwaiting`.

Here Ip(L)=integral_0^L exp(-p*s)ds, and Iout is the original full exit pressure integral including its 1/2 factor. The stable full normalized moments are

- Power: `P/KR^2 = (K/KR)^2*[I3(Ts-y)/2 + exp(-3*(Ts-y))*JafterPower]`.
- Exit: `P/KR^2 = exp(2*LQ+p*t)*[Iout_suffix(t)+exp(-3+k)*Jwaiting]`.
- Waiting: `P/KR^2 = exp(2*LT)*[Ip(wait-u)/2+exp(-p*(wait-u)-2*logone)*Ih]`.

Iout_suffix(t)=integral_t^1 exp(-3v+2kJ(v))/2 dv is enclosed directly with positive correlated cell lengths. Its integrand is replayed against the original transition AST. The reduction cancels source exponents before interval enclosure: it never divides a tiny amplitude cap or multiplies independently enclosed huge amplification and tiny attenuation. Prefix/suffix additivity identifies the same original forward absolute-pressure function, preserving P0/Pin and the complete future.

## Physical tensors and function joins

The checked actual general-K physical tensor assembler is replayed unchanged. Full current moment rows feed the original correlated stress operator; physical lambda/viscosity factors, moving cylindrical basis, source radius and `Ttheta_theta=r*partial_z(Trz)` are retained. The inverse-radius cap only encloses the exact 1/R source.

Each of the three joins consumes 60 current primitive mixed4 identities, original arbitrary-terminal stress/pressure AST join theorems and the exact ordinary-coordinate/radius proof. Current normalization and the common general-K physical operators transfer these identities to stress3, divergence2, diagonal2 and remainder2. Differently rounded valid endpoint envelopes are merged after the source-function equality; interval overlap or byte-equal bounds do not supply the function proof.

Waiting has a source-exact zero swirl axial-viscosity remainder because K is independent of logR and Z and the physical swirl is a pure radial power. Its angular/energy/pressure histories and full stress remain nonzero. This regional zero does not establish an independent global temporal-flat remainder.

## Scoped evidence and next tasks

Twelve views cover each entire chart, both endpoints and a fresh angle/time/viscosity point: 780 factored physical tensor/divergence/remainder rows, 694 with nonzero enclosures. Each new tensor interface checks 65 common rows; a fresh compact sector repeats each interface with different Z/time/viscosity. The original recurrence normalization has 90 identities and the three current primitive function joins have 180 identities. Invalid domains and a foreign current owner are rejected. The focused controller retains five tensor regions/four tensor joins and the separate 14/8 velocity-pressure atlas. Global cone/lift/NS/energy/points and temporal-recursion flags stay false.

Read-only review metadata: GPT-5.6 Luna / max. Scoped review passes with no material blocker: KR units, remaining-pressure orientation, restored absolute-pressure sign, nonzero histories and source tensor joins are consistent. The reviewer corrected an earlier sign note: the positive reduced remaining moment obeys P'=pP-K^2/2. Explicitly recording the heat endpoint exponent identity is optional provenance hardening for the next waiting-collar task, not a blocker for these three regions. The actual focused controller runs one shared producer/checker graph; compilation passes.

Raw working-file and Git-index SHA256 checks pass for all 656 source/receipt dependencies. Git whitespace checks pass. Existing unrelated experimental edits are preserved.

- [x] F57C5a-postpulse-power: actual full current steep-power tensor, native half-energy floor and absolute pressure.
- [x] F57C5b-entry-power-join: source-function completed-tensor equality and common physical bounds.
- [x] F57C5a-postpulse-exit-waiting: original sigma exit and full-history waiting tensors.
- [x] F57C5b-power-exit-join/F57C5b-exit-waiting-join: actual completed-tensor joins and common bounds.
- [ ] **Next: F57C5a-postpulse-collar/F57C5b-waiting-collar-join.** Reuse this checked owner. Bind current heat shape at offset zero to actual 1-epsilon, convert its full A/K/E/P and radial rows into common KR units analytically, retain the exact Ih/Ptail/full-Gamma source and current terminal identities, then assemble the original general-K tensor and prove the waiting-collar stress3/divergence2/diagonal2/remainder2 join. Do not divide amplitude caps or treat local heat K as K/KR. Keep original collar t in [0,3].
- [ ] F57C5a-postpulse-exterior/F57C5b-collar-exterior-join: convert the checked current zero-stress exterior theorem into the actual completed physical tensor and source-exact heat attachment. Distinguish zero similarity stress from the axial-viscosity remainder; preserve the canonical infinite Gamma function and physical operators.
- [ ] F57C5a-postpulse-flatten-outerpower: actual current flatten and outer-power tensors with source Jacobians, complete histories, remaining original pressure and the power-angular tensor join.
- [ ] F57C5a-pulse: actual entrance/main/exit/gap/end tensor from all five histories, selected C1/C2, meridional cross terms, complete future/2 and original Pin/P0. Local support differences are not actual tensor values.
- [ ] F57C5b/F57C6a-global: every actual tensor join and chart decomposition, including the eight internal support edges; keep actual tensor inventory distinct from 22 velocity-pressure traces.
- [ ] F57C4e-remaining/F57C6b/F57D/E/F: retained/core/incoming and axis attachment, independent global temporal remainder/cone/lift, prescribed-domain kinetic energy and resolved source-defined u,v,w,p.
- [ ] F58/F59/F60/F61: genuine n-dependent temporal equations/per-order repairs/summation, both oscillatory families and mean corrections, corrected NS and measured contraction/slenderness/winding.

For each bounded construction, save source evidence and the scoped receipt, mark exactly the implemented task complete and commit/push. Spatial moment ODEs, chart inventories and constant waiting K do not close temporal scale recursion.
