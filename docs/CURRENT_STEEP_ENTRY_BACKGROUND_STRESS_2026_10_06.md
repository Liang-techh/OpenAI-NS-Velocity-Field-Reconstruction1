# Current actual steep-entry tensor and angular-entry join (2026-10-06)

The original current steep-entry region now has full A/E/P/K stress mixed3, its symmetric completed physical tensor and the regional decomposition `residual=-div(T_B)+E_B`, including physical divergence/remainder mixed2. The actual current angular and entry tensors also have their first function-level completed-tensor join and 65 common physical contribution bounds. This completes F57C5a-postpulse-entry, its regional F57C5b/F57C6a tasks, and F57C5b-angular-entry-join. It does not close the global versions of those tasks.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_steep_entry_background_stress.py`, matching `_check.py`, producer `.json` and scoped `_check.json`. Focused controller stage: `currententrystress`. It consumes the checked [actual angular background tensor](CURRENT_ANGULAR_BACKGROUND_STRESS_2026_10_06.md), whose owner contains the checked [current22 velocity/pressure atlas](CURRENT_INTERFACE_ATLAS_2026_10_06.md). Inventory remains 14 affected adjacent / 8 internal velocity-pressure traces; the actual completed-tensor adjacent-join count is now one.

## Callable result and source domain

```python
from lei_ren_part1_paper_compliant_current_steep_entry_background_stress import CurrentSteepEntryBackgroundStress
stress = CurrentSteepEntryBackgroundStress(angular=checked_current_angular_tensor)
point = stress.entry(Z=".537", t=".417", log_tau="-2.6", theta=".41", viscosity=".8")
join = stress.interface(Z=("-.8", ".8"), log_tau=("-5", "-2"), viscosity=".2")
```

The source entry coordinate t belongs to [0,1] and is ordinary log-radius offset, not physical time. The heat-reference offset is q=-wait-Ts-2+t; physical time enters separately through log_tau. The source domain Z=[-1,1] has finite physical points R>0, |Z|<1 and tau>0. Z=+/-1 are infinity/source limits, not the axis. `theta=None` covers all angles; viscosity must be a finite constant nu>0. Compact positive-time bounds do not extend uniformly to tau=0.

Outputs contain the actual normalized full moments in the same exact KR/KR^2 units as the angular tensor, full stress mixed3, completed physical tensor, divergence/remainder mixed2, Cartesian residual decomposition, factored absolute-pressure mixed4 and common angular-entry tensor bounds. Every numerical row is a directed enclosure of a source function with exact logarithmic scale factors. No interval center or cap defines a resolved physical point value.

## Full current source construction

Let a=delta/2, r=1-mu, p=1+delta and po=1+2mu. The exact original sigmoid primitive is J(t)=integral_0^t sigma(v)dv. The full normalized shape is

`K/KR = exp((a-mu)*t-r*J(t))`.

The actual current native X history supplies A/KR=(K/KR)*X. The actual complete backward energy supplies

`E/KR^2 = exp(delta*t)*(after_entry + remaining_entry_energy(t))`.

The native energy packet has the half-energy normalization; the conversion is exactly E/KR^2=2*(K/KR)^2*energy_Taylor. The complete current preheat Gamma term and both epsilon atoms remain included in after_entry. Nonzero cumulative angular history remains in X.

The current original analytic Cp=0 identity already gives PR=-Ev2*thetaR^2*Jpost. Entry prefix pressure is the original kernel Ip(t), so

`P/KR^2 = exp(p*t)*(Jpost-Ip(t))`.

For stable evaluation, this uses the identical remaining expression: integral_t^1 original entry pressure plus the existing power/exit/waiting/full-heat future contributions. The suffix integrand is replayed against the production prefix AST, including its 1/2 factor and po rate. Original sigma/J and positive correlated cell lengths are retained. No new gauge, pressure constant or division by an amplitude cap is introduced.

The original entry signed-defect AST is replayed before removing 1/k, 1/delta and 1/(2p), giving 45 full-moment normalized mixed4 identities. A further 45 identities establish the actual native A=KX recurrence and exact equality between the original correlated entry inertial-stress AST and the general full-moment operator. Nonzero full E/P functions remain arbitrary in this proof. The previously checked 40 original stress baseline/homogeneity identities retain exact Qtheta*KR and Qz*KR^2 units.

## Physical tensor and actual interface

The general-K physical assembler is replayed directly from the checked current angular implementation AST. It preserves the physical lambda powers, viscosity factors, moving cylindrical basis, exact source radius and the completion Ttheta_theta=r*partial_z(Trz). General K_Z terms are allowed; here the entry K is source-independent of Z, while X/E/P retain their axial histories. The exact inverse-radius source is S_exact*exp(-q); its positive cap is only an enclosure. The leading remainder remains Er=Ez=0 and Etheta=-nu*partial_zz(utheta), with no global temporal-flatness claim.

The angular-entry join consumes the current postpulse source proof's 60 primitive mixed4 identities, including the current full-energy bridge, current XR/PR terminal assignments and current Cp=0 pressure theorem. It then instantiates the original arbitrary-terminal full-moment stress/pressure join and both current normalization/correlation proofs. The common ordinary logR/radius and unchanged general-K physical operators transfer the function equality to stress3, completed diagonal2, divergence2 and remainder2. Interval overlap is not used as the equality proof.

Five view calculations cover the whole entry interval, both boundaries, a middle spatial/time/viscosity sector and a fresh angle/time/viscosity point. They check 325 factored physical contributions, of which 310 have nonzero enclosures. The actual tensor interface checks all 65 common rows (62 nonzero enclosures), and a fresh compact sector checks another 65. These counts measure scoped coverage, not whole-project completion or mathematical blowup.

Read-only review metadata: GPT-5.6 Luna / max. Scoped review passes: current energy/pressure/moment rows resolve the prior old-terminal gap, and the original arbitrary-terminal AST join is instantiated with the current source functions. Direct suffix pressure uses the same original integrand and FTC additivity; a separate explicit t=0 suffix/prefix source identity is optional provenance hardening, not a blocker. The reviewer confirmed the physical inverse-radius rebase is an enclosure and the regional/global distinction remains explicit.

Exact working/Git-index source and receipt hashes pass for all 652 dependencies. The real focused producer/checker run preserves one graph, the 14/8 velocity-pressure inventory and one actual completed-tensor join. Compilation and Git whitespace checks pass; global/temporal gates remain false.

## Executable next construction

- [x] F57C5a-postpulse-entry: current complete-moment entry stress using original sigma/J and analytic pressure.
- [x] F57C5b-entry-completion/F57C6a-entry: actual regional completed physical tensor and residual/remainder decomposition.
- [x] F57C5b-angular-entry-join: current source-function tensor join and common physical bounds.
- [ ] F57C5a-postpulse-power: recover current steep-power full moments and actual tensor from the same complete history. Use the original source t=Ts*phase, distinguish phase derivatives from ordinary logR derivatives, and keep attenuated amplitude units in exact source logs. Preserve the native 1/4 half-energy floor and absolute pressure.
- [ ] F57C5b-entry-power-join: instantiate full-moment and completed tensor equality at entry t=1 / power phase=0, with the same actual XS/PS, full energy and radius. Bound stress3/divergence2/diagonal2/remainder2 without using interval overlap.
- [ ] F57C5a-postpulse-exit-waiting: recover actual steep-exit and waiting tensors from the current power/exit/heat histories, including original sigma derivatives, current Ts/wait and complete Gamma future. Then admit power-exit and exit-waiting tensor joins.
- [ ] F57C5a-postpulse-flatten-power-collar: attach current flatten/outer-power/heat-collar full stress to the same owner, with the complete actual exterior zero-stress theorem. Keep nonzero moments, original pressure and source phase Jacobians.
- [ ] F57C5a-pulse: actual entrance/main/exit/gap/end tensor from five histories, selected C1/C2, meridional cross terms, complete future/2 and Pin/P0; local support differences are not actual tensor values.
- [ ] F57C5b/F57C6a-global: finish all actual tensor interfaces and compose every chart's physical decomposition; count actual tensor joins separately from the current22 velocity-pressure atlas.
- [ ] F57C4e-remaining/F57C6b/F57D/E/F: other retained/core/incoming joins and axis attachment, independent temporal remainder/cone/lift/energy and resolved nonzero physical u,v,w,p.
- [ ] F58/F59/F60/F61: actual n-dependent recovery equations/per-order repairs/remainders/summation, both oscillatory families/mean corrections, corrected NS and measured contraction/slenderness/winding.

Save each next construction's current source evidence and scoped receipt, mark the precise task done and commit/push. The full reconstruction goal stays active.
