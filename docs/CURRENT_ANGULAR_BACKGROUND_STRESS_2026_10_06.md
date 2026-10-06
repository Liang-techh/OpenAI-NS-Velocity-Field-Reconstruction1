# Current actual angular background stress (2026-10-06)

The original angular region now has the actual current full-moment background stress, its prescribed symmetric tensor completion and the regional physical decomposition `residual = -div(T_B) + E_B`. This completes F57C5a-angular, F57C5b-angular-completion and F57C6a-angular. These outputs are source-bound directed enclosures with exact positive logarithmic factors. They are not resolved coefficient points or global tensor/cone/flatness/NS certificates.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_angular_background_stress.py`, matching `_check.py`, producer `.json` and scoped `_check.json`. Focused controller stage: `currentangularstress`. It consumes the checked current22 interface atlas; all 14 affected adjacent and 8 internal velocity/absolute-pressure source traces remain admitted.

## Callable result

```python
from lei_ren_part1_paper_compliant_current_angular_background_stress import CurrentAngularBackgroundStress
stress = CurrentAngularBackgroundStress(atlas=checked_current_interface_atlas)
packet = stress.angular(Z=".537", offset="-2.97", log_tau="-2.6", theta=".41", viscosity=".8")
```

The original source domain is s=[-4,0], Z=[-1,1]. Finite physical points require R>0, |Z|<1 and tau>0; Z=+/-1 are infinite-space/source limits, not axis points. Set `theta=None` for all angles. Any finite compact log(tau) interval and finite constant nu>0 are supported. There is no uniform bound at tau=0 or axis attachment in this regional result.

Outputs retain the full normalized A/KR, E/KR^2, P/KR^2 and K/KR rows through radial order4/axial order5; stress mixed3; physical stress divergence, tensor diagonal completion and axial-viscosity remainder mixed2; Cartesian tensor components, divergence, remainder and residual decomposition; and the stable factored absolute-pressure mixed4 companion. The pressure companion retains its signed negative coefficients and the positive source factor `exp(2*log(Bcurrent))`; it is not a numerical physical pressure value. Exact KR/KR^2 and amplitudes remain source factors, never selected values or divisions by amplitude caps.

## Full histories and pressure

The same checked selected repair, complete future energy, original analytic P0/Pin, current angular signed C5 controls, original beta functions and complete Gamma pressure tail are retained. No legacy angular-stress constructor supplies terminal constants. `angular_future_changes` contains both linear and quadratic terms in the two disjoint supports.

The two pressure rates are distinct: heat normalization p=1+delta; angular beta-kernel rate po=1+2mu. With original source amplitude ratios

`rS=exp(-2bp-rate)`, `rQ=rS*exp(-3Ts)`,
`rT=rQ*exp(-3+k)`, `rB=rT*exp(-2bh*wait-2logone)`,

the complete normalized post-angular pressure is

`Jpost = Iin + rS*I3(Ts)/2 + rQ*Iout + rT*Ip(wait)/2 + rB*remaining_heat_pressure(0)`.

The actual current PS/PQ/PT, waiting pressure and Cp assignment expressions are replayed symbolically. They give `Cp = PR + Ev2*thetaR^2*Jpost`. The accepted current original analytic-pressure function theorem has Cp=0, so the original terminal pressure is `PR = -Ev2*thetaR^2*Jpost`. Jpost replaces the datum through this identity; it is not added to an already absolute pressure. The original P0/Pin remains included once.

At angular offset s,

`P/KR^2 = exp(p*s) * [Jpost + (exp(-po*s)-1)/(2po) + remaining_beta_pressure(s)]`.

The original full angular defect-row AST is replayed on arbitrary nonzero terminal A/E/P functions, keeping the 1/k, 1/delta and 1/(2p) baselines. Its 45 normalized mixed4 identities match the new helper. The original stress AST then cancels those exact unit baselines before extracting KR and KR^2; 40 mixed3 cancellation/homogeneity identities pass. General K_Z and K_ZZ, full moment FTC and the moving physical-coordinate operator are retained.

## Tensor and bounds

The symmetric tensor has Trtheta=Ttheta, Trz=Tz and Ttheta_theta=r*partial_z(Tz), with the other independent components zero. Its radial divergence cancels exactly. In cylindrical components the leading remainder is Er=Ez=0, Etheta=-nu*partial_zz(utheta). The accepted general-K Cartesian/moving-basis/viscosity source theorem supplies the regional physical residual decomposition.

Positive source logs are reduced analytically before enclosure. Ts, waiting and log(1-epsilon) cancel against the exact KR normalization, keeping the physical Qtheta*KR, Qz*KR^2 and completed diagonal factors without materializing KR. The exact inverse radius is S_exact*exp(-q). Its directed `[0,cap]` enclosure enters the original stress operator; it does not define an exact inverse-radius field. All normalized and physical rows are enclosures of their source functions.

Seven view calculations cover the whole angular region, both original boundaries, both supports, a fresh point inside the first support and a fresh point between supports, with fresh positive-time/viscosity/angle parameters. Evidence: 455 factored physical tensor/divergence/remainder contribution rows, 434 nonzero contribution enclosures, 45 original full-moment normalization identities and 40 original stress baseline/homogeneity identities. Invalid domains, nonfinite time/angle and foreign current angular/pressure owners are rejected. The focused controller uses one actual producer/checker object under the retained checked graph; all global and temporal gates remain false.

Read-only review metadata: GPT-5.6 Luna / max. The reviewer found no sign, rate, tensor factor or physical viscosity error and confirmed the full-moment normalization. It emphasized the enclosure/factored-value distinction, which is explicit in the scoped receipt and this handoff.

Exact working/Git-index source and receipt hashes pass for all 648 dependencies. The actual focused producer/checker run preserves one shared graph, the 14/8 interface inventory and the regional gates; global/temporal gates remain false. Compilation and Git whitespace checks pass.

## Executable next tasks

- [x] F57C5a-angular: actual full A/E/P/K stress, exact KR/KR^2 units, current source pressure and complete future.
- [x] F57C5b-angular-completion: actual regional symmetric tensor and physical mixed derivative bounds.
- [x] F57C6a-angular: regional source-bound residual=-div(T_B)+E_B with the original axial-viscosity sign.
- [ ] F57C5a-postpulse-entry: build the current steep-entry full stress from the same angular terminal functions. Replay the original transition kernels and preserve exact normalization; prove the angular/entry tensor join on functions before quantitative bounds.
- [ ] F57C5a-postpulse-rest: build current steep-power/exit/waiting/collar and flatten/power full stresses from the checked history owner, preserving nonzero E/P and source phase Jacobians. Reuse the exact full Gamma exterior stress-zero theorem within its existing scope.
- [ ] F57C5a-pulse: recover actual entrance/main/exit/gap/end stress from the five cumulative histories and selected C1/C2, including m_actual cross terms, complete future/2, Pin/P0 and positive source logs. Do not substitute the checked local difference for the actual tensor.
- [ ] F57C5b-joins: prove actual completed tensor joins throughout these chart chains and all eight support edges. Bound common traces over the declared physical sectors; the current22 velocity/pressure atlas is a prerequisite, not a tensor-join certificate.
- [ ] F57C4e-remaining-joins: cover other retained/core/incoming chart boundaries and axis attachment before admitting the complete33 quantitative-interface gate.
- [ ] F57C6a-global: compose regional physical decompositions for every current chart with shared absolute pressure, basis and viscosity units.
- [ ] F57C6b/F57E/F: independently bound actual global temporal remainder, current cone margins/admissible lift and kinetic energy over the prescribed domain. Local support flatness and compact-time upper bounds do not establish these gates.
- [ ] F57D: solve the source coefficients to directed requested accuracy and provide actual nonzero u,v,w,p point evaluation with domain/error metadata; no interval-center point selection.
- [ ] F58: implement distinct n=1 and n>=2 defining equations, one common core, per-order five-moment/pressure repair, finite-order remainder and smooth summation.
- [ ] F59/F60/F61: both oscillatory families and mean corrections, averaged quadratic stress cancellation, independent corrected NS and measured contraction/slenderness/winding.

Complete the next bounded construction, save its current source evidence, mark its scoped task done and commit/push. The long-term reconstruction goal remains active.
