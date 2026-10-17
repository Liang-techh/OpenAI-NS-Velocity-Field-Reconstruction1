# Direct correlated scaled sources and exact Rc inlet identities

Tasks 14/15 remain **IN PROGRESS**. The long-term goal remains **ACTIVE / INCOMPLETE**. This stage replaces independent full-minus-background input subtraction with direct signed mean-history recovery, separates physical scale sectors, constructs the actual signed canonical collar primitive, proves 39 exact Rc inlet function identities, and resolves the true P0 pressure-source identity. It does not evaluate the original interior finite-N integrals or remove any aggregate cap.

## Direct same-history recovery

The earlier source graph expressed the correction as full physical U/Q minus its background. Evaluating those two inputs in independent interval boxes would lose the common history. The new `original_correlated_scaled_source.py` instead uses the actual signed normalized deltaV and delta-m history ports before applying any enclosure arithmetic.

With D=1-Z^2 and L=1-delta*Z^2:

```
delta-q = [2Z*deltaV - Z*(1-delta)*delta-m - D*(delta-m)_Z]/L,
delta-q_y = [Z*(1+delta)*deltaV + 2Z*deltaV_y - D*deltaV_Z]/L - delta-q.
```

Physical U=Pstar*V and Q=Pstar*q. The kernel `original_scaled_signed_density.py` returns the exact ordinary-C1 sectors:

```
Delta v = Pstar*v1 + Pstar^2*v2,
Delta F = Pstar*F1 + Pstar^2*F2,
Delta Q_y = Pstar*delta-q_y,
Delta boundary = R^2*(Pstar*F1+Pstar^2*F2)-R*Pstar*delta-q_y.
```

The linear and quadratic sectors retain every correlation and square from the accepted six-term pressure kernel. Neither Pstar nor physical radius nor the exact N is materialized. Ordinary C1 density output needs delta-m through ordinary Z3, because delta-q_ZZ differentiates the mean recovery. Missing derivatives and truncated records are rejected, not padded with zero.

| actual chart | signed deltaV Z0/Z1/Z2 | signed delta-m Z0/Z1/Z2/Z3 | local normalized graph |
| --- | --- | --- | --- |
| compact | 260375 / 260378 / 260382 | 261593 / 261594 / 261595 / 261596 | 332 nodes |
| Rm patch | 198182 / 198183 / 198184 | 198466 / 198467 / 198468 / 198469 | 350 nodes |

Fourteen actual row identities verify full V=background V+deltaV and full m=background m+delta-m before pressure recovery. Rm deltaV_y is formed directly as slow_B_y/N+B_over_Pstar_phi from its original total radial formula. This is the original radial N division, applied once to the slow term; the already scaled deltaV and histories are not divided by N again. Compact deltaV_y uses its saved signed profile derivative rows directly.

The normalized arithmetic has exact rational, SymPy and directed interval backends. New rational source-port controls check all six output roots on both real chart schemas. These fixtures are not values of the original functions, and their small synthetic N is not a surrogate for the original DoubleDyadicInteger. Actual signed Rc targets and interior histories remain unevaluated. Their saved covers cannot be used as function values.

## Actual signed canonical collar primitive

`original_collar_signed_primitive.py` binds the actual lower graph rather than leaving the singular-width remainder as a formula string. It creates 99 local graph nodes with ordinary Z0/Z1 pairs for:

```
A(s)=2*rho(s)^2*Z*(1-hb)*ratio(s)*Psi_core_rho(rho(s),Z)/L,
A0=32*Z*Psi_core_rho(4,Z)/L,
term4=epsilon^2*A0*sigma(s),
term5=epsilon^2*integral[t=0..s] (A(t)-A0)*sigma_prime(t) dt,
full singular sector = term4+term5.
```

The core radial source is node35715 for ordinary Z0 and node35833 for ordinary Z1. Both are radial order one. The new local substitutions use core rho variable35558 and constant4 node13. Nodes36065..36068 are radial-order-zero Psi(4) and are explicitly excluded. Exact original hb8, epsilon35565, rho35851, ratio35940/35941, Psi_core_rho35870/35871, L30382/30381 and sigma_prime35953 are pinned.

The local integration variable t is distinct from the outer query s. Epsilon squared occurs once outside the integral. A0 is independent of s. Source guards and independent rational polynomial integration/derivative controls pass. This constructs the canonical signed remainder function; its numerical value and the forward-to-canonical six-term source equality are still unproved. The saved width norm E remains a bound, not its integrand or signed value.

## Exact Rc inlet source normalization

`original_compact_inlet_identity.py` verifies that the complete outer graph is an identical source prefix of the compact graph. The compact x variable is introduced after that prefix and does not occur in the prefix. Source atoms from this verified prefix therefore remain fixed while the compact endpoint is substituted at x=1.

The bounded rewriter applies only exact rational arithmetic, zero-length integrals, log(1)=0, exp(0)=1 and the actual compact x-to-1 substitutions. It does not execute original functions or infer equal values from overlapping boxes. It reduces 38 local integrals to zero and proves:

```
compact.inlet.correction[key] = outer.O3_power(2).correction[key],
  key=m,h,k,e,p, ordinary Z orders0..5;
compact.inlet.integrals[S0,S1] = outer.O3_power(2).integrals[S0,S1],
  ordinary Z orders0..3;
compact.inlet.radius = outer.O3_power(2).radius.
```

There are 39 identities in total: 30 correction rows, 8 pressure-integral rows and one radius. All reduce to the identical actual source atom, with no unresolved rewrite dependency. Shared source prefix digest: `7ec4add5799967f01186af73b98d1c64ef45462b8a2c13ace2e24fc7ca3e7ad5`. Counterexamples keep distinct source atoms distinct and retain the unresolved inlet rather than setting it to zero. A symbolic retained inlet is not a numerical proof that its value is nonzero.

This is a new source-level bridge for cumulative histories and S0/S1. It does not yet prove the complete F/Q_y pressure-boundary pairing at Rc. Neither archive already contains a DeltaF/DeltaB theorem. The compact upper terminal zero was already consumed by the accepted numerical source and is not counted as new work. Aggregate endpoint and volume caps stay active.

## True core pressure source and all-order axial q coefficients

The new `build_core_pressure_identity` binds the saved native analytic pressure theorem to the current forward seed and canonical core constructor. Both true functions use all fourteen original pressure sources:

```
p(Z)=P0_physical/Pstar^2=-(M2/(1+Z^2)^2+M0+F_flat(Z)),
canonical: P0_physical=Pstar^2*p,
forward: P0_Z_taylor=epsilon_core*P0_physical_Taylor,
epsilon_core=exp(-logLambda).
```

The core radial epsilon is distinct from the pressure-datum epsilon parameter. The source and datum hashes, the analytic pressure-definition hash, source ASTs and previously accepted dimensional identities are checked. M2, M0 and F_flat are true source functions/scalars; their saved intervals remain enclosures. The forward point approximation does not become the true function. No mass endpoint, midpoint or flatten bound is selected as a value.

For q(Z)=(1+Z^2)^-2 at center Z0, write a=1+Z0^2, b=2Z0 and q(Z0+t)=(a+b*t+t^2)^-2. The independently checked identity

```
(a+b*t+t^2)*dq/dt + 2*(b+2t)*q = 0,
c0=a^-2, c[-1]=0,
a(n+1)c[n+1]+b(n+2)c[n]+(n+3)c[n-1]=0
```

matches the actual `q_jets` source for every n>=0. Since a>0 at real centers, local analytic uniqueness proves all Taylor orders, extending the inherited finite-prefix q-component proof without rerunning q_jets or the core recurrence. This extension is only for the rational q component; F_flat derivatives remain their source-owned Cauchy/enclosure rows, not newly evaluated exact coefficients. Ordinary derivatives are k!*c[k]; these are axial Taylor coefficients, not temporal Navier--Stokes coefficient recursion.

This resolves the P0/M2 *function-definition* gap recorded in the inherited parent adapter. The parent adapter remains unchanged as provenance; its older P0-unproved entry is superseded by the new pressure-source proof. Equality of the full two core numerical producers and the complete six-term source is still unproved.

## Actual run and preservation

New production exits0 in 6.363 s. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json` is 3,906,687 bytes with SHA256 `9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444`. It contains 1971 input pins and 275 implementation pins. Saved-only inspection validates all2246 pins, the 14 same-history row identities, the 39 Rc inlet function identities and packet invariants. Python compilation passes.

Exact controls cover the symbolic Pstar split, ordinary-C1 chain rule, pressure recovery linearity, missing-derivative rejection and signed interval arithmetic. Twelve new rational port-shaped controls cover delta-q, delta-q_y, v1/v2 and F1/F2 on compact and Rm. The collar controls cover the ordinary quotient rule and independently assembled polynomial signed split and primitive derivative. No accepted owner, producer, constructor, fixture, core recurrence or quadrature is replayed.

All three179-entry discrepancy views remain unchanged. The negative-pulse target and quiet pair were already consumed; do not subtract either again. Original raw213/reduced183 archives, all six signed core/lower contributions, exact N, four signed leading leaves19/20/41/42, pressure throughRv and all finite-N caps remain. Canonical joint slot182 is still -(inner+lower+active+compact). P1/m13 is distinct from m15.

## Next executable work

Attach the normalized delta-q/delta-q_y and DeltaF expressions to the Rc join. Prove actual full velocity/profile and radial continuity at this same source endpoint, then use the shared correction histories to pair the complete boundary R^2*DeltaF-R*DeltaQ_y before numerical boxing. The new cumulative-history identities alone do not authorize deleting the aggregate cap.

Build actual signed Rc target and Duhamel-history enclosures for the normalized evaluator; current control functions are exact source definitions, not numerical input records. Preserve the original phase, N, slow/fast radial terms, quadratic pressure history, physical weights and nonzero incoming pressure. The true P0 function identity is now bound; next use it to compare the full core solution/response operators, then complete the forward q-dependent collar profile/history bridge. The canonical signed primitive is ready for that comparison.

Global functional all-five/common-N closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion and smooth summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS1e-3 gate is not applied to this stage.
