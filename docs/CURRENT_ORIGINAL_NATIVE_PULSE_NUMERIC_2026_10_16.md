Latest continuation: [original repaired pulse P1 throughRv](CURRENT_ORIGINAL_NATIVE_PULSE_PRESSURE_2026_10_16.md). Use `OriginalNativePulsePressureVelocity` for pressure orders0/1. The base API documented below keeps its explicit leading-pressure-only scope.

# Original selected Rp-to-Rv pulse numerical velocity - 2026-10-16

`OriginalNativePulseVelocity.evaluate(Z,chart,coordinate,log_tau="-2",azimuth=".25")` now connects the actual quietRp numerical inlet to the original selected pulse and evaluates its entrance/main/gap/end velocity throughRv. It carries the original own moment/energy histories and leading pressure, with source-defined positive pulse selection and two narrow end corrections. The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Interior repaired P1 is explicitly unavailable in this pulse API, while the complete incoming P1 atRp is retained separately. The preceding compact/quiet APIs still provide their full P1 enclosures. Genuine temporal n>=2, full stress/heat matching and corrected NS remain open.

## Actual eight-function Rp interval callback

The reader consumes the accepted quiet numerical report at q=1 and the checked original interval-inlet, canonical native C5 source identity, original selected pulse, refined coefficient and fixed pulse-basis/energy archives. It does not construct another selected/future branch or replay the numerical ancestors. With Q=1+Z^2, the exact source functions are

```
u=U/Q
m1=C1*Z*Q, m2=C2*Z*Q
X=Xp
energy=C0+C_E*Z^2*Q^2
Mp=Pin/Q^2
pressure=P0+Mp.
```

The common constants U/C1/C2/C0/C_E/Xp/Pin are Z-independent source enclosures. u,m1,m2,X,energy andMp are reconstructed as actual C5 canonical functions, and their ordinary rows0..2 are intersected with the current independently continued quiet endpoint. The exact source identity permits this intersection; overlap alone is not treated as proof of source identity. Higher rows remain those of the same proved canonical functions, not zero padding or derivatives of interval endpoints.

The quiet packet already retains the six original independent P0/Pstar^2 rows, even though its simplified native frame retained only three. The new callback consumes all six, intersects the first three with the current frame and forms pressure=P0+Mp at C5. All eight incoming functions therefore have six genuine interval Taylor coefficients. No datum is reconstructed and no Pstar is applied during this normalized callback. The original selected ap and two scaled end coefficients retain their checked C5 enclosures on the four original whole-Z cells, preserving their common positive branch, incoming moments, full future energy and original source family. The coefficients remain bounded functions, not chosen interval endpoints or resolved signed point values.

The source weighted future equals mu*exp(-26)*complete_future/2. The reader converts it back by exactly2*exp(26)/mu for the existing positive terminal-energy formulas; it does not reinterpret a weighted future as an unweighted one. Kpulse is the fixed positive gp energy integral, distinct from the native moment constantK.

## Unchanged original segmented integrals

The pure original `CompliantAxialPulseField` methods are reused with this current data callback. No original selected cache or module object is overwritten. The profile is

```
gp(xi) = (integral_0^xi sigma(50v)dv)*sigma(11-xi)
Utheta = Pstar*u*exp(-t/2-xi), xi=mu*t
Uz = Utheta*ap(Z)*gp(xi)  # main
m_i(t) = exp(-lambda_i*t)*m_i(0)
         +ap*integral_0^t exp(-lambda_i*(t-v))*gp(mu*v)dv
lambda_i=1/2-i*mu, i=1,2.
```

The original forward main integrals retain the inherited moment terms and positive omitted-tail bounds beyond their finite integration window. The incoming and nonzero history terms are not reset in the inactive gap. Original backward bump primitives supply end moments with centers-3,-1 and half-width0.15. AtRv the future bump support is empty, giving the selected source linear terminal identities. This is not an independent signed residual evaluation of precisely resolved coefficients.

| chart | coordinate | original domain |
| --- | --- | --- |
| entrance | t=log(R/Rp) | t>=0, with mu*t in the main support |
| main | xi=mu*log(R/Rp) | 0<=xi<=11 |
| gap | xi=mu*log(R/Rp) | 11<=xi<13 and distance toRv>=4 |
| gap_end | s=log(R/Rv) | -1/mu<=s<=-4 |
| end | s=log(R/Rv) | -4<=s<=0 |

Only Z-independent scalar kernels are privately memoized. Copies of the original method namespaces route gp, gp energy, future gp energy and backward beta integrals through caches keyed by complete directed scalar endpoints, context and integration arguments. The defining code, source parameters, partitions and directed values are unchanged. The current Z-dependent callback, selected jets and future source are never cached as scalar kernels. This avoids recomputing identical radial integrals for each Z while preserving the source functions and their errors.

## Dependent recovery and original physical units

Let B=Uz/Utheta and m1=Mz/(R*Utheta). With L=1-delta*Z^2, d=1-Z^2 and Utheta_Z/Utheta=-2Z/(1+Z^2),

```
A = [2Z*B-(1-delta)Z*m1-d*(m1_Z-2Z*m1/(1+Z^2))]/L
Q0 = Utheta*A.
```

The positive n0 recovery weight is1-delta. The source-derived m1 derivative is retained. A numerical interval center is not differentiated as an exact global field. Current meridional recovery uses the same source streamfunction; physical Cartesian gradient/divergence diagnostics of interval centers are still unadmitted.

The physical terms are

```
utheta = lambda^(-1-delta)*Pstar*u*exp(-t/2-xi)
uz = lambda^(-1-delta)*Pstar*u*B*exp(-t/2-xi)
ur = lambda^(-1)*sqrt(Rp/2)*Pstar*u*A*exp(-xi)
psi = lambda^(1-delta)*Rp*Pstar*u*m1*exp(t/2-xi).
```

The radial sqrt(R) and exp(-t/2) cancel before interval boxing. log110 and the actual absoluteRp radius are applied once. In end coordinates, t=13/mu+s and xi=13+mu*s are represented by separate inverse-mu and finite atoms so the finite end offsets survive. The new physical reader uses the original Bhat/mhat end quantities with their common logE retained as separate atoms:

```
logE = -1/mu-3*(4mu)^(-2/3)-(4/3)*logmu+(2/3)*log4-log6/2.
```

It does not select a relaxed exp(-1000) cap as the actual end amplitude. The two selected end velocities retain their opposite signs. The n1 velocity/mass tail is the same source-admitted zero after2Rc. The physical velocity packet includes temporal orders0/1, where the n1 velocity terms are exact zeros on this region; the main/end pulses here belong to the leading spatial background, not temporal n>=2 recursion.

## Pressure scope

Leading pulse pressure is recovered from the original angular source, incoming Mp and the independent P0 once. It has physical factor lambda^(-2-2delta)*Pstar^2. Positive tiny decay caps remain directed errors, not exact zero sources. The full incoming repaired P1 atRp is retained as metadata for the next integral layer. It is not held constant, zeroed, or incorrectly appended as the interior P1. Accordingly the current physical pressure packet declares order0 only and `interior_P1=False`; global matching and relative pressure cancellation are still open.

## Focused evidence

The final focused run exited0 in130.531s. Saved-report inspection also exited0 and verified1854 scientific input hashes and78 implementation hashes. The18,378,400-byte report has SHA256 `08bf5e501b343ec41de0d1b41c6ce8fe956a7274a693b4c343bbdb4a5c29582c`.

- 11 independent original profile and Fubini convolution comparisons;
- 39 actual physical packets across Z=-0.3,0,0.3 and all original entrance/main/gap/end charts;
- 72 current native callback/same-N conditions and111 actualRp coefficient seam conditions;
- 24 checks that all8 callback functions retain genuine C5 rows;
- 9 nonzero main/opposite-sign end velocity conditions and48 gap-history/incoming-P1 conditions;
- 72 source-admittedRv terminal linear zero conditions;
- 186 physical pressure-scope/atomic-end-unit conditions,21 scaled-time conditions and5 rejected invalid requests.

The saved command is `python scripts/run_original_native_pulse.py --inspect-report`. Use `--write-report` only when a new implementation requires a new focused run. A request-stream disconnect is not a reason to replay this successful computation.

The existing read-only `/root/n1_atlas_map` reports GPT-5.6 Luna/max. It reviewed the actual native normalization, selected branch/future conversion, original integrator adaptation, physical radial/streamfunction weights, atomic end scale, preserved six-row P0 datum and private scalar cache ownership. It found no remaining material source/scaling issue within this scope. No new agent was spawned; scientific ancestors remain unchanged.

## Next construction work

1. Numerically integrate the original fixed leading N1p across each partial pulse chart to recover interior P1, starting from the full saved repaired P1 atRp. Reuse the current local U/m/Q ordinary rows, true physical measures and finite viscous endpoints. Full-pulse signed bounds or cancellations are not partial-prefix values.
2. Use the radial-free pressure identity to separate linear source integrals, quadratic error integrals and finite endpoint terms. Preserve q_Z, U_Z, q_ZZ as required for ordinary C1. For the main forward primitive, the same own-rate equation yields an index-derivative convolution for integral m dR; any implementation must independently verify its original factors and lower term.
3. Near the terminal end use a backward pressure representation only with the same actual full repaired pressure terminal theorem and the complete remaining forcing. Pressure zero is allowed by the source theorem only beyondRv; do not infer it from zero local velocity or an empty partial bump.
4. Then continue the original flatten/end/heat construction with all five histories, compatible analytic pressure and actual selected branch. Close physical energy/radial tail, high mixed/global joins and common-global-N admission, admissible stress and independently measured flat remainder.
5. Build genuine independent temporal n>=2 coefficient recursion/smooth sum, two-family oscillatory correction and independent corrected Cartesian NS validation. Keep all full-goal gates open until evidence proves them.
