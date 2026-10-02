# Quantitative pulse flat comparisons, both external joins and O5 flatten

The same compliant `.001` source now has quantitative support-endpoint
moment/velocity differences, independent O3 power inlet derivatives, both
external O4 joins and mixed derivatives throughout the original 100-unit
O5 flatten. This closes the **local leading O4 pulse spatial C4** gate.
Whole outer/core/axis C4, time derivatives, physical kinetic energy,
admissible stress, stress remainder and temporal coefficient recursion are
still incomplete.

## Same histories, quantitative differences

`lei_ren_part1_paper_compliant_pulse_flat_comparison.py` compares the actual
field to its zero-local-axial-input continuation with the **same exact
nonzero interface histories**. It never deletes incoming moments or the
selected positive terminal future energy.

For local log-radius distance `h`, the main support width is `mu*h` and
the end-bump variable satisfies `w <= 2*h/.15`. Original sigma/gp/beta
analytic flat bounds and the actual selected axial C5 functions bound the
forcing derivatives. The linear difference bounds are
`h*exp(lambda_i*h)*||B||`, valid in both orientations; the energy difference
is bounded by `h*exp(2*mu*h)*||B^2||`. The common `-.5` energy source
cancels in this difference only. Repeated derivatives use the original
ODEs, including every binomial derivative of `B^2`.

Paper (3.9), axial C5 input and the differentiated physical prefactors
produce all 15 velocity mixed orders through4 and the admitted cylindrical
physical r/z maps. Thirty endpoint/distance packets cover entrance, exit
and both edges of both end bumps. All upper bounds decrease toward the
endpoint, including the positive numerical tail caps.

The flat-limit assertion holds for each fixed actual source family and
fixed positive physical `tau`. It is not a uniform `tau -> 0` estimate,
and it is not the independently required stress flat remainder.

## Independent O3 inlet and exact source shapes

`lei_ren_part1_paper_compliant_power_inlet_C4.py` independently reconstructs
the left neighborhood `y=log(R/Rp) in[-1,0]` inside the actual long power
buffer. It rebuilds canonical incoming primitives from their constants,
evaluates the original analytic pressure datum, and transports the power
solution before recovering velocity/pressure mixed derivatives.

The checker translates the **production source expressions** through O2
slope, axial turnoff, O3 slope-mu and O3 power. These identities establish
the following exact whole-Z shapes with `q=1+Z^2`:

```
u=U/q,       Mz/R=M*Z,       mixed=K*Z/q,
angular=H/q, energy=E_Z*Z^2+E_Q/q^2, Mp=Pin/q^2.
```

All radial kernel call sites are independent of Z. Consequently the stored
`Z=.5` sample encloses the proved constant functions H and Pin; it does not
fit or extrapolate their axial dependence. The same-source canonical
incoming functions, flat pulse forcing and shared ODEs identify every
required mixed derivative across the inlet. Numeric interval overlap is
reported separately as a diagnostic, not as the proof.

Zero axial input does **not** generally imply zero radial velocity before
the pulse: the nonzero accumulated Mz history is preserved.

## Entire O5 flatten and terminal join

`lei_ren_part1_paper_compliant_flatten_mixed_C4.py` supplies the original
flatten for every `Z in[-1,1]` and `t=log(R/Rv) in[0,100]`. Set
`rho=log(q/2)` and `F=exp(rho*sigma(t/100))`. In Ev0 velocity units,

```
theta = F/q * exp(-(.5+mu)*t)
X_t   = 1 - (1-mu+rho*sigma_t)*X
e_t   = -.5 + (2*mu-2*rho*sigma_t)*e
P_t   = (Ev0^2/Pstar^2)*theta^2/2.
```

Angular/energy/pressure values use directed source integrals with exact
exponential weights. Whole-time boxes retain cell length `t/cells` and
remaining angular time `t*(cells-i-1)/cells`, avoiding invalid subtraction
of correlated intervals. Variable-rate ODE differentiation supplies every
mixed order<=4; sigma derivatives0..4 suffice because C5 here means
**axial** C5 input, not fifth y derivatives.

The terminal energy is half the **complete infinite future energy** in
`Rv*Utheta(Rv,Z)^2` units. Thus its full remaining integral in `Rv*Ev0^2`
units is `2*ev/q^2`. No q factor or factor2 is discarded. This energy stays
positive throughout flatten. The future-energy derivative boxes are bounds
on the same source function, not a polynomial fit to be evaluated at other Z.

Pressure remains original `P0+Mp`; formal Ev0 and the omitted exponentially
small pressure factor retain their positive defining logarithms and proved
caps. Their exact functions are not replaced by zero or by a cap.
The computed terminal pressure and scale values are **interval enclosures**
of those original exact sources. The pressure-decay log keeps its finite
`-26` offset separate from `-13/mu`.

At Rv the original sigma and its y derivatives vanish, so the primitive
ODEs match the actual pulse terminal ODEs. Linear and radial zeros are
inherited from empty future supports. At the flatten exit, sigma=1 with
flat derivatives, hence `F/q=1/2` identically and swirl is exactly Z
independent. The following power region still needs its independent
provider and joining certificate.

## Evidence and reproduction

- Flat comparisons: 360 independent closed-form derivative checks in both
  orientations, 4200 finite bounds, 840 exact-zero difference bounds and
  1920 monotone upper-bound checks pass.
- O3 inlet: 78 source/join identities, 840 actual mixed bounds and 60
  independently evaluated inlet overlap diagnostics pass.
- Flatten: 240 independent closed-form variable-rate derivative checks,
  900 actual mixed bounds, 450 inherited exact-zero bounds and 72 terminal
  overlap diagnostics pass. Functional identities separately prove the
  source units, rates, flat endpoint compatibility and original expressions.
- A reused read-only GPT-5.6 Luna/max worker confirmed the flatten rates,
  energy normalization and sigma derivative order. Its source-shape concern
  prompted the whole-Z production source proof described above.

Run from the repository with the experiment directory on Python's import
path (the script entry point handles this):

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage flatcomparison
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage externaljets
```

The complete ordered pipeline has 68 modules. The accepted checker JSON
is authoritative for joining/acceptance flags; the producer
report deliberately leaves those flags false until the checker passes.
The `full_pulse_C4_installed` flag refers only to the leading O4 pulse and its
two **local spatial** external interfaces at fixed positive tau, with
profile mixed order<=4 and the admitted physical cylindrical r/z map.
It does not certify time derivatives, the whole outer field, the core or
axis, Cartesian vector derivatives, physical energy, stress, recursion,
oscillatory correction or the complete forced NS residual.

Next: following power and both angular supports, steep/waiting/collar/Gamma
mixed derivatives and interfaces; then whole-field Cartesian derivatives,
physical energy and stress/flat remainder; then genuine n-dependent
coefficient recovery.
