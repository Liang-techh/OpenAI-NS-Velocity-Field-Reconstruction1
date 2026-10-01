# F18: larger-Md outer stress screening

The O.2 midpoint remains an actual stress-direction obstruction for Md=2,
3,4,5,6 on the local axial family [0.49,0.51], including Z=0.5. No fresh
core was generated for these parameters. The new results select the next
parameter route; they do not invalidate the underlying paper construction.

## New evidence

`experiments/root_st073/lei_ren_part1_paper_outer_tail_pressure_probe.json`
resolves the previous Md=2,3 midpoint noncertificates. Pressure is evaluated
by splitting the SAME fourteen-stage preheat integral at the current radius,
retaining all remaining atoms. No fitted or appended pressure is used.

`experiments/root_st073/lei_ren_part1_paper_large_Md_screen.json` rebuilds
the fourteen-stage datum and nominal waiting root for Md=4,5,6. It chooses
log(Pstar)=exp(Md)+11 and delta=min(1e-200,exp(-4log(Pstar)-30)), and checks
delta<c_delta*mu. Unknown paper constants and original parameter-rounding
errors remain unverified. The angular correction itself is not rebuilt.

At Z=0.5 and cutoff phase 0.5, define x=I_theta/(F R),
B=b^2, Q=b I_z/(F R), b=-S_z/F. The negative dot coefficient D=2x+Q is:

| Md | directed interval for D | leading strong margin/R^2 |
|---|---|---|
| 4 | [-4.13905,-4.03109] | [-44.98676,-42.35340] |
| 5 | [-2.68386,-2.55395] | [-28.36368,-26.02834] |
| 6 | [-1.85261,-1.71287] | [-19.77472,-17.79836] |

Display intervals above are rounded outward. Receipt exact endpoint tuples
are authoritative. Both phase 0.25 and 0.75 pass sampled asymptotic strong
cone checks; that does not establish a complete phase or axial-family proof.

## Algebra and numerical changes

The leading strong margin is evaluated as

    x * [(8-B^2/2)*x + (8+2B)*Q]

instead of subtracting the nearly equal rounded numbers kappa and 2. The
exact finite-placement identity is

    (x-2/R) * [(8-B^2/2)*x + (8+2B)*Q - (B+4)^2/R].

The finite negative-dot coefficient is D-(4+B)/R. Consequently D<0 rules
out every positive radial placement of the SAME normalized outer source.
A leading strong-margin failure alone is reported separately and is not
treated as that all-placement direction proof.

The cutoff integration now encloses only the monotone cutoff on exact
phase cells and integrates exp(s) exactly. This removes the repeated
exponential weight overestimate of rectangle quadrature. Bounded checks
confirm refinement nesting, overlap with the old directed rectangles,
containment of independent SciPy integral values, and twelve exact rational
finite-cone identities. These checks do not prove whole-stage admissibility.

The old CorrectedSourceProfile factory fails at Md>=6 because its separate
finite-float angular bump weights underflow. The pressure screening path
now invokes the same incoming angular ratio, nominal waiting equation,
and matched schedule directly, without constructing those unused bump
coefficients. This bypass is limited to pressure screening. Production
angular/axial corrections still need arbitrary precision weights.

## Next work, in dependency order

1. Implement the normalized propagation equations (6.24)-(6.26) directly
   for O.2. Retain inlet data, pressure terms, and explicit bounds on every
   dropped small term. Use the actual common preheat pressure datum.
2. Represent large-Md parameters by logarithms and normalized products.
   Avoid creating Decimal precision proportional to exp(Md) solely to
   subtract huge stage coordinates or form mu. Preserve tiny terms through
   signed logarithms or explicit interval bounds; never replace them by zero.
3. Screen the complete cutoff phase interval, splitting flat endpoints
   from the b>0 strong branch. First use Z=0.5 to select parameters; then
   certify axial cells. Local Z results are not whole-axis admissibility.
4. For a candidate passing direction AND margin, include finite 1/R shear
   terms and certify the whole O.2 stage before a fresh core run.
5. Replace finite-float bump weights on the production angular correction
   route with arbitrary-precision weights or a validated normalized solve.
   Recompute waiting, correction targets, and fourteen pressure atoms using
   the same final parameter set. Nominal waiting roots remain conditional.
6. Recompute complex pressure bounds, core contraction and C3 tail degree
   for that source. Only then regenerate finite core coefficients and the
   downstream five-moment/transition chain. Do not reuse Md1.1 solution rows.
7. Continue O.3/O.4 repair, flattening, heat collar, and exact exterior with
   the same source. Then implement actual temporal n-dependent recovery
   equations, oscillatory correction, and full Cartesian residual checks.

Full admissible background, temporal recursion, oscillatory correction,
and full corrected NS residual validation remain incomplete.
