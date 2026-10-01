# Compliant epsilon source and pressure/core transfer (F36)

The `c_epsilon=.001` source is now separate from the legacy `.01` source.
Its continuous waiting enclosure and fourteen pressure-atom envelopes are
regenerated. A correlated pressure perturbation bound proves existence and
closeness of its analytic core and supplies a callable enclosure of the
completed local 144-order coefficients. This avoids a full recurrence rerun
for this parameter change. It does not transfer the connecting cone or
finish the corrected outer field.

## Source and pressure comparison

`lei_ren_part1_paper_compliant_pressure_source.py` builds the new data without
changing the legacy parameter/datum modules or their hashes. The inherited
parameter constructor initializes common upstream data; the new class then
explicitly replaces epsilon, recomputes the continuous waiting equation and
updates the waiting edge. The source hash is independent of interval precision.

New implicit source:
`5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`.

Old implicit source:
`082d18b8f1ec2df84526fbdc190217d97dfb8825e823fa4c9e017c8a58d001bb`.

These sources define the **raw preheat profile, with H replaced by 1**. The
upstream swirl schedule is independent of epsilon and waiting. Its axial
flatten finishes at `y_f`, before the common start of waiting `y_t`. All
changes between these two preheat profiles occur after `y_t` and are then
Z-independent. Consequently `Delta P0` is an entire constant in Z and all
its positive-order axial derivatives are exactly zero. This correlation is
between the true sources, not between independently enlarged atom boxes.

The comparison fixes the same reference radius, `F0`, delta, j, Lambda and
Cstar. Pressure itself is independent of the reference radius because
`dR/R=dy` and the normalized preheat swirl is specified in y and Z.
Here `R_t` means the **start** of waiting; `R_tail=R_t*exp(tau)` is its end.
The normalization has the equivalent forms

`c_inf=Utheta(R_t)*R_t^((1+delta)/2)/(1-epsilon)`

and

`c_inf=Utheta(R_tail)*R_tail^((1+delta)/2)/(1-epsilon)`.

Thus epsilon changes the terminal normalization but does not rescale the
upstream profile. Waiting changes the collar location, not this normalization.

A coarse post-`y_d` pressure bound is too weak for core transfer after the
physical `Pstar^2` conversion. Instead the code proves
`y_t>=13/mu>cut`, with `cut=10*logPstar+1000`, and integrates the true
post-`y_t` envelope before using monotonicity:

`T_i <= exp(.6-cut)/(2*(1-epsilon_i)^2)`.

It follows that `|Delta P0/Pstar^2|<=T_old+T_new`; multiplying by `Pstar^2`
gives the physical and complex constant-source bound. Both caps remain
strictly positive. The code does not expand `exp(-13/mu)`, set that tail to
zero, or claim a pointwise profile approximation below `y_t`.

## Analytic core transfer

`lei_ren_part1_paper_compliant_core_transfer.py` rebuilds all twenty terms
of paper (8.50) on common old/new component envelopes. The new fixed-point
map has size and Lipschitz bounds below one half. `P0` is fixed within each
solution, while the transfer compares the two different data.

A constant pressure change still changes the core. From (8.16), (8.38):

`Delta g=2*(1+delta)*Z*Delta P0`,

`Delta Psi0=-(r/(2L))*Delta g`.

All remaining P0 dependence has been removed from the nonlinear operator;
the restored pressure contribution is `F0^2*V(Phi^2)`. With the same F0 and
other fixed data, the two operators therefore coincide. The common envelopes
contain the segment joining the two fixed-point balls. For their common
scaled Lipschitz upper q, the receipt proves

`||X_new-X_old||_Xh <= ||Delta X0||_Xh/(1-q) = D`.

It records physical `Uz` and pressure prefactors separately. Mixed C3
differences are bounded on `0<=r=Lambda*R<=4.1`, `-1<=Z<=1`. These are
normalized/scaled derivative bounds, not full physical Cartesian residuals.
The independent checker verifies six pressure, normalization and resolvent
identities, the physical pressure scale and the complete contraction gates.

## Finite coefficient enclosure view

`lei_ren_part1_paper_compliant_finite_core.py` binds the unchanged old state
to the new analytic sensitivity receipt. It exposes each retained coefficient
as the old interval plus a symmetric sensitivity cap. The Xh embedding is

`cap_nm = D_field*binomial(n+m,m)/(20^n*h^m*(n+1)^2*(m+1)^2)`.

The local center family is `[.49,.51]`; there are 33,060 retained A, Uz and P
coefficients through radial order 144. A is normalized Phi. Uz includes its
physical epsilon factor. The stored P rows are epsilon times physical P,
and their caps include that additional factor. The axis A and Uz rows and
positive-order pressure derivatives remain unchanged; only the axis pressure
value receives its constant-source cap. Positive caps are upper-bound data,
not proof that the actual coefficient differences are nonzero.

The old 20 MB state is preserved. This compact callable view provides new
coefficient enclosures; it does not recompute their point values or relabel
the old solution. The checker verifies 443 exact unchanged axis coefficients,
first axial recovery sensitivity, and early/middle/final coefficient samples.

## Next work and limits

1. Transfer the inner-exit, physical norm, Cstar/K1 and reference-join bounds
   using the new core differences. Regenerate source-bound receipts, not labels.
2. Re-establish the functional five-defect and five-bump correction for the
   compliant family; retain the same pressure correlation in all adapters.
3. Recompute selected-family outer waiting, actual angular/pressure repairs
   and the exact heat targets under `.001`.
4. Select the axial energy parameter from the complete future corrected tail;
   prove whole-outer cone margins and then the admissible lift/flat remainder.
5. Implement true order-dependent temporal recursion, oscillatory corrections
   and independent Cartesian residual validation.

The numeric epsilon threshold is satisfied. All Section 7 hypotheses, the
matched heat field, full background admission and temporal recursion remain
unproved. An actual H_delta pressure difference is generally Z-dependent;
it is excluded from this constant-source proof. The later angular/pressure
corrections must restore the prescribed preheat P0 rather than quietly
changing it to the uncorrected heat profile's pressure.

## Reproduction

Run the five Python files with the common prefix
`experiments/root_st073/lei_ren_part1_paper_` in this order:

1. `compliant_pressure_source.py`
2. `compliant_core_transfer.py`
3. `compliant_core_transfer_check.py`
4. `compliant_finite_core.py`
5. `compliant_finite_core_check.py`

Their adjacent JSON receipts contain the transitive input hashes, scope
flags, bounds and example enclosures. All five completed successfully for F36.
