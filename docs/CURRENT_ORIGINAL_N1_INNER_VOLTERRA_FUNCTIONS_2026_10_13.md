# True original n1 inner Volterra functions

The raw physical six-row n1 problem now has actual source, matrix and
Green-function bodies in the same graph as the true core and collar.
Its six ordinary C3 component functions denote the admitted complete
regular Volterra series. The zeroth and first words are stored separately,
and the fixed-point right-hand side retains every off-kinematic B0 entry
and every B1 axial-derivative entry.

This 2026-10-13 milestone installs function definitions. It does not claim
an evaluated original numerical n1 field or complete its radial cutoff,
five defects or terminal repair. The full original goal remains
**ACTIVE / INCOMPLETE**, including global heat/stress, genuine n>=2
temporal recursion, oscillatory cancellation and corrected NS validation.

## Restore and interface

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gate is
`current_original_n1_true_inner_Volterra_function_adapter_installed`.

```powershell
python scripts/run_original_n1_inner_volterra_functions.py
```

`CurrentOriginalN1InnerVolterraFunctions()` restores the checked core/collar
source and verifies the saved graph and all input hashes. `.systems` has
the true B0, B1, C0, g and N1 forcing functions for core and first collar.
`.functions` has the six raw physical C3 components; `.words[0]` and
`.words[1]` are actual finite integral bodies. `.fixed_point_rhs` exposes
the whole equation `W0+H(W)`. `.recovered` contains the n1 mean, physical
radial primitive, Q and V. `.finite_word(k)` extends actual finite words
on demand for the currently bounded indices 0..8; this API bound does
not truncate the all-word solution definition.

The saved common graph has 49,024 nodes. `.assert_graph()` preserves its
accepted prefixes and bindings. The CLI prints function handles and
scope; it does not select points from old enclosure centers. Reuse the
saved ancestors rather than constructing the old live hierarchy again.

## Physical system and complete series

```
x = sqrt(rho), rho=Lambda*R, epsilon_core=1/Lambda
W = (F1, Uz1, K1, P1, d_x_F1, d_x_Uz1)
D = (0,0,2,0,3,1)
d_x W + D*W/x = B0*W+B1*d_Z W+g
g4 = 2*epsilon_core*x*N1p
g5 = 2*epsilon_core*N1theta
g6 = 2*epsilon_core*N1z-4*epsilon_core^2*Z*rho*N1p/L
```

The angular and axial forcing use the original nested Z operators at
weights `-2-delta` and `-1-delta`. The pressure source uses the regular
rho expression `-S/2+Lambda*(2Q_rho+rho*Q_rhorho)`; it never divides
by R=0 at the axis. These are the actual leading core/collar fields,
with the same anchored F0, P0, width, scales and cumulative mean.

Remove only the three kinematic entries `(1,5)=1`, `(2,6)=1`,
`(3,6)=-1` to obtain C0. With G the diagonal regular inverse and
`J=G*Akin`, nilpotence gives `J^2=0`, and

```
W0 = (I+J)Gg
Hf = (I+J)G(C0*f+B1*d_Z f)
W = sum_{k>=0} H^k W0
W = W0+H(W)
```

The installed limit is the unique regular analytic solution. Existing
prefix/tail boxes remain enclosures, not selected F1/U1/K1/P1 functions.
The accepted raw preconditioned convergence packet allows an axial
derivative at every word, uses nested Cauchy radii, and bounds the entire
omitted sum. Its original geometric ratio is strictly below one.
No half-derivative-block simplification is used.

The normalized potential kernel Kchi and its different norm constants
are kept separate from this raw representation. They can only be used
after proving the full physical conjugation and source identity.

## Green path, axis and common-inner scope

For source x-coordinate t and endpoint X, the three source columns
4/5/6 have the six kernels

```
(t*(1-t^2/X^2)/2, t*log(X/t),
 -t*(1-t^2/X^2)/2, 1, (t/X)^3, t/X)
```

Each component uses its corresponding source row. The core path is
pulled back to the unit interval, keeping all endpoint powers. The first
collar adds the entire core segment and then its own phase integral;
the core contribution is never reset at rho=4. Its exact measure is
`dx=hb*x(v)/2 dv`, with `x(v)=2exp(hb*v/2)`.
Thin differences use `distance*exprel(-distance)` and exact logarithmic
distance, retaining positive nonzero widths instead of rounding to zero.
The log kernel has its integrable zero-endpoint continuation.

The regular axis powers are `(2,2,2,2,1,1)`, so all six values vanish
at the axis as axial functions. Core components use `0<=rho<=4`.
The selected n1 common-inner first interval is `0<=s<=sc/4`, ending
at `Rin=4/Lambda*exp(hb*sc/4)`. This is distinct from the leading
lower-pressure integration interval `sc/4<=s<=sc/2`. The same leading
collar body supplies both scopes; later modification increments vanish.

The accepted FTC jets of the anchored G/logF0/F0 are explicitly seeded.
Their higher derivatives follow the true rational gradient. The analytic
derivative engine differentiates solution limits and integral bodies,
never source bounds or numerical error envelopes.

## n1 radial recovery and evidence

```
mean1 = Uz1+K1
M1 = R*mean1
Q1 = [(1-delta)Z*Uz1-(1+delta)Z*K1-(1-Z^2)*mean1_Z]/L
V1 = R*Q1
```

Here M1 is the physical radial primitive, not the raw mean itself.
The n1 temporal recovery has the mean weight `1+delta`; the leading
background uses a different weight. These units must be retained when
the cutoff reconstructs its own M/Q/V.

The focused checker compares the raw matrix to the accepted mathematical
AST using independent symbolic fields and nonunit scales. It also checks
the nested Z forcing and physical-R viscous operator. Exact rational-field
arithmetic with rho=X^2 checks every row through Z3 without point sampling
or expanded-expression swell. A signed polynomial
fixture checks both Green paths, C3 zeroth/first words, genuine B1*dZ,
negative K1, exact core/collar seam and zero-axis data. Independent
recovery checks preserve mean/primitive/Q/V units and the n1 weight.
Actual-source checks consume the complete raw convergence admission and
require the true source-bound solution nodes and fixed-point bodies.

No original numerical solution oracle has been executed for this new
interface. Uniform original point accuracy, higher radial axis/join
contracts and complete matching remain separate open work.

## Next executable tasks

- [x] Bind the actual raw core/collar forcing and all six-row matrices.
- [x] Install the full Green path, complete regular limit and B1 feedback.
- [x] Store distinct finite W0/W1 words and the whole fixed-point RHS.
- [x] Recover physical n1 mean/M/Q/V with the correct temporal weight.
- [ ] Attach a same-source directed numerical evaluator for the complete
  functions. Reuse the existing regular-series encloser only after exact
  physical conjugation, parameter/datum and chart-coordinate guards.
  Keep complete tail and high-order Cauchy allowances; never select centers.
- [ ] Build first/second radial function jets from the actual six-row
  equation, preserving the finite even-component axis limits and each
  odd-component radial domain. Bind the required mixed axial derivatives.
- [ ] `N1-INNER-CUTOFF-FUNCTIONS`: apply the original flat cutoff to
  F1/Uz1 on first fraction q=s/(sc/4); recover M1 by its own integral,
  then K1/V1 and pressure. Keep cumulative history after the local
  cutoff becomes zero. Use the genuine five signed source integrands.
- [ ] Bind all signed inner endpoint sectors and remaining lower target
  functions; transfer original pressure/C5 bounds with nonzero errors.
- [ ] Compose full five defects, select actual linear alpha/beta repair,
  recover repaired own M/V/P and prove terminal five zero identities.
- [ ] Complete global matching, energy, heat and stress/remainder before
  true n>=2 temporal recursion, oscillatory cancellation and independent
  corrected Cartesian NS validation with fixed/restricted forcing.

The full objective remains in [PROJECT_GOAL.md](PROJECT_GOAL.md).
