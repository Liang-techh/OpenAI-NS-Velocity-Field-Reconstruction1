# Actual five-moment patch and implicit axial5 coefficient family

The actual core-to-Rsh-to-Rm history now feeds the original five-bump
repair. Its five coefficient functions have axial5 enclosures from one
smooth implicit solution, rather than independently fitted point values.
The patch supplies velocity, pressure and all five cumulative primitives
through axial5; mean-recovered radial velocity has axial4. The five
terminal reference identities are connected to this actual source as
functions of Z. Full mixed spatial derivatives and joins remain open.

## Bind the actual source before solving

Let G=log(Rz/Rsh), alpha=exp(-G-2)=Rsh/Rm and E=V110-4Z. At the actual Rsh
exit retain the centered shapes em,eh,ek,aa,eb,ep from the accepted reference
provider. The exact normalized source rows at Rh, in the paper's Rm,Am
units, are

```
d1=k1*E + alpha*(em-E)
d2=k2*E + alpha^(8/5)*(ek-(5/8)*E)
d3=alpha^(8/5)*eh
d4=Am^-2*[k4*E^2+alpha*(aa-E^2)] - .5*alpha^(6/5)*eb
d5=.5*alpha^(1/5)*ep.
```

The k1,k2,k4 are the SAME original full restoration integrals from the
accepted functional defect admission. The negative incoming E and E^2
terms and negative inherited swirl term in row4 remain explicit. Exact
transport algebra independently proves that these functions equal the
previous provider's actual five defects. Intersecting their two enclosure
representations sharpens these same functions; interval overlap itself
is not the identity proof.

Every incoming Taylor coefficient is combined with its original decay
log before applying the positive 10^-1000/Pstar^2 enclosure cap. The
transported tails' summed value/first-derivative bounds are each directly
checked against the older 10^-800/Pstar^2 C1 cap. Large dependent Rsh boxes
therefore do not need to satisfy a tiny norm themselves. The actual
source C2 bound on E-j is within the old admitted correlated C1 bound.
These facts bind the original complete source norm and anisotropic
coefficient scales to the actual functions. They do not assign that
Banach norm to arbitrary independently chosen Taylor-box functions.

## One implicit family through order5

Use the original normalized compact bumps, radius1/40 and centers5/4,
3/2,7/4. The coefficient order is h=(c1,c2,xi1,xi2,xi3). The same source
map is

```
L*h + Q(h,Am^-2) + d = 0.
```

The original strict self-map and uniform contraction certify h0 and its
C1 family with the unchanged anisotropic axial/angular boxes. For each
Taylor degree n=2..5, form Q using all known coefficients below n and
temporarily h_n=0. Its coefficient contains every lower-order quadratic
convolution and explicit Am^-2 derivative. The missing term is exactly
the same pointwise Jacobian J(h0)*h_n:

```
J(h0)*h_n = -d_n - Q(h_<n,Am^-2)_n.
```

A fixed preconditioner and a weighted interval inverse enclose this
linear solve. Midpoints choose only the preconditioner; actual defects,
weights, tails and amplitudes are not projected to midpoints. The actual
higher-jet inverse contraction is at most approximately1.138e-7.
The coefficient data remains enclosures of the unique implicit functions,
not newly recomputed point coefficients. This is the leading n=0 spatial
repair, not the future time/coefficient recursion.

## Partial fields and functional terminal closure

With x=R/Rm, f=sum xi_i*gamma_i and g=c1*gamma_1+c2*gamma_3,

```
u=Am*(x^.1+f), V=4Z+g.
```

All normalized moment changes use the original directed partial bump
integrals. Actual initial defects are added to those changes; P=P0+Mp
uses the original analytic datum and Ur is recovered from the SAME mean
primitive. The actual Rm inlet is bound by source transport and checked
for interval consistency in every retained axial derivative.

After x>=71/40 all bump integrals equal their full weights. The remaining
five defects are exactly Lh+Q(h)+d, hence identically zero for the unique
implicit family, including its derivatives through5. A terminal zero
enclosure is an analytic identity refinement after this proof. It does
not reset actual histories, choose midpoint coefficients or fit pressure.
The unreduced interval residuals are retained as diagnostics only.

The inherited C1 smallness theorem also supplies the LOCAL relaxed cone
bound for this same actual patch: margin>=approximately3.7999282,
|bw|<=approximately8.9743e-6 and kappa<1. This is a parameter/source-bound
local relaxed-cone theorem; it is not a complete divergence-form stress
lift, an independent flat remainder or a whole-background cone certificate.

## Reproduce and inspect

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage actualpatch`
from experiments/root_st073 with accepted prerequisites. The complete
ordered pipeline has92 modules. New source and checker are
`lei_ren_part1_paper_compliant_actual_moment_patch.py` and
`lei_ren_part1_paper_compliant_actual_moment_patch_check.py`, with separate
source-bound JSON receipts. No old accepted module is rewritten.
Ordinary callers require the current independent checker receipt by
default. Only producer/checker bootstrap bypasses that self-receipt gate.
Stale hashes or missing source/transport/implicit-jet acceptance prevent
downstream use. A read-only GPT-5.6 Luna/max audit verified the formulas and
requested these explicit cutoff and independent-check bindings.

Actual checks cover5 transported C1 tail bounds,36 factored source decay
caps,30 coefficient Taylor bounds and497 velocity/pressure/moment bounds
over the whole axial domain and complete patch x interval, plus selected
interior and terminal packets. Exact datum and implicit-map closure
checks preserve the actual source and limited scope.

Independent symbolic checks prove5 exact dominant/tail transport identities
and20 higher-order common-Jacobian convolution identities. A separate
finite fixture integrates21 physical bump weights and recovers all30
Taylor coefficients of five nonconstant manufactured smooth functions,
with a variable inverse amplitude. The fixture inverse contraction is
approximately.3777. It checks the nonlinear derivative algorithm and does
not admit the actual source. Its initial larger test box was outside the
map's contraction range; the actual paper parameter boxes were unchanged.

## Remaining dependencies

Restore bridge/switch/reshape/reference/patch radial-phase mixed4 and all
two-sided core/annular interfaces, using original flat support identities
and formal inverse-hb factors. Complete the inner/pre-O3 dispatcher and
bind its actual Rh reference terminal data to the accepted outer chart.
Then assemble whole physical Cartesian derivatives and time derivatives.

Terminal-energy domains, admissible divergence-form stress lift,
independently bounded flat remainder, actual n-dependent recursion and
oscillatory correction remain unfinished. The original unlocalized source
still has infinite whole-space kinetic energy; no new cutoff is inserted.
Full-field/stress-lift/global-energy/temporal gates remainfalse.
