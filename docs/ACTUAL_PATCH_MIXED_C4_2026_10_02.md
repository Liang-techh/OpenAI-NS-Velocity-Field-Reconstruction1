# Actual five-bump patch mixed spatial derivatives

The accepted actual five-moment patch now supplies all mixed derivatives
in x=R/Rm and Z, and y=logR and Z, through total order4. All four velocity
and pressure rows retain their physical radial and axial prefactors. The
five true physical primitive rows have x/Z mixed4 enclosures. The original
compact bump edges and the patch's Rm/Rh joins to reference branches have
functional smoothness proofs. This is local patch scope; other inner
charts/interfaces and the full Cartesian field remain unfinished.

## Source and units

The implicit coefficient functions are the SAME accepted axial5 family;
the current cumulative primitives and original pressure datum are copied
from `CompliantActualMomentPatch`. No old accepted module is rewritten.

V=Uz=4Z+g is the raw source axial velocity. It is NOT divided by Pstar.
Only the angular velocity is normalized by Pstar. With
am=Am/Pstar=exp(-.6)/(1+Z^2), H=x^.1+f,

```
Utheta/Pstar=am*H, Uz=V,
Ur/sqrt(Rm/2)=sqrt(x)*Q,
P/Pstar^2=P0/Pstar^2+am^2*S.
```

In particular the centered-energy RHS is g^2/Am^2-H^2/2, retaining
Am^-2=(1+Z^2)^2*exp(1.2-2logPstar). Replacing Am^-2 with am^-2 would lose
Pstar^-2 and change the actual field. Physical Rm and Pstar remain formal
positive common prefactors; the exact R derivative conversion is
partial_R^k=Rm^-k*partial_x^k.

## Differentiate original shapes and primitive equations

The compact bump is gamma_j(x)=beta((x-c_j)/r)/(r*N), with original
r=1/40, centers5/4,3/2,7/4 and the accepted normalization N. The generic
unscaled beta derivative tool is used with precisely these parameters,
not the outer angular pulse's different width. Ordinary derivatives are

```
gamma_j^(k)(x)=beta^(k)((x-c_j)/r)/(r^(k+1)*N).
```

Z derivatives act on the accepted axial5 coefficient family; x derivatives
act on these original bumps and x^.1. No capped numerical shape is
differentiated. The accepted shape0/shape1 covers refine the identical
functions, while higher derivatives use the original flat beta bounds.

Let M=Mz/Rm, T=Mtheta/(sqrt2*Rm^1.5*Am),
J=Mtheta_z/(sqrt2*Rm^1.5*Am),
E=(Mztheta-8ZMz+16Z^2R)/(Rm*Am^2), S=Mp/Am^2. Their exact equations are

```
M_x=V, T_x=sqrt(x)*H, J_x=sqrt(x)*H*V,
E_x=g^2/Am^2-H^2/2, S_x=H^2/(2x).
```

Repeated ordinary product rules provide x derivatives through4. The
zeroth primitive rows retain the actual current partial integrals. Their
true physical axial derivatives include all Am(Z) factors. With m=M/x,

```
Q=[2ZV-(1-delta)Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2).
```

The final sqrt(x) in Ur/sqrt(Rm/2) is differentiated before forming the
mixed grid. This is essential even though the parent packet stores only
Ur/sqrt(R/2)=Q. Axial5 primitive data supplies the extra Z derivative needed
for recovered Ur through mixed total4.

Log-radius derivatives use the exact Stirling transform
partial_y^k=sum_j S(k,j)*x^j*partial_x^j. Every reported mixed grid entry is
an ordinary derivative; stored axial Taylor jets still divide by n!.

## Functional joins

The six original support edges are49/40,51/40,59/40,61/40,69/40,71/40.
All beta derivatives through4 vanish at both endpoints. Positive
normalization and affine scaling preserve these zero limits. Actual
primitive constants remain continuous across every edge by their source
integrals.

Rm is strictly below all supports, so its field agrees with the incoming
unpatched reference branch on an open neighborhood, including all mixed
derivatives and inherited moment defects. Beyond71/40 every partial bump
integral is its full weight. The SAME accepted implicit equations make
the five terminal defect functions and their derivatives identically
zero. Rh=e*Rm is strictly beyond that endpoint and therefore matches the
exact reference field and all mixed4 derivatives. These are source
function identities, not equality inferred from overlapping intervals.

## Reproduction and remaining work

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage patchmixed`
from experiments/root_st073 with accepted prerequisites. The ordered
pipeline has94 modules. New producer/checker are
`lei_ren_part1_paper_compliant_actual_patch_mixed_C4.py` and
`lei_ren_part1_paper_compliant_actual_patch_mixed_C4_check.py`, each with a
source-bound JSON receipt.

Actual checks cover720 x/Z velocity-pressure bounds,720 y/Z bounds and900
physical primitive bounds over the entire patch/axial domain, original
support edges and interior packets. Independent closed physical primitive
integrals with nonzero histories compare135 x/Z derivatives and60 y/Z
derivatives. Five original normalized gamma derivatives are compared with
independent quadrature normalization and direct differentiation. Exact
checks prove4 coordinate identities,6 support edges and10 flat endpoint
limits. Finite fixtures do not admit the actual source.
A read-only GPT-5.6 Luna/max review found no material formula error and
confirmed the source units, primitive recovery, radial prefactor, Stirling
transform and conditional functional joins. Primitive grids currently use
x/Z; downstream callers needing explicit y/Z primitive grids can apply the
same exact Stirling transform to their retained derivative rows.

Remaining: actual bridge, both microscopic switches, long reshape and
reference/axial restoration mixed4; all their high-order core/annular
interfaces; missing inner/pre-O3 dispatcher and complete physical
Cartesian spatial/time assembly. Terminal-energy domains, admissible
stress lift, independent flat remainder, true n-dependent recursion and
oscillatory correction remain open. Original unlocalized whole-space
energy remains infinite. Full-field/stress/global-energy/temporal gates
remainfalse.
