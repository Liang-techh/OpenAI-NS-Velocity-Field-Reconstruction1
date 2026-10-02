# Actual reference continuation, axial restoration and five defects

The accepted actual Rsh exit now continues through the original reference
swirl branch to Rz, the original axial restoration to exact V=4Z, and the
unpatched reference branch through Rm and Rh. Callable source enclosures
retain axial derivatives through order5 for velocity, analytic pressure
datum and all five moment histories; mean-recovered radial velocity has
axial4. The actual normalized five defects at Rh are supplied through5.
These are enclosures of the original analytic functions, not recomputed
point coefficients or a connected five-bump repair.

## Source radii and inherited data

Use the original formal Rref=110*(Cstar*Pstar)^10 and Rsh=110*exp(T), with
the accepted selected T=400*Abar. The reference length is represented by
log(Rref/Rsh)=10*(logCstar+logPstar)-T. The offsets are exact source data:
Rz=e^-8*Rref, restoration exit=e^-7*Rref, Rm=e^-6*Rref, Rh=e^-5*Rref.
No enormous physical radius or exp(logCstar) is materialized. Small offsets
are never recovered by subtracting rounded huge absolute radius logs.

The angular field remains u=Pstar*(R/Rref)^.1/(1+Z^2). Before restoration
the actual axial field remains V110. Write E=V110-4Z; its source enclosure
is computed directly as j+epsilon*Psi(4,Z), plus the separately admitted
three bridge drive and three first-switch drive caps. The same actual
V110 passes unchanged through the second switch and long reshape.
The accepted C2 theorem refines only coefficients0..2 about j; higher
derivatives retain the analytic source and drive bounds.

The original analytic pressure datum P0 is inherited unchanged. Agreement
of angular or axial velocities never resets the accumulated moments.

## Center first, then enclose

The same current-amplitude-normalized shapes used by the reshape stage
satisfy

```
Mtheta=sqrt2*R^1.5*u*H, Mtheta_z=sqrt2*R^1.5*u*K,
Mz=R*m, Mztheta=R*A-R*u^2*b/2, Mp=u^2*p/2.
```

Here A is the normalized integral of V^2, separately from its swirl-energy
term. Center the shapes algebraically:

```
em=m-4Z, eh=H-5/8, ek=K-4Z*H,
aa=A-8Z*m+16Z^2, eb=b-5/6, ep=p-5.
```

With y=logR, their exact source equations are

```
em_y+em=E*alpha, eh_y+1.6*eh=0,
ek_y+1.6*ek=E*alpha, aa_y+aa=E^2*alpha^2,
eb_y+1.2*eb=0, ep_y+.2*ep=0.
```

Alpha is1 on the initial reference branch, then the original
1-sigma(t) for t=log(R/Rz) in[0,1], then0 after restoration. For constant
E before restoration, for example em=E+exp(-gap)*(em_Rsh-E).
All incoming terms, including the negative terms in these differences,
remain in the actual source expressions. A broad interval subtraction of
two reference baselines is not used to define the small defects.

The source-dependent history cap is positive exp(-1000*log10-2*logPstar).
Each incoming Taylor coefficient is combined with its exact source decay
in logarithms before it can be enclosed by this cap. Its Pstar^-2 factor
keeps negligible histories small enough for the fourth defect's units.
The cap is only an enclosure; it never replaces a source history by zero.

## Full restoration and normalized defects

The complete restoration kernels are

```
J(k,j;t)=integral_0^t exp[-k*(t-s)]*(1-sigma(s))^j ds,
(k,j)=(1,1),(1.6,1),(1,2).
```

Directed closed cells retain exact positive exponential weights. Interval
t includes the complete uncertain endpoint segment, with no midpoint
cutoff or omitted tail. At t=1 the bounds are intersected with the
previously admitted integrals of the identical original cutoff, multiplied
by exp(-k). This refines the same function and retains the incoming terms.

For x=log(R/Rm), Am=e^-.6*Pstar/(1+Z^2), the original (4.42) rows are

```
d1=exp(x)*em
d2=exp(1.6*x)*ek
d3=exp(1.6*x)*eh
d4=exp(x)*aa/Am^2-exp(1.2*x)*eb/2
d5=exp(.2*x)*ep/2.
```

In particular DeltaMztheta-8Z*DeltaMz=R*aa-R*u^2*eb/2. Both the reference
quadratic term and the swirl term are required in this identity. After
restoration all five physical defects are constant in R until a patch is
applied. Values at Rh normalize with Rm,Am, exactly as in the paper.
Zeroth-order d1,d2,d4 remain strictly positive for the actual source.

## Reproduction and evidence

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage restoreprofiles`
from experiments/root_st073 with accepted prerequisites. The complete
ordered pipeline has90 modules. The new producer and independent checker
are `lei_ren_part1_paper_compliant_reference_restore_profiles.py` and
`lei_ren_part1_paper_compliant_reference_restore_profiles_check.py`, with
separate source-bound JSON receipts.

The checker covers333 actual source transport-cap inequalities,350
velocity/pressure profile bounds,1080 centered moment bounds and30 actual
normalized defect derivative bounds. Independent finite fixtures compare
18 full restoration integrals and108 derivatives of separately integrated
physical moment primitives. Symbolic checks prove five physical primitive
RHSs, the fourth-defect reference/swirl cancellation, five normalization
factors and five unpatched radial constants. Finite fixtures do not admit
the actual source.

Exact checks preserve actual Rsh histories, Rz/restore inlet and P0, and
the exact restored4Z field. The source chain explicitly checks the core
formula, three bridge caps, three switch caps and unchanged Rsh velocity.
The original cutoff complement identity and its flat endpoints bind the
endpoint refinement to the same function. The E C2 check compares covers
directly: subtracting a second uncertain j would lose its correlation.
A read-only GPT-5.6 Luna/max audit found no material formula error and
requested these direct source/cutoff bindings.

## Remaining work

Connect the existing implicit five-bump solver to these actual defect
functions. Prove agreement with its correlated C1 source admission,
upgrade coefficient jets as required by the complete field, and transport
the actual inlet moments through the patch. Prove all five terminal
identities as functions of Z and retain the relaxed cone margin.

Actual bridge/switch/reshape/reference radial-phase mixed4 and all
two-sided high-order joins remain incomplete. Full inner/pre-O3/Cartesian
assembly, terminal-energy domains, admissible stress lift, independent
flat remainder, genuine n-dependent recursion and oscillatory correction
remain open. The original unlocalized whole-space kinetic energy remains
infinite. This stage adds no cutoff and claims no full-field, stress-cone,
global-energy or temporal-recursion completion.
