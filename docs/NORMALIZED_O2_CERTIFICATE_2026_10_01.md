# F19: normalized O.2 parameter-family certificate

Md=40 is now a viable **conditional local axial-family candidate** for the
complete slow-turnoff interval and its zero-axial buffer. This is a proof
using the prescribed common preheat pressure, not a newly assembled
background or a regenerated core. Md=48 and 64 also pass. The unsuccessful
Md=16,24,32 lower-bound tests are noncertificates, not direction obstructions.

## Scope and exact source conditions

The positive slab is Z in [0.49,0.51]. The parameter family is

    log(Pstar) >= exp(Md)+11
    delta = min(1e-200,exp(-4log(Pstar)-30))
    mu = .001*Pstar^(-4)
    epsilon = .01*delta
    Rref >= exp(10)

The pressure must equal the complete reference-plus-outer PREHEAT backward
integral, with H replaced by 1. Exact reference cumulative moments must
already hold at the repaired inlet. All remaining outer stages must use
the Section 6.1 angular ansatz and the same parameters. The waiting length
must be nonnegative. These are explicit conditions; the new receipt does
not silently attach this parameter family to the existing Md1.1 core.

The interval q in [0,1], where q=log(y)/Md and y=log(R/Rref), covers
1<=y<=exp(Md). The subsequent eleven-unit buffer has Uz=0, but retains its
accumulated moments. It is proved separately using b=0 and the weak
kappa=2 inequality. The inverse-u-squared bound used on the cutoff phase
is not applied to this buffer.

For Md=40, rounded downward lower bounds are:

| Quantity | Lower bound |
|---|---:|
| uniform normalized Q | 0.1905404332 |
| finite tau=Ntheta-2/R | 0.1905070298 |
| finite negative dot/R | 0.3644244556 |
| finite strong bracket | 0.3439670547 |

This certifies the relaxed cone condition (3.23): its strong branch for
kappa>2 and its weak branch at b=0. The admissible-stress lift itself has
not been constructed by this module.

Receipt exact endpoint tuples are authoritative. The complete construction
still needs whole-axis admissibility, later radial stages, actual heat
restoration, recursive coefficients and oscillatory correction.

## Derivation implemented

`lei_ren_part1_paper_normalized_slow_turnoff.py` implements bounds from
paper equations (6.16)-(6.18) and (6.24)-(6.26). It cancels the angular
amplitude analytically:

    Q = L*Itheta/(F R)
    J = L*Iz/(sqrt(R/2)*u^2)
    h = -2*Uz_y/u
    h*Nz = -2*f'*Z*J/L

The positive convolution of Q'+Q=g(B) yields
Q>=min(Q1,g(B(y))), since B decreases and g is affine increasing in B.
This supplies stronger cellwise Q floors than the global forcing minimum.

The complete preheat angular profile obeys

    u_pre(s,Z) <= u(y,Z)*exp(-(s-y)/2)/(1-epsilon), s>=y>=1.

Before the collar this follows from slope<=-1/2 and the decreasing factor
removing axial dependence. On the preheat collar the bracket is at most
one, while its normalization is 1/(1-epsilon). Also
|partial_Z log(u_pre)|<=2|Z|/(1+Z^2), with zero derivative after flattening.
Integrating these bounds over the SAME entire backward integral gives an
explicit pressure envelope inflated by (1-epsilon)^(-2). No tail integral
is deleted, added, or fitted. On positive Z, P<=0 and P_Z>=0 imply Pi<=0,
so |Z+Pi| is bounded by max(Z_upper,Pi_abs-Z_lower).

The J equation supplies |J|/y<=C+(J1_abs-C)/y. On 512 closed phase cells,
the code encloses sigma' and 1/y, then bounds the finite direction and
strong bracket. The cells include their full interiors; they are not
point samples. The first and last cells use a flat-endpoint derivative
bound and symmetry. Exact b=0 endpoints use the weak branch separately.

Checks compare Q and J/y against common-source Md=4,5,6 boxes. Older
pressure boxes overlap the sharper analytic envelope but can be wider;
the receipt does not claim they are contained in that envelope. Independent
scalar derivative checks are smoke checks, not the uniform proof.

## Logarithmic parameters and continuous waiting root

`lei_ren_part1_paper_logarithmic_outer_parameters.py` constructs directed
log(Pstar), log(mu), log(delta), log(epsilon). The huge pulse length remains
the symbolic term 13/mu. Stage distances are sums of shared edges, preserving
the one-unit joins, 100-unit flatten and three-unit collar exactly.

For the continuous waiting equation it uses X in [0,2]. The incoming ODE
rate before y_rel at Z=0 is at least 1-mu-32log(2)/100>1/2, and X begins
at 5/8. Positive integral bounds enclose the restore/drop transitions and
the collar deficit. The logarithmic waiting equation then proves a unique
positive continuous root and encloses it without nominal float quadrature.
The old nominal Md=4,5,6 roots lie inside these enclosures. For Md=40 the
waiting root lies approximately in

    [941541067348080059.5093, 941541067348080062.5094].

The exact root is not yet selected or transferred into a numerical pressure
datum. Angular and heat moment repair are not completed by this root bound.

## Ordered next actions

1. Use the Md=40 logarithmic parameter backend for the fourteen pressure
   atoms. Retain the continuous waiting-root enclosure and source dependence.
   Bound every positive atom, including flattening and preheat exterior;
   restore complex axial pressure jets from this same source.
2. Introduce normalized/logarithmic core quantities before rebuilding.
   Recompute the pressure majorants and choose j,Lambda,logC consistently.
   The Md1.1 contraction and degree-124 solution cannot be relabeled or
   assumed to dominate the vastly larger Pstar. Ordinary MP exponent values
   can also become infeasible after nested amplitude exponentiation; avoid
   materializing F0 when its logarithm suffices.
3. Generate the fresh core only after its contraction and C3 tail admission
   pass. Regenerate comparison, transition, uniform five-moment inverse and
   inlet for exactly the same source. Bind their hashes to this certificate.
4. Extend the O.2 proof to the whole required axial domain, treating Z=0,
   exact b=0 endpoints, and degenerate small-Q regions analytically.
5. Implement O.3/O.4, flatten/collar, corrected heat exterior and their stress
   conditions. Then implement genuine temporal n-dependent recursion,
   oscillatory stress correction and independent full residual validation.

The goal remains active. This checkpoint establishes an admissible O.2
parameter route under explicit source conditions, not the completed target.
