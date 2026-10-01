# Inner exit collar spatial width-coefficient bounds — 2026-09-30

The first and second width coefficients now have directed enclosures over every s in [0,2] in the Section9.25 inner exit collar. This is not the exterior heat collar. The axial coordinate remains Z=.3; the common conditional pressure-integral perturbation and all three available Z slots are retained.

The domain is covered by64 contiguous cells, with interval arguments for the whole cell. No sampled interpolation is used to infer the intervening values. Both the eight normalized states and all five unnormalized moment coefficients are enclosed in each cell. The finite pressure parameter is aggregated at1 only for diagnostic bounds; pressure powers above9 and width orders above2 remain excluded.

## Spatial primitives and formulas

Let chi0(s)=sigma(1-s) for s<=1, zero afterwards. Define J(s)=integral_0^s chi0, K(s)=integral_0^s J, and M(s)=sJ-K. A directed grid for integral sigma and integral x sigma gives, for s<=1:

J=s-integral_0^s sigma(x)dx;

K=s^2/2-s integral_0^s sigma(x)dx+integral_0^s x sigma(x)dx.

For s>=1, J=1/2 and K=integral_0^1 x sigma(x)dx+(s-1)/2. Since J and K are nondecreasing, their endpoint enclosures bound every interior point of each cell. Taylor-remainder integration covers interior primitive panels; positive flat-end bounds cover the singular-formula endpoints. Independent resolved quadrature checks J and K at s=.25,.5,.75,1,2.

Write D=D0, B=sqrt(Ra/2)Iz0, a=Ra F_R/F, and u=Uz0 from the same inlet. The first-width coefficients divided by h_b are -DJ/2 for g, -BJ for u, and the original six integral slopes multiplied by s. The dynamic second-width coefficients divided by h_b^2 are:

g2=-(Ra D_R M+D(s-J))/2;

u2=-sqrt(Ra/2)[Ra Iz_R M+Iz((1/2-a)M-DJ^2/4+s-J)].

The remaining coefficients are theta2=2s^2-DK, mz2=u s^2/2-BK, mixed2=2u s^2-uDK-2BK, axial2=u^2 s^2/2-2uBK, swirl2=s^2-DK, and p2=s^2/2-DK. At s=2 these reduce to the previously checked endpoint formulas. The leading switch's support lies in the direct-core comparison region, which is why the same inlet radial derivatives suffice.

## Whole-domain conditional pressure contribution

Maximum retained coefficient error bounds over the64 cells are:

| Coefficient / h_b^order | Value-slot bound |
| --- | ---: |
| u, first width | 2.25365546804558e-36 |
| u, second width | 8.17088079100805e-36 |
| g, first width | 1.54385303075742e-35 |
| g, second width | 6.71245613718127e-35 |

The second-width u bound is slightly wider than its endpoint bound because the entire spatial cell is enclosed and common primitive dependencies are relaxed into an interval box. Exact intervals, pressure contributions, Z derivative slots, and five-moment coefficients are recorded per cell.

## Replay and remaining scope

Run `python experiments/root_st073/lei_ren_part1_paper_collar_low_width_spatial.py`. The exact inlet reader restores the committed pressure-perturbation receipt, avoiding another finite-core regeneration. The spatial calculation evaluates analytic coefficients and directed primitives; it does not rerun the nonlinear spatial collar chain.

The profile radius expansion adapter now accepts general s, preserving the existing s=2 endpoint behavior. Its six independent mixed radial-flux derivative checks still pass.

These bounds concern retained coefficient functions, not the complete nonlinear collar solution. Original axis/adapter errors, omitted pressure/width orders, full spatial divergence, global axial coverage, stress-cone admissibility, relative flatness, exterior heat matching, and true temporal recursion remain open.
