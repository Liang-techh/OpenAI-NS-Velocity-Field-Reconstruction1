# Inner exit collar stress coefficients — 2026-09-30

The shared-pressure, five-moment collar now supplies inertial stress, shear,
and total stress through width order one at Z=.3. The new replay is
`experiments/root_st073/lei_ren_part1_paper_collar_low_width_stress.py`.
Its JSON retains directed nominal and conditional pressure-error intervals
at s=0,.25,.5,.75,1,2. These are sampled stress coefficients, not a uniform
stress-cone certificate or a full nonlinear collar solution.

Write W=h_b, chi=chi0(s), J=integral chi, a=Ra F_R/Fa,
D=Itheta0/Fa, and b=sqrt(Ra/2). Since y=ya+Ws, differentiation lowers
the width order by one. The recovered shear coefficients are

g1'=-D chi/2,

g2'=-(Ra D_R s chi+D(1-chi))/2,

u1'=-b Iz0 chi,

u2'=-b[Ra Iz_R s chi+Iz0((1/2-a)s chi-D J chi/2+1-chi)],

S_theta=2Fa[g1'+W(g2'+g1 g1')]+O(W^2),

S_z=sqrt(2/Ra)[u1'+W(u2'-s u1'/2)]+O(W^2).

The inertial terms use the same recovered F, Uz, P and all five actual
moments, including their axial derivatives, in Section 3.16–3.18.
No moment target or replacement pressure is substituted.

Leading terms satisfy S0=-chi I0 and T0=(1-chi)I0. Both identities were
checked with directed intervals at all six locations. The axial inertial
stress is approximately -9.14027665195e-6. At s=.5 the leading total axial
stress is approximately -4.57013832598e-6. These are unnormalized profile
stress values, not Cartesian momentum residuals.

For s>=1, chi0=0: S0=0, but the width-tied switch contributes S1=-I0.
Both shear-component identities were checked at s=1 and 2. Consequently
the endpoint cone cannot be assessed from leading shear alone. Interior
leading stress and shear are opposite in direction when 0<chi0<1, but
this does not by itself establish the kappa and angular-width inequalities.

The cached paper text, lines 2789–2824, confirms the strict cone:
T dot S < 0, kappa > 2, and (kappa-2)(T dot S_perp)^2 < 2(T dot S)^2,
where kappa=-|S|^2/(F S_theta). For the aligned leading pair, the cross
product vanishes and kappa=chi(D^2+E^2)/D, with E=Iz0/Fa. The relaxed
kappa<=2 branch reduces to D^2+E^2>2D. Stress-free points are excluded
from strict inequalities.

Importantly, lines 5891–5913 distinguish the admissible interval ending
at R_an from the later relaxed collar. After the finite-width switch,
chi=epsilon_b and the paper states kappa<1. The later portion should
therefore not be assigned the strict admissible-cone gate. Locate R_an
on the accepted schedule and apply each regional condition separately;
the low-width expansion does not determine that location or certify
the finite-width gate by itself.

Next work must enclose finite-width stress and its appropriate regional
directional margins, including vanishing-stress edges. Width-two stress needs width-three
velocity coefficients because radial differentiation divides by W.
Original source/axis errors, pressure orders above nine, omitted width
orders, global Z coverage, and the nonlinear collar remainder remain open.
This is spatial width recursion; true n-dependent temporal recursion is
still open.
