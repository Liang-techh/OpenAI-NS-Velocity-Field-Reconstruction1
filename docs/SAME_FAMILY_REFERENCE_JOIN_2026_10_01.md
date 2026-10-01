# F29: same-family long reshape and axial restoration

The F27/F28 selected analytic family now has a Section 9.5 log-coordinate
prescription and analytic relaxed-cone bounds from R110 to
Rh=exp(-5)Rref over the whole real axial interval [-1,1]. The R110 amplitude
budget is recomputed from the actual integrated exit rather than assumed.
All five moments continue from the axis and P=P0+Mp remains unchanged.
Terminal moment correction, the selected-radius corrected outer/heat field,
admissible stress lifting and temporal recursion are still unfinished.

## Reproduce and callable interfaces

```powershell
python experiments/root_st073/lei_ren_part1_paper_shared_reference_join_bounds.py
python experiments/root_st073/lei_ren_part1_paper_shared_reference_join_check.py
```

The first producer binds F28, F27, the same analytic preheat pressure,
core/seed/normalization receipts, and the existing cutoff. It supplies:

- `shape_log_jet`: log(u/Pstar) and the full angular-shear Z jet, given the
  admitted inlet B=log(Cstar*u1*(1+Z^2)) jet.
- `reference_log_jet`: reference swirl at an exact log(R/Rref) offset.
- `restore_axial_jet`: the (9.38) velocity and its log-radial derivative,
  given the admitted actual inlet v1 jet.

These are callable interval log interfaces. Their whole-axis examples
enclose the admitted inlet family; they are not point solutions or finite
evaluations of the complete physical field and its moments.

## Normalization and R110 amplitude

G is the implicit real integral of LH0/(H0^2+sigma^2), anchored at the
unique H0 root in [-j,0]. The executable receipt binds the accepted seed
branch and records directed L positivity, the two exterior sign regions,
root endpoint signs, and H0'>4 on [-j,0]. The identity

```
H0=(4Z+j)(1-Z^2)+(1-delta)Z/2
```

shows that H0<0 for Z<=-j and H0>0 for Z>=0. Thus G>=0 and the uniform
Phi<=2 bound gives Cstar*Fcore<=2. This family does not use a pole-log
primitive or choose a new integration constant.

The fixed Abar upper bound from F27 includes 10 plus the unweighted core
log and velocity norms. The actual integrated relative log, both short
switches, constant-power interval, sqrt220, and the summed C2 bound
log2+1+2 for log(1+Z^2) add less than 10. Therefore

```
||B||C2 <= Abar-10+extra < Abar <=2Abar.
```

The switches use the integrated bound
h*(.8+Cq*K^7*h), evaluated conservatively at K>=10^6. No h^-1 cutoff
derivative appears. This explicitly closes the amplitude item left implicit
by the F28 exit receipt.

## Conservative parameter choice and exact topology

The paper uses T=400 times the exact data size A. This companion chooses
T=400Abar, with Abar fixed once as a proved upper bound for A. This is a
conservative construction choice, not a claim to compute the exact A norm.
F27 already imposed the stronger radius restriction using this same Abar,
so the core, pressure, Cstar and Rref are not changed.

Keep Rref=110(Cstar*Pstar)^10 as its exact defining expression. Relative
interfaces are exact offsets from logRref:

| Interface | Offset |
|---|---:|
| Rz, start axial restoration | -8 |
| eRz, end axial restoration | -7 |
| Rm, start later moment repair | -6 |
| Rh, reference matching interface | -5 |

The stored radius margin gives Rsh=110exp(T)<Rz. The repair interval
[Rm,2Rm] lies strictly after restoration and before Rh. Absolute logarithms
are too large to retain these unit offsets in ordinary finite-precision
subtraction. The API uses the offsets directly; it never exponentiates
Abar, logCstar, or logRref.

## Angular cone and pressure

The unchanged cutoff derivative bound8 and ||B||C2<=2Abar imply
|a-.8|<=.08, safely a in [.7,.9]. Hence logu_y>=.05. The angular log
derivative interpolates the actual entry derivative and -2Z/(1+Z^2).
For V=4Z+error, the leading contribution of Hv*(-2Z/(1+Z^2)) is
nonpositive; the small velocity error bounds the remaining term.

The continuously inherited Mz/R is a positive Z-independent average of
the initial mean and subsequent velocities. Its error and the velocity
error stay below epsilon0, including the restoration, so |W-Wref|<=10epsilon0.
The direct (9.32) bounds give:

| Gate | Result |
|---|---:|
| Angular source SQ | >=1.6114999994723, required >1.4 |
| Damping 2-a/2 | in [31/20,33/20] |
| Q after one log-radial unit | >1/2 |
| D=3 crossing derivative | positive |
| Axial restoration \|bw\| | <=0.0012799288928393, required <0.01 |
| Axial restoration kappa | <1 |
| D(a-bw)-2a | positive |

For the early pressure, F<=4/Cstar and |logF_Z|<=Abar+1 give

```
||Mp(110)||C1 <=1760(3+2Abar)Cstar^-2
                 <=exp(4Abar)Cstar^-2 <=1.
```

The producer checks this in logarithms. On the shaping interval the C1
pressure increment is <=(5+30Abar)u_sh^2. The radius compatibility gives
u_sh^2/Pstar^2<=exp(-2)/(1+Abar)^2. The remaining reference pressure
increment has C1 bound7.5Pstar^2. Together with the unchanged P0 bound,
the total C1 pressure is <(Kp+100)Pstar^2.

The N=LIz/sqrt(R/2) normalization retains the actual inherited stress.
Each actual moment and P has C1 norm<=2K, and a conservative (9.13)
coefficient gives |N110|/K<20. The source of N_y+N is bounded by
(500+4Kp)Pstar^2. The admitted Rz restriction suppresses the inherited
110N110/R term without resetting it. Thus |N|<=KN Pstar^2 and
|J=N/u^2|<=20KN on [Rz,Rh]. This cancels the large amplitude in the
restoration product bw=2J*Vy/Q. The resulting relaxed cone covers both
the reshape and the restoration; full admissibility is not asserted.

## Independent checks and remaining work

The check producer differentiates the actual five primitive definitions and
the inertial stress formula (9.13), independently of the scalar bounds
producer. Both (9.32) and (9.35) reduce identically to zero, allowing arbitrary
Vy. The same check verifies the H0 factorization and the real-integral G
branch derivative/anchor identities. Seven callable log-interface checks
pass, including preservation of a unit relative offset with an unresolved
huge absolute base. The global cone proof uses analytic bounds, not these
fixtures. A read-only Luna/max review checked the dependency assumptions and
identified the normalization, full shear jet, and inherited-stress details
addressed in this implementation.
The final narrow review found no remaining mathematical blocker within this
implicit prescription/analytic relaxed-cone scope.

The next critical task is the actual normalized five-defect ledger and
Section 10 contraction. Use the whole-axis functions and inherited moments,
not a point fit or the historical hardcoded long-reshape family. Then rebuild
the corrected outer at this Rref, with its own corrections and exact heat
tail but the same preheat axis datum. Only after five functional identities
hold can pressure/radial velocity/stress and exterior matching be claimed.
Section 11 lifting, flat remainder and the n-dependent coefficient recursion
follow that matched leading profile; 144 spatial core orders are not these
time coefficient orders.
