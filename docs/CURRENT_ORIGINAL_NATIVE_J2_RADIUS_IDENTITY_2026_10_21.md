# Current original native J2 radius identity - 2026-10-21

Goal ACTIVE / INCOMPLETE. This stage proves the two actual J2 endpoint radius
functions agree on the inherited owned axial cell Z in [0,1/2], through ordinary
Z orders 0 and 1. It does not prove the J2 B jump is zero.

## Exact source identity

Write D=1000+4*(11+exp(40)) and A=4*exp(-D)>0. The frozen source definitions are

- frozen_macro: R(x)=A*exp(2*hbB+x*(log(100)-log(A)-2*hbB)); native coordinate node 70265, radius node 71448.
- first_switch: R(s)=100*exp(hbS*s); native coordinate node 73190, radius node 85905.

At the upper macro endpoint x=1, the identical node 66985 and its negative
66986 cancel. The positive real log domain permits
A*exp(log(100)-log(A))=A*(100/A)=100. At the lower switch endpoint s=0,
100*exp(hbS*0)=100. Therefore both endpoint functions equal the exact rational
100, and both ordinary Z derivatives equal zero. This is endpoint equality;
it is not an equality of the full two radial chart functions.

The source parameters hbB and hbS retain distinct node identities 8 and 9.
The original MacroFlow finite-log guard, the exact width recipe
hb=epsilon_b=cstar*K^-100, and the saved finite mpf log-enclosure tuples justify
the finite real parameter domain. The enclosure is used only to establish
that domain. No endpoint, midpoint, width cap, or finite-N error bound is
selected as an actual source function value, and h is not evaluated.

## Actual execution and retained scope

The frozen consumer was entered once and passed in 0.47408070001984015 seconds.
All 57 source checks passed; the saved proof graph has 10 steps. The source,
entered, and full result checkpoints are preserved; the result SHA256 is
5b7fb2939d3783d5f65ea26bb10e044cd7b2d6fb994ecdbfa67701045b9db520.
The saved-only independent review passed. No ancestor operator, DAG walker,
endpoint canonicalizer, quadrature, recurrence, or numerical exp/log evaluation
was called. Runtime model metadata is unavailable and recorded as null.

The R identity set is now J2,J3,J4,J6,J8,J9,J10 (seven). Inherited B zeros
J0=520 and J1=521 remain accepted. Every remaining B identity at nodes 522..535
is still open, including J2, whose U and physical Q joins are not proved here.
The native identity remains
I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy).
The parent F-T reduction is unchanged:
F-T=1/2*(I_N_Rc-I_O_active+Rc^2*integral_1^2 x*DeltaS dx).

The ledger remains cumulative 212, raw signed 213, reduced 183, with all six
signed channels, every old error, one reserved contribution, 17 K integrals,
and 2 outer-B ports. The 109 parent checkpoint hashes and five failed-v1
artifact pointers are retained; this leaf adds three checkpoints (112 total).
The publication authenticates 1834 source pins. Numerical completion, global
five-moment closure, physical construction completion, temporal recursion,
and the full NS residual goal remain unclaimed.

## Next executable task

Prove the actual J2 U C1 source join 140685/140686 versus 148085/148086.
Use common Pstar33572, inverse Pstar69245 and inverse N99411; discharge the
endpoint transport from frozen_macro nodes70152/70153,71706/71707 to
first_switch nodes72994/72995,87147/87148, respecting native s=0.
Then prove physical Q C1 (140780/140784 versus148180/148184), including
positive L25272, the Q0 recurrence, and m140749/148149. The bounded static
extraction is embedded in the new report as next_UQ_static_audit and is not
an accepted identity. Do not rerun this radius consumer or the broad seam pass.

Report: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J2_radius_identity.json

Full saved result: experiments/root_st073/original_native_J2_radius_identity_execution.json
