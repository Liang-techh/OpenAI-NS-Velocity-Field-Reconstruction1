# Actual two-micro comparison positivity and finite-N C0 radii

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The actual comparison direction is now proved positive on both micro charts, for every s in [0,2] and Z in [37/100,93/250]. A coarse outward enclosure is R*D_comparison in [1.563,21.314]. This uses the comparison's own complete histories; the second comparison is not identified with the core. Actual finite-N velocity-channel C0 errors and the cumulative mean-correction C0 error are below the positive radius 2^-1024 on this axial tube. Ordinary Z1/Z3, Cy_Z3, QZZ, the outer RQy and terminal closure remain unresolved.

## Comparison direction from its own source

Write L_Z=1-delta*Z^2, d=1-Z^2, g_axis=L_Z*H0/(H0^2+sigma_core^2), and use the original normalized comparison histories JH, JK and JM. The current C6 source gives

```
W  = 1-(1-delta)*Z*JM-d*JM_Z
A0 = (1-delta/2)*JH-(1-delta)*Z*JH_Z/2-d*JK_Z+(2delta-1)*Z*JK
B0 = (1-delta)*Z*JH/2+d*JK
R*D_comparison = rho/L_Z * [g_axis*B0/(2*barphi)
                            +epsilon*(-W+A0/(2*barphi))].
```

This exact identity follows by dressing the angular moments with F0, retaining F0_Z/F0=-Lambda*g_axis, and using R=epsilon*rho and epsilon*Lambda=1. The new symbolic control verifies the cancellation. The physical Pstar is absent from this direction formula; Lambda and L_Z are different quantities. Source DAG nodes bind these identities directly, alongside the unchanged actual source recipes and graph prefixes.

The core bounds are consumed from the accepted signed-core receipt. Because the prescribed comparison has alpha=1-sigma_flat(s-1) in [0,1] and the core Phi_rho is negative,

```
Phi_core(rho(s),Z) <= barphi(s,Z) <= Phi_core(4,Z).
```

This ordering does not assert equality after alpha begins turning off. The axial comparison retains its exact core inlet plus the alpha-weighted core radial derivative integral. Saved same-source U-U0 mixed bounds give barU in [1.4,2] and |barU_Z|<=5. The complete core inlet integrals and the positive Duhamel kernels then give

| own comparison history | enclosing interval |
|---|---:|
| JH, with rate 2 and source 2*barphi | [1/4,2] |
| JK, with rate 2 and source 2*barphi*barU | [7/20,4] |
| JM, with rate 1 and source barU | [7/5,2] |

Ordinary first Z derivatives are bounded from their defining differentiated integrands and the same positive kernels. In particular, the new provider uses the saved global Phi_Z and Phi_rhoZ bounds in

```
|(log barphi)_Z| <= C/Phi_min
 +2*hb_cap*rho_max*(C1/Phi_min+|Phi_rho|_max*C/Phi_min^2).
```

Here C and C1 are upper envelopes, never signed function values. No cap is differentiated. They give finite JH_Z and JK_Z bounds, retained separately in the report. The resulting epsilon error envelope has upper display `9.86539174966477362392472e-408906090034569537`, while the positive leading bracket has lower display `0.39099517774492621535805`. Directed arithmetic yields the strictly positive comparison-direction interval quoted above on both micro charts. No legacy C3-to-C5 cap transfer is used.

## Actual loop gates and phase normalization

The actual native equation is phi_y=-chi_R*(R*D_comparison)*phi/2, with y=log(R/Ra), dy=hb*ds and hb<=chi_R<=1. Thus the native phi remains positive and nonincreasing. The current generic source has E=F0*sqrt(2R)*phi/Pstar>0 and E_y=F0*sqrt(2R)*(phi/2+phi_y)/Pstar. The exact root identity is a=(E-2E_y)/E=chi_R*(R*D_comparison)>0. Neither a nor hb is rounded to zero. The original parameter logs prove 0<eta<1 and dstar>0 without materializing eta.

On the active source branch, kappa=a+b^2/a, Delta=kappa-2<eta, and define cut_value=sigma_flat(1-Delta/eta) in [0,1]. Then v=a*(1+t0^2+2q^2) is the convex combination kappa*(1-cut_value^2)+(2+2eta)*cut_value^2. The scalar helper calls cut_value "cut"; it does not bound the raw argument 1-Delta/eta by 1. Hence |b|<=kappa/2, v<=4, |B_phi|/E<2 and |A_phi|<=2. The source-owned flat branches give A=B=0 exactly.

The actual Poisson recipe has r=u/sqrt(1+u^2), h=1/sqrt(1+u^2) and 1-r^2=h^2>0. The Fourier series for (cos(psi)-r)/(1-2r*cos(psi)+r^2) has mean zero and squared mean 1/[2(1-r^2)]. Therefore the actual t has means t0 and t0^2+2q^2, proving Phi_phase(2*pi)=1 and a monotone inverse over phase [0,1]. This proof does not materialize q, u, r or the actual fast phase, and does not replace r by 1. The source-anchored values A(0)=B(0)=0 give the conservative primitive bounds |A|<=2 and |B|/E<=2 along a phase path of length at most 1. Endpoint-1 cancellation is not required for these estimates.

## Certified error radii and physical normalization

B in the following formulas denotes the current B_over_Pstar. It is already normalized by Pstar.

```
delta_E/E = expm1(A/N),        |delta_E|/E <= 4/N
delta_V   = B/N,              |delta_V|/E <= 2/N
delta_U   = Pstar*delta_V
|delta_m| <= 4*Pstar*hb*E_upper/N.
```

The last bound covers the joint first/second-micro history, starting only at the original zero correction inlet sc/2. The rate-1 kernel is at most 1 and the total native y-length is at most 2hb. The incoming first history is retained in the second chart. The actual background mean is not zeroed. The factor hb appears in the history integral; it is not inserted into a pointwise velocity error.

The exact anchored amplitude has F0=exp(-logC-Lambda*G). The original H0 root lies in [-j,0], H0 is increasing from that root through this target tube, G>=0 and the actual selected logC>0. The new bound F0<=1 is therefore sufficient for the physical E-channel envelope. The report retains log(Pstar*E_upper), exact directed log radii and the unchanged common integer N=2^(2^J). N, hb, F0 and the fast phase are not materialized.

| C0 error channel | directed upper-log display |
|---|---:|
| relative delta_E | -8.771843144176635750788e+408906090034569677 |
| physical Pstar*delta_V | -8.771843144176635750788e+408906090034569677 |
| cumulative physical delta_m | -8.91490221116304428462642e+408906090034569677 |

All five stored radii are strictly below 2^-1024. These are certified error radii, not actual signed corrections, chosen phases, or values of a cap. Zero is only the center of a numerical correction error box; the exact function graphs and the nonzero background are retained. This is a C0 reconstruction statement, not a Navier-Stokes residual test. Derivatives of these radii cannot be substituted for derivative bounds.

## Execution and next work

New production exits 0 in 3.701 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_direction_C0.json` is 16,440,058 bytes, SHA256 `45e41b113fa861934497c28cf5e2cd29a02e14874723d20fd13ab5ddcf6cef52`, with 1998 input and 355 implementation pins. New source/scale gates, exact comparison cancellation, positive-history interval arithmetic, phase normalization, finite-N inequalities, syntax/LF checks, inherited-packet checks and saved-only inspection pass. Fixed Luna/max workers performed bounded source scans, helper implementation and independent mathematical/directed-arithmetic review. Separate runtime model metadata is unavailable and is not claimed.

Every old signed contribution, all 212 forward terms, the six inner/lower signed sources, existing endpoint identities, all numerical derivative error caps and all physical completion limits remain unchanged. No accepted producer, owner, fixture, recurrence, norm computation, old quadrature or endpoint/phase seam proof was replayed. This receipt is new same-source arithmetic, not a replay or admission of older C3 target caps. The separate scaled-Poisson helper remains an unaccepted preserved draft.

NEXT bound actual normalized root ordinary Z1..3 and physical slow-y0..1, or construct a certified real slow-y tube. Use the new positive comparison direction and actual a/E gates to control normalized ratios before evaluating the implicit inverse. Then enclose Cy_Z3, joint mean-Z3/QZZ and outer RQy. The new C0 mean-correction radius does not bound any of these derivatives. Forward Rp/independent reserved-pulse functional terminal matching, global five-moment/common-N closure, exact heat/energy, stress/flat remainder, true n>=2 recursion/summation, realizable corrections and corrected Cartesian NS remain incomplete. Goal ACTIVE.
