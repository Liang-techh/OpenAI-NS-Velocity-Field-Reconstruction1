# Actual micro b/t0/p2 mixed bounds from the original background

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The actual
background b, t0 and generic p2 now have source-bound ordinary Z0..3 and
physical y0..1 logarithmic absolute norm caps on both micro charts,
s in [0,2], Z in [37/100,93/250]. The root q, finite-N derivative errors,
Cy_Z3, mean-Z3/QZZ and terminal/global closure remain open. These caps are
not signed function values or finite-N correction radii.

## Same-source pressure and amplitude

The normalized fourteen-stage pressure is p0=P0/Pstar^2:

```
p0(Z) = -M2*(1+Z^2)^(-2)-M0-F_flat(Z).
```

The saved M2/M0 intervals and the original nonzero flattening source are
retained. F_flat at order zero keeps its saved mass enclosure; higher
ordinary derivatives use k!*flatten_complex_mass_upper/rho_strip^k.
The new rational jets are evaluated over the entire axial interval, through
Z5. They are not copied from the saved narrow .371 pressure packet.
The following are rounded endpoint displays; exact mpf tuples in the report
are authoritative.

| ordinary Z order | p0 lower display | p0 upper display |
|---|---:|---:|
| 0 | -2.56442245601233104321723543 | -2.557740839917303945376079057 |
| 1 | 3.325289571074092607728672403 | 3.356373132682160781341583535 |
| 2 | 2.409371141845996946311475762 | 2.526458623607183970332789419 |
| 3 | -36.12073841855728018231057621 | -35.49028618261144223176139932 |
| 4 | 75.64384378116791071195798368 | 79.79958142209460769518388429 |
| 5 | 562.0382165012228367537099552 | 592.1535517738009301744138239 |

The physical scale is log(Pstar)=exp(40)+11. Physical P0=Pstar^2*p0, and that
Pstar^2 is restored explicitly in the reduced drives below. The exp argument
node 5 in the source DAG denotes its parameter handle.

For the same selected logC and exact anchored gauge,

```
F0=exp(-logC-Lambda*G),  0<=G<=Gbar,
-logC-Lambda*Gbar <= log(F0) <= -logC.
```

Gbar is the saved modulus on the common analytic tube; no norm producer is
rerun. This supplies a whole-tube positive F0 lower bound in logarithmic
form. F0 and its inverse are never materialized. Ordinary Z derivatives of
F0^2 through order four are proved to have the coarse absolute cap one from
the exact exponential derivative ratios and this amplitude envelope. This
cap encloses those derivatives; it does not define F0 or its values.

## Exact reductions before interval boxing

Set d=1-Z^2, L=1-delta*Z^2, rho=R/epsilon, Dbar=R*D_comparison,
and gaxis=-F0_Z/(Lambda*F0). For either specified history set define

```
S0 = -2delta*Z*JB+d*JB_Z+2(1+delta)*Z*JC-d*JC_Z
Pr = 2(1+delta)*Z*p0-d*p0_Z.
```

Using the comparison histories, the full signed axial drive is R*Q, where

```
Q = hydro+Pstar^2*Pr/(2L)
    +F0^2*(R*S0/(2L)-rho*gaxis*d*(JB-JC)/L)
hydro = [-Wbar*barU+(1-delta)*(JM-Z*JM_Z)/2
         +2delta*Z*JA-d*JA_Z]/(2L).
```

The actual root reductions are

```
b  = -sqrt(2R)*chi_R*Q/(F0*barphi)
t0 =  sqrt(2R)*Q/(F0*barphi*Dbar).
```

The native phi and Pstar cancel from b; chi_R cancels from t0. This is an
exact same-source cancellation, before any norm bound is taken.

The generic p2 instead uses the original native U and native histories:

```
T = (1-delta)*Z*JM+d*JM_Z
N = -U+(1-delta)*(JM-Z*JM_Z)/2+U*T+2delta*Z*JA-d*JA_Z
    +Pstar^2*Pr+F0^2*(R*S0-2rho*gaxis*d*(JB-JC))
p2 = sqrt(R/2)*N/(F0*phi*L).
```

The source algebra controls prove these reductions, including R*Lambda=rho.
Generic p2 is distinct from the comparison pressure drive and from b.
Physical y=log(R/Ra)=hb*s, hence R_y=R and rho_y=rho. The physical-y
prefactor derivative contributes 1/2 in the b/t0/p2 formulas. The first-chart
chi_y contributes at most 32/hb; on the second chart chi=hb and chi_y=0.
Those scales are retained in the logarithmic caps.

Let C denote the same selected positive logC. For each k=0,1,2,3 the saved
strict coarse envelopes certify the following. They are upper bounds on
log absolute norms, with the exact directed bounds retained separately.

| channel | strict log absolute norm envelope |
|---|---:|
| first b_Zk | < 2C |
| first b_yZk | < 102C |
| second b_Zk and b_yZk | < -98C |
| both t0_Zk and t0_yZk | < 2C |
| both p2_Zk and p2_yZk | < 2C |

## Complete inlets and native histories

Only two new all-coefficient core rows, (rho order,Z order)=(0,5),(1,5), are
formed. Existing source Z0..4 and direction Z0..3 intervals retain their
exact endpoints. New fifth source rows supply Dbar_Z4 and its mixed-y row;
no temporal coefficient recursion is substituted.

The incoming histories are the full original inlet integrals at rho=4.
Their derivative bounds use the global core Xh bounds on rho in [0,4.1],
not a point-only value of the core. M,A,C have inlet weight 1, while B has
weight t with integral 1/2 over t in [0,1]. Binomial product bounds and these
weights give their complete ordinary inlet norms. Positive Duhamel kernels
then give max(initial cap, source cap/rate), separately on the comparison
and native functions. No incoming memory is reset.

The native angular equation is enclosed first. Its positive phi floor and
ordinary jets through Z4 give the ratio phi/barphi. The exact prescribed
axial equation U_y=-chi_R*(phi/barphi)*R*Q is then integrated with dy=hb*ds
once, from the original core U inlet. The report retains these nonzero
absolute bounds for the increment from that inlet:

| ordinary Z order | upper display for native U minus core inlet |
|---|---:|
| 0 | 1.88331075230515922311608e-204453045017285003 |
| 1 | 2.8213682082544285878197e-204453045017284971 |
| 2 | 8.46307722170365733977452e-204453045017284939 |
| 3 | 3.80973215029831960463167e-204453045017284906 |
| 4 | 2.28704176870170693821618e-204453045017284873 |

The native M/A/B/C and their y jets then use their own sources, including
A_y=U^2-A, B_y=phi^2-2B, C_y=phi^2-C and M_y=U-M. The second comparison is
not identified with the core or with these native histories. The y jets are
derivatives of the exact source equations; no nonexistent p2_y DAG handle
is invented.

## Execution and continuation

New production exits 0 in 3.422 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_b_p2_mixed_bounds.json`
is 16,778,857 bytes, SHA256 `fc28f40dcf58e3a5123725361baca36005902f0733080c2fff6852293a703e5d`, with 2000 input and
362 implementation pins. Saved-only inspection,
source identity/AST gates, exact rational pressure controls, reduced algebra,
directed scale/envelope arithmetic and inherited-packet checks pass.
Independent fixed Luna/max workers reviewed bounded source/scale arithmetic;
separate runtime model metadata is unavailable and is not claimed.

All old 212 forward terms, six signed inner/lower sources, existing endpoint
identities and derivative-error packets remain unchanged. No accepted
producer, owner, norm computation, recurrence, quadrature or seam proof was
replayed. Existing dirty files and the unaccepted scaled-Poisson draft are
preserved.

The implementation is [the reduced source provider](../src/openai_ns_reconstruction/original_micro_b_p2_mixed_bounds.py),
[the fifth-order source extension](../src/openai_ns_reconstruction/original_micro_source_Z5.py),
[whole-tube pressure](../src/openai_ns_reconstruction/original_pressure_axial_tube.py), and
[the report runner](../scripts/run_original_micro_b_p2_mixed_bounds.py).

Inspect the saved report without replay:

```
python scripts/run_original_micro_b_p2_mixed_bounds.py --inspect-report
```

NEXT bound the actual q=cut*sqrt(gamma/(2a)) mixed derivatives using the
source-owned active/flat branch, positive gamma/eta and normalized root
bounds. Here cut=sigma_flat(1-Delta/eta), Delta=kappa-2 and
gamma=2+2eta-kappa; the flat source branch has q=0 exactly. Then control
actual loop/inverse-phase derivatives and finite-N
errors, Cy_Z3, joint mean-Z3/QZZ and outer RQy. Independent reserved-pulse
terminal functional matching and global five-moment closure remain open,
as do exact heat/energy, stress/flat remainder, true n>=2 recursion/summation,
realizable corrections and independent corrected Cartesian NS validation.
Goal ACTIVE.
