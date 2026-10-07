# Original whole-period Z supports tighten every terminal derivative bound

Successor: [CURRENT_NATIVE_PAIRED_C1_2026_10_07.md](CURRENT_NATIVE_PAIRED_C1_2026_10_07.md) adds exact paired Poisson derivative supports, improving all five terminal Z upper bounds again and refreshing weighted source contributions. Genuine linear q_Z and actual controls/global N/closure/recursion remain open.

Checked source: [96916dcb](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/96916dcb1a684253a5261aa224575c985a1815be). Same original source family, candidate N2048, full Z[-1,1] and24 continuous radial cells as the preceding [direct q² stage](CURRENT_NATIVE_COLLECTED_Q2_2026_10_07.md). Five original terminal Z absolute upper ranges now strictly improve. No strict C0 target reduction is claimed in this step; previous C0 support improvement remains.

## Original derivative theorem

The original primitive identity A=a/2*(phi-psi/(2pi)) and inverse phase give

```text
psi_Z = (2pi*phi*nu_Z - T2_Z|psi)/(1+t²)
|A_Z| <= |a_Z|/2 + a*L/(4pi)
L = 2pi*|nu_Z| + |T2_Z|psi
```

The existing `serial.periodic_parameter_theorem()` proves the original signed/small Poisson estimates |theta_u|<=2, |W1_u|<16pi and |(h^-2*W2)_u|<1000pi, uniformly over the complete period and either sign of r. Its exact Mobius denominator identity controls the angle derivative without taking large independent powers of a nearly vanishing Poisson denominator.

Use T>=|t0|, TZ>=|t0_Z|, Q>=|q|, QZ>=|q_Z|, UZ>=|u_Z|, R>=|q²| and RZ>=|(q²)_Z|. With q² and its Z row supplied by the original branch-local direct jet,

```text
nuZ_upper = 2*T*TZ + 2*RZ
firstZ = 4pi*(TZ*Q + T*QZ + T*Q*UZ) + 16pi*T*Q*UZ
T2Z_upper = 4pi*T*TZ + 4*firstZ + 400pi*RZ + 4000pi*R*UZ
```

The former4 q q_Z and800pi q q_Z terms become2 R_Z and400pi R_Z before absolute ranges are taken. This retains the original derivative function and its correlations; it is not a new cutoff. Genuine linear q_Z and u_Z terms remain in firstZ and the B derivative.

The original B=E*beta product rule uses independent proved C0 bounds |A|<=5/4 and |beta|<=3/2. Its moving-angle term is controlled by the original identity q*h^-1*w=(t-t0)/2 and |q*h^-1*w|/(1+t²)<=(1+|t0|)/2. Original E_Z, a_Z, b_Z, p2_Z and all cross terms remain. No derivative is inferred by differentiating a C0 support cap.

On an active original cutoff branch, a<=kappa<=2+eta<=5/2 and |b|<=v/2<=3/2. Only their C0 outer ranges are restricted, retaining any tighter original range. The original nonzero derivative rows stay intact; original positive a lower proofs remain the denominators. A source proved q=q_Z=0 gives exact zero primitive derivative supports. Otherwise the original active derivative expressions are retained, including cutoff-crossing sensitivity.

## Implementation and measured result

Files use `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_periodic_C1_support_transport.py`: `periodic_C1_support`, branch-local original consumer, density/oracle and complete24-cell transport.
- `current_native_periodic_C1_support_transport.json`: actual source/whole-period support decisions, continuous transport, terminal targets and baseline comparisons.
- `current_native_periodic_C1_support_transport_check.py` / `_check.json`: direct independent original scalar support comparisons, source/context/ledger/memory guards and actual target comparisons.

`NativePeriodicC1Density` evaluates the unchanged original cutoff/signed-u phase body with the accepted direct q² adapter. It intersects only A_Z/B_Z outer ranges before density. Original C0 A/B objects remain; the earlier C0 support layer then applies before nonlinear density as before. A tighter original derivative range is retained as the same object. The exact original source sign, basis, ledger, roots, global radius phase, E/V/cross terms and separate P0/P0_Z remain. Nonlinear evaluation precedes every overlapping branch hull.

There are138 actual branchwise derivative decisions: 80 use the new tighter support and 58 preserve a tighter original range. All24 continuous original cells and240 live integral C0/Z rows resolve. All five terminal Z upper bounds strictly reduce and no terminal upper bound worsens.

The table reports approximate **logarithms** of absolute upper ranges. Ratios are ratios of those logarithms, not ratios of actual target errors:

| Original Rc target | Previous log upper | New log upper | Previous/new logarithm |
|---|---:|---:|---:|
| M_Z | 1.916991e+408906090034569677 | 8.147214e+408906090034569676 | 2.352941 |
| I_Z | 1.763918e+408906090034569677 | 6.630788e+408906090034569676 | 2.660194 |
| S_Z | 1.916991e+408906090034569677 | 8.147214e+408906090034569676 | 2.352941 |
| Cp_Z | 1.783947e+408906090034569677 | 6.83107e+408906090034569676 | 2.611518 |
| D=(J-M)/mu_Z | 1.916991e+408906090034569677 | 8.147214e+408906090034569676 | 2.352941 |

The resulting bounds remain enormous and insufficient for useful repair contraction. They are conservative source-function bounds, not evaluated target values or controls. The original linear cutoff sensitivity from extremely small eta remains; source-range progress has not established a globally compatible N or coefficient recursion.

Focused acceptance:54 unchanged original scalar phase/primitive/y/Z/phi comparisons, including12 direct A_Z/B_Z bound comparisons (before taking the minimum with the original first-jet enclosure). Original negative/positive Mobius, small nonzero r, p2 crossing0, transition, flat and symmetry cases remain. Checked original q² derivative prerequisites are reused, not reconstructed. 138 actual derivative range decisions,240 live integral ranges and original quiet pressure memory pass. Producer 99.266s / checker 5.109s. Read-only review: **GPT-5.6 Luna / max**, scoped PASS.

## Detailed next-agent work

- [x] **SOURCE3c1-q²-Z/C1-support:** six original branch-local mixed q² rows, accepted five-site consumer, uniform whole-period ordinary Z primitive support and all five actual terminal Z upper improvements. Reuse these producers and source receipts.
- [ ] **SOURCE3c1-paired-q/u:** for P=h^-1*W1, H=h^-2*W2 and u=c*q, c=p2/dstar, retain `(q*P)_Z=q_Z*(P+u*P_u)+q²*c_Z*P_u` and `(q²*H)_Z=q*q_Z*(2H+u*H_u)+q³*c_Z*H_u` before hulling. Prove bounded original signed/small combinations, retaining the correct Fourier/Poisson tails. Compare local source supports and actual terminal bounds; do not zero genuinely linear q_Z terms.
- [ ] **SOURCE3c1-range-loss-refresh:** consume the current live24-cell objects to compute the largest weighted target C0/Z contributions after the new supports. Rebase original/common coordinates correctly and account for each downstream attenuation. Do not repeat the accepted source-coverage inventory. Use the result to choose the next refinement.
- [ ] **SOURCE3c2-joint-source:** form the full original joint integrand `Dk*Mk*f_k-Ac*Dm*Mm*f_m` before separate key hulls. Preserve distinct m rate1 and k rate3/2 masses, Ac_Z, background C_Z and division by mu*Ac². One common mass for their difference is invalid. Retain original V/F/E/B cross terms through the Z derivative.
- [ ] **SOURCE3c3-dominant-refinement:** split only the currently dominant original source cell, with exact selected sc/native endpoints, original global N*y phase, genuine microscopic widths and incoming memory. Continuous source covers, not midpoint values, must define the integral range.
- [ ] **SOURCE3d-compatible-N:** preserve original N^-1/N^-2 coefficient functions before fixing N. Bind the current source support dependence to an all-N route and prove finite log N meets the repair inequalities. Huge finite target bounds or fixed-N2048 scaling do not establish the needed coefficient functions.
- [ ] **CONTROL1a:** exact typed five bump matrix, Gram/quadratic functions and inverse action; complete first-Z Picard map with current original target functions. Preserve `B_exact*h+N*r+Q_exact/N=0`. Caps are sidecars and do not define h.
- [ ] **CONTROL1b/2:** actual convergent C1 controls on full Z, repaired original velocity and functional terminal moment identities, including endpoints/support joins, pressure and heat datum compatibility.
- [ ] **REC/HIGH/OUTER/WAVE/PHYS:** genuinely distinct n=1/n>=2 n-dependent radial recovery, common core definition interval and independent moment repair; higher jets/joins and full cone/pressure/heat/energy; smooth sum/flat remainder; two oscillatory families and averaged stress cancellation; corrected Cartesian NS and measured scale recursion/material winding.

Gate: `current_original_whole_period_Z_support_before_fixed_N_transport_executed`. Actual controls, globally compatible N, terminal closure, coefficient recursion, oscillatory corrections and corrected NS gates remain false. Keep the long-term goal active.
