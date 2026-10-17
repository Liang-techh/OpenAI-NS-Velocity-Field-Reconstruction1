# Signed actual core direction on a continuous axial tube

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. A new directed model calculation plus the already saved same-source fixed-point perturbation proves actual Phi_rho<0 and Dbar_core=-2rho*Phi_rho/Phi>0 on the entire rectangle rho in [4,41/10], Z in [37/100,93/250]. This is a quantitative prerequisite for the micro root provider. It does not establish the second-micro comparison direction, a complete root provider, numerical mean-Z3/QZZ, or global five-moment closure.

## New actual intervals

The following are deliberately coarse outward bounds, checked against the exact serialized directed endpoints. They hold for every point of the stated rectangle, rather than only sampled axial values.

| actual normalized core quantity | enclosing interval |
|---|---:|
| Phi | [0.2536,0.3005] |
| Phi_rho | [-0.1221,-0.1152] |
| Dbar_core | [3.069,3.946] |
| ordinary Dbar_core_Z | [1.04e-48,2.03e-48] |
| ordinary Dbar_core_ZZ | [-1.83e-47,-9.3e-48] |

At the exact retained query rho=4,Z=371/1000, the displayed approximations are PhiÃ¢â€°Ë†0.2829799868805425, Phi_rhoÃ¢â€°Ë†-0.1198820205377527 and Dbar_coreÃ¢â€°Ë†3.389130711589433. These displays are not exact scalar values or substitutes for the report's two exact MPF endpoints. The saved original query interval lies inside the new axial rectangle. No midpoint of an error cap is selected.

## Same actual source, model and error

The original defining functions give

```
H0=(1-delta)Z/2+(1-Z^2)(4Z+j)
sigma=j/500
chi_core=H0^2/(H0^2+sigma^2)
Phi_model=sum_n (-chi_core*rho/2)^n/[n!(n+1)!].
```

Here j is the admitted actual source parameter from core_transfer.required_j. The current source delta is exp(-4*(exp(40)+11)-30). The shared analytic tube's delta=1e-200 is a universal cap, not the actual selected parameter. The new calculation uses the actual exponential recipe, verifies it is below that cap, and consumes the accepted exact graph/oracle parameter identity. Matching boxes are not used as function identities. The core seed chi_core is distinct from the micro transport control chi_R.

The new helper sums 64 terms with 600-bit directed arithmetic. Its ordinary rho/chi derivatives use falling factorials, followed by the ordinary Z chain rule. The absolute omitted-term ratio is q_abs*(n+1)/[(n+1-i)*(n+1-ell)*(n+2)], decreasing beyond the retained terms and strictly below 1. The full infinite model series is enclosed by its finite directed sum plus a symmetric geometric tail. Ten ordinary mixed rows with i+k<=3 are retained. Small 1-chi_core is preserved by evaluating 1-sigma^2/(H0^2+sigma^2), without assigning chi_core=1.

Each actual mixed row is the model enclosure plus/minus the UPPER endpoint of its corresponding saved Phi_minus_model_absolute_bound from compliant_core_uniform_bounds.json. The lower endpoint of a norm-bound record is never treated as a positive lower bound on the actual error. No scalar cap is differentiated; each derivative uses its own saved row. For example, the retained radial derivative perturbation upper is about 3.6971e-204453045017284960. Its exact endpoint and every old error remain in the report.

The already accepted common-core identity gives Dbar_core=-2rho Phi_rho/Phi. In particular, at fixed rho,

```
Dbar_Z=-2rho*(Phi_rhoZ/Phi-Phi_rho*Phi_Z/Phi^2).
```

There is no derivative of rho with respect to Z. The new symbolic controls verify ordinary Z0..2, with no factorial normalization. This dimensionless logarithmic slope is R*D_over_R for the common core under rho=Lambda*R; F0 is independent of rho and cancels. No extra Pstar factor is inserted.

## Micro coverage and remaining cross-map

New arithmetic on the saved actual log(hb) proves 0<hb<1/1000. For 0<=s<=2, rho=4exp(hb*s)<2000/499<4.1, using exp(x)<=1/(1-x). This covers both micro charts' calls to the core function without materializing hb or the fast phase. It is a new domain proof; the older tiny collar's domain estimate is not extended by assertion.

The comparison profile uses alpha=1-sigma(s-1), with six continued comparison histories. On the first micro alpha=1, so the profile/core FTC and same-inlet history identities provide the route to a current C6/C5 direction cross-map. That current source cross-map is the next executable step. After the cutoff begins on the second micro, the comparison direction cannot simply be replaced by the core logarithmic slope. The new receipt explicitly leaves second comparison positivity and the full C5 mixed-root provider false. Ordinary core derivatives through total order 3 do not supply the required whole-micro Z3/slow-y mixed-root bounds.

The parent 212-term forward packet, six signed inner/lower contributions, nonreset means/pressures, all existing source correlations, every error cap and every completion limit are unchanged. This stage adds no m15 contribution and performs no terminal subtraction. P1 remains distinct from m15.

## Execution and handoff

New production exits0 in 3.14 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_signed_direction.json` has 16,360,725 bytes, SHA256 `18968d9762341c53e0f161193d2e998126abb4cd22d73e3df393a62d1cb4d12e`, 1997 input and 352 implementation pins. Focused syntax/LF checks, exact source/parameter guards, model derivative/tail controls, symbolic direction derivatives, parent-packet invariants and saved-only inspection pass. Bounded Luna/max workers implemented the model and independently reviewed the source mathematics and directed arithmetic. The fixed role is reported by the collaboration tool; separate runtime model metadata is unavailable and is not claimed.

No accepted owner, builder, fixture, recurrence, norm producer, old quadrature, 14-row endpoint proof or 8-row phase seam proof was replayed. Only new scalar series/arithmetic and saved receipts were used. The separate original_scaled_poisson_root.py helper remains an unaccepted untracked draft; it is not used in this accepted report.

NEXT bind the new signed core rows to the current first-micro comparison D_over_R source using alpha=1 and all six same-inlet histories. Then quantify the second-micro comparison deviation and supply the actual root Z0..3/slow-y0..1 provider, retaining chi-normalized root variables and the positive logarithmic defect 1-|r|. Enclose Cy_Z3, the joint nonreset mean/QZZ volume source and outer RQy before advancing the forward Rp/independent reserved-pulse terminal comparison. No old C3 cap transfer without a same-function quantitative proof.

Global five-moment/common-N closure, exact heat/finite energy, admissible stress/flat remainder, genuine n>=2 temporal recursion/summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. Goal ACTIVE.
