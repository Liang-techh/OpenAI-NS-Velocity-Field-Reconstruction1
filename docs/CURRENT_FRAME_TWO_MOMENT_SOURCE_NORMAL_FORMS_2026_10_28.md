# Source-defined pulse endpoint normal forms for both moment rows — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new frozen leaf was entered once and returned eight true symbolic check flags, including a duplicate generic endpoint check described below. It checked the byte hashes of the original kernel and own-field source against the saved own-field receipt, and located specific defining source statements. It did not construct, traverse or evaluate the archived expression graph, call an ancestor, compute an integral or remove an error term. The inherited receipt is not a new proof of the complete archive's semantics.

For each radial moment index `i=1,2`, the source defines `lambda_i=1/2-i*mu`, with `mu>0`. Write `alpha=logpref` and define, without evaluating it,

```
J_i = integral_0^11 exp(lambda_i*xi/mu) gp(xi) dxi
W_i = exp(-13*lambda_i/mu-alpha)
Pi_i = W_i*J_i/mu
H_i = exp(alpha-log(mu)+2*lambda_i/mu)
E_i = exp(-11*lambda_i/mu).
```

Here `J_i` is distinct from the beta-normalization integral called `J` locally in the original kernel source. The equality for `Pi_i` follows from its defining integrand. In the main chart, change variables `xi=mu*v` in its history integral at `xi=11`; the normalized moment is `E_i*(M_i+ap*J_i/mu)`. The gap formula at the same coordinate is `-H_i*sum_j(C_j*row_ij)` and `H_i*mu*W_i=E_i`. Therefore the two normal forms agree **if** the actual selected coefficients satisfy the source row identity `sum_j(C_j*row_ij)=-mu*(W_i*M_i+ap*Pi_i)` with the same inlet and kernel functions. These expressions are exact function definitions and conditional identities, not numerical values or certified enclosures.

The source uses `row_1=A` and `row_2=A+mu*D`. The second row is finite-mu data: from the original integral definition of `D`, `mu*(-s)*exprel(-mu*s)=exp(-mu*s)-1`, with the continuous `exprel(0)=1`, so the integrand of `A+mu*D` becomes the normalized beta integrand times `exp(lambda_2*s)`, since `lambda_2=lambda_1-mu`. The leaf checks an abstract exponential identity and linear distribution; the interpretation of `exprel`, integrals and the source program is a mathematical derivation, not an independently evaluated graph. The result keys named `first_row_endpoint_given_selector` and `second_row_endpoint_given_selector` run the identical generic symbolic expression: the second does not explicitly substitute the second-row source or its selector. The generic algebra can be applied to the second row under the stated premises, but that application and the actual second-row selector were not independently checked by this leaf. Neither flag is a numeric check.

The pinned source definitions also give the same main/gap `y=xi/mu`, `R=exp(logRp+y)` and `theta=Pstar*native_u*exp(-(1/2+mu)*y)` before inlet rewriting. The original flat `sigma` definition gives `gp(11)=gp_xi(11)=0`, hence `gp_y(11)=mu*gp_xi(11)=0`, under the stated source-function interpretation. The leaf records these source-derived formulas and anchors relevant definitions; it does **not** independently differentiate the flat function or check every archived expression edge.

Actual frame endpoint equality remains unproved. In particular, the original inlet rewriter copies integral binder names as metadata while recursively replacing expression operands; a general capture-free argument and ordinary Z derivative transport through order two are still missing. A substitution of a free `v` into `integral_v(x+v) dv` illustrates why a name-based rewrite alone cannot establish capture freedom: without alpha-renaming, the free `v` can become bound. This example is a semantic warning, not a claim that this actual graph contains such a capture. The `m` used by the physical `Q` formulas is the first moment; the second normalized history is separately retained.

No actual `Z`-jet, frame seam, m15 cancellation, numerical value, certified error improvement, global five-moment closure or later physical layer is claimed. The once-only reserved contribution, nonzero histories, independent reference offset and prior errors remain unchanged; `P1` and `m15` remain distinct. Runtime model metadata was unavailable (`null`).

A parallel reviewer reported an execution-rule violation: an unfrozen inline Python invocation decoded three gzip archives without first authenticating their raw hashes and traversed geometry-node DAG closures to look for `Z`. Its exact code and complete stdout are unavailable; its self-reported roots and scope are preserved in `experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_geometry_uncheckpointed_incident.json`. This was a source consumer despite the reviewer's initial description as static research. Its output is unaccepted and is not evidence for the formulas or checks above. Do not replay that geometry-closure family. The separately frozen leaf in this document did not use the reviewer's graph output.

Next: freeze a new bounded lexical-scope audit of the specific actual inlet substitutions and needed main/gap endpoint roles, using distinct binder identities and free-variable sets with fail-closed collision reporting; do not repeat the uncheckpointed geometry-closure family or the prohibited earlier graph-comparison or physical-root traversal families. Then prove ordinary Z transport and connect the saved source rows to these endpoint normal forms before claiming a seam identity.
