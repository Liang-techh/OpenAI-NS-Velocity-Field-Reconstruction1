# Compact bump normalization: closed form — 2026-10-29

Goal ACTIVE / INCOMPLETE. The compact raw bump integral has a direct special-function representation. A new independently frozen leaf statically checked the exact beta definition, its `[-1,1]` integration bounds and `normal=integrate(beta)` in the current Rc controls source. It did not import a historical constructor, traverse a source graph, invoke an ancestor, or perform numerical quadrature.

For the source-declared bump `beta(t)=exp(-1/(1-t^2))` for `|t|<1`, zero otherwise, write

```
J0 = integral_-1^1 beta(t) dt
   = exp(-1)*sqrt(pi)*U(1/2,0,1)
   = exp(-1/2)*(K_1(1/2)-K_0(1/2)).
```

The first identity follows by symmetry and `s=t^2/(1-t^2)`: `J0=exp(-1) integral_0^infinity exp(-s) s^(-1/2)(1+s)^(-3/2) ds`, the integral representation of Tricomi U. The Bessel form follows independently from `t=tanh(u)`, `x=2u`, integration by parts using `d[sinh(x)/(1+cosh(x))]/dx=1/(1+cosh(x))`, and the integral definition of K. Boundary terms vanish; the original bump is absolutely integrable and flat at both endpoints. These are analytic function identities for the stated definition. The static syntax check is not a proof that every archived provider graph node has the same runtime semantics.

The new callable `src/openai_ns_reconstruction/compact_bump_normalizer_closed_form.py` evaluates both forms in separate mpmath contexts at caller-selected precision, without numerical quadrature. The single frozen fixture returned successfully: all 12 checks passed, including the three source-AST bindings and agreement at 50 and 90 digits. Its 90-digit Bessel evaluation begins `0.4439938161680794378230489211705526637612`; the cross-form discrepancy at that precision was about `1.23e-91`. These are point approximations and diagnostics from the same mpmath library, not a certified enclosure or proof of independent provider implementations. No actual selected controls, amplitude, finite-N field, m15 value or pressure was evaluated. The source node's interpretation and whole-Z source/provider binding remain open.

Static review of `W/current_physical_measure_source_congruence_kernel.py` found caller-supplied same-source assumptions presented as unconditional flags; it is unaccepted and was not executed. A newly interrupted draft `src/openai_ns_reconstruction/pulse_signed_scaled_kernel.py` was also not executed or accepted after existing `pulse_correlated_density_delta.py` and `pulse_scaled_pair_quadrature_v2.py` were found to cover its intended functionality. Both drafts remain untracked. Preserve them as drafts; do not reenter them as prior proof or substitute their conclusions for source bindings.

For m15, the producer of historical F's original gap radius and finite endpoint tensor is identified, but the saved frame rewrite does not establish the rewritten local endpoint owner or capture-free ordinary Z transport. Current geometry metadata is not automatically the historical field radius. Retain `P-F=(P-G)+(G-F)`, all four physical correction roots, the independent reference offset, 212/213/183 inherited counts, six signed channels, the once-only reserved term and every old error. P1 and M remain distinct from m15.

Next: bind historical F's actual field radius and endpoint tensors through the inlet substitution, including source-owned `mu`, physical measure, orientation, Pstar and ordinary Z1; then use already existing correlated-density and normalized pair integration interfaces only with appropriately source-bound signed jets and separate error transport. Obtain actual compact amplitude and controls with provenance before using J0 to compute a physical source. Global five-moment closure, heat exterior, admissible stress/remainder, true time-order recursion and Cartesian NS validation remain incomplete. Runtime model metadata was unavailable (`null`).
