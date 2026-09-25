# Dyadic transfer of the complete bridge defect

`experiments/root_st073/adaptive_bridge_recursive_defect.py` evaluates the
unforced full Cartesian momentum residual on the *same* 15 bridge similarity
nodes at `k=11,13,15,17,19`. It compares the unmodified bridge, the shared
24-mode correction, and the two-sided cone-aware correction. Here
`tau = 2**(-k-1)`, and the diagnostic normalized defect is `tau**1.5 R`.
This normalization tests the observed approximate `tau**(-3/2)` growth;
it is not a proof that the exact solution has that similarity law.

| Bridge | Maximum at `k=11` | Maximum at `k=19` | Growth exponent per halving | Normalized 15-node vector L2 at `k=11 -> 19` |
| --- | ---: | ---: | ---: | ---: |
| Unmodified | `7.288e5` | `2.906e9` | `1.4952` | `5.377 -> 5.291` |
| Shared 24-mode | `3.909e5` | `1.581e9` | `1.4977` | `3.342 -> 3.296` |
| Two-sided cone-aware | `3.264e5` | `1.312e9` | `1.4967` | `3.238 -> 3.207` |

For the cone-aware field, successive two-halving ratios of the normalized
sampled vector L2 are `0.9960, 0.9970, 0.9981, 0.9991`. The ratios are
moving toward one, so this correction mostly lowers the prefactor. Across
two halvings, `tau**(-1.5)` multiplies an unchanged normalized defect by
eight. Even to keep the *sampled absolute vector norm* from growing, the
normalized defect would need a ratio at most `1/8`; the measured ratios
are nowhere close. The sampled maximum likewise rises by about `4,022`
between `k=11` and `k=19` in the cone-aware field.

This identifies a concrete gate for the next construction: its
scale-dependent background, radial moment repair, and wave/mean stress
update must reduce the normalized *complete* bridge defect at each
dyadic transfer, not merely refit a fixed compact basis at a few time
knots. The [OpenAI paper's Appendix A and Sections 5, 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
provide the structural ingredients for that update; the present bridge
does not implement them.

These are 15-node finite-difference results within the experimental
adapter's extended `k <= 20` window, beyond the original local core's
registered `k <= 6` interval. They are neither continuum bounds nor
spatial-volume L2 estimates. The cone-aware fit had not met its optimizer
termination condition and failed off-scale holdouts; it remains a diagnostic
candidate. No scale recursion or PDE acceptance is claimed.

The separate outer-moment audit and attempted coupled repair are in
`ST073_OUTER_MOMENT_REPAIR.md`. They show that the current mean modes also
leave a large residual-stress tail at the bridge's outer edge.
