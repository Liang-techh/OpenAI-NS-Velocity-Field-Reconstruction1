# Actual twenty-term core operator and radial extraction

The normalized core now has a callable production source operator. It consumes the original gauge's ten angular and ten axial expressions, uses the current admitted swirl source, restores the common pressure integral, and generates finite radial coefficients. This advances the core solver; temporal scale recursion remains open.

Run the scoped stage:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcoreoperator
```

Implementation and receipts are `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py`, its producer JSON, and `_check.py` / `_check.json`. Existing source, amplitude, radial-scaling and physical-field receipts are consumed by hash; they are not regenerated or relabelled.

## Implemented behavior

- The evaluator executes only the original pure symbolic gauge prefix, then interprets its actual expressions in an ordinary Taylor coefficient ring. Each expression is paired with its exact production majorant formula. Swapping source terms without swapping their bounds is rejected.
- `source_terms` restores `F0^2=S_scaled/epsilon_core^2`, `ell_scaled=-g=epsilon_core ell_phys`, and `P0_phys=P0_step/epsilon_core`. `M=integral Psi/rho` and `Pcal=F0^2 integral Phi^2` are computed from the same fields. `Pcal_Z` differentiates this same product and primitive, including the source derivative.
- `green_inverse` implements the full zero-axis inverses of `2D2+chi` and `2D1`. Raw J1/J2 have no one-half; the full equation contains it once. The angular inverse retains the axial convolution with chi. The Bessel Phi model and linear Psi model keep their original axis values.
- `correction_map` implements `T=(epsilon/2 R J2 Etheta, epsilon/2 J1 Ez)`. `next_rows` extracts the full Phi, normalized Uz and scaled pressure row. One extra axial coefficient is required and consumed, rather than filling an unavailable derivative with zero.
- The exact identity `Psi_axis=0` is imposed structurally after confirming the common U0 seed. Subtracting two copies of a nondegenerate U0 interval would lose correlation and introduce a spurious correction amplified by epsilon inverse.

The fresh current-source runtime at Z=.293 generates radial rows0..4. A checked adapter is also callable at a new Z. These are finite coefficient enclosures of the formal exact source, not a completed nonlinear point field.

## Evidence and precise scope

Twenty-one independent formal scalar coefficients at radial source orders0,1,2 agree exactly with unmodified `advance_scaled_one`. The pressure rows for this comparison are constructed by a separate convolution. Twenty-four independent Green equation identities pass. The Bessel and linear model equations pass. All-order coefficient extraction follows the production gauge identity, ordinary derivative/product/integral coefficient laws, and the prior accepted whole-production AST scaling.

The norm proof uses the paper8.30 **supremum** over radial order, axial order and Z. It does not substitute an l1 coefficient norm. Vandermonde and both complete convolution sums give product constant256. The mixed source `Jnu((M Psi)_Z rho f_rho)` is bounded by `20*4*16*16/h=20480/h`; no additional independent radial80 factor is required. Fixed Cauchy coefficients commute with Jnu and contribute their single summed factor `(1+r)/(1-r)^3`. The cross-swirl bound retains H through the exact partial-fraction identity for `H/(H^2+sigma^2)`.

Accepted gates:

- `current_core_actual_twenty_term_operator_and_radial_extraction_certified=true`.
- `current_core_actual_twenty_term_norm_incidence_certified=true`.

Still false: `current_core_actual_operator_contraction_and_tail_certified`, `finite_rows_and_tails_bound_to_same_nonlinear_fixed_point`, core/first and all-four bridge joins, complete nonlinear point evaluation, global stress/tensor admissibility, and temporal recursion. Norm incidence is not a blanket admission of these later claims.

## Next agent work

1. Admit the actual callable operator on the same Xh model-centered ball. Replay the stored twenty positive majorants, size and Lipschitz sums at their declared interval precision. Bind current parameters, pressure datum, selected amplitude, Cauchy disks, fixed coefficient moduli and commuting resolvent to those exact definitions. Keep analytic existence distinct from finite implementation enclosures.
2. Use the Banach solution's unique radial coefficient equations to identify the production recurrence coefficients. Prove finite output inclusion for all admitted parameter/axial domains; no midpoint or interval-overlap identity is permitted.
3. Bind the same solution's nonlinear norm tail and Bessel model tail to the actual `normalized_jets` and rooted-core consumers. Keep ordinary Taylor and ordinary derivative factorials explicit. Only then close the common finite/tail gate.
4. Connect the common pressure primitive to original bridge `4C=PD+PI`, then the original phase0 stress-free recovery and drive. Retain exact hb pullbacks through mixed4; admit core/first only after the equations match.
5. Continue global stress, flat remainder, required-domain energy and true n-dependent temporal recursion after the shared background is admitted.

Do not rerun unchanged121 swirl seed rows, prior33-region physical reports or the earlier36 radial-scaling identities. Reuse their hashes unless a defining source changes. Mark each completed task with its scoped receipt and commit.
