# Current core radial production source proof - 2026-10-05

The new `current_core_recurrence_source` stage proves that the two production radial recurrence implementations express the same equations after exact source scaling. It also proves the leading axis A/Uz identities and the compliant axis pressure seed units. It does not admit the actual scaled swirl-source jet or complete the core/first bridge join.

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcorerecurrence
```

Receipt: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_recurrence_source_check.json`.

## What is proved

The exact units are Atilde_n=Lambda^-n A_n, Utilde_n=Lambda^-n U_n, Ptilde_n=epsilon_core Lambda^-n P_n, elltilde=epsilon_core ell_phys and Stilde=epsilon_core^2 F0^2. Lambda is independent of Z and epsilon_core Lambda=1. The core scaling epsilon is distinct from the pressure-datum parameter epsilon=.001 delta.

Four derived RHS transformations, one inlined amplitude derivative and the original positive-epsilon domain guard turn the entire original advance_one AST into the actual advance_scaled_one AST. All other helpers, validations, Taylor convolutions, radial loops, cross terms, pressure terms and swirl terms are structurally identical. Thirteen homogeneous additive term families prove the scaling for arbitrary radial n and axial Taylor length, wherever the original L denominator is nonzero. The constant W0/H0 terms retain zero radial weight.

Both unmodified production functions also pass36 exact scalar identities with independent formal Taylor jets at n0/1/2/3. These bounded fixtures check the executed operators; they are not the arbitrary-n proof. No interval overlap or parameter midpoint is used.

The original core axis expressions and normalized field formulas give Atilde1=-(chi+epsilon_core beta)/4 and Utilde1=epsilon_core*(-g/(2L)) as exact Z-function identities. The seed's ell is already scaled and must not receive a second epsilon. The pressure seed is exactly epsilon_core Pstar^2 times the same datum's normalized pressure coefficients. The current bridge and rebuilt core parents retain the same actual-family, implicit-source and datum identities.

The algebraic pressure row Ptilde1=Stilde is proved only for a formal exact Stilde. The existing S_Z_taylor input is an enclosure; this stage deliberately does not relabel it as the exact nonlinear source.

## Next executable work

1. Bind exact S(Z)=epsilon_core^2 exp(-2 selected_logCstar-2 Lambda G(Z)) to the production anchored amplitude source, selected norm-family definition and the same complex tube.
2. Connect the existing S_bound consumer to transfer.F0_squared_Xh_norm_upper including its Cauchy_weight_supremum_upper. The current weight is1, but the new admission must reject or adjust a changed weight; do not silently omit it. Match the actual eta/2 Cauchy radius and ordinary Taylor coefficient convention, with no extra factorial.
3. Prove finite recurrence residual/tail enclosures belong to the same admitted nonlinear fixed point as normalized_jets. Shared hashes, object identity and the contraction record alone do not complete this proof.
4. Match pressure primitive 4C against PD+PI, then replay actual phase0 controls and hydro/pressure/swirl drives to obtain the core/first functional mixed4 join.

All33 regional physical source maps and three non-core bridge joins remain accepted at11a35cc9. Complete nonlinear points, full interface smoothness, uniform native pulse C4, global stress/cone/flat remainder/required-domain energy and genuine n-dependent temporal recursion remain open. Radial Taylor coefficient generation is distinct from temporal recursion.
