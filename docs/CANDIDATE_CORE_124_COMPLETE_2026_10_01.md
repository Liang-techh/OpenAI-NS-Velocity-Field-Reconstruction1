# Degree124 candidate core: combined normalized mixed C3 gate passed

The refined-pressure Lambda=1e48 coupled candidate recurrence has completed all124 radial steps, with128 initial axial Taylor coefficients and260-digit directed arithmetic. Four axial coefficients remain in the final row, sufficient for three axial derivatives. Pure recurrence arithmetic totaled approximately621.3 seconds; checkpoint serialization adds wall time.

The combined receipt includes propagated finite coefficient uncertainty and the infinite analytic radial Taylor remainder. Its equation-identification dependency verifies the exact angular and axial fixed-point substitutions; radial coefficient uniqueness identifies the analytic solution's coefficients with the directed recurrence relative to the accepted stored source.

For every mixed derivative of total order at most3, over0<=scaled_R<=4.1 at axial center Z=0.3, the largest combined midpoint error bounds are:

| Normalized component | Largest mixed C3 error bound |
| --- | --- |
| Phi=F/F0 | 1.84246630169076e-23 |
| Psi=Lambda*(Uz-U0) | 7.51073122007975e-13 |

All20 derivative/component budgets pass the internal1e-12 target. This is a normalized leading-core computation relative to the accepted stored datum. It is not the full corrected Navier–Stokes residual target, a whole-axis finite-field certificate, or a physical-coordinate error certificate.

At the provisional exit scaled_R=4 and the same axial center, the retained field's angular inlet trace Ttheta/F is enclosed in approximately[-1.75612e-63,1.75612e-63]. Two directed expressions from the integrated angular equation overlap. The trace has not been fitted or reset. Containment of zero is not an exact matching certificate; the analytic core has zero trace by its exact equation, while the retained polynomial is an approximation.

## Next dependency

Value positivity of the exact analytic candidate is already established throughout scaled_R in[0,4.1] and Z in[-1,1]. The next unresolved core condition is the whole-axis radial derivative sign, especially near the root of H0 where chi is small. The paper's shifted Bessel model gives the exact axis slope -(chi+epsilon*beta)/4. A weighted estimate for the first anchored correction and second-order remainder is needed; the existing generic correction norm does not prove this sign near the root.

After this condition is settled, regenerate the shared exit/collar data from the accepted candidate, restore five-moment terminal identities as functions of Z, and certify pressure compatibility and the full admissible stress cone. Spatial radial recurrence remains distinct from the later n-dependent temporal recursion, which is still open.

Artifacts:

- `lei_ren_part1_paper_candidate_gauge_core_lei_ren_part1_paper_candidate_pressure_axis_jets_refined_state.json`: all coupled exact interval rows.
- `lei_ren_part1_paper_candidate_combined_core_budget.json`: all20 combined derivative budgets and input hashes.
- `lei_ren_part1_paper_candidate_inlet_trace.json`: directed exit trace at Z=0.3.

Recompute the final diagnostics:

    python experiments/root_st073/lei_ren_part1_paper_candidate_combined_core_budget.py
    python experiments/root_st073/lei_ren_part1_paper_candidate_inlet_trace.py
