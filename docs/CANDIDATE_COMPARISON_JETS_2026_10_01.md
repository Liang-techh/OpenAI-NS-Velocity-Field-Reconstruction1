# Lambda120 candidate comparison with analytic axial jets

`experiments/root_st073/lei_ren_part1_paper_candidate_comparison_jets.py` now uses the completed 124-order candidate tensor directly in the Section 9.23 comparison ODE. It supplies center Z = 0.3 Taylor coefficients for the comparison field, six radial integrals, pressure, D = I_theta/F, and I_z. The initial amplitude is cancelled algebraically only in normalized ratios; its logarithmic derivative and squared coupling remain present. No amplitude fit or pressure replacement occurs.

The comparison agrees with the actual core before h_b = 0.005. From h_b to 2 h_b it integrates phi_y = alpha s phi_core,s / phi_core * phi and U_y = alpha s U_core,s, together with all radial moment increments. Beyond 2 h_b it propagates the frozen fields and exact polynomial moment increments. The generated receipt includes the frozen comparison at physical R = 110. This is an auxiliary comparison field, not the actual exit bridge or completed transition field at R = 110.

All three inlet D coefficients are contained in the independently generated analytic core inlet enclosures. At y = 0.01, the finite 32-step calculation gives phi approximately 0.2793810980 and D approximately 3.4388246955. The change in D after refining to 64 steps is approximately 6.75844e-12. This is a numerical refinement observation, not an error certificate.

The axial driver derivatives are Taylor coefficients of the same computation; no off-center Z queries occur. `exit_driver_jets(y, Z, chi)` supplies A, B, A_Z and B_Z with chi treated as a function only of radial log coordinate y.

Scope is deliberately explicit: intervals cover finite input arithmetic and rounding, not RK4 discretization error or infinite radial remainders. The tensor is used only at its actual center; other axial centers fail. Whole-axis transition, terminal five-moment repair, admissible stress, heat matching and temporal recursion remain incomplete.

Next: consume these driver coefficients in the exit bridge and switches, propagate the actual transition moments, then attach analytic radial tails and a controlled ODE remainder. The existing old Lambda1e36 build path cannot be used as candidate evidence.

## Analytic driver hooks and actual initial exit

`ExitTangents(..., analytic_driver_provider=P)` accepts `P.exit_driver(y,Z)` returning scalar A, B, A_Z, B_Z. `ExitSwitches(..., analytic_driver_provider=P)` accepts `P.switch_driver(x,Z,blend)` returning D, D_Z, I_z, I_z_Z, Fbar, Fbar_Z, F0, F0_Z and optional bar/barE. Both bypass axial stencils. Focused hook checks exercise the downstream state/tangent RHS and preserve the old fallback.

`lei_ren_part1_paper_candidate_exit_analytic.py` connects the actual candidate tensor and comparison to the existing initial exit integrator. It uses MP midpoint projections of directed inputs and retains nonzero amplitude, actual pressure primitive and all five moments. Its positive epsilon comes from the conservative uniform analytic upper bound on |G|. The experiment compares 8 and 16 exit steps with 32 comparison steps. It is a center numerical experiment, not a proof of cone or moment closure. The initial exit ends at y = 0.01; the full actual R = 100..110 transition still requires a candidate continuation/switch adapter.

The completed tensor accepts its center up to MP rounding of 1e-255 only; this permits passing a 260-digit representation of 0.3 and does not enable an off-center finite-difference stencil or a whole-axis atlas.
