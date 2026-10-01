# Resumable candidate core and analytic coefficient identification

The Lambda=1e48 candidate now has a resumable coupled radial recurrence. Its fixed target is radial degree124, initial axial Taylor length128, precision260, and axial center Z=0.3. Physical pressure is taken from the refined accepted-source receipt. The angular amplitude is factored as F=F0*A, while all F0-squared swirl and pressure couplings are retained.

Each completed radial order saves exact directed interval endpoint tuples atomically. Resuming checks the settings, accepted source, pressure data and nested dependencies, and source/driver hashes. Alternative pressure inputs receive independent state files. States created before the final loader changes were preserved locally under obsolete names; their old coefficients were not relabeled with new hashes.

The degree4 fixtures compare the reconstructed F=F0*A, Uz and pressure against the unfactored original recurrence, and A against the existing gauge recurrence. They also compare uninterrupted computation with JSON serialization and resume. These checks pass. They check algebra and storage, not matching or temporal recursion.

## Identification with the analytic candidate

The symbolic identity receipt verifies the complete angular and axial equations after substituting

    s=Lambda*R, epsilon=1/Lambda,
    F=F0*Phi, Uz=U0+epsilon*Psi,
    P=P0+epsilon*Pcal,
    F0_Z/F0=-Lambda*g, g=L*H0/(H0^2+sigma^2).

Here M=(integral_0^s Psi)/s and Pcal=F0^2*integral_0^s Phi^2. All10 angular and10 axial terms survive. The equations become

    2(s*Phi_ss+2*Phi_s)+chi*Phi=epsilon*Etheta,
    2(s*Psi_ss+Psi_s)=B+epsilon*Ez,

with chi=g*H0/L and B the exact axis axial equation divided by L. Both symbolic differences are exactly zero. The majorant uses 1+delta as a safe upper bound on |1-delta| where the exact angular transport contains 1-delta; this does not change the recurrence.

Given the analytic solution supplied by the accepted-datum contraction, its radial Taylor coefficients are uniquely determined by the axis data: order n determines the next angular and axial coefficients through nonzero divisors 2L(n+1)(n+2) and 2L(n+1)^2. Pressure coefficients follow by radial integration. L is zero-free on the previously admitted analytic tube. Thus the directed coupled recurrence encloses those exact Taylor coefficients, relative to the accepted stored source. This argument does not supply an independent derivation of the original construction parameters.

## Current result and next action

The refined-input state has reached radial order34 of124. Pure recurrence arithmetic took about177 seconds; file saving adds wall time. At this partial order, the largest finite coefficient uncertainty over 0<=s<=4.1 among mixed derivatives of total order at most3 is approximately5.52e-62 for Phi and4.34e-15 for Psi. Both maxima occur in the third axial derivative. The remaining infinite Taylor tail is still far above the1e-12 internal target at order34, so the total mixed derivative target has not passed.

The combined receipt covers only the axial Taylor center and the whole scaled radial interval. Continue the same immutable state to order124, recompute the combined error budget, then use the verified core data to regenerate exit traces and shared collar/moment data. Whole-axis functional matching, physical-coordinate error conversion, stress-cone closure and the n-dependent temporal recursion remain open.

Run one bounded segment:

    python experiments/root_st073/lei_ren_part1_paper_candidate_gauge_core.py --pressure-file lei_ren_part1_paper_candidate_pressure_axis_jets_refined.json --seconds 180

Recompute identities and combined bounds:

    python experiments/root_st073/lei_ren_part1_paper_gauge_fixed_point_identity.py
    python experiments/root_st073/lei_ren_part1_paper_candidate_combined_core_budget.py
