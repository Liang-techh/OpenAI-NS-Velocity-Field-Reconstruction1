# Shared candidate core moments and analytic inlet enclosures

The refined Lambda=1e48 degree124 core now provides all five core radial moment contributions and a common inlet dataset at scaled_R=4, Z=0.3. Every quantity is derived from the same accepted pressure anchor and coupled coefficients. This supplies the core contribution to future matching; it does not close the five terminal identities.

With s=Lambda R and epsilon=1/Lambda, the normalized moment functions are:

| Paper moment | Physical scale | Normalized function |
| --- | --- | --- |
| theta | F0 epsilon² | 2 integral_0^s sigma Phi d sigma |
| z | epsilon | integral_0^s Uz d sigma |
| theta_z | F0 epsilon² | 2 integral_0^s sigma Phi Uz d sigma |
| z_theta | epsilon | integral_0^s Uz² d sigma - epsilon F0² integral_0^s sigma Phi² d sigma |
| p | F0² epsilon | integral_0^s Phi² d sigma |

The receipt stores ordinary axial Taylor coefficients of these normalized functions through order3. Amplitude-dependent derivatives are restored algebraically where needed: for example, Mtheta_Z/(F0 epsilon²)=mtheta_Z+ell*mtheta, where ell=F0_Z/F0. This cancels the common amplitude before computing I_theta/F, rather than dividing two independently widened physical amplitude intervals.

Pressure is recovered as P=P0+epsilon F0² integral_0^s Phi². It is the prescribed radial primitive, so P_R=F² for the retained angular polynomial. This completes its pressure consistently beyond the stored radial pressure prefix; it does not add a fitted pressure term or change P0.

## Infinite radial remainder transfer

The combined core error receipt identifies the exact analytic candidate's Taylor coefficients and bounds the omitted radial orders. Its uniform radial remainder bounds are integrated into the moments. For a product, the kth derivative error is bounded by the Leibniz sum of E_f B_g+B_f E_g+E_f E_g, where B bounds the finite field and E bounds its omitted radial tail. The tiny F0² terms remain included.

The resulting analytic moment enclosures are transferred through the same original inlet formulas. The receipt separates the finite field from the analytic core enclosures. It reports normalized radial velocity, recovered pressure, I_theta/F0, I_theta/F, I_z/sqrt(epsilon), and their available axial jets.

Both total core stress components and their axial derivatives through order2 contain the exact zero required by the analytic core equation. In normalized units, the value intervals are approximately:

- Ttheta/F0: [-4.970e-64,4.970e-64].
- Tz/sqrt(epsilon): [-1.502e-16,1.502e-16].

Zero containment is a consistency check, not a terminal matching certificate. Six independent polynomial integral checks and five comparisons with the established physical stress evaluator pass. The independent five-moment angular trace and the direct angular-defect integral also agree through overlapping directed enclosures.

## Remaining integration

The receipt covers the axial center only. Whole-axis radial derivative sign, functional terminal five-moment closure, shared collar regeneration, full stress-cone margins and temporal coefficient recursion remain open. The next collar implementation must consume this source and its error receipt together; it must not mix the old Lambda=1e36 inlet with the new candidate.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_candidate_shared_inlet.py
