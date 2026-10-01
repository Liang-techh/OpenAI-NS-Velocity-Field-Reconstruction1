# Same-source pressure refinement for the candidate core

The previous beta=2 normalized true-mass enclosure had width about 4.44e-13, almost entirely from slope_transition_ref. This is a genuine limitation of the available enclosure, even though the source itself is unchanged. At Z=0.3 the pressure contribution to the first axial radial coefficient is

    U1_pressure = [d P0' - 2(1+delta) Z P0]/(2 L).

For P0=-exp(28) M2 (1+Z^2)^(-2), the common M2 dependence can be retained exactly. Its midpoint uncertainty alone gives a normalized Psi value error budget about 0.8867917147661 at scaled_R=4.1. This is an upper uncertainty budget, not a lower bound on the actual field error. It shows that a 1e-12 infinite Taylor tail does not by itself certify total computed-field accuracy.

An independent directed integration of slope_transition_ref with 128 panels and Taylor order24 now narrows that mass enclosure. Its intersection with the old enclosure is nonempty; all other accepted pressure stages are preserved. The first-coefficient source value budget becomes approximately 1.08859957801794e-16. This covers that contribution only, not all orders or mixed derivatives.

The refined physical axis-pressure receipt contains orders0 through160, uses the same accepted schedule and the same physical pressure scale, and preserves the flatten Cauchy enclosure. All161 refined physical coefficient intervals were checked to be contained in their original intervals. No old source receipt is overwritten. Original parameter derivation errors remain outside this claim.

The resumable candidate gauge recurrence must use a separate state when switching to these refined input intervals. A coarse-input state cannot simply retain its old rows and acquire the new input hash. After generating the refined coupled rows, add the finite coefficient error bounds to the analytic infinite Taylor tails before asserting any field error target.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_candidate_pressure_mass_refinement.py
    python experiments/root_st073/lei_ren_part1_paper_candidate_refined_pressure_jets.py

The independent scaled-coordinate evaluator now also bounds finite coefficient uncertainty over 0<=scaled_R<=a at the axial Taylor center. It keeps the axis U0 subtraction algebraic in Psi=Lambda(U-U0). Sixteen known polynomial derivative checks and three known coefficient-error checks pass. This remains a local axial-center calculation; functional matching and temporal recursion are open.
