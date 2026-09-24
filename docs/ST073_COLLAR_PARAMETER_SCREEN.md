# ST073 collar parameter screen

At `tau = 0.5/64` and `r = 0.0056890761915166545`, the current compact field passes the local wave-cone test at only one of five sampled axial positions (`z = 0.0025, 0.00275, 0.003, 0.00325, 0.003432627453438968`). This is a local diagnostic, not a spacetime certificate.

The reproducible screen in `experiments/root_st073/curl_wave_cone_parameter_screen.py` decomposes the current field into its base plus three existing correction groups, precomputes each group's kinematics, and recombines the full quadratic momentum residual for 63 amplitude combinations. The best combinations pass two of five points. The earliest improvement scales the incremental collar by 1.25 while leaving the other groups unchanged; the maximum residual over the five center points changes from 688,423 to 690,264. Increasing that collar scale further through 3 does not increase the pass count. The radial-swirl group has zero effect on this inner-radius sample because its support is farther out.

The three lower sampled heights remain outside the local cone for all screened combinations. Thus adjusting these existing group amplitudes does not provide a sufficiently wide axial neighborhood for the proposed wave packet. The next construction step is to change the spatial shape or placement of the inner/axial transition correction, or jointly redesign the wave envelope and its time scale. A pointwise cone pass alone is insufficient: the resulting wave must also control exact-curl remainders, viscosity, and the full spacetime momentum residual.

Machine-readable results: `experiments/root_st073/compact_potential/curl_wave_cone_parameter_screen.json`. This scan is sampled at one time and radius; it makes no global acceptance claim.
