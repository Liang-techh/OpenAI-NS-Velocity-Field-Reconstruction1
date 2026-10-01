# Actual-candidate repaired matching-region cone samples

The new Lambda120 candidate and its ten-update five-bump repair now have a same-source support-aware stress diagnostic at Z = 0.3. Run `python experiments/root_st073/lei_ren_part1_paper_candidate_repaired_cone_scan.py`.

All 39 samples pass the relaxed cone. Sampling includes the three compact supports (centers 1.25, 1.5, 1.75 and radius 0.025), their edges, eleven positions per support, the gaps and x = 1, 2. Prescribed P0 and its first axial derivative are unchanged at every sample.

Ratios S_theta/F, S_z/F, T_theta/F and T_z/F are formed in the retained pressure/width ring before nominal evaluation. The test uses the original stress and actual corrected cumulative moments, not substitute moment targets. For kappa <= 2 its normalized relaxed margin is -dot/(-S_theta/F) - (2-kappa), where dot = (T_theta*S_theta + T_z*S_z)/F^2. Direction requires dot < 0.

Kappa stays approximately 0.8. Therefore none of these samples satisfies the stronger kappa > 2 admissible condition. The relaxed connecting-region condition and the stronger final admissible-stress condition remain separately reported. This diagnostic establishes no uniform support/axial cone, finite Taylor remainder, quadrature enclosure, or final admissible reconstruction.

The next step is a directed integral backend for the actual exponential bump family, followed by a uniform repair/cone estimate and downstream outer-region construction. Nominal sampling should not be expanded indefinitely as a substitute for those bounds.
