# Latest completed radial core receipt

The current production state is degree124/124 for the fresh center family [0.49,0.51]. Its normalized mixed C3 Taylor-tail upper bound is 0.000000000000746728625300476742790181607464024148938587047970157512882532601774676557607048563664246015 and the target1e-12 passes. The normal comparison loader was exercised successfully against the same completed state hash. The full core inlet moment and physical pressure enclosures were regenerated at124; Phi value error at scaledR4 is 3.49220244372785849560678239641592394262774874981271772947258669333289939957461940145348396e-86. Full Phi lies in approximately[0.046903439,0.517711916] and D in[1.07178345,30.68073765]. These broad intervals are not reduced by the tiny radial truncation tail; axial family dependence must be treated separately. All claims remain relative to accepted construction data. No comparison ODE error, terminal functional closure or temporal recursion is certified.

The degree48 discussion below is an earlier recorded checkpoint, retained to explain the error-propagation implementation. The next production action is controlled comparison/exit propagation and axial partition refinement, rather than adding radial orders beyond124 without evidence.

# Analytic core tail propagated into fresh interval inlet

The center family Z in [0.49,0.51] has reached radial degree 48 of 124. The maximum normalized mixed C3 analytic tail bound is 6279526932563632723608550706467788319003.73337882298061274505924125439806528849974796204879; the 1e-12 target remains unmet. Exact source/state identities and pressure data are unchanged.

The new interval_core_inlet_bound module encloses the infinite analytic candidate's core moments at a fixed scaled radius using the accepted-data radial tail bounds. At scaled radius 4, the normalized Phi value error is at most 4.49373081115103962191250082671682004964781281284209235387009640205411625803595139070928125e-33 and the full Phi enclosure has lower endpoint 0.0469034390076638056724425470336681057587832308863214921339461877720222514063701947147145192. The denominator remains positive, so the module also returns bounded core-inlet D and physical I_z. The still-large higher axial derivative tails explain why this value-level success does not pass the full C3 gate.

For each axial order k up to 3, the derivative tail bound is divided by k! to form a symmetric Taylor error coefficient. The physical U correction uses Psi error divided by Lambda. Finite radial coefficients provide uniform absolute derivative envelopes on [0,4]. Products then bound the errors in Phi squared, Phi U and U squared by the Leibniz convolution. Integrating these uniform error jets with the actual positive radial weights supplies theta, z, theta_z, p, u_squared and weighted_phi_squared errors. Physical pressure is restored from the same P0 plus S times the pressure moment divided by Lambda, preserving S=F0 squared.

This error propagation is for the core inlet only. It does not advance the analytic errors through the comparison ODE, certify original construction parameters or prove a stress cone. The comparison diagnostic receipt remains finite-core-only even though a separate core-inlet enclosure is now available.

The independent polynomial fixture checks 32 axial coefficients against closed-form integrals of known full linear radial fields with nonzero axial error derivatives. It checks the physical U error conversion and six weighted moments, without using the production integral helper for expected values. This fixture supports the error arithmetic, not a certificate for the complete Navier-Stokes construction.

## Next tasks

1. Continue the frozen interval production core to 124 and pass the mixed C3 target.
2. Maintain separate receipts for finite-core comparison arithmetic and infinite-core inlet bounds.
3. Propagate bounded core uncertainty through comparison/exit ODEs with a discretization error enclosure; interval RK arithmetic alone is insufficient.
4. Narrow axial families where coefficient dependency makes inlet ratios unnecessarily wide.
5. Add pressure-row and fresh scalar-center consistency evidence before connecting the new driver to transition and functional moment repair.
6. Complete whole-axis matching, heat exterior, final admissible stress, temporal coefficient recursion and oscillatory corrections before full Cartesian residual validation.