# Analytic radial velocity jets for shared-field stress

velocity_radial_jets(profile,logR,Z) augments the existing shared values_with_jets output with Utheta_y and Uz_y, where y=logR. Utheta_y uses the installed angular logarithmic slope, including shared bump/heat point definitions. Uz_y uses the incoming cutoff derivative before Rp, the shared pulse product derivative on pulse support, the owned end-bump product derivative in its two compact supports, and zero after axial support.

For pulse Uz=A*a*gp(mu*(logR-logRp)), the jet is A*a*(slope*gp+mu*gp-prime). For end Uz=A*sum(c_i beta_i), the jet is A*(slope*sum(c_i beta_i)+sum(c_i beta_i-prime)). The live runtime and basis are reused; no amplitude, correction coefficient, pressure datum, moment or terminal mass is replaced.

The helper is ready to be consumed by Section 3 stress evaluation once the actual five-moment/P provider is accepted. It provides nominal numerical derivatives of the same definitions; inherited coefficient/input errors and full stress/energy/recursion remain open. Tiny heat corrections remain separately exposed by the underlying field, since full amplitude sums can round them away.

Shared-field check passed in incoming/pulse/end/heat: maximum local relative difference about 1.61e-16.
