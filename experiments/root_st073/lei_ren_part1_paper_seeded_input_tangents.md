# Actual seeded input derivative transport

Run `python experiments/root_st073/lei_ren_part1_paper_seeded_input_tangents.py`.

The provider reads actual corrected inner terminal moments and momentsZ from inner.evaluate_x(e,Z), differentiates the reference target using terminal F/FZ, then transports raw offset derivatives through the existing Rp/Ep normalization. It preserves unresolved inner entries [3,5] and inherited subtraction/quadrature uncertainty.

Before Rp the source has Uz proportional to Z and Utheta proportional to 1/(1+Z^2). Consequently logEp_Z=-2Z/(1+Z^2). Incoming row derivatives use these source factors directly; the original reference is recovered from dimensionless integrals, rather than subtracting a large seed from a rounded seeded sum. Reference swirl energy divided by Ep^2 is Z independent, so its derivative cancels analytically. At Z=0 direct reference linear factors are required; the actual factory obtains them from a nonzero reference quadrature, avoiding division by zero.

The nominal future angular-energy derivative remains a fourth-order float-backed Z difference at steps 1e-4 and 5e-5. It is explicitly sourced, not set to zero. At Z=.3 its log magnitude is about -83.15; refinement change has log magnitude -109.73. This is not an error enclosure. Actual coefficient derivatives replay the linear equations around 1e-201 under the supplied nominal inputs.

Next differentiate angular correction coefficients and heat defects analytically, replace float-backed caches, and provide shared continuous incoming primitives. Terminal mean/Z-mean closure, global finite energy and recursive/stress/oscillatory closure are not certified by this receipt.
