# Actual-candidate finite five-bump repair

The new Lambda120 R110 data now drives the Section 10.2 five-bump inverse. `experiments/root_st073/lei_ren_part1_paper_candidate_five_bump_inverse.py` builds the same-source centered defects, reference background and corrected field; it does not load the old source pickle or inverse receipts.

The correction solves A h + Q(h,h) = -d in paper coefficient order (c1, c2, xi1, xi2, xi3). Value and first axial derivative propagate together through AxialDual, including Am_Z. Each Picard increment is saved separately with exact MP atom tuples. The infinite inverse and the temporal n-dependent coefficient recursion are different constructions and remain unimplemented here.

Nominal materialized coefficients at Z = 0.3:

| Coefficient | Value |
| --- | --- |
| c1 | -9.9974682385e-15 |
| c2 | 7.7630972290e-15 |
| xi1 | -8.9574552416e-37 |
| xi2 | 2.0426723581e-36 |
| xi3 | -1.1341041624e-36 |

Ten nonlinear updates reduce the largest retained value/first-derivative polarization residual to approximately 6.23766e-210. This is finite-map algebra at 473 digits, using finite 64-order bump quadrature. It is not an integral-error certificate, a relative bound for every flat source component, a physical momentum residual, or a proof of terminal functional matching.

The corrected field is materialized at x = 1.25 using the actual same-source reference and cumulative moment corrections. The inherited P0 and P0_Z are preserved. Its stress evaluation succeeds, but no cone certificate is inferred from that fact. Very tiny response increments are retained independently even when a materialized coefficient sum cannot resolve them.

The endpoint source now keeps canonical delta as the decimal string 1e-200. This prevents a false mismatch when consumers compare MP values rounded in different precision contexts; the construction parameter is unchanged.

Next: independently integrate the actual physical bump moment differences and their axial tangents; quantify quadrature error; verify repaired stress over every bump support; then replace center source data by axial functions and prove a uniform inverse. Heat matching, infinite radial/source errors and temporal recursion remain open.

## Independent physical integral replay

Run `python experiments/root_st073/lei_ren_part1_paper_candidate_bump_integral_replay.py`. This reads and validates the persisted actual-candidate inverse controls, amplitude and defects, then independently integrates the physical polarized moment differences on the analytic reference profile over the bump supports. The fresh order-128 quadrature does not call the production map's apply method or reuse its product weights.

The largest retained value residual is 1.1390562734e-27; the largest first-Z residual is 2.5951210115e-50. These differ from the 6.23766e-210 finite-map algebra residual because the replay uses independent finite quadrature and rebuilds bump normalization at 180 digits. Neither number is a physical momentum residual or a rigorous quadrature error bound.

The helper's physical scaling is checked at Rm = 1 and Rm = 3. The physical theta and theta-z differences carry Rm squared; the swirl contribution to the z-theta difference also carries Rm squared. Testing only Rm = 1 would hide omitted radius factors. The fixture checks radius scaling and tangent propagation; its JSON is distinct from the actual-candidate replay receipt.

This completes an independent finite physical-moment consistency check at Z = 0.3. It does not replay the full source profile, bound its reconstruction errors, prove whole-axis terminal matching, or establish a uniform stress cone. Next: inspect the repaired stress across each bump support and develop controlled integral/functional inverse bounds.
