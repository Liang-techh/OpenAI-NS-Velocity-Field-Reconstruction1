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
