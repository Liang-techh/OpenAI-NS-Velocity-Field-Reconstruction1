# Common-source reference background with retained five defects

ReferenceDefectBackground evaluates the terminal reference interval x=R/Rm in [1,e], equivalently axial-restoration phase 2+log(x). It consumes CenteredComponentDefects from the same R110 source and preserves its explicit P0/P0_Z. Whole defect evaluations are cached per Z; the default delta is inherited from the canonical source rather than reparsed at a different precision.

The exact homogeneous reference moments are

    Mz = 4ZR
    Mtheta = (5/8)*sqrt(2)*R^(3/2)*u
    Mtheta_z = 4Z*Mtheta
    Mztheta = 16Z^2*R-(5/12)*R*u^2
    Mp = (5/2)*u^2

where u=Am*x^(1/10). The actual input/reference velocities agree throughout this interval, so their five original moment differences stay constant in R. The normalized centered defect parts determine these constants through Eq.(10.3). Every source part is separately converted into physical moment increments and first-Z data; reference_power is also retained in the part maps. Public sums may round away tiny contributions and do not replace the authoritative parts.

This reconstructs the supplied source integrals; it does not set moments to terminal targets or reset the pressure datum. Pressure is P0+Mp. F, Uz and their derivative/stress inputs follow the same reference velocity, and Ur follows the same cumulative axial moment and first-Z derivative. Individual raw axial/swirl energies cannot be inferred from d4 alone and remain unavailable in this provider.

## Resolved independent evidence

The source-lineage fixture checks all five radial moment densities against independent finite differences, with maximum relative error 3.20e-30. Defect parts remain constant over R, reference/source part sums recover every moment, P0 is preserved, and the whole source is evaluated only once for a fixed Z. This fixture uses an explicit resolved surrogate source and does not certify the actual core.

The actual check uses the authoritative local R100 cache through build_provider(resume=True), then recreates the live switches and complete centered defect calculation before joining a finite degree 3 bump response. That cache is restricted to Z=.3. It is not a uniform axial profile or a global source certificate.

## Remaining requirements

A finite same-source reference/bump join does not establish five-moment functional closure. The finite inverse response has a remainder; tiny defect responses must remain separate. Uniform axial smallness/C2 bounds, source/jet/quadrature errors, the relaxed cone throughout the supports, full outer energy/heat matching, genuine n-dependent temporal recursion and oscillatory corrections remain open.

## Actual common-source join

The actual cached-R100 run completed at Z=.3 and x=1.25, with pressure order 9 and width order 2, using a fresh complete centered defect evaluation. P0/P0_Z are unchanged, the shared stress recovery executes, and the joined part map contains 69 source labels, reference_power and five_bump_correction. The actual background and corrected field/moment receipts are saved in five_moment_reference_background_check.json. An internally generated ignored reference snapshot now lets subsequent diagnostics reuse this exact source state without recomputing the collar and flat defects.

This confirms the finite same-source join at one axial/radial point. It does not verify uniform moment closure, cone admissibility, individual raw energies, full outer matching or temporal recursion. Rounded public fields/moments can still hide tiny individual defect responses; the separate source and correction records remain necessary.
