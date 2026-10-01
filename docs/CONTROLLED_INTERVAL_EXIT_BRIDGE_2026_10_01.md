# Controlled initial actual exit bridge for the fresh axial family

The fresh degree124 core and controlled comparison now feed the actual Section9.25-9.26 exit equations over y in [0,0.01]. This is an initial exit-profile enclosure, not the frozen auxiliary comparison field and not a completed Navier-Stokes solution. The profile radial coordinate is R=s/Lambda; lab-space and time mapping remain a separate construction layer.

## Joint field and moment propagation

The comparison-cell provider supplies16 core cells and16 switched-comparison cells. Each packet encloses the complete radial cell, including partial-cell cumulative moments. It retains field order2 and driver order1 in the axial Taylor variable, over all centers in [0.49,0.51]. A separate endpoint packet supplies the terminal comparison field and drivers.

The actual exit integrates g=log(F/Fa), U and six normalized moments. Its equations are g'=A and U'=exp(g)B, with A=-chi D/2 and B=-chi sqrt(s/(2 Lambda)) Phi(4)/barPhi I_z. The cutoff multiplier is exactly chi=epsilon_exit+(1-epsilon_exit) sigma(1-y/hb). The positive terminal epsilon_exit uses the same conservative uniform analytic G bound as the nominal exit; it is different from the radial scaling epsilon_radial=1/Lambda. This policy is not a certificate of the Section9 contraction constants.

Within each cell, directed A ranges enclose g and exp(g), then directed exp(g)B ranges enclose U. Those full field ranges enclose all six moment integrands. The integration therefore encloses comparison uncertainty, analytic core tails and exit discretization together; no inputs are projected to midpoints.

At the endpoint, g lies approximately in[-0.04214379,-0.00125077] and U in[1.96000000000001,2.04000000000001]. The physical angular field F remains strictly positive. The tiny terminal prescribed shear is retained separately in the endpoint radial derivatives and late-shear increment, even where it cannot be resolved as an increment to a much larger field value.

The amplitude jet is recovered from the same F0 and ell. Physical moments use F0 epsilon_radial squared for angular/mixed moments, epsilon_radial for axial moments, and S=F0 squared for pressure/swirl quadratic terms. Physical pressure remains the supplied P0 plus the actual cumulative pressure moment. The radial velocity value is recovered from the axial moment and its first axial derivative, enforcing the original divergence-form relation.

## Evidence and limits

A read-only mathematical review confirmed the exit normalization, multiplier, successive cell integrations, moment scaling and radial velocity expression. An independent512-step RK diagnostic on a known polynomial source checks18 exit-field, moment and first-axial-derivative coefficients against the family enclosure. This numerical diagnostic supports the implementation; the cell-range integral argument supplies the enclosure.

The result retains first axial derivatives of F,U,P and physical moments. It does not enclose the axial derivative of Ur, which requires additional axial moment derivatives. Original construction-parameter errors, final stress-cone conditions, the physical post-collar exit through R=100, the R=100..110 switch, terminal functional moment repair and whole-axis coverage remain incomplete. No temporal n-dependent recursion or full NS residual certificate is claimed.

## Continue

1. Continue the actual exit from its enclosed endpoint to R=100, retaining the nonzero epsilon shear and bounding all field/moment deviations. Do not substitute the auxiliary frozen comparison for the actual exit.
2. Preserve interval errors through the R=100..110 switch and supply the actual analytic axial/radial derivatives for stress evaluation.
3. Treat broad axial-family dependence separately from tiny radial tails; refine the atlas or correlated algebra where cone and moment bounds are too wide.
4. Generate functional defects and moment repairs from these actual exit/transition data, retaining the same pressure/amplitude datum.
5. Complete heat exterior compatibility, final admissible stress, temporal coefficient recursion and oscillatory correction before Cartesian residual validation.