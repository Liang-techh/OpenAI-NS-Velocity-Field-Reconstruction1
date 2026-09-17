# Executable structure identities (CR010)

Run `python -m openai_ns_reconstruction.constrained_structure` after installing
`pip install -e ".[symbolic]"`. Result: `artifacts/constrained/structure_identities.json`.
The actual run used SymPy 1.14.0 and verified five exact symbolic identities.

For rho=x²+y², the candidate is curl(-y H,x H,0)+(-y S,x S,0).
The implemented H is amplitude*tau^(-A)*Z*b*q; S is the implemented swirl
factor. The v2 swirl polynomial only changes S. These functions depend on
rho,z,t, so the following results apply to both implemented families:

1. The poloidal component has zero divergence for H of class C2, by
   commutation of mixed partial derivatives. The script differentiates the
   vector potential, then its curl, rather than substituting zero divergence.
2. The swirl component has zero divergence for S of class C1. The script
   differentiates -y*S(rho,z) and x*S(rho,z) in Cartesian coordinates.
3. Every C1 axisymmetric pressure has -y*p_x+x*p_y=0. Away from the axis
   this is r times the azimuthal pressure gradient. Thus changing pressure
   cannot remove the measured azimuthal residual for fixed velocity/force.
4. The radial/swirl/axial representation commutes with every planar rotation.
   The script verifies preservation of rho and the vector transformation
   algebra for an arbitrary real angle and arbitrary scalar coefficients.
5. At fixed scaled coordinates, lr/lz*tau^(-A)=tau^(-1/2), using exact rational
   A=101/200 and D=99/200. Shape parameters and normalization must be constant
   in time for this statement. Future time-dependent extensions need a new check.

The global application assumes the smooth bump extension described in the
representation document. That extension's smoothness is not established by
these symbolic computations. The checks prove ansatz identities, not that
floating-point evaluation is exact, that NS momentum holds, that finite
sampling is a uniform error bound, or that a singularity exists. They are
symbolic checks, not a Lean kernel proof. The failed PDE results remain valid.
