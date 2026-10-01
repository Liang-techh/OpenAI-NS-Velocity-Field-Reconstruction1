# Same-pressure unnormalized collar endpoint through second width — 2026-09-30

The analytic normalized collar coefficients are now converted back to unnormalized profile fields and five radial moments. This precedes the similarity-to-Cartesian space/time transformation; it is not a new physical time evolution or temporal recursion.

Let W be the formal width variable with unit coefficient. Its coefficient of power w is the physical coefficient divided by h_b^w. Keep Fa=F(Ra,Z), Ra, and the original inlet pressure P0 common to every field. Normalized moment coefficients transform as:

| Profile moment | Normalized-state multiplier |
| --- | --- |
| theta | Fa Ra^2 theta |
| z | Ra mz |
| theta_z | Fa Ra^2 mixed |
| z_theta | Ra axial - Fa^2 Ra^2 swirl |
| p | Fa^2 Ra p |

This conversion is applied independently at both retained width orders over the directed second-Z atom ring. At the endpoint, R=Ra(1+2W+2W^2), F=Fa[1+Wg1+W^2(g2+g1^2/2)], Uz=Uz0+Wu1+W^2u2, and each moment is its same-core inlet value plus its converted first/second-width terms. Pressure is the original inlet P plus the p-moment increments; the analytic pressure datum is never reset or replaced.

Radial velocity is recovered from the same axial-flux primitive:

Ur=[2ZR Uz-(1-delta)Z mz-(1-Z^2)mz_Z]/[(1-delta Z^2)sqrt(2R)].

The first Z derivative is obtained by explicitly differentiating this formula, using the second derivative of mz. Ur_ZZ is left unavailable because it needs a third moment derivative. Six mixed width/Z derivatives of Ur and Ur_Z were checked independently using MP differentiation of a resolved exponential-radius model, including width powers zero through two. The pressure-datum object identity and unavailable Ur_ZZ are preserved.

`lei_ren_part1_paper_collar_width_physical_endpoint.py` exposes the conversion. Its fixture is a resolved mathematical check, not an actual-source enclosure. The accepted finite-core interval inlet adapter is being connected separately to obtain actual coefficients.

The endpoint polynomial does not establish full spatial collar continuity, full-field divergence, omitted-width-order bounds, original-source errors, or global five-moment terminal identities. These require the spatial and error propagation layers; no completed background or scale-recursion claim follows from this endpoint adapter.
