# Continuous partial axial-square energy

This receipt evaluates the actual inner-seeded axial-square primitive through `R_v` with the live continuous pulse, end basis, coefficients, and coefficient tangent.

## Normalization

For `xi = mu log(R/R_p)` and `t = log(R/R_v)`, the pulse contribution is

```text
R_p E_p^2 * a^2/mu * integral_0^xi exp(-2 s) gp(s)^2 ds.
```

The end contribution is `R_p E_p^2` times the sum of the positive partial bump atoms. Each atom contains `exp(-26 - 2 mu center)` and `integral exp(-2 mu s) beta(s)^2 ds`. The pulse support is `0 <= xi <= 11`; end support is `-3.15 <= t <= -.85`, so `xi = 13 + mu t > 11` on the end band and the cross term is zero by support.

The incoming prefix is `R_ref I_uz2(y,1) Z^2` plus the measured inner axial offset, with its analytic Z derivative. The offset is applied once and is never rounded to zero.

The normalized Z jet uses `2 a a_Z` for the pulse and `2 sum_i c_i c_i_Z K_i_partial` for the end atoms. The physical scale contributes `scale_Z = 2 LZ scale`, with `LZ = -2 Z/(1+Z^2)`.

## Numerical checks

- quadrature order: `96` (positive MP Gauss--Legendre; enclosure: `False`)
- samples: `incoming_reference, pulse_bulk_xi_5, end_bump_1_center, end_bump_2_center, terminal_Rv`
- maximum Rv boundary relative difference: `0.0`
- maximum checked radial/Z integrand relative difference: `4.6950879000110920453107151394143090748650067023113e-17`

The boundary compares this partial provider with the existing complete provider at `R_v` for axial square, axial-square Z jet, `z_theta`, and `z_theta` Z jet. Forward checks compare fourth-order log-radius differences with `R U_z^2`, `2 R U_z U_{z,Z}`, `R(U_z^2-U_theta^2/2)`, and `R(2 U_z U_{z,Z}-U_theta U_{theta,Z})`. At the two end centers, axial square and axial-square Z use the separately retained normalized end atoms; direct whole-sum stencils are retained as precision-loss diagnostics because the full accumulated pulse baseline overwhelms the local end atom at the working precision.

## Limits

The result is a nominal numerical primitive. Quadrature error, inherited inner/incoming uncertainty, coefficient-input uncertainty, full five-moment closure, heat continuation, finite energy, and scale recursion remain uncertified. The complete provider should be called separately for `R >= R_v`; terminal residuals remain materialized.
