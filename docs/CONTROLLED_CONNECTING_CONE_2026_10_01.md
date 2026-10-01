# C120-F15: complete R110-to-Rm relaxed-cone certificate

The fresh fixed-parameter connecting field now certifies the relaxed `kappa<=2` cone over the complete profile-radius interval `R110<=R<=Rm`, for the entire accepted axial family `Z in [0.49,0.51]`. Two whole logarithmic cells cover long reshape/reference continuation to Rz; 64 gap-free cells cover axial restoration/reference continuation from Rz to Rm. The composed certificate checks the common core/family, dependency hashes and exact boundary `rz_log+2=rm_log`.

This resolves the long-reshape margin/direction non-certificates reported in C120-F14. Velocities, cumulative moments and pressure have not changed. It remains a relaxed-cone result; strong admissibility, whole-axis coverage and original parameter remainders are not established.

## Cancelling common angular amplitude

Let `Ctheta=R*sqrt(2R)*u`, `zeta=u_Z/u`, `d=1-Z^2`, `L=1-delta*Z^2`, and `J=Jtheta`. The transport is

`H=-1+(1-delta)*Z*Mz/R+d*Mz_Z/R`.

The source angular atoms at `R0=110` are divided by the shared source angular scale, analytically cancelling the source amplitude `F0`. In terms of the normalized source packet,

`n_theta=eps^2*theta0/(2*R0^2*phi0)`,

`n_theta_Z=eps^2*(theta0_Z+ell*theta0)/(2*R0^2*phi0)`,

and likewise for the angular/axial mixed moment. Here `ell=F0_Z/F0`; it must appear in the physical axial derivative atom. These are not derivatives of separately divided interval ratios.

The source contribution is multiplied by the exact ratio `A=exp(-1.6*y+B*sigma(y/T))`. With `V=4Z+g`,

`C0=1-delta/2-d*V_Z+(2delta-1)*Z*V`,

`Cz=(1-delta)*Z/2+d*V`,

the angular inertial stress is

`Itheta/F=(R/L)*[H+A*source_combo+J*(C0-Cz*zeta)-Cz*J_Z]`.

The repeated `zeta*J` terms are grouped before interval evaluation. Axial mass is similarly centred as `Mz/R=V+(Mz0/R0-V)*exp(-y)`, avoiding a spurious zero lower bound from two independently evaluated complementary weights.

On long reshape/reference continuation, axial shear is exactly zero. The prescribed angular shear is `st=-a=-0.8-2B*sigma'/T`. Therefore `kappa=a` and the relaxed cone margin is `Ttheta/F-(2-a)`. The axial inertial stress is irrelevant to this branch test and is neither replaced by zero nor fitted. Pressure retains the existing axis datum.

## Restoration and evidence

Restoration cells retain both actual shear components and use the general cone evaluator. Positive cumulative restoration primitives are enclosed by their saved full-support upper bounds; after support completion, the saved full intervals are used. Every physical radial factor and field is evaluated over the full cell. All 64 cells pass; no unresolved cell remains.

An independent scalar stress fixture checks angular inertial/total normalized ratios against `evaluate_mp_stress`. Common angular-amplitude scales `1`, `1e-25`, and `1e25` preserve both ratios. The synthetic fixture intentionally need not pass the cone. An additional restoration adapter regression contains 88 existing point-packet coefficients in their full cells; it supports translation consistency, while the analytic primitive bounds establish whole-cell coverage.

New companions and receipts under `experiments/root_st073/` use the `lei_ren_part1_paper_` prefix:

- `interval_long_reshape_angular_cone` and `_check`;
- `interval_axial_restore_cone_atlas` and `_check`;
- `interval_connecting_cone_certificate`.

Next recover the paper's post-five-bump shear/flatten construction and exact heat exterior, preserving terminal moment identities and analytic pressure. Higher axial smoothness/Ur_Z, whole-axis coverage, complete finite energy, flat remainder, genuine temporal recursion, oscillatory corrections and independent full Cartesian NS residual remain incomplete. The existing 40-cell reference-annulus certificate is separate and remains valid.
