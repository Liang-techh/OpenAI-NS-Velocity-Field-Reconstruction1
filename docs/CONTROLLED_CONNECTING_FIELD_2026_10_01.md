# C120-F14: callable connecting field from R110 to Rm

The local fixed-parameter C1 background now has callable long swirl reshape, reference continuation, axial restoration and final reference continuation up to the five-bump annulus. The unified API is `IntervalAxialRestoreField.evaluate_log_offset(y)`, where `y=log(R/110)` and `0<=y<=log(Rm/110)`. It encloses the fresh axial family `Z in [0.49,0.51]`; scalar axial queries, whole-axis coverage and physical-time recursion are not installed by this API.

## Coordinates and shared data

The exact boundaries are

- `Tshape=4e152`;
- `rm_log=10*(5e151+14)-6`;
- `rz_log=rm_log-2`;
- axial restoration phase `t=log(R/Rz)` in `[0,1]`;
- final reference continuation phase `t in [1,2]`, ending at `Rm`.

These are profile-radius coordinates. The shaping interval parameter is not physical blow-up time. All layers retain the same accepted core/source, centered `g=Uz110-4Z`, shape amplitude `B` and its first axial derivative, five source moments and analytic axis pressure. Pressure remains `P0+Mp`; no fitted pressure datum is introduced.

## Controlled long reshape

With the increasing flat cutoff `sigma`,

`Utheta = exp(-logC-log(1+Z^2)+y/10+B*(1-sigma(y/Tshape)))`,

`Uz=4Z+g`, and `Utheta_y=Utheta*(0.1-B*sigma'/Tshape)`.

The source endpoint uses the actual R110 velocity packet. Completed shaping uses the exact reference amplitude. Common positive amplitude factors are cancelled analytically when retaining the prescribed normalized angular shear, `-0.8-2B*sigma'/Tshape`.

All five cumulative moments are built from the source moments and entire endpoint-normalized kernels

`J(k,m;y)=integral_0^y exp(-k*ell+m*B*(sigma(y/T)-sigma((y-ell)/T))) d ell`,

with `(k,m)=(1.6,1),(0.2,2),(1.2,2)` for angular, pressure and swirl-energy primitives. No finite quadrature window replaces these integrals.

For `B<=0`, positivity and the universal bound `0<=sigma'<=32` give

`E(k+32m|B|_max/T;y) <= J <= E(k;y)`,

where `E(k;y)=(1-exp(-ky))/k`. The first axial derivative is bounded by `32m|B_Z|_max/T * min(y^2/2,1/k^2)`. The derivative bound 32 follows from cutoff symmetry and, on the left half, `sigma<=exp(4-1/s^2)`, `(1-s)^-3<=8`, and monotonicity of `s^-3 exp(-1/s^2)` up to `s=1/2`. It does not rely on numerically locating the sharper derivative maximum.

For `y>=Tshape`, the exact identity `J=E(k;y)+exp(-k*(y-Tshape))*I_full` reuses the earlier whole flat-kernel defect enclosure and its axial derivative. Thus the reference continuation carries the entire inherited shaping defect, rather than resetting moments at the shaping endpoint.

## Axial restoration

The restoring field is `Uz=4Z+g*(1-sigma(t))`, with `Uz_y=-g*sigma'(t)`. Angular velocity keeps the reference power law. Three cutoff primitives enter axial mass, mixed angular/axial moment and axial-square moment. Partial primitives use directed cell integration over their full support. At and after the completed restoration, the frozen full weights from the five-defect receipt are reused with exact endpoint normalization.

All cumulative moments, pressure, radial velocity and radial derivative use these same fields and primitives. At the inlet, 22 field/moment value and first axial coefficients reduce exactly to the preceding layer. At restoration completion, `Uz=4Z` and `Uz_y=0` identically. Partial cutoff primitives are checked against independent scalar quadrature. The unified dispatcher preserves the terminal pressure packet.

## Current limitation and next work

Long-shape production samples at its inlet, one logarithmic unit and midpoint do not certify the stress margin/direction with current interval algebra. Its endpoint and the restoration inlet sample certify the relaxed cone; four restoration/reference samples also certify only the relaxed cone. These samples are not a whole-path cone proof. Wide common-amplitude and cumulative-moment dependencies must be cancelled analytically before claiming the full connecting region admissible.

Next retain endpoint-normalized inertial stress formulas and source contributions explicitly, resolve the whole long-shape/restoration cone, then implement the paper's shear modulation to final strong admissibility and flatten/heat collar/exact exterior. Higher axial smoothness, `Ur_Z`, original parameter remainders, whole-axis matching, full finite-energy verification, flat remainder, genuine n-dependent temporal recursion, oscillatory corrections and independent full Cartesian NS residual remain incomplete.

Sources and receipts are `interval_long_reshape_field`, `interval_long_reshape_field_check`, `interval_axial_restore_field`, and `interval_axial_restore_field_check` under `experiments/root_st073/`, with the `lei_ren_part1_paper_` prefix. Independent moderate-scale fixtures support implementation; the analytic bounds and retained source receipts support the production enclosures.
