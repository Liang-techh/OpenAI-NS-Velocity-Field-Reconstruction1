# Component-valued shear switches from R=100 to R=110

The short Section 9.4 shear switches now receive the derivative-aware prescribed continuation at R=100. The auxiliary drivers, pressure datum, endpoint velocities and all five moments use one common source. The axial and swirl quadratic integrals are carried separately from the collar through continuation, rather than reconstructed by treating their difference as either integral.

The first switch turns off the axial driver while retaining a=epsilon*D. The second changes a smoothly to 4/5 with zero axial driver. Both short switches use s=x/h_b, x=log(R/100), and retain h_b and epsilon as width atoms. The normalized eight-state ODE uses AxialDual arithmetic, differentiating its actual initial data and drivers without production finite Z differences. The final constant-a=4/5 segment to R=110 uses analytic power-law moment primitives.

Axis pressure P0 is inherited directly from the common component core. Total pressure is reconstructed from P0 and the cumulative pressure moment; total-minus-integral subtraction is not used to recover P0. Radial velocity and stress use the five moments and their automatically differentiated axial data. No targets are substituted for those moments.

The actual pressure9/width2 degree18 source completes R=100..110. Boundary field/moment scaled differences are below 1e-200; the leading F110/F100 ratio is 1.1^(-0.4), approximately 0.962593502656. Uz_Z retains a nonzero complete-preheat pressure-tail atom. At R=110, g_y=-0.4, Uz_y=0, and the raw axial-minus-swirl identity equals the corresponding cumulative moment at serialized working precision.

The independent resolved scalar physical-moment replay passes: its maximum scaled discrepancy is approximately 9.61e-13 at h_b=1e-6 and 1.02e-9 at h_b=1e-5. Those comparisons include finite ring truncation and scalar/switch RK errors; they are not remainder enclosures. The v2 Step 3 formulas and subsequent reshape inputs were checked directly in the supplied PDF, including equations (9.28)--(9.31).

At the actual Z=.3 probe, the leading D(110) is approximately 5.69e37 (>3), V-4Z is approximately 1.00e-14, and V_Z-4 approximately 1.91e-23. The reshape input B=log(Cstar*Utheta110*(1+Z^2)) is approximately -8.62e36, with B_Z approximately -8.05e35; both fit the local 2A budget for A=1e150. This is a local value/first-Z check, not the required uniform C2 source or stress-cone bound.

## Reproduction and limitations

- `python experiments/root_st073/lei_ren_part1_paper_pressure_width_switches_fixture.py`: independent scalar physical-moment ODE replay from the same R=100 input at resolved width.
- `python experiments/root_st073/lei_ren_part1_paper_pressure_width_switches_check.py`: actual degree18 pressure9/width2 source, including R=100 boundary continuity, R=110 fields/moments, pressure-tail retention and the constant-power slope.

The actual check writes an ignored local R100 cache before importing the switches. `--resume` reuses that same Z=.3 source for bounded switch repairs without recalculating the collar. The default command reconstructs the source and replaces the cache; use it again after changing upstream source inputs or algorithms.

The pressure/width/core Taylor and switch/collar RK errors are not enclosed. Automatic first-Z differentiation of a finite RK state is not independent radial/Cartesian divergence verification. Long reshape, functional terminal moment repair, finite global energy, exact heat, admissible stress and temporal coefficient recursion remain open.
