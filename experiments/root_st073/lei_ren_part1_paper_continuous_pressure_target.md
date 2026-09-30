# Actual terminal-pressure restoration target

ContinuousPressureMoments.terminal_pressure_jet integrates the actual inner-anchored pressure through preheat, compact angular bumps and the complete infinite heat tail. It retains baseline, current bump, heat correction and their Z jets separately. required_bump is the integral required to enforce the paper normalization P(infinity,Z)=0; additional_bump_required is the change from the current integral. This extracts a target and does not change the pressure datum or any coefficient.

At Z=.3 the current nominal P(infinity) is -0.01208803828438 and P(infinity)_Z is +0.01330793205620. The current signed compact pressure increment has log magnitude about -2.719157e28, whereas the required total increment is about +0.012088. Those are materialized numerical values, not a proof of incompatible analytic targets: the propagated pressure subtracts totals about 4e12, and pressure-stage nodes were initially supplied by binary64 leggauss. Inherited inner and schedule precision must also be assessed.

Run continuous_pressure_target.py --mp-node-audit to compare same-order MP pressure-stage nodes against the current nodes, keeping inner anchor, bump, heat and angular schedule inputs fixed. This comparison isolates one numerical source; it is not a complete error enclosure, a new inner construction or a coefficient feasibility proof.

Primary construction in the locally cached Lei-Ren Part I text gives P(R,Z)=Mp(R,Z)-Mp(infinity,Z) and restores P0+integral F^2=0. A Z-dependent datum reset alone changes the momentum equation. Once numerical contamination and actual seeded defects are separated, restore the pressure and angular-mean constraints together, regenerate dependent axial inputs, and replay stress.

Same-order MP pressure-stage node audit completed: P(infinity)=-0.01229506811512 and its Z derivative=+0.01353585480563. The shifts from binary64 nodes are -0.0002070298307370 and +0.0002279227494353 (about 1.7 percent of the original residual). This confirms sensitivity to node representation but does not explain the remaining mismatch. Next compare MP quadrature orders and trace the inherited inner pressure datum using separately normalized components. Do not use the residual as an angular-bump target until numerical contamination and actual target mismatch are separated.

The independent read-only worker confirms the pressure term attribution. Its suggestion to subtract P(infinity,Z) directly is not installed: although that produces the paper exterior normalization, its Z derivative changes the axial momentum equation and the existing inner stress-free datum. A coherent global update must also reconstruct the inner solution.

The next angular stress-target diagnostic should use the full combination

```text
Btheta = (1-delta/2) Mtheta - (1-delta) Z Mtheta_Z / 2
         - (1-Z^2) Mtheta_z_Z + (2delta-1) Z Mtheta_z
transport = -R + (1-delta) Z Mz + (1-Z^2) Mz_Z
Btheta_target = -2 L R [Utheta*transport/(L sqrt(2R)) + S_theta]
```

This is an equation for coherent profile/moment restoration, not permission to overwrite the integral. The simplified normalized angular mean 1/(1-delta/2) applies only after flattening with G=1 and a constant heat reference slope; it must not be used at arbitrary finite flattening or actual corrected heat points. Preserve derivative identities and mixed moment terms.

MP pressure-stage order comparison completed on one shared field: orders 96, 128, 192. The 128-to-192 terminal changes are 1.5568102589e-15 for P and 1.7139195511e-15 for PZ. Both converge near P=-0.01229506638282 and PZ=0.01353585289852. This nominal convergence does not enclose quadrature error, but shows the measured 0.0123 mismatch is not explained by the tested pressure-stage order. The core is currently built from a legacy angular pressure anchor before the continuous angular schedule is installed; an opt-in route now reconstructs the core from a continuous preflattening integral instead. Post-Rv tail jets remain omitted and must be restored coherently later.
