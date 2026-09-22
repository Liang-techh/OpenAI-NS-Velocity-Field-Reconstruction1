# ST064 — first scale-native prototype and a negative dynamical result

Task #1084. User-authorized implementation of the scale-recursion priority on #939. This round implements the clock, a physical field, independent derivative/volume checks, actual material winding, and MATLAB visualization. It does NOT complete the scale-recursive NS mechanism. No default, old candidate, publication or schedule changed.

## What actually runs

ST064-S uses the actual ST063-G2R spatial snapshot at t_ref=0.5, not its time polynomial extrapolated beyond0.75. q=tau/tau0=2^-k, tau0=.5,T=1,h=.005,k in[0,6]; ar=q^.5,az=q^(.5-h). The new physical factors are A=q^-1 A0(s/ar^2,z/az), B=q^(-1-h) B0, C=q^(-.5-h) C0, p=q^(-1-2h)p0. The poloidal representation preserves divergence. All time-dependent-coordinate derivatives and physical convection, pressure, viscosity and force are included.

The exact powers and rigid spatial profile are autonomous assumptions, not an NS scaling symmetry. Normalized profile collapse is imposed by construction, NOT evidence of dynamically generated recurrence. The source's coupled similarity core, annular momentum-flux correction and matched exterior are not implemented.

Viscosity remains.01 and the old independently specified force parameters [.10145917594737912,.15152417019028727] are fixed in physical coordinates. No force=R substitution. Velocity and pressure support shrink inside the old cylinder; force outside that shrinking support is retained in the full residual integral. The NEW initial time is.5 and energy is normalized ONCE there by mu=.9898451568522917, not at every frame and not under the old E(.25)=1 convention.

## Measured finite-ladder results

|k|time|core radial RMS|core axial RMS|axial/single-transverse aspect|N_nu|full spatial L2|scale-resolved sampled max|
|---|---|---|---|---|---|---|---|
|0|.5|.28252364|.34590264|1.73146643|.18403400|6.34112235|2.12398353|
|3|.9375|.09988719|.12357321|1.74956276|.18595742|29.73842602|45.83654046|
|6|.9921875|.03531546|.04414635|1.76784822|.18790095|143.14875080|1039.23726607|

Both physical lengths shrink. Radius becomes1/8, aspect and characteristic diffusion-time turns rise only2.10%. Moments use omega_z squared in a fixed scale-aligned reference cylinder r<=.35,|z|<=.6; they are observation-window-dependent, not detected intrinsic boundaries. Energy decreases1->.13186868 while characteristic velocity rises. This is not amplitude collapse or per-time normalization.

ALL checked physical momentum results fail the original absolute.001 targets by large margins. This is a rejected rigid-template NS ansatz, not a better candidate than ST063. Twelve independent scale checks comprise seven boundaries and five new interior log-times. Each includes40/64/96quadrature, separate space/time Cartesian derivative ladders and a73x141scale-resolved spatial peak scan. These are NEW-window audits, not an execution of the old fixed-window ST063 validator or a continuum supremum certificate.

At k0, the new time-derivative L2 is6.10435 versus convection.67374 and viscosity.27639. Inheriting a spatial snapshot does not balance its new dilation derivative. At reference r=.1,z=0 the material radial contraction rate is only4.17%of the imposed support contraction rate. A shrinking Eulerian support is not by itself a material cascade.

At k6 only5of4096uniform-box points hit the shrinking velocity support: their L2 estimate110.59668 differs markedly from mapped full-domain143.14875. The latter includes the unshrunk-force exterior. Fixed-box sampling alone must not certify a concentrating field.

## Actual particle winding, not streamline counts

The same12initial particles were integrated in physical time through k6, with coordinate-drift terms correctly included in the equivalent similarity ODE. They average only.06364325total turns byk6, and0/12remain inside the shrinking core observation cylinder. Independent direct physical-coordinate integration and tighter tolerance checks pass.

Separately reseeding12particles at identical similarity coordinates each scale gives per-interval mean turns.00623299->.00634194. Those are different particles and the modest increase follows the imposed scaling. Neither record demonstrates many recursively accumulating turns in the same core. Diffusion-time characteristic turns and actual finite-interval turns are different measurements.

## Deliverables and tests

The complete offline package contains the actual profile, source snapshot/raw coefficients, scale parameters, physical derivative evaluator, all12scale audits, trajectory arrays, native MATLAB model/reference MATs, GUI,13focused tests and Python numerical PNG/GIF previews. Root MATLAB entry: start_scale. It has three views: fixed physical coordinates, ONE isotropic magnification, and separately labeled anisotropic similarity coordinates. Streamline seeds are scale-aligned, not material. The surface uses an explicitly normalized omega_z threshold. Physical speed color units are fixed.

Native MATLAB was unavailable: NO native GUI execution, MATLAB screenshot or cloud success is claimed. Python reference and MAT checks passed. The existing ns_basis.m is reused; a native test_scale_export(false/true) entry is included for numerical/GUI checks.

13localtests passed4.46s. Independent clean-directory13tests passed; the actual k6quadrature replay exited1 and its entire integral-result dictionary matched exactly. Source-snapshot velocity/pressure agree with the independent old evaluator within5.58e-14/1.88e-14 after the declared one-time normalization. The time-FD checker initially moved x instead of t and failed3tests; the correction and original failure are retained. Plot permissions and bounded rendering timeouts are retained; completed frames were checked and reused, not invented.

Executed physical evaluator committed at c3b5bd0ab0cd77da485c1789dcb5f367d969a83e; Gitblob fd7c35babcb8b5663ca7405fc28a4cdcf3d1c7da matches local bytes. This branch contains that evaluator and this result record, NOT every dependency/MAT/profile array. The COMPLETE files are in the conversation ZIP. No PR/merge/default promotion.

Source G2R rawSHA b9505453786af0ae0b0fd6e6d784d14c6892827b72cf2a70ace2f9f87acd13a0. New profileNPZ SHA e1e57ead8f055f0c8a4aba242c1270392457cbc59955f06104297ad88121300e.

Next construction must solve the similarity-coordinate momentum balances and core/exterior compatibility rather than compress the old snapshot further. Kinematic observables implemented; dynamic recursion and full NS target remain unachieved. All PDE/source-correspondence/blow-up flags false.
