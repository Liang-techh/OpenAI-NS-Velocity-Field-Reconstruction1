# Original steep power: physical stress and retained remainder — 2026-10-04

The original Ts-long steep-power companion now maps the same complete angular/energy/absolute-pressure moments into physical stress mixed3, the completed tensor diagonal, divergence mixed2 and axial-viscosity remainder mixed2. Exact exponential-mode cancellation is performed before enclosure. Original velocities, coefficients, pressure datum and source units are retained.

The focused producer/checker passes with 330 current input hashes, 156 finite signed physical rows and 40 exact zeros. An independent native Cartesian fixture exercises both positive viscosities .01 and .7, the completed tensor and two nonzero remainders. Its six Cartesian decomposition errors are below 3.4e-53; 24 divergence and 12 remainder mixed2 comparisons are below the declared 1e-36 tolerance. This moderate fixture does not certify global NS accuracy or reproduce the actual Gamma history; that source history is consumed separately from the accepted common receipts.

## Original source map and physical units

```
y = Ts*phase, ell = Ts*(1-phase), phase in [0,1]
q = log(R/Rtail) = -wait-1-ell
K = KQ*exp(k*ell), K_y = -k*K
delta = 2*a, k = 1-a, p = 1+delta
lambda^2-lambda^(2*delta)*z^2/nu = tau
R = r^2/(2*nu*lambda^2), Z = z/(sqrt(nu)*lambda^(1-delta))
u_phys = sqrt(nu)*u_source, p_phys = nu*p_source
Ttheta_theta = r*partial_z(Tz)
Etheta = -nu*partial_zz(utheta), Er = Ez = 0
```

The general-K mapper is the unchanged CompliantCollarPhysicalC2.collar. Ephemeral adapters supply the actual power shape, common stress packet and original phase dispatch. The generic mapper receives q for source log factors and shape derivatives; native assembly receives phase for the original radius branch. All K rows are derivatives with respect to y=q plus a constant; no Ts^j phase-derivative factors are introduced.

Stress uses the original nu*lambda^(-2-delta) factor. Divergence mixed r^i z^j uses nu^(1/2-(i+j)/2)*lambda^(-3-delta-i+j*(delta-1)), with its actual radial and Qtheta/Qz factors. The remainder uses nu^(1/2-(i+j)/2)*lambda^(-1-delta-i+(j+2)*(delta-1)) and the original B source factor. Positive source factors remain logarithmic; interval caps are not chosen as fields.

## Complete-moment cancellation

With L=1-delta*Z^2, d=1-Z^2, and the complete pressure defect PdQ at the steep-exit inlet,

```
Dtheta = k*K/L + 4*S*exp(-q)*K
NpQ = -Z-2*p*Z*PdQ+p*Z*KQ^2/3+d*partial_Z(PdQ)
Dz = 2*k*Z*K^2/(3*L)+exp(-p*ell)*NpQ/L
```

The full theta divergence source factor has rate -bh, so the two original modes have total rates -1.5 and -2.5. The full axial factor has rate -p; its pressure homogeneous mode cancels on differentiation. The resulting axial ordinary rows are Dz, -2*k*Z*K^2/L and 6*k*Z*K^2/L. Arbitrary complete endpoint pressure functions are used in the symbolic identity; no pressure datum is reset.

The actual divergence-function AST is replayed and compared with the full-moment derivative equations through mixed2. The general velocity operator gives the retained source remainder coefficient

```
2*(1-2*(1-delta)*Z^2-delta^2*Z^4)*K_y/L^3
    -4*Z^2*K_yy/L^2
```

It is generally nonzero in steep power. The waiting region's pure-radial-power zero remainder is not propagated into this region. Global flatness, physical-volume bounds and energy remain unfinished.

## Actual physical right join

The accepted power/exit similarity source joins provide K derivatives through4, stress mixed3 and pressure mixed4. The new physical bridge also binds:

- Actual thetaS/thetaQ and power/exit theta/rate expressions, including the actual flat sigma endpoint proof.
- Original packet to flatten_mixed source calls and 45 native velocity mixed4 identities.
- Both actual pressure mixed[P]/fields[P] replacements and their packet routing, retaining the common absolute-pressure normalization and y derivative convention.
- The actual generic mapper's ebracket expression and six remainder mixed2 source operators, together with the same functional K4 join.
- The actual physical q at power phase1 and exit phase0, and the identical source log factors.

These are function-level source joins. Numerical interval overlap is not used as a substitute. The resulting completed stress, divergence, diagonal and generally nonzero remainder join the accepted steep-exit physical companion.

## Reproduction and next dependency

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steeppowerphysical
```

This ordered stage runs the stress producer/checker followed by physical producer/checker. The actual physical cone is admitted separately by STEEP_POWER_CONE_2026_10_04.md and its current-source receipt. Next transport the same complete moments and physical/cone construction through original steep entry, then the angular repair/power/flatten regions.

Upstream finite-width feedback, independent global flat/volume/energy bounds, true n-dependent coefficient recursion, oscillatory correction and final corrected Cartesian NS validation remain unfinished.
