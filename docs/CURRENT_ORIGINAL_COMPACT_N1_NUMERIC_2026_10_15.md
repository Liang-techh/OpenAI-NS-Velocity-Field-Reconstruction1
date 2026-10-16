Superseding pressure milestone: [CURRENT_ORIGINAL_COMPACT_N1_PRESSURE_2026_10_16.md](CURRENT_ORIGINAL_COMPACT_N1_PRESSURE_2026_10_16.md). The velocity-only API below remains valid; the new wrapper adds the original complete compact P1 C0/C1 enclosure.

# Original reserved n1 velocity and own mass numerical enclosure - 2026-10-15

`OriginalCompactN1Velocity.evaluate(Z,x,log_tau="-2",azimuth=".25")` now adds the independently repaired temporal n=1 velocity to the original compact background on x=R/Rc in[1,2]. `OriginalReservedN1Numeric` supplies C0/C1 axial/swirl coefficient enclosures, supported primitives, cumulative physical mass M1 through Z1 and R*Q1 as a value. The original selected coefficients remain signed functions with directed uncertainty. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Full n1 pressure and genuine temporal n>=2 are still open.

## Same actual coefficient functions and disjoint supports

The reader consumes the checked original reserved C1 coefficient archive, exact linear operator and legacy C1 terminal-integral theorem. Their ownership chain is the same checked `full.coefficients`, reserved own fields and physical response, rather than unrelated numeric ranges. The original source family, fixed original N and positive mu remain.

Put s=log(R/Rc), L=log2, ell=L/400 and centers c_j=jL/40 for j=12,13,14,15,16. The exact new supports are [(j/40-1/400)L,(j/40+1/400)L], disjoint from the old background repair supports. With beta(t)=exp(-1/(1-t^2)) on |t|<1 and zero outside, let

```
Nhalf = integral_-1^1 beta(t)*exp(ell*t/2) dt
shape_j(s) = beta((s-c_j)/ell)/(ell*exp(c_j/2)*Nhalf)
b_j(s) = sqrt(2/Rc)*shape_j(s)
Uz1 = alpha1*b_12+alpha2*b_16
Utheta1 = beta1*b_13+beta2*b_14+beta3*b_15
F1 = Utheta1/sqrt(2R)
```

The physical azimuthal velocity uses Utheta1 directly. F1 is a different source quantity and is not substituted for Utheta1 in the vector.

## Radius cancellation and source-sector refinement

The archive uses alpha_j=alpha_hat_j/(mu*sqrt(2Rc)). Combining the bump prefactor before boxing gives

```
Uz1 = (1/Rc)*sum_j (alpha_hat_j/mu)*shape_j(s)
R*Uz1 = exp(s)*sum_j (alpha_hat_j/mu)*shape_j(s).
```

The reader reconstructs alpha_hat/mu from the actual signed original m11/m14 source sectors, with the Astar quotient derivative already present. It does not amplify the archive's relaxed exp(-1000) normalized cover by1/mu. Representable inner sectors are evaluated with directed exponentials; every unmaterializable tiny sector is separately proved below the positive exp(-64logPstar) cap, retaining its allowed sign. C1 source rows include the original axial quotient/product rules. Fixed coefficient units are Z-independent and are not differentiated again.

Angular coefficients retain their actual normalized C0/C1 enclosures and original physical unit. The combined azimuthal log unit is

```
logUnit_theta = logPressureUnit_m15-logPstar-logSigma-logRc
logSigma = -logPstar/2+1/2.
```

Its exact original pressure-unit recipe is decomposed into separate logRc, logU0, logPstar and combined inverse-mu atoms. This preserves the37/(2mu) term, the negative (4mu)^(-2/3) term, the original reference factor and finite offsets. No physical hyper-exponential radius, coefficient or amplitude is formed. An independent symbolic lift of the original `pressure_unit_log(...,power=2)` proves equality with the implemented atoms.

## Actual cumulative mass, partial bumps and terminal source identity

The incoming pre-Rc record is the normalized inner-cutoff primitive J_norm=rho*mean_u. It is converted exactly once to physical mass M_in=J_norm/Lambda^2. A reserved source already named physical M1/m11 would not receive this factor again.

The physical dR=R ds integral cancels Rc analytically:

```
M1(s) = M_in+sum_j (alpha_hat_j/mu)*I_j(s)
I_j(s) = exp(c_j/2)/Nhalf * integral_-1^clip((s-c_j)/ell) beta(t)*exp(ell*t) dt.
```

A degree3 symmetric Taylor integral with the fourth Taylor coefficient as full-cell remainder encloses each partial weighted bump. The full-cell coefficient encloses the actual fourth derivative divided by4! at every point, so Taylor-Lagrange gives the entire signed remainder as that interval times width^5/80; it is not a fourth-degree polynomial with an omitted fifth-order remainder. The complete normalizer/weights come from the checked original operator. Evenness of beta evaluates the remaining suffix using a reflected primitive with exp(-ell*t), avoiding subtraction of almost equal full and partial integrals. Source-local Taylor coefficients use their factorial once.

The same original C1 theorem reduces the actual full bump integrals to the selected repair response and proves terminal physical m11=0. Therefore both forward M_in+partial and backward minus-the-suffix are valid C0/C1 enclosures of the same cumulative mass; their coefficientwise intersection is exported. Incoming history is not reset. After the final axial bump has ended, the suffix is exactly empty, so M1 and M1_Z are zero by the checked source identity. This is not a new independent numerical residual of resolved repair coefficients. Full leading background histories remain nonzero.

With the true mass derivative and local axial source,

```
R*Q1 = [2Z*(R*Uz1)-(1+delta)*Z*M1-(1-Z^2)*M1_Z]/(1-delta*Z^2).
```

The positive-order factor is1+delta. R*Q1 is a value; Q1_Z and physical gradients are not exported from C1 mass data.

## Physical temporal orders0 and1

The wrapper appends

```
psi1 = lambda^(1+delta)*M1
uz1 = lambda^(-1+delta)*Uz1
utheta1 = lambda^(-1+delta)*Utheta1
ur1 = lambda^(-1+2delta)*(R*Q1)/sqrt(2R).
```

These reserved coefficients are already in their physical source units. No additional Lambda^-2 is appended. The normalized cutoff mass conversion above is the one required conversion. Atomic radius/time factors, the factor110 and Cartesian rotation remain consistent with the background. Ordinary and scaled-time outputs retain the original delta cancellation before time comparison.

The bounded exact source has one shared meridional streamfunction. Interval centers are numerical proxies, not definitions of exact coefficient functions or a differentiated global point oracle. Finite differences of these independently summarized centers have not been certified divergence-free. The actual-source structural recovery, full error intervals and missing physical-gradient/global-join gates remain explicit. The current physical pressure is temporal-order-zero only; P1 is not silently replaced by a computed zero.

## Focused evidence

The focused run exits0 in86.979s; saved-report inspection exits0. It exports33 physical leading+n1 packets across3 axial positions, all five new bump centers, two partial axial bumps, Rc and the zero-velocity suffix, including ordinary/scaled time. Checks include10 independent defining integrals of weighted beta, an independent original physical beta-unit identity,90 source/N/scope conditions,90 mass-intersection/recovery conditions,42 actual support conditions,120 physical-unit/order/pressure-scope conditions,9 Rc normalized-to-physical mass seam conditions,54 source-zero-suffix/nonzero-background conditions,21 scaled-time conditions and4 rejected requests. Maximum tested raw weighted-beta integral width is2.811908355347347e-8.

The successful actual-source smoke preceded this run. The first driver attempt found a Python import-name collision between the inspection module and the saved-report inspector; that driver-only collision was fixed before the accepted run. Accepted ancestor source/report bytes are unchanged.

The existing read-only worker `/root/n1_atlas_map` reports GPT-5.6 Luna / max. It confirmed the exact legacy C1 coefficient/terminal ownership chain, the normalized inner mass conversion, beta physical unit, supported primitives and physical temporal weights. The full-cell C4 Taylor-Lagrange contract was checked explicitly; no remaining material enclosure or unit issue was identified. No new agent was spawned.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_n1_numeric.json`;24,123,543 bytes; SHA256 `182801c4ca21480bcc5d182e719d9fc35ec2ef6b763f0839afe2af977f8198d6`;1825 scientific input pins and72 implementation pins.

## Next original construction work and limits

1. Numerically assemble P1 from its incoming core/lower/active pressure, actual compact fixed leading N1p and local swirl-pressure integral. Preserve the same pressure source and physical dR measure. The pressure terminal radius is Rv; velocity/mass terminal closure at2Rc does not zero interior P1.
2. Continue the original nonzero-mu quiet background from2Rc toRp, with full leading histories and the source-admitted n1 zero-velocity tail. Use stable analytic exprel/own-rate kernels for the long quiet length Tw-2-log2; the short compact polynomial moment implementation is not a long-range point oracle.
3. Resolve signed coefficient functions/finite-N phases where viable, then establish higher mixed rows and global n1 joins. C0/C1 enclosures are not numerical C3 bounds, and current Q1 values do not supply Q1_Z.
4. Close signed S0/S1/total Q_y pressure, common global N/physical C2, compatible analytic preheat/exact heat exterior, physical energy/radial tail, admissible stress and flat remainder.
5. Construct genuine independent temporal n>=2 recursion, smooth sum and two-family oscillatory stress correction, then independent corrected Cartesian NS validation. One repaired n=1 velocity layer is not temporal scale recursion or a full corrected solution.
