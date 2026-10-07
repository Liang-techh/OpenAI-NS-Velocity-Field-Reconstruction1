# Uniform finite-N repair controls, cumulative defects and quiet primitive shear

Latest successor: [CURRENT_RADIAL_PRESSURE_ERROR_BOUNDS_2026_10_06.md](CURRENT_RADIAL_PRESSURE_ERROR_BOUNDS_2026_10_06.md), implementation [d446dc0e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d446dc0ecc8756be088225c162c18bd9b9dbc431). BOUND1b1c and BOUND1b2c native error bounds are implemented. Completed signed tensor and global physical/recursion gates remain open.

Implementation and focused receipt: [23dffae0](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/23dffae0e1efd146de3154904d1eb7677063f44a). This extends the accepted independent five-bump map to queries in N. The old signed defects and control boxes remain data only for their original N=10^12.

## Result and API

CurrentUniformRepairMajorants in experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_uniform_repair_majorants.py supplies query(N,logx=(0,1)). N must be a positive Python integer at least27,303,666; logx is a finite subset of[0,1]. No ancestor graph is rebuilt.

The same current feedback family, implicit source and analytic axis pressure datum are checked against the accepted independent-repair and repaired-history receipts. Fixed bump weights and the divided linear inverse are independent of N. The exact source equation is B*h_N+Q_N(h_N,h_N)=-d_N, where Q_N has its explicit1/N and the source-correlated D=(J-M)/mu is formed before enclosure. The saved N=10^12 signed defects and coefficients are never reused at other N.

Let CA=||B^-1||inf, CQ bound N*||Q_N||, D bound the actual normalized source defects and R=2*CA*D. The public query returns the uniform implicit h-ball, quadratic Jacobian cap2*CQ*R/N and inverse nonlinear Jacobian capCA/(1-kappa), kappa=2*CA*CQ*R/N. The accepted repair-only threshold gives strict contraction and inclusion. This is an existence/bound interface; it does not resolve arbitrary-N signed point coefficients.

Three disjoint normalized log bumps are fixed at1/5,1/2,4/5 with radius1/40. Raw beta derivatives are bounded by their original polynomial formulas and exp(-1/w)*w^(-2k) maxima. Product-rule log jets through4 and C(Z)=1/(1+Z^2) axial derivatives through5 retain common Pstar*Ad*exp(-1-3mu/2) outside. At most one bump is active, so the bound uses R*G_j, without multiplying by the number of bumps.

## Cumulative defect and support strips

The API returns prefix and remaining repair primitive bounds, followed by bounds for all five source cumulative defects in the original complete physical units. Axial/swirl cross terms and both square terms remain in the map. The cumulative factors are:

| Defect | Physical factor outside the normalized bound |
| --- | --- |
| M | A*R0 |
| J | sqrt(2)*A^2*R0^(3/2) |
| I | sqrt(2)*A*R0^(3/2) |
| S | A^2*R0 |
| Cp | A^2 |

A=Pstar*Ad*C(Z)*exp(-1-3mu/2), R0=Rw*exp(1). These factors are returned as logarithms and C powers; giant amplitudes or radii are not exponentiated.

Actual quiet q=1+y gives the original source offset t=1+q=2+y. This lies after the cutoff support ending at1/2. The scalar source integral has already stopped, and only the original homogeneous transport remains. Its exponent cancels the physical radial power when expressed in fixed A/R0 units. With the same implicit terminal equation, total defect equals negative remaining repair primitive. Thus the cap is min(incoming+prefix,remaining), retains shrinking remaining strips and is exactly zero beyond every bump. It does not zero the incoming moments or fit terminal samples.

## Repaired quiet primitive margin

The repair is supported in this quiet interval, disjoint from modulation. With alpha=1/2+mu and D_N=exp(-alpha)-mu*R*G0/N>0:

`|a-a0| <= 2*mu*R*(G1+alpha*G0)/(N*D_N)`

`|b| <= 2*sqrt(mu)*R*G1/(N*D_N)`.

The conservative sufficient threshold **N>=68,533,403** gives a-2>=3mu/2 and consequently a-2+b^2/(2+3mu)>=3mu/2 throughout the quiet repair band. The existing N=10^12 passes this estimate. This threshold certifies only repair and this primitive inequality. It does not select a common N for completed tensor cones.

Upstream flat tapers still use the separately accepted correlated primitive theorem M_N>=2mu*sigma+mu*chi^2/4 for N>=1. The new quiet bound does not replace that loop theorem. Full tensor admissibility requires diagonal, pressure, directional and alignment bounds as well.

Evidence:74 exact source/Jacobian/raw-beta/log-jet/axial/transport/terminal identities;12 independent Jacobian/inverse ball fixtures;45 raw derivative samples;30 independent bump/quiet quotient fixtures;184 continuous prefix/complement integral comparisons;7 current-source queries and6 invalid-query guards;915 working/index dependency hashes. Worker source/math review: **GPT-5.6 Luna / max**, read-only, no tests/constructors/edits/descendants. Run producer and checker directly using the two .py files; their .json receipts are committed.

## Tasks for subsequent agents

- [x] **BOUND1b1a:** actual all-N normalized terminal source caps, including correlated D before enclosure, exposed without fixed-N control substitutions.
- [x] **BOUND1b1b:** uniform implicit control ball and nonlinear Jacobian inverse for every N above the repair-only threshold.
- [x] **BOUND1b2a:** original raw bump ordinary4/axial5 profile caps, prefix and shrinking complement primitives, fixed physical units and exact terminal closure.
- [x] **BOUND1b2b (quiet primitive):** quotient perturbation estimates, positive swirl denominator and sufficient quiet primitive reserve.
- [x] **BOUND1b1c (native bounds):** expose upstream modulation signed partial histories and their mixed logR/Z jets for variable N. Keep common source correlation and separate positive kinetic terms; do not subtract two independent J/M boxes to recover D.
- [x] **BOUND1b2c (native bounds):** convert own M/Cp bounds to radial velocity and absolute pressure error jets. Bind increment_rows and its differentiation recurrences. Preserve original axis datum, radial powers, C versus C^2 derivative caps and Pstar sectors. Cover all partial supports and mixed logR4/Z5 rows; terminal function-level zero is insufficient by itself for these derivative bounds.
- [ ] **BOUND1a5c / BOUND2:** insert these bounds and the correlated upstream taper inequality into the actual completed signed tensor. Establish diagonal/residual/pressure, alignment and both directional vector margins throughout O2/O3. Preserve pointwise phase/cutoff correlations where a source reserve vanishes.
- [ ] **BOUND2/3 / COMMONN:** combine repair, radial, pressure, alignment and vector estimates to choose one common finite N. Publish every sufficient inequality. Do not promote27,303,666 or68,533,403 to a whole-cone threshold.
- [ ] **CONT4f2b:** add reference_restore/restore_buffer physical providers from the accepted endpoint germs and original radii. Then compose every original/current-modified adjacent and internal interface without relabelling the obsolete full33 packet.
- [ ] **ENERGY / REC:** bound physical kinetic cross terms and radial/time tails; implement actual n=1 and n>=2 recovery equations, a common core domain, independent per-order repairs, divergence-preserving truncation and smooth summation. Local inverse kernels or frequency N changes are not coefficient recursion in n.
- [ ] **WAVE / PHYS / DYNAMICS:** actual oscillatory/mean correction and averaged stress cancellation, resolved physical u/v/w/p, corrected Cartesian NS residual, measured contraction/relative axial elongation/material winding.

The successor supplies the remaining native upstream mixed-history and radial/pressure error bounds. Completed tensor/directional analysis remains open and is not implied by these source error bounds. Global interfaces, full modified admissible stress, common N, finite energy, actual coefficient recursion and corrected NS remain unadmitted. Preserve the persistent full goal and unrelated dirty files. The controller's last full runtime stage remains currentmodifiedphysicalvelocity. Previous report: [CURRENT_RSH_INTERFACE_AND_CORRELATED_SHEAR_2026_10_06.md](CURRENT_RSH_INTERFACE_AND_CORRELATED_SHEAR_2026_10_06.md).
