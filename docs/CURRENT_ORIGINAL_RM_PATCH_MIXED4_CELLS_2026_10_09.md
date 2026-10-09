# Successor: actual Rm full signed generic inputs

Checked source [a38f648b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a38f648b11045271de0427c8086c29d6ee2cf101); evidence/tasks [CURRENT_ORIGINAL_RM_GENERIC_INPUTS_2026_10_09.md](CURRENT_ORIGINAL_RM_GENERIC_INPUTS_2026_10_09.md). The new live leading-patch owner now supplies full signed inertia/shear and relaxed numerators in its five canonical bases. Actual quotient/root/phase/density/integration, whole-Z, finite-N Rc/all24/global N and real n-recursion remain open. Historical text below is unchanged.

---

# Actual leading Rm patch mixed4 and closed radial source cells

Checked source [e50bec23](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e50bec23330e72b21af3c2cfa6bfe1db87379de1). Full reconstruction **ACTIVE / INCOMPLETE**. The actual leading Rm patch now supplies physical velocity, pressure and all five primitives with mixed derivatives of total order<=4 in x=R/Rm, y=logR and R coordinates. Current-radius raw rows and closed radial cells include all three flat bump supports and the exact Rh endpoint. Native Z=0,.5 remain conditional source frames. Whole-axis, actual finite-N density integration/Rc repair/all24/global N and real n-recursion are still open.

Predecessor: [actual Rm defect/inverse/partial patch](CURRENT_ORIGINAL_RM_DEFECT_PATCH_INVERSE_2026_10_09.md). The same accepted coefficient functions, centered incoming Rm data, original analytic P0, canonical Rm/Pstar units and beta normalization are retained. Only accepted source callbacks run; ancestor producers/full checkers and old pressure/repair constructors are not rerun. Root implemented and checked the work. The existing **GPT-5.6 Luna / max** worker reviewed source indexing, prefactors, raw units, Q scope and provenance read-only.

## Original radial source equations

H=x^.1+f and V=4Z+g have ordinary x derivative rows0..4 using the original normalized compact beta. Zeroth partial primitives come from the actual CURRENT x value of the accepted patch, rather than reinitializing at Rm. Let G=x*m. The source equations are

    G_x=V; theta_x=sqrt(x)*H; mixed_x=sqrt(x)*H*V
    energy_x=g^2*invAm2-H^2/2
    pressure_primitive_x=H^2/(2x).

Primitive derivative k+1 uses source derivative k. This indexing is AST-bound to the original append statements. Ordinary derivatives of G/x recover m through4, retaining the full nonzero incoming partial history. Q is recovered by the original denominator formula for every radial row. Since Q differentiates m in Z, Q has ordinary Z Taylor coefficients0..4 only. A temporary internal zero padding enables arithmetic but is never exported as a Z5 certificate.

The operator uses directed beta ordinary derivatives through4 with radius1/40, original source normalization and centers5/4,3/2,7/4. Exact rational support edges49/40,51/40,59/40,61/40,69/40,71/40 use the exact flat endpoint identity. Closed cells crossing an edge use the accepted vanishing flat-tail envelopes; an interval crossing support never takes a reciprocal of a zero-containing denominator.

## Physical mixed derivatives and pressure memory

Axial amplitude products are differentiated before the final grid; shared Z-independent radius units are attached once after coordinate conversion:

    Utheta=Pstar*am*H; Uz=V; Ur=sqrt(Rm/2)*sqrt(x)*Q
    Mz=Rm*G
    Mtheta=sqrt2*Rm^1.5*Pstar*am*theta
    Mtheta_z=sqrt2*Rm^1.5*Pstar*am*mixed
    Mztheta=Rm*Pstar^2*am^2*energy+8Z*Rm*G-16Z^2*Rm*x
    Mp=Pstar^2*am^2*pressure_primitive
    P=Pstar^2*P0+Mp; am=exp(-.6)/(1+Z^2).

The centered energy is uncentered before reporting physical Mztheta. Its baseline term has its own unit, while the centered term carries Pstar^2. The am/am^2 Z products and sqrt(x) radial derivatives are retained. Rm powers are combined as Rm^(a-k) before attaching a physical unit, so a gigantic directed logRm offset is never independently subtracted from itself. Normalized x/y grids are exported alongside the physical grids. Original P0 appears only at radial order0, including its ordinary Z derivatives there. Positive radial derivatives contain only the pressure increment. Every point, cell and Rh path exposes separate P0, pressure-axis and increment rows; mismatched P0 objects reject.

    Dy^k=sum_j Stirling2(k,j)*x^j*Dx^j
    DR^k=Rm^-k*Dx^k.

Ordinary Z Taylor coefficients become derivatives via n! once, at final grids. Each x/y/R mixed4 grid has15 entries with k+n<=4. These grids are not a claim of Z0..5 for every radial order. Retained underlying primitive/velocity source rows are ordinary x derivatives with axial Taylor jets through5 where available; radial Q/Ur exports stop at Z4.

## Raw current-radius rows for later density compilation

The adapter exports the original current-radius source conventions:

    m=Mz/R
    h,k=am*x^-1.5*theta,mixed
    e=Mztheta/(R*Pstar^2)
    p=Mp/Pstar^2
    absolute_pressure=P0+p
    theta_velocity=Utheta/Pstar; axial_velocity=V
    radial_velocity_rows=(Dy+.5)^k Q.

All five histories are converted from x to y derivatives exactly once. The radial half shift includes the differentiated sqrt(R) prefactor. Absolute pressure is assembled from ALREADY y-converted p rows, with P0 at row0, and receives no second coordinate conversion. Independent closed-integral checks caught and resolved this duplicate pressure-conversion error before acceptance.

These are genuine conditional leading-background raw rows. They are not the missing finite-N Rc solution or complete five-density oracle. Exact phase, shear, both N-dependent levels, actual stress/raw/inertial density assembly and source integrals still have to be connected.

## Closed source cells and exact endpoint

`OriginalRmPatchMixed4Cells()` requires the accepted new receipt. `evaluate('0',(5,4))` returns a point; `cell('.5',(17,10),(9,5))` returns a directed enclosure on the whole closed cell; `evaluate('0','Rh')` uses exact x=exp(1); `cell('0',(71,40),'Rh')` encloses the complete terminal-support window. Rational cells must be strictly ordered and within[1,e], while Rh is accepted only as an exact right endpoint. Whole-Z requests reject.

Partial beta integrals on a cell retain monotone endpoint enclosures, actual signed control coefficients and original incoming defects. A cell spanning the active support does not reset terminal moments. Only a cell entirely beyond71/40 uses the full-weight identity of the SAME unique leading implicit solution. Its original defects and unrefined residual diagnostic remain recorded. The cell is not an interpolation of samples.

On the exact complete-support window, the same unique leading-map identities give G=4Z*x, m=4Z and V=4Z. The radial formula is refined algebraically to Q=4*((2+delta)*Z^2-1)/(1-delta*Z^2), with positive x derivatives of m/Q exactly zero. This avoids spurious cancellation of leading terms in tiny high-Z coefficients. Active and crossing-support cells retain their genuine partial primitives. Rm has an open zero-correction neighbourhood. At Rh every bump and its derivatives vanish and every partial weight is complete, so the same implicit identities join to the exact reference branch through mixed4 at each admitted source frame. This does not establish all inner/outer interfaces or whole-axis functional closure.

## Evidence and limits

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells.py, .json.gz, _check.py, _check.json. Producer 52.782s; checker 59.579s; terminal exit0. **361 staged exact dependency hashes PASS**. The gzip report preserves exact source records and is hashed as bytes.

15 exact closed primitive/Stirling/radial-shift identities pass. Two independent finite polynomial-basis fixtures derive nonzero current primitives by CLOSED integration, then symbolically differentiate the physical expressions and raw source units. 810 physical x/y/R mixed4 and 270 raw current-unit comparisons pass. These are diagnostic finite fixtures, not native beta parameters or source-value selections. Numerical reference tolerance remains separate from production directed interval certification.

Native evidence preserves 9720 factored physical mixed4 records, 10 current primitive object joins, 36 separate P0 checks, 180 exact edge-derivative zero checks and 4 endpoint source joins. 3240 cell-versus-native-point enclosure comparisons pass in shared physical units. Algebraic prefactor checks bind normalized and physical rows without expanding a giant exponential; pressure inclusion compares the increment and then translates by the exact same P0. 15 invalid requests reject. Raw radial rows never promote a padded Z5 coefficient.

The new gate closes full leading-patch radial mixed4/raw-unit/cell availability at two native frames. Whole-Z analytic source providers, finite-N density compilation/Rc coefficients/all24/global N, completed heat/stress/flat layers and real n-recursion remain false. This is progress in the matched-background dependency chain, not completion of recursive scaling or corrected NS.

## Detailed next production tasks

- [x] **FULL LEADING PATCH x/y/R MIXED4.** Original beta derivatives, inherited current primitives, k+1 RHS indexing, Q scope, physical prefactors and pressure memory. No finite differences of enclosure caps.
- [x] **RAW CURRENT UNITS / CLOSED PATCH CELLS / Rh.** Original history normalization, radial half shift and single pressure conversion; directed whole-cell partial integrals, support crossings, exact edges and full terminal window.
- [ ] **NEXT CONNECT PATCH DENSITIES.** Feed new raw rows and canonical units into original generic-shear/fixed-phase Z density calculus. Preserve exact all-N phase, current actual shear/denominators, both N-dependent levels and full inertial terms. Integrate active[1,71/40] and terminal[71/40,e] with the same solved leading coefficients; leading cells alone are not finite-N density integrals.
- [ ] **NEXT REFERENCE/RESTORATION/POSTRESTORE SOURCE CELLS.** Extend actual Rsh->Rm functions to closed radial cells with true sigma prefix changes, exact amplitude/radius factors and signed incoming tails. Keep pressure P0 and high derivatives separate; compile the actual five-density windows. Accepted Rh_reference is a different window.
- [ ] **WHOLE-Z SOURCE-CORRELATED PROVIDERS.** Replace cached frames0,.5 by axial cells from analytic core/fixed-point, anchored pole primitive and fourteen pressure atoms, with derivative/error, midplane/pole and uniform operator controls. Two-point leading jets do not close functional terminal identities on the axis.
- [ ] **CONTINUOUS CORE/MICRO/SWITCH/LONG RESHAPE.** Complete both original microscopic functions/cells with dy=hb*ds and true current modes/inlets; full long-reshape kernels for0<=y<=T retain actual anchored B, incoming/body/tail source and P0.
- [ ] **ALL24 TRUE COMPLETE TARGET / FINITE-N Rc REPAIR.** Assemble all original source windows and Rc_E/Rc_E_Z, then solve B*h+N*r+Q/N=0 with exact same owner/N, unique Jacobian/derivative bounds and outer feedback. Do not substitute the leading map controls for this repair.
- [ ] **GLOBAL N / FIVE AXIAL TERMINAL FUNCTIONS.** Reconcile repair, cone, interface, centered-target N^-2/Picard/limit/tail inequalities with one N; close functions of Z rather than sampled leading values.
- [ ] **MATCHED HEAT / FINITE ENERGY / STRESS / FLAT.** Complete annular/flatten/pressure interfaces, exact heat exterior, axis regularity and energy tail, regional admissible stress cone and flat remainder. No claim follows solely from these patch derivatives.
- [ ] **REAL n-RECURSION / MEAN CORRECTIONS / TWO PULSES / CORRECTED UVW.** Actual n-dependent recovery/repairs/smooth summation with divergence-preserving cutoffs, averaged stress cancellation, independent corrected Cartesian residual and separate contraction/elongation/recursive-scale/material-winding diagnostics.

Mark DONE with defining code, scoped report/receipt and a commit; continue production after necessary changed-scope checks. Full objective remains active.
