# Actual signed Rm defect and pressure inlet - 2026-10-15

Successor: [original leading implicit Rm patch reader](CURRENT_ORIGINAL_RM_PATCH_NUMERIC_2026_10_15.md)
now uses these signed unmodified defects to certify control Z2 and recover
leading patch fields, moments, pressure and radial coefficient. Finite-N
transport and fixed-time physical patch assembly are next; their incoming
histories remain separate and nonzero.


`OriginalRmDefectNumeric.evaluate(Z)` now exports the actual unpatched Rm inlet for arbitrary interior decimal Z. It supplies all five signed UNMODIFIED defect enclosures through Taylor Z2, the five finite-N own moment enclosures in generic S=Pstar units, canonical-frame and current-amplitude finite-N diagnostic defects, and the same-source analytic pressure datum through Taylor Z5. It preserves the original selected N, incoming correction histories and nonzero positive-order n1 primitive.

## Center the source before interval arithmetic

Use am=exp(-3/5)/(1+Z^2), epsilon=1/Lambda, and the actual reference/restoration histories. The signed background centers are

```
d0 = j*gamma1 + epsilon*Mpsi
d1 = (j*Kalpha_over_Phi + epsilon*Kpsi_over_Phi)/2
d2 = (H/Phi-5/4)/2
S  = j^2*gamma2 + 2*j*epsilon*Dpsi + epsilon^2*Apsi
d3 = S/(Pstar^2*am^2) - (B/Phi^2-5/6)/2
d4 = (C/Phi^2-5)/2
```

These are the source identities M-4Z, K/(2Phi)-4Z*H/(2Phi), H/(2Phi)-5/8, the centered energy defect and pressure defect. Computing raw order-one histories first and subtracting them would introduce rounding error much larger than the true defect.

The three centered angular memories are transported from the actual LONG TERMINAL inlet across L+1+1 with rates1.6,1.2,.2. Their exact finite memories are not reset to zero. Their numerical centers are zero proxies, with symmetric error enclosures proved from the whole incoming-coefficient logarithm minus its true finite decay distance. The numerical error scale exp(-64*logP) is used only after this inequality is checked. It is not a chosen field value and does not weaken an old accepted cap.

## Finite-N own moments and amplitude

The incoming C3 rows are unsigned ordinary-derivative log bounds. Retain their external1/N, physical chart measures and nonzero memories; convert coefficient k once by1/k!. They bound unknown SIGNED increments in

```
m = Mz/(R*Pstar)
h = Mtheta/(sqrt(2)*R^(3/2)*Pstar)
k = Mtheta_z/(sqrt(2)*R^(3/2)*Pstar^2)
e = Mztheta/(R*Pstar^2)
p = Mp/Pstar^2
```

The raw identities are m=M/Pstar, h=am*H/(2Phi), k=am*K/(2Phi*Pstar), e=A/Pstar^2-am^2*B/(2Phi^2), and p=am^2*C/(2Phi^2). In particular e is a combined energy row; do not assign its correction to A alone or use pressure p as both B and C.

`actual_finite_N_canonical_frame_defect_enclosures` keeps the original am normalization. `actual_finite_N_defect_enclosures` uses the actual modified am_N=am*exp(A/N). The latter retains exp(-a)-1 and exp(-2a)-1 separately, including the energy -5/12 and pressure5/2 cross terms. Small exponential increments are expanded before adding an order-one unit, so interval rounding cannot erase their scale. These are diagnostic enclosures of the actual finite-N inlet; the signed finite-N integral centers and phase are unresolved.

## One analytic pressure datum in the correct units

The core fixed P0 Taylor rows have core-ODE units epsilon*physical_P0. Multiply every row by exp(logLambda-2logP) to export physical_P0/Pstar^2. The same normalized datum is added once to generic p and excluded from all five defects. Its Z0..5 entries remain Taylor coefficients, not ordinary derivatives; the underlying analytic datum and every accepted source/error atom are retained.

## Focused receipt

```
python scripts/run_original_Rm_defect_numeric.py --write-report
python scripts/run_original_Rm_defect_numeric.py --inspect-report
```

The run exits0 in144.940s and saved-report inspection exits0. Four actual packets use Z=-.3,0,.3,.17. There are90 independent modified-radial-primitive/exact-five-identity checks,180 finite centered-moment coefficients,9 original-integer/incoming-history checks,17 pressure-datum/once-only checks and24 angular-error/control-unit comparisons. The independent fixture directly integrates the five raw physical density primitives under R=s^10; it does not integrate the implementation's centered recurrence. Source-current pressure is also compared against the independent scalar core export.

At all sampled points the d2/d4 errors lie below the original positive angular control unit10^10*t^2*Pstar^-2. This is a pointwise numerical prerequisite for the original implicit patch, not proof of whole-axis closure. Read-only Luna/max algebra/unit review found no material blocker.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_defect_numeric.json`,1032780 bytes, SHA256 `162d657ce701d95586fe0fa6daf18b1044dba66e366e45ce3e67eae51695698f`. It pins1758 source inputs and45 implementations. Accepted ancestor source/report bytes remain unchanged.

## Next source-preserving patch step

Evaluate/certify the ORIGINAL LEADING implicit map L*h+Q(h,invAm2)+d=0 using these actual signed unmodified defects and original bump weights. Recover the available Z2 controls using the same Jacobian; do not pad missing Z3..5 rows with zeros. Continue leading patch fields with exact x=R/Rm and dx/x. Separately transport the original finite-N five corrections from this nonzero inlet through the active patch toward Rh and the later reserved repair.

The source `current_original_true_Rm_patch_finite_N_C5_functions.py` explicitly keeps leading controls and finite-N corrections separate. Re-solving the leading map on the finite-N diagnostic defects would change that construction and is NOT the next acceptance condition. A separate finite-N perturbed-branch claim, if ever needed, would require its own uniform contraction certificate and amplitude variation; the existing leading-only certificate cannot establish it.

Current work does not export patch controls/fields, higher-Z defects, physical inverse/gradients, full Cartesian C2 or a global field. Global common N, pressure joins/energy/tails, exact heat, admissible stress/flat remainder, genuine n>=2 recursion/smooth sum, oscillatory cancellation and corrected NS remain open. Tasks14/15 stay IN PROGRESS; the full goal remains ACTIVE / INCOMPLETE.
