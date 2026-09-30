# Continuous incoming axial field installation

`ContinuousIncomingOuterField` extends the continuous seeded exterior jet
adapter. Before Rp, axial point values and cumulative means come from the
same continuous incoming cutoff. The actual corrected inner mass and its Z
derivative are transported once; no finite terminal mass is replaced by zero.

The installer regenerates the reference `I_z`, `I_theta_z`, and `I_uz2`, then
reapplies the raw measured inner moment offsets to the normalized inputs.
The energy target changes by old prior minus new prior, preserving the same
future angular target. The matrix, pulse and energy atoms are reused to solve
for new coefficients, and those coefficients feed both values and means.
The mixed integral retains the original angular amplitude primitive. Its
float-backed provenance, and the inherited swirl energy, remain explicit.

The resulting reference factors are carried into the zero-Z input tangent
path rather than reconstructing them from a different legacy integral.
Incoming mean Z derivatives use the exact axial factorization plus actual
inner moment derivatives. Future angular target derivatives still use the
existing declared finite difference.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_incoming_outer.py`
to regenerate the installed shared-candidate receipt. It compares cumulative
mass and its Z derivative on both sides of Rp and samples incoming velocity
and mean at four locations. This is numerical matching evidence, not an
integral error enclosure or a proof of global finite energy, full NS stress
closure, or scale recursion. Angular jets, all exterior moments, and exact
terminal mean closure remain outstanding.
