# Continuous incoming axial field installation

`ContinuousIncomingOuterField` extends the continuous seeded exterior jet
adapter. Before Rp, axial point values and cumulative means come from the
same continuous incoming cutoff. The actual corrected inner mass and its Z
derivative are transported once; no finite terminal mass is replaced by zero.

The installer regenerates the reference `I_z`, `I_theta_z`, and `I_uz2`, then
reapplies the raw measured inner moment offsets to the normalized inputs.
The energy target is regenerated from the new prior and the live angular
tail nominal contribution, including the shared bump atoms. The matrix,
pulse and energy atoms are reused to solve
for new coefficients, and those coefficients feed both values and means.
The mixed integral uses the continuous MP angular amplitude primitive.
Swirl energy, measured inner offsets and the live tail model
still retain their inherited uncertainties.

The shared angular correction now owns continuous bump atoms for point
values, coefficient equations, pressure and energy integrals. Its coefficient
Z derivatives come from the implicit two-equation solve. The actual swirl
jet includes `amplitude * h_Z`, which the previous adapter omitted after
flattening. `h`, `h_Z` and `h_t` remain separate MP values: forming `1+h`
alone can round a tiny nonzero correction to one. Pressure exposes the
backward bump contribution and its analytic Z jet separately from baseline
pressure; no full pressure derivative is claimed.

The resulting reference factors are carried into the zero-Z input tangent
path rather than reconstructing them from a different legacy integral.
Incoming mean Z derivatives use the exact axial factorization plus actual
inner moment derivatives. Future angular target derivatives still use the
existing declared finite difference.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_incoming_outer.py`
to regenerate the installed shared-candidate receipt. It compares cumulative
mass and its Z derivative on both sides of Rp, samples incoming velocity
and mean, checks installed angular jets on both bump supports, and compares
the bump pressure Z jet with a fourth-order independent difference. This is
numerical matching evidence, not an
integral error enclosure or a proof of global finite energy, full NS stress
closure, or scale recursion. Heat/baseline input bounds, all exterior moments,
complete pressure jets and exact terminal mean closure remain outstanding.
