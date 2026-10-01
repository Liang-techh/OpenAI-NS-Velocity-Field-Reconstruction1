# Pressure-integral error propagated to the collar endpoint — 2026-09-30

The accepted coherent pressure-integral error envelope now enters the same finite radial recurrence and reaches the analytic collar endpoint. It remains a conditional pressure-error contribution; original axis/source parameters, adapter-generation roundoff, infinite radial series, omitted pressure powers, and omitted width orders are not certified.

The new calculation uses the committed exact axis input at Z=.3, radial degree18, pressure order9, width order2, and 473-digit directed arithmetic. The 22 initial axis Taylor terms require pressure bounds through derivative order21, within the existing accepted coherent order24 receipt. The axis data are shared with exactly zero perturbation.

## Pressure grouping and propagation

The stored P0 is a polynomial whose pressure power zero contains the pre-Rv prefix and power one contains the post-Rv tail. Both groups contain mathematical pressure integrals; the prefix is not simply the reference constant. For every Taylor term k, the combined normalized bound Bk is converted to exp(28)Bk and conservatively applied to both groups. Each subset's error is bounded by the positive combined triangle bound, so this is valid but may double the seed bound at pressure=1. Nonzero error atoms are retained even when the corresponding nominal coefficient is zero.

The finite radial recurrence keeps nominal and perturbation coefficients separate. Regenerated rows enter the same interval inlet and first/second-width formulas used by the previous endpoint calculation. P0 remains the shared nominal pressure datum with its propagated integral-error enclosure; no pressure tail is appended afterward.

The output reports each width coefficient separately. Pressure powers zero through nine are summed only to bound the retained pressure=1 perturbation; the tiny width h_b is not materialized. Thus the width-w output is a bound on that coefficient divided by h_b^w, and the physical contribution requires multiplication by h_b^w. These are not full-field error or residual bounds.

## Representative conditional error bounds

| Quantity | Width zero | Width one coefficient / h_b | Width two coefficient / h_b^2 |
| --- | ---: | ---: | ---: |
| Uz value | 4.507310936091168e-36 | 2.253655468045584e-36 | 7.910684621753725e-36 |
| Ur value | 2.415265927005230e-53 | 1.150264261201741e-52 | 2.379668446876275e-52 |
| Ur_Z | 2.096641867791618e-52 | 9.908799217651555e-52 | 1.914554982256747e-40 |

The larger Ur_Z second-width bound is retained as computed; it is not replaced with an expected smaller number. The pressure-value width-zero bound is0.6417014430653714 in the unnormalized profile pressure units. This is an absolute pressure-integral error contribution, not a momentum residual. Exact intervals and every retained atom are in the receipt, including all five endpoint moment contributions.

## Nominal regeneration and scope

The regenerated nominal finite recurrence is also compared with the older rounded coefficients:760 comparisons per field, including the first four axial Taylor entries and ten pressure powers at each radial degree. Outside counts are396 for F,432 for Uz, and360 for P. These differences are preserved and do not get folded into the pressure-integral error. They require a separate generation-roundoff comparison to make a claim about the earlier stored field.

Run `python experiments/root_st073/lei_ren_part1_paper_interval_pressure_perturbed_core_check.py`. It uses only committed exact input and existing accepted error receipts, regenerates the finite coefficient recurrence, and writes an exact endpoint receipt. Runtime is about100 seconds on the current host; the expensive spatial collar/source chain is not rebuilt.

Next: account for original axis/adapter-generation errors, bound retained-versus-omitted pressure and width orders, and propagate the endpoint bounds through the whole collar/annulus and axial domain. Global functional moment closure, admissible stress, relative flatness, and true temporal n-dependent recursion remain open.
