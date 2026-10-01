# Candidate-specific centered five-row defects

The actual refined R110 candidate endpoint is now an input to `CenteredComponentDefects`, replacing the old stored Lambda1e36 reference source for this computation. Run `experiments/root_st073/lei_ren_part1_paper_candidate_component_defects.py` after generating `candidate_transition_analytic.json`.

The endpoint reader validates transition dependencies, restores exact MP tuples, exposes the original P0/P0_Z, and accepts only R = 110 at the actual center Z = 0.3. It never reconstructs tiny angular amplitude from a shortened logarithm. The five rows and first axial derivatives use the existing component/dual algebra and are stored with exact MP atom tuples and separately labelled contributions.

Nominal values at this center:

| Paper row | Nominal defect |
| --- | --- |
| d1 | 2.2343710096e-15 |
| d2 | 5.6898620665e-16 |
| d3 | negative, nonzero, extraordinarily small; see exact atom tuple |
| d4 | 5.9150131214e-41 |
| d5 | negative, nonzero, extraordinarily small; see exact atom tuple |

These are normalized matching defects, not momentum residuals or proof that terminal moments vanish. The source is a fixed-parameter MP midpoint projection. Atom storage does not establish a source pressure/exit-width parameter family. Quadrature and source parameter remainders, the infinite radial tail, and the axial domain remain unenclosed. Formal assembly is not a solved repair.

Next concrete tasks:

1. Feed these new d rows to the actual `FiveBumpMomentMap`, using correction target `-d`, canonical Am and the same Z. Do not load the old `reference_defect_Z03.pkl` or old inverse receipts.
2. Retain nonzero d3/d5 and first axial tangents when solving controls. Report the five changed moment rows independently of the Newton stopping residual.
3. Compute actual bump product integrals and their error bounds before claiming a directed inverse certificate. Check the repaired field and stress throughout bump supports.
4. Replace point inputs by axial functional bounds and solve/verify repair functions on the complete axial interval. Enclose source and quadrature errors in that same problem.
5. Complete outer/preheat/heat compatibility and verify the five terminal identities before claiming a matched leading background.
