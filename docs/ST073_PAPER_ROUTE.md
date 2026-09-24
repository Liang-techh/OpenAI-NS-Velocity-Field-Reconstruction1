# ST073 continuation guided by the OpenAI paper

Primary source: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
Read against the live official PDF on2026-09-24. Statements below describe the
paper construction; they do not certify our implementation or its theorem.

Source map: Section2 describes anisotropic contraction and annular residual
cancellation by oscillatory momentum flux. Section4/Proposition4.2 identifies
leading tangential stress; its operator removes axial viscosity, so it must
not be equated with our full corrected ST073 residual. AppendixA constructs
radial moment repair and a heat exterior. Sections7-9 address oscillatory stress
realization, mean corrections and residual improvement; Section10 handles
compact forcing and whole-space construction.

Implementation decisions: reuse the bundled heat evaluator, but recompute
matching data from corrected ST073 jets. Do not reuse old leading moment/cone
values without a derivation. Prioritize actual radial moment/stress compatibility
before treating material winding as evidence of the whole construction.

New direct-join experiment: heat_join_screen.py freezes one heat amplitude
matched to inner midplane swirl at k0. At k6 midplane, swirl value jump is only
.0013394 but rtheta stress jump is-10.19394; off-midplane axial velocity jumps
reach2.10343. Thus value matching alone cannot attach this exterior. No fitted
force was introduced. The pure heat-swirl exterior is not globally finite-energy
in z and does not carry the inner axial/radial flux.

Next solve a transition carrying these actual stress and transport mismatches,
then derive moment constraints and test a realizable correction stress. Preserve
the local recurrence as a boundary-data provider, not a replacement for the
paper's outer correction mechanism. Source-inspired modifications and autonomous
choices must be labeled separately. Full momentum max/L2 gates stay1e-3.
