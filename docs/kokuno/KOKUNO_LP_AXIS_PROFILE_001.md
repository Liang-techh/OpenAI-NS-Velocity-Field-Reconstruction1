# KOKUNO-LP-AXIS-PROFILE-001

## Scope

This increment turns the explicit leading-axis data in KokunoYumeto's public Navier--Stokes reconstruction workbench into executable, vectorized profile code. It is intentionally smaller than a full leading-profile or 3-D velocity reconstruction.

## Public source and release provenance

- corrected public release: Zenodo record `22678406`, corrected 2026-09-09;
- corrected PDF name: `released_ns_reader_corrected_20260909.pdf`;
- executable formula evidence used here: `KokunoYumeto/yang-mills-interacting-workbench` commit `143f6773feb424ad9ed3a8d116653200f20346b7`, file `navier-stokes/navier_stokes_workbench.tex`.

The corrected PDF is release provenance for this increment, but it has **not** been independently parsed in this implementation pass. Formula claims below are pinned to the public workbench commit rather than being presented as a fresh PDF transcription.

## Source formulas implemented

The source sets

- `A = 1/2 + h`, `D = 1/2 - h`, with small positive `h`;
- `d(eta) = 1 - eta^2`;
- `0 < j0 <= 0.05` and `U_*(eta) = 4 eta + j0`;
- `H_*(eta) = D eta + d(eta) U_*(eta)`;
- `W_*(eta) = 1 - d(eta) U_*'(eta) - 2 D eta U_*(eta)`.

The workbench also records the simplified identity

`-W_*(eta) = 3 - 8 h eta^2 + (1 - 2 h) j0 eta`.

The executable API additionally exposes the analytic eta derivatives

- `U_*' = 4`;
- `H_*' = D - 2 eta U_* + 4 d`;
- `W_*' = 16 h eta - (1 - 2 h) j0`.

`KokunoLeadingAxisProfile(h=0.005, j0=0.025)` uses autonomous demo values inside the public source bounds. They are not recovered OpenAI parameters. No amplitude parameter is exposed at this layer, so there is no `u -> 0` scaling knob.

## Coordinate mapping status

The public workbench uses the native similarity coordinate

`eta = z / q^(1/2-h)`

and writes

`tau = 1-t = q(1-eta^2)`.

Eliminating `eta` gives

`1-t = q - z^2 q^(2h)`.

The current repository/user Eq. (4.5) convention elsewhere has been written as

`1-t = q - z^2 q^(-2h)`.

These are not silently identified here. The serialized profile records the mapping state as `pending_explicit_convention_reconciliation`. A future full 3-D Kokuno candidate must first establish the correct convention bridge/native q-solver rather than swapping the exponent by assumption.

## Executable API

`src/openai_ns_reconstruction/kokuno_leading_axis_profiles.py` provides:

- vectorized `U_*`, `H_*`, `W_*` evaluation;
- analytic eta derivatives;
- deterministic SHA identity;
- fail-closed JSON save/load with pinned source, formula, coordinate, parameter-origin and truth metadata.

The tests check formulas, the simplified `W_*` identity, centered-finite-difference derivative agreement, source bounds, positive-`j0` asymmetry/nontriviality, absence of an amplitude-collapse parameter, round-trip identity, and provenance/truth mutation rejection.

## Truth boundary and next integration seam

This is `public_reconstruction_source=true`, but `full_leading_profile_reconstructed=false`, `full_3d_velocity_candidate=false`, `velocity_api_compatible=false`, `pde_validated=false`, `paper_exact=false`, and `openai_field_identified=false`.

The next representation increment should reconcile the native `q/eta` convention (or implement a separately named native Kokuno q-solver) and then assemble only the minimum radial/base-profile information required to lift these axis data into a regular 3-D candidate. Replay against these source formulas will remain a provenance check, not independent Navier--Stokes validation.
