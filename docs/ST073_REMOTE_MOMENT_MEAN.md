# Remote radial moment correction for the ST073 mean field

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), Appendix A, imposes radial moment identities on its leading mean profile before stress realization. This experiment applies **analogous normalized slice identities** to the existing corrected `JoinedField`: `I_1=∫(U²-E²/2)dX=0` and `I_2=∫(√(2X)E-H_power)dX=0`, using the bundled heat exterior's analytic infinite-tail quadrature. The full ST073 profile depends on `tau`, so this is a finite-window construction rather than the paper's exact leading profile.

The joined field has its inner/bridge interfaces at `X=3/64` and `X=3/16`; it is pure heat swirl outside the second interface. A fifth-order bump `B(y)=1024y^5(1-y)^5` preserves several radial endpoint jets. The correction is `delta E=a_E(eta)B(X)` and `delta U=a_U(eta)B_X(X)`. The latter has zero radial integral. On each slice, `I_2` determines `a_E` linearly, then `I_1` gives a quadratic equation for `a_U`; the continuous positive root is selected.

Putting the bump immediately beside the core (`X∈[3/64,3/16]`) requires `a_E≈29.6` at `eta=0`, much larger than the original normalized swirl peak near 1.4. A reserved remote heat patch reduces this demand:

| Patch in `X` | `a_E` at `eta=0`, `tau=1/128` | Positive `a_U` root |
| --- | ---: | ---: |
| `[3/64,3/16]` | 29.56 | 0.697 |
| `[1,4]` | 0.300 | 1.322 |
| `[4,16]` | 0.0375 | 2.617 |

The remote `a_U` coefficient grows because its bump is wider, while the derivative `B_X` used in `delta U` shrinks with patch width. The chosen `[4,16]` coefficients vary smoothly with `eta` and negligibly across the three sampled times. Degree-four Chebyshev functions of `eta/0.5` approximate them; this is a fitted finite-window profile, not exact all-scale data.

`MomentMatchedJoinedField` turns the coefficients into a physical axisymmetric velocity. It adds `delta u_theta=√nu q^(-A)a_E B` and a meridional streamfunction `delta psi=nu^(3/2)q^(1-A)a_U B`, with `delta u_r=-psi_z/r`, `delta u_z=psi_r/r`. The implementation differentiates the moving `q`, `eta`, and `X` coordinates when computing `u_r`. Thus the correction is solenoidal by construction. Its pressure remains the joined field's pressure.

At disjoint `eta=±0.15,±0.35` and `tau=.012,.064`, direct integration of the **implemented physical field** gives maximum normalized moment defects `1.15e-5` for `I_2` and `7.77e-4` for `I_1`. Sampled finite-difference divergence in the patch is at most `1.63e-8`. These are moment and incompressibility diagnostics, **not** the requested full momentum threshold. The added mean correction creates physical full-momentum residuals of order `1e3–1e4` at sampled patch points, whereas the unmodified heat exterior there has residual near numerical zero. This residual is a concrete stress target for a later supported oscillatory correction. The field still has no axial localization of the heat exterior and therefore no finite total energy claim; critical-time smooth forcing is also open.

Reproduce in order:

```powershell
python experiments/root_st073/paper_moment_bridge.py
python experiments/root_st073/fit_paper_moment_coefficients.py
python experiments/root_st073/moment_matched_joined_screen.py
```

The scripts write `paper_moment_bridge.json`, `paper_moment_coefficients.json`, and `moment_matched_joined_screen.json` beside them. All retain `accepted: false`. The next step is to evaluate the resulting stress cone and transport equation on this moment-corrected mean field, then realize the required stress without losing the radial identities or the heat exterior.

The subsequent stress-target screen is in `ST073_REMOTE_STRESS_TARGET.md`: the current moment-corrected mean does not yet have an admissible sampled stress cone.
