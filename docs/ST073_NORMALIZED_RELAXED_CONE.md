# Normalized relaxed-cone gate for the remote mean fields

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) requires a strict **relaxed** cone on the joined leading profile before Appendix C can modulate the shear into the admissible cone. We now evaluate its equations (4.7), (4.11), (4.15), (4.16), (4.20), and (4.21) on snapshots of our current fields. The profile pressure `Pi` is reconstructed from the inner axis pressure and `Pi_X=E^2/(2X)`; this is not yet the pressure returned by the patched physical field.

At `tau=.0084`, each of the joined baseline, remote two-moment field, and selected two-meridional-shape field was sampled at `eta=0,.2,.35` and `X=6,10,14`. **Each field passed 0 of 9 relaxed-cone points**, hence also 0 of 9 admissible-cone points. The central `X=10,eta=0` values illustrate three distinct failures:

| Field | `a` | `v_s` | `P_c` | upper bound `U(P_c,J_c)` | Failure |
| --- | ---: | ---: | ---: | ---: | --- |
| Joined | 2.009 | 2.009 | 0.180 | undefined | `P_c <= 2` |
| Two-moment | 1.692 | 8774 | 775 | 35.1 | `v_s` exceeds upper bound |
| Two-shape (`odd=-3,even=12`) | 1.692 | 67.2 | -319 | undefined | `P_c <= 2` |

The two-shape streamfunction reduces the normalized shear strength at this point by more than two orders of magnitude relative to the two-moment field, but sends the inviscid vector to the wrong side of the relaxed cone. The baseline's small `P_c` is consistent with the five-moment and axis-pressure mismatch already exposed by `ST073_FIVE_MOMENT_AND_SHEAR_LOOP.md`; that mismatch alone does not prove the numerical `P_c` value. These observations are sampled, not a global impossibility theorem.

This diagnostic uses a fixed-time snapshot of fields built from finite-window coefficients. It does not establish that their normalized profiles are independent of `tau`, meet the exact leading equations, or are smooth on the full `eta in [-1,1]` domain. The reconstructed `Pi` satisfies the leading radial relation by definition, while the present physical patched fields retain the old pressure. Finite differences in `X` and `eta` and piecewise Gaussian radial quadrature are used, without interval certification.

The next constructive step is to choose a pressure datum and a joined `U,E` profile that satisfy all five outgoing moment conditions **and** have a positive relaxed-cone margin on a radial interval. Only then does the paper's high-frequency periodic shear loop and five-moment restoration apply. The current simple remote-bump family is not such an input. Physical finite energy, forcing, complete momentum and volume-L2 acceptance remain open.

Run `python experiments/root_st073/normalized_relaxed_cone.py`; results are in `experiments/root_st073/normalized_relaxed_cone.json` with `accepted: false`.

The later `ST073_REMOTE_PATCH_GEOMETRY.md` screen moves and widens the two-moment patch. Large width narrows one sampled cone gap but does not produce a full relaxed-cone pass.
