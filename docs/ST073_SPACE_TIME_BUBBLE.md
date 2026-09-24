# Two time elements and an interior temporal correction

The compact patch was marched to three nodes with `tau` step `2.5e-6`, creating two Hermite time elements. Run from the repository root:

```powershell
python experiments/root_st073/curl_wave_patch_evolution.py --time-step 2.5e-6 --stages 3 --output-name curl_wave_patch_evolution_3stage.json
python experiments/root_st073/curl_wave_patch_trajectory.py --evolution-name curl_wave_patch_evolution_3stage.json --interval-index 1 --output-name curl_wave_patch_trajectory_second_interval.json
python experiments/root_st073/curl_wave_temporal_bubble.py
```

The Hermite potential and pressure coefficients agree at the shared interface to `2.1e-11` in the recorded calculation. Its endpoint slopes match by construction. This gives an exactly divergence-free callable velocity with a continuous potential time derivative and pressure across the two elements. It does not imply small momentum residual.

On the second element, the 16-node held-out midpoint full-momentum maximum is `1,748,646`, compared with `727,061` for the original field. The third node's instantaneous derivative/pressure projection is `573,667`, showing again that node fits do not control the interval interior. The first harmonic's fitted slope norms at the three nodes are `4,015`, `13,445`, and `45,178`.

To correct the interior without disturbing either endpoint, `curl_wave_temporal_bubble.py` adds a compact odd time polynomial to the exact-curl potential and an even polynomial to pressure. Both have zero value and first time derivative at the element endpoints. At the midpoint, the potential bubble has zero value but a freely fitted time derivative; the pressure bubble equals its fitted pressure. A 25-node spatial fit reduces the **direct** 16-node held-out midpoint maximum from `1,748,646` to `463,797`, below the original field's `727,061`. The endpoint coefficients are unchanged exactly in the report.

This single-time correction fails away from its fit time. On four spatial nodes held out from fitting, the quarter-point maxima are `1,565,617` and `1,762,787` after the bubble, against original-field values `579,224` and `579,660`. Therefore the second element remains rejected. The result identifies a concrete next discretization: fit temporal basis functions at several interior times **simultaneously**, with the full nonlinear residual and pressure, while preserving the interface conditions. That is closer to the coupled amplitude/pressure and mean-flow evolution required in Sections 7--9 of the [OpenAI Navier--Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), but the paper's estimates have not been established for this reconstruction.
