# Actual directed finite-core collar endpoint — 2026-09-30

The accepted finite core at Z=.3 has now been evaluated with directed 473-digit interval arithmetic. All eight first- and second-width endpoint coefficients, all five unnormalized profile moments, F/Uz/P, and Ur/Ur_Z are recorded with exact interval endpoint tuples. Original P0 is preserved. No source core, transition RK chain, or R100/R110 continuation is regenerated.

This encloses arithmetic on the stored finite core input. It does not enclose the original core-series approximation, source-parameter errors, pressure truncation, omitted width orders, or the full nonlinear collar solution. It is a point-Z endpoint result, not a global axial-function or stress-cone certificate.

## Actual second-width coefficients

For display, the pressure-power-zero value slot of each normalized coefficient, divided by h_b^2, is:

| State | Display value |
| --- | ---: |
| theta | 5.0391693023633759 |
| mz | 2.40000000000002 |
| mixed | 6.04700316283610148 |
| axial | 2.880000000000048 |
| swirl | 1.0391693023633759 |
| p | -0.960830697636624096 |
| g | -2.94333839527565701 |
| u | 2.10531208355682245e-23 |

These are individual formal coefficients, not complete physical velocities, residuals, or blow-up measurements. The receipt retains every pressure atom and all three axial derivative slots, and performs the unnormalized moment/field conversion before exposing Ur and Ur_Z. The actual h_b remains separate.

## Independent radial derivatives and legacy comparison

Independent MP differentiation of the same finite scalar core at pressure power zero verifies R*d/dR of I_theta, Iz, and D inside their directed intervals. Review corrected the I_theta radial transport differentiation to F_R*transport+F*transport_R; the Utheta factor includes sqrt(2R) and must first be canceled against the original denominator. Final results use the corrected formula.

The old 260-digit inlet calculation is audited in 180 comparisons: six fields, three axial derivative slots, ten pressure powers. 174 stored scalar outputs lie outside the much narrower new intervals. The receipt records all scalar differences rather than enlarging the intervals to accept them. This comparison does not prove a physical field error of that count or establish which approximation errors dominate. It confirms that the older rounded output is not a substitute for the new arithmetic enclosure. The original finite source data and old caches remain unchanged.

## Portable replay

`lei_ren_part1_paper_interval_collar_core_inlet_input.json` contains exact binary tuples for Lambda, delta, and every cached F/Uz/P coefficient, plus the original trusted pickle SHA256. Every atom survives the JSON loader exactly even when called at the default MP precision. It is a portable copy of stored finite input, not an independently certified source construction.

Run:

```
python experiments/root_st073/lei_ren_part1_paper_interval_collar_core_inlet_check.py
```

The committed JSON input is selected by default. No pickle or long source rebuild is required; a replay takes roughly35 seconds on the current host. The output is `lei_ren_part1_paper_interval_collar_core_inlet_check.json`. The export script is only for regenerating the portable input from the trusted local accepted pickle.

Next dependencies are original pressure/core/source error propagation into these inlet atoms, omitted-pressure/width bounds, extension from the endpoint to the whole collar and axial domain, and functional moment matching/stress admissibility. Ur_ZZ still needs third-Z moment data. True temporal n-dependent recursion remains a later distinct layer.
