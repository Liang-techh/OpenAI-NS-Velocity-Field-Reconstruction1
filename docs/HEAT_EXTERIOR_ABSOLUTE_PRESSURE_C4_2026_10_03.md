# Conditional pressure-tail candidate in the Gamma exterior

The new pressure-tail evaluator computes axial Taylor coefficients through five and log-radius/Z mixed derivatives through total order four on Z in [-1,1], t=log(R/Rtail)>=3. Integral and derivative checks pass. Transfer of the admitted absolute pressure closure into this C4 provider is NOT yet proved.

Use CompliantHeatPressureC4(allow_conditional=True).exterior(Z,t) only for candidate work. Default construction rejects the unresolved transfer. Original velocities, pressure datum, repair coefficients and forward-pressure histories are retained. The original dispatcher and collar remain unchanged.

## Implemented algebra

The earlier CompliantAbsoluteMomentClosure proves its own original infinity offset is zero from P0=-Mp_pre(infinity) and Delta_bump=Delta_heat. Matching family/source hashes ensure current dependencies; neither hashes nor overlapping intervals establish the transfer to another provider.

For a=delta/2, p=1+delta, S_current=S*exp(-t), K=H_delta, the candidate uses the full infinite Gamma integral:

```
xi = 2*(1-Z^2)*S_current
A_p = integral_0^infinity exp(-p*v)*H_delta(xi*exp(-v))^2/2 dv
C = (Ev0^2/Pstar^2)*theta_base^2
P_candidate/Pstar^2 = -C*exp(-p*t)*A_p
P_candidate/Utheta^2 = -A_p/K^2 < 0
d_t(P_candidate/Pstar^2) = C*exp(-p*t)*K^2/2
```

Orders one through four follow from FTC and the exact Leibniz rule; axial data use the full Gamma Taylor composition. No convergent inverse-radius series, pressure fit or cap endpoint defines the candidate function.

The partial production bridge checks 33 source-expression and symbolic identities, including theta_base=exp(Ltail-log(1-epsilon)) and Prel*Ptail*tailmult=theta_base^2. These prove formal unit algebra. The current P3/infinity-offset identities still use placeholders and do not establish the actual retained-history transfer.

## Outstanding source bindings

Independent read-only review identified four required bindings:

- Bind both inlet buffer objects to the same SharedOuterBuffer callable, beyond matching power() call syntax.
- Bind exact S=exp(-(logRref+13/mu+tail_finite))=1/Rtail and its Gamma argument independently of all numerical caps.
- Bind actual C4 Ptail and forward_pressure() integrals to the closure's original P0, Mp and Prv histories.
- Explicitly bind the C4 Gamma derivative evaluator to the canonical full positive Gamma expectation.

Until these hold, source_history_transfer_conditional=True, absolute_pressure_same_source_mixed4_available=False and checker all_passed=False. The original closure's zero infinity offset is not transferred merely by replacing C4 pressure with the negative tail.

## Focused evidence

The checker recomputes the earlier absolute source proof, partial source/unit identities and four pressure derivative identities. A finite-parameter fixture compares the complete hyperu function with direct positive Gamma quadrature, integrates the whole pressure tail and its first axial derivative, and checks seven mixed derivatives independently. These exercise candidate algebra and units; they do not certify the actual source transfer.

Five actual packets cover three points, the whole-Z t=3 inlet and the unbounded exterior. Their 75 forward/tail mixed-bound overlaps and earlier C1 agreement are diagnostics only. At t=3 the original whole-Z forward pressure interval is about [-2.48627,2.48717], while the candidate tail enclosure is approximately [-3.708755e-18400774051555611784,0]. A broad forward sum does not establish a physical mismatch; the zero cap endpoint also does not define zero physical pressure.

```
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4.py
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4_check.py
```

A successful candidate run reports candidate_checks_passed=True and all_passed=False.

## Next work

Complete the four source/history bindings first. Then derive the independent exact heat PDE and moment-normalized exterior stress identity, tighten collar pressure and integrate the admitted companion into the main dispatcher. All-region stress/cone/flat remainder, physical energy, coupled n=1, n>=2 recursion, oscillatory corrections and the independent Cartesian residual remain open.
