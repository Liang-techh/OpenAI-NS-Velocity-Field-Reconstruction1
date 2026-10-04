# Source-bound absolute pressure in the Gamma exterior

The pressure evaluator computes axial Taylor coefficients through five and log-radius/Z mixed derivatives through total order four on Z in [-1,1], t=log(R/Rtail)>=3. The defining-function and actual retained-history proofs now transfer the admitted original absolute pressure closure to this C4 companion. Producer and independent checker pass.

Use CompliantHeatPressureC4().exterior(Z,t). Original velocities, pressure datum, repair coefficients and forward-pressure histories are retained. The original dispatcher and inward collar have not yet adopted this companion.

## Implemented algebra

The earlier CompliantAbsoluteMomentClosure proves its own original infinity offset is zero from P0=-Mp_pre(infinity) and Delta_bump=Delta_heat. Matching family/source hashes ensure current dependencies; neither hashes nor overlapping intervals establish the transfer to another provider.

For a=delta/2, p=1+delta, exact S=1/Rtail, S_current=S*exp(-t), K=H_delta, the pressure uses the full infinite Gamma integral:

```
xi = 2*(1-Z^2)*S_current
A_p = integral_0^infinity exp(-p*v)*H_delta(xi*exp(-v))^2/2 dv
C = (Ev0^2/Pstar^2)*theta_base^2
P/Pstar^2 = -C*exp(-p*t)*A_p
P/Utheta^2 = -A_p/K^2 < 0
d_t(P/Pstar^2) = C*exp(-p*t)*K^2/2
```

Orders one through four follow from FTC and the exact Leibniz rule; axial data use the full Gamma Taylor composition. Numerical caps enclose the exact function; they do not define its radius or normalization.

The production bridge retains the amplitude/unit identities, including theta_base=exp(Ltail-log(1-epsilon)) and Prel*Ptail*tailmult=theta_base^2. The new retained-history proof verifies 76 actual source/algebra bindings. It follows the original pulse, flattening, power, switch, waiting and collar integrals, with their original datum and units. Identical full integrals remain exact symbolic atoms in additive identities; the previous P3 placeholders have been removed.

## Completed source bindings

All four links identified by independent review are now verified:

- Both inlet objects resolve to the same SharedOuterBuffer.power callable, through their actual import/constructor paths and C4 constants chain.
- Exact S=exp(-(logRref+13/mu+tail_finite))=1/Rtail and its Gamma argument are bound independently of numerical caps.
- Actual C4 Ptail and forward_pressure() integrals are bound to the closure's original P0, Mp and Prv histories, including the full collar/exterior future integral split at t=3.
- The C4 Gamma derivatives and infinite-tail evaluator are bound to the canonical full positive Gamma expectation and the admitted angular high-jet source proof.

The combined bridge has no unresolved bindings. Receipts now report source_history_transfer_conditional=False, absolute_pressure_same_source_mixed4_available=True and all_passed=True. The original zero infinity offset transfers through these exact source/history identities. Runtime callable paths are normalized; receipts store portable source filenames.

## Focused evidence

The checker recomputes the earlier absolute source proof, both new bridges, source/unit identities and four pressure derivative identities. A finite-parameter fixture compares the complete hyperu function with direct positive Gamma quadrature, integrates the whole pressure tail and its first axial derivative, and checks seven mixed derivatives independently. The exact bridges establish source transfer; finite fixtures supplement them.

Five actual packets cover three points, the whole-Z t=3 inlet and the unbounded exterior. Their 75 forward/tail mixed-bound overlaps and earlier C1 agreement are diagnostics only. At t=3 the original whole-Z forward pressure interval is about [-2.48627,2.48717], while the candidate tail enclosure is approximately [-3.708755e-18400774051555611784,0]. A broad forward sum does not establish a physical mismatch; the zero cap endpoint also does not define zero physical pressure.

```
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4.py
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_pressure_C4_check.py
```

The current run reports candidate_checks_passed=True and all_passed=True. Actual exterior stress, global stress and temporal recursion flags remain false.

## Next work

The exact heat PDE and conditional terminal-normalization theorem are implemented; see HEAT_EXTERIOR_STRESS_EQUATIONS_2026_10_03.md. Next transfer actual angular/energy terminal histories to eliminate both homogeneous stress constants, tighten inward collar pressure and integrate the companion into the main dispatcher. Finite-width bridge feedback/second switch, all-region stress/cone/flat remainder, physical energy, coupled n=1, n>=2 recursion, oscillatory corrections and the independent Cartesian residual remain open.
