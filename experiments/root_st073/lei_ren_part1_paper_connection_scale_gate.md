# Necessary connection input rejection

`connection_scale_gate.py` applies necessary conditions from source
Sections 9.4, equations 9.4,9.16,9.17. From the supplied axis values,
log(Cstar F0)=-Lambda G, so A>=10+Lambda*max(|G(-1)|,|G(1)|).
Also logK>=logCstar+Lambda*max(|G(-1)|,|G(1)|), using the inverse-F
norm at the same axis endpoints. These lower bounds can reject a
candidate. They cannot certify its C3 norms or the whole connection.

The current fixture has logCstar=165.786, logRref=1802.562, but
A>=2.829e36. The long-shape restriction alone requires
logCstar>=1.132e38. The fixture also fails the axial-radius condition.
From Kp>=1, source epsilon0<=1/(1e6*2001), hence its prescribed
j<=epsilon0/32<=1.562e-11; fixture j=.02 fails this necessary condition.
The short-collar condition hb=cstar*K^-100 also fails with hb=.005.
These explicit failures retain the fixture as an algorithm-development
candidate; they forbid claiming a supplied source-compatible pair.

The complex bound A_Omega=max(-Re G) over the complex neighborhood
is distinct from real max(G). No substitution between them is made.
The actual C3 bounds, pressure bound Kp, moment constants and the
equality hb=epsilon remain unresolved inputs.

`build_source_core` now accepts keyword j,logC,logPstar,delta as one
shared configuration. It recomputes Rref=110*(Cstar Pstar)^10, the
corrected outer profile, axis-pressure anchor, axis data and core factory.
It returns shared_parameters and pressure_recomputed_for_shared_parameters.
The nondefault wiring receipt uses j=.01,logC=170,logPstar=13,delta=1e-34
and checks matching axis/core delta and a positive decreasing local F.
That test does not claim the changed candidate satisfies source gates.

The wiring check can be replayed with `experiments/root_st073` on
PYTHONPATH:

```python
from lei_ren_part1_paper_core_adapter import build_source_core
b = build_source_core(precision=100, degree=3, Lambda="1e36",
                      j=".01", logC="170", logPstar="13", delta="1e-34")
assert b["axis"].delta == b["core"].delta
assert b["pressure_recomputed_for_shared_parameters"]
```

Before accepting a replacement, supply shared constants, control A/K,
check all radius inequalities, and rerun the pressure/core/connection
construction. Old receipts belong to their recorded fixtures and must
not be attached to a different candidate.

## Replacement work order

1. Bound the same outer profile's pressure C3 jets and establish its Kp
   input; choose the Section 10 moment constants before choosing j.
2. Choose j from the source tolerance. Recompute H0,G and the axis root;
   choose Lambda with both j^-2 and observed pressure/narrow-layer scale
   restrictions included. Preserve historical failure probes.
3. Bound A in the mixed core coordinates and the frozen R<=110 interval,
   retaining all derivatives up to order three. A real endpoint lower
   bound cannot stand in for this upper bound or the complex A_Omega.
4. Choose logCstar and logRref simultaneously with source Eq9.17; choose
   delta consistently with the outer restrictions. Use the shared builder
   to recompute pressure and core, then repeat the bounds for that pair.
5. Bound K, choose one common hb=epsilon from its restriction, and check
   the full radial ordering in log coordinates. If a collar is below
   physical-coordinate precision, preserve its switch phase explicitly;
   do not round it to a zero-width switch or set epsilon to zero.
6. Replay the core, exit and short switches for the accepted shared pair,
   keeping quadrature, derivative and approximation errors separate.
7. Only then append the long angular reshape, axial restoration and
   actual five-moment repair. Measure stress and full-field diagnostics
   with the same pressure/velocity identity.
