# Coherent continuous-pressure inner rebuild

The pilot constructs the nonlinear inner core with the continuous angular pressure prefix installed before deriving its pressure Taylor coefficients. It then reconstructs exit/reshape/axial restore/five-bump matching and regenerates the continuous exterior axial inputs from that same source. It does not shift the pressure of an already completed field.

ContinuousAxisPressureJets uses K/Pstar^2 = -2.5 - integral_Rref^Rv Utheta^2/(2 Pstar^2) dlogR at Z=0, with exact pre-Rv (1+Z^2)^-2 separation. Variable angular stages use MP Gauss quadrature; constant stages use exact exponential atoms. The profile may solve a waiting root and replace its initial schedule, so installation targets the actual profile.schedule object.

This is an opt-in source reconstruction: build_source_core and build_candidate accept continuous_pressure=True and pressure_order=192. Existing callers retain the legacy default. The pilot evaluates actual P(infinity), P(infinity)_Z, the pressure transport combination and the rebuilt inner join stress. The full post-Rv pressure tail is deliberately omitted from the dominant axis Taylor jet, and five-bump unresolved input entries remain recorded.

A variable-stage exponential integral agrees with its independent analytic primitive to relative 6.3e-81 at 80 digits. Full source replay, not this narrow check, determines the reported pressure improvement. No claim of exact pressure closure, finite energy, admissible stress or scale recursion is made.

Initial actual replay at Z=.3 gives P(infinity)=2.15130102858891e-69 and P(infinity)_Z=-2.36840480211623e-69, with pressure-transport coefficient 3.44602898707912e-69. The absolute pressure defect is 1.7796940893e-67 times the legacy value. This is a nominal single-Z improvement; unresolved inner input entries 3 and 5 remain and the omitted post-Rv axis jet is not restored.

Use `build_joined_field(continuous_pressure=True)` and pass that source to ContinuousIncomingOuterField. The shared profile detects this source provenance and automatically selects MP pressure nodes with order at least the axis anchor order. The legacy default remains available for controlled comparisons. The pilot now replays the coupled five-moment stress at flattening and heat with this public path.

The public source/quadrature route replay completes successfully and preserves the same terminal pressure improvement. The coupled total axial stress is about 9.8244e-97 of the legacy value at Rv+50 and 1.7797e-67 at Rtail+4. Total angular stress is unchanged at saved precision (ratio 1). These compare the actual rebuilt source with the legacy binary64-96 pressure evaluation; they are not stress-cone or PDE residual certificates. This leaves angular/mixed-moment matching as the next principal stress task.
