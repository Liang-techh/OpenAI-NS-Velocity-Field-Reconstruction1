# Candidate angular inlet trace from the original equation

The candidate inlet diagnostic now uses the same refined coupled coefficients and computes the angular stress trace through the exact integrated angular equation. It does not prescribe or fit a zero trace on the retained polynomial.

Write s=Lambda R, F=F0 Phi, g=-F0_Z/(Lambda F0), and let

    Q = epsilon/L * [W(s Phi_s+Phi)+H Phi_Z
                     +delta/2*(Phi-2Z Phi Uz)] - g/L * H Phi.
    D = 2(s Phi_ss+2 Phi_s)-Q.

Then the finite retained field has

    Ttheta/F0 = (1/s) integral_0^s sigma D(sigma) d sigma
              = 2s Phi_s - (1/s) integral_0^s sigma Q(sigma) d sigma.

Both expressions are evaluated with directed interval arithmetic and their enclosures must overlap. For polynomial D=sum D_k s^k, the first expression is sum D_k s^(k+1)/(k+2). W and H include the full coupled axial polynomial, and Phi_Z is the gauge derivative before adding the amplitude logarithmic derivative already represented by g.

The default exit is s=4, within the certified analytic radial rectangle up to4.1. The diagnostic covers Z=0.3 only. At radial order56 it gave a normalized finite trace Ttheta/F interval approximately[-1.75612e-63,1.75612e-63]. This is consistent with the exact analytic core equation, whose angular trace is zero. It is not a matching certificate: containment of zero does not show exact cancellation in the retained polynomial, and the collar has not been regenerated on the new candidate data.

Reproduce after a bounded radial segment has stopped:

    python experiments/root_st073/lei_ren_part1_paper_candidate_inlet_trace.py

After order124, combine the full coefficient/tail budget, update this trace, and regenerate shared exit, pressure and collar data. The remaining whole-axis bottleneck includes the weighted sign of Phi_s near the axis anchor where chi is small; the generic analytic correction norm alone does not settle that sign.
