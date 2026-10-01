# F34: later angular candidate and actual preheat waiting source

The same-family angular candidate now continues from the Rp inlet through the pulse interval, exact 100-unit flatten, post-flatten power buffer, steep transitions/power and waiting interval. Its angular cumulative primitive and formal logarithmic amplitude are callable. This does not assemble the complete velocity/moment field after Rv: axial coefficients, angular/pressure corrections and the exact heat exterior are still pending.

## Transported angular moment

The API in `experiments/root_st073/lei_ren_part1_paper_shared_outer_angular_candidate.py` uses X=Mtheta/(sqrt(2)*R^1.5*Utheta). The physical moment definition gives X_y=1-(3/2+partial_y log Utheta)X. All stages derive from this equation, using the same F32 inlet. The axial pulse changes no angular moment or swirl, so its unselected amplitude does not obstruct this candidate calculation.

On the pulse interval, Xv=1/(1-mu)+(Xp-1/(1-mu))*exp(-13*(1-mu)/mu). The tiny actual memory has a proved nonzero exp(-1000) magnitude cap; its exact source definition is retained. Xv is independent of Z because swirl and angular primitive share the same (1+Z^2)^-1 factor before flattening.

On flattening, t=log(R/Rv) in [0,100], F=(q/2)^sigma(t/100), q=1+Z^2, and

    X(t,Z) = [Xv + integral_0^t exp((1-mu)s) F(s,Z) ds]
             * exp(-(1-mu)t) / F(t,Z).

Closed cells integrate the exponential weights exactly and retain directed axial jets of F. The post-flatten buffer retains the inherited axial dependence of X even though swirl is Z-independent. The steep-in integrating factor has rate (1-mu)*(1-sigma); the steep-power rate is zero; steep-out has rate (1-delta/2)*sigma. The waiting rate is k=1-delta/2. No moment history is reset.

## Original continuous preheat waiting equation

The actual raw X_t(0) and directed collar integral J=integral_0^3 exp(k*t)[1-sigma(t)+sigma(t)*f((3-t)/2)]dt determine the same source root

    tau = [log(X_t(0)-1/k) + log(1-epsilon)
           - log(epsilon) - log(1/k+J)]/k.

The exact equation is unchanged, and the root enclosure lies inside the earlier Md40 source enclosure. Its width improves from 3 to approximately 0.00630990. The small log(1-epsilon) correction and positive departure X-1/k at the end of waiting are stored separately. The preheat axis pressure source/datum hashes remain unchanged; this refinement does not replace H_delta by one as an exact heat field.

## Formal amplitude, evidence and remaining work

The angular log origin at Rv is separated into log(Utheta(Rp,0)/Pstar), the inverse-mu term -13/(2mu), and the finite offset -13. Relative stage profiles and coordinates are stored separately. Consumers must not collapse these into rounded absolute log origins and subtract them later; that would lose the unit transitions and heat collar.

The companion `shared_outer_angular_candidate_check.py/.json` independently verifies seven physical/angular/waiting identities, twenty interfaces in value and retained first axial derivative, whole-axis C1 calls, the strict waiting refinement and the separately retained positive waiting departure. All-order smoothness follows from the original flat cutoff; all-order jets are not numerically evaluated here.

Next: build the selected-radius H_delta heat factor/collar/exterior with positive deficit and tail data; close angular and pressure moments; insert the complete corrected swirl-energy tail into F33 and select ap; then install the actual coupled pulse and full matched background. Exact heat matching, whole outer cone, global admissible stress/flat remainder and n-dependent temporal recursion remain unfinished.
