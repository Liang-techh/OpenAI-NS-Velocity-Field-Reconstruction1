# Source exit after the initial collar

The Section 9.23 comparison is frozen in R after y=2hb. Its velocity
values and Z jets remain fixed; all five moments and pressure continue
with their actual integral powers. Queries beyond the finite core no
longer ask the core polynomial to extrapolate.

The frozen driver has exact algebraic coefficients:

```
D(R) = d1 R + d0 + dm/R
I_z(R) = c3 R^(3/2) + c1 R^(1/2) + cm R^(-1/2).
```

`ExitContinuation` uses these coefficients to integrate the post-collar
prescribed shear's leading changes and bound the changes in velocity.
Writing Re for the collar endpoint, its log-F primitive is
`d1*(R-Re)+d0*log(R/Re)+dm*(1/Re-1/R)`.
The axial primitive of sqrt(R/2)*I_z is
`c3*(R^2-Re^2)/(2*sqrt(2))+c1*(R-Re)/sqrt(2)+cm*log(R/Re)/sqrt(2)`.
The tiny positive epsilon remains explicit in these increments and in
the actual prescribed shear passed to the stress evaluator.

For values/moments at finite precision the implementation uses analytic
frozen continuation only when conservative conditional velocity bounds
are below 10^(-precision+20), including a relative Uz guard. It also
reports corresponding bounds for the five moment values. If the guard
fails, it requires integration of the full ODE rather than silently
accepting the approximation. Bounds are conditional on the supplied
comparison and endpoint; endpoint RK, pressure approximation and Z-jet
uncertainties are not covered. No corrected-tail derivative certificate
or full source constant choice follows from them.

Run `python experiments/root_st073/lei_ren_part1_paper_exit_continuation.py`.
The receipt samples the exact H0 root and Z=.3 at R=1,10,100, recording
velocity/moment bounds and actual relaxed/admissible cone tests. The
comparison's independent frozen-extension receipt covers R=110 too.

The actual Step 2 continuation stops at R=100. The subsequent short
shear changes on 100..110, restoration of the reference power law,
Section 10 moment repair, heat/exterior joining and time-scale recursion
remain necessary. This enlarged exit is not a localized finite-energy
field. It can supply `LocalCoreExitField` as its radial provider, while
retaining those domain and approximation limits.
