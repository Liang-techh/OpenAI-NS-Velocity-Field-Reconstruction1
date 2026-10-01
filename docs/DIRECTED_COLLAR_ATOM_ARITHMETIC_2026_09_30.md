# Directed collar atom arithmetic — 2026-09-30

The finite collar equations retain independent pressure and width coefficients. `IntervalPressureWidthJet` now provides this rectangular finite ring with directed `IntervalDifference` coefficients. Each atom retains a nominal interval and a separate perturbation interval; operations do not evaluate the formal parameters at (1,1).

Supported operations are addition, subtraction, multiplication, reciprocal, division, integer powers, exp, log, and sqrt. Domains require the constant coefficient to satisfy the applicable nonzero or positive conditions. Context and truncation orders must agree. The normalized convention is the coefficient of pressure^p width^w, including factorial normalization.

For nonlinear series the constant atom is removed structurally before forming the nilpotent remainder. Subtracting an interval from itself is insufficient: it can leave a constant interval, for which a finite nilpotent expansion would not be valid. The retained nonconstant polynomial has total nilpotence order at most pressure_order+width_order+1.

Shared perturbations use cancellation-aware formulas: exp difference=exp(nominal)*expm1(difference), log difference=log1p(difference/nominal), and sqrt difference=difference/(sqrt(total)+sqrt(nominal)). Exact zero perturbations remain zero.

The bounded fixture checks coefficient containment against the legacy resolved finite ring, independent resolved coefficient perturbations, mixed/tiny atoms, wide nominal intervals, exact zero perturbations, and invalid domains. Its JSON reports actual check counts. These are arithmetic certificates for the retained finite rectangle, not certificates for omitted pressure/width orders, ODE integration, original source coefficients, or the physical velocity field.

The companion second-Z interval ring is the next adapter required before collar equations can propagate the value, first Z derivative, and second Z derivative using this arithmetic. The legacy second-Z ring explicitly coerces to ordinary MP pressure/width jets, so it cannot silently stand in for this adapter.
