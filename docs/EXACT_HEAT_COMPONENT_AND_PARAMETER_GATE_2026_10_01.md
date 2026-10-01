# Exact heat component and unresolved parameter admission (F35)

The selected-family outer candidate now has a specified exact heat component,
a flat collar interpolation, and directed angular-moment, pressure-moment and
swirl-energy tail defects. This closes the heat-component calculation only;
it does not match the incoming five moments or certify a complete background.

## Implemented component

For `a=delta/2`, the exterior uses the actual Gamma integral

`H_delta(xi)=Gamma(1+a)^(-1) integral_0^infinity exp(-v) v^a (1+xi*v)^(-a) dv`.

The construction retains its positive deficit instead of setting `H_delta=1`.
With `S=1/Rtail`, `t=log(R/Rtail)` and `xi=2(1-Z^2) S exp(-t)`, it records
`1-H_delta=a*S*Dhat` and whole-axis first-derivative enclosures. The collar
equals the preheat candidate at its inner endpoint and the exact heat profile
at its outer endpoint. Its definition retains the formal selected radius and
amplitude because the physical scales cannot be expanded numerically.

The three returned heat differences are **scaled** quantities. Their physical
factors are respectively `sqrt(2)*c_infinity*Rtail^(1-a)*S`,
`c_infinity^2*Rtail^(-1-delta)*a*S`, and
`c_infinity^2*Rtail^(-delta)*a*S`. The positive bound `exp(-1000)` is an upper
bound on `S`, never its value. Consumers must preserve these exact/formal
factors. Backward targets include the separate epsilon and epsilon-squared
preheat atoms with their original signs.

The companion checker verifies nine independent Gamma/PDE/moment identities,
positive defects, exact axis endpoint values, collar endpoint sources and
three finite Gamma-integral reference examples. The scope is recorded in
`experiments/root_st073/lei_ren_part1_paper_shared_exact_heat_component.json`
and its `_check.json` companion; the corresponding Python files reproduce
the calculations.

## Parameter discrepancy requiring a new source

The supplied Lei--Ren Part I v2, Section 7.1, requires `c_epsilon<=0.001`.
The existing logarithmic source uses `epsilon=0.01*delta`. No existing receipt
proves an extension of Section 7 to that value. The isolated heat component
remains valid under its own weaker epsilon conditions; it cannot certify
the full Section 7 construction.

Preserve the old source, hashes and completed finite coefficients. The next
step is a separate `c_epsilon=0.001` source with its own waiting equation,
fourteen pressure atoms and fingerprint. Quantify its pressure difference
from the old source, then establish fixed-point/finite-core transfer or
regenerate the affected core. A new fingerprint alone is not a transfer proof.

## Remaining dependencies

1. Admit the compliant source and its core/transition family.
2. Solve actual angular and pressure corrections, using the exact heat targets.
3. Select the axial energy parameter from the complete corrected future tail.
4. Assemble the matched outer field and prove its cone margins.
5. Construct admissible stress/flat remainder, genuine order-dependent temporal
   recursion, oscillatory corrections and independent Cartesian residuals.

None of these remaining dependencies is claimed complete by F35. In
particular, the finite 144-order radial core calculation is distinct from
the requested temporal coefficient recursion.
