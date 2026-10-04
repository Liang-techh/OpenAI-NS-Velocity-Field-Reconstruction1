# Original pulse-end / flatten functional interface — 2026-10-04

The actual pulse `s=0` / flatten `t=0` similarity interface now connects the full five-moment data through total mixed order four, theta/axial stress through total mixed order three, and the same absolute analytic pressure through total mixed order four. Ordinary derivatives on both sides are derivatives in `log R`; the proof covers the original `Z ∈ [-1,1]` as functions, rather than selecting sample values or relying on overlapping intervals.

## Current implementation

`experiments/root_st073/lei_ren_part1_paper_compliant_pulse_end_flatten_join.py` consumes the accepted pulse-end stress, flatten stress/cone, and original inlet sources with current hashes and a common source family. The new controller stage is `pulseendjoin`.

The local source bridge replays the actual pulse `data()` and power `incoming()` jet producers. It compares the `U/(1+Z²)` and `Pin/(1+Z²)²` coefficients through axial order five, including the actual reciprocal-square recurrence. It also replays the live buffer endpoint at `Z=0` and the serialized final endpoint at `Z=1/2`, the two constructors' pressure-constant assignments, and the fifth provider's serialization of the live fourth-provider constants. The original analytic datum's source and enclosure hashes agree across the five-moment repair, outer buffer and inlet routes.

At the interface, the actual beta routines and backward integrals have empty future support. All meridional pulse values and the required derivatives therefore vanish. The signed angular memory is retained as `1/(1-mu)+(Xp-1/(1-mu))*H`, with the exact logarithmic definition of positive `H`; it is not replaced by its equilibrium value. Both energy routes retain half of the same complete future integral, including selected angular corrections and the entire heat tail.

The pressure bridge replays both original forward histories and the pulse stress companion's actual extraction of the flatten canonical pressure. The amplitude bridge and production-radius equality are consumed before comparing the different reference rates `1/2+mu` and `1/2+a`. No original radius, support, source parameter, velocity, pressure datum or selected coefficient is changed.

The functional proof replays the actual flatten shape, derivative recurrences, full pulse stress sectors, downstream stress operator and pressure derivative operator. Exact binary source halves are represented as rational halves during symbolic algebra to prevent floating-point polynomial expansion from introducing cancellation noise. This changes only their exact algebraic representation. It rationalizes no parameter, enclosure or field value.

The focused checker recomputes the source proof and consumes fresh hashes. Its independent fixture evaluates the original full moment and stress formulas with unequal reference rates, nonzero signed incoming memory and an analytic datum. It checks 435 physical similarity derivatives on both sides, through the interface's stated orders, with tolerance `1e-55`. This fixture is independent interface evidence, not a project-wide physical residual or cone certificate.

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulseendjoin
```

## Next work

1. Build the pulse-end physical adapter with all meridional terms, original viscosity/time units, the completed tensor and diagonal, and the actual three-component remainder. Compose the physical interface using the admitted similarity moment/stress/pressure connection.
2. Establish continuous pulse-end admissibility through both original supports, preserving signed correlations and exact log factors. Complete the required derivative-interface bounds at the four internal beta support endpoints.
3. Extend the companions through the original gap, main pulse and entrance with their actual source histories and both-sided joins.
4. Resolve upstream finite-width bridge feedback and complete the global tensor cone, flat remainder, volume norms and required-domain energy bounds.
5. Implement true `n=1` and `n>=2` coefficient recovery with moment repair and smooth summation, followed by oscillatory correction and independent corrected Cartesian residual/dynamics.

This closes one original similarity interface. Pulse-end physical decomposition, its cone, the whole pulse, completed global tensor admissibility, required-domain energy, temporal coefficient recursion and oscillatory correction remain open.
