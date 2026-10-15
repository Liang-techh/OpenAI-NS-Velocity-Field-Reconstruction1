# Actual first-micro modified velocity enclosures — 2026-10-15

## What is now computed

`OriginalFirstMicroNumeric` extends numerical leading+n=1 velocity coverage
from r_minus through the first actual modified microscopic chart, ending at
native phase s=1. It encloses the original modified function using a numerical
background proxy and separately certified signed-source correction errors.
The fast phase and the signs of the corrections are not numerically resolved.
Neither a magnitude cap nor a freely selected angle defines the field.

```python
from openai_ns_reconstruction.original_first_micro_numeric import OriginalFirstMicroNumeric
field = OriginalFirstMicroNumeric()
left = field.selected('.3', '.5', log_tau='-2', azimuth='.25')  # s=sc/2
collar = field.selected('.3', '.75', log_tau='-2')             # exact flat collar
active = field.evaluate('.3', '.625', log_tau='-2')           # actual native s
exit = field.evaluate('.3', '1', log_tau='-2')
```

`evaluate` takes actual native s in [original_sc/2,1]. `selected` takes an
exact source multiple of sc in [1/2,1], avoiding a decimal approximation to
the original 798-bit sc singleton. These are different from the preceding
n=1 cutoff fraction phase/(sc/4). The output label
`source_cover.whole_original_phase_cover=[0,1]` denotes the normalized fast
periodic phase cover. It does not enlarge the native source domain below sc/2.

The returned normalized leading Phi/Psi/own mean have native Taylor Z0..2;
Qpsi has Z0..1 with the leading recovery weight1-delta. The physical packet
returns cylindrical and Cartesian leading+n=1 velocity plus its shared
streamfunction and factored coordinates. It requires native coordinates;
no new arbitrary Cartesian inverse or physical spatial gradient is supplied.

## One exact original N, no smaller substitute

The exact selected source integer is

`N = 2^(2^1358356628656378313)`.

The adapter checks the canonical repair-band integer node, its selected J,
its N>=N0 proof and same-family range/source receipts. It computes logN by
shifting the binary exponent of a directed log2; neither 2^J nor N is
constructed. The accepted range theorem is uniform for all integer N>=N0.
The original common N is retained across the field and correction histories.
Its repair-only selection does not admit missing global derivative, cone,
heat or interface frequency requirements.

The unchanged Section11 source is

```
Phi_N=Phi*exp(A/N)
U_N=U+Pstar*B_over_Pstar/N
phase=fractional_part(N*hb*(s-sc/2))
E=sqrt(2R)*F0*Phi/Pstar, V=U/Pstar
```

Consequently normalized `delta_Psi=Pstar*B_over_Pstar/(N*epsilon)`.
No extra Pstar or F0 dresses this correction. A/N and B/N caps use the
actual complete original phase cover, retaining the exact defining signed
C3/C5 source bindings. This gives a certified approximation to the original
modified function without resolving an impossibly large spatial frequency.
It does not provide the signed loop value at a resolved fast phase.

For each ordinary axial row0..2, the actual source-log bound is proved less
than log(1e-1000) before a finite numerical error bound is used. That number
is an error cap, not a selected correction value. Taylor coefficients divide
ordinary derivative bounds by k! exactly once. For |A/N|<=1, exp(A/N)<=3
bounds the full value/first/second exponential increment, including the
squared first-derivative term. Phi's actual background jets multiply this
increment before the error is added.

## Histories and Q

The first-micro correction starts only at the true r_minus source inlet.
Background own histories, including the original core mean, remain nonzero.
For the rate-one correction mean,

`delta_m = hb/N * integral_sc/2^s exp(-hb*(s-t))*B(t,Z,phase) dt`.

Each traversed partial source cell keeps its actual partial length and one
physical hb measure. The positive Duhamel kernel is <=1, so the complete
cellwise N-scaled density bounds give a prefix enclosure. A terminal history
cap is never substituted for an arbitrary prefix. Analogous bounds retain
all five common-unit correction histories m/h/k/e/p, including rate-zero
pressure memory. These are bounds for signed source integrals, not computed
signed integral centers. The normalized mean correction multiplies delta_m
by Pstar/epsilon.

The actual leading Qpsi is recovered from those same corrected fields:

`Qpsi=[2Z*Psi_N-(1-delta)Z*Mpsi_N-(1-Z^2)*(Mpsi_N)_Z]/(1-delta*Z^2)`.

The complete own mean_Z2 is used for Qpsi_Z. The positive-order f/u have
already vanished before r_minus. Their terminal primitive remains unchanged
through this chart; mean1 equals that primitive divided by actual rho, and
q1 uses its distinct weight1+delta. It is not reset to zero. The physical
radial value therefore retains this nonzero n=1 history term.

## Quantitative third-derivative limitation

The current conservative source bounds reveal a specific remaining gap.
Across the four original axial cells, the largest logs of the ordinary
A/N derivative caps, divided by logN, for orders0..3 are approximately

`[-1, -0.768006, -0.385481, -0.002956]`.

For normalized B/N, the corresponding leading ratios are

`[-0.958331, -0.726337, -0.343812, +0.038714]`.

Finite log(Pstar/epsilon) terms do not change these displayed digits.
The extra first-micro log(hb)/logN is about -0.0163; it does not make the
current mean_Z3 cap small. Thus U_Z3/mean_Z3 numerical accuracy and Qpsi_ZZ
are not asserted. Their actual positive log-error bounds remain exported,
without exponentiation, silent zero padding, or changing N. This shows that
the present caps are insufficient; it does not show that the true third
derivative is large. Local correlated bounds or a coherently propagated
new frequency admission are separate next work.

## Evidence

Driver: `scripts/run_original_first_micro_numeric.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_micro_numeric.json`.

The focused run exits 0 in **167.644 seconds**, producing twelve source-bound
packets: three axial points Z=-0.3,0,0.3 at native s=1/4,5/8,1, plus three
exact selected-sc coordinates. Ninety-nine independent finite-width
partial-prefix/exponential comparisons pass, alongside 81 original-N small
correction error proofs, 27 exact-flat/zero-swirl checks, eight r_minus
physical velocity coefficient seams and 19 precision60/85 native/physical
coefficient overlaps. Seven unsupported requests are rejected. The same
n=1 terminal primitive is retained at the first-micro exit.

The independent prefix fixture directly integrates an oscillating signed
source with a finite positive width/frequency, for five own histories
and three axial rows. Production bounds use traversed partial cells; the
expected integral does not use the production prefix formula. A separate
direct differentiation checks the full exponential increment errors.

At Z=0.3,s=1, Phi is approximately 0.28297998688054 with relative error
below 1.902e-44. Psi is approximately
`5.6206211351e+204453045017284567`; its absolute error is below
`4.321e+204453045017284547`, giving relative error below **7.688e-21**.
Mpsi and Qpsi value/first derivative have comparable sampled relative
errors. Psi is normalized by epsilon: physical profile recovery retains
`U=U0+epsilon*Psi`, so its enormous normalized absolute error is not an
enormous physical velocity error. These relative accuracy figures apply to
the displayed sample, not a uniform global field claim. Original F0 and
physical temporal exponentials remain factored.

The 4002142-byte report pins 1750 source inputs and 28 implementation files.
Fresh saved-report inspection exits 0 in 2.294 seconds. The existing
read-only Luna/max worker reviews source units, exact selected N, the
partial-prefix bounds, correction/proxy distinction and retained n=1 own
history; no material underbound or normalization error is found.

## Next dependencies

Continue the actual second_micro background through native phase1..2 with
its distinct coupled controls, nonzero first-chart fields and histories.
Transfer original selected-N correction enclosures and own histories using
that chart's physical measure. Establish physical total radial derivatives
with the A_phi/B_phi terms and original correlated Poisson conditioning;
Z derivatives at fixed phase do not supply them. Tighten the unresolved
third-derivative error rather than treating the existing cap as a point.

Later modified charts, exact signed phase/integral resolution where needed,
complete five-defect/selected reserved repair and repaired own P/Q/V, global
m13/m15 pressure tails/common global N, Cartesian C2/axis joins, preheat/heat,
energy/stress/flat remainder, genuine temporal n>=2 recursion, smooth sum,
oscillatory cancellation and corrected Cartesian NS remain unfinished.
The full goal remains ACTIVE / INCOMPLETE.
