> Successor: [actual original-N macro numeric values](CURRENT_ORIGINAL_MACRO_NUMERIC_2026_10_15.md) now reach R100 with true radius and own-mean decay. This second-micro document remains the accepted prerequisite; next work starts at first_switch.

# Actual selected-N second-micro velocity and local axial Z3 — 2026-10-15

## Evaluated domain and source

`OriginalSecondMicroNumeric` continues actual modified leading+n=1 physical
velocity enclosures through the second microscopic window, native s in[1,2].
It consumes the complete first-window source and the original second-window
comparison direction; the field is not restarted at s=1.

```python
from openai_ns_reconstruction.original_second_micro_numeric import OriginalSecondMicroNumeric
field = OriginalSecondMicroNumeric()
inlet = field.evaluate('.3', '1', log_tau='-2', azimuth='.25')
inside = field.evaluate('.3', '1.5', log_tau='-2', azimuth='.25')
exit = field.evaluate('.3', '2', log_tau='-2', azimuth='.25')
```

The source radius remains R=Ra*exp(original_hb*s), rho=4*exp(original_hb*s).
The terminal radius is R0=Ra*exp(2hb), the actual following macro inlet.
These native coordinates differ from the selected n=1 inner/lower fractions.
Requests outside[1,2], including first-chart selected-sc coordinates, are
rejected. No raw uncut n=1 solver domain is enlarged.

The physical packet supplies cylindrical/Cartesian leading+n=1 values,
shared streamfunction and factored coordinates. Original F0, time and
geometry logarithms remain unmaterialized. Numerical background proxies
have full actual source/correction errors; signed loop values and the
high-frequency phase are not resolved.

## Two coupled controls, one physical measure

The accepted bridge producer retains these cumulative control masses:

```
first window: hb*wc(1) + hb^2*ws(1), wc(1)=ws(1)=1/2
second window: hb^2*(s-1)
```

Each mass is combined with its own known comparison direction and quotient.
Both the angular and axial equations keep their controls, and all three
axial drive parts (hydro, original pressure, swirl) remain. The finite-width
nonlinear prefix error is propagated after the actual source factors are
cancelled or bounded in their original units.

The original physical log-radius measure is dy=hb*ds. The extra hb in the
second control produces hb^2 in the integrated field source. It does not
justify applying the physical Jacobian a second time. Actual background mean
errors are bounded cumulatively from the nonzero core inlet through the two
windows, using the true history equation.

The canonical selected repair-only integer remains

`N=2^(2^1358356628656378313)`.

No smaller finite candidate replaces it. The same original signed Section11
functions and phase frac(N*hb*(s-sc/2)) define the modified field. Full-phase
source caps enclose their numerical proxies, while the shared source N,
pressure datum, family and original width remain pinned. This selection
still does not admit missing global derivative/stress/heat requirements.

## Incoming history and derivative scope

Correction history bounds include **both** first_micro sc/2..1 and
second_micro1..s. Every traversed partial cell contributes its own length,
its complete N-scaled density cap and one hb/N factor. Positive Duhamel
kernels are bounded by1. This is a conservative enclosure of the defining
signed integral; it is not an evaluated signed integral center. Only the
original r_minus correction inlet is zero.

At s=1 the complete five incoming correction bounds match the accepted
first-chart endpoint. Their nonzero value/derivative memories are preserved
at later second-chart queries. Rate-zero pressure correction memory is also
retained. The positive-order n=1 terminal primitive is unchanged after its
f/u support; mean1 uses primitive/actual_rho and q1 retains weight1+delta.

Local second-window A/N and normalized B/N now have small ordinary derivative
caps through Z3. Their largest log-cap/logN ratios for orders0..3 are roughly

```
A/N: [-1, -0.873688, -0.597333, -0.320979]
B/N: [-0.974640, -0.832018, -0.555664, -0.279310]
```

Finite original Pstar/epsilon normalization terms do not change those digits.
The actual background already supplies genuine Z3 jets. The third exponential
Taylor error uses

`3*(X3+X1*X2+X1^3/6)`

for positive bounds Xk on the normalized Taylor coefficients of A/N, after
proving |A/N|<=1. The full background Phi jet multiplies this increment;
B_over_Pstar is converted to normalized Psi with Pstar/epsilon exactly once.
Ordinary cap rows divide by k! once.

The numerical interface therefore returns **Phi/Psi Taylor rows0..3**,
**own mean rows0..2**, and **Qpsi rows0..1**. It does not fill the absent mean_Z3
or Qpsi_ZZ with zero. The first-window incoming mean_Z3 cap remains too broad
and decays only by an almost-unit kernel in the second window. Small local
velocity derivatives do not erase that memory. Formal third-history error
logs remain exported. This is a cap limitation, not proof that the true
mean derivative is large.

These are native fixed-phase Z derivatives. They do not provide total
physical radial derivatives, where the fast A_phi/B_phi terms survive after
N cancels1/N. Physical gradients and vorticity beyond r_minus remain open.

## Executed evidence

Driver: `scripts/run_original_second_micro_numeric.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_second_micro_numeric.json`.

The focused run exits 0 in **191.175 seconds**, evaluating ten original
second-window packets, including Z=-0.3,0,0.3 at s=1,3/2,2 and an off-grid
Z=0.17,s=1.37 query. Thirty-seven independent coupled-control measure,
third-exponential and nonzero incoming-memory comparisons pass. The driver
also passes 153 selected-N correction-error proofs, 57 first/second native
and physical coefficient seam overlaps, 63 incoming/terminal-history checks,
and 21 precision60/85 native/physical coefficient comparisons. Seven
unsupported requests are rejected.

The control fixture integrates the first smooth control and the second
constant control independently at finite width. It compares the angular
source and all three axial source measures. A separate direct third
exponential derivative checks the new Taylor error formula. Finite-frequency
signed prefix integrals start at the original initial radius even when the
query is in the second chart; their first-chart incoming memories are
nonzero at s=1.

At Z=0.3,s=2, the normalized Taylor coefficient Phi[3] is approximately
-4.74699504548e-48 with error below 5.562e-85, a relative bound below1.172e-37.
Psi[3] is approximately `1.13962018904e+204453045017284568`, with relative
error below7.688e-21. Multiply both Taylor coefficients and their errors by6
for ordinary third derivatives. The huge normalized Psi exponent must retain
the original epsilon factor when recovering U; it is not a large physical
velocity error. Sampled relative bounds do not establish uniform global
relative accuracy or measured blow-up dynamics.

The saved 3491005-byte report pins 1750 source inputs and 30 implementation
files. Fresh saved-report inspection exits 0 in 2.287 seconds. The existing
read-only Luna/max worker independently reviews cumulative controls, partial
prefixes, width measures, local Z3 calculus and incoming mean/primitive
retention, and finds no material mathematical underbound.

## Next dependency

Continue at the actual R0 macro inlet with the original frozen-comparison
background and its length log(100/Ra)-2hb, retaining actual micro endpoint
fields and own histories. Transfer the same canonical-N correction sources
and first/second incoming histories. Derive total radial/phase rows with the
original correlated Poisson conditioning; tighten the incoming third-mean
bound without changing or zeroing its source.

Global modified chart coverage, exact signed phase/integral resolution where
needed, complete five-defect/selected reserved repair/own pressure recovery,
global m13/m15 tails/common global N, Cartesian C2/axis joins, heat/preheat,
energy/stress/flat remainder, true n>=2 temporal recursion, smooth summation,
oscillatory cancellation and corrected Cartesian NS remain unfinished.
Full goal ACTIVE / INCOMPLETE.
