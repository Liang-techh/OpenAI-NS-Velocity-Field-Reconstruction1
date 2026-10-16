# Original finite long reshape: logarithms and six histories — 2026-10-15

The physical successor is now [CURRENT_ORIGINAL_LONG_VELOCITY_2026_10_15.md](CURRENT_ORIGINAL_LONG_VELOCITY_2026_10_15.md).
It consumes this numerical reader and assembles physical Cartesian value terms,
with original time/geometry factors and n1 primitive cancellation. Direct physical
inversion, derivatives and later charts remain open.

`OriginalLongNumeric` continues the actual unmodified R110 inlet through
the original finite long reshape. It evaluates the selected-N logarithmic
angular profile, normalized axial velocity, own mean and recovered Q,
and separately transports all six unmodified leading-background histories.
The selected repair-only integer and all incoming correction memories remain.

```python
from openai_ns_reconstruction.original_long_numeric import OriginalLongNumeric
field = OriginalLongNumeric()
middle = field.evaluate('.3', '.5')       # q in [0,1], y=T*q
near_inlet = field.at_log_radius('.3', '1')
terminal = field.evaluate('.3', '1')
```

## Exact source, logarithmic geometry

T is decoded from the accepted C6 geometry certificate's equal exact mpf
endpoint tuples. It is the original finite singleton T=400*A, not a bound,
shortened numerical window or infinite-time replacement. The actual geometry is

```
y = T*q
R = 110*exp(y)
rho = 110*Lambda*exp(y)
dy = T*dq
B_log = -Lambda*G + log(Phi110) + log(220)/2 + log(1+Z^2)
logPhi = log(Phi110) - B_log*sigma(q) - 2*y/5
logUtheta = B_log*(1-sigma(q)) - log(1+Z^2) + y/10 - logCstar
```

No R, Phi, F0 or exponential of T is materialized. The physical-similarity
swirl logarithm includes the original logCstar once; fixed-time physical
Cartesian vector assembly is a separate next step.

G0 uses the accepted six-pole anchored primitive with its full directed
enclosure. Higher G Taylor coefficients use the independently checked exact
rational derivative L*H/(H^2+(j/500)^2), with G[k]=gradient[k-1]/k. This works
at arbitrary interior decimal Z and does not interpolate the old 0/.5 frames.
The positive real H0-squared denominator is preserved before inversion.

At q=1 the inlet logPhi and G contributions to logUtheta cancel analytically,
giving T/10-logCstar-log(1+Z^2). At very small q the B_log*sigma product is
bounded in logarithms before applying an error cap. A clipped sigma value
is never selected as the source function.

## Finite kernels and nonzero histories

The four angular histories are stored in CURRENT Phi units, with Psi held
at the actual unmodified R110 value. Define

```
D(k,m) = exp(-k*y + m*B_log*sigma(q))
K(k,m) = integral_0^y exp(-k*t + m*B_log*(sigma(q)-sigma(q-t/T))) dt
H/Phi = (H110/Phi110)*D(1.6,1) + 2*K(1.6,1)
Kpsi/Phi = (Kpsi110/Phi110)*D(1.6,1) + 2*Psi*K(1.6,1)
B/Phi^2 = (B110/Phi110^2)*D(1.2,2) + K(1.2,2)
C/Phi^2 = (C110/Phi110^2)*D(.2,2) + K(.2,2)
Mpsi = Psi + (Mpsi110-Psi)*exp(-y)
Apsi = Psi^2 + (Apsi110-Psi^2)*exp(-y)
```

Full finite kernels are enclosed using local positive rates
k +/- 8*m*abs(B_log[0])/T and derivative Bell bounds. Numerical kernel centers
are constant-rate proxies with complete directed errors. A center equal to
1/k at huge y is not a claim that the defining finite kernel is an infinite
integral. Inherited coefficient/derivative factors enter the logarithmic decay
bound before an exponential error cap. All nonzero R110 inlets survive.

These are leading-background histories before finite-N modification. Raw units
remain M=U0+epsilon*Mpsi, K=U0*H+epsilon*Kpsi, and
A=U0^2+2*U0*epsilon*Mpsi+epsilon^2*Apsi. C remains undressed; later pressure
recovery applies F0-squared ratios once. This does not export pressure.

## Original finite-N and positive-order n1 data

The exact integer remains N=2^(2^1358356628656378313). Local A/N is added to
logPhi/logUtheta; the normalized axial B_loop/N correction and own mean
correction retain the complete preceding prefix plus the current T*dq measure
once. All five signed source histories remain present as logarithmic bounds.
Their centers and the original fast phase are unresolved, rather than chosen.

Positive-order n1 velocity coefficients vanish here, but the actual terminal
primitive remains nonzero and is retained. Its mean is tagged as the actual
primitive divided by (110*Lambda), times exp(-y); recovered q retains the same
source. This factored history has not yet been merged into a physical vector.

## Focused run and persisted report

The 172.601 s focused run exits0 with11 original packets,81 independent finite-kernel/memory/cutoff/G checks,165 selected-N/flat-tail small-error proofs,93 R110 field/six-history seam coefficients,86 history/integer/primitive checks,12 exact terminal cancellations,145 finite-log/history checks,11 exact-T measure relations and44 precision comparisons. Five invalid requests are rejected. Saved-report inspection exits0; input/implementation pins total1754/38. Existing read-only GPT-5.6 Luna/max review confirms the full perturbed-rate enclosure, finite source measures and retained memories.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_numeric.json`,
1,476,990 bytes, SHA256 `7190258a960c80fe738d7d1c612519355eea4411c4896cd70a78028a17a9d997`.
The proof intervals serialize as exact endpoint tuples. The prior serialization
failure did not invalidate the completed numerical checks; the repaired run and
saved-report inspector both finish successfully.

At Z=.3 the finite-y history proxies evolve as follows, with complete directed
errors in the report. They are not substituted equilibrium source values.

| y=log(R/110) | H/Phi proxy | C/Phi^2 proxy |
| --- | --- | --- |
| 0 | 1.03535947970 | 1.07552601709 |
| 1 | 1.20666482633 | 1.78691246054 |
| 10 | 1.24999997585 | 4.46888020197 |

## Limits and next work

logPhi/Psi have Taylor Z0..3, own mean Z0..2, Qpsi Z0..1, and the six
unmodified history ratios Z0..2. The packet provides source enclosures and
numeric proxies with complete errors, rather than exact signed point values.

Tasks14/15 remain IN PROGRESS. Next assemble the long-chart fixed-time physical
vector with its retained n1 decay and original time/geometry factors. Continue
the reference/restoration charts using this actual endpoint. Incoming meanZ3,
total radial/fast-phase derivatives, merged modified six moments, own pressure,
signed five-defect/selected repair, global common N, C2/heat/energy/stress and
flat remainder, genuine temporal n>=2/smooth sum, oscillatory cancellation and
independent corrected NS remain open. Full long-term goal ACTIVE / INCOMPLETE.
