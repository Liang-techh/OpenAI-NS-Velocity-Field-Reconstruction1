# Correlated micro/macro m15 source and signed micro I1 enclosures

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The new
implementation combines the accepted first-plus-second micro pressure source
with one macro local source in a single graph namespace. It also computes
signed interval enclosures for all six full finite-N micro I1 integrals and
their ordinary Z1 derivatives on Z in [37/100,93/250]. These are interval
enclosures of source functions, not selected cap values or point samples.

## One graph and one copy of each source

The common Q4 base has 183700 nodes. The accepted 624-node correlated
suffix is retained, and the independent 102-node RQy suffix is relocated
after it. Relocation follows operation-specific graph references: rational
numerators/denominators are literals, integral binder names and metadata
integers remain unchanged. The complete serialized graph has
184687 nodes; all output handles
are included in the graph reconstruction.

Only the macro local six-I1 source is added to the whole micro pair. The
second micro contribution is not added again through the RQy seam packet.
The macro measure is R^2*length*dm exactly once on m in [0,1]. Q4 U, Q and
Q_y are already physical; no additional Pstar factor is applied.

The signed pressure reduction retains

```
B = endpoint_difference[R^2*q]/L + endpoint_difference[R^2*q^2]/2
    - 2*Z^2*endpoint_difference[R^2*U^2]/L^2
D = -I[q]/L + h*Z*I[q_Z]/(2*L) - I[q^2]/2 + d*I[U*q_Z]/L
    + 2*h*Z^2*I[U^2]/L^2 + 2*Z*d*I[U*U_Z]/L^2
h=1-delta, d=1-Z^2, L=1-delta*Z^2.
```

The micro2/macro0 RQy contraction consumes the accepted per-order endpoint
identities. The separate reduced-B contraction uses the common physical
U/Q source and radius identities, with the actual endpoint graph ports.
The first lower and macro outer boundaries remain. The pressure increment
and the shifted tail have separate definitions: the increment contains the
first-lower RQy term, whereas a tail added to the accepted shifted H anchor
omits that term. The original first-lower M15/H source is retained once.

## Signed micro integral computation

The actual core axis source is U0=4*Z+j. The physical core field and its
normalized mean are U0+epsilon*Psi and U0+epsilon*mean_Psi. The saved current
core/inlet bounds control these corrections. Native micro U uses the saved
actual physical provider, distinct from comparison profiles. Native M is
recovered by the nonreset equation M_y=U-M with y=h_b*s. Its positive Duhamel
kernel bounds the deviation from U0 while preserving the incoming history.

Ordinary M0..M2 enclosures and the saved physical M3 norm radius supply the
linear Q recovery through Q_ZZ. The M3 radius is used as an error enclosure;
its center is not asserted to be the actual derivative. All products use
ordinary derivatives. No missing derivative is padded by zero.

Write h_b=exp(loghb) and W=Ra^2*h_b, where R=Ra*exp(h_b*s) and
Ra=4*exp(-(4*logPstar+1000)). The saved scale-log field named `hb` contains
loghb. W is one exact positive physical source scale,
retained by its logarithmic decomposition. The normalized integrals use
exp(2*h_b*s)*ds on [sc/2,1] and [1,2]. The exact source sc is retained, and
the source log-hb bound controls the geometric remainder.

The accepted joint finite-N defect caps are divided by this same W using
directed logarithmic arithmetic. Each ordinary Z0/Z1 error is below
2^-1024 in normalized units. That symmetric error radius is added once to
the joint background enclosure. It is never used as a signed center.

The table shows outward-rounded enclosures for I/W; exact binary interval
tuples in the report are authoritative. Ordinary Z1 commutes with W because
the physical scale and native endpoints are Z-independent.

| full finite-N micro integral / W | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-5.821472, -5.773983] | [11.808000, 11.936001] |
| q_Z | [11.808000, 11.936001] | [31.999999, 32.000001] |
| q_squared | [16.669445, 16.944769] | [-69.485090, -68.179203] |
| U_q_Z | [17.475840, 17.760769] | [94.592000, 95.360001] |
| U_squared | [4.380800, 4.428289] | [23.680000, 23.808001] |
| U_U_Z | [11.840000, 11.904001] | [31.999999, 32.000001] |

Applying the actual signed ordinary C1 pressure coefficients to these
integrals gives D/W in [23.292141, 23.905884] and D_Z/W in
[149.485006, 151.984964]. The half-volume contribution is
(D/2)/W in [11.646070, 11.952942]. These are signed volume enclosures;
the missing outer pressure terms are not absorbed into them.

## Checks, limits and next work

New production exits 0 in 7.793 seconds. Report
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_macro_correlated_m15.json` is 18,788,901 bytes, SHA256 `8354af18754496450fb41ff0335e2c82f503c0d64852496da11ee1546fb23577`, with
2007 input and 380
implementation pins. Independent bounded Luna/max mathematical review,
saved-only inspection, graph reconstruction and scoped source-byte checks
pass. Separate runtime worker model metadata is unavailable. The new
production executes no accepted ancestor. During independent review, one
reviewer mistakenly executed the accepted
`original_n1_micro_correlated_m15.exact_algebra()` six-identity check once.
That result is excluded from acceptance; old data were not modified.
Accepted producers, owners, norms, recurrences and quadrature were not
replayed. The report records this execution exception explicitly.

The old 212 forward terms, six signed inner/lower contributions and all
existing error packets are unchanged. The macro I1 integrals are exact
source definitions at this stage; their signed numerical enclosures remain
open. The outer macro RQy is retained. These local micro enclosures do not
establish the complete pressure increment, terminal five-moment closure or
an independent NS residual estimate. P1 and m15 retain their separate data
and source definitions.

Implementation: [combined source](../src/openai_ns_reconstruction/original_micro_macro_correlated_m15.py),
[signed micro I1](../src/openai_ns_reconstruction/original_micro_background_I1_numeric.py),
[runner](../scripts/run_original_micro_macro_correlated_m15.py).

```
python scripts/run_original_micro_macro_correlated_m15.py --inspect-report
```

NEXT obtain signed macro six-I1 enclosures with the same native histories
and finite-N defect control. Carry its outer RQy to the next source before
separate boxing, then combine with the reserved-pulse terminal functional
comparison. Exact heat/energy, admissible stress/flat remainder, genuine
time-order recursion/summation, realizable stress correction and independent
Cartesian NS validation remain incomplete. Goal ACTIVE.
