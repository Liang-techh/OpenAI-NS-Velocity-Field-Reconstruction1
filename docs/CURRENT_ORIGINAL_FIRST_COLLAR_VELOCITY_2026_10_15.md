# Physical leading+n1 first-collar velocity — 2026-10-15

Subsequent progress: [actual native radial derivatives and core source seam](CURRENT_ORIGINAL_FIRST_COLLAR_DERIVATIVES_2026_10_15.md)
are now evaluated on the same admitted first-half collar. Complete physical
Cartesian C2 joins and higher n1 axial numerical data remain open.

## Delivered interface and domain

`OriginalFirstCollarVelocity.evaluate(Z,s,...)` combines the evaluated
actual leading collar, its own mean/Q, and the complete-error n1 transport
into physical Cartesian u,v,w, cylindrical velocity, coordinates and one
meridional streamfunction. Its original domain is `|Z|<1, 0<=s<=1/2`.
Time uses the original log(tau) or source-scaled time; azimuth is arbitrary.

```python
from openai_ns_reconstruction.original_first_collar_velocity import OriginalFirstCollarVelocity
field = OriginalFirstCollarVelocity()
packet = field.evaluate('.3', '.5', log_tau='-2', azimuth='.25')
u_terms = packet['Cartesian']['u']
v_terms = packet['Cartesian']['v']
w_terms = packet['Cartesian']['w']

late = field.evaluate('.3', '.5', log_tau=None,
                      scaled_time='1e18', azimuth='.25')
```

The output is a signed exponential sum with directed coefficient and log
enclosures. It preserves the original unmaterializable amplitude and
positive collar width, and does not provide an ordinary float vector.
This is a native chart for physical components at an explicitly represented
physical location, not a new direct Cartesian inverse. The accepted
`OriginalCoreVelocity.cartesian/at_time` API keeps its original core domain.

## Source and geometry factors

Let `rho=4exp(hb*s)`, `R=epsilon*rho`, `d=1-Z^2`,
`lambda=sqrt(tau/d)` and `L=1-delta*Z^2`. The accepted core physical
weights at rho=4 are reused algebraically; original collar fields and
geometry are then expanded separately. This does not evaluate a core
continuation at the collar point.

All radial/swirl velocity terms and the physical radius receive the
additional log `hb*s/2`. Every streamfunction term receives `hb*s`.
Axial velocity, physical z and the time weights receive no separate
collar geometry factor. Cartesian components are rebuilt after these
changes using the actual common azimuth.

In the natural n1 frame of the transport interface,

```
uz = lambda^(-1-delta)*(U0 + epsilon*Psi_s)
   + epsilon*lambda^(-1+delta)*u_s

ur = lambda^(-1)*sqrt(epsilon*rho/2)*(Qaxis + epsilon*Qpsi_s)
   + epsilon^(3/2)*lambda^(-1+2delta)*sqrt(rho/2)*q_s

utheta = F0*lambda^(-1-delta)*sqrt(2*epsilon*rho)*Phi_s
       + F0*lambda^(-1+delta)*sqrt(2*rho/epsilon)*f_s

psi = epsilon*rho*lambda^(1-delta)*(U0 + epsilon*Mpsi_s)
    + epsilon^2*rho*lambda^(1+delta)*mean_u_s
```

Leading source packets export `inlet + hb*increment_per_hb`; they receive
exactly one separate original log(hb). The n1 value packets f,u,k instead
export complete actual finite-width increments and receive no extra hb.
The same own mean `mean_u=u+k` and its Z derivative recover q with the
n1 `(1+delta)` weight. No independent fitted radial component is added.
The phase is already present in each source integral; no extra factor s
is applied to a field increment.

The physical u and v each retain nine terms and w retains five terms.
Small leading corrections, collar increments and temporal orders stay
separate. The exact source uses one streamfunction; symbolic divergence
identities do not claim a separately measured Cartesian residual bound
for the rounded coefficient centers.

## Original width and log-scale schema

The width remains the original defined positive scalar
`hb=cstar*(actual physical norm sum K in 9.16)^(-100)`. The accepted ledger
provides a log enclosure, not the numerical value of that functional norm.
No endpoint, midpoint or upper cap is selected as the physical field width.

Each term is `coefficient*exp(sum(log_scale_terms))`. Log rows have three
forms, and a consumer must recognize all three:

1. Ordinary directed summaries contain `approximation` and
   `source_enclosure_exact_mpf_tuples`. Their approximation is a numerical
   proxy and their entire enclosure must accompany it.
2. `negative_of_exact_mpf_tuple` is the exact original log(Cstar) anchor
   contribution. It must remain separate from the relative F0 log.
3. `kind=original_log_hb` defines log(hb) and supplies exact directed endpoint
   tuples directly. It deliberately has no selected `approximation`.

The geometry row `kind=original_hb_linear_log` is an ordinary bounded log
summary augmented by its exact definition `multiplier*hb*phase`. Its zero
numerical proxy is explicitly marked and is not a choice hb=0. It preserves
the positive-width source and applies the original width-log smallness
proof before constructing a numerical upper cap.

Readers must use the endpoint tuples for a tagged width row, or reject
an unsupported tag. Missing `approximation` must never default to zero.
Summing all huge logs can lose useful scale information even when their
interval sum remains valid; keep the separate terms for comparisons.

## Temporal comparison

The actual finite-width swirl fields have the exact ratio

`log(abs(utheta_n1/utheta_n0)) = logLambda + 2delta*log(lambda) + log(abs(f_s/Phi_s))`.

The original F0 and full collar radius cancel. The same shared radius
geometry tag is removed from both returned physical terms in the
independent comparison. Leading width corrections and complete n1
increments remain in the coefficient enclosures. This compares the first
two temporal orders, and does not admit higher orders or their sum.

## Executed evidence

The focused run exits 0 in **162.954 seconds**. It computes ten native
physical packets, including the offgrid point Z=0.17,s=0.37, and three
source-scaled-time packets. Seven independent symbolic geometry,
streamfunction and rotation identities pass, together with 190 emitted-term
geometry/width-factor checks and 13 independent order-ratio comparisons.

A separately constructed accepted core evaluator agrees at the phase-zero
seam in 60 coefficient enclosure comparisons, including physical Cartesian
components, streamfunction and coordinate coefficients. All added phase-zero
field increments and geometry logs are exactly zero. Six unsupported
requests are rejected. This checks the value seam; it does not certify
first/second spatial derivative joins.

The 5842867-byte report pins 1478 source inputs and 13 implementation files.
A fresh-process saved-report inspection exits 0 without rebuilding source
owners. The existing read-only Luna/max worker confirms the geometry,
normalization, width/temporal split and common-streamfunction recovery.

At Z=0.3,s=0.5, ordinary log(tau)=-2 classifies the first n1 swirl as larger
than the leading swirl. Original source-scaled time S=1e18 classifies it
as smaller. The comparisons retain original delta*time products and
cancel F0 and the actual collar radius analytically. They do not establish
admission of the full temporal series or blow-up diagnostics.

Driver: `scripts/run_original_first_collar_velocity.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_collar_velocity.json`.

## Remaining work

Task 14 has physical components on the core and the admitted first-half
collar; it remains IN PROGRESS. Task 15 also remains IN PROGRESS. The full
goal remains ACTIVE / INCOMPLETE.

1. Establish quantitative first/second spatial derivative joins, beyond
   the checked value seam and available native axial source derivatives.
2. Admit and evaluate remaining collar domains; apply source-owned cutoff
   and independent own five-moment repair to numerical fields.
3. Supply a safe generic log-term reader before downstream materialization
   or plotting. A new direct Cartesian collar inverse needs a separately
   resolved/factored phase; the current core inverse must not be reused by
   silently rounding rho to 4.
4. Complete global matching, five-moment numerical closure, exact heat
   exterior, finite-energy diagnostics, admissible stress and flat remainder.
5. Construct genuine temporal n>=2 coefficients, their own repair and
   smooth summation, oscillatory corrections and independent corrected
   Cartesian NS validation.

The interface does not establish a global physical velocity field, useful
combined physical relative error, full recursion or a blow-up theorem.
