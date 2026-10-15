# Original fully coupled n1 radial coefficients — 2026-10-15

## Delivered and scope

`OriginalN1CoreNumeric` computes finite radial Taylor prefixes of the
**temporal n=1** original regular core. It obtains the leading coefficients
from the accepted original coupled production recurrence, retains all raw
B0/B1 entries, axial derivatives, pressure and nonzero feedback, and returns
numerical approximations with directed original-source coefficient bounds.
It uses the current leading/n1 source family and no interval midpoint.

The new implementation is independent of the older three-row `axis_solve`
shortcut. Radial degree is an approximation index within temporal n=1;
it is not temporal order and does not implement n>=2 recursion.

```python
from openai_ns_reconstruction.original_n1_core_numeric import OriginalN1CoreNumeric
solver = OriginalN1CoreNumeric(precision=70, quadrature_nodes=192)
packet = solver.coefficients('0.3', degree=16)
prefix = solver.evaluate_prefix(packet, '2')
prefix['normalized_prefix']  # f,u,k,p,mean_u,q, finite-prefix errors only
```

The finite coefficient bounds do **not** include the infinite radial
remainder. Consequently this is not an original full-n1 point oracle,
a certified nonzero-radius n1 solution, n0+n1 physical assembly, a global
field or a full NS solution. The full goal remains ACTIVE / INCOMPLETE.

## Exact normalization and recurrence

Let rho=Lambda R, epsilon=1/Lambda, d=1-Z^2, L=1-delta Z^2.
The local normalization is

- f=epsilon F1/F0(Z);
- u=Uz1/epsilon, k=K1/epsilon;
- p=P1/(epsilon Pstar^2).

The F0 derivative F0_Z/F0=-g/epsilon enters once, including in the
second weighted axial forcing. The original gauge is not replaced by a
constant amplitude. The leading recovery uses **1-delta** in its mean
coefficient; the n1 recovery uses **1+delta**:

`q=Q1/epsilon=[(1-delta)Zu-(1+delta)Zk-d(u_Z+k_Z)]/L`.

For the full normalized right sides Rf,Ru,Rp obtained from original raw
rows 5,6,4, respectively,

`f_(n+1)=[rho^n]Rf/[4(n+1)(n+2)]`,

`u_(n+1)=[rho^n]Ru/[4(n+1)^2]`,

`p_(n+1)=[rho^n]Rp/(n+1)`,

`k_n=-n*u_n/(n+1)`.

All axis constants vanish. The recurrence keeps two extra actual leading
Z orders and one extra radial order before constructing Q, Q_Z and the
pressure head `(n+1)(n+2)Q_(n+1)`. Missing derivatives are rejected rather
than padded with zero. Exact radial constants are promoted explicitly.

The numerical proxy omits the unmaterializable F0-squared contribution;
the separate directed branch retains its actual nonzero Cauchy enclosure,
all original pressure atoms in aggregates and their source uncertainties.
The proxy is never selected as a replacement original source. Outputs
retain enormous scales through arbitrary exponents and normalized terms.

## Important older axis limitation

The older `current_original_n1_regular_system.axis_solve` omits the direct
B0[6,1] contribution to the axial recurrence. In raw units this coefficient
is `-8 epsilon^2 Z rho F/L`; in the new normalized axial RHS it is
`-8 Z rho F0^2 A f/L`. It first affects the rho^3 axial coefficient. It is
not absorbed into the pressure term. Comparisons against that older row
would therefore certify the wrong equation.

The new implementation includes this term and tests it with a nonzero,
Z-dependent F0 fixture. Older archived three-row axis data are not used as
numerical truth by this solver. Historical archives are preserved; the
scope of any downstream axis claim must be checked before reuse.

The consumer trace isolates this defect to the legacy axis packet and its
checker. Its Uz1_rho3 and K1_rho3 are stale; its stored P1_rho3 and Q through
rho2 do not acquire this omitted term at that order. Current Cartesian
axis export uses equation-derived rho1/rho2 limits and its transverse
curvature specifically uses Uz1_rho1. It does not consume the legacy axis
archive. The positive-radius high source dispatcher also differentiates
the actual six-row equation and rejects the axis. Those source paths and
the first axis curvature are unaffected by this legacy rho3 defect.

## Executed evidence

The actual-source run exits 0 in 89.469 s. It computes degree16 prefixes
at Z=-0.3,0,0.3, rho=0,1,2,4 with 50/70 decimal digits and independently
changed 96/192 pressure quadrature panels. The focused checks cover:

- 30 exact normalized raw-matrix state-coefficient identities;
- 36 exact mixed coefficient comparisons against direct raw equations;
- nonconstant F0 and g, their gauge derivative, and nonzero direct F0-squared axial feedback;
- 72 cross-precision finite-prefix/error comparisons;
- 2040 directed mixed coefficient slots, with 1824 nonzero point proxies;
- 54 degree8/12/16 refinement diagnostics.

At Z=0.3,rho=2, the degree16 normalized finite prefixes include
f=-0.18941311938863450413, u=6.33600000000000000000,
p=-1.97950082489533307016 and q=-4.13760000000000000000.
Errors for the full serialized finite-prefix decimals are approximately
3.04e-71, 2.20e-70, 5.48e-27 and 1.62e-70, respectively. These are
**finite-prefix errors**, not full n1 point errors or NS residuals.

For f at that point, the degree8->12 change is 3.63e-12 and the
12->16 change is 1.28e-20. Some u/p scalar differences round to zero
because later coefficients are extremely suppressed. Their nonzero
coefficients remain in the saved report; zero displayed differences do
not establish an exactly vanishing correction or an infinite-series tail.
A separate saved-report inspection exits 0 in 2.135 s without owners.

## Reproduction and next work

Run `python scripts/run_original_n1_core_numeric.py --write-report` for
exact raw-equation tests and original-source computations. Run
`python scripts/run_original_n1_core_numeric.py --inspect-report` to verify
the saved source/implementation pins without rebuilding historical owners.
The JSON report saves all mixed coefficient enclosures, 12 sampled finite
prefixes at two precisions, and radial degree8/12/16 refinement diagnostics.
Refinement differences are observations, not a radial remainder proof.

Next obtain an original full regular n1 evaluation with a valid radial
remainder (or actual full Volterra quadrature/iteration plus its own word
tail). Do not apply an existing Volterra/Picard word tail to a radial Taylor
polynomial. Then propagate actual nonzero histories, apply the source-owned
cutoff/moment repair and assemble n0+n1 through one streamfunction. Higher
inner interfaces, inverse physical chart selection, global physical/heat/
stress layers, common global N, n>=2 recursion and corrected NS stay open.


The tail route was checked against existing estimates. The radius-5 analytic
source bounds in `current_original_n1_inner_analytic_domain.py` are leading/
source estimates, not a complex-rho sup norm of the solved n1 field. The
accepted weighted feedback estimates use Z Cauchy loss and Volterra word
length, with a uniform Fmax basis, not a local F0(Z) radial circle norm.
A Cauchy radial tail would need a new joint-holomorphy and circle-bound
lemma for the full coupled solution. The next computational route should
therefore evaluate Wbase=(I+J)Gg and actual H iterations with all B1 d_Z
terms, propagate original-source/quadrature errors and then apply the
accepted word remainder in its own original component frame. The current
finite radial prefixes provide independent coefficient diagnostics.
