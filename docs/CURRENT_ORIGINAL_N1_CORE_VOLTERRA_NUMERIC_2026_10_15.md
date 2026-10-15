# Original full n1 core numerical values and errors — 2026-10-15

## Delivered

`OriginalN1CoreVolterraNumeric` evaluates the original regular temporal-n1
core in its natural component frame, including a directed absolute error
relative to the **full infinite coupled solution**. This advances beyond
the finite unknown-coefficient prefixes in the preceding checkpoint.

The evaluator analytically integrates finite forcing polynomials with the
same accepted Kchi principal inverse, executes a nominal C0+B1*dZ feedback
iterate, and bounds all omitted feedback using the accepted original
index-zero word remainder. Separate budgets cover the admitted leading
source tail and finite potential-kernel tail. No radial Taylor remainder
of the unknown n1 solution is assumed. No interval midpoint is selected.

```python
from openai_ns_reconstruction.original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric
solver = OriginalN1CoreVolterraNumeric(precision=70, quadrature_nodes=192)
result = solver.evaluate('0.3', '2', degree=18, kernel_degree=24)
f = result['normalized_original_n1']['f'][0]  # value, original absolute error
fZ = result['normalized_original_n1']['f'][1] # first ordinary Z derivative
q = result['normalized_original_n1']['q']
```

Domain: native 0<=rho<=4, -1<Z<1, with decimal-string inputs and the same
original source family. Axis values and ordinary Z derivatives of all six
state components, mean and q are exactly zero. The API is not restricted
to stored sample locations. It restores archives directly without replaying
historical owner constructors.

## Normalization and physical scope

With epsilon=Lambda^-1 and x=sqrt(rho), the six numerical components are

`f=epsilon F1/F0(Z), u=Uz1/epsilon, k=K1/epsilon,`

`p=P1/(epsilon Pstar^2), fx=partial_x f, ux=partial_x u`.

The same natural frame is used by the accepted regular-series bound:

`Wphysical=diag(F0*Lambda,epsilon,epsilon,epsilon*Pstar^2,F0*Lambda,epsilon) Wnatural`.

This is a full **regular n1 core coefficient** numerical approximation.
It is not an exact source point oracle, an ordinary physical velocity
vector, a completed n0+n1 physical assembly, a whole matched inner/outer
field, independent n>=2 temporal recursion or a full NS solution. In
particular, the values below must not be interpreted as physical speeds.
Absolute F0 remains in the previously installed logarithmic representation.
The full goal stays ACTIVE / INCOMPLETE.

## Why the full error applies

The exact accepted decomposition is

`W=sum_j>=0 V^j W0, W0=Gprincipal*g, V=Gprincipal*(C0+B1*dZ)`.

Gprincipal contains the infinite swirl kernel

`Kchi=sum_j>=0 (-2chi)^j T^j K0`,

`K0=s*(1-s^2/x^2)/2`, `chi=H0^2/(H0^2+sigma^2)`.

The other kernels are exact: axial `s*log(x/s)`, kinematic
`-s*(1-s^2/x^2)/2`, pressure `1`, and axial derivative `s/x`.
In rho, the zero-potential monomial inverse divides source row n by
4(n+1)(n+2) for f, 4(n+1)^2 for u and n+1 for p, then raises its degree.
Repeated analytic integration produces the signed Kchi polynomial response.
This is integration of a finite forcing approximation, not evaluation of a
finite radial Taylor polynomial of the unknown as a presumed full solution.

The principal cancellation in the computed feedback is performed through
the exact identity H-H0=d(U-U0). It does not subtract rounded near-equal
O(1) quantities. Every remaining raw B0/B1 coefficient, pressure term,
ordinary Z derivative and original amplitude gauge is retained once.

The directed forcing approximation includes the actual nonzero implicit
swirl and all original pressure atoms through their accepted enclosures.
The point proxy alone omits unmaterializable F0-squared source values;
that is not a selection of a different original source.

Three distinct errors are combined:

1. The accepted original leading fixed-point norm supplies the nonlinear
   source-tail bound. The known Bessel model has its own factorial radial
   tail; its Z derivatives follow Cauchy after checking |chi|<=1 on the
   same original complex Z tube. An exact difference-jet product rule
   propagates those omitted leading tails through the full forcing map.
2. The finite Kchi series has its own kernel and derivative remainder,
   multiplied by a bound for the entire finite source integral.
3. The saved regular all-word remainder T0 bounds every true feedback
   word after W0. It is multiplied by rho for f,u,k,p and sqrt(rho) for
   fx,ux; its accepted Cauchy derivative bound covers ordinary Z1.

The rigorous enclosure stays centered on the directed principal polynomial
integral. Distance from the returned nominal W0+V W0 value to that interval
charges the computed nominal feedback correction conservatively. No claim
that an index-one word remainder certifies that nominal iterate is made.
Fmax in the accepted majorants only bounds F0-squared coefficient norms;
it does not replace local F0(Z) in the exact natural-frame conjugation.

## Executed results

`python scripts/run_original_n1_core_volterra_numeric.py --write-report`
exits 0 in 34.672 s. The run changes precision 50->70, pressure quadrature
96->192, leading source degree12->18 and Kchi degree16->24. Twelve original
core points per level cover Z=-0.3,0,0.3 and rho=0,1,2,4.

- 12 exact source-difference/value/derivative comparisons against the
  unnormalized physical weighted operators pass, with nonconstant F0/g.
- 24 exact principal-equation and kinematic comparisons pass.
- 168 cross-precision comparisons lie within the combined full-n1 errors.
- 90 computed nominal feedback components are nonzero and stored separately,
  retaining tiny terms that would disappear in an ordinary scalar sum.
- Five unsupported input/domain requests are rejected.

At Z=0.3,rho=2, the refined value/error pairs for full serialized decimals are:

| Natural component | Approximation, shortened | Full original absolute error |
| --- | ---: | ---: |
| f | -0.18941311938863450413 | 2.48e-36 |
| u | 6.33600000000000000000 | 2.20e-70 |
| p | -1.97950082489533307016 | 5.48e-27 |
| q | -4.13760000000000000000 | 1.62e-70 |

Across the sampled nonzero values, maximum relative errors are 2.75e-29
for f, 1.53e-28 for fx and 2.77e-27 for p. These are original normalized
coefficient errors, not NS residuals. Derivative precision is separate:
for example f_Z at this point has an error of approximately 5.66e-10,
so value accuracy must not be reused as a derivative error.

The 764951-byte JSON report pins 1476 original inputs, the current module
and its supporting numerical implementations. It saves all approximation,
full error, source-enclosure tuples, nominal feedback increments and tail
ledgers. A fresh `--inspect-report` process exits 0 in 2.077 s without
historical owners. A read-only Luna/max review found no material error in
the natural-frame conversion, leading-source difference map or tail factors.

## Next work

Task 15 is still in progress: extend the evaluated core into the actual
first/lower collar with its nonzero histories, cutoff and moment repair,
and combine the leading core and n1 through the same streamfunction and
original physical/time weights. Keep logarithmic amplitudes and tiny
corrections factored; do not silently discard a scale in their sum.

Tasks 13/14 still require higher inner interfaces, inverse physical chart
selection and global physical assembly. Common global N, quantitative
five-moment/energy/support checks, compatible preheat/exact heat evaluation,
stress/flat remainder, genuine n>=2 recursion and oscillatory corrected
NS remain open. The preceding finite-radial-prefix API retains its original
finite-only error semantics; this new full-error interface uses a different,
explicitly admitted Volterra remainder route.
