# Actual long reshape physical mixed derivatives through four

The original full long reshape now supplies logR/Z velocity, pressure and
all five physical primitive derivatives of total order at most four. Its
actual R110 histories, selected T=400*Abar, full backward kernels and P0
remain unchanged. A local original post-switch power provider on R109..R110
supplies the left side of the R110 interface. This does not complete the
microscopic switches or prescribed-shear bridge.

## Corrected full kernel source

Read-only Luna/max review found that the previously printed full-kernel
definition and its independent fixture used the reversed cutoff difference.
The correct endpoint-normalized kernel is

```text
integral_0^y exp[-k*t + m*B*(sigma(y/T)-sigma((y-t)/T))]dt
```

It follows from the physical ratio of amplitudes at y-t and y, because
logu=y/10+B*(1-sigma(y/T)) plus the common axial datum. The old absolute
rate/Bell enclosures were symmetric and already enclosed both signs;
their numerical bounds did not change. The claimed exact source and
independent full-integral fixture needed correction. A new symbolic
integrating-factor identity now catches that sign error.

The original reshape, reference/restore, actual implicit patch and both
accepted mixed stages were regenerated and checked after this source
correction, so their source hashes remain current. Earlier source-bound
records are superseded by the new receipts, not silently treated as current.

## Physical radial derivatives

Write q=partial_y logu=.1-B*sigma'(y/T)/T. All original inverse-T powers
and sigma derivatives are retained. Ordinary exponential Bell derivatives
include q_y, q_yy and q_yyy; powers of q alone would be incorrect.

For radial order k>=1, physical Mtheta and Mtheta_z rows follow by
differentiating their exact RHS sqrt(2)*R^(3/2)*u and its raw-V multiple.
Their rates are 3/2+q. Physical Mz has derivative/R=V at every positive
radial order; derivatives of its normalized mean are different. The
Mztheta swirl contribution has rate 1+2q and Mp has rate 2q. The axial
V^2 term retains Pstar^-2. Pressure retains the same analytic P0.

Ur uses (D+1/2)^k Q after differentiating sqrt(R/2). Axial Taylor
coefficients become ordinary derivatives by multiplication by n!.
All exported physical normalizations are fixed at the current basepoint;
they are not differentiated as moving denominators.

Large derivative factors are combined with log(u/Pstar)^2 BEFORE numerical
positive exponential capping. The exact positive amplitude remains formal;
a cap defines only an enclosure, never the velocity source or a zero value.

## Shared actual axial source and interfaces

The new source graph retains two signed formal functions. The bridge term
is the original negative integral of chi*(phi_actual/phi_bar) times the
current-radius hydro, Pstar^2 pressure and F0^2 swirl contributions. The
first-switch term is the negative hb^2 integral with weight1-sigma(t),
the same actual/comparison quotient and radius100*exp(hb*t). Source hashes
bind their original formulas, cutoff and core. These exact functions are
distinct from their directed numerical caps.

The same two integral references occur in

```text
V100 = 4Z+j+epsilon_core*Psi(4,Z)+I_bridge
V110 = V100+I_first_switch
E    = j+epsilon_core*Psi(4,Z)+I_bridge+I_first_switch
```

Three identities at each Z Taylor order0..5 prove V110=4Z+E without choosing
independent representatives of the cap intervals. Existing actual-source
checks bind the enclosures to switch.v_cover and its unchanged Rsh copy.
The second switch has no new axial increment. The signed source integrals
are retained formally; their point values have not been numerically rebuilt.

At R110, both providers inherit exactly the same raw phi/V/five moments/P0
from switch.post(110). At Rsh, the reference inherits the identical actual
reshape parent, with the exact centered-moment change of coordinates.
Original cutoff endpoint flatness gives q=.1 and higher radial log-u
derivatives zero at both interfaces. Identical physical primitive RHSs
then give the functional mixed joins through four. Interval overlap is
only an additional consistency diagnostic. The earlier standalone
reference receipt's conditional Rsh gate is superseded by this combined
reshape/reference interface certificate.

## Evidence and reproduction

- 480 actual velocity/pressure and 600 five-primitive mixed bounds over
  eight whole-interval/endpoint/interior packets.
- 135 independent derivatives from full finite physical integrals with
  nonconstant log velocity and nonzero histories. This fixture checks
  formulas and does not admit the actual source.
- Five exponential Bell identities, six actual centered-history identities,
  eight flat endpoint log-velocity identities and eighteen shared signed
source Taylor identities.
- 480 factored positive-source cap inequalities and ten original cutoff
  endpoint checks, with source/family hashes and exact actual histories.
- Independent canonical source-tree equality checks every integral bound,
  sign, weight, quotient, current radius and physical drive scale against
  the original definitions. All R110 two-sided physical mixed grids also
  have identical exact directed endpoints; printed widths are not compared.

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage reshapemixed`
under experiments/root_st073, using accepted prerequisites. The complete
ordered pipeline has 98 modules. Producer/checker are
`lei_ren_part1_paper_compliant_long_reshape_mixed_C4.py` and its `_check.py`.

Remaining: microscopic switch and bridge mixed4, core/bridge interfaces,
inner/pre-O3 dispatcher and complete physical Cartesian assembly. Required
energy domains, admissible stress, independent flat remainder, actual
n-dependent recursion and oscillatory correction remain unfinished. Full
field/stress/global-energy/temporal gates remain false; the original
unlocalized whole-space energy remains infinite.
