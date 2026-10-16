# Actual original-N second/post-power values and six histories to R110 — 2026-10-15

The finite-long successor is now [CURRENT_ORIGINAL_LONG_NUMERIC_2026_10_15.md](CURRENT_ORIGINAL_LONG_NUMERIC_2026_10_15.md).
It consumes this actual unmodified inlet and all six histories, preserving exact T,
canonical N and complete correction memories in a logarithmic reader. The long-chart
physical Cartesian vector is still a separate open task.

`OriginalSecondSwitchNumeric` continues the accepted actual first-switch
endpoint through second_switch and post_power to the exact R110 endpoint.
One source owner supplies both charts; chart-sensitive numerical caches
are separate. Neither the incoming field nor any history is restarted.

```python
from openai_ns_reconstruction.original_r110_numeric import OriginalSecondSwitchNumeric
field = OriginalSecondSwitchNumeric()
second = field.evaluate('.3', '.5', log_tau='-2', azimuth='.25')
post = field.post_power('.3', '.5', log_tau='-2', azimuth='.25')
exit = field.post_power('.3', '1', log_tau='-2', azimuth='.25')
```

`OriginalPostPowerNumeric` is also independently constructible. The helper
`field.post_power_field()` shares source/inlet caches and resets the local
numeric cache; equal native strings on different charts cannot collide.

## Two original charts, distinct measures

Second switch uses native s in[0,1]:

```
R = 100*exp(hbS*(1+s))
rho = 4*exp(Y+hbS*(1+s))
(logPhi)_s = -hbS^2*(1-sigma(s))*R*D_over_R/2 - (2/5)*hbS*sigma(s)
U_s = 0
dy = hbS*ds
```

The comparison shear term uses its complement weight and two width
factors. The explicit power term uses sigma and one width factor. Their
prefix masses remain separate, including the exact half-masses at s=1.
The comparison direction continues at the actual radius, with its original
R0 modes converted once to R100 modes. The actual unmodified Psi is held
at its first-switch value; the local finite-N B/N correction still varies.

Post-power uses native q in[0,1]:

```
Lp = log(11/10)-2hbS
y = Lp*q
R = 100*exp(2hbS+y)
rho = 4*exp(Y+log(11/10)*q+2hbS*(1-q))
Phi = Phi_R2*exp(-2y/5)
U = U_R2
dy = Lp*dq
R(1) = 110 exactly
```

The physical reader keeps the complete log radius, including the ordinary
log(11/10)*q term. Only the unmaterialized original hb contribution is an
error-bound component. The point radius omits this tiny component and has
its full directed error. At R110 the width component vanishes exactly.
Bridge hbB and switch hbS remain distinct source nodes under the accepted
exact-reboxing certificate and identical K1 ledger.

The canonical selected integer remains
`N=2^(2^1358356628656378313)`. Its phase is the original log(R/r_minus):

```
second: frac(N*(Y+hbS*(1+s)-hbB*sc/2))
post:   frac(N*(Y+log(11/10)*q+2hbS*(1-q)-hbB*sc/2))
```

The origin shift never enters the physical radius. Local A/B error caps
and complete cumulative five-history density caps use their actual chart
and one native physical measure. Macro length is bounded above by Y and
post length by log(11/10) only inside conservative error bounds. No extra
hb is applied to the post-power density. Signed loop points and signed
correction-integral centers remain unresolved.

## Six actual leading-background histories

The packet now also supplies `unmodified_own_six` for H, Mpsi, Kpsi, Apsi,
B and C, with Taylor Z0..2 and directed absolute errors. These are
propagated point proxies with complete source errors; their centers are
not evaluations of exact signed histories. They describe the leading
background BEFORE finite-N modification. They are separate from the five
selected-N correction histories and from positive-order n1 own history.

The normalization is

```
M = U0+epsilon*Mpsi
K = U0*H+epsilon*Kpsi
A = U0^2+2*U0*epsilon*Mpsi+epsilon^2*Apsi
```

H, B and C are unchanged. C is undressed: pressure reconstruction must
apply the original F0-squared ratios once, without absorbing Pstar/P0.

The respective rates are2,1,2,1,2,1 and equilibrium targets are
Phi, Psi, Phi*Psi, Psi^2, Phi^2/2, Phi^2. R100 transport starts with the
nonzero actual core atoms and exact total bridge length Y. It uses the
positive uniform angular and axial source bounds over EVERY micro and
macro prefix. A signed endpoint integral is not a uniform source bound.

Both short switches retain their actual incoming moments and bound the
change by rate*hbS*s times the difference from the true initial target
plus its complete source variation. No comparison moments are substituted.
Post-power uses exact finite positive-kernel formulas:

```
H/Kpsi source weight: 5/4*(exp(-2y/5)-exp(-2y))
Mpsi/Apsi:           1-exp(-y)
B:                   5/6*(exp(-4y/5)-exp(-2y))
C:                   5*(exp(-4y/5)-exp(-y))
```

Incoming histories keep their own rate-dependent decay. The mean agrees
with the physical velocity recovery; n1 f/u remain zero and its terminal
primitive/current_rho and one-plus-delta q1 persist.

Phi/Psi retain Taylor Z0..3, modified mean Z0..2, and Qpsi Z0..1. The
inherited first-micro mean Z3 error cap remains unresolved. Six unmodified
Z2 histories do not close the five corrected terminal identities or
establish physical gradients/vorticity beyond r_minus.

## Focused run

Run `python scripts/run_original_r110_numeric.py --write-report` followed by
`python scripts/run_original_r110_numeric.py --inspect-report`. The report is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_r110_numeric.json`.

The driver independently integrates the full second angular control and
all six finite post-power ODEs, including genuine axial derivatives. It
checks native measures, physical geometry, both source joins, all incoming
correction/terminal memories, six-history joins, own-mean agreement,
R110=epsilon*rho, nonconstant post-power swirl and cross-precision overlaps.

The accepted run exits0 in225.136s with20 original packets,54 independent
control/six-ODE/geometry/measure checks,306 selected-N small-error proofs,
126 first/second/post native and physical coefficient seams,54 six-history
join coefficients,138 correction/terminal-memory checks,54 own-mean
overlaps,6 R110/nonconstant-swirl checks and39 precision60/85 native,
physical and six-history comparisons. Twelve unsupported requests are
rejected. Saved-report inspection exits0 in2.422s. The report pins1753
source inputs and36 implementation files.

At Z=.3 the numerical Phi proxy changes from0.282979986880542503 at R2 to
0.272394696752937809 at R110, a ratio0.962593502656167763. The H, B and C
history proxies have respective endpoint ratios about0.9966303081,
0.9933462316 and0.9965676204. Their different evolution follows their own
finite positive-kernel ODEs; the histories are not reset to local field
products. Every reported number retains its directed absolute error.

## Next action and full scope

Tasks14/15 remain IN PROGRESS; full goal ACTIVE / INCOMPLETE. Next begin
the genuine finite long reshape from the actual R110 inlet. Use the
original finite T=400*A singleton and

```
B(Z) = -Lambda*G(Z)+log(Phi_R110)+log(220)/2+log(1+Z^2)
log(Phi(q)) = log(Phi_R110)-B(Z)*sigma(q)-(2/5)*T*q
```

Use the unmodified actual Phi_R110 for B; keep local finite-N corrections
separate. G must be the actual anchored source, not Gbar. Implement a log/factored
profile reader before trying to exponentiate this expression: original T
and G terms can make an ordinary Phi scalar unmaterializable. Retain the
finite window and all six incoming histories; do not shorten it or replace
it with an infinite window. Transfer the same local selected-N source
caps, physical T*dq measure, global phase and nonzero correction memory.

Still open: later long/reference/restoration/outer charts; total radial
and fast-phase derivatives; tighter incoming mean Z3; complete signed
five-defect/selected repair and own pressure; global pressure/common N/C2
joins; exact heat exterior, finite energy and tail; admissible stress and
flat remainder; genuine n>=2 temporal recursion and smooth sum; mean and
two-family oscillatory cancellation; and independent corrected Cartesian
NS validation. This milestone is a local numerical reconstruction, not
the completed global field or a blowup proof.
