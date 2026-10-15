# Actual original-N first-switch velocity and retained history — 2026-10-15

`OriginalFirstSwitchNumeric.evaluate(Z,s,...)` continues actual modified
leading+n1 physical velocity enclosures from R100 through the original
first switch, native s in[0,1], ending at R1=100*exp(hbS).

```python
from openai_ns_reconstruction.original_first_switch_numeric import OriginalFirstSwitchNumeric
field = OriginalFirstSwitchNumeric()
inlet = field.evaluate('.3', '0', log_tau='-2', azimuth='.25')
inside = field.evaluate('.3', '.5', log_tau='-2', azimuth='.25')
exit = field.evaluate('.3', '1', log_tau='-2', azimuth='.25')
```

## Original unequal controls and source scales

The actual source equations remain

```
R = 100*exp(hbS*s), rho = 4*exp(Y+hbS*s), dy = hbS*ds
Y = log(25)+logLambda
(logPhi)_s = -hbS^2*R*D_over_R/2
U_s = -hbS^2*(1-sigma(s))*(Phi/barPhi)
      *(R*drive_hydro+R*Pstar^2*drive_pressure+R^2*F0^2*drive_swirl)
Psi_s = U_s/epsilon
```

The angular weight is1; only the axial source uses the original flat
1-sigma. Its prefix mass is evaluated by the accepted original pulse
integrator and equals1/2 at s=1. The angular prefix mass is the full s.
Replacing both controls by a half-mass changes the construction.

Comparison moments continue at the current R. Their R0 modes are converted
once to R100 modes via R0^j*100^(p-j); exp((p-j)hbS*s) retains their current
radius variation. Positive sums of modal norms bound every source prefix,
including signed hydro, pressure and swirl contributions before epsilon
normalization. The actual Phi/barPhi quotient keeps its complete macro
and switch exponential error. No comparison field replaces the actual
incoming field, and no endpoint signed integral is assumed to bound all
earlier prefixes.

Two width factors are retained in each local velocity source. Angular
uses hbS^2 directly in its norm bound; axial pressure/F0/epsilon factors
are normalized before the original width-log proof. The physical Jacobian
is applied once. Neither hb nor F0 is assigned a representable bound or
midpoint as its source value.

The accepted C6 geometry certificate checks that bridge hbB and switch
hbS are distinct native source nodes but exact reboxings of the same K1
width source. The certificate, full receipt/input hashes and shared K1
ledger are consumed before bridge log-width bounds can be used for hbS.
Their distinct roles remain explicit in geometry and selected-N phase.

## Actual incoming mean and finite-N history

The numerical background starts from the UNMODIFIED actual directed macro
endpoint. Local selected-N corrections are added separately; incoming
selected-N errors are therefore not added twice.

The original mean equation is

```
Mpsi(s) = exp(-hbS*s)*Mpsi100
          + hbS*integral_0^s exp(-hbS*(s-t))*Psi(t) dt
```

The nonzero actual Mpsi100 persists. Subtracting this common inlet gives
the conservative change bound

```
hbS*s*(norm(Psi100-Mpsi100)+uniform_norm(Psi-Psi100)).
```

This retains the true positive history kernel while bounding it by1;
the mean is not restarted or replaced by a comparison moment. The core
proxies have complete actual macro and local source errors. Phi/Psi have
Taylor Z0..3; own mean stays Z0..2 and Qpsi Z0..1. The inherited first-micro
mean Z3 modification error remains unresolved.

The canonical N is still `2^(2^1358356628656378313)`. The original phase is
`frac(N*(Y+hbS*s-hbB*sc/2))`; the origin shift belongs only to the periodic
correction, never to R or rho. Local A/B and five history-density caps are
the accepted full-phase magnitude enclosures, not signed point values.
All earlier micro and macro correction memories and the current switch
prefix enter once, using their respective physical measures. The
rate-zero pressure correction memory remains a bound, not an exported
completed pressure field.

The n1 f/u support is already finished. Its terminal primitive persists;
the directed mean uses primitive/current_rho, and q1 retains weight1+delta.
Its point coefficient uses the R100 radius proxy with the full original
hbS*s error enclosed. The physical reader retains the true Y+hbS*s log
weight. This approximation does not select a new width or assert exact
signed point evaluation.

## Focused run

Run `python scripts/run_original_first_switch_numeric.py --write-report`,
then `python scripts/run_original_first_switch_numeric.py --inspect-report`.
The report is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_switch_numeric.json`.

The driver checks representable finite-h source quadrature with the true
flat sigma, genuine axial derivatives through Z3, unequal angular/axial
mass, actual geometry and nonzero mean memory. Original packets check
the accepted R100 seam, cumulative histories, terminal primitive/current
rho, endpoint semantics and precision60/85 overlaps. The finite fixture
parameters test the numerical bounds; they never replace original scales.

The accepted run exits0 in280.202s with10 original packets,24 independent
unequal-control/full-source/Z3/geometry/mean checks,153 selected-N small-error
proofs,63 R100 native/physical coefficient seams,69 retained correction and
terminal-primitive memory checks,15 current-rho/endpoint checks and26
precision60/85 coefficient/log comparisons. Seven unsupported requests are
rejected. Saved-report inspection exits0 in2.307s. The report pins1752
source inputs and34 implementation files.

At s=.5 the original axial prefix mass lies in
[0.4793661109803435,0.4793661478776142], while the angular mass is.5.
At s=1 the axial mass is exactly.5 and the angular mass is1. These distinct
weights are retained before the two native width factors and physical
normalization are applied.

## Remaining work

Tasks14/15 remain IN PROGRESS; full goal ACTIVE / INCOMPLETE. Next continue
second_switch with its original angular blend of hbS*Dbar and 4/5*sigma,
and frozen unmodified actual U from this endpoint. The local selected-N
B/N correction still evolves. Carry all actual mean and five correction
histories. Then post_power must reach R110 using the exact
length log(11/10)-2hbS and angular exponent -2/5, without resetting memory.

Also open: total radial/fast-phase derivatives and physical vorticity
beyond r_minus; tighter incoming mean Z3; later charts and complete signed
five-defect/selected repair/own pressure; global pressure/common N and
Cartesian C2 joins; exact heat exterior, finite energy and tail;
admissible stress/flat remainder; genuine n>=2 temporal recursion and
smooth sum; oscillatory stress cancellation; and independent corrected
Cartesian NS validation. These velocity enclosures are not a completed
global field or a proof of the full blowup construction.
