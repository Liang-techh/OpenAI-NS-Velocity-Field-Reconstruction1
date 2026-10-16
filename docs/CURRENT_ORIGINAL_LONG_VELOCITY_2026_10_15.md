# Physical finite-long Cartesian velocity — 2026-10-15

`OriginalLongVelocity` assembles the original finite-long leading+n1 field
as physical cylindrical and Cartesian exponential sums. It consumes the
accepted actual finite-long numerical reader, its six unmodified histories,
same canonical selected N, local/cumulative corrections and terminal n1
primitive. No source profile, incoming history or original scale is reset.

```python
from openai_ns_reconstruction.original_long_velocity import OriginalLongVelocity
field = OriginalLongVelocity()
middle = field.evaluate('.3', '.5', log_tau='-2', azimuth='.25')
near_inlet = field.at_log_radius('.3', '1', log_tau='-3', azimuth='1.1')
late = field.evaluate('.3', '.5', log_tau=None,
                      scaled_time='1e18', azimuth='1.7')
u, v, w = (middle['physical']['Cartesian'][name] for name in ('u', 'v', 'w'))
```

Each component is a list of terms. Its meaning is the sum of
`coefficient * exp(sum(log_scale_terms))`, with complete directed source
enclosures and serialized-number error bounds. Original huge/tiny factors
remain separate. The result is not an ordinary float vector. Native Z,q or
Z,y is required; direct physical Cartesian inversion on this chart is open.

## Physical weights and coordinates

Let y=T*q, R=110*exp(y), rho=Lambda*R, epsilon=1/Lambda,
lambda=sqrt(tau/(1-Z^2)), tau=1-t, and U0=4Z+j. Leading fields are

```
uz = lambda^(-1-delta) * (U0 + epsilon*Psi)
ur = lambda^(-1)*sqrt(R/2) * (Qaxis + epsilon*Qpsi)
utheta = lambda^(-1-delta) * exp(actual_similarity_log_Utheta)
Qaxis = ((1+delta)*Z*U0 - 4*(1-Z^2))/(1-delta*Z^2)
Qpsi = (2*Z*Psi - (1-delta)*Z*Mpsi - (1-Z^2)*Mpsi_Z)/(1-delta*Z^2)
psi0 = lambda^(1-delta)*R*(U0 + epsilon*Mpsi)
```

Psi, Mpsi and Qpsi include the actual selected-N modification with its
complete cumulative source errors. The current swirl source log already
contains F0, current Phi_N and sqrt(2R) once. No additional Pstar, F0 or
radius factor is applied to this term.

Physical coordinates are r=lambda*sqrt(220)*exp(y/2),
x=r*cos(azimuth), y_physical=r*sin(azimuth), and
z_physical=Z*lambda^(1-delta). Cartesian velocity uses
u=ur*cos(azimuth)-utheta*sin(azimuth),
v=ur*sin(azimuth)+utheta*cos(azimuth), w=uz.

Time accepts finite decimal log(tau)<=0, or the source-scaled coordinate
S=-delta*log(tau)/2>=0 with log_tau=None. The accepted core time map
cancels delta exactly before forming its scaled-time products.

## Actual n1 primitive and cancellation before enclosure

The actual terminal primitive J(Z) persists after the local n1 f/u profiles
vanish. Its three Taylor Z rows determine Q_J through Z1:

```
mean_u = J/(110*Lambda) * exp(-y)
Q_J = (-(1+delta)*Z*J - (1-Z^2)*J_Z)/(1-delta*Z^2)
ur1 = Q_J * lambda^(-1+2delta) * Lambda^(-2)/sqrt(220) * exp(-y/2)
psi1 = J * lambda^(1+delta) * Lambda^(-2)
uz1 = utheta1 = 0
```

The radial memory exp(-y) combines with sqrt(R) before interval arithmetic.
The streamfunction memory cancels R completely. Consequently the n1 mass term
contains no y log, and the n1 radial term contains exactly one -y/2 log.
Large opposite logs are never boxed separately and subtracted later. Both
the nonzero J and its source enclosure remain in the physical packet.

The meridional field uses the shared actual mean/streamfunction recovery.
Numerical centers remain bounded proxies of the defining source functions;
this value reader does not certify derivatives of those proxies or a full
Cartesian C2 join.

## Focused physical evidence

A169.944s focused run exits0 with11 actual physical packets,44 independent fixed-time streamfunction/Cartesian/scaled-time checks,81 geometry/cancellation checks,38 R110 physical term-log seams,25 history/integer/coordinate checks and102 precision comparisons. Five invalid physical time/angle requests are rejected. Saved-report inspection exits0; source/implementation pins total1755/40. Existing read-only GPT-5.6 Luna/max review confirms all Lambda/time/radius and primitive weights, with no material algebra issue.

The independent check differentiates a manufactured physical streamfunction at
fixed r,z,t, using an independently solved axial inverse. It checks all cylindrical
and rotated Cartesian values, the streamfunction and coordinates under ordinary
and source-scaled time. Actual original-scale checks compare term source
log magnitudes across R110, canceling the common huge anchor before the swirl
comparison, and compare coefficient/log enclosures at two precisions. These
checks cover the assembly and the stated value seam; they do not prove full
field gradients, global moment closure or useful combined relative accuracy.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_velocity.json`,
3,521,643 bytes, SHA256 `f57ba9eca3f456d0814ac7149a2ca5a48067a533e5a29e2cdc52bdfcf3c2a0c6`.

## Limits and next work

This delivers physical value terms on the full original finite-long native
window, including its endpoints. It preserves the original finite T and
selected N and the leading/n1 temporal orders0/1. It does not add genuine
temporal n>=2 recursion or an admitted infinite sum.

Next continue reference/restoration from the actual long terminal profile,
with its six history ratios and current original scales. Direct physical
inversion, total radial/fast-phase derivatives, tighter incoming meanZ3,
merged modified six moments, own pressure, signed five-defect/selected repair,
global common N/C2/heat/energy/stress/flat remainder, genuine n>=2/smooth sum,
oscillatory cancellation and independent corrected NS remain open.
Tasks14/15 are IN PROGRESS; full long-term goal ACTIVE / INCOMPLETE.
