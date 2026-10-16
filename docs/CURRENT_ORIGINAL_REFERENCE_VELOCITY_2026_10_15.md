# Actual reference/restoration values to unpatched Rm — 2026-10-15

`OriginalReferenceNumeric` and `OriginalReferenceVelocity` continue the actual
finite-long terminal profile through reference, axial restoration and the
post-restoration interval to Rm. They retain all nonzero leading histories,
original selected-N local/cumulative corrections and the positive-order n1
terminal primitive. Rm here is the unpatched inlet, not a five-moment closure.

```python
from openai_ns_reconstruction.original_reference_velocity import OriginalReferenceVelocity
field = OriginalReferenceVelocity()
Rz = field.reference('.3', '1', log_tau='-2', azimuth='.25')
restore = field.restoration('.3', '.5', log_tau='-2', azimuth='.25')
Rm = field.postrestore('.3', '1', log_tau='-2', azimuth='.25')
near_Rsh = field.at_reference_gap('.3', '1', log_tau='-3', azimuth='1.1')
```

The four methods return normalized profiles/history sectors and physical
cylindrical/Cartesian value sums with full directed coefficient/log errors.
Coordinates remain native; direct physical inversion and gradients are open.

## Original finite geometry and angular profile

Let C=original logCstar, P=original logPstar, T=400*A the original finite
singleton, logref=10*(C+P), and y=log(R/110). The source intervals are

| Chart | Native coordinate | y |
| --- | --- | --- |
| reference | q in[0,1] | T*(1-q)+10*C*q+10*P*q-8*q |
| restoration | t in[0,1] | 10*C+10*P-8+t |
| postrestore | p in[0,1] | 10*C+10*P-7+p |

Reference length is the same source expression L=10*(C+P)-T-8>0.
Its physical measure is L*dq. Restoration and postrestore have unit length
and their native measures occur once. No interval endpoint defines L or T.
The physical reference radius is Rref=110*exp(logref), with Rz/Rrestore/Rm
offsets -8/-7/-6. Small offsets are separate log terms, so they are not
discarded inside the enclosure of the huge reference radius.

After the long endpoint, the angular profile has constant logarithmic slope
-2/5 for Phi and +1/10 for physical-similarity Utheta. Its defining logs are

```
logPhi = Lambda*G - log(1+Z^2) - log(220)/2 - 2*y/5
reference logUtheta = -log(1+Z^2) + (T/10-C)*(1-q) + (P-.8)*q
restoration/postrestore logUtheta = P + offset/10 - log(1+Z^2)
```

At reference q=1, C and T cancel before enclosure. The remaining swirl
shape relative to Pstar is exp(offset/10)/(1+Z^2), which is nonzero and
finite. The actual G, exact source parameters and all finite-N corrections
remain in their original units; no amplitude exponential is materialized.

## Separate j and epsilon histories

During restoration alpha=1-sigma(t), and V=4Z+(j+epsilon*Psi110)*alpha.
The original j and epsilon contributions are stored separately. Combining
Lambda*j with Psi110 in one normalized number would lose the much smaller
epsilon correction numerically.

The additional histories satisfy

```
gamma1_y + gamma1 = alpha
gamma2_y + gamma2 = alpha^2
(Kalpha/Phi)_y + 1.6*(Kalpha/Phi) = 2*alpha
Dpsi_y + Dpsi = Psi110*alpha^2
```

At Rz, gamma1=gamma2=1, Kalpha/Phi=H/Phi and Dpsi=Mpsi. These equalities
follow from the accepted incoming identities K=U0*H+epsilon*Kpsi and
A=U0^2+2*U0*epsilon*Mpsi+epsilon^2*Apsi, with U0=4Z+j and alpha=1 throughout
the preceding actual source. They are not independently chosen initial data.

Each restoration history uses the full finite cutoff kernel
J(k,j;t)=integral_0^t exp(-k*(t-s))*(1-sigma(s))^j ds. Directed positive-cell
enclosures use the original cutoff; scalar centers are quadrature proxies
with complete errors. The mixed source is 2*alpha; its constant-source
equilibrium is 2/1.6=1.25. A source density is not its equilibrium value.

The complete raw leading histories recover as

```
H = Phi*(H_over_Phi)
M = 4Z + j*gamma1 + epsilon*Mpsi
K = 4Z*H + j*Phi*Kalpha_over_Phi + epsilon*Phi*Kpsi_over_Phi
A = 16Z^2 + 8Z*j*gamma1 + j^2*gamma2
    + 8Z*epsilon*Mpsi + 2j*epsilon*Dpsi + epsilon^2*Apsi
B = Phi^2*B_over_Phi2
C = Phi^2*C_over_Phi2  # undressed
```

H/Phi, B/Phi^2 and C/Phi^2 relax toward 5/4, 5/6 and5 with rates1.6,
1.2 and.2. The actual incoming residuals survive. At postrestore alpha=0,
the axial/mixed additional histories decay from their actual restoration
inlets; the three angular histories keep their constant forcing. C remains
undressed and pressure is not exported here.

## Physical values and retained primitive

The physical assembler uses base velocity 4Z+j*alpha, base mean
4Z+j*gamma1 and their common radial recovery. The epsilon*Psi/Mpsi/Qpsi
terms remain separate. The original fixed-time weights and Cartesian rotation
are inherited from the accepted physical long reader.

The nonzero n1 primitive J has physical radial factor exp(-global_y/2),
while its streamfunction has no global_y dependence after analytic radius/
history cancellation. This cancellation is performed before atomic geometry
is enclosed. Local n1 axial/swirl coefficients remain exact zero. No extra
F0/Pstar factor is added to the already assembled angular source log.

## Focused saved evidence

A191.525s focused run exits0 with30 original packets,72 independent raw H/M/K/A/B/C ODE coefficients,279 profile/six-history seam coefficients,207 signed physical coefficient/scale seam checks,63 selected-N/history checks,51 finite-offset/swirl/primitive checks and78 precision comparisons. Five invalid requests are rejected. Saved-report inspection exits0; source/implementation pins total1756/43. Existing read-only GPT-5.6 Luna/max review confirms normalized source2/equilibrium1.25, separated raw moment identities, physical baseline and atomic geometry.

The independent test integrates all six raw defining ODEs under a manufactured
nonconstant axial source and reconstructs M/K/A from the separated sectors.
This checks the binomial cross terms and finite reference/restoration/postrestore
memories, rather than repeating their normalized closed formulas. Signed
coefficient and scale overlap is used at zero crossings; log(abs(coefficient))
is not a valid comparison when its source error interval crosses zero.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_velocity.json`,
9,918,864 bytes, SHA256 `2a162e801d76d25dc3a69dbad3df6dc0236c7e17a76855597b38e28063159144`.
Ordinary reference distances1 and10000, full chart endpoints and source-scaled
time are included. This evidence covers the delivered profiles/physical values;
it does not prove active patch closure or global NS validity.

The normalized scalar history proxies below retain full directed cutoff-integral
errors in the report. They demonstrate the surviving mean/square memory after
alpha has reached zero; they are not substituted source coefficients.

| Stage | alpha | gamma1 | gamma2 |
| --- | --- | --- | --- |
| restoration0 | 1 | 1 | 1 |
| restoration.5 | .5 | .980034117915 | .965881680930 |
| restoration1 | 0 | .607365011340 | .589507470988 |
| postrestore1 (Rm) | 0 | .223437100959 | .216867678994 |

## Remaining work

An additional predecessor gap is now reproduced with the actual saved R110
enclosures: the old long reader rejects finite y=10000 at its
`Mpsi_inlet_memory_Z0` absolute error cap. An ordinary large-looking distance
does not by itself make a huge normalized coefficient absolutely tiny. The
new reference transport evaluates representable ordinary decay distances
below1e30 directly, and exercises gap10000. Extend the earlier long reader
with the same source-preserving treatment before claiming complete numerical
coverage between R110 and Rm; do not weaken or remove its cap proof.

Next connect the actual signed five-defect and implicit active patch at this
unpatched Rm, preserving these sectors and histories. The current original
selected N is unchanged; global high-order stress/heat admission remains open.
Direct physical inversion, total radial/fast-phase rows, incoming meanZ3,
merged modified six moments/own pressure, full Cartesian C2, exact heat,
energy/tails, admissible stress/flat remainder, genuine n>=2/smooth sum,
oscillatory cancellation and independent corrected NS are not completed.
Tasks14/15 remain IN PROGRESS; full goal ACTIVE / INCOMPLETE.
