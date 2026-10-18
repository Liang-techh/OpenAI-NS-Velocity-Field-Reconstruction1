# Correlated fast-phase derivative and closed post-restore phase

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. A new signed-source
derivative interface and closed phase primitive are implemented. The new
Rm ordinary-C1 normalized endpoint diagnostic is **FALSE**.
This remains separate from signed pressure joining, independent terminal
m15 matching, global five moments and the final full-NS residual test.

## Same-source derivative

The checked C5 source gives B_phi/Pstar=-a*E*t/(4*pi*lambda), b=t0=0,
lambda=K*(1+t^2), K=1/[2*pi*(1+2*q^2)]. Let e=E/E0,
C=a*(1+2*q^2)/2 and bphi=B_phi/(Pstar*E0). Frozen Pstar/E0 and constant
a/q/K give bphi=-C*e*f(t_hat), f=t/(1+t^2). At fixed target phase,
Psi_Z=-Phi_Z/lambda and

    bphi_Z=-C*(e_Z*f(t_hat)+e*fprime(t_hat)*(t_Z+t_psi*Psi_Z)).

The signed interface keeps lambda*t_Z-t_psi*Phi_Z as one source expression
before taking an absolute value. An independent symbolic chain-rule proof
has exact residual zero. The C5 source byte and function AST hashes bind
the formula to its B_over_Pstar_phi, direction_at_inverse and
phase_lambda_at_inverse ports. Pstar*E0 is applied once.

The generic formula constructor accepts exact expressions for E/p2 and
requires frozen source parameters, a=4/5 and t0=0. It rejects log-cap
mappings, binary/SymPy floats and dependent frozen scales. It is not a
certificate that supplied rational expressions are the actual native data.
No actual signed native p2 or inverse-phase value has been installed.

## Closed phase function

For finite signed inputs, a separate interval evaluator implements
u=p2*q/dstar, r=u/sqrt(1+u^2), m=1-r^2, D=1-2r*cos(psi)+r^2 and
t=2q*sqrt(m)*(cos(psi)-r)/D. Here u is the loop parameter. For r!=0,

    omega=atan(r*sin(psi)/(1-r*cos(psi)))
    I=2*q^2*(psi+m*sin(psi)/(r*D)+(3*r^2-1)*omega/r^2)
    Phi=K*(psi+I), lambda=K*(1+t^2).

I is the closed source primitive of t^2. The evaluator retains I_u,
t_u, t_psi and the signed contraction
(1+t^2)*t_u-t_psi*I_u, with the exact r=0 removable limits. This provides
a computable phase unit for finite signed parameters without rebuilding
the old integral graph. The actual native parameter values remain a
separate dependency. A non-singleton r interval crossing zero is explicitly
rejected; the direct r^-2/r^-3 representation is not a uniform near-zero
interval implementation. Positivity checks must resolve on the supplied
finite interval before its result is returned.

Independent symbolic residuals for I_psi, I_r, t_u and t_psi are exactly
zero, using cos^2(psi)+sin^2(psi)=1 and the stated positive-denominator
domain. The r=0 formulas are separately recorded. These exact identities
are distinguished from the finite-input numerical consistency fixtures.

Only the explicit psi="2*pi" token invokes the exact endpoint identity
Phi=1 and I_u=0. An ordinary interval 2*pi remains an interval and is not
collapsed to a selected value. Global-precision-53 regression checks retain
high-precision input tuples. The r=0 and nonzero-r fresh fixtures compare
the primitive/derivatives and fixed-phase contraction independently.

## New actual error transport

The derivative identity cancels the compatible lambda_Z terms before
bounding. Independent rational inequalities |f|<=1/2 and |fprime|<=1,
plus C<=2, give

    logcap(bphi_Z)=ln(2)+logsum(logcap(e_Z)-ln(2),
                             logcap(e)+logcap(t_hat_Z)).

The adapter consumes saved y0/Z1 of t_hat at fixed phase, and saved e/e_Z.
No inverse/phase producer is rerun and no log cap is used as a signed value.
Only the fast Z1 row is changed; the value row, slow-y/N term, nonreset
mean memory and all other mixed rows remain exact copies. The fast term
still has no 1/N factor.

The new endpoint Z1 numerator uses universal coefficient bounds on the
already source-bound delta/Z domain, its positive L lower floor, and the
retained denominator-derivative correction. Native signed Qy, old Z0,
full directed logRm and external handles [183696,183697] with coefficient
-1 remain. Qy/Rm is normalized once by Wpost=Rm^2. The old large cap is
retained alongside the new proven bound; it is not silently deleted.

The old combined Z1 log-cap upper endpoint was 7.83963687085518704816048e+408906090034569676; the new upper
endpoint is 5.77958630625090433148534e+408906090034569676. The exact-tuple local C1 result is false.
A symmetric 2^-1024 enclosure is emitted only if both rows pass, and is
explicitly an enclosure rather than an actual signed value or zero.

The successful new saved-data calculation took 0.016000
seconds at 1600 interval bits. Saved-only
inspection passed with 76 packet interval
records, 2033 input pins and
455 implementation pins. Fresh fixtures
include exact symbolic cancellation, actual runtime cap propagation,
synthetic serialization/endpoint integration, safe interval coercion and
the closed primitive. The approximate signed-cancellation witness uses new
synthetic rational data and is explicitly not a certified native value.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_fast_phi_correlated_C1.json`:
497555 bytes; SHA256 `8b77d380fcfd4837c4bcfff10a1e4fe512b37d47b30340c5093552bb3b1ca108`.
Parent `lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json` remains immutable with SHA256
`78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb`. Its native, loop, finite-N and old
signed graph data remain preserved. The 212 cumulative terms, six signed
core/lower terms and all errors remain. Reserved n1 own_swirl stays
registered once against negative full-pulse m15; P1 is different data.

## Remaining source work

The saved native root/phase packets contain absolute caps for p2, u/r/alpha,
E/E0, t/Phi and the inverse. The saved signed six-I1 pressure-Q rows cannot
serve as loop-root q or p2. Obtain signed native history/root function ports
from the actual C6 source and its incoming histories on the strict Z and
chart domains. Then apply the signed contraction or closed phase with
those actual inputs, retain internal endpoint source identities, join the
pressure sources and compare reserved/full-pulse m15 as an ordinary-C1
function of Z. A small endpoint enclosure cannot establish that identity
by interval overlap.

Global five-moment closure, exact heat exterior, stress/flat remainder,
true n-dependent recursion, realizable corrections and independent
Cartesian corrected-NS validation remain unfinished. No 10^-3 full-NS
threshold is applied to this leading-background calculation.

```
python scripts/run_original_fast_phi_correlated_C1.py --inspect-report
```
