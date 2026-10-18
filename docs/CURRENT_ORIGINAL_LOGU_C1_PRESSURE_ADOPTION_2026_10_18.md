# Uniform log-u derivative and correlated pressure source adoption

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. The new ordinary-C1
normalized Rm endpoint diagnostic is **TRUE**: both Z0 and Z1 log bounds are
below -1024*ln(2). The saved correlated restoration pressure source has now
been adopted: **TRUE**. These statements do not establish
independent reserved/full-pulse m15 matching or global five-moment closure.

## Same-source logarithmic derivative

The post source has bphi=B_phi/(Pstar*E0)=-C*e*f(t_hat),
f=t/(1+t^2), C=a*(1+2q^2)/2<=2, a=4/5. The frozen positive E0 is a
normalization constant chosen from a whole-chart bound, not a selected
native field value. The source q and dstar are Z-independent, so
u_Z/u=p2_Z/p2. The latter is consumed as the actual signed function
enclosure saved in the preceding signed-native packet.

At fixed target phase, with I=integral_0^psi t^2 and Phi=K*(psi+I),

    W=u*fprime(t)*(t_u-t_psi*I_u/(1+t^2))
    bphi_Z=-C*(e_Z*f(t_hat)+e*W*(u_Z/u)).

The full combination is bounded before assigning separate inverse-phase
caps. Set m=1/(1+u^2), r=u/sqrt(1+u^2), D=1-2r*cos(psi)+r^2 and

    v=sqrt(m)*(cos(psi)-r)/D, h=sqrt(m)*sin(psi)/D, k=m/D
    k=1+2uv, h^2=k-v^2, (v-u)^2+h^2=1+u^2.

Exact source identities are

    t=2qv
    t_u=2qm*(2v^2-3uv-1), t_psi=-2qkh
    I_u=4q^2*(omega/u^3-mh*(3-k)/(2u^2)),

where omega=atan2(r*sin(psi),1-r*cos(psi)) has a positive denominator and
|omega|<pi/2. This removes the apparent singular factors before bounding.
For |u|>=1 and q in [1/2,3/2],

    |u*fprime(t)*t_u|<=4
    P=mk*h^2*(3-k)/(2u)
     =m*(1-v^2)/u+3mv-(4u^2+1)mv^3+2mu*v^4.

After division by (1+v^2)^2, the five coefficient bounds for P sum to
15/2, and the omega part is at most 12. The phase prefactor is at most 27.
Thus |W|<=4+27*(12+15/2)=1061/2<600, uniformly in the entire angle domain.
This is an analytic inequality, not a sampled maximum or an asymptotic
replacement of the actual finite-u function.

The actual source domain is certified without exponentiating logR:
q^2=3/4+5eta/4, eta<1 and the positive square-root source imply the q range;
the signed negative p2 coefficient Kp2 gives

    log|u| >= logR+logP+offset/10-log(1+Z^2)+log(-Kp2)-2ln2-logdstar > 0.

One ln2 is the physical p2 scale R*Utheta/2, the other uses the explicit
lower bound q>=1/2. The lower bound is not a selected q value. The full
post domain Z=[37/100,93/250], p=[0,1] is retained.

The new numerical row is

    logcap(bphi_Z)=ln2+logsum(logcap(e_Z)-ln2,
                             logcap(e)+ln600+logcap(p2_Z/p2)).

Its log upper endpoint is 8.30467165741843545916337.
The previous failed normalized endpoint Z1 upper log was
5.77958630625090433148534e+408906090034569676; the new one is
-1.43059066986408533838418e+408906090034569675.
Old Z0, nonreset M, slow-y/N, native signed Qy, the denominator correction
and every old large bound remain. Pstar*E0 is applied once, and the full
saved logRm is subtracted once. The fast term has no 1/N factor.

Only after both exact endpoint comparisons pass is the symmetric 2^-1024
ordinary-C1 enclosure emitted. It is an error enclosure, not an actual zero,
signed endpoint value or a terminal moment identity. The external Rm RQy
handles [183696,183697] retain coefficient -1.

## Reuse of the existing seam proof

The older pressure source retains uncontracted roots [186084,186085] and
its historical null contraction field. The same immutable source report
also contains actual_restoration_correlated_pressure, with accepted seam
proofs and contracted roots [186096,186099]. Its internal targets are
[171869,171870,178093,178094]; its external Rm pair remains [183696,183697].

The new consumer adopts the saved correlated source under its authenticated
graph namespace and source pins. It keeps the original uncontracted source
alongside that reference and attaches the new Rm C1 certificate only in its
own normalization Wpost=Rm^2. It does not overwrite old evidence, turn a cap
into a signed center, or delete the external atom. No source owner, graph
builder, old reducer or archive reconstruction is replayed. The accepted
source proof is reused as a premise; this is not a fresh archive-level audit.

The complete joined pressure numerical B/D/S1 values and the independent
terminal moment comparison remain unfinished. The adoption packet makes
the proved source dependency available to the next consumer; the small Rm
bound cannot substitute for the remaining functional comparison.

## Checks and saved artifacts

Fresh symbolic coordinate/derivative/W/polynomial residuals are exactly zero;
new finite interval fixtures include negative u=-10^80. A focused boundary
regression repaired the finite evaluator's |u|>=1 check at global precision
53. The original full serialized fixture was not retained; its captured
check summary/hash is saved, other mathematical function ASTs are unchanged,
and only the new boundary regression was run after that precondition patch.
No captured numerical case interval was reconstructed or invented.

The new actual endpoint calculation ran once in 0.203000
seconds at 1600 interval bits and was checkpointed before packaging.
Saved-only inspection passes with 202
packet interval records. The source-adoption calculation is separately
checkpointed and inspected. Old mathematical providers are not imported by
the saved inspectors.

Endpoint report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_post_restore_logu_C1.json`:
657774 bytes, SHA256 `1d778c10c272acafee88bd42178ad4ae8198d231de1f8b12d185a09da75f8072`.
Adoption report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_restoration_pressure_adoption.json`:
674832 bytes, SHA256 `2f10516f66168e513cca3ea16ce25d82bc78a4850c37c30ed867ffda0f558234`.

All 212 cumulative terms, six signed core/lower terms, prior errors and the
reserved n1 own_swirl registration count of one remain. P1 differs from m15;
pressure-Q differs from loop-root-q. No global legacy cap recertification,
full phase value installation or full-NS residual claim is made.

NEXT: consume the adopted pressure source in ordinary-C1 B/D/S1 transport
and the independent reserved/full-pulse m15 comparison, preserving the
external Rm term and its physical scale. Evaluate remaining source centers
or prove their same-source cancellations before boxing. Do not rerun the
accepted log-u endpoint calculation or earlier gamma/history/p2 providers.

Global five moments, exact heat exterior, stress/flat remainder, true
n-dependent recursion, realizable corrections and independent Cartesian
corrected-NS validation remain unfinished. The full-NS 1e-3 target is not
applied to this leading-background boundary diagnostic.

```
python scripts/run_original_restoration_pressure_adoption.py --inspect-report
python scripts/run_original_post_restore_logu_C1.py --inspect-report
```
