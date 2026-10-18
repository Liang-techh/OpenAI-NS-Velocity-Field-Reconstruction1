# Signed native histories and pressure source

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. New source-bound signed
enclosures for prescribed native A, b=B/phi^2 and c=C/phi^2 are installed
through ordinary Z2 on the post-restore chart. The full native pressure
source p2 now has signed scaled ordinary-Z0/Z1 enclosures. Its value is
strictly negative over Z in [37/100,93/250], p in [0,1]. These are function
enclosures over a domain, not a selection of sampled values or midpoints.

The phase/inverse, finite-N modified field and terminal pressure join have
not been evaluated in this packet. The inherited Rm C1 endpoint diagnostic
remains FALSE and was not recomputed. No independent m15/global-five-moment
functional closure is claimed.

## Histories and retained incoming data

The pinned bridge densities/rates give A_y=U^2-A, B_y=phi^2-2B and
C_y=phi^2-C. The C6 source has constant U=U_R110 before restoration,
Lref=logref-T-8, phi_y=-2phi/5 on reference/restoration/post, and
phi=Utheta/(F0*sqrt(2R)). Incoming A_R110/B_R110/C_R110 are never set to
zero or replaced by a cap value. A_R110 itself is not numerically selected.

On the long reshape, the saved bound |Bshape|/T<=1/200 and sigma'<=8 imply
phi_y/phi<=-9/25. Thus for k=0,1,2,

    |(phi^2)_Zk(y)| <= phi110_upper^2 * Dk * exp(-18*y/25)
    Dk = sum_j binomial(k,j) relphi_j relphi_(k-j).

The new endpoint evaluation, with lambda=2 for B and lambda=1 for C, is

    |J_Zk(T)| <= exp(-lambda*T)*incoming_cap_k
               + phi110_upper^2*Dk*exp(-18*T/25)/(lambda-18/25).

This follows from the finite Duhamel integral and 1-exp(-(lambda-18/25)T)<=1.
Its inlet term is retained. Whole-chart phi derivative caps are used only
inside the density bound. They are not endpoint values. A fresh rational
decaying-density proof verifies the discarded positive factor and bound.

Reference transport keeps the exact source decomposition

    A_Rz = U110^2 + exp(-(T+Lref))*(A_R110-U110^2)
    b_Rz = 5/6 + exp(-2Lref)*B_Rsh/phi_Rz^2 - (5/6)*exp(-6Lref/5)
    c_Rz = 5   + exp(-Lref)*C_Rsh/phi_Rz^2 - 5*exp(-Lref/5).

Ordinary binomial differentiation of invphi_Rz^2=2Rz F0^2/Utheta_Rz^2
provides its error multipliers. The saved restoration inverse-Utheta bound
and upper radius are conservative over Rz. All memory/subtraction log caps
are kept in the packet. Each resulting Z0/Z1/Z2 deviation is proved below
2^-1024 before a finite symmetric enclosure is emitted. This radius is
an error enclosure, not an actual deviation, endpoint value or zero.

At Rz the new B normalized deviation log-cap upper endpoint is approximately
-8.0113e+408906090034569674; the C counterpart is approximately
-8.5835e+408906090034569673. No exponential of these enormous logarithms or
of logR is materialized. The old large bounds remain in the immutable
parent and are not deleted by an unsupported endpoint cancellation.

## New restoration scalar

Let h=1-sigma, E=U110-4Z. The source sigma uses reciprocal squares:
sigma=exp(1/(1-s)^2-1/s^2)/(1+exp(1/(1-s)^2-1/s^2)) on the left half and
sigma(s)=1-sigma(1-s) on the right, with flat endpoints.

    gamma_j(1)=exp(-1)+integral_0^1 exp(s-1)*h(s)^j ds.

The saved gamma1 interval is reused intact, with no old quadrature replay.
Only gamma2 was newly integrated. Certified Simpson panels use the whole
panel fourth-derivative interval and remainder length^5*sup|f4|/2880.
The true endpoint tails 0<=integrand<=1 are retained. At 100 interval decimal
digits, 204 panels give

    gamma2(1) in [0.58950747057120941112361358,
                  0.58950747171186115484106431],

with certified width 1.1406517437174507281998325e-9 < 1e-8. Exact directed
tuples, tail/quad errors and panel budget are serialized. C6/Gamma/stable
sigma/flat derivative source bytes and ASTs are pinned. No library quadrature
estimate is presented as a certificate.

The newly evaluated post-restore history is

    A(p)=16Z^2+exp(-p)*(8ZE*gamma1+E^2*gamma2)
                    +exp(-(p+1))*epsilon_A_Rz.

Its ordinary Z0/Z1/Z2 rows include the full signed E jet, saved gamma1
uncertainty and faded A_R110 memory. b and c retain their certified memory
radii over the whole post chart. M is read from saved corrected endpoints;
their hull is valid because each ordinary-Z row is affine in exp(-p).
The old M producer is not rerun. At Rm the new A value enclosure is
approximately [2.1850613316725477,2.2194826683274523].

## Source-correct signed p2

The exact bridge generic_fields and inlet generic_functions are the source
of the signed formula. Older textual majorant expressions are not treated
as the signed numerator. With d=1-Z^2, L=1-delta*Z^2,

    H=-U+(1-delta)*(M-Z*M_Z)/2+U*((1-delta)*Z*M+d*M_Z)
      +2delta*Z*A-d*A_Z
    P=2(1+delta)*Z*p0-d*p0_Z
    g=(Utheta^2)_Z/Utheta^2=-4Z/(1+Z^2)
    S=-2delta*Z*b+d*(b_Z+g*b)
      +2(1+delta)*Z*c-d*(c_Z+g*c)
    K=p2/(R*Utheta/2)=(S+2H/Utheta^2+2(Pstar^2/Utheta^2)*P)/L.

The F0_Z terms cancel against the SAME phi normalization inside the
differentiated F0^2 B/C products. There is no uncancelled Lambda*G_Z swirl
term. The independent normalized p0 and its ordinary Z0/Z1/Z2 rows are
copied from actual_micro_b_p2_mixed_bounds.whole_tube_pressure.p0_jets.
Its flattening uncertainty is retained. Pstar^2 is applied once through
exp(2(logPstar-logP)-offset/5)*(1+Z^2)^2. Distinct saved logP/logPstar
enclosures are retained; no identity is inferred from interval overlap.
The finite hydro term 2H/Utheta^2 remains even when extremely small.

The ordinary derivative K_Z and the physical p2_Z scaled coefficient
K_Z+(logUtheta)_Z*K are both recorded. At Rm,

    K in [-32.986884507109012656713005,-32.400839355543820451024280]
    K_Z in [-59.124770766415126624938811,-56.398457064415970209710082]
    p2_Z/(R*Utheta/2) in [-38.062790693695879750634563,
                         -34.811473096354482465322295]
    p2_Z/p2 in [1.0553125466835842655635048,1.1747470575074251275050473].

The physical outer R is present once through the positive scale R*Utheta/2;
R is not materialized. Over the full post chart K lies approximately in
[-42.1411,-32.4008]. Fresh independent rational source expansion,
factorization, nonzero F0 derivative, ordinary Z derivatives, directed IV
containment, two-row hydro and no-double-Pstar fixtures pass. No unknown
hydro Z2 row is filled with zero to obtain a claimed derivative.

## Execution and next dependency

Three new calculations ran once and were checkpointed before packaging:
gamma2, signed histories, signed p2. Total new mathematical runtime was
2.860000 seconds. Saved-only inspection passed
with 140 interval records, 2037
input pins and 460 implementation pins;
it imports no mathematical provider. No accepted ancestor computation ran.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_signed_native_sources.json`:
566530 bytes; SHA256 `9887a353aeab1b24f301a5b6c240f83f58a5d54d001eaac6d16f932b87979561`.
Native immutable source SHA256 `10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a`.
The 212 cumulative terms, six signed core/lower terms, all errors and
external Rm handles [183696,183697] with coefficient -1 remain. Reserved n1
own_swirl is registered once against negative full-pulse m15; P1 is different
data and pressure-Q is different from loop-root-q.

NEXT: feed these signed scaled p2 and p2_Z/p2 ports into the same-source
closed phase/contraction. Keep the enormous positive u scale symbolic or
logarithmic; direct r close to -1 subtraction is not a reliable numerical
input. Obtain an interval-safe scaled/limiting contracted derivative before
reevaluating the fast C1 endpoint. Then join pressure sources and compare
reserved/full-pulse m15 as ordinary-C1 functions of Z. Neither sign nor a
small endpoint enclosure establishes that terminal identity.

Global five moments, exact heat exterior, stress/flat remainder, true
n-dependent recursion, realizable corrections and independent Cartesian
corrected-NS validation remain unfinished. No full-NS 1e-3 threshold is
applied to this leading-background source calculation.

```
python scripts/run_original_signed_native_sources.py --inspect-report
```
