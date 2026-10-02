# Original microscopic switches and complete post-switch power mixed4

The SAME actual core/bridge source now supplies physical phase/Z derivatives
through total order four on BOTH original microscopic switches. Each physical
row also has the exact formal conversion to logR derivatives. The complete
original R2..110 power interval supplies ordinary logR/Z mixed4 directly.
Original P0, raw axial V, all actual incoming moment histories and positive
source radii remain unchanged. These are directed enclosures of source
functions, not newly computed point coefficients or a complete NS field.

## Coordinates and width

```text
hb = cstar*K^-100 > 0
s  = log(R/100)/hb
R  = 100*exp(hb*s)
D_y^k = hb^-k*D_s^k, y=logR
R2 = 100*exp(2hb)
```

Neither hb nor hb^-1 is materialized as a source value. A numerical interval
[0,exp(-1000)] encloses the positive width only after a directed log proof.
The k-th logR derivative cancels hb^k against each factored source term
before computing a triangle bound in logs. The looser bound
log|capped D_s^k|-k*loghb is retained as a diagnostic; no finite
phase derivative is silently relabeled a logR derivative. The original
Ra-to-100 prescribed-shear bridge and its inlet/outlet joins remain open.

## Comparison moment histories

The true smoothed comparison Fbar and Vbar become radially constant after
Ra*exp(2hb), well before R100. Their own moment histories continue to vary.
They are enclosed by the existing positive-kernel cumulative bounds, not
replaced with unsmoothed moments, exact point values, or actual-field moments.

The original alpha=1-sigma((y-hb)/hb) gives the IBP representation
alpha(y)*g(y)+integral_0^y[-alpha'(t)]*g(t)dt, for g=logPhi_core and V_core.
All weights are nonnegative, independent of Z, and sum to one. The directed
inequality 2hb<log(4.1/4) keeps the full smoothing path in the accepted
core continuation box. Its axial0..6 jets therefore enclose the comparison
history, including the frozen endpoint, by interval Taylor composition.

For cumulative moments, the exact integrals of those covered histories have
positive masses (1-theta^2), (1-theta), and (1-theta^2)/2 as appropriate,
theta=Ra/R. These masses bound the integrated histories; they do not assert
that the comparison field was constant during smoothing. In the frozen
region the true moment ODE is Y_y=rY-lambdaY*Y, with rates (2,1,2,1,2,1).
Its higher y derivatives are enclosed using the actual moment bounds and
the now constant comparison phi/V. One Z derivative consumes axial order6
to supply direction order5, with all true F0 and F0^2 derivatives retained.

## Source factors retained until the final derivative

The derivative algebra carries a finite sum of source products as exponent
vectors against (loghb, 2logPstar, 2logF0base, 2logu-base). All bases are fixed
at the current Z basepoint; the coefficient jets carry their TRUE Z
derivatives. Bell and Leibniz derivatives, quotient factors, V^2 terms,
radial prefactors and n! factors are combined before each final physical
coefficient is numerically bounded. The final sum ledger records every
source exponent and the combined log.

This resolves the read-only Luna/max finding that width caps were taken
before multiplying some large amplitude/quotient factors. Independent
extreme-factor fixtures include hb^2*Pstar^2=exp(1), which would be damaged
by premature separate capping. Existing parent history enclosures remain
enclosures of their original functions; this stage does not recompute their
signed integral values or claim to restore every lost interval correlation.

The first switch retains

```text
a = hb*Dbar
V_s = -hb^2*(1-sigma(s))*(phi_actual/barphi)
      *(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)
```

with the full multinomial derivative rule. The second switch has
a=hb*(1-sigma(s-1))*Dbar+.8*sigma(s-1), and V_s=0. Dbar is distinct from
the axial drive. All current-radius factors are retained. Physical Ur
uses (D_s+hb/2)^k Q; physical primitives use their original differentiated
RHSs. The axial-quadratic contribution retains Pstar^-2.

## Signed axial source and local functional joins

The accepted common source graph is retained, with a first-switch partial
integral V(s)=V100+I_first_switch(0,s). Its integrand, sign, bound and
derivative recurrence share the original source namespace. The second
switch and full postpower use the SAME V110=V100+I_first_switch(0,1).
Caps do not define any of these functions.

At phase1, identical actual histories/P0 and original sigma endpoint jets
give identical physical derivatives. At phase2, a=.8 and all higher control
derivatives vanish. Every physical phase row is hb^k times the corresponding
postpower logR row, including the primitive RHSs and radial velocity
prefactor. This proves the functional R2 join, rather than inferring it
from overlapping derivative intervals.

On the complete following power region, use the formal coordinates

```text
zeta = fraction*(log(110/100)-2hb), 0<=fraction<=1
theta = exp(-zeta)=R2/R
logR = log100+2hb+fraction*(log(110/100)-2hb)
```

All six moment histories use the original positive power transport from
the actual R2 inlet. The exact angular identity is
log(F/F100)=-(2/5)*log(R/100)+.6hb-.5hb^2*JD.
It agrees with the accepted original post(110) source, and therefore with
the accepted local R110-to-reshape join. The checker proves positive-width
containment and the source log-radius ordering, including fraction0=R2 and
fraction1=110. Numeric radius/ratio intervals are explicitly enclosures.

## Evidence and reproduction

- 480 phase velocity/pressure and 600 phase five-primitive mixed bounds,
  across whole branches, endpoints and interior packets.
- 1080 exact formal phase-to-logR conversion rows and 1080 final source-sum
  ledgers, checked against their factored exponent vectors and coefficients.
- 180 velocity/pressure and 225 primitive mixed bounds on the complete
  postpower region and its two endpoints.
- Independent finite fixtures check 96 comparison-direction derivatives,
  192 switch-control derivatives, 135 full physical primitive/velocity
  derivatives and 18 extreme-factor coefficients. Fixtures validate formulas;
  actual-source admission comes from current family/source prerequisites.
- 62 symbolic source, endpoint, physical prefactor, FTC/Leibniz, IBP and
  positive kernel identities; 4431 directed positive-source cap proofs;
  exact actual R2 histories and original R110 V/P0.

Run `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py
--stage switchmixed` from this checkout using Python with mpmath and sympy.
The ordered pipeline now contains 100 modules. The independent receipt is
`experiments/root_st073/lei_ren_part1_paper_compliant_microswitch_mixed_C4_check.json`.

## Remaining work

Prescribed-shear bridge mixed4, core/bridge and bridge/R100 joins; complete
inner/pre-O3 dispatcher and physical Cartesian assembly; required-domain
and terminal energy; admissible divergence-form stress and independently
bounded flat remainder; genuine n-dependent recursion; mean/oscillatory
correction and full forced NS residual. Their completion gates remain false.
The original unlocalized whole-space energy remains infinite. This stage
does not establish blowup, dynamic scale measurements, or a recursion step.
