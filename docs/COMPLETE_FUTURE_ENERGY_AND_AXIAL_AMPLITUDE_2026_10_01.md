# Complete future swirl energy and actual axial amplitude (F38-C5)

The compliant angular/pressure repair now feeds the complete corrected
future swirl energy. This selects a unique smooth positive axial amplitude
on the entire axial interval Z in [-1,1], with

`1.0086895652225 <= ap(Z) <= 1.0114075201817`.

The actual end coefficients c1(Z), c2(Z) are defined by the same two linear
moment equations and this selected amplitude; c1 is negative and c2 is
positive. Their formal nonzero exponential scale is retained. The partial
axial pulse field, full five-primitive assembly and whole-outer cone are
subsequent work. No temporal recursion is claimed.

## Full corrected future energy

`CompliantFutureSwirlEnergy.future(Z)` returns a C1 enclosure of

`integral_Rv^infinity Utheta_corrected(R,Z)^2 dR / (Rv*Utheta(Rv,Z)^2)`.

It includes the entire following sequence:

1. The actual 100-unit flatten, with factor ((1+Z^2)/2)^sigma.
2. The -30log(mu) post-flatten power buffer.
3. Both actual angular bumps, including their signed linear and positive
   quadratic energy changes from the solved coefficients.
4. Both steep unit transitions and the Ts-long exact R^-3/2 segment.
5. The selected same-source waiting interval.
6. The epsilon collar and the infinite exterior, retaining separate
   epsilon/epsilon^2 atoms and the positive Gamma heat energy deficit.

All pieces use relative radii. Neither Rref nor exp(-1/mu) is materialized.
The power/waiting integrals use positive cancellation-free forms. The heat
deficit retains its exact formal factor a*S with S=1/Rtail; the strong finite
S cap only supplies a directed enclosure. It is not substituted for S.

The whole-Z normalized future energy lies between approximately
7.0615580051106e18 and 2.8246232020442e19. Its Section 7.34 weight is
`mu*exp(-26)/2`, because

`Rv*Utheta(Rv,Z)^2/(Rp*Utheta(Rp,Z)^2)=exp(-26)`.

The resulting weighted future energy is between approximately
4.5614516127117e-408906090034569130 and
1.8245806450847e-408906090034569129. These are positive arbitrary-exponent
quantities; they have not been replaced by zero. The integral is finite
because the actual delta is positive. This proves finiteness of the swirl
radial tail in these units, not the full three-dimensional physical kinetic
energy or time integral.

## Actual selected amplitude and affine end corrections

Write the end coefficients in the common pulse scale as

```text
cj(a,Z)=exp(log_end_scale)*(uj(Z)+vj*a),
nu_j=mu*Kj*exp(2log_end_scale)>0.
```

The pulse slopes vj depend only on the fixed parameters and shapes. The
incoming functions uj use the true m1, m2 jets from the repaired outer
inlet. Their row normalization is essential:

```text
Qi(Z)=mi(Z)*exp(-13lambda_i/mu-common_logpref),
lambda_i=.5-i*mu.
```

The same fixed row factor multiplies mi and its first Z derivative. The
implementation explicitly carries both actual source jets through a
proved finite positive factor cap over all Z. It does not use unscaled mi
in the row-normalized inverse, and the boxes are enclosures of the true
functions rather than definitions of replacement functions.

The energy equation is the actual Section 7.34 equation:

```text
Kpulse*ap^2+sum_j nu_j*(uj+vj*ap)^2
  =(1-exp(-26))/4-mu*incoming_energy+weighted_full_future_energy.
```

Kpulse is approximately .245; it is distinct from the pressure bound Kp=17.
Expanding gives A2*a^2+A1(Z)*a+A0(Z)=0, with A2>0 and A0<0 on the whole
axial interval. The completed-square positive root avoids repeated interval
dependence in sqrt(K)/K. The denominator for the implicit derivative is
bounded between approximately .49303207686278 and .49702831036357.

The residual at a=.9 is negative (between -.052042176289733 and
-.050973928231822), and at a=1.2 is positive (between .10192501992837 and
.10382412758688). This proves the intended positive branch, rather than
choosing a nominal a=1. The actual first derivative is enclosed by

`ap_Z=-(A1_Z*ap+A0_Z)/(2*A2*ap+A1)`.

The end-energy factors are super-small but positive. Each is admitted by a
directed log margin against a finite cap, without exponentiating the
unrepresentable formal log_end_scale. Their dependence on ap and the actual
incoming source remains in the quadratic equation. No evenness assumption
is imposed on ap: cross terms involving the incoming axial moments can
produce a small asymmetric contribution.

## Evidence and reproduction

New paired modules, all prefixed `lei_ren_part1_paper_`, are:

- `compliant_future_swirl_energy.py/.json` and its `_check.py/.json`.
- `compliant_axial_amplitude_selection.py/.json` and its `_check.py/.json`.

The energy checker verifies eleven independent physical normalization and
tail identities, whole-Z C1 positivity, and retained signed, epsilon and
Gamma contributions. The amplitude checker verifies seven independent
affine/quadratic/C1 identities, whole-Z positive-root and denominator gates,
end-coefficient signs, positive formal caps and independent source-jet row
scaling. Every receipt binds the same compliant pressure and five-moment
family; existing trial/legacy receipts are preserved.

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage energy
```

This ordered stage regenerates the energy pair, then the selected-amplitude
pair. The full `--stage all` pipeline also includes the preceding source,
inner, outer and angular stages.

## Next work

Install the actual compactly supported axial pulse with these ap/c1/c2
functions and its partial moments. At Rv, prove Mz=Mtheta_z=0 and
Mztheta=one half of the positive future swirl energy; zero Mztheta is only
the infinity target. Then compose the five primitives through flatten,
angular repair, steep/waiting and exact heat, recover Ur and pressure, and
prove all whole-Z terminal identities and interface smoothness.

The C4 repair/pulse/cone bounds, global admissible stress, flat remainder,
genuine order-dependent temporal recursion, oscillatory correction and
full Cartesian NS residual remain unfinished.
