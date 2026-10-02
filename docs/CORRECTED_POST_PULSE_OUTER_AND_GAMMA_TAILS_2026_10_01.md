# Corrected post-pulse outer field and Gamma tails (F38-D1)

The selected pulse now feeds one corrected post-Rv provider covering
flatten, post-flatten power, both actual angular bumps, steep entry/power/exit,
waiting, the epsilon collar and the exact Gamma exterior. Every region
returns the same five cumulative primitives, the inherited analytic P0,
and C1 axial enclosures. This completes the callable post-Rv assembly layer;
absolute global terminal identities, C4/cone and time recursion remain open.

All axial and mixed histories are inherited from the solved pulse.
`Uz=Ur=Mz=Mtheta_z=0` holds throughout this outer region by the pulse's
terminal identities and paper (3.9). It is not a moment reset. In this
region Ur and every derivative of Ur are exactly zero; the pulse/interior
still lacks the higher axial jets needed for a full global stress evaluator.

## Fields and coordinate charts

`CompliantCorrectedOuterField` supplies:

- `flatten(Z,phase)`, for the entire100-unit flatten.
- `power_buffer(Z,phase)`, up to four log units before Rrel.
- `angular_end(Z,s)`, with `-4<=s=log(R/Rrel)<=0`, including both actual
  implicit angular/pressure corrections.
- `steep_in(Z,t)`, `steep_power(Z,phase)`, `steep_out(Z,t)`, and
  `waiting(Z,phase)`, with the original unit transitions, Ts and waiting root.
- `collar_heat(Z,t)`, for every `t=log(R/Rtail)>=0`. The exact Gamma
  velocity applies from t=3; offsets beyond3 use current-radius normalizations
  and formal log factors instead of exponentiating an unrepresentable radius.

The angular bump evaluator handles interval boxes that intersect a support
even when their lower endpoint lies outside it. Fixed raw coordinates and
direct remaining bump integrals avoid one-minus-CDF cancellation.

The velocity is represented relative to `Ev0=Utheta(Rv,0)`. The Ev0/Pstar
scale retains its separate inlet, inverse-mu and finite log terms. Gamma
packets likewise retain the waiting endpoint, epsilon normalization, current
power and log-bracket terms separately. Formal `S_current=exp(-t)/Rtail`
is never replaced by its finite positive cap.

## Energy and pressure integration

The remaining corrected swirl energy and pressure integral are computed
directly from the current radius through the entire infinite exterior.
Each region uses its remaining transition/power/waiting pieces and the
same Gamma tail. Tiny epsilon and epsilon-squared atoms and the Gamma
deficit are retained separately. Energy corrections from both angular
bumps retain their signed linear and quadratic terms.

Normalized by the current velocity, the energy moment is

```text
Mztheta/(R*Utheta^2)=.5*integral_R^infinity Utheta_corrected^2 dR /(R*Utheta^2)>0.
```

Its initial value agrees with the pulse's selected positive terminal target.
The change from Ev(Z) units to Ev0 units uses the full C1 jet of `(1+Z^2)^2`,
not only its value. Since delta>0, the radial swirl tail is finite and this
energy moment has the zero cumulative infinity target. This is not a
certificate of the full three-dimensional or time-integrated kinetic energy.

The pressure prefactors and energy prefactors differ:

```text
Ns=exp(-1-mu);       Ps=exp(-2-mu)
Nq=Ns*exp(-2Ts);     Pq=Ps*exp(-3Ts)
Nt=Nq*exp(-1-a);     Pt=Pq*exp(-2-a)
Ntail=Nt*exp(-delta*waiting)
Ptail=Pt*exp(-(1+delta)*waiting), a=delta/2.
```

Pressure remains the forward cumulative quantity `P=P0+Mp`. With pressure
integrals in Ev0-squared units,

```text
Mp(R)=Mp(Rv)+Ev0^2*(Prv-Pr(R))
Pforward-Pbackward=P0+Mp(Rv)+Ev0^2*Prv,
Pbackward=-Ev0^2*Pr(R).
```

The constant difference is explicitly retained with its original source
atoms. The backward integral is an independent diagnostic until the
absolute preheat-source equality has been analytically composed. At Rv the
empty pressure increment is exactly zero, including its axial derivative.
Inconsistent negative-upper prefix enclosures fail rather than being clamped
to zero. No pressure datum is fitted or changed.

## Angular history and outstanding constants

The angular moment is propagated forward from the actual pulse through all
regions and the solved bumps. In the heat region it is represented as the
renormalized Gamma target plus the inherited constant defect:

```text
Xforward=Xheat_target+Ctheta*exp(-(1-a)*t)/K(t,Z)
Ctheta=(1-epsilon)*X(Rtail,Z)-Aheat_target(Rtail,Z)
X=Mtheta/(sqrt(2)*R^(3/2)*Utheta).
```

This retains an actual cumulative moment rather than replacing it by a
target. Proving Ctheta=0 and the pressure constant above equal to zero is
the next dependency. The matching equations are already defined; an
independent equality chain must bind their exact source normalizations to
the assembled preheat datum and the inner reference moment.

The Gamma targets use the true positive integral

```text
H_delta(xi)=Gamma(1+a)^-1*integral_0^infinity exp(-v)*v^a*(1+xi*v)^-a dv.
```

Current-radius finite remainder bounds enclose H, H_Z and the three local
tail targets. No convergent infinite Taylor series is assumed. The angular
target also has the exact positive Gamma expectation
`Aheat=E_Gamma[(1+xi*v)^(1-a)]/(1-a)`, supplying an independent finite-parameter
check of the normalization.

## Evidence and reproduction

The paired files use prefix `lei_ren_part1_paper_`:

- `compliant_corrected_outer_field.py/.json`.
- `compliant_corrected_outer_field_check.py/.json`.

Checks cover21 independent physical normalization/ODE/constant-offset
identities, seven directed region interfaces, whole-Z inherited zero
histories and positive energy, unchanged P0+Mp, exact empty pressure history
at Rv, selected pulse energy source agreement, and eight direct Gamma value/
axial-derivative fixtures with controlled infinite-tail remainders.
Interval overlap supplements the defining exact identities; it is not
alone a proof that the two unverified global constants vanish.

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage postpulse
```

The full all-stage pipeline includes this producer/checker after pulse.
Next establish the absolute pressure/angular identities, then obtain higher
axial derivatives and whole-region C4/cone bounds. Full matched-background
acceptance, full physical energy, global admissible stress/flat remainder,
genuine n-dependent recursion, oscillatory cancellation and corrected
Cartesian residual remain incomplete.
