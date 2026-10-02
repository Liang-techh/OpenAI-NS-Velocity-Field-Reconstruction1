# Original three-chart prescribed-shear bridge mixed4

The actual Ra..100 prescribed-shear bridge now has directed physical
velocity/pressure and five-primitive mixed derivatives through total four.
Both original microscopic smoothing charts retain varying comparison fields.
The following macro chart uses the frozen comparison fields with their
still-varying own moments. The actual core exit, R100 histories, original
analytic P0 and signed bridge source remain unchanged.

These are enclosures of the same analytic construction, with separate exact
function references. They are not newly reconstructed point coefficients,
dynamical measurements, a complete Cartesian field or a stress lift.

## Three charts and actual source

Write y=log(R/Ra), Ra=4/Lambda and hb=cstar*K^-100>0. The original controls
are alpha=1-sigma((y-hb)/hb) and chi=1-(1-hb)*sigma(y/hb).

| Chart | Coverage | Comparison alpha | Actual chi | Derivative coordinate |
| --- | --- | --- | --- | --- |
| first | y/hb in[0,1] | 1 | 1-(1-hb)*sigma(s) | s=y/hb |
| second | y/hb in[1,2] | 1-sigma(s-1) | hb | s=y/hb |
| macro | y in[2hb,log(100/Ra)] | 0 | hb | y=logR |

The macro fraction only selects coverage. It is not differentiated as y.
All radii stay formal: R=Ra*exp(hb*s) in smoothing, and
y=2hb+fraction*(log(100/Ra)-2hb) thereafter. Numerical width, radius and
theta=Ra/R intervals enclose these sources; they do not define them.

The actual angular field obeys logF_y=-chi*Dbar/2. The actual raw axial
field obeys

```text
V_y = -chi*(phi_actual/barphi)
      *(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl).
```

Full cutoff, quotient, comparison, moment, current-radius and three-scale
Leibniz derivatives are retained. In particular, barphi is radially varying
in the two smoothing charts, so its derivatives must enter the quotient.
The old frozen-comparison recurrence cannot be used there.

## Same analytic core: rectangular orders

Comparison direction axial5 requires one extra Z derivative. Its radial3
also requires core radial3, so the needed core rectangle is rho0..3/Z0..6.
The old actual core packet supplies total order<=5. The ten missing order
pairs are (0,6), (1,5..6), (2,4..6), (3,3..6).

For those pairs the new wrapper uses the admitted FULL analytic norms:

```text
Phi_norm <= Phi_model_Xh_norm_upper + correction_norm
Psi_norm <= fresh_Psi_model_Xh_norm_upper + correction_norm
|D_rho^i D_Z^k Phi| <= Phi_norm*E(i,k)
|D_rho^i D_Z^k V|   <= epsilon_core*Psi_norm*E(i,k)
E(i,k)=(i+k)!/[20^i*h_X^k*(k+1)^2*(1-rho_max/20)^(i+k+1)].
```

The missing V rows have zero derivatives of the baseline4Z+j. The extra
epsilon_core factor is therefore essential. h_X is the analytic norm
parameter, distinct from the microscopic width hb. These are full-norm
bounds containing both model and nonlinear correction, not correction-only
bounds or arbitrary new coefficient boxes. Existing sharper total5 rows
are reused unchanged.

The Phi model norm comes from the full Bessel coefficient/Cauchy majorant.
The Psi model norm comes from the fresh compliant transfer proof with its
pressure contribution. The older linear record's uncertified Psi-model
flag is not used as admission. Transfer, tube, linear and actual core
source hashes are all bound; the new checker verifies every requested
high-order embedding row rather than treating old C3 rows as high-order
evidence.

Euler derivatives use (rho*D_rho)^k, with Stirling coefficients, not powers
of rho alone. Core logPhi derivatives come from the positive core Phi
quotient series. The original continuation remains inside rho<=4.1.

## Exact comparison functions and moment enclosures

Separate formal functions retain

```text
logbarphi(y)=logPhi_core(4,Z)
            +integral_0^y alpha(t)*D_t logPhi_core(4exp(t),Z)dt
barV(y)=V_core(4,Z)+integral_0^y alpha(t)*D_t V_core(4exp(t),Z)dt.
```

The six own comparison moment shapes retain their actual core integral
initial values, positive integrating-factor kernels and the ODEs

```text
H_y+2H=2barphi        M_y+M=barV
K_y+2K=2barphi*barV   A_y+A=barV^2
B_y+2B=barphi^2       C_y+C=barphi^2.
```

Their numerical values are positive-kernel enclosures from the SAME smooth
histories. They are not exact moment point values, unsmoothed frozen
histories, or substitutes for the actual field's five moments. One common
comparison source namespace and argument binds these functions across
all charts. Alpha's IBP weights preserve pure axial orders through6 without
an inverse-width factor; mixed radial orders keep the original source
width factors explicitly.

## Physical derivatives and local joins

The accepted factored algebra carries hb, Pstar, F0 and physical current
amplitude factors through Bell, quotient, quadratic and axial derivatives.
Only final ordinary physical coefficients are numerically capped. Phase
derivatives are converted by hb^-k, with source powers cancelled BEFORE
forming the logR triangle bound. No inverse width or enormous derivative
magnitude is materialized.

Ur differentiates sqrt(R/2) before forming the grid. All five primitives
use their physical RHSs with fixed current-basepoint normalizations; raw V
is not divided by Pstar. The axial-quadratic term retains Pstar^-2, and
pressure retains the original P0.

At the core exit, alpha=chi=1 with flat positive-order cutoff jets. The
comparison agrees with the same stress-free analytic core on y<=hb, and
the actual controls therefore agree with the core ODEs at y=0 to all
required orders. Actual core primitive initial values/P0 are identical.
This argument relies on the accepted analytic fixed-point construction,
not on comparing independent interval representatives.

At y=hb and2hb, original flat alpha/chi jets and common exact comparison
and actual histories give the functional joins. At R100, chi=hb and the
same frozen comparison and actual histories feed the first switch;
D_s=hb*D_y identifies its derivatives. Exact directed phase1 physical
grids and identical actual smoothing-exit histories provide additional
consistency checks. Original R100 physical zeroth-radial rows/P0 agree
with the accepted microswitch inlet.

## Evidence and reproduction

- 80 requested rectangular core extension rows checked against full source
  norms and the general mixed embedding formula.
- 480 microscopic velocity/pressure and600 five-primitive bounds;
  240 macro velocity/pressure and300 primitive bounds.
- 1620 final factored source sums and logR ledgers,12435 positive-source
  cap inequalities, and current family/source/hash prerequisites.
- 42 independent Euler/log/core derivatives,96 varying-comparison
  direction derivatives from closed physical moment integrals, and192
  original alpha/chi/quotient/three-drive derivatives using actual finite
  partial integrals and FTC. Finite fixtures validate formulas, not actual
  parameter admission.
- 63 symbolic Euler, embedding, IBP, kernel, cutoff and functional source
  identities; exact actual core/R100 histories and original pressure datum.

Run `experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py
--stage bridgemixed`. The ordered pipeline now has102 modules. The receipt
is `experiments/root_st073/lei_ren_part1_paper_compliant_bridge_mixed_C4_check.json`.

The next dependency is complete source dispatch and the Rh-to-pre-O3
continuation/outer source binding, then whole physical Cartesian assembly.
Required-domain and terminal energy, admissible divergence-form stress,
independent flat remainder, genuine n-dependent recursion and oscillatory
correction remain unfinished. Full-field/stress/global-energy/temporal
gates remain false; original unlocalized whole-space energy remains infinite.
