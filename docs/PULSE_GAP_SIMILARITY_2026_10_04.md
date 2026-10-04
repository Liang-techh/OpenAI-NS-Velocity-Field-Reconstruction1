# Whole original inactive-gap similarity companion — 2026-10-04

The entire original inactive gap now has a same-source companion for all five
raw cumulative moments through mixed4, velocity through mixed4, full
meridional similarity stress through mixed3, and the same absolute pressure
through mixed4. It covers xi in [11,13-4mu], with d=13-xi in [4mu,2].
The old main gap xi[11,12] and end gap s[-1/mu,-4] are two charts of this
single source; ordinary differentiation is d_s=d_logR=-mu*d_d.

The axial input is exactly zero here. Nonzero cumulative linear moments,
radial velocity, mixed angular/axial history and selected energy loss remain.
No primitive is reset at a support edge.

## Uncapped sources and stable factors

Use the same selected C5 control functions Cj(Z), original normalized beta
and complete future energy F(Z). Full beta weights yield

`Mi=-sum_j Cj*exp(lambda_i*center_j)*Wi`,
`lambda_i=.5-i*mu`.

The actual normalized moment factor has the source logarithm

`logDi=(d/2-1)/mu+finite-i*d`,

where `finite=-3*L+2*log(u0)-log(6)/2-2*log(mu)`. The producer replays
the original logscale and both native gap recipes to establish
`logD0=-1/mu+finite`. Backward growing exponentials are not materialized,
and inverse-mu terms are grouped before enclosure. Direct replay binds both charts' leading factors and the actual xi=12/s=-1/mu switch. Full beta weights use the original fixed 256 integration cells.

The exact original histories are

`e0=F*exp(-2d)-expm1(-2d)/(4mu)`,
`J=J0*exp(-2d)`,
`X=1/(1-mu)+(Xp-1/(1-mu))*exp(-(1-mu)*(13-d)/mu)`.

The stable expm1 expression preserves the radial energy integral at the
right boundary d=4mu even when mu is far below working precision.

Absolute pressure remains

`P=-C^2/(2p)+exp(-p*d/mu)*(P0+C^2/(2p))`, `p=1+2mu`.

P0 is the exact current canonical flatten pressure getter in the original
pulse units. Its whole-end source enclosure also bounds its terminal value
at s=0; that enclosure is used only as a bound for the same function.
The existing analytic datum/FTC source identity is directly consumed.
The original cumulative pressure moment Mp is separately retained with its
same Pin, U and forward time decay. Thus absolute pressure and the fifth
cumulative moment are distinct retained quantities. The Mp inlet, swirl limit and decaying tail are separate sectors; the source U^2 is retained explicitly as 2*actual_log_inlet_U.

All raw linear moments and the signed angular-memory contribution are
constant on the inactive gap after exact source-factor cancellations.
The raw energy moment keeps the complete future, actual radial integral and
selected end loss separately. Their factors are never reconstructed by
subtracting enormous rounded numeric logs. The production logRp is evaluated
directly from its original logC/logP/Tw recipe.

## Source joins and evidence

Actual gap_shapes and the original full pulse_coefficients source are replayed
with arbitrary axial functions. All ordinary stress rows0..3 agree
identically, including full radial transport and signed histories. The
right boundary d=4mu is the same end s=-4 source through
velocity4/moment4/stress3/pressure4. Original main/gap and chart-switch
functional source certificates remain, using uncapped defining sources.

Focused checker PASS: 423 current hashes, 45 source-function identities and
945 finite signed coefficient rows across the whole gap, right boundary
and coordinate switch. It retains 75 structural zero pressure/velocity rows.
Three independent fixtures use direct complete original-beta integrals and
the original full stress evaluator. They supply489 comparisons:
stress60, velocity135, pressure45, constant raw moments24 and all five raw
moments225. All three radial histories are nonzero. Tolerance1e-55,
maximum positive enclosure miss3.3585e-71. Fixtures check formulas and units;
they do not certify a cone or the final corrected NS residual.

Run only the bounded new stage:

`python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulsegapsimilarity`

Next recover the full physical tensor, diagonal/divergence and actual
three-component remainder on this unchanged whole gap; then prove its
continuous admissible cone and compose the accepted end/tail. Main pulse,
exit, entrance and upstream finite-width bridge feedback remain. Completed
global tensor admissibility, global temporal flatness/volume/energy,
n-dependent recursion, oscillatory correction and corrected Cartesian
residual/dynamics remain unfinished.
