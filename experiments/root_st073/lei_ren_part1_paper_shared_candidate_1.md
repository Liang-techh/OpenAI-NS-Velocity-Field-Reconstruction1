# Shared candidate 1: actual core and short connection

Run `python experiments/root_st073/lei_ren_part1_paper_shared_candidate_1.py`.
The script recomputes the shared pressure anchor, degree 18 core, comparison,
actual exit moments and Z tangents, and short switches at 260 decimal digits.
Parameters: j=1e-14, Lambda=1e36, logCstar=5e151, logPstar=14,
delta=1e-200. The schedule uses 233 Decimal digits to preserve finite offsets
in its enormous log reference radius. The five necessary input gates pass;
these are lower-bound rejection tests, not sufficient source conditions.

Six core probes at the actual axis root and Z=.3, Lambda R=1,4,4.1,
have F>0 and F_R<0. The real-axis norm upper bound is approximately
1.348e108. This does not cover mixed core derivatives or the complex domain.
The budgets A=1e150 and logK_upper=1e152 remain provisional. The common
positive hb=epsilon is defined by loghb=-100-100*1e152. Its compliance
with the full source constants is unproved.

At Z=.3, stage 1 midpoint, stage 2 midpoint/end and R=110 pass the sampled
relaxed stress test. Stage 2 midpoint has a=.4,b=0; the endpoint and R=110
have a=.8,b=0. The stricter admissible test does not pass at these probes.
Explicit switch phase coordinates retain the distinct shears although all
three switch probes round to physical R=100. No physical mesh resolution
or uniform cone estimate is claimed. The inherited frozen continuation
also carries conditional truncation bounds, not a global ODE certificate.

Next work, in dependency order:

- [ ] Bound the mixed (R/Ra,Z) core C3 quantities, reciprocal F, pressure
  derivatives, moment/driver quantities and inverse D defining A and K.
  Retain signed/log quantities and distinguish bounds from finite samples.
- [ ] Bound A_Omega on the source complex neighborhood; do not replace it
  by a real-axis maximum. Check Cstar >= Lambda^2 exp(Lambda A_Omega).
- [ ] If either provisional budget fails, rebuild all components with the
  new parameters using this script rather than rescaling old receipts.
- [x] Implement the numerical long angular reshaping interval from R=110, transporting
  actual five moments and their Z jets. Integrate in shifted log radius,
  preserving offsets; match both values and shears at the short endpoint.
  See long_reshape.py/md/json. Conditional tail bounds and finite quadrature
  are implemented; global source constants and interval certification remain open.
- [x] Restore the numerical axial reference profile through Rh. See
  axial_restore.py/md/json for actual moment transport and defect diagnostics.
- [x] Implement the numerical representative five-bump correction, carrying
  pressure and Ur from the same fields and continuity formula. See
  inner_corrected_field.py/md/json; unresolved entries retain explicit intervals.
- [ ] Establish actual uniform moment closure with source norm and quadrature
  enclosures before claiming the complete source moment repair.
- [ ] Match that repaired background to the heat exterior and establish
  finite energy before claiming a global velocity field.
- [ ] Replay multiple times in the physical chart, measuring core radii,
  axial/radial ratio, velocity growth and interscale matching. A sequence
  of separately scaled local profiles is not established scale recursion.

Full stress/remainder closure and oscillatory correction remain unfinished.
The full momentum residual target is not evaluated by this receipt.
