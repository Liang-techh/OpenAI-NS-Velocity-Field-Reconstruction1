# Next constructive tasks and acceptance criteria

The next priority is **actual profile data**, not more tests of the toy exponent.
Preserve the existing numerical regression suite while implementing these dependencies.

## P0: Instantiate one admissible leading-profile construction

Translate the constructive choices behind Theorem 4.6 and Appendices A/B/C.
Produce the outer heat profile, inner analytic profile, matching moments,
shear modification and cone constraints. The pointwise moment solver is a
reusable subroutine, not a substitute for that construction.

The current executable parameter-selection chain already contains the actual
outgoing pressure/sigma-side work and the theorem-faithful Lambda/C algebra.
The ST073 exploratory bridge now has separated smooth directions that close
two **physical** tangential outer moments on two sampled scales. A third
radial window and constrained optimization reduce the moment-closure
momentum cost to about `1.43–1.50` times the pre-repair peak, but the
absolute residual still grows under dyadic refinement; see
`ST073_SEPARATED_MOMENT_CONSTRUCTION.md`.
The three-knot extension closes these sampled physical moments also at
`k=15` and reduces off-knot moment leakage by roughly sixtyfold, with
little change in the momentum growth. The next necessary advance is a
genuine residual-canceling update across scales, coupled to the radial
moments and nonaxisymmetric stress realization.
The first exact-curl wave on the moment-closed bridge matches its local
stress covariance but creates a much larger viscous cutoff residual;
one time slope per wave does not cure the spatial/temporal holdouts.
See `ST073_WAVE_ON_MOMENT_CLOSED_BRIDGE.md` before reusing that frozen
wave ansatz. The next wave implementation needs a transported spatial
amplitude, its pressure and mean corrections, and a wider admissible
support or a quantitative cutoff budget.
A first spatial potential/pressure/mean Taylor step on the late bridge
does reduce direct held-out midpoint momentum by about `2.14x` but leaves
an error of `4.86e9`; see `ST073_SPATIAL_WAVE_SLOPE.md`. Continue with
stable multi-stage amplitude transport, wider supported cone geometry,
and coupled radial moment repair rather than accepting the local fit.
Two `1e-10` time steps retain the local factor-of-two held-out momentum
gain but cannot span the pulse window; see
`ST073_TWO_STAGE_WAVE_EVOLUTION.md`. Next establish a stable supported
amplitude inverse or longer-step evolution with uniform residual and
moment control, then test its transfer between actual dyadic scales.
The first actual new-knot trial at `k=21` closes four sampled physical
moments to `8.53e-14`, yet lowers the full momentum peak by only `0.198%`;
the `k=22` holdout still grows with roughly the same `2^1.49` per-step
rate. See `ST073_NEXT_SCALE_TRANSFER.md`. Further moment-only scale knots
are therefore insufficient: couple the wave-amplitude inverse, stress
update, compact mean correction, and moment repair in one residual cycle.
The late-wave explicit coefficient march is locally effective at `dt=1e-9`
but unstable at `dt=1e-8` on direct interior-time momentum; see
`ST073_WAVE_STEP_HORIZON.md`. Implement a supported amplitude inverse
or stable interval solve over an appreciable part of the pulse, with
full residual checks between collocation nodes, before scale transfer.
At `dt=1e-8`, a four-update damped trapezoid coefficient correction
has relative fixed-point defect `7.80` and direct midpoint momentum
`5.93e11`; stronger ridge regularization also leaves the explicit
second-interval error near `3.27e11`. See
`ST073_INTERVAL_CORRECTOR_LIMIT.md`. Stop tuning this physical-time
Taylor fit and implement the pulse-coordinate transverse-amplitude ODE,
normal pressure identity, and support/cutoff controls of Proposition 7.2.
The first actual-source local transverse inverse now solves its frozen
principal ODE, but needs amplitudes `30–120` times the local background
speed; its narrow radial support spans about `247` diffusion times per
pulse half-width. See `ST073_LOCAL_PULSE_INVERSE.md`. Before repeating
the inverse across scales, open an admissible stress-cone support and
establish carrier/cutoff/viscosity balance, then use a moving phase and
complete wave–mean–moment residual cycle.
The same moment-closed bridge's sampled cone band at `k=11,19` cannot
contain the diffusion-balanced width around the existing source patch:
even an optimistic failure-bracket half-width is short by factors `5.91`
and `6.57`. See `ST073_CONE_SUPPORT_SCALE_GAP.md`. Search a different
radial–axial cone region or redesign the moment-matched mean profile
instead of simply widening this wave's cutoff.
A two-scale coarse search found a second cone-positive point at
`(y,eta)=(0.35,0)`. Its refined radial band nearly accommodates the
chosen pulse time scale, but at nearby `y=0.325` the best-centered
sampled axial half-width falls short by factors `5.11` and `4.97` at
`k=11,19`; see `ST073_MIDPLANE_CONE_CANDIDATE.md`. Focus the mean-profile
redesign on widening this axial cone or obtain a joint space–time
carrier/cutoff balance, then check a connected 3D cone on both scales.
The midplane wave trial now uses a pair with **positive actual exact-curl
covariance weights** and matches the center stress at both scales to
`~5.6e-15`; nevertheless its multiplier-`0.1` momentum rises by factors
`5,241` and `5,125` over the background, dominated by the axial-cutoff
viscous term. See `ST073_MIDPLANE_WAVE_TWO_SCALE.md`. Next require a
spatially supported amplitude/pressure solve and mean correction on a
wider axial cone, with complete interior-time residual and cross-scale
checks. Do not treat center stress matching as residual improvement.
The new one-time cross-scale projection does show a reusable candidate:
after `tau^1.5` normalization, the wave-induced defects at `k=11,19`
have `0.999939` cosine similarity. A curl-potential/pressure/mean slope
fitted only at `k=11` and transferred by `tau^-1` lowers the `k=19`
held-out momentum max from `1.44e13` to `5.84e12`. This is not an
evolved correction; see the same report. Next solve the supported
moving-normal amplitude equation throughout a pulse, insert its
derivative, pressure and mean into the full velocity field, and check
nonlinear momentum plus radial moments at interior times on at least
two adjacent dyadic scales. Require absolute residual improvement and
non-growing scale behavior before marking recursion established.
Direct `k=19` interior-time testing now rejects the constant transferred
slope: its corrected/frozen momentum maximum ratios are `52.5, 3.27,
0.404, 1.97, 12.4` at pulse fractions `-0.5, -0.2, 0, +0.2, +0.5`.
Integrating the pulse bump into the slope still gives `34.2, 2.97,
0.406, 1.92, 9.20`. Next integrate the supported potential-amplitude
coefficient ODE from `fit_slope` over a short pulse interval with a
stiff/adaptive solver, recomputing pressure algebraically at each state;
directly check full residual on held-out nodes at interior times before
attempting `k=11→19` transfer. The paper's Proposition 7.2 supplies
the transverse pulse inverse and Section 9 the full correction cycle;
the current collocation ODE is only an exploratory numerical proxy.
An actual two-step `k=19` explicit potential-coefficient march improves
the first held-out midpoint (`0.778` of frozen-wave maximum) but fails
at the second (`1.960`). Four half-size steps give midpoint ratios
`1.108, 0.958, 2.613, 6.788`, with growing coefficient slopes; see
`ST073_MIDPLANE_PULSE_MARCH.md`. Do not extend this explicit Euler
scheme by further step-size tuning alone. Build a supported pulse
inverse with stable time integration, pressure recovery, and a
continuum-in-time residual budget, then couple the mean/stress/moment
operations before claiming recursive contraction.
Implicit midpoint integration in the transferred two-harmonic-plus-mean
template space also fails: even with nine training nodes and solved
three-coordinate midpoint equations, the held-out midpoint maxima
are `12.8` and `315.6` times the frozen wave. See
`ST073_MIDPLANE_IMPLICIT_TEMPLATE.md`. Do not keep tuning the integrator
inside this fixed three-template space. Implement the paper's
spatially varying transverse amplitude along moving pulse paths,
including normal pressure recovery, supported time cutoff errors,
mean flow and moment restoration, before reattempting a dyadic
recursive residual cycle.
The first **spatial** frozen-principal inverse now solves equation
`(7.13)`-inspired paths at five nodes on both `k=11,19`. The
amplitude/background map has only `4.23%` relative L2 drift across
eight halvings, and midpoint complex amplitudes have `0.978%` shape
error after one scale factor. However mode-1 amplitudes reach
`23.4×` and `22.4×` background speed near the axial support edges,
and the endpoint amplitudes remain nonzero; see
`ST073_SPATIAL_PULSE_INVERSE_TWO_SCALE.md`. The next implementation
must move from independent frozen-normal paths to a supported
two-dimensional amplitude with the phase normal transported along
pulse paths, quantify spatial derivatives and endpoint cutoff errors,
then build its exact-curl velocity and normal pressure and measure
the complete momentum on the two scales. The similar normalized
inverse map alone is not a residual contraction.
The five-node inverse was reconstructed into a compact **exact-curl**
field and screened with full nonlinear momentum at `k=11,19`.
Sparse normal-pressure interpolation made the held-out maxima
`15.57×/15.16×` the frozen wave even after physical-coefficient
regularization. Fitting compact pressure gradients to the full
momentum reduced those factors to `3.663×/3.666×`, but the peak is
then viscosity-dominated (`1.82e10/6.73e13`); see
`ST073_INVERSE_CURL_RECONSTRUCTION.md`. The next wave design must
control spatial derivatives of the amplitude and curl/cutoff
remainders on a wider axial cone, not just match pointwise inverse
values or add pressure degrees of freedom. Require direct complete
momentum reduction at both scales before extending the recursion.
An axial mean-geometry change now opens the `eta=+0.05` cone sample at
`y=0.325` on both `k=11,19`, while a 34-variable compensation keeps the
twelve sampled outer moments at `k=11,15,19` below `2.41e-15` in normalized
units. Its cone margin is narrow, the negative side remains closed, and no
full-momentum gain has been shown; see `ST073_AXIAL_CONE_MOMENT_REPAIR.md`.
Next optimize a wider connected cone with moment constraints, then test the
supported correction's complete residual on adjacent scales and pulse times.
The integer-scale holdout found the two-knot cone fails at `k=12..16`.
Adding the same coefficient change at the middle `k=15` knot repairs this
single-point cone at all nine integer scales `k=11..19`, with normalized
moment defect `1.50e-14`; however its weakest cone margin is `0.001819`
and mean-only local momentum still grows `3655×` over the interval.
Do not infer continuous cone support or recursive contraction from this
discrete result. The next gate is a wider connected space-time cone and
absolute full-residual decrease after a supported wave/mean correction.
The widened-support frozen exact-curl wave reduces its midpoint full
momentum peak by about `21%` at both `k=11,19`, but its remaining maxima
are `2.17e9/8.07e12` and grow `3715×` across eight halvings; see
`ST073_WIDER_CONE_WAVE_SCREEN.md`. Continue with transported spatial
amplitude, normal pressure, and mean/stress repair rather than additional
width-only tuning.
The wider-cone compact time-slope/pressure projection lowers held-out
full momentum at pulse center to `18.6%/18.3%` of frozen at `k=11/19`,
and direct exact-curl replay matches that center result. A linear-in-time
realization fails at just `+0.1` pulse half-width: its held-out maxima
are `3.083×/3.074×` the same-time frozen wave. See
`ST073_WIDER_CONE_SLOPE_REPLAY.md`. Use the fitted slope only as an
initial condition or diagnostic for a nonlinear moving-normal pulse
amplitude solve; require direct residual reduction throughout the pulse.
Use those directions as a conditioning prototype, not as the paper's five
leading-profile moments or a recursive correction.
The next Stage-1 inputs are therefore the actual coefficient-family
`remainderBound/remainderLip`, the certified complex compact-set `realPartSup`,
and the resulting coefficient-space fixed point `phi/u/average/pressure`.
Do not reintroduce the superseded Lambda/C-selection task as if that algebra
were still absent.

Deliver a deterministic constructor and a manifest containing every free
parameter, equation reference, cutoff, support radius, quadrature order,
coefficient/table hash and truncation tolerance. Save the actual moment
systems B(eta), Q_eta, d(eta) and account for invertibility and smallness
uniformly on the required eta domain. Establish parity, axis regularity,
pressure balance, matching and support requirements. Do not set
`paper_exact=True` solely because a numerical plot or sampled constraints
look right.

## P1: Solve the background recursion from that profile

Use the landed Eqs. (5.2)-(5.7), Lemma 5.2 compact repair, Eq. (5.15) forward
reconstruction, and finite-prefix recursive cutoff scheduler only after
generating their genuine profile-derived inputs. The remaining work is the
actual converged coefficient hierarchy, its eta-dependent support/stress
closure, true uniform `C[j,m]` bounds, the infinite cutoff/local-finiteness
argument, and arbitrary-order residual decay. Record why each cutoff schedule
is admissible. Compare residuals and derivatives while increasing both
coefficient order and numerical resolution; separate truncation, quadrature
and roundoff errors.

## P2: Instantiate Sections 6-9

The dyadic charts, slow partitions, auxiliary-slot witness, physical
slow-support-to-adjacency bridge, and pointwise phase adapters are already
landed. Next instantiate them with the paper-exact Proposition 5.5 base field
and prove the uniform `LocalBaseBounds/C1-C2` needed for Eqs. (7.9)-(7.11),
then implement the stress cone and amplitudes, curl remainders, compact mean
corrections, and residual-improvement cycle. Each stage needs independently
checked conservation, support, moment and residual identities. A manually
tuned wave-frequency list or finite plot cannot establish the required
iteration and all-order convergence.

## P3: Final compact field, pressure and smooth force

The Section 10 spatial/time localization, closed-past branch, endpoint-jet
adapter and Borel scale arithmetic are already present as formal structure.
The next endpoint task is to materialize the actual closed-past residual's
full spacetime derivative family, prove its locally uniform limits as
t -> 1-, derive genuine analytic compact-template majorants, and only then
complete smooth force gluing through t=1. Check the paper's support, energy
and singular-path requirements. Only then replace the top-level incomplete
status with an evidence-backed completion state. Numerical residuals
supplement the analytic/formal evidence; they do not replace it, and a
force defined from the same residual is not independent verification.

## Upstream verification and handoff

Pin and review the relevant official Lean source definitions and theorem
statements; record source paths, theorem names and hashes. Run the actual
Lean build before claiming a successful formal cross-check. The source pin
currently present is observed commit metadata only.

Every new change must be integrated from the latest `main`, preserving newer
parallel work and the fail-closed provenance gates. Historical PR #5 is an
integration baseline, not a snapshot to restore. Never force-reset or
force-push over another agent's work.
