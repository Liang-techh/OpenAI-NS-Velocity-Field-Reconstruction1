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
