# ST073 implicit axis jets and wider axial core

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
uses a low-residual inner core, a radially matched exterior (Appendix A),
and oscillatory momentum flux to cancel the annular residual (Sections 2
and 7–9). Our ST073 inner recurrence is an autonomous local prototype,
not the paper's full profile or oscillatory correction.

The previous ST073 axis-jet implementation used a Lagrange series in an
auxiliary index. Its numerical implementation was restricted to
`|eta|<=0.5`. `extended_axis_jet.py` computes the same bivariate Taylor
coefficients directly from the implicit source-coordinate relation

`Q - (eta+zeta)^2 Q^(2h) = 1-eta^2+theta`, `Q(0,0)=1`.

It solves the coefficients of `Q` degree by degree, then uses
`Q*d(Q^beta) = beta*Q^beta*dQ` for the five axis powers in the original
recurrence. This changes the numerical axis-jet evaluator, not the
recurrence equations or the selected ST073 axis traces. Within the old
strip, tested velocities agree with the old evaluator to floating-point
precision at `X=.01,.03,.06`; independent fourth-order physical-coordinate
differences at `eta=.7,.8` gave momentum residuals around `1e-8` to
`2e-7` and divergence around `1e-9` or smaller at the two tested steps.
These are sample checks, not a certified continuum bound. The new local
analytic residual at `X=.06` stayed below `1.4e-7` on the sampled
`eta=0,.3,.5,.65,.75,.85` nodes at `k=5.5`; see
`experiments/root_st073/extended_axis_screen.json`.

`extended_compact_join.py` selects `eta_max=.97`, a flat inner segment
through `eta=.7`, and an axial cutoff ending at `eta=.94`. It uses the
same old radial bridge and pressure fit, then the inward heat-swirl
attachment. At `k=5.5`, `X=.03`, the old field's residual norms at
`eta=.3,.48` were about `18342,121`; the extended candidate gave about
`6.4e-7,1.9e-6`. Its new collar still gave about `3144,750,372` at
`eta=.75,.85,.9`. These values are in
`experiments/root_st073/extended_compact_join.json`.

At the same time and `X=.5859375`, the selected field's radial-bridge
residuals were about `2953,6600,4982` at `eta=.3,.65,.8`; see
`extended_radial_bridge_screen.json`. The pure heat exterior starts at
`r=0.04487` at the earliest registered time. With heat age `2048`, its
analytic-form spatial/time residual supremum is about `4.62e-4`, based on
numerically integrated radial pressure; the axial Gaussian energy factor
is about `11.34`. The finite-energy and divergence-free construction is
still **not** globally momentum-accepted; the radial bridge and axial
collar exceed `1e-3` by orders of magnitude, and no whole-volume L2 gate
has passed.

Moving the cutoff closer to `eta=1` reduces selected inner-collar samples
further, but enlarges the physical support and leaves radial-bridge
matching unresolved. Those endpoint variants are exploratory only.
Following the paper's actual mechanism, the next work should derive the
radial moment/stress mismatch for this extended background and realize a
dynamical, solenoidal correction in the bridge rather than regard the
wider low-residual core as a complete solution.
