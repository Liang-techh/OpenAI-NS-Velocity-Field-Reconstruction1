# Current global exterior gap

The current enriched candidate is a locally registered field. It does not
yet define a finite-energy velocity on all of R3, regardless of the
improvements in the local patch momentum residual.

`broad_meridional_constrained.load_saved_field()` inherits the dynamic
controlled mean built through `outer_feedback_evolution.build_current()`.
Its `JoinedField` / `GroupedJoinedField` base only accepts the registered
axial slab abs(eta) <= 0.5 and the registered time range down to k=20.
At fixed tau, arbitrary exterior z values are therefore not defined by the
current field API.

Within the slab, the radial exterior calls
`src/openai_ns_reconstruction/heat_exterior.py:physical`. That evaluator is
independent of z and explicitly has no global axial localization. Extending
that same nonzero radial heat field unchanged through all z would give an
infinite 3D energy integral, even though its radial tail decays. This is a
statement about that natural extension; the registered current API itself
does not define a global energy integral.

The compact Fourier wave and tangent corrections are supported in both
radial and axial coordinates. They cannot repair the non-global base.
The radial bridge pressure also does not supply an axial exterior pressure
closure.

The missing construction is an axially localized, globally defined
divergence-free exterior with an explicit pressure extension and matching
to the local field. A streamfunction/vector-potential cutoff can be a
useful intermediate candidate if it preserves the present patch exactly,
but all new cutoff/transition momentum terms must then be included and
corrected. Multiplying the meridional velocity by an axial cutoff would
generally destroy divergence freedom and is not an acceptable substitute.

The immediate code investigation is locating an analytic streamfunction
or vector-potential representation through the current mean wrappers.
No global closure or finite-energy certificate is claimed by this audit.

## Paper alignment and terminal-time limitation

Reference: OpenAI, Section 10.1, equations (10.1)--(10.5),
https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

The paper localizes a vector potential before taking its curl, with a
separate axisymmetric swirl term and explicit pressure cutoff. Its spatial
cutoff is fixed in physical coordinates. Its smooth force extension across
the terminal time additionally uses local residual flatness and compatible
limits of derivatives (Section 10.2).

Our registered slab only reaches k=20. A cutoff tied to eta can supply a
spatially global construction on that finite time interval, but its physical
support shrinks with tau. It does not establish agreement on a fixed
neighborhood through the critical time, or the paper's terminal force
regularity. Those are separate missing requirements after implementing the
finite-interval extension.

## Implemented finite-interval prototype

`global_axial_extension.py` reconstructs the current mean's analytic Stokes
streamfunction through the wrapper chain, applies axial and radial cutoffs
to it, and adds the compact curl wave corrections. The pressure extension is
explicit. Unknown wrappers raise rather than silently discarding components.
The implementation introduces no numerical radial primitive at field calls.

At the reference time the full wave support has eta range
[-0.278382, 0.255987], contained in the chosen plateau |eta| <= 0.30.
The axial cutoff ends at |eta|=0.49. The mean radial support ends at
r=0.01013845; the compact wave is inside that support at this time.
Eleven plateau probes agree exactly in velocity and pressure. Off-support
probes return zero. This is a spatially localized prototype, not a completed
critical-time construction.

Finite-difference divergence decreases by about four when halving the
central-difference spacing; the finest full-field sampled error is 8.068.
This supports the expected truncation trend but does not meet a 1e-3
numerical divergence tolerance. Analytic incompressibility depends on the
streamfunction adapters matching all meridional velocity contributions.
Additional checks at nonzero time offset will exercise the slope wrappers.

The sampled full momentum maximum in the new cutoff collars is
1.564028779e9, including pressure and cutoff derivatives. This is not a
spatial supremum or L2 integral, and no compensating force is declared.
The new transition requires dynamical correction before acceptance.

## Whole-support kinetic energy

`global_support_energy.py` builds a bounding cylinder from the union of the
moving mean support and the fixed compact wave support. This union matters
away from the reference time, when the wave may extend outside the shrinking
mean support. It then integrates one-half of |u|^2 with cylindrical volume
weights over the entire bounding cylinder at the reference time.

The order-6 panel rule (5,760 points) gives 0.0250475853; order 10
(16,000 points) gives 0.0250030141, a relative difference of about 0.1783%.
The finite positive sampled energy confirms that the constructed field is
nonzero; compact spatial support and bounded fixed-time ingredients provide
the finite-energy rationale. These two numerical rules do not certify an
integration error or a uniform bound at the critical time. Full-domain
momentum and time evolution remain unresolved.

The streamfunction reconstruction has now also been checked at k0+1e-3,
where slope contributions are nonzero, using five-point spatial derivatives.
Maximum sampled meridional-velocity discrepancy across reference and shifted
times is 3.66e-9. The actual geometric endpoint k0+1e-6 separately preserves
the whole wave support in the plateau and gives zero sampled velocity and
pressure differences. A generic BroadShearSlope adapter was corrected to
include its time multiplier; its concrete zero meridional base makes this
fix numerically neutral for the present candidate.
