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
