# Reuse checkpoint: executable legacy bridge

`CompactCandidate.as_legacy_local_field()` connects the new candidate to the
existing `local_field.LocalField` implementation. It supplies an explicit
vector potential; the old `curl_numeric` and `verify.jacobian_numeric` execute
its poloidal curl. The direct velocity is used only to extract the azimuthal
component, not to replace the legacy curl calculation.

Actual verification: `python -m pytest -q -W error tests/test_constrained_legacy_bridge.py tests/test_constrained_candidate.py`
returned 4 passed in 0.26 s. Random points, the axis and the exterior agree
with the direct implementation. `sum_vector_fields` assembles two potential
parts and reproduces the field. This creates an executable route for adding
legacy-compatible potential corrections without losing the curl structure.

The inspected `local_field.py` is composition machinery: its module contract
explicitly says individual corrections are caller-supplied. `background.py`
likewise accepts caller-supplied coefficient profiles and labels its radial
flux identity as kinematic, not a solution of recursive profile equations.
These files therefore cannot be counted as ready-made residual-cancelling
solutions. No broader claim about every legacy module is made here.

Useful reuse now: LocalField, curl_numeric, jacobian_numeric,
sum_vector_fields. Remaining work: supply and fit actual correction potentials
while maintaining nontriviality, support and the configured force restriction.
The scalar legacy route is a composition/reference path, not a fast replacement
for batched optimization. No improved PDE residual is claimed by this bridge.
