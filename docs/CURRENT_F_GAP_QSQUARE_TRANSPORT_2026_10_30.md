# F gap Q-squared ordinary-Z transport — 2026-10-30

Status: **ACTIVE / INCOMPLETE**.  This note records the independent signed
Leibniz derivation for the gap `I1_q_squared` role and one frozen source-only
structural result.  The run returned `status=returned`; its result SHA-256 is
`b516e5a72290062d788114f25237ab5019abc217838a68d125de24caa1323d95`, and its
authenticated source-manifest SHA-256 is
`f25ab1b64c7c236b419c46033713d1ef2d5a5b3ddfe6e334a788d1f1923e424c`.  The
synthetic fixture reported `12` checks and `7` fixtures, all passing.  This is
structural source-role evidence only, not an F source identity, numerical
result, error bound, NS certificate, or global closure.

The reproducible receipt is
`experiments/root_st073/f_gap_qsquare_transport_v1_receipt.json`, SHA-256
`963421b711d19f4c0b5fd26a415be6a9a7138d5ba8cf51bf1e758058f6f92123`.  It
contains the frozen source, runner, source/entry/result payloads, and the two
prior-result raw byte payloads.

## Role and coefficient

The saved F pressure-role leaf identifies the two selected roots as the gap
`I1_q_squared` channel at ordinary Z orders 0 and 1, with formal integrand
`Q^2` and formal m15 coefficient `-1/4`:

```text
F.pressure.m15[Z0] -> gap / I1_q_squared[0]
F.pressure.m15[Z1] -> gap / I1_q_squared[1]
c = -1/4
```

This position comes from the source pressure reduction

```text
S1 = ... - 1/2 * I1_q_squared + ...
m15 = 1/2 * S1 - (upper R*Q_y - lower R*Q_y),
```

so the q-squared integral enters m15 with `c = -1/4`.  The coefficient is a
formal constant, hence `c_Z = 0`, but the selected Z 1 row is not thereby
proved to be the ordinary-Z derivative of the Z 0 integral.

The frozen structural result matched the actual `I1_q_squared` expression and
the same F role correspondence; the prior role packet also matched.  Its
bounded direct graph reads were exactly:

```text
279381 = Q0 * Q0
279475 = Q1 * Q0
279476 = Q0 * Q1
279477 = (Q1 * Q0) + (Q0 * Q1)
```

The two factors were retained as opaque source references `Q0=277042` and
`Q1=277055`.  These refs were not indexed and no source value was selected;
there was no numerical evaluation, source-graph evaluation, or ancestor
replay.  The four direct reads establish the operand structure of the binomial
sum only.

## Independent signed derivation

Let

```text
I0(Z) = integral from l(Z) to u(Z) of M(Z,xi) * Q0(Z,xi)^2 dxi
I1(Z) = integral from l(Z) to u(Z) of
         M(Z,xi) * (Q1(Z,xi)*Q0(Z,xi) + Q0(Z,xi)*Q1(Z,xi)) dxi
F(Z,xi) = M(Z,xi) * Q0(Z,xi)^2
E_Q(Z,xi) = Q1(Z,xi) - partial_Z Q0(Z,xi) at fixed xi.
```

Ordinary Leibniz differentiation gives

```text
D_Z(c*I0) - c*I1
  = c_Z*I0
    + c * integral_l^u M_Z*Q0^2 dxi
    - 2*c * integral_l^u M*Q0*E_Q dxi
    + c*F(Z,u)*u_Z
    - c*F(Z,l)*l_Z.
```

For `c = -1/4`, the complete signed residual is therefore

```text
D_Z(c*I0) - c*I1
  = -1/4 * integral_l^u M_Z*Q0^2 dxi
    +1/2 * integral_l^u M*Q0*E_Q dxi
    -1/4 * F(Z,u)*u_Z
    +1/4 * F(Z,l)*l_Z.
```

The second term is required whenever the stored/source row `Q1` is not
authenticated as the fixed-coordinate derivative `partial_Z Q0|xi`.  It is
not a sign-free magnitude correction.

## Geometry and moving bounds

The gap chart candidate is

```text
M = R^2/mu
R = exp(logRp + xi/mu)
l = 11
u = 13 - 4*mu.
```

At fixed `xi`, with `r = logRp`,

```text
R_Z|xi = R * (r_Z - xi*mu_Z/mu^2)
M_Z|xi = 2*R*R_Z/mu - R^2*mu_Z/mu^2
        = M * (2*r_Z - mu_Z*(2*xi/mu^2 + 1/mu)).
```

The endpoint substitutions are only valid after their source owners and
ordinary-Z rows are authenticated.  Under the separately stated bound-role
conditions that the lower endpoint is the literal source-owned `11` and the
upper endpoint is the source-owned `13-4*mu`, one may set

```text
l_Z = 0
u_Z = -4*mu_Z.
```

Only under those conditions does the residual specialize to

```text
D_Z(c*I0) - c*I1
  = -1/4 * integral_11^(13-4*mu) M_Z*Q0^2 dxi
    +1/2 * integral_11^(13-4*mu) M*Q0*E_Q dxi
    + mu_Z * F(Z,13-4*mu).
```

Until the lower/upper bound roles are source-authenticated, retain the general
`-1/4*F_u*u_Z + 1/4*F_l*l_Z` form.  The endpoint term uses the boundary value
of the integrand; it does not replace the fixed-coordinate `M_Z` term.  If an
endpoint field such as `R*Q_y` is differentiated separately, its own total
chain rule also requires its fixed-`xi` derivative, its `xi` derivative, and
`xi_Z`.

## Implementable expression-DAG shape

The source-bound DAG should retain the following nodes without simplification:

```text
M              = quotient(R^2, mu)
R_Z_fixed      = R*(logRp_Z - xi*mu_Z/mu^2)
M_Z_fixed      = 2*R*R_Z_fixed/mu - R^2*mu_Z/mu^2
Q2_Z_candidate = Q1*Q0 + Q0*Q1
E_Q            = Q1 - fixed_xi_Z(Q0)
upper_term     = -(1/4) * at(F, u) * u_Z
lower_term     = +(1/4) * at(F, l) * l_Z
residual       = -(1/4)*Integral(M_Z_fixed*Q0^2)
                 +(1/2)*Integral(M*Q0*E_Q)
                 +upper_term+lower_term.
```

Required source-role bindings include the F gap `I1_q_squared` rows, the F
gap `Q` ordinary rows, `mu`, `logRp`, the two endpoint coordinates, the F
radius owner, and the gap measure `(R^2/mu) dxi`.  The existing saved role
mapping proves only bounded role reachability and metadata correspondence.  It
does not prove that the Z 1 row equals the complete DAG derivative containing
`M_Z` and the moving-boundary terms.

The frozen direct-read result strengthens this only to the structural statement
that the saved F Z 1 `I1_q_squared` row has the expected two-product shape and
the same F source-role correspondence.  It does not prove

```text
Q1 = partial_Z Q0 at fixed xi
```

and does not identify the opaque Q rows with a common source function.

## Existing source assumptions and open gates

The current source path contains explicit fixed-geometry assumptions:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields.py` builds `I1_q_squared` with the product rule `Q_Z*Q + Q*Q_Z`, but exposes only the value/Z integrand rows.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields.py` marks the radius `radius_Z_independent=True`; its defining Z 1 integral row differentiates the retained integrand row while keeping the radius/Jacobian factor unchanged.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py` declares `endpoint_independent_of_Z=True` and assigns ordinary derivative orders 0 and 1 to the separately supplied value/Z rows.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py` uses fixed-radius endpoint construction, so its existing endpoint C1 rows do not establish an F moving-radius derivative.

For the F rewrite, the following remain unproved and must remain explicit:

```text
F ordinary-Z transport of I1_q_squared       = false
Q1 = partial_Z Q0 at fixed xi                = unresolved
F radius owner and capture-free substitution  = false
mu_Z and logRp_Z source rows                  = unresolved
lower/upper bound ownership and transport    = unresolved
common physical measure/orientation/units    = unresolved
signed error transport                       = false
numeric-use admission                        = false
```

In particular, the successful structural match does not prove Q ordinary-Z
semantics, same-source ownership, the physical measure owner, finite endpoint
ownership, capture-free transport, or signed-error transport.

The formal candidate must not infer `radius_Z = 0` from a caller convention or
from equal opaque node identifiers.  It also must not infer endpoint
independence from a metadata declaration alone.

## Preserved history and non-claims

This local residual does not cancel or rewrite any historical terms.  Preserve
the `212` cumulative, `213` raw-signed, and `183` reduced ledgers; all six
signed channels; the once-only reserved contribution; all four physical pulse
correction roots; the independent reference offset; and the formal difference
`(P-G) + (G-F)`.  Preserve `P1 != m15` and every old historical, measure,
orientation, remainder, and error packet.

This note does not establish F/P/G function equality, m15 closure, a physical
Navier–Stokes solution, a global closure, numeric admission, error reduction,
or removal of any prior error.

## Source-only next step

The next executable leaf should verify the F `Q` ordinary-Z semantics and then
the F measure and finite-endpoint owners.  It should be frozen before entry,
authenticate its source manifest, and be limited to the directly required
source rows and role metadata.  It must not index the opaque Q rows unless that
direct read is explicitly included in the new frozen scope, must not evaluate
source values, and must not rerun an old graph, ancestor producer, or checker.
