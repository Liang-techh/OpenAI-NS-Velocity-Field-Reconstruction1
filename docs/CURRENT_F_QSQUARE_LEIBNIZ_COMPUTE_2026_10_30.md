# F Q-squared Leibniz callback calculator — 2026-10-30

Status: **SYNTHETIC RUN RECORDED / F SOURCE AUTHENTICATION PENDING**. This note records the
independent callback calculator interface and the source conditions that must
hold before it can receive real F rows. It is deliberately a synthetic
arithmetic component: it does not read an F archive, a saved source packet, an
old graph, or a producer module; it does not authenticate a source identity;
and it does not produce a certified bound or remove any prior error.

## Calculator interface

The new independent module is
`src/openai_ns_reconstruction/f_gap_qsquare_leibniz.py`. Its public callable is

```text
qsquare_leibniz_transport(
    M, M_Z, Q0, Q1, dZ_Q0,
    l, u, l_Z, u_Z,
    quadrature,
) -> QSquareLeibnizTerms
```

The compatibility aliases
`compute_qsquare_leibniz_transport` and `compute_qsquare_leibniz` point to the
same callable. The callback contract is:

- `M`, `M_Z`, `Q0`, `Q1`, and `dZ_Q0` are scalar callbacks of `xi`;
- `M_Z` and `dZ_Q0` are caller-supplied fixed-`xi` ordinary-Z derivatives;
- `l`, `u`, `l_Z`, and `u_Z` are finite endpoint data;
- `quadrature(f, a, b)` returns an oriented `d xi` integral and must not add a
  second physical-measure weight, since the measure is already in `M`.

The result dataclass has five finite arithmetic fields:
`term_mz_q0_squared`, `term_m_q0_q1_minus_dz_q0`, `term_upper_endpoint`,
`term_lower_endpoint`, and `total`. These are ordinary floating-point
arithmetic results, not intervals, uncertainty objects, or certificates. The
interface and no-double-weight rule are documented in the module docstring and
callable at `f_gap_qsquare_leibniz.py:1-29,135-172`.

## Four synthetic terms

For caller callbacks, the calculator evaluates exactly these four signed
terms:

```text
T_measure = -1/4 * integral_l^u M_Z * Q0^2 dxi
T_Q       = +1/2 * integral_l^u M * Q0 * (Q1 - dZ_Q0) dxi
T_upper   = -1/4 * (M*Q0^2)(u) * u_Z
T_lower   = +1/4 * (M*Q0^2)(l) * l_Z
total     = T_measure + T_Q + T_upper + T_lower
```

The second term intentionally retains the Q-derivative defect. It must not be
replaced by zero unless a separate source proof establishes
`Q1 = dZ_Q0`. The endpoint order is caller-controlled; the implementation
does not sort or take absolute values (`f_gap_qsquare_leibniz.py:185-234`).

## Actual original Q producer and fixed-delta condition

The original leading gap producer is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py:32-41`:

```text
d = 1 - Z^2
L = 1 - delta*Z^2
N = 2*Z*U - (1-delta)*Z*m - d*m_Z
Q = N/L
```

The producer constructs the gap coordinate, `U`, `m`, `Q`, `Q_y`, and radius at
`:82-87,114-141`. `Q` is formed as a C2 value/Z1/Z2 object, not as a sampled
number. `C2Algebra.mul` and `C2Algebra.div` implement the product and quotient
rows at
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C2_phase_primitives.py:53-71`.

In `radial_Q`, `delta` is passed through `a.fixed(delta)`
(`original_n1_pulse_own_fields.py:34-41`); `C2Algebra.fixed` assigns zero
ordinary-Z rows (`C2_phase_primitives.py:56-58`). Thus the producer-side Q1
formula assumes the source-owned scalar `delta` is Z-independent and that
`L != 0` (in practice the quotient certificate requires the intended positive
denominator). With `delta_Z=0`, the formula represented by the C2 quotient is

```text
Q1 = (N_Z*L - N*L_Z) / L^2,
L_Z = -2*delta*Z,
N_Z = 2*U + 2*Z*U_Z
      - (1-delta)*(m + Z*m_Z)
      + 2*Z*m_Z - (1-Z^2)*m_ZZ.
```

This also requires genuine U/m ordinary jets, a fixed-`xi` interpretation, and
valid differentiation under the defining integrals. The ordinary derivative
engine fails closed unless source parameters are marked Z-independent,
integral endpoints have zero Z derivative, and non-Z function substitutions have
Z-independent values (`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions.py:81-153`).

## F row-wise rewrite is not a derivative-commutation proof

`PulseSourceSubstitution` recursively copies supported operation rows and
replaces bound children at
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py:39-66`.
The F adapter builds its explicit input map from old native rows and rewrites
`selected_pulse`, `leading_fields`, and `pressure` at `:88-102`.

Therefore F's saved leading-gap Q0/Q1 rows are row-wise rewritten source
references. This code does not prove that:

```text
rewrite(Q1_old) = dZ(rewrite(Q0_old)) at fixed xi;
```

The proof would need a capture-free binding closure for every Q dependency and
an explicit proof that substitution commutes with the ordinary-Z derivative.
The separately generated `pulse_coordinates` radius at `:103-113` is only
current-geometry metadata; it is not shown to be the radius used by the
rewritten F field or pressure graph.

The structural Q-squared audit can therefore establish only the operation shape
`Q0*Q0` and `Q1*Q0 + Q0*Q1`. Its formal transport must retain
`Q1 - dZ_Q0`; it cannot infer that defect is zero from row order, metadata, or
opaque node identifiers.

## Geometry and owner gates remain open

The calculator's `M` and `M_Z` are caller inputs. For the proposed F gap
interpretation, the candidate geometry is

```text
R_F = exp(logRp_F + xi/mu_F)
M_F = R_F^2 / mu_F
R_F,Z|xi = R_F*(logRp_F,Z - xi*mu_F,Z/mu_F^2)
M_F,Z|xi = M_F*(2*logRp_F,Z - 2*xi*mu_F,Z/mu_F^2 - mu_F,Z/mu_F)
```

The lower/upper coordinates are candidates `xi_l=11` and
`xi_u=13-4*mu_F`, with a potential moving-boundary term
`xi_u,Z=-4*mu_F,Z`. The actual F radius owner, finite endpoint owners,
`mu_F`, `logRp_F`, their ordinary-Z rows, and the current-radius metadata's
identity with the rewritten field radius are not authenticated. See
`src/openai_ns_reconstruction/pulse_gap_actual_source_program.py:231-250,344-371`,
`src/openai_ns_reconstruction/f_gap_qsquare_transport.py:526-599`, and
`docs/CURRENT_F_GAP_EMITTED_DAG_2026_10_30.md:21-28,53-60`.

The canonical source expression for `mu` being Z-free does not by itself prove
P/F branch equality, common graph namespace, capture-free substitution, or
`mu_F,Z=0`; those remain separate owner gates. No radius, measure, endpoint,
orientation, unit, or Pstar identity is inferred here.

## Preserved ledger and non-claims

This calculator does not include, cancel, reinterpret, or reduce any inherited
source/error ledger. Preserve all of the following:

- `212` cumulative terms, `213` raw signed terms, and `183` reduced terms;
- all six signed channels;
- the once-only reserved contribution;
- all four physical pulse corrections;
- the independent reference offset;
- the formal `(P-G)+(G-F)` separation and every existing error packet;
- `P1` and `m15` as distinct objects.

The calculator is not a real F source evaluator, not a source-identity proof,
not a certified error transport, and not m15 closure.

## Execution and receipt — synthetic run only

```text
status: FROZEN SYNTHETIC RUN RETURNED; F source and geometry owners remain unproved
source_sha256: d0605e8246040f71e37730b1488dc268d7388481597426a2b1fd10a48640afdd
result_sha256: 95572731144c1d8a1514ece11e074dd34bb292459216aa2951452d452e141dcc
function_calls: 3
rational_table_checks: 15/15 passed
forward: (-15.75, 76.5, -126, 9, total -56.25)
reverse: (15.75, -76.5, -9, 126, total 56.25)
coincident_endpoints: (0, 0, 9, -4.5, total 4.5)
comparison: finite floating-point comparison only; not a certified error
receipt: experiments/root_st073/f_gap_qsquare_leibniz_v1_receipt.json
receipt_sha256: 4a8e28c54c97d5ef6078e8c0ca2985c15a5d184bf254dceb36d5ac227f2007cc
F source values: NOT SELECTED
F identity: NOT PROVED
geometry owner: NOT PROVED
errors removed or reduced: FALSE
```

The returned run is limited to the independent synthetic callback and its
hand-computed rational table. It is not a real F calculation, source-identity
proof, certified error bound, or error reduction. The F Q and geometry owners
must still be bound to the same frozen source graph by an independent
source-only proof before any real F rows can be connected.

## Next step

Freeze a narrowly scoped source manifest that binds F Q0/Q1, fixed-`xi` Q
derivatives, `mu/logRp`, radius, measure, and both endpoint owners to one source
namespace. Independently prove the row-wise rewrite/ordinary-Z commutation and
the moving-boundary chain rule. Until those gates pass, do not connect real F
rows or real F geometry to `qsquare_leibniz_transport`, and do not claim a
Leibniz residual, error reduction, or m15 closure.
