# Actual micro six-I1 C1 defect and pressure-error bounds

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The six
same-source finite-N I1 integrand defects now have ordinary Z0/Z1 norm caps
on the original micro charts for Z in [37/100,93/250]. Their RÃƒâ€šÃ‚Â²*hb-weighted
joint integrals, reduced volume D, and half-S1 pressure error excluding the
outer RQy term are certified below 2^-1024. The complete micro pressure
error is not certified small; its outer RQy majorant is retained.

## Actual background and increment units

The saved native_U and native M derivative intervals already enclose the
physical axial velocity U and physical normalized radial mean. C6 stores
V=U/Pstar and m=M/Pstar; Q4 restores physical U=Pstar*V and m=Pstar*m.
Therefore this stage uses the saved U/M intervals directly, with no extra
Pstar multiplication. M is the normalized radial mean, not R*M. All its
source inlets and histories remain present. Native and comparison histories
are not identified.

The physical background Q0..2 norm provider is new arithmetic on the saved
native U0..2 and M0..3 intervals using the same linear Q operator. These are
norm bounds, not selected signed U/M/Q values. The actual first and second
loops use the same native bridge.micro source, verified by its source AST.
The previous mean/ODE identification makes this background exactly the
U0=Ufull-delta_U and Q0=Qfull-delta_Q functions in the accepted defect graph.

## Six ordinary C1 differences

For full=background+increment, the source differences are

```
delta(q)       = delta_Q
delta(q_Z)     = delta_Q_Z
delta(qÃƒâ€šÃ‚Â²)      = 2Q*delta_Q + delta_QÃƒâ€šÃ‚Â²
delta(U*q_Z)   = U*delta_Q_Z + delta_U*Q_Z + delta_U*delta_Q_Z
delta(UÃƒâ€šÃ‚Â²)      = 2U*delta_U + delta_UÃƒâ€šÃ‚Â²
delta(U*U_Z)   = U*delta_U_Z + delta_U*U_Z + delta_U*delta_U_Z.
```

The equations and their ordinary Z1 derivatives are formed before taking
norms. Square terms and mixed increment terms are retained. Product rules
use ordinary binomial coefficients. The C1 products consume the available
U/Q/increment rows through order2; mean Z3 is retained in Q recovery.
No missing top derivative is filled by zero.

## Exact physical weights and joined endpoints

Write h=1-delta, d=1-ZÃƒâ€šÃ‚Â², L=1-delta*ZÃƒâ€šÃ‚Â². The actual physical I1 measure is
RÃƒâ€šÃ‚Â²*hb*ds, once, on [sc/2,1] and [1,2]. The saved upper graph handle 48 is
the exact number 2; it is not coordinate 48. There is no additional Pstar in
this measure. Each native width<=1 gives a local cap hb*R_upperÃƒâ€šÃ‚Â² times the
source sup norm, and the joint cap sums both local bounds. R, hb and the
chart endpoints are Z-independent, so ordinary Z differentiation introduces
no omitted measure terms.

The reduced pressure coefficients are retained exactly:

```
B = endpoint_difference[RÃƒâ€šÃ‚Â²*q]/L
    + endpoint_difference[RÃƒâ€šÃ‚Â²*qÃƒâ€šÃ‚Â²]/2
    - 2ZÃƒâ€šÃ‚Â²*endpoint_difference[RÃƒâ€šÃ‚Â²*UÃƒâ€šÃ‚Â²]/LÃƒâ€šÃ‚Â²
D = -I[q]/L + hZ*I[q_Z]/(2L) - I[qÃƒâ€šÃ‚Â²]/2 + d*I[U*q_Z]/L
    + 2hZÃƒâ€šÃ‚Â²*I[UÃƒâ€šÃ‚Â²]/LÃƒâ€šÃ‚Â² + 2Zd*I[U*U_Z]/LÃƒâ€šÃ‚Â²
micro_pressure_increment = (B+D)/2 - endpoint_difference[R*Q_y].
```

There is no extra Z on the U*q_Z coefficient. The factors 2 on the last two
volume coefficients and the boundary UÃƒâ€šÃ‚Â² coefficient are source-backed;
ordinary Z1 coefficient rows are differentiated as signed interval functions
before conversion to log norms.

The accepted free-Z fourteen-field seam theorem cancels the internal B and
RQy traces as functions. Since the background is one common native source,
the corresponding finite-N error traces cancel too. That theorem is
consumed as a receipt; its fingerprint checker is not rerun. The first
lower endpoint lies in the original source-owned flat collar: the finite-N
loops and zero-inlet history give exact zero increment there. This source
identity justifies a zero lower error. It is not inferred from a small cap.
The outer RQy term is never discarded.

## Directed error results and limits

Write Csel for selected positive logC. Every six-defect C1 log norm is below
-2600*Csel; the joint I1, reduced D and half-S1 error log norms are below
-2700*Csel. The local 2^-1024 threshold applies to these correction errors,
not to the leading-background or complete corrected NS residual.

The following are rounded displays of directed log caps divided by Csel;
exact tuples in the report remain authoritative.

| error expression | Z0 log-cap/Csel | Z1 log-cap/Csel |
|---|---:|---:|
| reduced_volume_D | -4769.62334199342517903509 | -2749.62334199342517903509 |
| outer_reduced_boundary_B | -5977.62334199342517903509 | -4968.62334199342517903509 |
| half_S1_pressure_without_RQy | -4769.62334199342517903509 | -2749.62334199342517903509 |
| retained_outer_RQy | 254.0 | 2075.0 |
| complete_micro_pressure_error | 254.0 | 2075.0 |

The complete pressure majorant remains positive in these units because of
the retained outer RQy bound. This does not establish a positive lower bound
on actual error. No cap, midpoint or overlap is used as a signed integral or
function value. The six signed integrals remain exact source functions and
have not been numerically evaluated here. All old 212 forward terms, six
signed inner/lower contributions and prior error packets remain unchanged.

Global five-moment closure, reserved-pulse terminal functional matching,
exact heat/energy, stress/flat remainder, true n-dependent temporal
recursion/summation, realizable stress corrections and independent corrected
Cartesian NS validation remain incomplete.

## Execution and next work

New production exits 0 in 2.966 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_I1_defect_bounds.json` is
18,006,425 bytes, SHA256 `88df21cdcf020204fa0f651dae58ed0ffb98cbc1bfe9adb360c990d770ab5492`, with 2005 input and
375 implementation pins. Saved-only inspection,
source-scale and joined-endpoint binding, exact C1 defect controls,
directed pressure-error checks and independent bounded Luna/max math review
pass. A pre-acceptance writer gate found an inconsistent nested-chart
schema check; the runner type check was corrected without changing the
mathematical module or enclosures. No report was accepted before that fix.
A later format-only terminal-newline normalization refreshed the runner
source pin; the report records the unchanged mathematical-packet digest and
no accepted calculation was replayed. Separate runtime model metadata is
unavailable. No accepted producer,
owner, norm, recurrence, quadrature or seam checker was replayed; the new
background Q bound uses new inputs to the pure arithmetic operator.

Implementation: [I1 defect provider](../src/openai_ns_reconstruction/original_micro_I1_defect_bounds.py),
[runner](../scripts/run_original_micro_I1_defect_bounds.py).

```
python scripts/run_original_micro_I1_defect_bounds.py --inspect-report
```

NEXT resolve the external RQy contribution using its original correlated
phase/source endpoint or an exact adjacent-source identity. Combine the
certified I1 error with actual signed background integrals, preserving the
common mean, pressure and incoming histories, then advance the terminal
reserved-pulse functional comparison. Goal ACTIVE.
