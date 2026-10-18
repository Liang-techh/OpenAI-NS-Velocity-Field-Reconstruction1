# Actual micro-to-macro RQy endpoint identity

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The new
proof connects the original second_micro upper endpoint s=2 to the original
frozen_macro lower endpoint m=0. Their physical R*Q_y traces agree as
functions of Z through ordinary derivative orders 0, 1, 2 and 3.
This is an exact source identity composed with accepted source-seam
receipts. No phase, cap or interval midpoint is selected as a function value.

## Exact endpoint and scale conventions

The second_micro upper RQy handles are 133745..133748. The macro lower
trace is the substitution m=0 in finite_boundary handles 141027..141030.
Handles 141114..141117 describe the macro upper endpoint and are not used
as its lower endpoint. Endpoint substitutions are Z-independent.

The native coordinates and Jacobians remain distinct: micro uses s with
dy=hb*ds; macro uses m with dy=length*dm. Physical y differentiation is
used on both sides. The radius recipes are

```
R_micro(s) = Ra*exp(hb*s)
R_macro(m) = Ra*exp(2*hb + length*m).
```

The accepted source theorem supplies the common physical Q_y jets,
original global phase and single N, and complete incoming mean and
pressure histories. The new calculation verifies the actual Q4 boundary
wrappers and the endpoint radius substitution. It does not replay that
source theorem or replace either boundary output by an arbitrary value.
Physical Pstar scaling is retained once in Q4; native Jacobians are never
identified with each other.

## Scope of the cancellation

For adjacent signed pressure increments, the boundary contribution is

```
second micro: -RQy(s=2) + RQy(s=1)
macro: -RQy(m=1) + RQy(m=0).
```

Only the middle pair -RQy(s=2)+RQy(m=0) cancels in the joined expression.
The second-micro lower and macro upper boundaries, signed S1 integrals,
and all predecessor histories remain. This local adjacent increment does
not replace the earlier first-micro contribution. A cumulative macro pressure already
contains the micro contribution; summing two cumulative outputs would
double count it. The new representation uses source differences or local
increments to retain each source once.

The previous standalone micro report and its outer RQy cap are unchanged.
No old error term is erased, and the new identity does not certify a small
complete pressure error. Signed background integrals have not been
numerically evaluated by this endpoint proof. The old 212 forward terms,
six signed inner/lower contributions and all prior error packets remain.

## Execution, checks and next work

New production exits 0 in 6.436 seconds. Report
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_macro_RQy_seam.json` is 18,041,834 bytes, SHA256 `1ee3ff1234cad6bf79fb24a0a558be50a3c10058a7746c176ea5d4b54d260c3d`, with
2006 input and 377
implementation pins. Independent bounded Luna/max mathematical review,
source/graph binding and saved-only report inspection pass. Separate
runtime worker model metadata is unavailable; the workers have the fixed
Luna/max role. No accepted producer, owner, norm, recurrence, quadrature or
old seam fingerprint calculation was replayed.

Two pre-acceptance writer attempts failed on a nested contract-type check
and an identity-label mismatch. Those runner guards were corrected while
the independently reviewed mathematical module stayed unchanged. No report
was accepted before these fixes.

The new graph suffix and the inherited correlated-I1 suffix each name their
own copy of the common Q4 base. Equal integer suffix handles across packets
must not be identified. A later combined graph must relocate appended nodes
and preserve their source bindings.

Implementation: [endpoint proof](../src/openai_ns_reconstruction/original_micro_macro_RQy_seam.py),
[runner](../scripts/run_original_micro_macro_RQy_seam.py).

```
python scripts/run_original_micro_macro_RQy_seam.py --inspect-report
```

NEXT combine the source-linked second-micro/macro pressure representation
with the accepted first-micro representation, actual signed background
six-I1 integrals and a controlled macro contribution.
Carry the macro upper boundary forward to the next adjacent source before
boxing it separately. Then advance the reserved-pulse terminal functional
comparison, preserving ordinary Z1 and physical scale. Global five-moment
closure, exact heat/energy, admissible stress and flat remainder, genuine
order-dependent time recursion, realizable stress correction and independent
Cartesian NS validation remain incomplete. Goal ACTIVE.
