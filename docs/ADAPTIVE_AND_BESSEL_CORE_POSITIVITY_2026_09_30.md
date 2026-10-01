# Adaptive finite-core positivity and Bessel reference — 2026-09-30

The finite-core positivity calculation now maintains an exact partition of
the entire closed axial interval [-1,1]. Each accepted interval has a
positive directed Bernstein lower bound over every normalized radius
x=R/R_a in [0,1]. All other intervals remain pending or explicitly
unresolved. The union is checked for gaps and overlaps before every save.

Run `python experiments/root_st073/lei_ren_part1_paper_adaptive_core_positivity.py --seconds 180`
to continue the saved partition. Stored binary MP endpoints are exact;
dependency hashes prevent reusing accepted proofs after input changes.
The radial degree remains18 and the coherent preheat-pressure error remains
included. The initial geometry resolves the H0-root scale, but its estimated
root is only used to place cells; positivity relies on interval arithmetic.
Stopping at a run limit leaves a partial proof, not a counterexample.

The latest completed run has142 accepted intervals and216 pending intervals,
with no refinement-limit failures. Accepted cells cover approximately
52.9 percent of the axial interval length; this is only coverage of this
finite-model positivity calculation, not project completion or progress
toward temporal scale recursion. The exact fraction and all endpoints are
in the saved receipt. Initial breadth-first traversal found no positive
coarse cells; its pending partition was retained when traversal changed
to depth-first, with the migration recorded and no accepted proofs reused.

The recurrence now routes delta-dependent factors through the supplied
scalar converter, so interval users evaluate them with directed arithmetic.
The independent coupled/unfactored coefficient fixture still passes.

## Paper-centered route

Cached Lei–Ren text (8.20)–(8.26), (8.47)–(8.59), and Lemma9.1 give

F=F0 Phi, scaled_R=Lambda R, epsilon=1/Lambda,

Phi_model=B(q), q=scaled_R*(chi+epsilon beta)/2,

chi=H0^2/(H0^2+sigma^2), beta=(3-delta/2+jZ)/(1-delta Z^2),

B(q)=sum_n (-q)^n/[n!(n+1)!].

`lei_ren_part1_paper_bessel_core_reference_bounds.py` computes directed
reference bounds for the actual parameters and all |Z|<=1,
0<=scaled_R<=4.1. The structural bound 0<=chi<=1 avoids multiplying
independent intervals for the same rational factors. It yields q<2.1.
The cubic alternating truncation

1-q/2+q^2/12-q^3/144

is a lower bound, is decreasing on this interval, and gives a uniform
reference lower bound approximately0.265381076388888889. The degree18
reference-series tail is also bounded by the first omitted alternating term.
An independent resolved Bessel evaluation checks the endpoint convention.

At H0=0, the epsilon beta term remains: dropping it would erase the
paper's nonzero radial slope there. This reference is not substituted
for the coupled nonlinear core.

The normalized model-comparison error budget for retaining Phi>=1/8 is
approximately0.140381076388888889. If an actual bound of the form
|Phi-Phi_model|<=K8/Lambda is established, the corresponding sufficient
K8 threshold is recorded in the receipt. K8 is a core-comparison constant;
it must not be confused with the later global data-size K containing
reciprocal amplitudes. No actual K8 bound is asserted here.

## Next dependencies

Complete the exact axial partition or center the normalized recurrence
on the Bessel model and bound its coupled defect. For every unresolved
cell, preserve the distinction between interval dependence and an actual
zero. Then combine global positive bounds and derivative bounds into
reciprocal-core C3 estimates, complete frozen-profile K contributions,
and enclose the infinite radial remainder.

The Bessel reference is certified positive; the actual whole-axis finite
core remains uncertified until every partition cell passes. The infinite
core, finite-width stress cone, time recursion and oscillatory correction
remain open.
