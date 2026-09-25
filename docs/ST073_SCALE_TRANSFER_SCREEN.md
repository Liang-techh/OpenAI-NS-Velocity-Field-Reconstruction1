# Sequential scale-transfer screen

`experiments/root_st073/scale_transfer_recurrence.py` starts with the
width-six joint correction and fits the next amplitude vector in each of
the bands `k=1..3` and `k=3..5`, freezing previous knots. A quintic
smoothstep makes the time-dependent amplitudes C2 at the joining scale.
Each correction uses the existing divergence-free, endpoint-preserving
streamfunction and swirl modes. The Cartesian finite-difference momentum
evaluator sees the changing amplitudes at every time stencil.

The bounded fit finds increments of only about `1.3e-10` in each band.
At off-grid `k=2.25,4.25`, the largest sampled momentum peak remains
`1115.4877104` and the larger physical-volume L2 remains `2.5520260`
to numerical precision. Interface velocity and pressure changes are zero
at the sampled endpoints. This particular sequential transfer map gives
no nontrivial improvement over its fixed-scale starting state.

This is an experimental recursion algorithm, not a derived scale law or
interscale error bound. It does not establish a critical-time limit,
whole-space finite energy, the paper's oscillatory stress correction,
or either `1e-3` momentum gate. The next meaningful route is to solve
the evolving shear/phase/amplitude and mean-correction dynamics from
the official paper, with matched inner and outer data, before claiming
a physical scale recursion.
