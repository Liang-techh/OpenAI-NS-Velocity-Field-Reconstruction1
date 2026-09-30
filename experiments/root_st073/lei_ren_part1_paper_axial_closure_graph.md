# Shared-atom closure of the actual quadrature rows

Run `python experiments/root_st073/lei_ren_part1_paper_axial_closure_graph.py`.
The solver now exports its actual canonical 192-node matrix, scaled incoming
rows, pulse rows, and mu. The graph uses that matrix, not an invented fixture.
Incoming and pulse quantities remain shared atoms; their huge exponents must
not be converted into literal powers-of-ten integer denominators.

For A c + b + a p = 0, Cramer coefficients are exact rational linear forms
of b1, b2, a*p1 and a*p2. Expanding the two terminal row expressions cancels
all terms exactly. Differentiating these forms also cancels all terms when
the matrix is Z-independent and the shared atoms are differentiable in Z.
The current receipt verifies both formal row identities. This is a conditional
functional statement, not a measurement of the continuous Z derivatives.
Row 2 is the mixed angular/axial moment; it is NOT the derivative of row 1.

The graph deliberately separates the exact dependency identity from its
rounded coefficient replay (about 1.20e-416 and 2.94e-417). The coefficient
drift from the existing solve is about 1.26e-416. Energy replay using graph
evaluations is about 5.36e-443; the exact algebraic energy root is not certified.
A changed independent basis integral (matrix entry perturbed by 1e-100) gives
nonempty residual terms instead of being silently replaced by zero.

The physical mean has a different requirement. If continuous weights differ
from the quadrature weights, its remaining numerator is

    (A_cont-A_quad)c + (b_cont-b_quad) + a(p_cont-p_quad).

Nonzero error here still produces a 1/r radial tail, regardless of tolerance.
Moreover Mz_Z must close in the same continuous representation throughout a
Z neighborhood. The existing actual nonzero tail has NOT been overwritten.
No finite-energy claim or exact continuous integral claim is made.

Next implement one shared continuous/conservative basis abstraction for:
pointwise bump/pulse values, partial primitives, full primitives, and Z jets.
The row solve and radial transport must consume those same objects. Record
whether each full atom is a continuous integral, a certified enclosure, or a
conservative numerical basis definition. A conservative model must also
recover its velocity by differentiating its primitive and independently
compare it with the source profile; do not keep a different pointwise velocity
and merely force the full primitive to zero. Resolve incoming core/annular
primitive provenance and both linear rows, then the shared energy root.
