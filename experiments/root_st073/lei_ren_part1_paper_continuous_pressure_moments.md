# Continuous pressure moments

`lei_ren_part1_paper_continuous_pressure_moments.py` provides the bounded fifth
moment and the propagated pressure jet for the prepared `Z = .3` outer field.

The source datum is read once at `Rh` from
`source.inner.evaluate_x(mp.e, Z)`.  For every later radius the provider adds
the same increments to both `Mp` and the actual pressure `P`:

* `d Mp / d log R = Utheta^2 / 2` and
  `d Mp_Z / d log R = Utheta Utheta_Z`;
* preheat reference stages use the shared angular schedule, exact exponential
  atoms for constant stages, and order 96 Gauss quadrature only for variable
  stages;
* the two signed compact bump atoms remain separate and use
  `exp(-(1 + 2 mu) t) (d h + d h^2 / 2)` with its analytic Z derivative;
* the heat region uses `ContinuousHeatMoments.pressure_increments`, preserving
  separate reference and correction terms, order 96 collar quadrature, and
  analytic exterior polynomial atoms.

The public aliases are `moments_jet`, `moments`, and `pressure_jet`.  Each row
contains `p`, `p_Z`, `P`, `PZ`, separate reference/bump/heat increments, the
local radial integrand, and explicit scope flags.  Full five-moment closure,
finite-energy certification, stress/PDE certification, and quadrature error
enclosures remain false by design.

Run the prepared checks with:

```powershell
C:\Users\z5242\AppData\Local\Programs\Python\Python311\python.exe `
  experiments/root_st073/lei_ren_part1_paper_continuous_pressure_moments.py
```

The run writes the companion JSON receipt.  It checks the `Rh` anchor, radial
integrands at `Rref`, `Rp`, flattening, bump, and both heat regions, and checks
the signed bump and heat corrections independently of their much larger
reference sums.
