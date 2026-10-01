# Core entrance trace and reciprocal derivative bounds

## Results

The adaptive degree-18 finite-core positivity receipt now contains 209 accepted axial cells, 208 pending cells and no refinement-limit failures. Each accepted cell covers the entire normalized radial domain x=R/Ra in [0,1]. These counts are not a project completion percentage.

The reciprocal mixed C3 computation has finished all 142 cells of an immutable earlier positivity snapshot. It bounds all derivatives with total order at most three in x and Z, retaining the inverse-axis amplitude as a logarithmic prefactor. The aggregate upper log is approximately 5e151. This is a finite-polynomial result on the accepted snapshot domain. It does not certify whole-axis positivity, the infinite core, or the full paper constant K. Physical R derivatives require the corresponding Ra inverse powers.

The snapshot SHA256 is 6d595f95ab05c41b0db572ff1e7a86c75548ac743f2613252f6b8c88460e5cb0. The driver checks the snapshot and dependency hashes before resume. Exact MP endpoints, rather than displayed decimals, are authoritative.

## Entrance compatibility diagnosis

At Z=.3, the retained degree-18 / pressure-parameter-order-9 inlet gives

    (I_theta + 2 R F_R) / F
      in [2.378912224157204345383e-28,
          2.378913714738176832892e-28].

This interval excludes zero after including the conditional pressure-integral perturbation currently propagated into the inlet. At chi=1 the bridge instead prescribes angular shear -I_theta. Consequently its log-F radial trace differs from the finite core trace by minus one half of the displayed normalized stress. A small finite residual does not prove C1 matching at this entrance.

The axial stress interval contains zero; this receipt does not establish an axial trace mismatch. Missing pressure-parameter orders and the infinite radial remainder are not enclosed. This is not a counterexample to the paper construction. Nor does the gap-to-collar-width ratio alone establish a new cone condition.

The leading bridge identity T=(1-chi)I follows from its prescribed shear. It cannot be used as evidence that the independently truncated core has exactly zero entrance stress.

Resolving the retained stress polynomial by pressure-parameter degree localizes the dominant defect to degree zero: its coefficient divided by the pressure-free inlet F is approximately [2.3789125968e-28,2.3789133421e-28]. Degree one is enclosed within +/-3.7264524313e-35; degree two within +/-9.335001e-72. These are coefficients of stress divided by the pressure-free amplitude, not Taylor coefficients of the full stress-to-F ratio. Increasing pressure order alone cannot change the nonzero degree-zero coefficient; the radial recurrence and entrance identity need attention. No bound on missing orders follows from this breakdown.

## Next dependencies, in order

1. Identify the exact stress-free core identity and separate the radial truncation defect from pressure-parameter truncation. Use the same accepted analytic preheat datum and coupled core recurrence.
2. Compute the retained stress defect by order, including pressure degrees above nine, without changing the underlying pressure datum or fitting a compensating tail.
3. Establish a nonlinear infinite-core remainder enclosure, potentially using the paper's Bessel-centered comparison and its fixed-point estimates. A Bessel reference bound alone is insufficient; preserve all swirl and pressure coupling terms.
4. Recover a compatible stress-free entrance from the analytic construction. Do not silently set the measured finite trace defect to zero.
5. Finish remaining axial positivity cells and extend reciprocal mixed C3 bounds to the full required domain. Add frozen D/E and inverse-D contributions before certifying K.
6. Bound the bridge comparison error with its vanishing edge factor, then apply the paper's distinct strict and relaxed cone inequalities.
7. Continue functional five-moment matching, exact heat exterior, flat remainder, and genuine n-dependent temporal coefficient recovery. Oscillatory correction and full Cartesian residual validation remain subsequent stages.

## Reproduce

Run from the repository root with Python and mpmath:

    python experiments/root_st073/lei_ren_part1_paper_core_boundary_stress_gap.py
    python experiments/root_st073/lei_ren_part1_paper_positive_cell_reciprocal_C3.py --seconds 50
    python experiments/root_st073/lei_ren_part1_paper_adaptive_core_positivity.py --seconds 180

The reciprocal driver processes the immutable 142-cell snapshot, not newly accepted cells from the evolving adaptive receipt. Its completed receipt should remain paired with that snapshot. No global cone, temporal recursion, or full corrected residual completion is claimed.
