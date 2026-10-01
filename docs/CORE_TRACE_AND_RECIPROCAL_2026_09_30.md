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

Here "pressure-free" means the coefficient of pressure-parameter degree zero: the accepted physical analytic preheat pressure datum is still present. It does not mean zero physical pressure.

## Angular defect integral implemented

Paper equations (3.13), (3.19), and (8.5) give, at fixed Z,

    d_R [R T_theta] = R D_theta / L,
    D_theta = 2L(R F_RR + 2F_R) - B_theta,
    T_theta(R) = (1/(LR)) integral_0^R s D_theta(s) ds.

The last formula uses axis regularity and L independent of R. The diagnostic `lei_ren_part1_paper_angular_defect_integral.py` computes this integral from the full retained pressure-degree-zero radial polynomial, including nonlinear products. It uses the stored accepted axis data without replacing P0.

At R=4/Lambda and Z=.3 the integral divided by axis F is approximately 6.731847608842599e-29 and lies inside the independent directed inlet stress interval. The absolute sum of contributions from radial defect orders below18 is approximately 2.601107373415409e-261. Thus those computed recurrence equations cancel at the stored source precision; the observed trace is carried by the remaining finite-polynomial defect orders. The first such equation needs radial row19, which the degree18 input does not supply.

This computation is MP arithmetic with stored coefficients, not a directed source-generation error certificate or infinite-tail bound. A total stress enclosure merely containing zero would not establish the exact identity. Completion requires a controlled analytic core solving the equations with regular axis data, plus compatible entrance evaluation.

## Coupled radial continuation to degree20

The new `lei_ren_part1_paper_radial_trace_continuation.py` replays the coupled radial recurrence from the same22 axis Taylor rows, including the same physical pressure datum. It reaches degree20, the highest degree these rows support with the required first axial derivative. The pressure-parameter degree-zero angular entrance trace divided by axis F is:

| Radial degree | Signed trace / axis F | Absolute ratio to previous degree |
| --- | --- | --- |
| 18 | 6.731847608842599e-29 | — |
| 19 | -3.374359703680501e-31 | 0.00501253132832 |
| 20 | 1.533799865309319e-33 | 0.00454545454545 |

The magnitude decreases by a factor of approximately43890 from18 to20. The retained normalized degree18 prefix, including its first axial derivatives, is unchanged in both continuations. An independently evaluated moment trace agrees with the angular defect integral, and the regenerated degree18 trace agrees with the committed source diagnostic.

This extends the actual coupled finite core, rather than resetting its trace or modifying P0. It remains a point-Z, pressure-parameter degree-zero MP computation; no full pressure-parameter sum or infinite radial tail certificate is implied. These finite-order ratios are not extrapolated into a geometric tail bound. Higher orders need higher analytic axis jets and a validated analytic majorant or contraction argument. This spatial radial recurrence is distinct from the requested n-dependent temporal recursion.

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
    python experiments/root_st073/lei_ren_part1_paper_angular_defect_integral.py
    python experiments/root_st073/lei_ren_part1_paper_radial_trace_continuation.py
    python experiments/root_st073/lei_ren_part1_paper_positive_cell_reciprocal_C3.py --seconds 50
    python experiments/root_st073/lei_ren_part1_paper_adaptive_core_positivity.py --seconds 180

The reciprocal driver processes the immutable 142-cell snapshot, not newly accepted cells from the evolving adaptive receipt. Its completed receipt should remain paired with that snapshot. No global cone, temporal recursion, or full corrected residual completion is claimed.
