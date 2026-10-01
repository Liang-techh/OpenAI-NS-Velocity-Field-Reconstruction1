# Same-source general axial-center core factory

`experiments/root_st073/lei_ren_part1_paper_candidate_general_center_factory.py` now recomputes finite coupled core rows at a supplied real axial center or center interval. It does not read the completed Z = 0.3 radial tensor.

For each request, `GeneralCenterCoreFactory.build(center, degree=8, required_depth=3)` performs the following linked computation:

1. Evaluate the exact anchored six-pole amplitude primitive for Lambda = 1e120 at the requested center. The fixed data remain j = 1e-14, delta = 1e-200, logC = 5e151 and sigma = j/500.
2. Generate center-dependent gradient and U0 jets from the analytic axis data. Rebuild ell = -Lambda*g and the jets of S = F0 squared, retaining the nonzero swirl terms.
3. Generate analytic physical-pressure jets from the same accepted fourteen-stage pressure datum. Its accepted schedule SHA remains 736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4. No pressure fit or new pressure tail is introduced.
4. Call the coupled functional radial recurrence with the explicit center and all four jet families. An input length of degree + required_depth + 1 leaves the requested axial depth in the final radial row.

Run `python experiments/root_st073/lei_ren_part1_paper_candidate_general_center_factory.py`. The saved receipt contains two genuine new computations: Z = 0.5, and every local center in [0.49,0.51], with radial degree 8 and retained axial depth 3. It saves exact interval endpoints for seeds and all coefficient rows. All 216 scalar coefficient intervals are contained in the interval-family coefficient bounds. Amplitude branch/root checks pass and its enclosure lies inside the global analytic amplitude enclosure.

A center interval represents a family of local Taylor series about each center in that interval. It is not a single Taylor series about an interval, nor an extrapolation of the old center. The same analytic pressure datum is shared; the core coefficients themselves are recomputed.

## Remaining construction work

- This finite degree-8 factory does not replace the completed degree-124 Z = 0.3 core. Bound the infinite radial remainder for each new interval before promoting it to a completed leading core.
- Provide a resumable version of the explicit-center coupled recurrence for production depth. Do not modify the frozen center-bound legacy recurrence by monkeypatching its center.
- Build a center-aware comparison adapter consuming the new rows, their exact amplitude and the same pressure data. The existing CandidateComparisonJets rejects non-0.3 centers and should retain that guard.
- Propagate actual comparison and transition equations, including analytic axial jets, then generate functional defects and independent moment repairs from those new source data.
- Cover the axial domain by controlled interval families, with particular treatment near the narrow amplitude-gradient regions. A sample grid alone does not prove whole-axis matching.

Original construction-parameter derivation errors and infinite radial tails remain outside this receipt. No new comparison/transition field, terminal functional closure, heat exterior or temporal n-dependent recursion is claimed.
