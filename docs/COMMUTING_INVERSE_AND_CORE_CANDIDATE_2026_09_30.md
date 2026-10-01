# Commuting inverse and a contraction-compatible core candidate

## Tighter angular inverse

The paper's operator A=(chi/2)J2 has a fixed coefficient depending only on Z. It commutes with radial inversion, giving the exact identity

    (A^n f)_(k+n) = chi^n f_k / [2^n (k+1)_n (k+2)_n].

Cauchy bounds chi^n by the nth power of its complex modulus Cchi, without a generic algebra constant at each iteration. For r=h/(eta/2), axial convolution is bounded by

    S(r)=sum_(p>=0)(p+1)^2 r^p=(1+r)/(1-r)^3.

The binomial ratio in the output norm is at most1, and (m+1)/(m-p+1)<=p+1. The radial factorial ratio is at most (n+1)^2/[n!(n+1)!]. Therefore

    ||A^n|| <= S(r) (10 Cchi)^n (n+1)^2/[n!(n+1)!], n>=1,
    ||R|| <= 1 + S(r) sum_(n>=1) (10 Cchi)^n (n+1)^2/[n!(n+1)!].

The identity term n=0 is included separately. `lei_ren_part1_paper_commuting_resolvent_bound.py` sums through n=73 and encloses the infinite tail with directed arithmetic. The resulting inverse norm upper bound is approximately624734.786491445, with logarithm13.34508249704647, compared with the prior generic logarithm390.6324362373922.

Independent weighted-Bessel evaluation contains the full series;1760 coefficient checks support the implementation. The whole-order inequalities above, rather than these finite checks, establish the bound. The same Cauchy convolution also replaces the outer generic multiplier factor256 for fixed analytic coefficients by S(r), approximately80/27. Products of the unknown fields retain the original256 constant.

## Current and candidate contraction gates

The full20-term nonlinear majorant is recomputed with the tighter inverse and fixed multipliers.

| Parameter | Scaled correction size upper | Scaled Lipschitz upper | Half-radius/half-Lipschitz gates |
| --- | --- | --- | --- |
| Existing Lambda=1e36 | log upper25.70930153954 | log upper14.94475260719 | Unproved |
| Candidate Lambda=1e48 | 0.146355076889422 | 3.093311033865688e-6 | Both pass |

The candidate keeps the same accepted analytic P0, j=1e-14, delta=1e-200, h, and logCstar=5e151. It recomputes the complex F0-squared amplitude bound using the candidate Lambda and the certified G modulus; it does not reuse the old amplitude enclosure. The complex Cstar guard also passes at this candidate, and Lambda>=500 and Lambda>=j^-2 hold.

Under the accepted stored-parameter inputs and the paper's complete analytic-space estimates, the candidate map preserves the unit ball about the leading pair and is a contraction there. This supplies an analytic fixed point for the held-fixed datum. It is not a claim that the current degree20 numerical field or its collar has already been rebuilt at this parameter, nor an enclosure of its infinite radial tail.

The current Lambda=1e36 receipt is preserved separately from `lei_ren_part1_paper_nonlinear_candidate_Lambda48.json`. The candidate is a concrete sufficient choice under these bounds; it is not asserted to be minimal or necessary.

The pressure adapter convention has been inspected: `AxisPressureJets.R_a` is continuation metadata, while `ContinuousPreheatPressure` constructs the accepted reference branch contribution5/2 on [0,Rref] and its post-reference stages independently of Ra. Changing the candidate core exit therefore does not change this adapter's P0 formula. This verifies the adapter convention; it does not prove that the rebuilt inner connection satisfies all terminal matching conditions.

## Integrate the candidate coherently

1. Generate candidate axis jets, normalized coupled radial coefficients, velocity, pressure and all five moments from Lambda1e48 and the same accepted datum. Preserve amplitude-factor arithmetic.
2. Derive a candidate analytic norm enclosure from the fixed-point ball and apply the conditional radial tail formulas. Use enough radial and axial orders for required derivative bounds, not only a small endpoint trace.
3. Verify the source-pressure adapter and its inner reference conventions for the new Ra=4/Lambda. Propagate any actual datum change through the accepted schedule, pressure receipts, model norms and contraction gates.
4. Rebuild candidate inner exit and connection data together, with a controlled analytic stress-free trace. Avoid reusing old Lambda1e36 inlet receipts.
5. Complete functional five-moment/heat matching and stress-cone conditions before true temporal recursion.

The original parameter derivation, full reconstructed source errors, candidate radial tail, shared matching layers, global cone, temporal recursion and full corrected residual remain open.

Reproduce from the repository root:

    python experiments/root_st073/lei_ren_part1_paper_commuting_resolvent_bound.py
    python experiments/root_st073/lei_ren_part1_paper_nonlinear_map_majorant.py
    python experiments/root_st073/lei_ren_part1_paper_nonlinear_map_majorant.py --Lambda 1e48 --output-name lei_ren_part1_paper_nonlinear_candidate_Lambda48.json
