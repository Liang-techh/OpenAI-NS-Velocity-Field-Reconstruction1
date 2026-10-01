# Separate Lambda120 regeneration

The Lambda120 candidate is selected because the weighted angular derivative
estimate passes on the full analytic axial interval. The previous degree124
Lambda48 field remains a separate calculation.

The new core driver imports the existing gauge recurrence without modifying
its source or completed states. It uses the same accepted fourteen-stage
physical pressure datum, and a separately enclosed anchored amplitude.
Source hashes, target parameters, pressure inputs, row counts and axial depths
are checked before resuming. Exact interval endpoint tuples remain the
authoritative restart data.

## Anchored physical amplitude

For `g=L*H/(H^2+sigma^2)`, the six simple poles are the roots of
`H(r)=+/-i*sigma`. Their residues are `L(r)/(2*H'(r))`. The amplitude module
evaluates the anchored primitive

`G(Z)=sum_r L(r)/(2*H'(r))*log((Z-r)/(Z0-r))`,

then preserves `F0=exp(-5e151-Lambda*G)` without rescaling. All six poles and
the real anchor have strict Rouche disk enclosures of radius `1e-240`.
Root separation, logarithm branch consistency, the real primitive, derivative
identity and positivity must pass before the API returns a usable amplitude.
At `Z=.3`, `G=8.620576621613044845639278153395...`; an independent high-precision
integral is enclosed by the primitive interval. That integral is a diagnostic,
not the construction used to bound the amplitude.

The true amplitude interval is contained in the earlier coarse global
amplitude bound. The new core uses this narrower interval and retains its
squared coupling in every applicable recurrence term.

## Bounded continuation

From the repository root:

```powershell
python experiments/root_st073/lei_ren_part1_paper_candidate_Lambda120_pipeline.py --seconds 180
```

Repeat the command until `finite_core_complete` is true. Each radial order is
saved atomically. After order124, the command generates separate Lambda120
combined error, angular trace and shared five-moment inlet receipts. The
pipeline status distinguishes these from terminal matching and recursion.
The seconds option bounds recurrence time between complete radial steps;
initialization, serialization and downstream diagnostics add wall time.

The normalized mixed C3 budgets must pass before the downstream receipts are
marked regenerated. These budgets apply at `Z=.3`, throughout `s in [0,4.1]`.
The analytic weighted sign estimate has whole-axis scope; the finite core
does not yet have that scope.

The weighted sign receipt records the original paper-text cache hash under
`work_paper_cache/lei_ren_part1.txt`. Replaying its source provenance currently
requires that local reference cache as well as the tracked numerical inputs.

## Still open

- Recover all five terminal moments as axial functions on a shared domain.
- Construct and validate the coherent transition/flatten/heat collar.
- Establish admissible stress-cone margins in every matching region.
- Implement the distinct n-dependent temporal recovery equations and repairs.
- Construct oscillatory corrections and independently validate the full
  corrected Cartesian Navier–Stokes residual.

None of those items follows from degree124 completion or stress intervals
containing zero at a single axial center.
