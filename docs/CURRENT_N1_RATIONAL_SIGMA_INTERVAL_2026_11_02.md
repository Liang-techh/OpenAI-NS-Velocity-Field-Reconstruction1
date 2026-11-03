# N1 flat sigma: independent rational interval leaf (2026-11-02)

The long-term physical reconstruction remains ACTIVE and incomplete. This entry adds an independent exact-rational interval implementation of the flat cutoff formula that the N1 source describes. It is a mathematical leaf, not a live source-owner binding, numeric evaluation of a physical graph, terminal match, or Navier–Stokes validation.

## Implementation and evidence

`src/openai_ns_reconstruction/n1_sigma_interval_rational_v1.py` encloses both the cutoff and its derivative with respect to its **input** `x` on a closed interval with `Fraction` endpoints, including cells crossing `0`, `1/2`, or `1`. The implementation splits support and symmetry branches; derives endpoint derivative bounds rather than dividing by a zero endpoint; and encloses exponentials with alternating rational Taylor bounds or a conservative positive-series fallback. It uses exact integer/rational arithmetic and explicit resource limits. A returned upper bound is a bound, not a selected function value. The implementation does not import or execute the old point-sigma family or physical construction.

After code, fixture, runner and the two source-definition files were frozen, one synthetic fixture call returned **35 / 35 checks passed**. The tests include independent rational Bernoulli brackets, exact support and midpoint identities, reflection, cells crossing the support boundaries and midpoint, input/control rejection, and resource guards. For the fixture's test interval `[-1,2]`, the returned enclosure was `sigma ∈ [0,1]`, `d sigma/dx ∈ [0,36]`; the latter is a conservative bound. These checks support the implementation but do not establish graph ownership or equivalence of a particular live object. The fixture uses synthetic arguments only; no actual N1 input interval was evaluated.

The one-entry runner and complete saved result are in `scripts/run_n1_sigma_interval_rational_v1_once.py` and `experiments/root_st073/n1_sigma_interval_rational_v1_receipt.json.gz`. The raw archive SHA-256 is `ec96aef18d6456f0b30528c933ca83cc284f8745667b0311d296f4b6e3f1d6c1`; the source checkpoint SHA-256 is `6f33933c4d144d7e60be36e142805e66544f55a3e5fdef8a654689eba5a5764b`. The archive preserves frozen code bytes, input hashes, entry, fixture result and trace if any. Authenticate the raw gzip hash before decoding. Do not re-enter the saved one-entry run.

## Physical and provenance limits

`d sigma/dx` is not ordinary `Z1`: any later composition must use the actual input's ordinary Z derivative and the chain rule. `source_identity_verified=false`, `owner_binding=false`; the source files were hashed but not executed or used to certify an actual owner. Seven saved N1 integral occurrences still have no verified binder/parent object, and the actual analytic core's whole-cell mixed-jet/tail oracle is unavailable. The new interval result cannot replace that oracle.

Read-only tracing found the native lower port `[129567,129568]` belongs to a Q4 first-micro boundary, whereas corrected `[53450,53462]` is relocated from the corrected patch. Each owner retains its own source/history chain, but a shared live `R/Pstar/delta/N`, physical-coordinate and incoming-history object has not been established across them. Native coordinate metadata `66988` and corrected metadata `29` are distinct identifiers, not by themselves proof that physical coordinates differ. Keep both sides and all error terms; do not infer cancellation from matching formulas or overlapping bounds. Four-cell directed repairs and prior conditional identities were already computed and must not be repeated.

The 212 accumulated m15 terms, six signed inner/lower contributions and inherited errors remain unchanged; P1 and m15 remain distinct. No global five-moment closure is claimed.

## Next executable task

Connect the interval leaf to the **actual** N1 cutoff only after obtaining a source-construction-time witness of the shared live graph, binder and caller/owner chain, and a real interval for the cutoff argument and its ordinary Z derivative. Such a witness cannot be fabricated from matching names or caller-supplied tokens. In parallel, identify source-owned, whole-cell coefficient and radial-tail data for the analytic core, and trace the cross-owner lower history and physical-scale connection without re-entering earlier computations. Until these exist, keep the new interval leaf separate from claims about actual `g5`, corrected pre-Rp or m15.
