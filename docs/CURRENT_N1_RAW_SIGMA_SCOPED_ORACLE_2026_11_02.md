# N1 raw sigma scoped-oracle connection (2026-11-02)

The reconstruction goal remains ACTIVE and incomplete. This entry connects two raw first-chart sigma leaves to the existing scoped graph evaluator under an explicitly supplied radial interval. It does not establish a live owner, execute a complete `g5` root, compute physical `m15`, or establish terminal five-moment closure.

## Real graph leaf execution

The authenticated raw graph archive has SHA-256 `1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe`. Its checked rows are:

```text
35848: bound variable original_lower_collar_s
35949: exact_original_source_flat_sigma(argument=35848)
35953: original_flat_sigma_prime(argument=35848)
```

The runner checked these rows and their source-module SHA records, supplied binder `35848` with the exact interval `[0,2^-128]`, and invoked the existing scoped graph evaluator once for each root `35949` and `35953`. The new oracle takes the respective exact rational bounds from the authenticated previous lower-collar receipt; it checks the lexical node ID and interval inclusion, the graph row and target graph hash. It does not invoke the old interval calculation again. The `35949` result uses the saved sigma interval and `35953` uses the distinct saved input-derivative interval. Both graph evaluations returned `conditional`, `certified=false`, with `source_identity_verified=false` and `owner_binding=false`. The saved interval upper bounds remain bounds, not point values.

This is the direct raw `sigma(s)` and `d sigma/ds` path. It is not the separate normalized inner cutoff `sigma(4s/(sc/4)-2)` or node `52490`; the input derivative is not an ordinary `Z` derivative. The prior source-domain argument places the relevant first-chart range inside the supplied interval, subject to the source-construction assumptions recorded in the previous handoff. Authentication of saved bytes and matching graph records do not prove live graph owner identity.

The new fixture uses synthetic interval data and passed **19 / 19 checks** for the fresh v2 entry. The actual saved graph and previous receipt were raw-hash authenticated before gzip decoding; the real leaf evaluations were each called once after the fixture passed. No integral, complete raw `g5`, prior physical source or old fixture was executed. The real evaluator results are not certified physical-function identities. The two new oracle-supported leaves are among the previously recorded 33 opaque instructions; all other source, integral, substitution and core-provider obligations remain.

## Failed v1 and preserved evidence

The first component entry failed its synthetic fixture after 3 / 19 checks, because the fixture used an undefined local `F` alias. Before that failure it had called the new oracle on synthetic sigma and sigma-prime data; the real graph and old receipt had not been decoded and no real root was evaluated. Its unchanged frozen code and full failure are archived in `experiments/root_st073/n1_raw_sigma_scoped_oracle_v1_failure.json.gz`, raw SHA-256 `247d1c6abed8d30321293a885feb302985c246b4f6245bc91ff43725e9965189`. Do not re-enter the failed v1 run.

The v2 runner and fixture are distinct; v2 reuses the frozen v1 oracle implementation and uses different synthetic bounds. Its successful receipt is `experiments/root_st073/n1_raw_sigma_scoped_oracle_v2_receipt.json.gz`, raw SHA-256 `698f558d6ae5d179ea7700a48fc0af2d31d630be64fe1c24b639cb04e68d507a`. The v2 source checkpoint SHA-256 is `26a39326f41684759d08c38200548a300dc0e506832d9c4c0b9c605a7eacf933`, and entry SHA-256 `86dc51da7d791cd94e73218a42b5f5ba1f38b82b874c467fec5647addb624d6d`. Archive contents retain exact bounds, per-node results, frozen code and both input hashes. Authenticate raw gzip bytes before decoding; do not re-enter v2.

## Corrected raw g5 formula

`CURRENT_N1_RAW_G5_GRAPH_2026_11_02.md` previously printed `-2 + delta` as the raw axial weight. The actual source uses `-2 - delta`:

```text
Nt = -axial_second(g, F, rho, delta, -2 - delta)
g5 = 2*epsilon*Nt
```

The document is corrected. The separate feedback coefficient `b1=-2+delta` is a different quantity and was not changed. This correction is from source inspection; it is not a new numerical evaluation.

## Remaining work

An actual `g5` C1 oracle still requires scoped integral binders, source-owned analytic-core mixed jets and full tails, and common coordinate/scale/owner proof. Cross-owner native and corrected lower endpoints remain distinct. Preserve the 212 forward terms, six signed inner/lower contributions, reserved contribution once and all inherited errors; `P1` is not `m15`. Next, prove the live integral binder/source identity without replaying prior calculations and provide actual whole-Z core C1 coefficients and tails before attempting the independent negative complete-pulse m15 comparison.
