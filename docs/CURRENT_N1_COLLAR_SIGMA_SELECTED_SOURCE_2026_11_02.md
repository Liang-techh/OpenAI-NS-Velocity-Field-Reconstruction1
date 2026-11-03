# N1 lower-collar sigma on the selected source domain (2026-11-02)

The full reconstruction remains ACTIVE and incomplete. This entry evaluates a previously uncomputed rational enclosure on a domain deduced from the saved selected `sc` source record. It is a conditional source-domain result, not a live graph binding, an m15 value or terminal five-moment matching.

## Actual selected input and direct lower collar

The saved ordinary JSON `current_inner_exit_strict_collar.json` contains identical lower and upper exact binary tuples `[0, M, -679178314328189951, 798]`, with `M.bit_length() == 798`. The lower-collar source uses that selected value for `sc`, constructs `sc/4` and `sc/2`, and defines its direct shear factor with `sigma(s)`. The new consumer checks tuple equality, sign and bit length and proves, without materializing `sc`,

```text
0 < sc/2 < 2^(-679178314328189154) <= 2^-128.
```

It applies the prior rational interval leaf once to the new enclosure `[0, 2^-128]`. After preserving the leaf result, it intersects the derivative upper bound with the independently documented endpoint majorant `8!*(9/8)^8*(56/27)*(2^-128)^13`. The saved exact rational outputs imply:

```text
0 <= sigma(s) < 2^-2032
0 <= d sigma(s)/ds < 2^-1646
```

These are enclosures over the stated source domain `[sc/4,sc/2]` and the subintegral variables within `[0,s]`, conditional on the saved source construction and fixed source definition. They are bounds, not sampled values or an evaluation of `chi=1+(hb-1)*sigma(s)` as a whole. The selected endpoint is not a finite-N error bound or a midpoint substituted for a function.

## A different normalized cutoff

The first-chart N1 cutoff instead has `q=s/(sc/4)` and argument `x=4q-2`; it must not use the tiny direct `s` interval as its sigma argument. For `q∈[0,1]`, its argument spans `[-2,2]`, with exact flat branches `sigma=0` for `q≤1/2` and `sigma=1` for `q≥3/4`. The consumer separately called the rational leaf on the active argument `[0,1]` and saved the resulting interval `[0,1]` and `d sigma/dq∈[0,144]`. This separately labeled calculation is retained in the receipt but explicitly excluded from the runner's accepted direct-lower status; it does not bind normalized cutoff node `52490` or establish live owner identity.

For ordinary slow `Z`, the runner explicitly supplies fixed-`s` and fixed-`s/sc` conditions, under which the respective sigma factors have derivative zero. Those are conditional coordinate assumptions consistent with the source geometry, not a proof that an arbitrary moving coordinate has zero derivative. They do not set the ordinary Z1 derivative of the full N1 fields, history, pressure or m15 to zero.

## Verification and retained evidence

A single frozen run passed **25 / 25 synthetic checks**, then called the new consumer once on the saved source mapping and the published rational evaluator twice on the two different intervals. It did not re-enter an old fixture, execute source graphs, decode gzip, or evaluate the physical N1 field. Source files were hash-pinned, not imported as live owners. The archived result contains exact Fraction endpoints, tuple, conditional Z flags and all frozen hashes. Archive: `experiments/root_st073/n1_collar_sigma_actual_domain_v1_receipt.json.gz`; raw SHA-256 `b42ad7a699de3d9bc501b6f8f9d8040919079cafd214ce7bb9274c639b90d209`. Source checkpoint SHA-256 `6ccb4455bfc803531dc89123c82aae59c4bdcec7ec845042e4e6bbb5d7c54549`; entry SHA-256 `91a8124a4654b00f16929a9c0b3a7412457aab75074229249b2f8dadac402b85`. Authenticate raw gzip before decoding; do not re-enter the saved run.

The shared live graph/binder and actual owner chain remain unproved, as do the whole-Z analytic core oracle, corrected pre-Rp signed correlation and physical terminal matching. Native and corrected lower endpoints remain separate. Keep the 212 forward terms, six signed inner/lower contributions, reserved contribution once and every inherited error; P1 is not m15.

## Next executable task

Trace the exact source-construction chain and live binder for the first-chart cutoff and both lower ports; require actual owner/history/physical-scale correspondence before attaching either enclosure to the saved physical graph or cancelling endpoint terms. In parallel identify a source-owned whole-cell analytic-core coefficient and radial-tail provider, keeping the n0 core data separate from the n1 temporal response. Continue toward the independent negative complete-pulse m15 comparison without repeating ancestor calculations.
