# Partial core second-word K0 exit — 2026-11-02

This v2 pass constructs a finite K0 second-feedback exit at `x=2` from the authenticated saved eight-cell first-feedback profile and raw `C0/B1` rows with a common saved source lineage. It is a conditional coefficient enclosure, not a complete physical core state or a five-moment closure. `P1` is not `m15`.

## Row-5 algebra and missing terms

For row 5, the saved zeroth ordinary-Z order is

```text
Z0 = C0[0] W[0] + B1[0] W[1].
```

The known first order is

```text
Z1 = C0[1] W[0] + C0[0] W[1] + B1[1] W[1].
```

`B1[0] W[2]` is still missing in columns 1, 2 and 3: 24 unknown terms. The row-5 `C0` row is nonzero in columns 1, 2, 3 and 5; `B1` is nonzero in columns 1, 2 and 3. Columns 4 and 6 are structural zeros.

## Kernel partition and retained output

At `x=2`, with `t=x/2`, cell `i` is `[a,b]=[i/8,(i+1)/8]`. The exact kernel masses used were

```text
M_F  = (b^2-a^2) - (b^4-a^4)/2
M_Fx = (b^4-a^4)/2.
```

The `dx` and derivative scales are included. Each of the `F1` and `d_x_F1` fields has 4,824 directed-interval/formal-power records for `Z0` and 8,568 for `Z1`. These are enclosures, not function values; formal powers are retained, and no records or uncertainties are canceled. The fixture passed 27 synthetic checks; the leaf was then called once on saved inputs.

## Receipts

The v1 attempt failed in its synthetic stage with `KeyError: '0'` before reading any actual receipt. Its frozen failure archive has SHA-256 `a229abba26f655ad15bf5e02e8bad2974c858bf9f6b9a67273b33a73fbfea709`.

The successful v2 receipt is [`core_secondword_partial_v2_receipt.json.gz`](../experiments/root_st073/core_secondword_partial_v2_receipt.json.gz), SHA-256 `d7de2dbcd1fbc90bf039d0b8527cc82be6e7826e475c2c0c1fdfa6ed3c8585fe`. Its source manifest SHA-256 is `ebf53e8dc6eca94921e1212e2121d3b33dfa66d0d191a621d8b7e0c7b9415864`, entry SHA-256 is `edc0124a755f007470e86202429b3659d0fdb018f8ffc91880c371d58e976def`, and result SHA-256 is `ddd42a3d51413d4dab54cb1dd82d61f51e58b4d49f528839e2ffe9ac09664419`.

The archived source, entry, full result, and frozen executable bytes were checked against the saved files. Staged executable hashes match both manifests. Frozen files retain their original trailing blank lines and the v2 runner's CRLF line endings; the staged whitespace check accepts CR at end of line and disables only blank-at-EOF checking.

The saved family/log and source chain are jointly traceable. They do not prove a complete physical owner, Kchi or gauge compensation, inlet or axis corrections, the all-word tail, or the global construction. This pass was not added to the earlier first-collar terms and did not change the 212-item ledger; five-moment closure must not be claimed.

## Next independent work

Build the first-word `Z2` radial profile independently in every core cell, retaining higher-order dependencies and certified errors. Before connecting it onward, verify the common owner from the core exit to the first collar, the scale factors, and absence of double counting.
