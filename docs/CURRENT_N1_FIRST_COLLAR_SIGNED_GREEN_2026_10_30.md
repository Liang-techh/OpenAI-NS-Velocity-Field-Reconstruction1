# Signed first-collar Green computation interface, 2026-10-30

This checkpoint adds a computable signed matrix-jet and Green-integral path toward the actual first-collar `W1_Z2`. It does not yet evaluate the original source, bind the formal source sectors, complete the all-word series, compute m15, or prove five-moment closure. The earlier `n1_word1_z2_first_collar_v1` result, its separate local error radius and all earlier source/error records remain in force.

## Implemented computation

`src/openai_ns_reconstruction/n1_first_collar_matrix_jets.py` independently evaluates the original raw system's dynamic `C0` and `B1` entries with ordinary Z jets of order zero through two. Its inputs are the caller's same-source `F`, `Fr`, `FZ`, `U`, `Ur`, `UZ`, `Q` jets and the original fixed-Z geometric parameters. It removes exactly the three fixed kinematic `B0` entries; the resulting 10 `C0` and 6 `B1` entries retain their signs. The caller must provide coherent actual field and derivative jets, including `FZ` and `UZ`; the leaf does not establish their origin.

`src/openai_ns_reconstruction/n1_first_collar_signed_green.py` supplies the order-major matrix adapter and a directed-interval Green integrator. It transports all six supplied core-exit `W0` Z0–Z3 jets into the collar and integrates the signed local `g` Z0–Z3 sources. At each outer cell it combines those values with the local matrix jets using

```
h_ZZ = C0_ZZ W0 + 2 C0_Z W0_Z + C0 W0_ZZ
     + B1_ZZ W0_Z + 2 B1_Z W0_ZZ + B1 W0_ZZZ.
```

It then integrates this feedback with the original six signed kernels, using `u=2 exp(hb*v/2)`, `du=hb*u*dv/2` and `v=s*t`. It adds the supplied core-exit `W1_Z2` history after the homogeneous transport. `Z` is differentiated at fixed phase and equivalently at fixed physical radius when the same-source scales are Z-independent. Cell range integration requires the callbacks to enclose their complete input phase intervals; the cell intervals and final finite-word sums are returned. The source callbacks, core data and family identifier are caller assertions, not certified by this leaf.

The current integrator materializes `hb` as an interval and combines source terms into interval sums; it **does not preserve formal `hb`, `F0`, or `Lambda` source-sector identities**. It must not replace the source-factored checkpoint or support source-correlation or endpoint-cancellation claims. It also does not supply or subtract any all-word tail. For original-source application, a new adapter must establish exact same-source identity, retain each formal sector and its physical factors, and provide certified signed first-collar `g_Z0` through `g_Z3`, leading jets and complete core inlet history.

## Checks and actual status

The one-entry synthetic fixture returned successfully with 60 checks: 48 signed matrix entries/ordinary derivatives against independently differentiated symbolic expressions, six Green results against a high-precision direct integral using separately derived analytic W0 source masses, and six inlet identities. The comparison used nonzero synthetic data; it is not an evaluation of the paper source. The source/entry/result receipt is `experiments/root_st073/n1_first_collar_signed_fixture_v1_receipt.json.gz`, SHA-256 `2e8e8a799efb30ce371a3a46992d8449a3fdc526127c2368be858e14d44b1f71`. The one-entry result SHA-256 is `a78ba42a259c0baa3ad1744fa735b10c18484b4a53c35751092b266ccc433daa`.

The first-collar source graph defines the functions, but the available numerical inputs do not furnish a complete signed first-collar `g_Z3` or independently bound all the required leading jets. A local finite-word enclosure is neither the complete n1 field nor an NS residual. `P1` remains distinct from the m15 moment.

## Next executable step

Build source-factored interval adapters for actual collar leading jets and forcing through ordinary `Z^3`, including the full core-exit `W0` history and original `hb/F0/Lambda` factors. Use them to produce signed local Green sectors with a separate certified integration remainder; retain the earlier finite-word and all-word error records. Only after those bindings are established connect the independently defined own-swirl and pulse m15 terms, without duplicate terms or unsupported endpoint cancellation.
