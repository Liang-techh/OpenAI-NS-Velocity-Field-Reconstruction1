# Core g_Z3 and W0_Z3 Green enclosure

The saved actual analytic-core archive contains directed enclosures of
`Phi`, axial `U`, and radial mean `M` ordinary mixed derivatives through
total order five on `rho in [0,4.1]`, `Z in [-1,1]`. This step uses those
same-source enclosures to construct the genuine raw `g4_Z3`, `g5_Z3`, and
`g6_Z3` source sectors on eight core radial integration cells, and applies
the canonical Green kernels to enclose the six components of `W0_Z3` at
the core exit `rho=4` for real `Z in [-1,1]`.

The source intervals used on each cell are the **same whole-core enclosures**,
not independent cell evaluations. They retain the saved nonlinear-correction
and factorial-tail uncertainty and are generally very wide. Signed interval
coefficients are preserved; a cell interval, error bound or midpoint is not
reported as a function value. The Green calculation uses directed whole-cell
inclusion rather than point quadrature. The first cell uses the safe bound
`0 <= t log(1/t) <= 1/2`; it does not evaluate `log(1/0)`.

The highest derivative in the pressure forcing can be formed without an
axis singularity. With `M` the radial mean, `d=1-Z^2`, `L=1-delta*Z^2`,
and `rho*M_rr=U_r-2*M_r`, the exact identity is

```
rho*Q_rr = [2*Z*rho*U_rr - (1-delta)*Z*(U_r-2*M_r)
            - d*(U_rZ-2*M_rZ)] / L.
```

This uses only mixed derivatives of total order at most five and never
divides by `rho`. The source adapter implements the original axial operators
and ordinary Z derivatives. `F0(Z)` remains a positive formal factor with
its saved whole-Z log enclosure, and the outer powers of `Lambda` are kept
separately. Powers of `Lambda` inside the relative-F0 derivative coefficients
are already enclosed numerically; `F0(Z)` is not a fixed global point anchor.
No source exponential is materialized.

The independent synthetic check used direct symbolic differentiation with
a non-unit F0 factor, an axis check, and exact constant and rho-dependent
Green integrals. All 18 checks passed. The actual-data run returned eight
source cells and respectively 8, 16, 16, 16, 8, 16 unmerged Green terms for
`F1`, `Uz1`, `K1`, `P1`, `d_x_F1`, `d_x_Uz1`. These are formal interval
enclosures, not useful relative-accuracy or sign claims. Exact inputs,
source pins, entered record, coefficients and logs are in
`experiments/root_st073/n1_core_g3_word0_v1_receipt.json.gz`.
The staged whitespace check reports one extra blank line at the end of the
frozen `n1_core_forcing_z3.py`; it is retained so the executed source remains
byte-identical to the one-entry manifest. The other staged files pass.

## Same-source inlet join

The new core-exit `W0_Z3` enclosure was joined to the previously computed
first-inlet `B1` factors at the selected narrow Z interval near 0.371. The
join checks the same source family, the exact selected `Cstar` record and
physical norm source hash, the same `Lambda` transfer source, the actual
amplitude definition, and the core/first inlet coordinate. It does not infer
an identity from interval overlap. B1's source factors and the core W0's
`F0(Z)`/`Lambda` factors remain separate; no giant exponential is evaluated.

Thus the missing `B1 W0_ZZZ` term in the previous **seam feedback source**
now has an enclosure: rows 4, 5, 6 respectively add 0, 40, 48 unmerged
product terms to the previously saved 3, 28, 39 known terms. The independent
Cauchy fallback and all-word tail remain separate. The same-source result is
in `experiments/root_st073/n1_word2_inlet_complete_v1_receipt.json.gz`.

This is still a source enclosure at one seam, not the radial Green integral
for `W1_Z2`, an NS residual, or five-moment terminal closure. `P1` here is a
state component and is not `m15`. No cumulative m15 term or prior error was
removed or appended; the 212-term ledger, six signed inner/lower terms and
reserved contribution remain as before.

Next, construct the actual signed `W1_Z2` Green-source enclosures on the
radial cells, including all six Leibniz terms and complete core-to-first
history. This requires W0 jets and matrix jets over those cells, not merely
their seam values; retain the whole-core interval uncertainty, the Cauchy
fallback and the raw omitted-word tail independently. Subsequent physical
mean/pressure and m15 propagation still requires same-source history and
endpoint identities.
