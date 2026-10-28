# Quiet X quotient jets and the higher-row gate — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new frozen leaf authenticated the raw Rp source archive before decoding, the quiet-flow source and an existing Rp receipt. It bound the saved incoming `E` and `h` ports and the current quiet `X` ports through ordinary Z order five. It did not traverse expression edges, evaluate source functions, call an ancestor, or remove old errors. The old receipt's `X_Z0`–`X_Z5` labels and passed flags are inherited evidence; they were not rerun by this leaf.

Let `A(Z)` be the incoming modified `E`, `h(Z)` the incoming modified angular history, and let `s` be the quiet coordinate or its Rp endpoint. The source's exact finite-mu formulas give, with `lambda=1-mu`,

```
E = A exp(-(1/2+mu)s)
H = exp(-3s/2) [h + A s exprel(lambda s)]
X = H/E = T W + K
W = h/A;  T = exp(-lambda s);  K = T s exprel(lambda s).
```

The simplification preserves the entire `exprel` continuation, including `lambda=0`; it does not evaluate an integral or substitute a bound. It assumes nonzero `A`, same-source ordinary C5 input jets for `A,h`, and `mu` and `s` (including the endpoint length) independent of Z. Those source/geometry premises remain to be proved for the actual port.

The new callable `quotient_jets` computes ordinary derivative jets `W_n` from supplied `A_n,h_n`:

```
W_0 = h_0/A_0
W_n = [h_n - sum_(k=1)^n binom(n,k) A_k W_(n-k)] / A_0.
```

Consequently, under the stated assumptions, `X_0=T W_0+K` and `X_n=T W_n` for `n>=1`. Fifteen exact symbolic checks passed: the source-ratio reduction, six Leibniz identities, direct quotient differentiation through order two, and propagation checks at orders 0, 1, 2, 4 and 5. The symbolic coefficient checks do not independently differentiate or evaluate the archived source graph. This is a **spatial Z-jet quotient recurrence**, not the order-dependent temporal coefficient recursion required later in the project.

The saved incoming `E` references are `[262708,262709,262710,262711,262712,262713]`, the incoming `h` references are `[262726,262727,262728,262729,262730,262731]`, and both saved quiet and current-frame `X` references are `[265772,265776,265782,265789,265797,265806]`. Reference agreement does not establish the functions' derivative identities. The inherited Rp receipt lists all six X rows, but its checker compares against the old canonical expression with additional source assumptions; it is not a new proof of current inlet-rewrite transport.

The formula gives a sharper all-Z compatibility criterion: for real finite parameters `T` is nonzero, so the current X4 and X5 vanish exactly when the corresponding ordinary derivatives of `W=h/A` vanish, provided the premises above hold. On an interval, `W'''' identically zero` implies `W` is a polynomial of degree at most three, hence `W'''''` also vanishes. Proving `W` constant would suffice but is stronger than required. None of these all-Z conditions has yet been proved for the actual modified incoming source; a node ID or a finite collection of samples cannot establish it.

Actual current X C2/C5 jet admission, geometry independence, inlet substitution and axis transport, Gamma improper-integral differentiation, X4/X5 compatibility, physical seams and m15 cancellation remain open. No source numerical value, certified bound, global five-moment closure or later physical layer is claimed. All earlier errors, nonzero histories, the independent reference and once-only reserved contribution remain; `P1` remains distinct from `m15`. Runtime model metadata was unavailable (`null`).

Next: prove the actual modified compact-exit `h/A` source function and its same-source ordinary Z jets, using bounded source interfaces and preserving the geometry and substitution gates. Apply the formula to current X1/X2 and prove or disprove the all-Z fourth-order condition before accepting the old constant-X C5 continuation.
