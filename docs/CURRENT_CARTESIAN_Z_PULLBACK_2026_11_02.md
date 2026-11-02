# Cartesian axial pullback for scalar profiles — 2026-11-02

The long-term Navier–Stokes reconstruction remains **active and incomplete**. This checkpoint supplies and tests a conditional coordinate derivative operator for a caller-supplied scalar profile; it does not construct the missing complete velocity, pressure, forcing or NS residual. The 212-term m15 ledger, its six signed inner/lower contributions and all inherited errors remain untouched. P1 is not m15.

## Coordinate distinction

The row-5 Green variable `x=sqrt(rho)` is a normalized radial coordinate, not Cartesian `X`. Earlier references to a physical `dx` in the collar propagator mean the integration measure in that radial x-coordinate, as opposed to its phase parameter. They do not by themselves establish Cartesian `dX`, a physical cylindrical-volume measure, or physical-x cell ownership.

The pinned Cartesian source defines, with `d=1-Z²`,

```text
lambda² = tau/d,   R=(X²+Y²)/(2 lambda²),   rho=R/epsilon,
z_phys=lambda^(1-delta) Z,   t=1-tau.
```

Consequently, at fixed Cartesian `X,Y,tau` and Z-independent `epsilon,delta`, `rho` is generally Z-dependent. For a C² scalar profile `f(rho,Z)`, set `D = d/dZ` along that fixed-Cartesian-transverse path. Then

```text
rho_Z  = -2 Z rho/d,          rho_ZZ = -2 rho/d,
D f    = f_Z + rho_Z f_rho,
D² f   = f_ZZ + 2 rho_Z f_rhoZ + rho_Z² f_rhorho + rho_ZZ f_rho.
```

Thus even at `Z=0`, where `rho_Z=0`, the second derivative generally contains `-2 rho f_rho`. A fixed-normalized-radius derivative cannot simply be substituted for a physical fixed-Cartesian derivative.

Let `L=1-delta Z²`, `B=d/(lambda^(1-delta) L)`, and `C=-(3-delta)Z/d+2 delta Z/L`. Wherever the coordinate map is locally invertible,

```text
partial_z f  = B D f,
partial_zz f = B² (D² f + C D f).
```

The operator also handles the full product rule for `F=lambda^p f` when `p` is Z-constant: with `q=p Z/d` and `q_Z=p(1+Z²)/d²`, `D F=lambda^p(Df+q f)` and `D² F=lambda^p(D²f+2q Df+(q_Z+q²)f)`. This is a scalar-profile operator, not yet a full Cartesian vector or NS operator.

## Implementation and validation

The new leaf requires `f`, `f_rho`, `f_Z`, `f_rhoZ`, `f_rhorho` and `f_ZZ`, plus explicit fixed-coordinate declarations. It rejects missing jets, invalid domains, Python floats and undirected, malformed, reversed or nonfinite exact-MPF interval records. Interval arithmetic is directed; no midpoint is selected. The caller must establish that `rho`, `lambda`, `X`, `Y`, `tau`, `epsilon`, the profile jets and the chosen source/gauge satisfy the stated geometry. The leaf does not independently certify those relationships or a live source owner.

The frozen component ran once and passed **37 synthetic checks**. An independent rational polynomial example checks unweighted and `lambda²`-weighted first and second physical axial derivatives; additional checks cover `Z=0`, missing radial jets, fixed-coordinate declarations, domain guards and exact endpoint records. The pinned source text was matched after the exclusive entry. This is source-definition evidence, not execution of the original source graph or proof that the saved row-5 sectors already provide the required radial jets. The prior 23 MB row-5 receipt was authenticated by raw SHA only, without decompression or re-entry.

The full [receipt](../experiments/root_st073/cartesian_z_pullback_v1_receipt.json.gz), including frozen source bytes and fixture result, has raw SHA-256 `97cd8dcf64d256ad29b84ef45a57c377bdab5b4df96790d1c7a9c4dceda8ba88`. Its source manifest is `839dfa46d44fcbcbd7fd4cbb63deb6e6ce1b3eb5bf89828cf43e21ec13785255`, entry `c8ca00b7d800de5e70a53cee0e98564c95f11eefd277a1be7a2b2a1ea802037d`, and result `b5b80597615565e098c72f760facd33dd6c172b47ef389eed3d4abfd570a345f`.

The fixture retains a trailing blank line in its frozen bytes. The staged whitespace check disables only `blank-at-eof`; staged leaf, fixture and runner hashes match the manifest.

## Remaining work

The current row-5 finite first-word products still omit the baseline `W^(0)` as a separate state contribution: the unknown `R=W-W_finite` includes that baseline as well as further feedback, and must not be called only a higher-word tail. Establish actual source-linked `rho` derivatives for the row-5 expression and state, plus common owner, gauge and physical-cell identity, before applying this pullback to those sectors. A saved whole-domain remainder with the needed Z2 is not yet bound. The independent reserved own-swirl versus negative pulse comparison, functional five-moment matching, exact heat exterior, stress/remainder, time-order recursion, corrections and independent Cartesian NS validation remain open.

Workers were requested using the fixed Luna/max role; runtime model metadata was not exposed for independent verification.
