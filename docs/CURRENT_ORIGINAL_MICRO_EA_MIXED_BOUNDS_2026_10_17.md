# Actual micro E/a mixed derivative bounds

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. This
milestone provides a partial actual-background provider on both micro charts,
for slow coordinate s in [0,2] and Z in [37/100,93/250] = [.370,.372]. It
encloses ordinary Z orders 0 through 3 and physical y orders 0 through 1 for
E and a. It does not provide the full mixed-root or terminal/global closure.

## What this receipt establishes

The comparison direction is differentiated from its own source histories. With
d = 1-Z^2, the exact leading-history term is

```
B0 = (1-delta)*Z*JH/2 + d*JK
Dbar = R*D_comparison = rho/L * [gaxis*B0/(2*barphi)
                                  + epsilon*(-W + A0/(2*barphi))].
```

The d*JK coefficient is retained; it is not delta*JK. The actual loop
variables are bound by

```
a = chi_R*Dbar,
E_y = (1-a)*E/2.
```

Here y=log(R/Ra), R=Ra*exp(hb*s), and dy=hb*ds. Ordinary Z
differentiation holds chi_R fixed. On the first chart, the
physical-y derivative of chi_R contributes the explicit 32/hb factor, so the
provider records the first-chart bound

```
|a_yZ^k|  <= 32*|Dbar_Z^k|/hb + |Dbar_yZ^k|.
```

On the second chart chi_y is zero and the corresponding bound retains the
physical factor hb:

```
|a_yZ^k|  <= hb*|Dbar_yZ^k|.
```

The E derivative entries are relative ratios E_Z/E and E_yZ/E. They do not
choose E=1, and they are not signed field values or differentiated C0 caps. The
provider uses direct source graph ports for the comparison histories, native
fields, geometry and loop branches. It computes only the three fresh Xh rows
(0,4), (1,3), (1,4); the existing rows are consumed by their saved receipts and
no accepted ancestor producer, recurrence, norm computation or quadrature is
replayed.

The stored exact Xh and derivative entries are outward norm enclosures. They are not signed values. In particular, this receipt does not
bound a finite-N modified E_N, and it does not differentiate any C0 cap into a
derivative estimate.

The all-coefficient derivative embedding used for the three new rows is

```
(i+k)! / [20^i * h^k * (k+1)^2 * (1-rho_max/20)^(i+k+1)],
```

with rho_max=4.1 and the same saved Xh parameter. The actual Phi norm is
bounded by the model norm plus the uniform self-map size, and U-U0 uses
epsilon times the actual Psi norm, once. The saved map is uniform on the
unit ball and includes all twenty terms; at its fixed point X*=X0+N(X*),
the correction is bounded directly by that uniform self-map size. No
fixed-point norm or coefficient recurrence is rerun.

Ordinary source jets through Z4 are propagated with binomial product and
reciprocal rules. Core Phi/Phi_rho and U/U_rho bounds feed the prescribed
comparison integral; the complete inlet jets and positive Duhamel kernels
bound JH, JK and M. Their exact physical-y ODEs give the mixed jets.
Dbar_Z3 consumes fourth source jets because its defining expression already
contains JH_Z, JK_Z and M_Z. The comparison is never replaced by core data
after its cutoff changes, and a missing physical-y jet fails the build.

The following are rounded displays of the serialized directed upper caps from the inner
`actual_Ea_mixed_bounds` packet. Exact lower/upper mpf tuples remain in the
report and are authoritative; displayed entries describe norm caps, never signed
function values.

| serialized norm cap | Z0 | Z1 | Z2 | Z3 |
|---|---:|---:|---:|---:|
| E_Z/E | 1.0 | 5.37250623360879124912809e+408906090034569567 | 2.88638232301653198504437e+817812180069139135 | 1.55071070229845417515852e+1226718270103708703 |
| E_yZ/E | 11.1566361934974848740407 | 5.99390974956706940399137e+408906090034569568 | 3.22023174932375892280872e+817812180069139136 | 1.7300715146906837256932e+1226718270103708704 |
| a_Z/a | 1.0 | 1.202360937163858586777e+33 | 3.56534353445527495492644e+65 | 1.58583521106322095795331e+98 |

The implementation is [the source-jet provider](../src/openai_ns_reconstruction/original_micro_Ea_mixed_bounds.py),
[ordinary-jet arithmetic](../src/openai_ns_reconstruction/original_ordinary_jet_bounds.py),
and [the report runner](../scripts/run_original_micro_Ea_mixed_bounds.py).
For continuation, inspect the saved receipt only:

```
python scripts/run_original_micro_Ea_mixed_bounds.py --inspect-report
```

The accepted writer has already run and should not be replayed.

## Exact scope and open work

The current scope is the actual background E/a derivative provider only:
ordinary Z0..3, physical y0..1, s in [0,2], and Z in [.370,.372]. The narrow
saved .371 packet is not promoted to a whole-Z result. The source ports are
directly provenance-bound, while the comparison second chart remains its own
history and is not identified with the core.

The following remain open: actual b, t0, q and p2 bounds; Cy_Z3; mean-Z3 and
QZZ; finite-N modified E_N derivatives; the terminal/global closure; and any
full five-moment or corrected-field claim. The next source step is a whole-Z
log F0 lower bound followed by pressure jets and b/p2 control. The remaining
phase and mean source terms must retain their actual graph ports and physical
scales.

Production exits 0 in 4.007 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_Ea_mixed_bounds.json` is
16,577,711 bytes with SHA256 `a48e8ba9de2f2865598d78b6ba7cfcc94f95ec4b94cdda45b7acb2c168413576`, 1999 input pins and
358 implementation pins. The report's independent review
and saved-source checks are retained here as receipts only; separate runtime
model metadata is unavailable and is not claimed.

The milestone is partial progress toward the active goal. It makes no claim of
full mixed-root closure beyond this partial E/a provider, C0 derivative
closure, full first/second-micro construction, terminal matching, or global
incompressible-NS completion. The old 212 forward terms, six signed sources,
and inherited derivative-error packets remain unchanged.
