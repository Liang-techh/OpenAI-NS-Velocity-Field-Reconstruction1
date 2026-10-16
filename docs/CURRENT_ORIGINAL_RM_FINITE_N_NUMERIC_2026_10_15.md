Latest continuation: [Physical original finite-N Rm/Rh velocity](CURRENT_ORIGINAL_RM_PATCH_VELOCITY_2026_10_15.md) now assembles the same source into cylindrical/Cartesian leading+n1 value sums. The finite-N numerical source and its limits below are unchanged.

# Actual original finite-N Rm-to-Rh patch values - 2026-10-15

`OriginalRmFiniteNNumeric.evaluate(Z,x)` now continues the actual selected-N value enclosures through the original active patch1<=x=R/Rm<=e. Use x="Rh" for the exact source endpoint. It supplies generic modified E/V and all five own moments through TaylorZ2, unchanged-P0 pressure, and similarity axial/swirl/mass/radial values, with radial recovery through Z1. Original leading controls remain fixed; no finite-N leading-control re-solve is introduced.

The numerical center is the accepted leading point proxy, paired with a complete directed enclosure of the unknown SIGNED finite-N correction. The actual fast phase and signed correction integral centers are not resolved. Magnitude bounds are not selected as field values. The same original DoubleDyadic N, true source family and nonzero positive-order n1 primitive remain.

## One actual correction and one physical measure

The incoming `finite_N_incoming_history_errors` already include1/N. For a native cell, use width=log(x_right/x_left), rates m=1,h=k=3/2,e=1,p=0, and

```
D_out = exp(-rate*width)*D_in + mass*unknown_signed_driver
mass = -expm1(-rate*width)/rate   (rate>0)
mass = width                    (rate=0)
|driver| <= original_N_scaled_density_cap/N
```

Pressure memory is exactly1. New drivers divide by the original N once; inherited corrections are not divided again. The original C3 archive includes every first-order and quadratic finite-N density term and is uniform for all N>=N0. Source B is B_over_Pstar, so local fields are E_N=E*exp(A/N), V_N=V+B_over_Pstar/N in generic S=Pstar units.

The eleven original native cells are1,49/40,5/4,51/40,59/40,3/2,61/40,69/40,7/4,71/40,2,Rh. Coincident Z covers are merged into a single physical cell; they cannot duplicate width or decay the inlet twice. Per-cell driver caps are consumed directly, so no already-applied suffix is applied again. All ordinary derivative bounds are converted once to Taylor coefficients by1/k!.

Exact rational native endpoints and decimal queries determine traversal before interval rounding. The code asserts a contiguous original partition. Queries49/40 +/-1e-140 correctly traverse one/two cells even when rounded interval coordinates straddle the same boundary. A decimal at the uncertain upper e boundary is rejected unless certified interior; exact Rh is always available.

## Dependent recovery and canonical defect units

The generic frame is V=Uz/Pstar, m=(Mz/R)/Pstar. With M=Mz/R=Pstar*m,

```
Q = Pstar*(2Z*V-(1-delta)*Z*m-(1-Z^2)*m_Z)/(1-delta*Z^2)
Ur = sqrt(R/2)*Q
```

The full modified mass derivative is retained. Native y/fast-phase derivatives are not required for these velocity values; they remain necessary for later physical gradients/vorticity.

Canonical partial-defect changes are

```
Delta d0 = x*Pstar*Delta m
Delta d1 = x^(3/2)*(Pstar*Delta k-4Z*Delta h)/am
Delta d2 = x^(3/2)*Delta h/am
Delta d3 = x*(Delta e-8Z*Delta m/Pstar)/am^2
Delta d4 = Delta p/am^2
```

There is no second inner x in the energy cross term: Delta M=Pstar*Delta m, while Delta d0=x*Delta M. The independent raw m/e inversion check covers all five rows and Z0..2 at several x. At Rh the original LEADING five identities remain, while finite-N defects retain unknown signed error intervals; their actual values are not asserted zero or nonzero. The later reserved repair is still required.

Pressure uses the same source-normalized analytic P0/Pstar^2, added once to the modified p row. Source geometry remains atomic: global y=10logCstar+10logPstar-6+log(x). No huge physical radius or phase is materialized.

## Focused receipt

```
python scripts/run_original_Rm_finite_N_numeric.py --write-report
python scripts/run_original_Rm_finite_N_numeric.py --inspect-report
```

The accepted run exits0 in74.861s; saved-report inspection exits0. Seventeen actual packets cover five axial positions and original inlet/interior/Rh points, a new x=1.4 query and the two extremely close cell-boundary queries. Checks include48 independent signed finite-N density/dx-over-x/memory coefficients,45 independently inverted canonical-defect coefficients,2 exact boundary classifications,75 Rm inlet coefficients,126 fixed-control/N/primitive/history checks,104 unique-measure/pressure-memory checks,34 independent generic radial recoveries,51 unchanged-P0 checks and25 terminal nonzero-width error intervals.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_finite_N_numeric.json`,5026004 bytes, SHA256 `c5424697bb9d4043ed872fc5543689fbcd8ed1e9f67b833c8b6cb252478a3e5f`. It pins1762 source inputs and49 implementations. Accepted ancestors are unchanged.

- [x] Original selected-N local E/V and five-history value/error transport through the fixed leading patch.
- [x] Same-source pressure and own-mass radial recovery on the admitted native patch.
- [ ] Fixed-time physical cylindrical/Cartesian patch assembly and outer continuation.
- [ ] Later reserved repair and global modified five-moment closure.

## Remaining full objective

Next assemble physical patch values with the original lambda/time factors and atomic radius geometry, retaining the n1 terminal primitive. Then continue outer/Rh matching and reserved repair. Signed fast-phase centers, total native radial derivatives, higher-Z numerical rows and physical C2 remain open, as do new whole-axis numerical admission, global common N, pressure joins/energy/tails, exact heat, admissible stress/flat remainder, genuine n>=2 recursion/smooth sum, oscillatory cancellation and independent corrected NS. Tasks14/15 remain IN PROGRESS; full goal ACTIVE / INCOMPLETE.
