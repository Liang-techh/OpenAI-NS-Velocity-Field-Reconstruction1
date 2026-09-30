# Step 3 finite exit switches

`lei_ren_part1_paper_exit_switches.py` implements the finite Section 9.4
Step 3 continuation from the actual `ExitContinuation` endpoint at `R = 100`
to `R = 110`. It keeps the source switch locations exact in
`x = log(R / 100)`:

* `0 <= x <= h_b`: `a = epsilon * Dbar` and
  `b = -epsilon * Ebar * (1 - sigma(x / h_b))`;
* `h_b <= x <= 2 h_b`: `a` is the flat-step blend from
  `epsilon * Dbar` to `4 / 5`, with `b = 0`;
* `2 h_b <= x <= log(1.1)`: `a = 4 / 5` and `b = 0`.

The class `ExitSwitches` exposes `evaluate(y, Z)` and `evaluate_R(R, Z)`.
Its state carries `g = log(F / F_100)`, `U^z`, the five cumulative moments,
and their analytic `Z` tangents. The first two stages use fourth-order
centered arbitrary-precision driver tangents. The final constant-power
stage uses exact exponential primitives for all five moments. The starting
quadratic axial and swirl integrals must be supplied by the provider as
`raw_quadratic_integrals`; omitting the collar is treated as an error.

Run the focused receipt with:

```powershell
$env:PYTHONPATH = 'src;experiments/root_st073'
python experiments/root_st073/lei_ren_part1_paper_exit_switches.py
```

The generated JSON compares 16 and 32 switch steps while holding the
provider tangent at 32 steps, checks all five moments and their `Z` jets at
`R = 100`, and samples both switch boundaries plus `R = 110`. This is
finite numerical evidence for the candidate only. It does not certify the
source parameter gate, the relaxed cone globally, pressure closure, moment
repair, PDE, finite energy, or the paper's scale recursion.

Recorded result at Z=.3: the maximum R=100 field/moment/jet/slope
matching relative error is 1.60e-160. All six switch/endpoint samples
pass the relaxed cone. The 16/32-step maximum relative difference at
R=110 is 7.03e-12; the largest interior difference is 2.44e-9. At R=110,
a=.8,b=0. Recovering D from the recorded zero-b relaxed margin D-2 gives
D approximately 5.68e37, above the source endpoint threshold three.
This is an algebraic reading of the same cone diagnostic, not an
independent extra validation.

The current fixture fails necessary source parameter gates recorded in
connection_scale_gate.md/json. Rebuild the shared outer pressure and
core for a replacement before using the long angular reshape. No old
receipt should be reused after that change.
