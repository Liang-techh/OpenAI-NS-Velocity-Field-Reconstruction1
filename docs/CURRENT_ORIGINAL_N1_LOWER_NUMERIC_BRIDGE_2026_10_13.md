# True n1 lower-collar pressure C1 numerical bridge

The true signed pressure functions on the entire original lower collar
`Rin=sc/4` to `rminus=sc/2` now have a same-source directed numerical
interface. Their defining objects remain the accepted true C1 integral
FunctionRefs. The numerical prefix supplies enclosures of those functions
and their ordinary first Z derivatives with the nonlinear core, own mean
history, native pressure datum and all nonzero collar errors retained.

The full original goal remains **ACTIVE / INCOMPLETE**. Complete pre-Rp and
pulse pressure transfer, quantitative reserved coefficients, repaired own
fields, five terminal zeros, heat/stress, higher-order recursion,
oscillatory cancellation and corrected NS validation remain open.

## Saved interface and recovery

The quartet under `experiments/root_st073` is
`lei_ren_part1_paper_compliant_current_original_n1_lower_numeric_bridge`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gate is
`current_original_true_lower_collar_pressure_C1_numeric_bridge_installed`.

```powershell
python scripts/run_original_n1_lower_numeric_bridge.py
python scripts/run_original_n1_lower_numeric_bridge.py --evaluate --Z .371 --output tmp/n1_lower_pressure.json
python scripts/run_original_n1_lower_numeric_bridge.py --evaluate --Z .371 --include-inner --output tmp/n1_pressure_to_rminus.json
```

The first command checks and restores saved results without executing
ancestor mathematical suites. `--evaluate` restores numerical owners
once and executes a fresh original query. For repeated queries, retain
`CurrentOriginalN1LowerNumericBridge` and its accepted owners and use
`evaluate()/report()`. A live same-Z true Rin packet can be passed as
`Rin_packet` to compose the inner and lower prefix without re-evaluating
the inner solution.

Saved cases use Z=.371, independent Z=.211 and the continuous interval
Z=[-1,1]. Each contains four signed pressure rows: m13 and m15, each with
ordinary Z0 and Z1. The Z=.371 case also composes all five actual Rin
moments through rminus. The three velocity moments already cover their
full support; both pressure moments still require the remaining exterior.

| Query | Passing pointwise value/Z1 rows | Largest relative lower-pressure interval width |
| --- | --- | --- |
| Z=.371 | 4/4 | 9.53021667222e-29 |
| Z=.211 | 4/4 | 3.39924911599e-29 |
| Z=[-1,1] | Range enclosure; relative target inapplicable | Crosses zeros |

Whole-Z range widths include genuine axial variation and zero crossings.
They are not pointwise relative errors, terminal-zero errors, or an NS
residual. No uniform relative accuracy at zeros is claimed.

## Exact function and numerical containment

The physical path and source equations are

```text
rho(s) = 4 exp(hb s), R = epsilon rho
chi = 1 + (hb-1) sigma(s)
ratio(0) = 1
d_rho log(ratio) = (chi-1) Phi_core,rho / Phi_core
Psi_rho = ratio chi Psi_core,rho
rho meanPsi = 4 meanPsi_core(4) + integral_4^rho Psi_actual drho
U = Uaxis + epsilon Psi
m = Uaxis + epsilon meanPsi
Q = [2ZU-(1-delta)Zm-(1-Z^2)m_Z] / (1-delta Z^2)
```

This is the leading source's `(1-delta)` recovery equation. The n1
recovered field elsewhere continues to use `(1+delta)` as required.
Incoming histories start at the same true core rho=4; the comparison
core and its interval center do not replace the actual collar.

The direct pressure source is

```text
Q_y = rho Q_rho
Q_yy = rho Q_rho + rho^2 Q_rhorho
L = 1-delta Z^2
S = [Q+(1-delta)Z Q_Z/2+Q_y]/L
    + Q(Q/2+Q_y)
    + U[-2ZQ+(1-Z^2)Q_Z-2ZQ_y]/L
N1p = -S/2 + (Q_yy+Q_y)/R
m13_lower = integral_sc/4^sc/2 N1p epsilon hb rho ds
m15_lower = -integral_sc/4^sc/2 N1p epsilon^2 hb rho^2 ds
```

Both viscous radial derivatives, the positive hb measure and the single
m15 minus sign remain. Exact production-AST checks compare the true
differential N1p with the legacy forcing formula, including its ordinary
first Z derivative, for arbitrary analytic U and Q. They also compare
all nine leading U/m/Q radial rows and four source-correlated collar
difference equations.

The numerical prefix uses the actual core plus complex comparison and
Cauchy derivative errors. Sigma and sigma-prime have nonzero bounds
factored with two powers of the original hb. Ratio, transported Psi and
own mean retain their own additional hb powers and errors. No error or
source parameter is set to zero.

For range integration, a whole-cell enclosure of the genuine density is
multiplied by the positive phase length and physical measure. This is
**containment of the true integral**, not equality of a range product
with its defining FunctionRef integral. All physical terms stay in the
original `(hb,F0_anchor,Pstar_squared,Lambda)` log frame.

## Source identity and scope guards

The accepted common-inner bridge binds the unique twenty-term nonlinear
core, directed model-plus-error callbacks, native P0, the unique F0
anchor and original amplitude recurrence. This bridge additionally binds:

- the same live core/mean/source-frame owners;
- exact sparse binary sc/4 and sc/2 endpoints;
- hb from the same native global norm and `source_log_hb_enclosure`;
- original stable-sigma and ordinary sigma-jets provider ASTs;
- the actual ratio/Psi/mean rho=4 histories;
- exact absence of finite-N source modifications throughout sc/2;
- genuine signed C1 integral nodes with Z-independent endpoints;
- same-Z same-frame Rin packets for optional prefix composition.

Wrong-domain, foreign, wrong-Z and changed-binding packets are rejected.
The output does not supply full external pressure C1 bounds or promote
the repaired-family and recursion gates.

## Detailed next actions

- [x] Bind lower physical source and exact original endpoints to the true graph.
- [x] Compare true and numerical pressure differential recipes and ordinary Z1.
- [x] Retain actual core/ratio/Psi/mean histories and all Cauchy comparison errors.
- [x] Execute two original points and the full real-Z interval enclosure.
- [x] Compose a same-source actual Rin packet through rminus for all five moments.
- [x] Bind pre-Rp active/compact/quiet histories and derivative enclosures to the true C1 functions using live exact function replay and all 228 original primitive handles. See [the pre-Rp numerical bridge](CURRENT_ORIGINAL_N1_PRERP_NUMERIC_BRIDGE_2026_10_13.md).
- [x] For active pre-Rp source, bind all 17 true charts, 68 numerical chart views and 228 native cells; compare the seven S0/nine S1 physical terms and both finite endpoint derivatives. At O3 retain q_y=-q and q_yZ=-q_Z.
- [x] For compact pre-Rp source, bind the actual repair-band C3 rows, selected N and repair integer, Pstar applied once, full incoming histories and `recover_Q_increment`; use the accepted `normalized_first_Z_upper` only as a symmetric derivative enclosure.
- [x] For quiet pre-Rp source, compare the exact K=Rp Pstar kappa recipe and signed value formulas; preserve their exact Z1=0 without zeroing other fields.
- [x] Convert pre-Rp magnitude bounds with `preRp_pressure_source.Rc_unit_log`, canceling common Rc first. Use physical power 1 for m13 and power 2 for m15. Magnitude caps do not supply signed centers.
- [ ] Bind pulse main/gap/end rows, quadratic and viscous errors, selected finite N, incoming mean and energy histories to the true pulse FunctionRefs.
- [ ] Compose lower/pre-Rp/pulse contributions exactly once; apply beyond-Rv zero support only to the pressure forcing integrals it covers.
- [ ] Enclose all five complete true moment functions and their genuine first Z derivatives in common physical units.
- [ ] Transfer those ranges through the actual physical reserved inverse, retaining Astar_Z, nonzero mu and original radius factors.
- [ ] Extend required true coefficient/pressure C3/C5 derivative bounds.
- [ ] Build repaired own M/K/Q/V/P using the actual compact patches and n1 recovery equation.
- [ ] Prove all five actual terminal zeros and joins before closing n1 repair.
- [ ] Continue the full heat, admissible stress, flat remainder, n>=2 recursion and oscillatory correction objective.

Mark an item complete only with a checked source/receipt or a verified
implementation. Preserve saved owners and packets instead of rebuilding
accepted ancestors after a stream observation timeout.
