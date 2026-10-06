# Whole current gap and reciprocal gap-end signed cone

Current successor: [CURRENT_FLATTEN_POWER_CONE_2026_10_06.md](CURRENT_FLATTEN_POWER_CONE_2026_10_06.md), commit [caa9f74f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/caa9f74f851391267adb0ffaacfdf98f30e5a8d0). F57C-cone1b-flatten is now complete; the next work is current angular/tail cones.

Implementation and scoped whole current gap/gap-end cone receipt: commit [b592adf4](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b592adf478cfa4adc63011a30f3c80693e91ad3f).

## Constructed result

The current pulse gap xi[11,12] and reciprocal gap-end phase[0,1], all Z[-1,1], now satisfy the original strict two-vector cone. This closes **F57C-cone1b-gap**. Together with the earlier main, exit and pulse-end proofs, five consecutive current pulse regions have whole-domain cone bounds on the same admitted velocity/pressure/history graph.

This is not a completion of the global admissible tensor, actual oscillatory correction or coefficient recursion. The remaining 28 registry regions, global completed-tensor condition, edge weights, actual homogeneous pulses/covariance, genuine n-dependent recovery, flat completed remainder, corrected NS, prescribed-domain energy and measured dynamics remain separate construction requirements.

## Source formulas and signs

The same ordinary radial source coordinate is `s=logR`, with geometric distance d=13-xi and s=-d/mu. Gap uses d[1,2]; reciprocal gap-end uses `d=4mu+(1-4mu)*(1-phase)`, d[4mu,1]. Exact phase0 is d=1,s=-1/mu; phase1 is d=4mu,s=-4. These are exact source identities, not recovery by dividing independently rounded log intervals.

Let finite=-3*saddle_L+2*log(saddle_u0)-log(6)/2-2*log(mu). Keep each original moment scale:

`logD0=finite-1/mu`,
`logD1=finite+(d/2-1)/mu-d`,
`logD2=finite+(d/2-1)/mu-2d`.

D0/D1/D2 are positive source factors from the same shared pulse map. They are distinct from distance d; none is replaced by the main/exit local D=1 convention. The radial velocity and meridional remainder use D1, while signed stress sectors retain their own selected D recipe. All source coefficients remain directed enclosures of the original functions.

The original gap Bh and axial velocity jets are structurally zero. Thus actual shear is a=2+2mu, b=0 and kappa-2=2mu, retained without rounded subtraction. This does not zero the radial velocity, five cumulative histories, absolute pressure, or the four axial stress sectors.

Normalize by the positive common sqrt(R/2)*B. The actual signed theta equilibrium has lower bound g/2, g=(mu-delta/2)/(1-mu)>0. Each of four theta corrections is bounded by g/64, giving the conservative total theta lower g/4. Each of four axial corrections is bounded by `(g/4)*exp(-1000)/4`. Therefore the full signed stress obeys

`2*theta^2-2mu*axial^2 > 0`,

with normalized lower bound `1-mu*exp(-2000)>0` after division by `2*(g/4)^2`. The stored producer retains the actual tiny positive theta and quadratic margins; ordinary double precision underflow is not interpreted as zero stress.

All nine signed sectors per region remain: theta equilibrium plus four theta corrections, and four axial histories. Absolute coefficient bounds apply only to errors around the signed equilibrium. Full diagonal, divergence, Cartesian tensor, radial velocity and remainder records remain attached.

Pressure memory has exact grouped log

`logPstar+logU+(-13-d)/(2mu)-13-d`.

Its signed distance derivative is `-(1/(2mu)+1)`, so its maximum is at d=4mu. All other retained correction logs increase with d and reach their uniform maxima at d=2. The producer stores signed derivatives separately from positive monotonicity magnitudes. The current Rv pressure datum and raw native forward pressure are both retained; they are different functions.

## Evidence and API

Use the `experiments/root_st073/lei_ren_part1_paper_compliant_` prefix:

- `current_pulse_gap_cone_operator.py`: 36 original generic signed-factor/shear/cone identities, exactly three allowlisted historical receipt statements removed, 17 live current source-log/endpoint identities, and whole-current bounds. Mathematical statements and AST checks remain unchanged; old admissions are not loaded or promoted.
- `current_pulse_gap_cone.py`: producer and `CurrentPulseGapCone(maincone=checked_current_main_exit_cone)`. Both complete whole views are obtained directly from its checked registry. `native(...)` scopes gap and gap-end flags separately; other regions retain nested inherited admissions.
- `current_pulse_gap_cone.json` / `_check.json`: exact source composition, 33 positive directed bounds and checked input hashes. The checker rejects foreign source, wrong D2 selection, nonzero axial source and omitted pressure memory.
- `current_pulse_gap_cone_views.json.gz`: deterministic gzip containing both full unpruned signed whole views. Read with `json.loads(gzip.decompress(path.read_bytes()))`; its hash is bound by the receipt.

Working/index source hashes passed for all 787 dependency files. The focused controller and checked per-region API passed after the pressure-rate correction.

Focused controller stage: `currentgapcone`. Reuse the checked warm graph. The unchanged source inventory is 33 regions / 32 adjacent / 14 internal tensor traces and primitive atlas 14 / 8. Existing current complete source theorems supply exit-gap, gap-coordinate and gap-end attachments, including the same full beta-square defining integrals and cell-partition additivity. This cone adds sign bounds to that same graph; it does not recompute or replace the matching data.

Read-only review uses GPT-5.6 Luna / max. Its pressure-rate metadata finding was corrected and the producer/check/controller/scoped API rerun after that change. No remaining material issue was found in the scoped review.

## Next executable tasks

- [x] **F57C-cone1b-gap:** original D0/D1/D2 and Q recipes, actual zero axial shear, all signed histories, both whole domains and inherited three current attachments.
- [ ] **F57C-cone1b-flatten:** add same-source shear/cone adapters for `flatten` offset[0,100] and `outer_power` phase[0,1]. Use owner `flatten`, its actual velocity jets and source amplitude/radius recipes. Keep complete forward pressure and selected energy. Bind the pulse-end/flatten join and flatten/power join. Recompute current whole-domain sign margins; historical cone checks alone cannot admit the current graph.
- [ ] **F57C-cone1b-angular:** `outer_angular` s[-4,0], all Z[-1,1], with all four support traces. Bind current selected angular controls, full raw moment/absolute pressure history and ordinary s derivatives. Preserve both correction supports, positive/negative signed contributions and actual source shear. Add uniform edge directions where stress vanishes.
- [ ] **F57C-cone1b-tail:** current `steep_entry`, `steep_power`, `steep_exit`, `waiting`, `heat_collar` and `heat_exterior`. Consume exact Gamma/higher-pressure identities and current source maps. Prove regional strict margins where stress is nonzero; establish the stress-free exterior and its limiting edge behavior separately. Retain positive constant-viscosity physical scaling and full Cartesian remainder.
- [ ] **F57C-cone1b-entrance:** `pulse_entrance` xi[0,.02], `O3_power`, `O3_slope_mu`, and O2/reference regions. Include actual incoming moments and raw pressure in the normalization; bind variable ordinary log-amplitude slopes. Keep source-proven zero/nonzero terms. Publish whole domains and source attachments, distinguishing failed signs from inconclusive enclosures.
- [ ] **F57C-cone1b-inner:** remaining actual patch/restore/reshape/switch/bridge/core. If the current inner profile only meets the relaxed cone, construct the original periodic shear loop, finite uniform N and independent five-moment repair on reserved intervals. Preserve the analytic collar and absolute pressure. Any changed source family requires dependent tensors, joins and receipts to be rebuilt.
- [ ] **F57C-cone1c:** actual two homogeneous pulse families, covariance integrals with cosine factor1/2 and torus area, finite column errors and determinant, uniformly positive squared amplitudes, flat edge weights/smooth square-root extension, fixed-base signed linear lift, mean correction and full averaged quadratic cancellation. Reference matrix algebra is not the actual pulse construction.
- [ ] **F57C-recursion:** actual n=1 and distinct n>=2 recovery and moment repair; cancel/absorb the nonflat leading origin E; finite-order estimates, smooth sum, resolved physical fields, corrected Cartesian NS residual, prescribed-domain energy and scale-recursion/material-line diagnostics.

The full long-term goal remains active. Do not mark it complete at this regional cone milestone.
