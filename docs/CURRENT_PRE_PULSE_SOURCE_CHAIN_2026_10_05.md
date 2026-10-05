# Current Rh-to-Rp source chain — 2026-10-05

CurrentPrePulseSourceDispatcher extends the accepted current Rh owner through the original O2_slope, O2_axial, O2_buffer, O3_slope_mu and O3_power charts. Fourteen downstream chart owners now share the current construction family. Every pre-pulse chart uses the same CompliantPrePulseMixedC4 object already admitted at Rh.

Use branch codex/st073-transition-next:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentprepulse
```

```python
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import CurrentPrePulseSourceDispatcher
field = CurrentPrePulseSourceDispatcher()
packet = field.evaluate("O3_power", Z=".3", coordinate=1)
```

The API returns source-field enclosures. O2_axial's phase=log(y)/Md and O3_power's Tw phase select coverage; their returned derivatives remain ordinary logR/Z derivatives.

## One current history and pressure owner

The loader first consumes actual_Rh_source_join_check.json, including its compliant .001 analytic pressure function/projection and current unique implicit patch closure. It reuses that exact Rh reference provider for all five subsequent charts; it does not construct independent nominally matching outer fields.

The original source methods are unchanged. New AST bindings require each actual parent call, inlet cache, five-history dictionary, kernel arguments, ordinary radial log-amplitude jets and packet velocity arguments. The common packet retains the same normalized_jets pressure callable and physical_mixed operator.

The existing independent closed-integral and cutoff-coordinate fixtures are consumed by current hashes. Existing universal comparisons of original source equations are reused for arbitrary input histories; their historical Rh admission is not reused as current evidence. Current Rh admission comes from the new dedicated receipt.

## Functional interfaces and retained tails

source_histories interprets the bound original hist AST in arbitrary exact source symbols. Neutral integrating factors/kernel values establish25 five-history identities at Rref, O2 slope exit, turnoff-to-buffer, Rd and Rw. At turnoff/buffer the same slope parent, y=exp(Md) and turnoff-kernel arguments are used.

Original sigma endpoint derivatives1..4 vanish. The two sides retain the same log-amplitude/axial-velocity source jets, parameters, primitive initial values and five ODEs:

```text
m_y = V-m
h_y = u-3h/2
k_y = uV-3k/2
e_y = V²/Pstar²-u²/2-e
p_y = u²/2
```

The common physical_mixed callable differentiates physical radial prefactors before producing grids. These source identities imply675 mixed4 interface rows across five interfaces.

After V vanishes, m and k continue as their nonzero decaying histories. No history is reset to zero. The O3 power source keeps its exact positive mu and slope -1/2-mu even when ordinary enclosures cannot resolve mu relative to1/2. Tw remains the original positive logarithmic source length.

## Acceptance and remaining interface

Focused acceptance passed:345 current dependencies,14 downstream owners,25 exact history source identities,675 implied interface mixed rows and675 whole-domain mixed rows across the five new charts. A fresh lazy O3_power call loaded the new receipt and confirmed identity of its provider with Rh_reference.

The producer records source_chain_proved=True and source_chain_certified=False. Runtime certification becomes true only after consuming current_pre_pulse_source_dispatcher_check.json. Legacy dispatchers/receipts stay unchanged.

The Rp-to-pulse external join remains false. The source scan located a concrete remaining bridge: pulse inlet data currently pass through SharedOuterBuffer.power and saved fifth-order/canonical constants. Their absolute-history to canonical-unit function definitions must be source-bound to this current chain, including common P0/Pstar/mu/Tw and the actual Rp radius relation. Interval overlap at Rp is insufficient.

Current physical/Cartesian composition, complete shared leading inputs/remainders, production point fields, completed global tensor/flatness/required-domain energy, true n-dependent recursion and oscillatory/corrected dynamics remain open.
