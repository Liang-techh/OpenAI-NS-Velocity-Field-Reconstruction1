# Corrected own-swirl C1 function source and cached point provider

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. The corrected complete
m15 source is now propagated through the existing five-moment coefficient
formulas into a new append-only C1 field-source provider. It also exposes
the previously accepted corrected reserved-response enclosure at exact
Z=371/1000, full-support x=2, without repeating its computation. The legacy
default forward chain is not reconfigured. The current Rm-to-Rp physical
germ join and functional terminal matching remain open.

## New implementation and actual result

The new module is
`src/openai_ns_reconstruction/original_corrected_own_swirl_provider.py`.
Its `rebind` API consumes the saved 53,500-node graph, old coefficient,
response-matrix, residual and repair-velocity roots. It replaces only the
old physical full-m15 C1 pair53175/53176 with corrected53474/53475. The
immutable graph prefix is reassembled from authenticated saved rows and
accepted deltas; no ancestor graph constructor or numerical provider runs.

Only declared arithmetic operands are traversed: sum/product arguments,
negative argument, positive-quotient numerator/denominator and analytic
unary argument. Literal integers, derivative orders and other metadata
are preserved. Replacement ports stop traversal. Only changed downstream
rows are copied after the complete prefix. Cycles, forward references,
unsupported visited operations and incorrect source digests are rejected.

The new dependency propagation ran ONCE. It appends80 rows and yields a
53,580-node graph with canonical SHA256
232138925df3c1e8cceca3edd591d88e06e7da2a6c01b137c5e51a3776be30ff.
Its ordinary-C1 coefficient and field roots are:

```json
{
  "coefficients": {
    "alpha1": [
      53228,
      53250
    ],
    "alpha2": [
      53233,
      53255
    ],
    "beta1": [
      53503,
      53511
    ],
    "beta2": [
      53513,
      53515
    ],
    "beta3": [
      53517,
      53519
    ]
  },
  "repair_velocity": {
    "Uz1": [
      53379,
      53380
    ],
    "Utheta1": [
      53572,
      53577
    ],
    "F1": [
      53578,
      53579
    ],
    "physical_radius": 53392,
    "log_R_over_Rc": 24978
  }
}
```

All three beta coefficients change. Alpha1/alpha2, Uz1, the physical
radius, log(R/Rc) and every response-matrix port remain the identical old
roots. The m12/m13/m15 linear-response residual expressions are rebound
with the new beta and m15 operands. Their new graph roots are expressions,
not newly evaluated zero residuals or fresh numerical certificates.

The established field formulas are retained:

    Uz1 = sum(alpha_i * axial_bump_i)
    Utheta1 = sum(beta_i * angular_bump_i)
    R = exp(logRc+s)
    F1 = Utheta1 / sqrt(2*R).

Thus the m15 change is propagated through the angular coefficient system,
not inserted as an unrelated terminal number. The other four physical
moment source pairs remain unchanged. The preRp m15 segment is separately
bound to53472/53473 from the accepted corrected packet; it is never inferred
from replacement of the full moment. No repaired own M/V/P C3 graph or
Cartesian velocity evaluation is claimed by these C1 source bindings.

## Callable data-only provider and one reserved source slot

`CorrectedOwnSwirlProvider(binding_packet, cached_point_packet)` exposes
deep-copy `field_source_bindings()` and `coefficient_source_bindings()`.
Its `evaluate(Fraction(371,1000),2)` returns the accepted parent directed
binary tuples in the response packet's `coefficient_C0_C1` shape. Other
queries are rejected. This method performs no arithmetic on the cached
enclosure and does not instantiate the old response or forward provider.
Numerical Utheta1/F1 and a whole-Z numerical coefficient provider remain
unimplemented. This class is a new source/data view, not a drop-in instance
of the old provider's context/interval/evaluation chain.

The full-support response identity is retained through ordinary orders0,1:

    response(x=2) = [0,0,1]
    reserved_corrected = -full_corrected
    target = -complete_pulse
    reserved_corrected-target = -(inner+lower+preRp_corrected).

The response uses the original physical measure -R*theta0*Utheta1 ds.
Physical reserved roots53476/53477 and normalized roots53492/53493 refer
to that same source. The independent negative-pulse target retains roots
53494/53495. U=reference*A*Rv^2 remains positive and Z-independent; ordinary
Z1 divides by the same unit53427. P1 is not m15; pressure Q is not loop q.

The canonical saved ledger is `macro_endpoint_m15_numeric.json`, with213
raw terms and183 reduced terms. Its raw slot212 and reduced slot182 are
the same `reserved_n1_own_swirl_once` record in all three saved Z rows,
canonical term SHA256
36a0581ecc4e55fea62f548bdd8859a1813c1ea1793487d5f8d700fbcde98fc7.
The older base repaired-Rp report's201-entry reduced list is not this183
entry source namespace.

The new overlay replaces that existing source leaf ONCE. It does not add a
214th raw term. The180-entry joint source view retains179 old source
references; only its reduced index182 is overridden. All212 forward terms,
their exact raw prefixes, six structural channels (five signed contributors
plus core-swirl error-only), all finite-N/width errors and prior endpoint
reductions remain referenced unchanged. Saved raw/reduced digests and
source-record locators are recorded for each Z. Old terminal transport
handles remain provenance in their separate namespace and are never
relabeled as corrected terminal roots. No ledger is numerically resummed.

## Validation and limits

Each actual component was checkpointed before subsequent inspection:
source input,80-row rebinding, then data-only provider/overlay. Recorded
component seconds are {"data_only_provider_installation": 0.001127199997426942, "new_C1_dependency_rebinding": 0.2947610999981407}.
Production command wall time was4.5507986s including authentication.
The saved-only inspector passed: new graph hash/count, unchanged parent
packets, source pins, exact cached tuples and the one-slot override.
It does not rerun the rebinder, source producers, recurrences, quadratures,
norms, old checks or completed fixtures. New Python passed LF/single-EOF/AST
preflight before production.

The worker's new moderate toy DAG exercises independent before/after
arithmetic, metadata preservation, guards and cache isolation. Its observed
inline harness stdout reported `fixture_harness_passed=True`; no raw harness
receipt was saved. This delegated observation is recorded as such, not an
independently retained execution log. The worker summary quoted an incorrect
fixture hash; root preflight stopped before production, pinned the actual
6,806-byte artifact SHA256
91b61cc485259c6c100ad6b9b5681628e0be194dc9f16a68f57c08126953a448,
and copied those exact bytes. The original artifact and its completed
fixture were not modified or replayed.

All three child tasks used fixed `luna_max_worker` with `fork_turns:none`.
The tool returned role/task metadata, but no independent runtime model
attestation. Read-only provider and response-identity audits are separate
from the worker's bounded new module/fixture work.

The accepted parent point enclosures and their inherited oracle/error
premises are unchanged. A tiny single-point discrepancy enclosure is not
an all-Z zero function. No cap, midpoint, bound overlap or unproved endpoint
cancellation is used to declare closure.

NEXT establish a current physical-germ source join from Rm into the
patch/transition/compact/quiet chain on a common Z domain. Current absolute
Rm data covers[0.37,0.372]; old cumulative tables use Z=-0.3,0,0.3. There is
still no saved endpoint identity joining current RmRQy ports183696/183697
to the old lower endpoint. Shared source family or N does not supply that
identity. Use new source consumers and genuinely new evaluations; do not
replay this turn's rebinding, accepted parent repair, or any ancestor.
Then propagate the new provider through the admitted current forward chain.
Independent five-moment functional closure, exact heat exterior,
admissible stress/flat remainder, true n-dependent recursion, realizable
mean/oscillatory corrections and Cartesian NS validation remain incomplete.

Report experiments/root_st073/lei_ren_part1_paper_compliant_current_original_corrected_own_swirl_provider.json: 2297305 bytes,
SHA256 bf4ab120edd8dc02925977dc37614215a2be00a83b8393d1d5763c9ad50cba42.
