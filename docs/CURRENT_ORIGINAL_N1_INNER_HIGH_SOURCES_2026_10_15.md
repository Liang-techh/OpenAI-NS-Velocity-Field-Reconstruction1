# Original n1 inner high source functions — 2026-10-15

## Delivered

The same regular six-component Volterra solution now supplies positive-radius
inner-core, inner-first and lower-collar state sources through physical y4,
with P1 and all five cumulative moments through y3. Ordinary Z is retained
through 3, except Q1/V1 through 2. The 414 previously installed lower scalar
rows are unchanged. Complete nonzero mass, pressure and signed histories are
retained; no new solution, model, enclosure center or common-N value is chosen.

```python
from openai_ns_reconstruction.original_n1_source_dispatcher import OriginalN1SourceDispatcher
source = OriginalN1SourceDispatcher().extend_inner_high_sources()
query = source.query('inner_first', 'F1', fraction='5/8',
                     physical_y_order=4, ordinary_Z_order=3)
pressure = source.query('lower', 'P1', fraction='1/2',
                        physical_y_order=3, ordinary_Z_order=3)
```

`extend_inner_high_sources()` restores a checked 47986-node suffix of the
accepted 2232772-node graph. The combined graph has **2280758 nodes**. The
small suffix is stored in a compressed archive (about 0.4 MB); the original
large archive is reused. No live historical owner or process is needed.
The default base dispatcher remains available with its weaker inner contract.
Call the extension before issuing queries, since it preserves the exact
archived graph numbering. This is source-only recovery, not a full owner or
numerical solver restore.

## Mathematical source construction

The source equation is the original raw six-row equation

`d_x W + D W/x = B0 W + B1 dZ W + g`,

where `x=sqrt(rho)`, `D=(0,0,2,0,3,1)` and
`W=(F1,Uz1,K1,P1,d_x_F1,d_x_Uz1)`. Hence

`d_y W = x/2*(B0 W+B1 dZ W+g)-D W/2`.

Every original B0/B1/g entry is used. On the core, `partial_rho=d_y/rho`;
on the first chart, `partial_s=hb*d_y`. The derivative engine expands these
actual equations. It never adds a `radial_order` field to a regular n1
component. It applies exact product, quotient, capture-aware composition and
full moving-endpoint Leibniz rules to the saved source expressions.

The actual finite dependency graph uses **regular Z through 5**, and the
same admitted nonlinear analytic core through **rho4 / Z6**. The Z6 core
leaves are new source-function derivative bindings in this milestone, not
previously evaluated data or a previously certified C6 bound. Their original
analytic function definitions, source family, twenty-term operator, Banach
solution and convergence/Green providers remain unchanged. Finite guards
reject higher requests, changed definitions and invented regular radial tags.
The sparse kinematic/B1 matrix determines the dependencies; a generic
radial-order-plus-Z upper estimate is not the actual dependency graph.

The cutoff keeps its actual argument and original sigma provider for that
argument. New sigma3/sigma4 derivatives use the original `_sigma_left` pin;
interior rational formulas are evaluated only on their active interval.
Flat branches retain exact zero. All active-cutoff products and both signed
pressure/moment densities are differentiated, rather than extending the
pre-existing flat-collar formula over the active region.

## Focused evidence

```powershell
python scripts/run_original_n1_inner_high_sources.py --check --write-report
python scripts/run_original_n1_inner_high_sources.py --restore-source
```

The focused construction/check run exited 0 in **12.750 s**:

- 414 unchanged lower scalar rows;
- 12 actual six-row equation identities with original matrix entries;
- 16 comparisons between differentiated existing first rows and retained
  original second rows, on both core and first chart;
- five independent partial/Leibniz/capture fixtures;
- four original logistic sigma derivative identities;
- 2214 mixed source queries across the three complete native intervals;
- 104 regular radial requests expanded through the actual equation;
- three rejected unadmitted source requests;
- lossless compressed archive roundtrip.

A separate fresh process restored the source archive in **8.890 s**, exited 0
and issued **84** high-derivative source queries without constructing a
Volterra/derivative owner. `full_owner_restored=false` and
`numerical_oracle_available=false` remain explicit.

The base regional dispatcher was also checked after its small API change:
6318 scalar slots / 18954 source queries, 18 invalid requests and five
operator identities pass in 9.422 s. No ancestor mathematical suite was rerun.
These are finite source/operator and archive checks, not numerical NS,
quantitative analytic bounds, all-order smoothness or blow-up evidence.

## Remaining work: task 13 is still IN PROGRESS

1. Compare actual new higher germs at core rho4 -> first s0, across the
   inner-flat overlap, Rin -> lower, and rminus -> first_micro. Preserve the
   common analytic solution and incoming histories; do not infer a new global
   high continuity gate from regional query tests.
2. Construct regular axis and Cartesian component limits separately. The
   current core query still rejects fraction=0: cancellation of rho factors
   in an interior expression is not an axis proof.
3. Recover physical-coordinate velocity/pressure through the original
   similarity map, and implement certified physical-radius chart selection.
4. Restore original-scale numerical solvers/evaluators and error budgets,
   then physical energy/tail, compatible analytic preheat/exact heat,
   admissible stress/flat remainder and common N.
5. Continue genuine independent n>=2 coefficients and moment repairs,
   smooth summation, oscillatory cancellation and corrected Cartesian NS.

The regional positive-radius derivative slots are available. Higher inner
interface germs, global physical completion, numerical values/bounds,
axis/Cartesian recovery, temporal recursion and complete corrected NS remain
**OPEN**. The full long-term objective stays **ACTIVE / INCOMPLETE**.

## Artifacts

- `src/openai_ns_reconstruction/original_n1_inner_high_sources.py`
- `scripts/run_original_n1_inner_high_sources.py`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_high_source_functions.json.gz`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_high_source_functions_check.json`

The archive records exact prefix/suffix identities, original equation roots,
regular radial equation requests, new analytic derivative nodes, sigma pin,
all original input hashes and implementation hashes. Resume from these files;
the old exec session 84468 remains unavailable and is not a dependency.
