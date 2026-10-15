# Original n1 regional source dispatcher — 2026-10-15

## Delivered

`src/openai_ns_reconstruction/original_n1_source_dispatcher.py` provides
`OriginalN1SourceDispatcher`. It restores the accepted 2232772-node source
archive directly, checks archive/receipt/input identity and shared graph
prefixes, and installs one regional query API. No ancestor owner is rebuilt
and no historical mathematical suite is replayed.

The 27 views comprise all 17 pre-Rc charts, reserved/quiet/main/gap/end/tail,
inner core, inner first chart, lower collar and an overlapping inner-flat
high-derivative view. These are 27 views, not 27 disjoint physical intervals.
The axis is explicitly excluded from the logarithmic-radial query contract.

Each view exposes Uz1, Utheta1, F1, M1, mean1, K1, Q1, V1, P1 and the five
cumulative moments m11..m15. Incoming pressure and mass histories are retained.
V1 and Q1 remain auxiliary recovery quantities; no Cartesian velocity export
or numerical value is inferred from their names.

## Query

```python
from openai_ns_reconstruction.original_n1_source_dispatcher import OriginalN1SourceDispatcher
source = OriginalN1SourceDispatcher()
source.contract('inner_first')
query = source.query('O2_axial', 'P1', fraction='1/2',
                     physical_y_order=3, ordinary_Z_order=3)
query.as_dict()
```

```powershell
python scripts/run_original_n1_source_dispatcher.py --query O2_axial P1 --fraction 1/2 --dy 3 --dZ 3
python scripts/run_original_n1_source_dispatcher.py --query tail m15 --group cumulative --tail-log-offset 2 --dZ 3
```

Fractions interpolate the complete original **native** interval, including
symbolic scale-dependent endpoints. Tail offsets mean log(R/Rv). Returned
source/coordinate/radius handles use the archived graph plus an append-only
query suffix. The CLI includes that suffix with a query. The manifest pins
the source archives and preserves the small new core-construction suffix;
it does not duplicate or replace the large accepted source graph.

This dispatcher selects already constructed physical y derivative rows; it
does not mistake native derivatives for physical ones or numerically invert
the symbolic radius maps. Source expressions still depend on original axial
variables, parameters and source solvers.

## Derivative limits and unfinished task 13

| Views | State physical y | P1 / moments physical y | Ordinary Z |
|---|---|---|---|
| 17 pre-Rc, six outer and inner_flat | 0..4 | 0..3 | 0..3, Q1/V1 only 0..2 |
| inner_core, inner_first, lower | 0..2 | 0..1 | 0..3, Q1/V1 only 0..2 |

Core rows are lifted from the actual saved rho rows by D_y=rho D_rho.
Utheta1=sqrt(2R)F1 and M1=R mean1 use the corresponding weighted product
rules. All five cumulative first derivatives use their saved signed densities.
The original core endpoint is rho=4; the first chart is s in [0,sc/4]; lower
is s in [sc/4,sc/2]. The high inner-flat view is restricted to
[3sc/16,sc/4]. No local high-row formula is extended over the missing interval.

**Task 13 remains IN PROGRESS.** There is now a working regional dispatcher,
but no accepted uniform global y4/pressure-y3 interval contract. Next:

1. Construct genuine inner third/fourth radial rows outside the flat cutoff
   collar from the same Volterra problem. Retain C3 state and C2 Q/V axial
   source rows and prove their derivative source wiring.
2. Extend actual pressure and all five cumulative moments through physical
   y3 on those intervals; do not infer them from matching traces.
3. Complete and check the lower-collar high own-field interval contract,
   keeping the original nonzero incoming mass and signed pressure history.
4. Add regular axis/Cartesian recovery separately, including parity and
   component limits, then resolve automatic physical-radius chart selection
   with original-scale numerical or interval services.
5. Compose the original similarity/physical map into actual velocity and
   pressure recovery. Original numerical evaluation, common N and certified
   bounds remain separate from the source API.

Global n1 numerical fields, complete physical velocity export, compatible
preheat/exact heat exterior, admissible stress/flat remainder, genuine n>=2
recursion, smooth sum, oscillatory cancellation and full corrected NS remain
open. The long-term objective stays **ACTIVE / INCOMPLETE**.

## Focused evidence

```powershell
python scripts/run_original_n1_source_dispatcher.py --check --write-report
```

The single focused run exited 0 in **11.437 s**. It addressed all **6318**
admitted scalar derivative slots at three native parameters (**18954**
source queries), checking exact saved source/coordinate/radius wiring.
**18** invalid domain/derivative/type requests were rejected without graph
mutation. **Five** independent symbolic operator identities checked the core
physical-radial lift and weighted products. These are adapter checks, not new
NS residual, smoothness, numerical enclosure or blow-up evidence.

Reports:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_source_dispatcher.json`
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_source_dispatcher_check.json`

## Stream-disconnection recovery

A disconnected reply does not determine the saved work state. Before retrying,
use this checkpoint and source archives. The old exec session 84468 is gone;
its in-memory owners are unavailable. The current entry point requires neither
that session nor historical replay. The prior recovery record found this
thread's server-closed stream/remote-compaction failures; the server-side
cause of request c985edb3-6c8a-43d6-9580-ded1325339ae remains unestablished.

Resume only the five unfinished items above. Do not recreate the completed
seven-interface milestone or relabel regional source queries as numerical
velocity evaluation or temporal recursion.
