# Current flatten -> power -> angular source and physical chain - 2026-10-05

The current reconstruction now has twenty-three accepted downstream profile source owners and twenty-three downstream physical source owners. This update adds `outer_power` and `outer_angular` after the checked current flatten. Earlier twenty-one-owner physical results are retained by their source hashes, registry and datum rather than regenerated.

This is regional derivative-source progress. It does not establish complete nonlinear point fields, global stress/cone/flatness/energy, or temporal scale recursion.

## Current source acquisition

`CurrentPowerAngularSourceAssembly` reuses one checked `CurrentPulseFlattenSourceAssembly`. `CurrentPowerAngularC4` inherits the original `power()` and `angular()` algorithms but replaces their legacy constructor and serialized input lookup with current defining objects.

The defining raw future is `pulse.fifth.fourth.energy.base`. This is the object used inside both the admitted C4 complete-future prefix and its C5 extension. It is a different instance from `pulse.high.base.future`; the previous F46 acquisition instruction pointing to that latter object is superseded. Independent instances are related by typed, pinned Md40/precision160 parameter constructors, defining-source assignments and common datum/source hashes, rather than object-identity assertions or endpoint equality alone.

At each real `Z` argument:

- Actual physical angular coefficients and Gamma energy jets come from `pulse.fifth.angular(Z)`.
- AST-bound `append_fifth` and `angular_fifth` preserve the first five C4 coefficients, including the original physical scale and Gamma prefix, while adding the fifth coefficient.
- The complete future-energy callable used by flatten retains the same admitted C4 prefix and C5 extension.
- Angular history, positive energy and absolute pressure come from `current_flatten.flatten(Z,100)`; pressure remains `P0+Mp` plus the original integrating-factor increment.
- Actual Rv five histories come from `terminal_histories(Z)`. Zero linear histories follow from empty native future supports. The nonzero angular/energy/pressure histories remain present; Rp reference data are not relabeled as Rv terminal histories.

The post-angular energy retains all terms: steep entry, steep power, steep exit, waiting, the infinite heat tail, `1/delta - 2*epsilon*W + epsilon^2*W_squared`, and the Gamma energy deficit with its complete source factor. Twelve new symbolic identities replay those actual current/C4/C5 production expressions. The existing power/angular functional theorem is consumed only after this current source bridge; interval overlap is not the join proof.

Original radius laws remain:

```text
flatten:       logR = logRp + 13/mu + t
outer_power:   logR = logRp + 13/mu + 100 + (Lrel-4)*phase
outer_angular: logR = logRp + 13/mu + 100 + Lrel + s
```

The source assignments prove flatten t=100 equals power phase=0 in radius, and power phase=1 equals angular s=-4. The common primitive ODEs, original flat supports and complete future decomposition identify the source functions and their mixed interface derivatives.

Source acceptance covers whole real Z in [-1,1], power phase in [0,1], angular s in [-4,0], both endpoints of each chart and four support-edge crossing intervals. Ten packets contain 600 velocity/pressure ordinary mixed logR/Z rows through total order four. A fresh unsaved Z exercises live coefficient/future acquisition. Run:

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentpowerangular
```

## Physical conversion

`CurrentPowerAngularPhysicalAssembly` dispatches both provider and evaluation through the same checked source assembly. It forwards the existing core, pre-pulse and native pulse graph through `source.before.dispatch`. Both new charts use the same current outer provider; their fixed unit comes from exactly the current flatten object used by that provider.

The unchanged mapper consumes the already differentiated mixed source grid and restores fixed Ev0 for all velocity components and Pstar squared for absolute pressure. It retains the original radius branches, moving cylindrical basis derivatives, Cartesian spatial derivatives through total order four and first physical-time derivatives at fixed position. Normalization is not differentiated twice.

Six new whole-Z physical packets cover power/angular inlet, full domain and exit. Each contains 35 Cartesian spatial multiindices for ux/uy/uz/pressure plus first fixed-position time derivatives: 216 source contributions per packet, 1,296 newly checked contributions including endpoint views. Adding the two whole-domain packets to the retained 4,506 regular contributions gives 4,938 regular contributions across twenty-three owners; the earlier 216 supplemental gap contributions remain separately retained.

An independent finite full-field fixture checks both original POST radius branches and 120 full ordinary mixed4 rows with nonzero synthetic velocity components, differentiating before fixed-unit normalization. The previous independent Cartesian/time and flatten-unit fixtures remain unchanged and are retained by hash. These fixtures do not admit Md40 production points.

Run:

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentpowerangularphysical
```

Default constructors consume their acceptance receipts before publishing the new certified flags. APIs are `CurrentPowerAngularSourceAssembly.evaluate(chart,Z,coordinate)` and `CurrentPowerAngularPhysicalAssembly.evaluate(chart,Z,coordinate,log_tau,theta)`; the physical adapter accepts the twenty-three admitted downstream charts only.

## Next work

First extend current ownership through steep entry/power/exit, waiting, heat collar and exact heat exterior using the same future/pressure defining objects. Admit each interface before extending the physical adapter, and retain completed packets rather than rebuilding them. Complete quantitative native pulse mixed4 interfaces and full common leading/remainder/point admission in their own scopes. Global tensor/cone, flat remainder and energy on the required physical domain remain separate gates.

True temporal recursion still requires distinct n=1 and n>=2 recovery equations, per-order moment repairs, finite-order remainder estimates and smooth summation. Coordinate scaling or adding these two spatial charts does not complete that recursion. Mean/oscillatory stress cancellation and independently measured corrected NS residual/dynamics follow afterward.
