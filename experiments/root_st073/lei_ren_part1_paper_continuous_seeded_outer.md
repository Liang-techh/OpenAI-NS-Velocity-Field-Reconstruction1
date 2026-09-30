# Installed continuous seeded axial exterior

ContinuousSeededAxialProfile replaces the pulse and translated end-bump axial values and cumulative mean together. It retains the original schedule, angular correction, heat tail and pressure objects. ContinuousSeededOuterField uses the existing joined physical velocity wrapper and recovers Ur from the same seeded mean and its Z derivative.

The complete seeded incoming record is retained; pre-Rp transport is inherited and its mass offset is applied once by the existing joined adapter. After Rp, normalized seeded base rows enter the continuous component directly. Before Rv, M/R is E exp(-lambda end) times the normalized cumulative row. After Rv, it is Ev exp(-end) times the terminal row, so cumulative mass remains constant even through later heat stages. A terminal row is evaluated, not reset to zero.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_seeded_outer.py`.

The receipt samples Rp, startup, the plateau pulse, cutoff, first end bump, Rv and after-Rv at Z=.3. It separately records object identity, post-Rv mass conservation and the recovered radial velocity/divergence diagnostic. Prepared .3 seeded/continuous receipts are the actual prior shared candidate; neighboring Z inputs are solved on demand from the retained actual inner terminal.

This is a provisional physical installation. The shared pulse is now used in values and primitive definitions, but numerical primitives retain bounded omitted pieces and uncertified quadrature. Inherited incoming quadrature, unknown inner entries and float-Z caching/derivatives remain explicit. A small arithmetic terminal mean is still nonzero and cannot certify global finite energy. Outer jets, full five moments, uniform stress/remainder and recursive/oscillatory closure remain open.

## Installed radial diagnostic

The actual adapter recovered Ur from its own continuous seeded mean. Independent radial differences at the pulse plateau give mapped-divergence relative cancellation 6.24e-24 and 3.90e-25 as radial step decreases from 1e-5 to 5e-6. This is a local diagnostic with the existing float-backed axial stencil, not a uniform divergence bound. Post-Rv mass propagation replays to 2.10e-443. All four schedule/angular/tail/pressure object identities are retained.

The startup primitive uses integration by parts; truncating the subtracted sigma integral produces a NEGATIVE omitted correction. The component and physical mean receipt now propagate an absolute bound and correction sign, with no positive-bound label for that branch. Startup .015 at actual candidate mu has been replayed through the component. Quadrature error and incoming uncertainty remain separate and unenclosed.

Installed startup point xi=.015 also has nonzero Uz and propagates omission sign -1. Seven named profile samples are now in the receipt.
