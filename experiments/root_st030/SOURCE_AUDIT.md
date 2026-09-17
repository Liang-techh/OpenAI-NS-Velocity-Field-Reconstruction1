# Primary-source audit, 17 September 2026

Source statements, autonomous implementation choices, and numerical observations
are kept distinct. None of the sources below supplies an independently verified
full candidate satisfying this project's exact fixed-force/support/core contract.

## Official OpenAI
https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
https://openai.com/index/navier-stokes-solution/

The theorem and introductory construction, oscillation, heat/exterior and
localization pages were inspected including PDF images. Their background-plus-
oscillatory-plus-mean-correction construction is not the same as a low-degree
axisymmetric field. Their general smooth forcing is not this project's fixed
solenoidal two-parameter curl family. Potential-level localization motivates
ST031 complete curls; its specific annular modes are autonomous, not extracted
paper coefficients. No arbitrary residual-defined force is used for acceptance.

## Kokuno corrected reconstruction
Pinned source: KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7, navier-stokes/RESEARCH_STATE.md.
Corrected edition: https://doi.org/10.5281/zenodo.22678406

The state explicitly leaves imported profile/stress existence, all-stage
realization, final force for the same assembled witness and independent endpoint
validation unfinished. Component indices/check counts are not an accepted NS
candidate. GitHub index/state were read; bulk TeX/ZIP download failed DNS.
No external replay program was run in this round.

## Duraiswami / swirl-collapse, submitted15September2026
https://arxiv.org/abs/2609.17642
https://arxiv.org/html/2609.17642v1
https://arxiv.org/src/2609.17642v1/anc/README_anc.txt
https://gitlab.umiacs.umd.edu/ramanid/swirl-collapse

Pressure gauge, Newton/continuation, axis Cauchy versus unstable Dirichlet and
asymmetric matching inform the optimizer-versus-ansatz experiment design.
The reported0.2% profile moment match is not a fullCartesian NS max/L2 gate.
Paper and ancillary index were readable. Raw Python requests failed web content-
type handling; direct downloads failed DNS, and no GitLab connector was found.
Thus no core_solver/axis_march or other external numerical routine was executed.
ST032 reuses this repository's existing asymmetric basis, not that paper's solver.

## Pressure robustness
https://arxiv.org/abs/1906.03009 — Linke, Merdon, Neilan.
https://arxiv.org/abs/2401.04456 — Di Pietro, Droniou, Qian.

Motivation: distinguish pressure gradients and solenoidal velocity dynamics.
Our compact-pressure fixed-column QR variable projection is an independent
finite-dimensional least-squares construction, not their FEM discretization or
a continuum Helmholtz projector. Exact derivative formulas were tested directly.

## Compact representations / Euler comparison
https://arxiv.org/abs/1810.08020 — Gavrilov.
https://arxiv.org/html/1810.08020v1
https://arxiv.org/abs/cs/0502092 — Deriaz/Perrier.

Compact divergence-free representations are meaningful; compact support alone
is not an impossibility argument. A compact steady Euler field is not this
viscous prescribed-force/core-scaling problem. No Euler benchmark is substituted.

## Residual-adaptive sampling
https://arxiv.org/abs/2207.10289 — Wu, Zhu, Tan, Kartha, Lu.
https://arxiv.org/abs/2210.00279 — Gao, Yan, Zhou.

Motivation: inspect residual peaks in a separate TRAINING pool, not just one
fixed collocation mean. ST033 is an independent polynomial-field sampler with
quartic peak penalty, not a copied PINN or an exact RAD/RAR implementation.
Earlier holdouts informed the new experiment direction; fresh ST033 holdout
seed9172912 was not used for training or coefficient selection.

## Parallel repository work, not merged
https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/pull/233
Exact reviewed head:4b77580b2bedbbc9b047471cb77b2229549192f2.

Source-native finite leading CORE plus complete-curl oscillatory probe, unchanged
core pressure and f=0 raw diagnostic. Its24point/3time rawRMS comparison improves
about1.74% but remains O(1)-O(10). No global outer/heat matching, remeasured mean
correction on that core, or final compatible force is attached. That evidence
cannot validate an ST030-series field or the complete original-domain target.
