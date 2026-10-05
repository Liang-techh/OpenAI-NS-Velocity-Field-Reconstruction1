# Original main and exit: complete similarity source companion — 2026-10-04

The original main xi in [.02,10] and exit xi in [10,11], with Z in [-1,1],
now have all five raw cumulative moments through mixed4, actual similarity
velocity through mixed4, the same absolute pressure through mixed4, and full
meridional stress through mixed3. Ordinary derivatives use y=log(R/Rp),
d_y=mu*d_xi. The angular normalization is Mtheta/(sqrt(2)*R^(3/2)*Utheta).

## Defining sources and numerical bounds

The new uncapped_selected_pulse_source callable uses the exact incoming row
factors exp(-13lambda_i/mu-logG), the original divided-difference inverse,
nu_j=mu*Kj*exp(2(logG-log(mu))), and the positive quadratic branch. Both scalar
and Taylor modes check the positive branch; Taylor mode uses the exact
implicit recurrence through every supplied axial order. logG is the common
pulse-integral prefactor, not logE. Original full gp, beta, Gram, incoming and
corrected-future integral definitions and their current source paths remain.

Production coefficients remain directed C5 enclosures of these functions.
No point parameters, root, finite quadrature row, pressure box, or positive
factor cap is selected as an exact production value. Exact production point
evaluation remains false. This increment adds a source companion, not a
replacement physical field.

The full source retains local gp forcing, separate original incoming Mz and
Mtheta_z histories, signed angular memory, partial main squared-energy
integral, selected terminal energy loss, and canonical absolute pressure
memory. Forward linear integrals keep their positive omitted tail; terminal
end histories are not substituted over the main pulse. The energy uses the
same selected identity to restore the backward representation while keeping
its local gp integral. Raw Mp and the analytic datum remain separate.

## Stress, derivatives and interfaces

All original transport and shear terms remain. The stress split has local
terms, distinct incoming histories, the pressure memory, and the selected
energy loss; the nonzero axial input and full Uz_y are retained. Exact
ordinary rows 0 through 3 of the split equal the original full stress formula
for arbitrary axial functions.

Raw incoming Mz and Mtheta_z, angular memory, absolute-pressure memory and
selected-loss radial rates cancel exactly. The mixed rows therefore retain
their histories without giant exponential materialization. The same function
covers xi=10. The current original functional xi=11/gap certificate supplies
the selected linear and energy joins, Z derivatives through 5, common
pressure data and source ODE compatibility.

37 source identities and 27 actual source AST bindings are installed. The
current admitted gap source receipt and C5 source receipts are consumed
without rerunning unchanged upstream or downstream chains.

## Focused verification

Focused checker PASS: 429 current hashes, 37 source identities, 27 source AST bindings and 2025 finite signed rows (stress 750, velocity 300, five raw moments 825, absolute pressure 150). Three independent original-integral/full-stress fixtures give 169 checks, including two nonzero-shear cases. Normalized tolerance 1e-50, maximum positive enclosure miss 9.83375766025783602734756116449e-92. C5 positive-root implicit derivatives and invalid-branch rejection pass; this is formula/unit evidence, not corrected NS accuracy.

The independent moderate-parameter fixture uses the original startup/exit gp
integrals, both complete beta weights and the original full-stress evaluator.
It checks the uncapped positive root, both terminal linear equations, the
quadratic equation, C5 implicit derivatives and rejection of an invalid
positive branch. Its forward energy is independent of the companion's
backward energy. Physical comparisons are normalized before testing, so tiny
amplitudes cannot make a zero result pass.

Unchanged independent quadrature is saved with the generator hash, precision
and mpmath version and reused. This fixture cache is not production data.

Run the bounded stage:

    python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulsemainexitsimilarity

## Next implementation

1. Lift the multiple velocity/history sectors to the completed physical
   tensor, diagonal and divergence, and all three remainder components.
   Include every local/incoming cross product in radial convection.
2. Close xi=10 and xi=11/gap physical joins through stress3,
   pressure4, velocity4 and remainder2 in the unchanged nu/time/radius units.
3. Derive continuous main/exit cone bounds with
   sigma=2Uz_y/((2+2mu)Utheta) and
   kappa-2=2mu+(2+2mu)sigma^2. Preserve gp/moment correlations; zero-shear
   or tiny terminal estimates do not apply.
4. Complete entrance/upstream feedback and global tensor admissibility,
   temporal flatness, volume norms and required-domain energy.
5. Implement actual n-dependent recursion, oscillatory correction and the
   final corrected Cartesian residual and measured vortex dynamics.

Physical main/exit completion, main/exit cone, global admissibility, global
flatness, energy, recursion and corrected NS validation remain false.
