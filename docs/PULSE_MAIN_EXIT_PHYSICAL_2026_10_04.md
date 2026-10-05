# Original main/exit completed physical tensor — 2026-10-04

The original main xi in [.02,10] and exit xi in [10,11], Z in [-1,1],
now have completed physical stress through mixed order 3, diagonal and
divergence through mixed order 2, and all three remainder components through
mixed order 2. The source is the accepted main/exit five-moment, velocity
and absolute-pressure companion. This increment does not certify a main/exit
cone, a globally flat remainder, global energy, recursion or corrected NS
accuracy.

## Physical operators and nonlinear histories

The accepted physical lift is replayed with its stress, diagonal, divergence
and Cartesian tensor loops unchanged. Its remainder inputs contain nine
radial sectors: local/incoming time terms, local/incoming radial viscosity,
local nonlinear transport, terms with one incoming history, the incoming
history square, and local/incoming axial viscosity. Theta and axial errors
retain their full axial-viscosity terms.

Ur*Ur_r retains local-local, local-incoming, incoming-local and
incoming-incoming products. Uz*Ur_z retains both local and incoming radial
histories. Exported ordinary source rows already differentiate the source
exponentials. Frozen exact source logs are restored as factors at the
basepoint and are not differentiated a second time.

The tensor completes Ttheta_theta=r*partial_z(Tz); radial tensor divergence
and velocity divergence remain exactly zero. The Cartesian tensor,
divergence and remainder use the same fixed positive viscosity, original
implicit lambda, physical time, radius and pressure units as the admitted
gap/end companions. Local main/exit D is 1; incoming histories and terminal
energy loss carry their separate exact logs.

## Source-functional interfaces

The same original main/exit source function covers xi=10. At xi=11 only the
local gp input and its required jets vanish. Its forward cumulative
convolutions generally remain nonzero. For both linear moment histories,

    local_i(11,Z)+exp(-11*lambda_i/mu)*incoming_i(Z)
        = D_i(2)*selected_gap_M_i(Z).

The actual row/velocity algorithms are replayed using the current selected
two-row inverse certificate, the defining forward integral and original
flat-endpoint receipt. All required mixed velocity rows agree; the entire
stress ordinary rows through 3 and energy/loss/absolute-pressure rows
through 4 agree as axial source functions. The full radial history square
is retained. Nonzero Er and Etheta are not reset at the interface; Ez
vanishes there because the local axial input is flat.

Actual main/gap logR, logB, logH, pressure-memory Q and selected moment
factors agree. The native pressure calls share inlet, swirl, canonical datum
and time 11/mu. The physical calls pass identical log_tau, theta and
viscosity arguments. Production pulse_exit(11), pulse_gap(11) and
pulse_gap_end(-2/mu) all give logRp+11/mu. Chart-local coordinate labels are
not treated as physical radii.

## Focused verification

Focused checker PASS: 437 current hashes, 189 source/operator/interface identities, 24 AST bindings, 1346 finite signed physical rows and 94 structural zero rows. Consume 28 admitted full physical operator identities. The independent nonzero-Uz Cartesian fixture at nu=.01/.7 gives 114 checks through stress3 and diagonal/divergence/error2, with all three errors nonzero. Tolerance 1.0e-55, maximum positive enclosure miss 2.19168267783394594220548535231e-89. This checks formulas/units, not corrected NS accuracy.

The current source proof has 189 identities and 24 AST bindings; the
accepted 28 full meridional physical operator identities are consumed
without rerunning unchanged upstream or tail stages.

The independent Cartesian oracle uses the original gp startup integrals
and exact plateau primitives, nonzero local axial input, nonzero incoming
radial histories, and moderate smooth axial fixtures. It directly computes
time derivative, convection, pressure gradient, the full vector Laplacian
and completed tensor divergence from the original implicit lambda relation.
Both viscosity values .01 and .7 exercise stress mixed3 and diagonal,
divergence and three-component error mixed2. All three fixture errors must
be nonzero. Fixture parameters and integral values are not production
point selections or a cone certificate.

Use only the bounded increment:

    python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulsemainexitphysical

## Next work

1. Derive original full-shear main/exit stress-cone expressions with
   sigma=2Uz_y/((2+2mu)Utheta) and
   kappa-2=2mu+(2+2mu)*sigma^2.
2. Bound the complete main/exit domains continuously, retaining gp,
   gp_xi, amplitude and moment correlations and all signed stress terms.
3. Compose the original xi11 physical join with the accepted gap/end/tail.
4. Restore entrance and upstream finite-width feedback, completed global
   tensor admissibility, global temporal flatness, volume norms and
   required-domain kinetic energy.
5. Implement actual n-dependent coefficient recovery, independent repairs,
   oscillatory stress correction and final corrected Cartesian residual
   and measured vortex/particle dynamics.

Exact production point parameter selection remains false. The present
artifact is a source-bound regional physical decomposition and a generic
formula/unit test, not a final globally evaluable corrected velocity field.
