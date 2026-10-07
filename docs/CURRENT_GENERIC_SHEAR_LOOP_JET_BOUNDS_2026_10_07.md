# Actual generic-loop inverse and phase-held primitive derivative bounds

**Successor:** signed source derivative expressions are now exposed in [CURRENT_GENERIC_SHEAR_SIGNED_JETS_2026_10_07.md](CURRENT_GENERIC_SHEAR_SIGNED_JETS_2026_10_07.md). Actual point-loop/higher-order rows, changed transport/new repair/N and true recursion remain open.

Checked implementation: [6acd9f09](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/6acd9f091b800f149ef8bd88094242c825211b94). The actual original17-chart source cover now supplies uniform logarithmic bounds for the Section11 inverse phase and original loop primitives through y,Z,yy,yZ. The adapter also supplies phi,phi_phi,y_phi,Z_phi bounds needed to insert one common fast phase. These are bounds of the loop defined by the original analytic source expressions; source caps never define point field values. Signed point-loop evaluation, higher mixed rows, changed histories/new repair and common finite N remain open. The long-term goal remains active.

## Implemented API and units

`CurrentLoopJetBounds().chart(chart)` consumes the checked whole actual generic input/scales receipt, the original positive denominator theorem, full ordinary source quotient norms and the original checked sigma derivative caps. `run()` assembles17 charts. Producer/checker use `lei_ren_part1_paper_compliant_current_generic_shear_loop_jet_bounds`.

Available fixed-phase orders are (0,0),(1,0),(0,1),(2,0),(1,1). The six source y2/Z1 rows do not by themselves justify y2Z for the nonlinear loop, Z2 or full mixed4. A is dimensionless; returned B bounds are for the original B/Pstar, because E=Utheta/Pstar. Pstar is a fixed global constant under slow derivatives and must be restored for raw physical B.

The module forms new executable logarithmic majorants. It does not rebuild source graph ancestors, exponentiate actual amplitudes/eta/radius/width, or create an installed physical point field. Absolute caps may bound derivatives of the existing signed source expressions; they may not be substituted for signed a,b,p1,p2,E values in a later point solver.

## Branch and conditioning construction

Keep Delta=kappa-2 and eta separate. If Delta>=eta, q and every q derivative vanish; the exact loop is Phi=psi/(2pi), psi=2pi*phi and A=B=0, including derivatives. On active support,

    x=1-Delta/eta
    gamma=2eta-Delta>=eta
    q=sigma(x)*sqrt(gamma/(2a)), gamma<3.

The checked original flat cutoff supplies global sigma derivative caps32 and1792. q slow derivatives use sigma composition and the logarithmic square-root derivatives, with the same actual a_min and global eta. The cap q_star is inherited from the current source scale selection.

The inverse remains well conditioned in formal logs through the exact identity

    v/a=1+t0^2+2q^2>=1
    K=a/(2pi*v)=1/[2pi*(1+t0^2+2q^2)]
    lambda=Phi_psi=K*(1+t^2)
    1/lambda<=2pi*(1+t0^2+2q^2).

This works on both active and zero-q branches and avoids using an enormous independent kappa upper. For u=p2*q/d_star,h=sqrt(1+u^2),r=u/h, preserve

    1-|r|=1/[h*(h+|u|)]>=1/(2h^2)
    1-r^2=1/h^2
    D=1-2r*cos(psi)+r^2>=1/(4h^4).

No subtraction of two rounded unit values is used. Both positive lambda and Poisson D lower bounds are recorded as finite logs for every chart.

## Fixed-phase inverse and primitives

First bound t=t0+2q/h*w(psi,r), its slow derivatives and its psi derivative using the correlated Poisson denominator. Integral derivatives at fixed psi follow

    T1=integral_0^psi t, T2=integral_0^psi t^2
    Phi=K*(psi+T2)
    psi_i=-Phi_i/lambda
    psi_ij=-(Phi_ij+Phi_psi_i*psi_j+Phi_psi_j*psi_i+Phi_psi_psi*psi_i*psi_j)/lambda.

The primitive definitions remain

    A=a/2*(phi-psi/(2pi))
    M=-a*T1/(2pi)-b*phi
    B/Pstar=E*M/2.

The T1 chain includes the inverse angle:

    T1hat_i=T1_i+t*psi_i
    T1hat_ij=T1_ij+t_i*psi_j+t_j*psi_i+t_psi*psi_i*psi_j+t*psi_ij.

The implementation bounds these exact recurrences and supplies phase/mixed derivatives by implicit differentiation of lambda. All sums/products use directed logarithmic absolute majorants with exact-zero sentinels.

For the later common phase phi=N*log(R/r_minus), the total derivative is

    D_y f=f_y_at_phi+N*f_phi
    D_y^2 f=f_yy_at_phi+2N*f_y_phi+N^2*f_phi_phi
    D_Z D_y f=f_yZ_at_phi+N*f_Z_phi.

The velocity `modulate()` interface consumes slow A_y/B_y. Passing the total derivative there would double count the fast term. This adapter supplies bounds for those terms; it does not choose N or evaluate a changed velocity.

## Evidence

Passed: 7 independent exact second-order chain/normalization identities, 1054 finite derivative log caps across17 actual source covers, positive phase/Poisson conditioning on each, the invalid eta guard, and 96 comparisons against the independent scalar loop on modest active/cutoff-transition/flat artificial sources. Actual production scales were not exponentiated. Working/index audit matched 1003 hashes. Read-only mathematical reviewer: GPT-5.6 Luna / max.

## Next production tasks

- [x] **LEFT4c1-H/scales:** whole current H0-2 margin, complete original weak/strong branches, full11-unit O2 buffer and actual logarithmic scales.
- [x] **LEFT4c2-bounds:** actual17-chart inverse/A/B y,Z,yy,yZ and phase/first slow-phase log bounds, correlated positive phase/Poisson denominators and exact flat branch. This completes the bound subtask only.
- [ ] **LEFT4c2-signed-input-jets:** expose signed/factored ordinary jets of a,b,p1,p2,E using `SourceQuotient`, the raw numerator rows, formal R and actual positive E/C/a certificates. Preserve the existing four fixed log bases, full inertial pressure/energy and original history rows. Never evaluate a denominator box containing zero or use absolute norm caps as field values.
- [ ] **LEFT4c2-point-loop:** attach an actual-source analytic/factored loop and monotone inverse provider, using Delta and eta without rounded2+eta. Retain exact zero branch and one global slow-constant scale object. Use the new conditioning and derivative bounds for its domain contract.
- [ ] **LEFT4c2-common-phase:** rebase all original charts to phi=N*log(R/r_minus), including the formal microscopic left offset. Restore raw Pstar units only at the physical consumer and keep slow/total derivative interfaces separate.
- [ ] **LEFT4c2-required-orders:** produce sufficient original signed source rows and higher implicit recurrences for physical velocity mixed4 and full signed stress mixed3. Current y,Z,yy,yZ admission must not be relabeled full mixed4.
- [ ] **LEFT4c3-own-transport:** express/integrate all five changed moment defects over each formal source interval. Keep unchanged analytic P0, incoming original moments and quiet-gap memory. Derive radial and absolute pressure recovery from those changed histories, including all pressure/energy/meridional cross terms.
- [ ] **LEFT4d-new-repair:** solve the new functional five-bump problem on(Rc,2Rc) with the actual Ac/mu/P0 and newly transported defects; certify inverse/uniqueness and terminal identities as functions of Z.
- [ ] **LEFT4c1-outer:** attach unchanged original admissibility from2Rc through Rb, all later interfaces and exact-zero heat exterior. Modification-box input is already admitted; post-repair outer admission remains separate.
- [ ] **LEFT4e-new-N/cone:** assemble all actual loop/transport/repair constants before selecting one new finite whole-source N and proving the modified stress cone. Historical scoped N is not the new common threshold.
- [ ] **CONT / ENERGY:** changed source interface identities, high derivatives and whole physical energy/tail.
- [ ] **REC:** genuine n=1 and n>=2 recovery equations, independent moment repairs, finite-order remainder and smooth summation.
- [ ] **WAVE / PHYS:** two pulse families, mean corrections/averaged quadratic stress cancellation, corrected uvw/p/f, independent Cartesian residual and measured recursion/contraction/elongation/material winding.

True new gate: `current_actual_phase_held_inverse_and_loop_y2_yZ_log_bounds_certified`, with `phase_held_loop_primitive_derivative_bounds_certified=True` scoped to the explicitly listed orders. Signed installed point-loop, higher derivative admission, changed moments/new repair/common N, global cone, actual coefficient recursion and full corrected NS remain false. Original cone inventory remains15 strict nonzero regions/17 open plus exact zero exterior.
