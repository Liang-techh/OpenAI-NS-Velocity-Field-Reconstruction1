# F26: completed fresh finite core, local exit chain, and global frozen profile

Date: 2026-10-01. This checkpoint supersedes the pending computation and frozen-profile gates in F25. The reconstruction objective remains incomplete.

## Completed results

The fresh inlet-tolerance core has completed radial order 144/144, with 148 initial axial jets and retained axial depth three. The final resumed batch performed 125 updates in 1994.579 seconds (this is the last batch's compute time, not total project time). Its terminal state SHA256 is:

    16477c4ce787db62ce2f28bd345ce011c77da9500a2e6ce4c64370233edbdb9b

Seed acceptance passed after termination. The finite-plus-same-source-tail exit loader produced 20 mixed C3 enclosure rows. The new-parameter comparison produced 32 phase cells; the actual bridge produced 48 phase cells; continuation reached physical R=100; the controlled switch used 32 cells and exact constant-power moments to reach R=110. The local-axis relaxed (3.23) direction condition at R=110 passed. All five downstream receipts share the same parameter family, and 57 phase/derivative checks passed.

These finite results apply to Z in [0.49,0.51]. They use the new j approximately 1.7478350439606e-22, not the old-j branch. The shared width and exit epsilon are the identical symbolic scalar h_b=epsilon_b=cstar*K^-100. The physical K and fixed K1 bounds needed to admit this scalar numerically remain unresolved.

Separately, analytic companions now establish every frozen-profile condition in (9.14) on the full real axial domain Z in [-1,1], over physical Ra=4/Lambda through R=110:

| Gate | Result |
| --- | --- |
| Df>0 on Ra..110 | Passed |
| Df>=4 on 100..110 | Lower bound approximately 150 |
| Hf=(Df^2+Ef^2)/Df>=2+4gamma=2.04 on Ra..110 | Passed by two analytic regions |
| high-chi region chi>=0.99 | Hf>=Df>=2.8913664589602 |
| low-chi region chi<=0.99 | Positive pressure-dominated bound far above 2.04 |

Rounded values are summaries. Exact interval tuples in the JSON receipts are authoritative.

## Angular positivity proof

The same-source analytic fixed point and uniform core bounds from F25 supply the angular equation

    2r Phi_rr + (4-epsilon*W*r/L) Phi_r = forcing.

Square completion bounds the mixed chi and sqrt(chi) terms and gives forcing<=-epsilon*margin, with margin>0. Explicit guards establish a positive ODE coefficient below five and an initial slope strictly above the barrier. Thus -Phi_r>=epsilon*margin/5. At scaled r=4 this yields Qentry>=0.31845729166667, where Qf=L*Df/R.

With theta=Ra/R, the frozen identity is

    Qf=theta^2*Qentry+(S0/2)*(1-theta^2)
       +[A(mu_z)/Ra]*theta*(1-theta).

Here mu_z=Mz(Ra)-Ra*v. The moment-error ledger bounds the pointwise normalized quantity |A(mu_z)/Ra|, not the unnormalized moment. Using theta*(1-theta)<=0.5*(1-theta^2) gives a convex lower bound between positive Qentry and Qreference. Since 0<L<=1, Df>=R*Qf. This proves positivity and the terminal lower bound.

The angular receipt alone deliberately leaves the Hf test false; the companion H receipt completes it.

## Full Hf proof

For chi>=0.99, the leading-model polynomial for 3B+2qB' is negative at the entry, including the analytic correction. Consequently Dentry>=3. The signed logarithmic-gradient decomposition gives a positive frozen source coefficient greater than 1.9. For x=R/Ra>=1,

    Df>=3/x+1.9*(x-1/x)>=2*sqrt(1.9*1.1)>2.04.

For chi<=0.99, the axial polynomial localizes Z to [-j/3,-j/5]. The full core/frozen velocity and its axial derivative obey |V|<j and |V_Z|<5. Their cumulative integral averages obey the same sup bounds. The energy-average and pressure-increment ledgers retain all physical amplitude terms and their axial derivatives.

All fourteen preheat pressure atoms have nonnegative contributions to the pressure operator in this region. The six beta-two atoms give a lower bound 0.25*j*Pstar^2. Subtracting the velocity, cumulative energy, and pressure-increment bounds leaves N>=0.125*j*Pstar^2>0. The implicit physical amplitude is strictly positive and bounded above by epsilon^2*exp(-1000); it is not replaced by zero. Thus |Ef|>=N_lower*sqrt(2*epsilon)/f_upper and Hf=Df+Ef^2/Df>=2|Ef|>2.04.

Read-only mathematical review checked the two-region argument and requested explicit normalized moment units and cumulative average bounds. Both are now recorded in the producers and receipts. These global analytic proofs do not consume a local finite coefficient tensor.

## Evidence and identities

Under experiments/root_st073/, prefix lei_ren_part1_paper_:

- shared_interval_core_Z049_Z051.json and paired state: terminal 144-order finite core.
- shared_core_seed_check.json: fresh terminal seed acceptance.
- shared_tolerance_core_exit.json: finite-plus-tail exit enclosures.
- shared_tolerance_comparison.json, shared_tolerance_exit_bridge.json, shared_tolerance_exit_continuation.json, shared_tolerance_exit_switch.json, shared_tolerance_R110_cone.json: completed local exit chain.
- shared_tolerance_phase_check.json: 57 checks and shared parameter identity.
- shared_frozen_angular_bounds.py/.json: whole-axis Df positivity and terminal Df gate.
- shared_frozen_H_bounds.py/.json: complete frozen (9.14) input test.

Analytic core family SHA256:

    a4056322a3946c0ec0d43cc14d6f82d9a022410a5c077bab3577e87f8e9d7bdf

Local exit parameter family SHA256:

    933017f45cc5cb9d5b090a5de8afa8985e12b504f053fc8349242c8521ac17f7

The historical core-generation receipt still reports finite_plus_tail_exit_validation_completed=false because that producer only generates the core. The companion exit receipt records the now completed exit validation. Similarly, historical F25 parameter metadata is not rewritten to claim full admission.

## Ordered remaining work

- [x] Complete fresh 144-order recurrence and terminal seed acceptance.
- [x] Execute the fresh local comparison, actual bridge, continuation, switch and R110 direction chain.
- [x] Prove all full-axis frozen-profile inputs (9.14).
- [ ] Construct a fixed Section 9 K1 ledger covering cutoff interpolation, products, reciprocals, cumulative moments, D/E/H perturbations, and actual bridge direction estimates. K1 must be independent of j, profiles, K, Pstar, Lambda and Cstar; do not reuse unrelated constants named K1.
- [ ] Bound every physical C3 term defining K on the required core/frozen coordinates. Normalized derivative bounds and compact finiteness alone do not constitute a numerical K certificate.
- [ ] Establish radius compatibility (9.17), including the long-reshape length 400A. Check whether the current minimum Cstar must increase; do not assume the amplitude guard alone supplies radius admission.
- [ ] Transfer frozen margins to the actual connecting profile across the full domain, including bridge and switch. The local R110 result does not certify the whole transition cone.
- [ ] Build the long reshape and solve the five functional terminal moment defects, with controlled axial derivatives.
- [ ] Attach the exact heat exterior and restore its corrected pressure datum. Preheat pressure alone does not complete this matching.
- [ ] Construct the admissible stress lift and bound the flat remainder for the assembled background.
- [ ] Implement actual temporal n-dependent scale recursion and the two-family oscillatory corrections with averaged quadratic stress cancellation.
- [ ] Assemble u(x,y,z,t), v(x,y,z,t), w(x,y,z,t) for the completed field and evaluate the independent full Cartesian momentum residual and blow-up diagnostics.

Radial order 144 is spatial series generation. True temporal scale recursion, whole-axis finite evaluation, global actual stress admission, the full heat-matched background, and the final forced NS reconstruction remain incomplete.
