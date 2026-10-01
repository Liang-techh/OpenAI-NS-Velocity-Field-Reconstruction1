# F31: selected-family reference and initial outer field

The actual F30 implicit terminal five-moment identities now feed a new same-family O.1 reference continuation and O.2 outer candidate: slope transition, axial turnoff, and the eleven-unit zero-axial buffer through Rd. This construction uses the selected F27 radius, original Md40 schedule, actual delta and analytic preheat pressure P0. It does not reuse the old hardcoded A/logC physical field.

## Local units avoid materializing the radius

`experiments/root_st073/lei_ren_part1_paper_shared_outer_initial.py` exposes `SharedOuterInitial.reference(Z,offset)`, `.slope(Z,y)` and `.axial(Z,phase=...)` or `.axial(Z,buffer_offset=...)`. The coordinate y is the exact relative logarithm log(R/Rref); Rref = 110*(Cstar*Pstar)^10 stays formal. The code is callable for real axial boxes inside [-1,1] and retains first axial derivatives. It reports utheta/Pstar, Uz, P/Pstar^2 and Ur/sqrt(R/2); a full physical-coordinate evaluator and Ur_Z are still absent.

Let m=Mz/R, h=Mtheta/(sqrt(2)*R^1.5*Pstar), k=Mtheta_z/(sqrt(2)*R^1.5*Pstar), e=Mztheta/(R*Pstar^2), p=Mp/Pstar^2, and u=Utheta/Pstar. The original physical primitives imply the five local ODEs

    m_y = Uz - m
    h_y = u - 3h/2
    k_y = u Uz - 3k/2
    e_y = Uz^2/Pstar^2 - u^2/2 - e
    p_y = u^2/2.

The exact functional F30 closure supplies the reference moments after the bump support. These identities are used only outside the completed repair, not as a reset inside the repaired patch. O.2 then transports every inherited moment and computes P=P0+Mp, radial velocity and Q/N from those same primitives. The axial and mixed moments remain accumulated after Uz becomes zero.

## Stable turnoff integrals

The turnoff uses the original B=1-sigma(log(y)/Md), with y from 1 to exp(Md). Its normalized accumulated moments use integrals of exp(s-y) B(s) and exp(s-y) B(s)^2. Substitution r=y-s gives a decaying endpoint kernel. Each closed cell integrates its exponential weight exactly and bounds B by monotone endpoints. Beyond the chosen 800-unit kernel window, a positive enclosing tail is retained. The symmetric flat-cutoff formula evaluates the small complement directly to avoid losing it in subtraction from one.

This prevents the exp(exp(Md)) accumulated mass from overwhelming a numerical grid and avoids subtracting absolute log radii of the far larger selected Rref. The physical pressure datum continues to include all fourteen original preheat atoms, including positive late tails; no pressure fitting is performed.

## Evidence and scope

The producer receipt contains ten callable examples through the Rd buffer. `lei_ren_part1_paper_shared_outer_initial_check.py/.json` independently verifies twenty identities: five physical-to-local normalization identities and five cumulative ODE identities for each of reference, slope and turnoff. Twelve interfaces across four axial positions preserve values and retained first derivative enclosures. A whole-axis C1 call and nonzero axial history after Uz=0 are also checked. Full smooth matching follows analytically from the original flat cutoff; high-order interface jets are not numerically evaluated here.

This is the initial outer candidate, not the complete corrected outer profile. O.3 slope-mu/power buffer, O.4 pulse and axial moment corrections, flatten/angular correction, steep restoration/waiting, selected-radius heat collar/exact heat coefficients, global cone/admissible stress lift, flat remainder and temporal recursion remain unfinished. Finite example diagnostics do not certify the whole outer cone.
