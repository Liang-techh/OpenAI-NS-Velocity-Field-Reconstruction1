# F33: actual trial-amplitude axial correction map

The Section 7.31 two linear moment equations now have a callable same-family trial map for a in [.9,1.2], using the actual F32 incoming moments. `experiments/root_st073/lei_ren_part1_paper_shared_outer_pulse_map.py/.json` preserves the unique smooth functions c1(a,Z), c2(a,Z) implicitly and encloses their finite factors and fixed-a axial derivatives. This is not the energy-selected ap or a completed O.4 velocity field.

## Stable inverse and exact source normalization

The two rows differ by microscopic mu. Direct subtraction loses this difference at ordinary precision. For x=d-s, d=3 or 1, the new matrix uses row1 = integral exp(-x/2+mu*x) beta(s) ds and the exact divided difference D=(row2-row1)/mu. Its quotient (exp(mu*x)-1)/mu is the positive integral from 0 to x of exp(mu*v), enclosed between x and x exp(mu*x). The divided-difference determinant lies in [-0.271848,-0.269975], so the inverse has a finite nonzero gap.

Incoming functions are m1=Mz/(Rp Utheta(Rp)) and m2=Mtheta_z/(sqrt(2) Rp^1.5 Utheta(Rp)^2). Both use the actual F32 cumulative primitives and physical Pstar factors. After common pulse rescaling, their C1 magnitude is bounded by exp(-1000). Symmetric value/derivative boxes enclose these true sign-indefinite functions; the positive cap is a magnitude bound, not a positivity claim or a replacement by zero.

## Common pulse saddle and retained tails

The full weighted pulse integrals share k0=1/(2mu), u0=(4mu)^(1/3), L=u0^-2 and w=u0/sqrt(6). With v=1+xw, the exact phase difference is phi=x^2(2v+1)/(6v^2), with phi''=v^-4>0. Both rows use the same formal logarithmic prefactor and finite row factor exp(i(2+u)). Directed closed cells enclose the centered integrals; convex tangents bound the Gaussian tails without a spurious large L multiplier. Earlier pulse support has a separately proved, nonzero exp(-1000) magnitude cap.

Coefficients are represented as c_j=exp(common_logpref-log(mu))*C_j. The exponential is not materialized or rounded to zero. The finite C_j inverse solves row1*C=-mu*R1 and D*C=-(R2-R1), which are exactly equivalent to both original equations. The resulting signs c1<0<c2 are separated for the recorded trial range.

The independently enclosed fixed pulse-energy constant is Kpulse in [0.2443923749493,0.2457111997123], inside the paper's (0.24,0.246). This is distinct from the preheat-pressure bound Kp=17. `shared_outer_pulse_map_check.py/.json` independently verifies twelve physical-row, normalization, divided-difference, saddle and full future-energy identities plus five interval/source examples.

## Remaining dependency

Section 7 determines waiting tau, then angular (d1,d2), then axial (ap,c1,c2). The energy equation includes the entire future corrected swirl, including angular bumps and the exact heat tail. The target at Rv is half of the positive swirl-energy tail beyond Rv, not zero; Mztheta vanishes at infinity. The new map leaves this term explicit and ap unselected. Next build the later angular candidate/heat data, solve the angular pressure/moment corrections, then evaluate the actual energy root and install the coupled pulse field. Whole outer cone, exact heat matching, global stress lift and temporal recursion remain incomplete.
