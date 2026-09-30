# Heat functional with separate tiny deficit

The new provider evaluates H(xi) = E[(1+xi V)^(-h)] with V distributed as Gamma(1+h,1), and first/second derivatives from the same integral. It retains deficit=1-H, logH, and the correction to H-prime at zero separately; rounded full H is not a substitute for these atoms. Positive quadrature is normalized by h before integration, avoiding absolute tolerance loss of a whole tiny integral.

For xi<=1e-8, one quadratic Taylor polynomial supplies all three jets. With c_m=(h)_m(1+h)_m, analytic truncation bounds are c3*xi^3/6, c3*xi^2/2 and c3*xi. These cover analytic truncation only, not finite precision arithmetic. The branch transition is an approximation within these bounds, not an exact globally smooth heat solution.

At h=5e-201, finite-difference functional checks at xi=.2 and 1 have maximum relative difference 4.28e-19. At xi=exp(-1e6), nonzero deficit, logH and derivative correction remain represented. Quadrature uses 80 working digits; working precision is not a certification of accuracy.

This provider is NOT installed in the shared schedule, angular moments, pressure or live correction targets yet. Next install it using logH and separate deficits (never subtract H-1 after rounding), propagate all target/jet changes coherently, and add heat cumulative/tail moments. Complete heat matching, finite energy and recursion remain uncertified.
