# Complete nonlinear core map majorant

The new `lei_ren_part1_paper_nonlinear_map_majorant.py` implements the coupled map from paper(8.50)–(8.51), including all pressure and swirl terms. It uses the existing accepted complex-axis, linear inverse and pressure mass receipts.

## Map and coefficient estimates

The map is

    X = X0 + epsilon N(X),
    N(X) = (1/2)(R J2 Etheta(X), J1 Ez(X)).

The affine unit ball is centered at the paper's leading pair X0=(Phi0,Psi0), with the sum norm. The current separate norm upper bounds are Phi<=47313.40216951903 and Psi<=1.917518615499939e15 on this ball. These are conservative bounds relative to the accepted stored parameters.

The coefficient verification in(8.32) supplies explicit constants:

| Operation | Norm upper constant |
| --- | --- |
| Product | 256 |
| Radial average M | 1 |
| Radial primitive V or multiplication by scaled_R | 80 |
| J1 | 80 |
| J2 | 40 |
| Jnu(f scaled_R d_R g) | 20480 times norm(f) norm(g) |
| Jnu((d_Z f) g), including an extra scaled_R d_R g | (20480/h) times norm(f) norm(g) |

Fixed axial coefficients commute with radial inversion and are bounded by the existing Cauchy multiplier estimate. The angular non-derivative fixed terms combine exactly as -beta Phi. The term d H0/(H0^2+sigma^2) Psi Phi remains present, using the sharper denominator-factor identity

    H0/(H0^2+sigma^2)
      = (1/2)[1/(H0-i sigma) + 1/(H0+i sigma)].

The pressure correction is Pcal=F0^2 V(Phi^2). Its size and Lipschitz estimates use the same analytic amplitude F0 and all coupling factors. The restored datum is not replaced. Extremely small amplitude factors are kept in arbitrary-precision arithmetic.

The receipt lists20 split terms, their polynomial powers in Phi/Psi, coefficient bounds, size bounds, and Lipschitz bounds. For a positive monomial majorant c Phi^p Psi^q, derivative bounds in both variables give a conservative Lipschitz estimate in the coupled sum norm. Angular terms are multiplied by the full resolvent bound; axial terms use J1 directly.

## Actual parameter gate remains unproved

At epsilon=1e-36 the resulting logarithmic upper bounds are:

| Quantity | Upper logarithm |
| --- | --- |
| epsilon times size of N | 407.4556429016299 |
| epsilon times Lipschitz constant of N | 396.6910940233444 |

Neither upper bound proves the paper's sufficient half-radius/half-Lipschitz gate. This is not a proof that the actual map is noncontractive, that the analytic core does not exist, or that the paper fails. It quantifies the loss in the present conservative estimates. The angular cross-swirl term is the largest inverted-remainder size bound; multiplying it by the coarse linear resolvent estimate dominates the coupled gate.

Independent algebra checks compare the complete original(8.50) expressions with the split forms at30 parameter pairs, including nonzero delta, epsilon, pressure and swirl inputs. These checks verify the split identities; they do not prove contraction by sampling.

## Next work

1. Replace the coarse angular resolvent estimate by a tighter validated inverse or a coefficient-preconditioned bound preserving its factorial structure.
2. Bound the residual at the actual leading pair separately from the Lipschitz action on its unit ball; exploit analytic identities without discarding couplings.
3. Sharpen fixed-data multipliers near the H0 root, keeping the pressure-induced Psi profile and its correlation with H0. Do not treat these as independent fitted fields.
4. Recompute both map gates. If proof requires different Lambda/Cstar, propagate that change through the shared source schedule and pressure construction before claiming a compatible core.
5. Only after contraction is certified, apply the conditional infinite radial tail factors and recover the stress-free core entrance and collar matching.

Original parameter derivation, full core source errors, infinite core remainder, smooth matching, global cone, temporal recursion, oscillatory correction and full residual validation remain open.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_nonlinear_map_majorant.py
