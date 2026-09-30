# Separated five-defect response for Section10 moment repair

The fixed five-bump map now has a finite formal inverse response to five independent centered defects. This implements the paper's local contraction expansion without merging tiny d3/d5 into larger cancelling scalar coefficients.

Write h=sum h_n, with n denoting degree in the five defect variables, not a temporal coefficient index. Then

    h_1 = -A^(-1) d
    h_n = -A^(-1) sum_(i+j=n) Q_Z(h_i,h_j), n>=2.

The response stores each five-tuple monomial separately. At degree 3 there are 55 monomials and each has a five-vector of bump-coefficient contributions. Coefficients may themselves be pressure/width jets or first-Z duals. Original centered source labels remain separate in the actual receipt's linear response; nonlinear terms currently use the serialized aggregate defect jets, rather than expanding every source-stage crossproduct.

## Evidence

The independent fixture applies the actual five-bump matrix and bilinear map coefficient by coefficient, then compares evaluated coefficients against an independent scalar fixed-point iteration on a resolved small-defect input. At degree 3 the coefficientwise residual is 2.01e-87, scalar coefficient relative error 1.67e-22, and independently differenced first-Z relative error 1.48e-52.

The actual check consumes the committed 120-digit centered defect values/tangents at Z=.3, pressure order 9 and width order 2; it does not rebuild or change the core or pressure datum. It requires separately nonzero d3 and d5 response vectors, and retains all 69 original labels in the linear response.

## Limits

Finite response degree is not exact five-moment closure, a convergent-series certificate, or the n-dependent temporal recursion. Source and quadrature errors remain unenclosed; higher defect-degree and finite pressure/width remainders remain open. A rounded aggregate h or aggregate moment residual cannot establish closure relative to the actual tiny rows. Next steps are a source-aware remainder/convergence treatment, corrected-field partial-moment recovery, uniform analytic input bounds, and the relaxed cone audit before joining the supplied outer construction.

The actual serialized-input run completed with 55 separated monomial responses and 69 labeled linear source responses. Both direct tiny-row responses remain nonzero before any aggregate coefficient sum. No corrected global field has been installed by this receipt.
