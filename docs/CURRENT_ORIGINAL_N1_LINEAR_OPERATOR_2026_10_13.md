# Original n1 physical linear repair operator — 2026-10-13

The positive-order repair now has a concrete same-family physical angular amplitude, five newly normalized compact bumps, the axial two-by-two and angular three-by-three matrices, and a source-function C1 inverse-action API. Their original radius factors and finite positive mu are retained. The new local operator gate is checked; the full n1 repair and actual temporal recursion remain incomplete.

This consumes the [complete pre-repair five-moment values](CURRENT_ORIGINAL_N1_FIVE_DEFECT_VALUES_2026_10_13.md). Those values remain enclosures. They are not selected as affine midpoint functions or passed to the coefficient API as scalars. A complete genuine signed defect-function adapter is still required before this operator can install the actual alpha/beta functions and repaired own M/V/P.

## Actual background and support

The accepted C3 repair-band graph is an immutable prefix. Its amplitude leaves are the original ordinary O3 endpoint E rows, including the genuine first, second and third Z derivatives. They are matched to the checked power-to-Rp source binding. The original common-source convention is E=Utheta/Pstar. Consequently,

```text
alpha = 1/2 + mu
Astar(Z) = Pstar*A_rc(Z)*(2Rc)^alpha
Utheta0(R,Z) = Astar(Z)*(2R)^(-alpha)
```

The physical constant is applied once. Astar positivity is inherited from the actual positive amplitude enclosure on each of the four whole axial cells. Log Rc is restored from the original microscopic origin and checked independently against the accepted native Rp-frame expression. Numerical range endpoints do not define the amplitude functions.

The old compact corrections occupy log(R/Rc)/log(2) in [7,9]/40, [19,21]/40 and [31,33]/40. Their actual raw-beta graph arguments are checked, rather than inferred from a geometric drawing. The new supports are [119,121]/400, [129,131]/400, [139,141]/400, [149,151]/400 and [159,161]/400. They are strictly disjoint and lie within the reserved gap [11,17]/40. The axial pair uses the first and last supports; the angular triple uses the middle three.

These narrower supports are an autonomous reconstruction choice, not the paper's original prescribed bump set. Both old local corrections and all their local derivatives vanish on the new supports. The leading axial velocity there is zero and the leading swirl retains its original power law. Nonlocal histories, pressure and density are retained.

## Normalization and original matrices

Let s=log(R/Rc), ell=log(2)/400, c_j=j*log(2)/40 for j=12,...,16 and beta(t)=exp(-1/(1-t^2)) for |t|<1, zero otherwise. Define

```text
J(k) = integral_-1^1 beta(t)*exp(k*ell*t) dt
J_j = ell*exp(c_j/2)*J(1/2)
b_j(R) = sqrt(2/Rc)*beta((s-c_j)/ell)/J_j
integral b_j(R)*dR/sqrt(2R) = 1
H(k) = J(k+1/2)/J(1/2)
M_j(k) = exp(k*c_j)*H(k)
```

The physical axial rows are sqrt(2Rc)*M_j(1/2) and sqrt(2Rc)*(2Rc)^(-mu)*M_j(1/2-mu). The physical angular rows have factors 2Rc, (2Rc)^(-1-mu), (2Rc)^(-mu), multiplying M_j(1), M_j(-1-mu), M_j(-mu), respectively.

The axial second row is evaluated in divided form D_j=(M_j(1/2-mu)-M_j(1/2))/mu using the entire exprel kernel. No subtraction of independently boxed near-equal rows is used. Its determinant has an exact translated-bump formula with a strictly positive negative-of-determinant. The equally spaced angular supports give a fresh translated exponential Vandermonde determinant. Both identities are independently checked.

The actual original normalized determinant bounds are approximately:

```text
-det(axial divided matrix) = 0.088345841987096972
 det(angular normalized matrix) = 1.0407806602793233e-5
```

The positive normalizer J(1/2) is approximately 0.443993842518479277. Directed Taylor integration retains whole-cell omitted-derivative errors and nonzero endpoint-strip allowances. Matrix and determinant interval widths are about 1e-18 or smaller. The original axial determinant still includes its true mu factor; the divided form does not replace mu with zero.

## Coefficient function API and scope

The physical equations are unchanged:

```text
Bz*alpha = -[m11, m14/Astar]
Btheta*beta = -[m12, m13/(2*Astar), -m15/Astar]
```

`CurrentOriginalN1ReservedLinearOperator.coefficient_functions(moments)` requires five same-graph C1Function handles. It applies the physical row factors, preserves the exact signed axial cancellation before division by mu, includes Astar_Z in the quotient derivative and applies the fresh inverse matrices. Scalars, magnitude caps and interval midpoint values are rejected as functions. The API is a defined linear operator; actual complete defect handles have not yet been bound to it.

The installed gate is only `current_original_n1_reserved_normalized_linear_repair_operator_installed`. Full pressure C1 accuracy, actual selected n1 coefficient functions, repaired own M/V/P, five terminal zeros, temporal recursion, global stress and final NS residual validation remain false.

## Reproduction and focused checks

Run `python scripts/run_original_n1_reserved_linear_operator.py`. It reads checked archives, preserves the accepted background prefix and never reruns ancestor producers/checkers. The source quartet is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator` with .py, .json.gz, _check.py and _check.json.

The checker independently integrates the physical bump measures with Rc=3.7, Astar=2.3 and mu=.08; checks all five normalizations and both physical matrices; verifies the fresh determinant identities; and exercises the production C1 inverse action against an independent solve with nonconstant amplitude and nonconstant signed moments. It checks the factor two in m13, the negative physical m15 change and the amplitude quotient derivative. Synthetic fixture functions are never archived or installed as original fields.

## Next agent tasks

- [x] Bind actual physical Astar ordinary Z0-Z3 handles and positivity on four whole axial cells.
- [x] Check original log Rc against the accepted native-frame identity.
- [x] Install fresh disjoint n1 supports in the exact old-correction gap.
- [x] Build physical dR/sqrt(2R) normalizers with nonzero quadrature errors.
- [x] Construct the correct positive-order axial pair and angular triple matrices.
- [x] Prove fresh translated-bump determinant identities and positive normalized bounds.
- [x] Implement the stable divided axial row and physical C1 coefficient inverse API.
- [x] Independently check nonunit physical weights, five-moment signs and C1 derivatives.
- [ ] Bind complete signed defect functions and genuine Z derivatives; do not select current interval values as functions.
- [ ] Complete full pressure C1 source enclosures. The existing pre-Rp q-Z2 magnitude archive is too wide to close this gate; use correlated signed derivative integrals or a sharper actual q-Z2 source.
- [x] Consume the checked [outer-pulse C1 source and local derivative grid](CURRENT_ORIGINAL_N1_PULSE_C1_2026_10_13.md); full pressure derivative composition remains open.
- [ ] Extend finite lower endpoint recovery using genuine U_yZ and U_ZZ rows.
- [ ] Install actual alpha(Z)/beta(Z), add the corrections to F1/Uz1 and recover that family's own M/V/P.
- [ ] Prove all five terminal identities as functions and check repaired support/joins/axis/heat compatibility.
- [ ] Continue genuine n2 and temporal recursion, admissible stress, oscillatory cancellation and independent corrected NS validation.

Full original goal **ACTIVE / INCOMPLETE**.
