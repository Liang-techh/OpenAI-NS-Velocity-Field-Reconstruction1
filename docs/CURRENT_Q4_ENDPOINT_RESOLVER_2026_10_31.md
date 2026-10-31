# Authenticated Q4 endpoint partial resolver — 2026-10-31

A new resolver consumes the saved physical U_y endpoint ports and the actual Q4 gzip archive. It verifies the raw archive SHA-256 before decoding, checks that the selected parameter is the archived Z-independent sc source parameter, binds the chosen exact sc and sc/2 dyadics, and follows the source expressions without calling any earlier builder or evaluator. The port receipt is separately authenticated by the single-entry runner; the general resolver requires its caller to authenticate the ports.

The resolver implements exact rational and dyadic arithmetic, sums, products, positive quotients with a known positive denominator, exact comparisons, three-valued logical OR and lazy flat branches. It represents scalars as a rational coefficient times `2**exponent`, without expanding the extremely large exponent. Its environment accepts exact dyadics for declared source-variable nodes only; attempts to replace an arithmetic root are rejected. Operations it cannot resolve remain explicit opaque source references; they are not interpreted as zero, interval midpoints, bounds or field values. Source expressions with unresolved factors retain their original domain assumptions.

The actual authenticated run returned **unresolved** for all three physical U_y ordinary Z0/Z1/Z2 rows. In those respective rows it did, however, select both original flat branches: B_phi nodes [99405,99406,99407] and native slow-B derivative nodes [128202,128239,128295]. The corresponding active bodies [99399,99400,99401] and [128201,128238,128294] were not followed through those branch decisions. It retained the physical source expressions and unresolved background, rather than treating the whole velocity as zero. Traversal visited 111, 201 and 306 distinct nodes, respectively.

The remaining opaque operations include analytic unary expressions (including Pstar), source-variable Z, source parameters and the common N expression, negative, function substitution, definite integrals and the leading-function recipe. Thus this entry provides a real authenticated **partial** source-expression resolver, not a numerical Q4 velocity or a certified error bound. No free Z value was supplied. The selected sc is a new explicit instance; the historical assignment and corrected source owner remain unproved. No old errors, signed terms or histories were removed.

## Validation and next work

The new implementation and all inputs were frozen before its one permitted entry. A fresh synthetic fixture passed 27 checks, including exact evaluation of a small closed expression as 19/64, lazy-branch short circuit, missing or unsupported source operations, invalid endpoint bindings and rejection of an environment that overrides an arithmetic root. The actual result was saved in full before inspection. Its raw source hash was verified before gzip decoding.

Next, implement fresh versioned support for negative and capture-safe function substitution and preserve analytic unary nodes as exact symbolic operations, while retaining unresolved source parameters, integrals and domain obligations. Use a fresh fixture, freeze and enter that new consumer once; do not reenter this one or any ancestor family. Then supply legitimate actual source parameters and source functions, with independent numerical/error certification, and bind corrected physical U_y and leading m-history at the same endpoint. Preserve the 212 cumulative terms, six signed contributions, the reserved contribution once, the independent negative complete pulse, ordinary Z1 and all prior errors. Global five-moment closure and the full NS goal remain open.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | ee8ac7e6a410ab565f89f3e85b2d7040389ffceb2b97072328f10455918e6f17 |
| Frozen source manifest | aaa7b1c1bbff5612797bc9fa97f4229c31960f3a8fe547e0f8fadf9c053f8d13 |
| Exclusive entry | 74e5990edcbecb7507bd6c8331e5c2a21eba551c572983b2dca2d13ddef47be7 |
| Complete result | 573e607e6812943de6d5fae06e0082c14faf0f10fc1e630692b413f7bb6c33db |
| Resolver | e8ff5870edda026b61ebd9908eb453663eed86c8b0d6ce895dcae9be39ec437c |
| Fixture | cc2d315e78d4a055dd1e68096645b9e97fdfa0b506df9a9b57899613107fb442 |
| Runner | 9338534f9a433d67c09e4643c7784b4e27399b743f959f47be5a65131f946feb |
