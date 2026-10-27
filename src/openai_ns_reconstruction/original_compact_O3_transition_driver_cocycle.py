"""Standalone exact algebra for the finite-``N`` exponential driver cocycle.

This file is intentionally a new, unexecuted algebra kernel.  It has no
project imports, no source-file reads, and no calls into any earlier proof or
checker family.  The only symbolic variables used by the sparse witness are
the ordinary jet names ``En0..En3``, ``Ec0..Ec3``, ``Ac0..Ac3``, ``D0..D3``,
``eps``, ``expAc``, and ``expD``.

The certificate is algebraic.  It records the exact cocycle identity and its
ordinary ``Z``-derivative recipes, while source binding remains a caller-owned
gate.  In particular, this module never evaluates a source, a numerical
exponential, a finite ``N`` value, or a source gap.  It also does not assert a
quantitative bound or a global closure statement.
"""

__all__ = ("derive_exponential_cocycle",)


class _ExponentialCocyclePolynomial:
    """Small immutable-by-convention sparse polynomial over integer scalars.

    A monomial is a sorted tuple of symbolic atom names.  No numerical field
    values are introduced: integer coefficients only record exact Leibniz-rule
    multiplicities.  Every operation creates a fresh instance.
    """

    __slots__ = ("terms",)

    def __init__(self, terms=None):
        normalized = {}
        if terms:
            for monomial, coefficient in terms.items():
                if coefficient:
                    key = tuple(sorted(monomial))
                    normalized[key] = normalized.get(key, 0) + coefficient
        self.terms = {
            monomial: coefficient
            for monomial, coefficient in normalized.items()
            if coefficient
        }

    @classmethod
    def atom(cls, name):
        return cls({(name,): 1})

    @classmethod
    def scalar(cls, value):
        return cls({(): value}) if value else cls()

    def _coerce(self, other):
        if isinstance(other, _ExponentialCocyclePolynomial):
            return other
        if isinstance(other, int):
            return self.scalar(other)
        raise TypeError("exponential-cocycle algebra accepts integer scalars only")

    def __add__(self, other):
        right = self._coerce(other)
        merged = dict(self.terms)
        for monomial, coefficient in right.terms.items():
            merged[monomial] = merged.get(monomial, 0) + coefficient
        return self.__class__(merged)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return self.__class__(
            {monomial: -coefficient for monomial, coefficient in self.terms.items()}
        )

    def __sub__(self, other):
        return self.__add__(-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other).__sub__(self)

    def __mul__(self, other):
        right = self._coerce(other)
        product = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in right.terms.items():
                monomial = tuple(sorted(left_monomial + right_monomial))
                product[monomial] = product.get(monomial, 0) + (
                    left_coefficient * right_coefficient
                )
        return self.__class__(product)

    def __rmul__(self, other):
        return self.__mul__(other)

    def is_zero(self):
        return not self.terms

    def as_rows(self):
        """Return deterministic coefficient rows for an audit-friendly witness."""

        rows = []
        for monomial in sorted(self.terms):
            rows.append(
                {
                    "monomial": monomial,
                    "coefficient": self.terms[monomial],
                }
            )
        return tuple(rows)


def _exponential_cocycle_derivative_atom(atom):
    """Return the exact first ``Z`` derivative of one symbolic atom.

    The exponential atoms use the chain rule
    ``(exp(eps X))_Z = eps X_Z exp(eps X)``.  Jet atoms stop at ordinary
    order three; this is deliberate, because no fourth derivative is part of
    the contract.
    """

    atom_poly = _ExponentialCocyclePolynomial.atom
    if atom in ("En0", "En1", "En2"):
        return atom_poly("En" + str(int(atom[-1]) + 1))
    if atom in ("Ec0", "Ec1", "Ec2"):
        return atom_poly("Ec" + str(int(atom[-1]) + 1))
    if atom in ("Ac0", "Ac1", "Ac2"):
        return atom_poly("Ac" + str(int(atom[-1]) + 1))
    if atom in ("D0", "D1", "D2"):
        return atom_poly("D" + str(int(atom[-1]) + 1))
    if atom == "expAc":
        return atom_poly("eps") * atom_poly("Ac1") * atom_poly("expAc")
    if atom == "expD":
        return atom_poly("eps") * atom_poly("D1") * atom_poly("expD")
    return _ExponentialCocyclePolynomial()


def _exponential_cocycle_z_derivative(polynomial):
    """Apply one ordinary ``Z`` derivative by the exact product rule."""

    result = _ExponentialCocyclePolynomial()
    for monomial, coefficient in polynomial.terms.items():
        for position, atom in enumerate(monomial):
            atom_derivative = _exponential_cocycle_derivative_atom(atom)
            if atom_derivative.is_zero():
                continue
            remainder = monomial[:position] + monomial[position + 1 :]
            remainder_poly = _ExponentialCocyclePolynomial({remainder: 1})
            result = result + coefficient * remainder_poly * atom_derivative
    return result


def _exponential_cocycle_gate_value(contract, names):
    """Read a source-gate field without evaluating source objects.

    The root wrapper may place gates at the top level or under
    ``source_gates``.  Only plain mapping lookup is performed; callable
    source values are treated as opaque presence evidence and never called.
    """

    if not isinstance(contract, dict):
        return False, None, None
    containers = [(contract, "")]
    nested = contract.get("source_gates")
    if isinstance(nested, dict):
        containers.append((nested, "source_gates."))
    for container, prefix in containers:
        for name in names:
            if name in container:
                return True, container[name], prefix + name
        shared_n = container.get("shared_symbolic_N")
        if isinstance(shared_n, dict):
            if any(name.endswith("_node") for name in names) and "node" in shared_n:
                return True, shared_n["node"], prefix + "shared_symbolic_N.node"
            if any(name.endswith("_descriptor") for name in names) and "descriptor" in shared_n:
                return (
                    True,
                    shared_n["descriptor"],
                    prefix + "shared_symbolic_N.descriptor",
                )
    return False, None, None


def _exponential_cocycle_gate_record(contract, label, names, expectation):
    present, value, key = _exponential_cocycle_gate_value(contract, names)
    if expectation == "present":
        passed = present and value is not None
    elif expectation == "true":
        passed = present and value is True
    elif expectation == "node":
        passed = present and type(value) is int and value == 18
    elif expectation == "descriptor":
        passed = present and type(value) is str and value == "J1358356628656378313"
    else:
        passed = False
    return {
        "label": label,
        "required": True,
        "passed": passed,
        "received_key": key,
        "reason": (
            "opaque source binding present"
            if expectation == "present" and passed
            else "exact gate value present"
            if passed
            else "caller must provide the required source gate; no source value was evaluated"
        ),
    }


def _exponential_cocycle_source_gate(contract):
    """Build the conditional source-binding gate record."""

    records = (
        _exponential_cocycle_gate_record(
            contract,
            "actual_source_bindings",
            (
                "actual_source_bindings",
                "actual_source_binding",
                "source_bindings",
                "source_binding",
            ),
            "present",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "shared_symbolic_N_node",
            ("shared_symbolic_N_node", "shared_N_node"),
            "node",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "shared_symbolic_N_descriptor",
            ("shared_symbolic_N_descriptor", "shared_N_descriptor"),
            "descriptor",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "ordinary_Z_independent_N",
            ("ordinary_Z_independent_N", "ordinary_z_independent_N"),
            "true",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "ordinary_derivative_vs_Taylor_correct",
            (
                "ordinary_derivative_vs_Taylor_correct",
                "ordinary_derivative_is_not_Taylor_coefficient",
            ),
            "true",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "parent_normalization_data_only",
            ("parent_normalization_data_only", "normalization_data_only"),
            "true",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "finite_N_producer_C3_exprel_native_expm1_bound",
            (
                "finite_N_producer_C3_exprel_native_expm1_bound",
                "finite_N_producer_semantics_C3_exprel_native_expm1_bound",
                "finite_N_producer_C3_exprel_native_expm1_bound_validated",
                "finite_N_producer_semantics",
            ),
            "true",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "actual_C3_exprel_producer_clauses",
            (
                "actual_C3_exprel_producer_clauses",
                "actual_producer_clauses",
                "c3_exprel_producer_clauses",
            ),
            "present",
        ),
        _exponential_cocycle_gate_record(
            contract,
            "actual_C3_exprel_producer_inputs",
            (
                "actual_C3_exprel_producer_inputs",
                "actual_producer_inputs",
                "c3_exprel_producer_inputs",
                "actual_input_ref_arrays",
                "input_ref_arrays",
            ),
            "present",
        ),
    )
    return {
        "records": records,
        "all_passed": all(record["passed"] for record in records),
        "source_evaluated": False,
        "source_gap_zero_proven": False,
        "status": "passed" if all(record["passed"] for record in records) else "pending",
    }


def _exponential_cocycle_polynomial_witness():
    """Construct the four exact scaled residual rows.

    Write ``P_X = X*exprel(eps*X)`` and ``E_X = exp(eps*X)`` only in the
    accompanying recipe.  The polynomial witness uses the division-free
    relations ``eps*P_X = E_X - 1`` and
    ``eps*P_(Ac+D) = expAc*expD - 1``.  Thus the base residual is

    ``En*(expAc*expD - 1) - Ec*(expAc - 1)
    - (En - Ec)*(expAc - 1) - En*expAc*(expD - 1)``.

    ``eps`` is a constant under the ordinary ``Z`` operator.  Four successive
    applications of the fresh product/chain-rule operator give the rows for
    orders zero through three.
    """

    atom = _ExponentialCocyclePolynomial.atom
    en = atom("En0")
    ec = atom("Ec0")
    exp_ac = atom("expAc")
    exp_d = atom("expD")
    h = en - ec
    sum_increment = exp_ac * exp_d - 1
    ac_increment = exp_ac - 1
    d_increment = exp_d - 1
    lhs = en * sum_increment - ec * ac_increment
    rhs = h * ac_increment + en * exp_ac * d_increment

    residual_rows = []
    lhs_rows = []
    rhs_rows = []
    current_lhs = lhs
    current_rhs = rhs
    for order in (0, 1, 2, 3):
        residual = current_lhs - current_rhs
        lhs_rows.append(
            {
                "order": order,
                "operator": "D_Z^" + str(order),
                "coefficient_rows": current_lhs.as_rows(),
                "is_nonzero": not current_lhs.is_zero(),
            }
        )
        rhs_rows.append(
            {
                "order": order,
                "operator": "D_Z^" + str(order),
                "coefficient_rows": current_rhs.as_rows(),
                "is_nonzero": not current_rhs.is_zero(),
            }
        )
        residual_rows.append(
            {
                "order": order,
                "operator": "D_Z^" + str(order),
                "residual_polynomial": residual.as_rows(),
                "is_zero": residual.is_zero(),
            }
        )
        if order != 3:
            current_lhs = _exponential_cocycle_z_derivative(current_lhs)
            current_rhs = _exponential_cocycle_z_derivative(current_rhs)
    return {
        "variables": (
            "En0",
            "En1",
            "En2",
            "En3",
            "Ec0",
            "Ec1",
            "Ec2",
            "Ec3",
            "Ac0",
            "Ac1",
            "Ac2",
            "Ac3",
            "D0",
            "D1",
            "D2",
            "D3",
            "eps",
            "expAc",
            "expD",
        ),
        "residual_definition": (
            "eps*R = En0*(expAc*expD - 1) - Ec0*(expAc - 1) "
            "- (En0 - Ec0)*(expAc - 1) - En0*expAc*(expD - 1)"
        ),
        "lhs_rows": tuple(lhs_rows),
        "rhs_rows": tuple(rhs_rows),
        "lhs_coefficient_rows": tuple(row["coefficient_rows"] for row in lhs_rows),
        "rhs_coefficient_rows": tuple(row["coefficient_rows"] for row in rhs_rows),
        "residual_rows": tuple(residual_rows),
        "rows": tuple(residual_rows),
        "row_count": len(residual_rows),
        "lhs_rows_nonzero": all(row["is_nonzero"] for row in lhs_rows),
        "rhs_rows_nonzero": all(row["is_nonzero"] for row in rhs_rows),
        "all_rows_zero": all(row["is_zero"] for row in residual_rows),
        "operator": {
            "kind": "ordinary_Z_differential_algebra",
            "product_rule": "D_Z(U*V)=D_Z(U)*V+U*D_Z(V)",
            "expAc_rule": "D_Z(expAc)=eps*Ac1*expAc",
            "expD_rule": "D_Z(expD)=eps*D1*expD",
            "eps_rule": "D_Z(eps)=0",
            "highest_order": 3,
            "ordinary_derivative_not_Taylor_coefficient": True,
        },
    }


def _exponential_cocycle_c3_producer_intertwiner():
    """Check the C3 ``exprel`` producer against an independent derivative.

    The producer side is written explicitly for orders zero through three.
    Multiplying by ``eps`` replaces every ``A0*exprel(eps*A0)`` by the
    division-free increment ``expAc - 1``.  The other side is obtained by
    differentiating the independent base polynomial
    ``E0*(expAc - 1)``.  They are compared only as formal sparse polynomials;
    no source producer or input reference is evaluated here.
    """

    atom = _ExponentialCocyclePolynomial.atom
    eps = atom("eps")
    exp_ac = atom("expAc")
    increment = exp_ac - 1
    e0 = atom("En0")
    e1 = atom("En1")
    e2 = atom("En2")
    e3 = atom("En3")
    a1 = atom("Ac1")
    a2 = atom("Ac2")
    a3 = atom("Ac3")

    explicit_scaled = (
        e0 * increment,
        e1 * increment + eps * e0 * a1 * exp_ac,
        (
            e2 * increment
            + 2 * eps * e1 * a1 * exp_ac
            + eps * e0 * exp_ac * (a2 + eps * a1 * a1)
        ),
        (
            e3 * increment
            + 3 * eps * e2 * a1 * exp_ac
            + 3 * eps * e1 * exp_ac * (a2 + eps * a1 * a1)
            + eps
            * e0
            * exp_ac
            * (a3 + 3 * eps * a1 * a2 + eps * eps * a1 * a1 * a1)
        ),
    )

    independent_rows = []
    independent = e0 * increment
    for order in (0, 1, 2, 3):
        independent_rows.append(independent)
        if order != 3:
            independent = _exponential_cocycle_z_derivative(independent)

    checks = []
    for order in (0, 1, 2, 3):
        explicit = explicit_scaled[order]
        independent_row = independent_rows[order]
        residual = explicit - independent_row
        checks.append(
            {
                "order": order,
                "scaled_explicit_source_producer": explicit.as_rows(),
                "independent_D_Z_row": independent_row.as_rows(),
                "residual_polynomial": residual.as_rows(),
                "is_zero": residual.is_zero(),
            }
        )

    return {
        "kind": "C3_exprel_producer_intertwiner",
        "orders": (0, 1, 2, 3),
        "producer_formulas": (
            "F0 = E*Ac*exprel(eps*Ac)",
            "F1 = E1*Ac*exprel(eps*Ac) + E*A1*expAc",
            "F2 = E2*Ac*exprel(eps*Ac) + 2*E1*A1*expAc "
            "+ E*expAc*(A2 + eps*A1^2)",
            "F3 = E3*Ac*exprel(eps*Ac) + 3*E2*A1*expAc "
            "+ 3*E1*expAc*(A2 + eps*A1^2) "
            "+ E*expAc*(A3 + 3*eps*A1*A2 + eps^2*A1^3)",
        ),
        "indexed_producer_formulas": (
            "F0 = E0*A0*exprel(eps*A0)",
            "F1 = E1*A0*exprel(eps*A0) + E0*A1*exp(eps*A0)",
            "F2 = E2*A0*exprel(eps*A0) + 2*E1*A1*exp(eps*A0) "
            "+ E0*exp(eps*A0)*(A2 + eps*A1^2)",
            "F3 = E3*A0*exprel(eps*A0) + 3*E2*A1*exp(eps*A0) "
            "+ 3*E1*exp(eps*A0)*(A2 + eps*A1^2) "
            "+ E0*exp(eps*A0)*(A3 + 3*eps*A1*A2 + eps^2*A1^3)",
        ),
        "scaled_source_producer_formulas": (
            "eps*F0 = E0*(expAc - 1)",
            "eps*F1 = E1*(expAc - 1) + eps*E0*A1*expAc",
            "eps*F2 = E2*(expAc - 1) + 2*eps*E1*A1*expAc "
            "+ eps*E0*expAc*(A2 + eps*A1^2)",
            "eps*F3 = E3*(expAc - 1) + 3*eps*E2*A1*expAc "
            "+ 3*eps*E1*expAc*(A2 + eps*A1^2) "
            "+ eps*E0*expAc*(A3 + 3*eps*A1*A2 + eps^2*A1^3)",
        ),
        "scaled_source_producer_polynomials": tuple(
            polynomial.as_rows() for polynomial in explicit_scaled
        ),
        "independent_derivative_polynomials": tuple(
            polynomial.as_rows() for polynomial in independent_rows
        ),
        "base_independent_expression": "D_Z^j[E0*(expAc - 1)]",
        "checks": tuple(checks),
        "check_count": len(checks),
        "all_checks_zero": all(check["is_zero"] for check in checks),
        "ordinary_derivative_not_Taylor_coefficient": True,
        "source_evaluated": False,
        "division_free_scaling": True,
    }


def _exponential_cocycle_recipes():
    """Return source and derivative recipes without source evaluation."""

    return {
        "source_formula": {
            "epsilon": "eps = 1/N > 0; D_Z(eps)=0",
            "exprel": "exprel(q) = expm1(q)/q, continued analytically with exprel(0)=1",
            "F": "F(E,A) = E*A*exprel(eps*A) = E*expm1(eps*A)/eps",
            "delta_F": "DeltaF = F(En,An) - F(Ec,Ac)",
            "cocycle": (
                "DeltaF = H*Ac*exprel(eps*Ac) "
                "+ En*exp(eps*Ac)*D*exprel(eps*D)"
            ),
            "jet_definitions": (
                "D = An - Ac",
                "H = En - Ec",
                "En, Ec, An, Ac are independent C3 ordinary Z jets",
            ),
            "domains": {
                "source_x": "[0,1]",
                "ordinary_Z": "[-1,1]",
                "regularity": "arbitrary C3 functions",
            },
            "assumptions_not_used": (
                "no E equality",
                "no A equality",
                "no phase equality",
                "no P2 equality",
                "no zero assumption",
            ),
        },
        "exponential_split": {
            "identity": "exp(eps*(Ac+D)) = exp(eps*Ac)*exp(eps*D)",
            "expm1_relation": "expm1(q) = exp(q) - 1",
            "division_free": (
                "eps*(X*exprel(eps*X)) = exp(eps*X) - 1"
            ),
            "scaled_cocycle": (
                "eps*DeltaF = (En - Ec)*(expAc - 1) "
                "+ En*expAc*(expD - 1)"
            ),
        },
        "ordinary_Z_expression_recipes": {
            "derivative_kind": "ordinary D_Z^j, not a Taylor coefficient",
            "orders": (0, 1, 2, 3),
            "exp_chain_rule": (
                "E_X[0] = exp(eps*X)",
                "E_X[j+1] = D_Z^j(eps*X_Z*E_X[0])",
                "E_X[1] = eps*X1*E_X[0]",
                "E_X[2] = (eps*X2 + eps^2*X1^2)*E_X[0]",
                "E_X[3] = (eps*X3 + 3*eps^2*X1*X2 + eps^3*X1^3)*E_X[0]",
            ),
            "entire_increment_chain_rule": (
                "P_X = X*exprel(eps*X)",
                "eps*P_X = E_X[0] - 1",
                "P_X[0] = X*exprel(eps*X)",
                "P_X[j] = sum_{k=0}^{j-1} binom(j-1,k)*X[k+1]*E_X[j-1-k] for j>=1",
            ),
            "product_rule": (
                "D_Z^j(U*V) = sum_{k=0}^j binom(j,k)*U[k]*V[j-k]"
            ),
            "delta_F_recipe": (
                "DeltaF[j] = sum_{k=0}^j binom(j,k)*En[k]*P_An[j-k] "
                "- sum_{k=0}^j binom(j,k)*Ec[k]*P_Ac[j-k]"
            ),
            "cocycle_recipe": (
                "DeltaF[j] = sum_{k=0}^j binom(j,k)*H[k]*P_Ac[j-k] "
                "+ sum_{a+b+c=j} j!/(a!*b!*c!)*En[a]*E_Ac[b]*P_D[c]"
            ),
            "ordinary_derivative_normalization": "no division by j!",
        },
        "removable_zero_gap": {
            "statement": (
                "exprel is entire at zero; P_X=X*exprel(eps*X) is used "
                "without dividing by X, D, or Ac"
            ),
            "no_gap_division": True,
            "no_ac_division": True,
            "native_expm1_semantics_left_to_source_gate": True,
        },
    }


def _exponential_cocycle_operator_recipes():
    """Describe four reference-array outputs without touching source arrays."""

    recipes = []
    for order in (0, 1, 2, 3):
        recipes.append(
            {
                "order": order,
                "output": "DeltaF_" + str(order),
                "expression": (
                    "D_Z^" + str(order)
                    + "[F(En,An)-F(Ec,Ac)] evaluated on actual input references"
                ),
                "rhs_expression": (
                    "D_Z^" + str(order)
                    + "[H*Ac*exprel(eps*Ac) + En*exp(eps*Ac)*D*exprel(eps*D)]"
                ),
                "input_reference_convention": (
                    "En[" + str(order) + "], Ec[" + str(order) + "], "
                    "An[" + str(order) + "], Ac[" + str(order) + "], "
                    "D[" + str(order) + "], H[" + str(order) + "]"
                ),
                "source_evaluated": False,
            }
        )
    return {
        "kind": "operator_output_recipe",
        "array_binding": "actual_input_ref_arrays supplied by the root wrapper",
        "recipes": tuple(recipes),
        "count": len(recipes),
        "no_source_evaluation": True,
    }


def derive_exponential_cocycle(contract, journal=None):
    """Return the conditional exact finite-``N`` exponential cocycle certificate.

    ``contract`` is only a source-gate mapping.  The returned certificate is
    exact at the formal-algebra level even when the caller has not supplied
    source bindings.  Once the root wrapper supplies the required gate fields,
    ``source_gate.status`` reports whether all binding metadata is present;
    this function still never evaluates a source or proves a source gap is
    zero.

    ``journal`` is intentionally caller-owned and is not mutated or read as a
    source.  It is recorded only as a presence flag.
    """

    source_gate = _exponential_cocycle_source_gate(contract)
    witness = _exponential_cocycle_polynomial_witness()
    c3_producer = _exponential_cocycle_c3_producer_intertwiner()
    recipes = _exponential_cocycle_recipes()
    operator = _exponential_cocycle_operator_recipes()
    return {
        "family": "NS O3_transition finite-N exponential driver cocycle",
        "certificate_kind": "conditional_exact_algebra_certificate",
        "status": "conditional_exact_certificate",
        "source_gate": source_gate,
        "shared_symbolic_N": {
            "node": 18,
            "descriptor": "J1358356628656378313",
            "ordinary_Z_independent": True,
        },
        "formal_assumptions": (
            "eps=1/N>0",
            "Z-independent eps",
            "En, Ec, An, Ac independent C3 ordinary jets through order 3",
            "D=An-Ac",
            "H=En-Ec",
        ),
        "algebra": witness,
        "cocycle_primary4": witness,
        "c3_exprel_producer4checks": c3_producer,
        "recipes": recipes,
        "operator_output": operator,
        "journal": {
            "provided": journal is not None,
            "mutated": False,
            "source_evaluated": False,
        },
        "scope_exclusions": (
            "no actual source binding or source-gap-zero proof",
            "no numerical field values",
            "no numerical exponential or approximation",
            "no N-sized materialization",
            "no quantitative cap or global closure",
            "no old generic five-density proof",
            "no Duhamel kernel",
            "no temporal-n recursion",
        ),
    }

"""Static source wrapper v2 for the finite-N O3 transition driver.

The file is concatenated after the fresh exponential cocycle draft.  It is
deliberately unexecuted while being prepared: only source text, ASTs, hashes,
and explicit literal rows are inspected.  No project module is imported, no
source provider is called, no recursive graph walk is performed, and no
numeric or source-gap value is evaluated.
"""

import ast as _od2_ast
import hashlib as _od2_hashlib


OD2_COMPONENT = "compact_O3_transition_driver_cocycle"
OD2_INPUT_SCHEMA = "compact_O3_transition_driver_cocycle_input_v1"
OD2_RESULT_SCHEMA = "original_compact_O3_transition_driver_cocycle"
OD2_MATH_DRAFT_SHA256 = (
    "655fd80395cc216445f3871e89f2fe480faed2ea7d8adbb41f1350fd1e219a35"
)
OD2_NATIVE_PACKET_SHA256 = (
    "81c49a108764002f2e705aea62e58a5df3c3a52ee99fb57c748b1cc215ae27ad"
)
OD2_C3_PACKET_SHA256 = (
    "acc24cacb22c4040fa1cc906fd55cd8b9ddb2858227509106c3d729dfffb5741"
)
OD2_C3_PRODUCER_PACKET_SHA256 = (
    "3bd8f5e091f73bdce891a3ac567e6bd9cb6f9a8ca1b971fcfee3b4b27c148e2e"
)
OD2_C3_SCALAR_SHA256 = (
    "44c6c4431594f8e25a564a8ef06ac2423cd18b1385627d4012065680620e4e16"
)
OD2_CACHED_RESULT_SHA256 = (
    "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
)
OD2_NATIVE_ARCHIVE_SHA256 = (
    "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
)
OD2_C3_ARCHIVE_SHA256 = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
OD2_N_NODE = 18
OD2_NATIVE_INVN_NODE = 99411
OD2_C3_INVN_NODE = 10287
OD2_J = "1358356628656378313"
OD2_ORDERS = (0, 1, 2, 3)
OD2_FIELDS = ("m", "h", "k", "e", "p")
OD2_FLAGS = (
    "O3_transition_finite_N_driver_source_bound",
    "O3_transition_exponential_driver_cocycle_ordinary_Z3_proved",
    "O3_transition_twenty_signed_density_gap_functions_defined",
)

OD2_JETS_NATIVE = {
    "E": (241552, 241553, 241554, 241555),
    "V": (0, 0, 0, 0),
    "A": (244098, 244099, 244100, 244101),
    "B": (244136, 244137, 244138, 244139),
    "F_raw_r": (244352, 244355, 244359, 244364),
}
OD2_JETS_C3 = {
    "E": (4590, 4591, 9988, 21725),
    "V": (4592, 4593, 9989, 21726),
    "A": (4718, 4725, 10273, 22254),
    "B": (4727, 4732, 10279, 22255),
    "F_N": (4735, 4739, 11937, 24061),
}

# Exact basenames are intentional.  Heuristic substring role matching caused
# the previous wrapper to bind the wrong coefficient source.
OD2_CREATOR_RECORDS = {
    "coefficient_pairs": (
        "lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls.py",
        "527c2f1e7aefa4ff653ceb13d5ed0a3814d52e43a0c9d3b637c8d5a36d89e4bd",
    ),
    "C2_density": (
        "lei_ren_part1_paper_compliant_current_original_C2_density_transport.py",
        "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    ),
    "C3_source": (
        "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
        "de41858de44ff9cea9b5d2b17aeca55eeca2c114969ca4df383ec1345b6badeb",
    ),
    "loop_source": (
        "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
        "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    ),
    "outer": (
        "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
        "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    ),
    "control_connection": (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d",
    ),
    "limitcontrols": (
        "lei_ren_part1_paper_compliant_current_original_C3_limit_controls.py",
        "25176df81492100f6bc9f8a1ada93213a1814fb21a29eb295ca81e1385b06a71",
    ),
    "repair": (
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
        "62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4",
    ),
}

OD2_FAMILY_HASHES = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def _od2_sha_text(value):
    return _od2_hashlib.sha256(value.encode("utf-8")).hexdigest()


def _od2_ast_sha(tree):
    return _od2_hashlib.sha256(
        _od2_ast.dump(tree, include_attributes=False).encode("utf-8")
    ).hexdigest()


def _od2_emit(journal, event, **fields):
    if journal is None:
        return
    payload = {"event": event, **fields}
    if callable(journal):
        journal(payload)
    elif isinstance(journal, list):
        journal.append(payload)


def _od2_spec(source):
    if not isinstance(source, dict) or not isinstance(source.get("input"), dict):
        return {}
    return source["input"]


def _od2_packet(spec, key, expected, blockers):
    packet = spec.get(key)
    if not isinstance(packet, dict):
        blockers.append("input.%s packet missing" % key)
        return {}, {}
    if packet.get("sha256") != expected:
        blockers.append("input.%s SHA changed" % key)
    payload = packet.get("payload")
    if not isinstance(payload, dict):
        blockers.append("input.%s payload is not parsed JSON" % key)
        payload = {}
    return packet, payload


def _od2_row(nodes, node, label, blockers):
    row = nodes.get(str(node)) if isinstance(nodes, dict) else None
    if not isinstance(row, dict):
        blockers.append("%s literal row %s missing" % (label, node))
        return {}
    return row


def _od2_exact_row(nodes, node, expected, label, blockers):
    row = _od2_row(nodes, node, label, blockers)
    for key, value in expected.items():
        if row.get(key) != value:
            blockers.append("%s row %s field %s changed" % (label, node, key))
    return row


def _od2_clause(text, clause, label, blockers):
    if not isinstance(text, str) or clause not in text:
        blockers.append("source text clause missing: %s" % label)


def _od2_creator_bindings(spec, blockers):
    records = spec.get("creator_source_records")
    if not isinstance(records, dict):
        blockers.append("input.creator_source_records missing")
        return {}, {}
    meta = {}
    selected = {}
    for role, (basename, expected_sha) in OD2_CREATOR_RECORDS.items():
        record = records.get(basename)
        if not isinstance(record, dict):
            blockers.append("creator basename missing for %s: %s" % (role, basename))
            continue
        text = record.get("source_text")
        if not isinstance(text, str) or record.get("path") is None:
            blockers.append("creator source text/path missing for %s" % role)
            continue
        actual_sha = _od2_sha_text(text)
        if actual_sha != expected_sha or record.get("sha256") != expected_sha:
            blockers.append("creator SHA changed for %s" % role)
        try:
            tree = _od2_ast.parse(text, filename=str(record.get("path")), mode="exec")
            ast_sha = _od2_ast_sha(tree)
        except (SyntaxError, ValueError, TypeError) as error:
            blockers.append("creator AST parse failed for %s: %s" % (role, error))
            continue
        clauses = {
            "coefficient_pairs": (
                "def coefficient_pairs",
                "g.unary('exprel',arg)",
                "g.unary('exp',arg)",
                "E.value,A.value",
            ),
            "C2_density": (
                "def coefficient_C2",
                "g.unary('exprel',arg)",
                "N_scaled_density=True",
                "original dlogR: Jacobian exactly once",
            ),
            "C3_source": (
                "def coefficient_C3",
                "CurrentC3SourceTargets.source_packet",
                "ordinary_slow_Z_derivative_order=3",
                "Taylor_coefficient_factorial=6",
            ),
            "loop_source": (
                "def finite_N_functions",
                "g.unary('expm1',x.value)",
                "r=a.mul(roots['E'],exm1);s=a.scale(loop['B_over_Pstar'],invN)",
                "all_five_signed_finite_N_density_cross_terms_and_expm1_retained=True",
            ),
            "outer": (
                "O3_transition",
                "physical_Jacobian_applied_once=True",
                "exact_original_mu_and_exponential_O2_log_radius_Jacobian=True",
            ),
            "control_connection": (
                "O3_transition=Rm+17+sy.exp(40)+x",
                "measure='original dlogR: Jacobian exactly once'",
                "N_scaled_density=True",
            ),
            "limitcontrols": (
                "signed=target.CurrentC3SourceTargets(target_field=self.lower.ranges.target)",
                "same_original_C2_limit_branch=True",
            ),
            "repair": (
                "field.control.ranges.target.functions['actual_terminal_amplitude']",
                "field.control.ranges.target.functions['terminal_N_scaled_histories']",
                "ordinary_slow_Z_derivative_order=3",
            ),
        }
        for clause in clauses[role]:
            _od2_clause(text, clause, role, blockers)
        meta[basename] = {
            "path": record.get("path"),
            "source_sha256": actual_sha,
            "ast_sha256": ast_sha,
            "bytes": len(text.encode("utf-8")),
            "source_text_not_executed": True,
        }
        selected[role] = meta[basename]
    return meta, selected


def _od2_slices(spec, blockers):
    slices = spec.get("literal_slices")
    if not isinstance(slices, dict):
        blockers.append("input.literal_slices missing")
        return {}
    out = {}
    for name, archive_sha in (
        ("native", OD2_NATIVE_ARCHIVE_SHA256),
        ("C3", OD2_C3_ARCHIVE_SHA256),
    ):
        item = slices.get(name)
        if not isinstance(item, dict):
            blockers.append("literal slice %s missing" % name)
            out[name] = {}
            continue
        if item.get("archive_sha256") != archive_sha:
            blockers.append("literal slice %s archive SHA changed" % name)
        if not isinstance(item.get("nodes"), dict):
            blockers.append("literal slice %s nodes missing" % name)
        if not isinstance(item.get("chart_ports"), dict):
            blockers.append("literal slice %s chart_ports missing" % name)
        out[name] = item
    return out


def _od2_validate_geometry(native, c3, blockers, validate_native=True, validate_c3=True):
    nt = native.get("native_O3_transition", {})
    if validate_native:
        if nt.get("chart") != "O3_transition":
            blockers.append("native O3 chart changed")
        if nt.get("domain", {}).get("semantic") != ["0", "1"]:
            blockers.append("native x domain changed")
        geom = nt.get("log_radius_geometry", {})
        if geom.get("offset_literal") != "17 + exp(40) + x":
            blockers.append("native exact transition offset changed")
        if geom.get("physical_measure_literal") != "1" or geom.get("jacobian_literal") != "1":
            blockers.append("native physical Jacobian literal changed")
        if geom.get("physical_dlogR_jacobian_applied_once") is not True:
            blockers.append("native physical Jacobian-once marker missing")
        branch = nt.get("source_branch_variable", {})
        if branch.get("name") != "sigma" or branch.get("J_graph_node") != 241546:
            blockers.append("native sigma branch descriptor changed")
        if branch.get("J_formula") != "J_sigma(x)=integral_0^x sigma(t) dt":
            blockers.append("native J_sigma clause changed")
        for key, digest in OD2_FAMILY_HASHES.items():
            if native.get("inputs", {}).get("archive_source_family", {}).get(key) != digest:
                blockers.append("native source-family hash changed: " + key)
    cj = c3.get("rebuilt_O3_transition_source_jets", {})
    if validate_c3:
        for key, digest in OD2_FAMILY_HASHES.items():
            if cj.get("source_family_literals", {}).get(key) != digest:
                blockers.append("C3 source-family hash changed: " + key)
        if c3.get("rebuilt_archive", {}).get("raw_sha256_before_gzip_decode") != OD2_C3_ARCHIVE_SHA256:
            blockers.append("C3 rebuilt archive SHA changed")
    return nt, cj


def _od2_native_binding(native, slices, blockers):
    nt, _ = _od2_validate_geometry(native, {}, blockers, validate_c3=False)
    finite = nt.get("finite_N", {})
    nrow = finite.get("N", {})
    if nrow.get("node") != OD2_N_NODE or nrow.get("J") != OD2_J:
        blockers.append("native N descriptor changed")
    if nrow.get("definition") != "N=2^(2^J)" or nrow.get("Z_independent") is not True:
        blockers.append("native N positivity/Z-independence changed")
    inv = finite.get("invN", {})
    if inv.get("node") != OD2_NATIVE_INVN_NODE or inv.get("numerator") != 1 or inv.get("denominator") != 18:
        blockers.append("native invN descriptor changed")
    roots = nt.get("root_rows", {}).get("ordinary_Z0_Z3", {})
    loop_rows = nt.get("loop_rows", {})
    loops = loop_rows.get("ordinary_Z0_Z3", {})
    fn = finite.get("ordinary_Z0_Z3", {})
    jets = {
        "E": roots.get("E"),
        "V": roots.get("V"),
        "A": loops.get("A"),
        "B": loops.get("B_over_Pstar"),
        "F_raw_r": fn.get("delta_E_r"),
    }
    for key, expected in OD2_JETS_NATIVE.items():
        if jets.get(key) != list(expected):
            blockers.append("native actual %s jet array changed" % key)
    background = nt.get("background_recipe", {})
    expected_background = {
        "E": "E_O2_buffer_endpoint * exp(-x/2 - mu*J_sigma(x))",
        "V": "0",
        "E_y": "-(1/2 + mu*sigma(x))*E",
        "V_y": "0",
    }
    for key, expected in expected_background.items():
        if background.get(key) != expected:
            blockers.append("native background %s recipe changed" % key)
    expected_background_jets = {
        "E": list(OD2_JETS_NATIVE["E"]),
        "E_y": [241563, 241564, 241565, 241566],
        "V": [0, 0, 0, 0],
        "V_y": [0, 0, 0, 0],
    }
    if background.get("ordinary_Z0_Z3_node_ids") != expected_background_jets:
        blockers.append("native E/E_y/V/V_y ordinary source rows changed")
    if roots.get("Delta") != [242546, 242480, 242486, 242493]:
        blockers.append("native Delta ordinary source rows changed")
    lazy_rows = loop_rows.get("lazy_rows", {})
    for name in ("A", "B_over_Pstar"):
        lazy = lazy_rows.get(name, {})
        expected_nodes = loop_rows.get("full_saved_C5_jet_rows", {}).get(name)
        for key, expected in {
            "node_ids": expected_nodes,
            "flat_value": 0,
            "Delta": 242546,
            "eta": 21,
            "flat_predicate": 242552,
            "active_body_not_evaluated": True,
            "ordinary_orders": [0, 1, 2, 3, 4, 5],
        }.items():
            if lazy.get(key) != expected:
                blockers.append("native lazy %s field %s changed" % (name, key))
    source_formula = finite.get("source_formula_binding", {})
    expected_formula = {
        "r": "E*expm1(A/N)",
        "s": "B_over_Pstar/N",
        "m": "s",
        "h": "r",
        "k": "V*r+E*s+r*s",
        "e": "2*V*s+s^2-(E*r+r^2/2)",
        "p": "E*r+r^2/2",
        "no_numeric_values": True,
    }
    if any(source_formula.get(key) != value for key, value in expected_formula.items()):
        blockers.append("native finite-N source formula binding changed")
    cached_scaling = finite.get("cached_parent_N_scaling_binding", {})
    if cached_scaling != {
        "native_N_scaling_applied_once_bound": True,
        "h_driver_label": "N*r",
        "underlying_h_row": "r",
        "reused_as_data_only": True,
        "new_scaling_proof": False,
    }:
        blockers.append("native cached N-scaling binding changed")
    expected_density = {
        "m": [244378, 244379, 244380, 244381],
        "h": [244352, 244355, 244359, 244364],
        "k": [244436, 244437, 244438, 244439],
        "e": [244538, 244539, 244540, 244541],
        "p": [244526, 244527, 244528, 244529],
    }
    ports = slices.get("native", {}).get("chart_ports", {})
    if ports.get("jets") != jets:
        blockers.append("native chart_ports jets do not match packet arrays")
    if ports.get("coordinate") != 241536 or ports.get("Z_variable") != 20:
        blockers.append("native chart_ports coordinate changed")
    if ports.get("domain") != ["0", "1"] or ports.get("physical_measure_node") != 1:
        blockers.append("native chart_ports geometry changed")
    if ports.get("density") != expected_density:
        blockers.append("native chart_ports density rows changed")
    if ports.get("normalization") != "raw correction density; N factor applied once by definition":
        blockers.append("native chart_ports normalization changed")
    for field in OD2_FIELDS:
        values = fn.get("five_signed_density_increments", {}).get(field)
        if values != expected_density[field]:
            blockers.append("native density rows missing: " + field)
    expected_port_lengths = {
        "background_incoming_from_O2_buffer_endpoint": (OD2_FIELDS, 7),
        "native_chart_incoming_correction": (OD2_FIELDS, 6),
        "background_local_source_history_contributions": (OD2_FIELDS, 7),
        "background_outlet_history": (OD2_FIELDS, 7),
        "finite_N_local_source_history_contributions": (OD2_FIELDS, 6),
        "finite_N_outlet_history": (OD2_FIELDS, 6),
        "native_O3_power_endpoint_correction": (OD2_FIELDS, 6),
        "pressure_incoming_integrals": (("S0", "S1"), 4),
        "pressure_local_source_integrals": (("S0", "S1"), 4),
        "pressure_cumulative_integrals": (("S0", "S1"), 4),
        "pressure_ports": (("m13", "m15"), 4),
    }
    native_ports = nt.get("ports", {})
    for port_name, (port_fields, expected_length) in expected_port_lengths.items():
        port = native_ports.get(port_name, {})
        for field in port_fields:
            if not isinstance(port.get(field), list) or len(port[field]) != expected_length:
                blockers.append("native port %s.%s changed" % (port_name, field))
    for name, expected in {
        "background_representative": (241633, False),
        "finite_N_representative": (244575, True),
        "pressure_S0_representative": (249951, None),
        "pressure_S1_representative": (249968, None),
    }.items():
        row = nt.get("integral_port_metadata", {}).get(name, {})
        if row.get("node") != expected[0] or row.get("physical_Jacobian") != 1 or row.get("physical_Jacobian_applied_once") is not True:
            blockers.append("native integral metadata %s changed" % name)
        if row.get("numerical_value_not_installed") is not True:
            blockers.append("native integral metadata %s is numeric" % name)
        if expected[1] is not None and row.get("finite_N_correction") is not expected[1]:
            blockers.append("native integral metadata %s finite-N marker changed" % name)
    shared = nt.get("shared_mu_eta_descriptors", {})
    mu = shared.get("mu", {})
    if {
        "logmu_node": mu.get("logmu_node"),
        "mu_node": mu.get("mu_node"),
        "Tw_node": mu.get("Tw_node"),
        "recipe": mu.get("recipe"),
    } != {
        "logmu_node": 16,
        "mu_node": 17,
        "Tw_node": 205147,
        "recipe": "mu=1/1000*exp(-4*(exp(40)+11))",
    }:
        blockers.append("native shared mu descriptor changed")
    eta = shared.get("eta", {})
    if {
        "node": eta.get("node"),
        "name": eta.get("name"),
        "source_attribute": eta.get("source_attribute"),
        "quantity_is_exact_original_parameter_not_range_endpoint": eta.get("quantity_is_exact_original_parameter_not_range_endpoint"),
        "Z_independent": eta.get("Z_independent"),
    } != {
        "node": 21,
        "name": "eta",
        "source_attribute": "eta_log",
        "quantity_is_exact_original_parameter_not_range_endpoint": True,
        "Z_independent": True,
    }:
        blockers.append("native shared eta descriptor changed")
    lazy_parameter = shared.get("lazy_branch_parameter", {})
    if {
        "eta_node": lazy_parameter.get("eta_node"),
        "flat_value": lazy_parameter.get("flat_value"),
        "predicate_node": lazy_parameter.get("predicate_node"),
        "predicate_literal": lazy_parameter.get("predicate_literal"),
        "active_body_not_evaluated_on_flat_branch": lazy_parameter.get("active_body_not_evaluated_on_flat_branch"),
    } != {
        "eta_node": 21,
        "flat_value": 0,
        "predicate_node": 242552,
        "predicate_literal": "Delta >= eta",
        "active_body_not_evaluated_on_flat_branch": True,
    }:
        blockers.append("native shared lazy-branch descriptor changed")
    nodes = slices.get("native", {}).get("nodes", {})
    _od2_exact_row(nodes, 244325, {"operation": "product", "arguments": [244098, 99411]}, "native", blockers)
    _od2_exact_row(nodes, 244331, {"operation": "analytic_unary", "name": "exp", "argument": 244325}, "native", blockers)
    _od2_exact_row(nodes, 244351, {"operation": "analytic_unary", "name": "expm1", "argument": 244325}, "native", blockers)
    _od2_exact_row(nodes, 244352, {"operation": "product", "arguments": [241552, 244351]}, "native", blockers)
    _od2_exact_row(nodes, 244378, {"operation": "product", "arguments": [244136, 99411]}, "native", blockers)
    _od2_exact_row(nodes, 99411, {"operation": "positive_quotient", "numerator": 1, "denominator": 18}, "native", blockers)
    _od2_row(nodes, 244098, "native A0", blockers)
    for row in (244098, 244099, 244100, 244101, 244136, 244137, 244138, 244139):
        _od2_row(nodes, row, "native higher jet", blockers)
    exp_rows = finite.get("exp_expm1_jet_constructor_rows", {})
    if exp_rows.get("exp_A_over_N_rows") != [244331, 244332, 244335, 244339, 244344, 244350] or exp_rows.get("expm1_A_over_N_rows") != [244351, 244332, 244335, 244339, 244344, 244350]:
        blockers.append("native exp/expm1 higher-jet constructor rows changed")
    if exp_rows.get("constructor_literal") != "exm1=FunctionJet((g.unary('expm1',x.value),)+rows(ex)[1:])" or exp_rows.get("higher_rows_reuse_full_exp_recurrence") is not True:
        blockers.append("native exp/expm1 constructor semantics changed")
    return {
        "namespace": "native_true_outer_C6_C5",
        "archive_sha256": OD2_NATIVE_ARCHIVE_SHA256,
        "chart": nt,
        "jets": jets,
        "raw_density": fn.get("five_signed_density_increments"),
        "ports": nt.get("ports"),
        "chart_ports": ports,
        "finite_N": finite,
        "higher_jet_semantics": {
            "constructor_source": OD2_CREATOR_RECORDS["loop_source"][0],
            "constructor_source_sha256": OD2_CREATOR_RECORDS["loop_source"][1],
            "expm1_rows": finite.get("exp_expm1_jet_constructor_rows"),
            "source_ports_are_actual": True,
            "no_old_helper_calls": True,
        },
    }


def _od2_c3_binding(c3, producer, scalar, slices, blockers):
    cj = c3.get("rebuilt_O3_transition_source_jets", {})
    _, _ = _od2_validate_geometry({}, c3, blockers, validate_native=False)
    if producer.get("schema") != "transition_C3_source_literals_v1":
        blockers.append("C3 producer source packet schema changed")
    producer_sources = {
        row.get("file"): row.get("sha256")
        for row in producer.get("authenticated_source_files", [])
        if isinstance(row, dict)
    }
    for role in ("C3_source", "repair", "C2_density", "coefficient_pairs", "control_connection"):
        basename, expected_sha = OD2_CREATOR_RECORDS[role]
        if producer_sources.get(basename) != expected_sha:
            blockers.append("C3 producer authenticated source changed: %s" % basename)
    if cj.get("chart") != "O3_transition" or cj.get("ordinary_Z_orders") != list(OD2_ORDERS):
        blockers.append("C3 O3 transition chart/order descriptor changed")
    jets = {key: cj.get(key if key != "F" else "F_N") for key in ("E", "V", "A", "B", "F")}
    for key, expected in OD2_JETS_C3.items():
        actual_key = "F" if key == "F_N" else key
        if jets.get(actual_key) != list(expected):
            blockers.append("C3 actual %s jet array changed" % key)
    ports = slices.get("C3", {}).get("chart_ports", {})
    expected_port_jets = {key: jets["F"] if key == "F_N" else jets[key] for key in ("E", "V", "A", "B", "F_N")}
    if ports.get("jets") != expected_port_jets:
        blockers.append("C3 chart_ports jets do not match rebuilt arrays")
    if ports.get("coordinate") != 4582 or ports.get("Z_variable") != 20:
        blockers.append("C3 chart_ports coordinate changed")
    if ports.get("domain") != ["0", "1"] or ports.get("physical_measure_node") != 1:
        blockers.append("C3 chart_ports geometry changed")
    nodes = slices.get("C3", {}).get("nodes", {})
    for jet_name, jet_rows in jets.items():
        for node in jet_rows if isinstance(jet_rows, list) else ():
            _od2_row(nodes, node, "C3 %s jet" % jet_name, blockers)
    scalar_rows = scalar.get("selected_scalar_nodes", {}).get("literal_rows", {})
    if scalar.get("rebuilt_archive", {}).get("raw_sha256_before_gzip_decode") != OD2_C3_ARCHIVE_SHA256:
        blockers.append("scalar addendum archive SHA changed")
    for node, expected in (
        (0, {"operation": "exact_rational", "numerator": 0, "denominator": 1}),
        (1, {"operation": "exact_rational", "numerator": 1, "denominator": 1}),
        (48, {"operation": "exact_rational", "numerator": 2, "denominator": 1}),
        (13766, {"operation": "exact_rational", "numerator": 3, "denominator": 1}),
        (10287, {"operation": "positive_quotient", "numerator": 1, "denominator": 18}),
    ):
        row = scalar_rows.get(str(node))
        if not isinstance(row, dict) or any(row.get(k) != v for k, v in expected.items()):
            blockers.append("C3 scalar row %s changed" % node)
    for node in (13766, 10287):
        scalar_row = scalar_rows.get(str(node), {})
        _od2_exact_row(
            nodes,
            node,
            {
                "operation": scalar_row.get("operation"),
                "numerator": scalar_row.get("numerator"),
                "denominator": scalar_row.get("denominator"),
            },
            "C3 scalar",
            blockers,
        )
    frows = c3.get("rebuilt_F_N_immediate_literal_rows", {})
    expected_rows = {
        4735: {"operation": "product", "arguments": [4590, 4718, 4734]},
        4739: {"operation": "sum", "arguments": [4736, 4738]},
        11937: {"operation": "sum", "arguments": [11932, 11933, 11936]},
        24061: {"operation": "sum", "arguments": [24054, 24055, 24056, 24060]},
    }
    for node, expected in expected_rows.items():
        _od2_exact_row(nodes, node, expected, "C3 F producer", blockers)
    _od2_exact_row(nodes, 4734, {"operation": "analytic_unary", "name": "exprel", "argument": 4733}, "C3 F0", blockers)
    _od2_exact_row(nodes, 4733, {"operation": "positive_quotient", "numerator": 4718, "denominator": 18}, "C3 F0", blockers)
    _od2_exact_row(nodes, 4737, {"operation": "analytic_unary", "name": "exp", "argument": 4733}, "C3 F1", blockers)
    _od2_exact_row(nodes, 11930, {"operation": "analytic_unary", "name": "exprel", "argument": 11929}, "C3 F2", blockers)
    _od2_exact_row(nodes, 11931, {"operation": "analytic_unary", "name": "exp", "argument": 11929}, "C3 F2", blockers)
    _od2_exact_row(nodes, 11929, {"operation": "product", "arguments": [4718, 10287]}, "C3 F2", blockers)
    _od2_exact_row(nodes, 24059, {"operation": "sum", "arguments": [22254, 24057, 24058]}, "C3 F3", blockers)
    for node, expected in {
        4736: {"operation": "product", "arguments": [4591, 4718, 4734]},
        4738: {"operation": "product", "arguments": [4590, 4725, 4737]},
        11932: {"operation": "product", "arguments": [9988, 4718, 11930]},
        11933: {"operation": "product", "arguments": [48, 4591, 4725, 11931]},
        11934: {"operation": "product", "arguments": [4725, 4725, 10287]},
        11935: {"operation": "sum", "arguments": [10273, 11934]},
        11936: {"operation": "product", "arguments": [4590, 11931, 11935]},
        24054: {"operation": "product", "arguments": [21725, 4718, 11930]},
        24055: {"operation": "product", "arguments": [13766, 9988, 4725, 11931]},
        24056: {"operation": "product", "arguments": [13766, 4591, 11931, 11935]},
        24057: {"operation": "product", "arguments": [13766, 4725, 10273, 10287]},
        24058: {"operation": "product", "arguments": [4725, 4725, 4725, 10287, 10287]},
        24060: {"operation": "product", "arguments": [4590, 11931, 24059]},
    }.items():
        _od2_exact_row(nodes, node, expected, "C3 F producer support", blockers)
    _od2_exact_row(nodes, 48, {"operation": "exact_rational", "numerator": 2, "denominator": 1}, "C3 F producer scalar", blockers)
    _od2_row(nodes, 22254, "C3 A3 support", blockers)
    if frows.get("F3", {}).get("A_ZZZ_support_row") != 22254:
        blockers.append("C3 F3 A_ZZZ support row changed")
    if c3.get("rebuilt_leading_function_recipe_metadata_rows", {}).get("rows", {}).get("21725", {}).get("Taylor_coefficient_factorial") != 6:
        blockers.append("C3 E3 factor-6 metadata changed")
    e3 = scalar.get("actual_x_and_leading_recipe_rows", {}).get("E3_full_literal_row", {})
    if {
        "operation": e3.get("operation"),
        "native_chart": e3.get("native_chart"),
        "quantity": e3.get("quantity"),
        "Z_order": e3.get("Z_order"),
        "coordinate": e3.get("coordinate"),
        "Z_variable": e3.get("Z_variable"),
        "ordinary_slow_Z_derivative_order": e3.get("ordinary_slow_Z_derivative_order"),
        "Taylor_coefficient_factorial": e3.get("Taylor_coefficient_factorial"),
        "phase_independent_leading_source": e3.get("phase_independent_leading_source"),
        "defining_quantity_not_a_range_value": e3.get("defining_quantity_not_a_range_value"),
        "derivative_of_range_endpoint": e3.get("derivative_of_range_endpoint"),
    } != {
        "operation": "current_original_leading_function_recipe",
        "native_chart": "O3_transition",
        "quantity": "E",
        "Z_order": 3,
        "coordinate": 4582,
        "Z_variable": 20,
        "ordinary_slow_Z_derivative_order": 3,
        "Taylor_coefficient_factorial": 6,
        "phase_independent_leading_source": True,
        "defining_quantity_not_a_range_value": True,
        "derivative_of_range_endpoint": False,
    }:
        blockers.append("C3 E3 full source literal row changed")
    if e3.get("source_family") != OD2_FAMILY_HASHES:
        blockers.append("C3 E3 source-family binding changed")
    e3_recipe = e3.get("recipe", {})
    if {
        "module": e3_recipe.get("module"),
        "module_sha256": e3_recipe.get("module_sha256"),
        "function": e3_recipe.get("function"),
        "native_chart": e3_recipe.get("native_chart"),
    } != {
        "module": "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py",
        "module_sha256": "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7",
        "function": "WholeZAllNOuterRcFunctions.leading_packet",
        "native_chart": "O3_transition",
    }:
        blockers.append("C3 E3 leading recipe source binding changed")
    e3_projection = e3.get("source_projection_binding", {})
    if {
        "module": e3_projection.get("module"),
        "module_sha256": e3_projection.get("module_sha256"),
        "function": e3_projection.get("function"),
    } != {
        "module": "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
        "module_sha256": OD2_CREATOR_RECORDS["C3_source"][1],
        "function": "CurrentC3SourceTargets.source_packet",
    }:
        blockers.append("C3 E3 projection source binding changed")
    if frows.get("F3", {}).get("row") != 24061:
        blockers.append("C3 F3 row 24061 not bound")
    selected_ab = c3.get("rebuilt_A_B_primitive_rows", {}).get("selected_literal_rows", {})
    for node in (4718, 4725, 10273, 22254, 4727, 4732, 10279, 22255):
        selected = selected_ab.get(str(node))
        literal = _od2_row(nodes, node, "C3 A/B primitive", blockers)
        if not isinstance(selected, dict) or selected.get("operation") != "original_lazy_flat_branch":
            blockers.append("C3 A/B primitive source row %s changed" % node)
        if literal.get("operation") != "original_lazy_flat_branch":
            blockers.append("C3 A/B literal row %s changed" % node)
    if selected_ab.get("22254", {}).get("ordinary_slow_Z_derivative_order") != 3:
        blockers.append("C3 A3 ordinary order metadata changed")
    if selected_ab.get("22255", {}).get("ordinary_slow_Z_derivative_order") != 3:
        blockers.append("C3 B3 ordinary order metadata changed")
    routes = c3.get("rebuilt_direct_upstream_integral_routes", {})
    density = {}
    groups = routes.get("groups", {})
    for field in OD2_FIELDS:
        density[field] = []
        for order in OD2_ORDERS:
            entry = groups.get(field, {}).get("Z" + str(order), {})
            if not isinstance(entry.get("density_ref"), int):
                blockers.append("C3 density route missing %s Z%s" % (field, order))
            density[field].append(entry.get("density_ref"))
    return {
        "namespace": "rebuilt_C3_repair_band_archive",
        "archive_sha256": OD2_C3_ARCHIVE_SHA256,
        "jets": jets,
        "density_refs": density,
        "routes": routes,
        "ports": ports,
        "F_literal_rows": frows,
        "source_creator_chain": c3.get("source_creator_chain"),
        "producer_source_hashes": producer_sources,
        "literal_rows": c3.get("rebuilt_A_B_primitive_rows"),
        "excluded_primary_F3": 22419,
        "excluded_quiet_inlet_m3": 24145,
        "repair_local_bump_F_G_used": False,
    }


def _od2_integrals(c3, scalar, slices, blockers):
    routes = c3.get("rebuilt_direct_upstream_integral_routes", {})
    groups = routes.get("groups", {})
    full = scalar.get("all_80_integral_rows", {}).get("full_literal_bodies", {})
    body_nodes = slices.get("C3", {}).get("nodes", {})
    cells = {
        0: ([0, 1], [1, 4], 0, 247),
        1: ([1, 4], [1, 2], 247, 28),
        2: ([1, 2], [3, 4], 28, 49),
        3: ([3, 4], [1, 1], 49, 1),
    }
    rows = []
    for field in OD2_FIELDS:
        for order in OD2_ORDERS:
            route = groups.get(field, {}).get("Z" + str(order), {})
            integral_ids = route.get("integral_rows", [])
            if len(integral_ids) != 4:
                blockers.append("C3 %s Z%s does not have four integral rows" % (field, order))
            for cell_index, node in enumerate(integral_ids):
                expected_left, expected_right, lower, upper = cells.get(cell_index, (None, None, None, None))
                if expected_left is None:
                    blockers.append("C3 integral partition cell index changed %s Z%s" % (field, order))
                    continue
                body = full.get(str(node))
                if not isinstance(body, dict):
                    blockers.append("scalar addendum integral body missing %s" % node)
                    continue
                literal = _od2_row(body_nodes, node, "C3 integral", blockers)
                common_fields = {
                    "operation": "definite_integral",
                    "variable": "native_O3_transition",
                    "actual_current_chart": "O3_transition",
                    "measure": "original dlogR: Jacobian exactly once",
                    "lower": lower,
                    "upper": upper,
                    "exact_native_left": expected_left,
                    "exact_native_right": expected_right,
                    "N_scaled_density": True,
                    "exact_function_integral": True,
                    "numerical_value_not_installed": True,
                }
                for target_name, target in (("body", body), ("literal", literal)):
                    for key, value in common_fields.items():
                        if target.get(key) != value:
                            blockers.append("C3 %s integral %s field %s changed" % (target_name, node, key))
                if order in (0, 1):
                    for target_name, target in (("body", body), ("literal", literal)):
                        if target.get("Z_derivative_order_held_phase") is not True:
                            blockers.append("C3 %s integral %s missing held-phase flag" % (target_name, node))
                        if "ordinary_slow_Z_derivative_order" in target:
                            blockers.append("C3 %s integral %s has an unexpected ordinary-order field" % (target_name, node))
                else:
                    for target_name, target in (("body", body), ("literal", literal)):
                        if target.get("ordinary_slow_Z_derivative_order") != order:
                            blockers.append("C3 %s integral %s ordinary order changed" % (target_name, node))
                        if target.get("integration_endpoints_kernel_radius_and_global_phase_Z_independent") is not True:
                            blockers.append("C3 %s integral %s endpoint-independence flag changed" % (target_name, node))
                rows.append({"field": field, "ordinary_Z_order": order, "partition_cell_index": cell_index, "node": node, "body": body})
            product_id = route.get("integrand")
            product = routes.get("direct_product_literal_rows", {}).get(str(product_id), {})
            expected_args = route.get("product_arguments")
            if product.get("operation") != "product" or product.get("arguments") != expected_args:
                blockers.append("C3 route product changed %s Z%s" % (field, order))
            literal_product = _od2_row(body_nodes, product_id, "C3 route product", blockers)
            if literal_product.get("operation") != "product" or literal_product.get("arguments") != expected_args:
                blockers.append("C3 literal route product changed %s Z%s" % (field, order))
            if not isinstance(expected_args, list) or len(expected_args) != 2:
                blockers.append("C3 route product arguments malformed %s Z%s" % (field, order))
            elif route.get("kernel_ref") != expected_args[0] or route.get("density_ref") != expected_args[1]:
                blockers.append("C3 density_ref/product argument mismatch %s Z%s" % (field, order))
    if len(rows) != 80:
        blockers.append("C3 integral binding count is not 80")
    return {
        "count": len(rows),
        "rows": rows,
        "cells": cells,
        "physical_measure": "original dlogR: Jacobian exactly once",
        "source_values_evaluated": False,
    }


def _od2_cached(spec, blockers):
    packet = spec.get("cached_signed_result")
    if not isinstance(packet, dict) or packet.get("sha256") != OD2_CACHED_RESULT_SHA256:
        blockers.append("cached signed result packet/SHA changed")
        return {}, {}, {}
    payload = packet.get("payload")
    if not isinstance(payload, dict):
        blockers.append("cached signed result payload missing")
        return packet, {}, {}
    result = payload.get("signed_density_algebra_result", {})
    defects = result.get("factorized_signed_defects", {})
    rows = defects.get("rows")
    formulas = defects.get("formulas")
    if not isinstance(rows, list) or len(rows) != 20 or not isinstance(formulas, dict):
        blockers.append("cached factorized signed defects rows/formulas malformed")
        return packet, {}, {}
    if set(formulas) != {"e2", "h", "k", "m", "p2"}:
        blockers.append("cached factorized formula channel set changed")
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or row.get("residual_zero") is not True:
            blockers.append("cached residual_zero data row changed")
            continue
        seen.add((row.get("channel"), row.get("ordinary_Z_order")))
        if row.get("channel") not in ("m", "h", "k", "e2", "p2") or row.get("ordinary_Z_order") not in OD2_ORDERS:
            blockers.append("cached factorized row channel/order changed")
    if len(seen) != 20:
        blockers.append("cached factorized signed defects do not cover 20 channel/order rows")
    return packet, defects, result


def _od2_atom(name, namespace, archive_sha, node, order, **extra):
    value = {
        "atom": name,
        "namespace": namespace,
        "archive_sha256": archive_sha,
        "node": node,
        "ordinary_Z_order": order,
    }
    value.update(extra)
    return value


def _od2_atomrefs(native, c3):
    n = native["jets"]
    c = c3["jets"]
    refs = {}
    for order in OD2_ORDERS:
        refs["En_Z%s" % order] = _od2_atom("En_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, n["E"][order], order)
        refs["Ec_Z%s" % order] = _od2_atom("Ec_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c["E"][order], order)
        refs["Vn_Z%s" % order] = _od2_atom("Vn_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, n["V"][order], order)
        refs["Vc_Z%s" % order] = _od2_atom("Vc_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c["V"][order], order)
        refs["Bn_Z%s" % order] = _od2_atom("Bn_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, n["B"][order], order, semantic="B_over_Pstar", formal_quantity="normalized_loop_B", physical_B_not_identified=True)
        refs["Bc_Z%s" % order] = _od2_atom("Bc_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c["B"][order], order, semantic="primitive B", formal_quantity="normalized_loop_B", physical_B_not_identified=True)
        refs["Fc_Z%s" % order] = _od2_atom("Fc_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c["F"][order], order, semantic="exact C3 F_N exprel producer", N_scaled_source=True, no_extra_N_scale=True)
        refs["Fn_Z%s" % order] = _od2_atom(
            "Fn_Z%s" % order,
            "native_true_outer_C6_C5",
            OD2_NATIVE_ARCHIVE_SHA256,
            n["F_raw_r"][order],
            order,
            operation="scale_by_common_N",
            factor={"node": OD2_N_NODE, "descriptor": "N=2^(2^J)", "J": OD2_J, "ordinary_Z_independent": True},
            operand={
                "atom": "native_raw_delta_E_r_Z%s" % order,
                "namespace": "native_true_outer_C6_C5",
                "archive_sha256": OD2_NATIVE_ARCHIVE_SHA256,
                "node": n["F_raw_r"][order],
                "ordinary_Z_order": order,
            },
            scaled_once=True,
        )
    refs["eps"] = {
        "atom": "eps",
        "operation": "shared_inverse_N",
        "descriptor": "eps=1/N",
        "J": OD2_J,
        "ordinary_Z_independent": True,
        "archive_sha256": [OD2_NATIVE_ARCHIVE_SHA256, OD2_C3_ARCHIVE_SHA256],
        "namespace": "shared_native_C3_parameter",
        "node": OD2_N_NODE,
        "ordinary_Z_order": None,
        "native_invN_ref": _od2_atom("native_invN", "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, OD2_NATIVE_INVN_NODE, None, numerator=1, denominator=18),
        "C3_invN_ref": _od2_atom("C3_invN", "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, OD2_C3_INVN_NODE, None, numerator=1, denominator=18),
    }
    return refs


def _od2_instantiate_terms(rows, refs, blockers):
    out = []
    for row in rows:
        channel = row.get("channel")
        physical_channel = {"e2": "e", "p2": "p"}.get(channel, channel)
        coefficient_denominator = 2 if channel in ("e2", "p2") else 1
        def terms(key):
            result = []
            term_rows = row.get(key, [])
            if not isinstance(term_rows, list):
                blockers.append("cached factorized %s terms are not a literal list" % key)
                return result
            for term in term_rows:
                if not isinstance(term, dict):
                    blockers.append("cached factorized %s term is not a literal mapping" % key)
                    continue
                exponents = []
                exponent_rows = term.get("exponents", [])
                if not isinstance(exponent_rows, list):
                    blockers.append("cached factorized %s exponents are not a literal list" % key)
                    exponent_rows = []
                for exponent in exponent_rows:
                    if not isinstance(exponent, dict):
                        blockers.append("cached factorized %s exponent is not a literal mapping" % key)
                        continue
                    atom_name = exponent.get("atom")
                    if atom_name not in refs:
                        blockers.append("cached factorized atom is not source-bound: %s" % atom_name)
                    exponents.append({"power": exponent.get("power"), "atom": refs.get(atom_name)})
                coefficient = term.get("coefficient")
                result.append({
                    "coefficient": term.get("coefficient"),
                    "physical_channel": physical_channel,
                    "physical_coefficient": {"numerator": coefficient, "denominator": coefficient_denominator},
                    "exponents": exponents,
                })
            return result
        out.append({
            "channel": channel,
            "physical_channel": physical_channel,
            "coefficient_denominator": coefficient_denominator,
            "ordinary_Z_order": row.get("ordinary_Z_order"),
            "residual_zero_data_only": row.get("residual_zero") is True,
            "direct_signed_difference_source_terms": terms("direct_signed_difference"),
            "factorized_gap_source_terms": terms("factorized_gap"),
            "source_values_evaluated": False,
        })
    return out


def _od2_density_gaps(native, c3):
    rows = []
    for field in OD2_FIELDS:
        for order in OD2_ORDERS:
            raw_node = native["raw_density"][field][order]
            c3_node = c3["density_refs"][field][order]
            raw_ref = _od2_atom("native_raw_%s_Z%s" % (field, order), "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, raw_node, order, raw_density=True)
            left = {
                "operation": "scale_by_common_N",
                "factor": {"node": OD2_N_NODE, "descriptor": "N=2^(2^J)", "J": OD2_J, "archive_sha256": OD2_NATIVE_ARCHIVE_SHA256},
                "operand": raw_ref,
                "scaled_once": True,
            }
            right = _od2_atom("C3_Nscaled_%s_Z%s" % (field, order), "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c3_node, order, N_scaled_density=True, no_extra_N_scale=True, route_namespace="rebuilt_direct_upstream_integral_routes")
            rows.append({
                "id": "%s_Z%s" % (field, order),
                "field": field,
                "ordinary_Z_order": order,
                "formula": "N*native_raw_density_%s^%s - actual_C3_Nscaled_density_%s^%s" % (field, order, field, order),
                "left": left,
                "right": right,
                "residual_zero_claim": False,
                "source_values_evaluated": False,
                "new_density_polynomial_check": False,
                "new_differentiation": False,
            })
    return rows


def _od2_f_gap_recipes(native, c3):
    refs = _od2_atomrefs(native, c3)
    rows = []
    for order in OD2_ORDERS:
        rows.append({
            "ordinary_Z_order": order,
            "native_E_ref": _od2_atom("En_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, native["jets"]["E"][order], order),
            "C3_E_ref": _od2_atom("Ec_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c3["jets"]["E"][order], order),
            "native_A_ref": _od2_atom("An_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, native["jets"]["A"][order], order),
            "C3_A_ref": _od2_atom("Ac_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c3["jets"]["A"][order], order),
            "native_F_N_ref": refs["Fn_Z%s" % order],
            "C3_F_N_ref": refs["Fc_Z%s" % order],
            "D_definition": "D=An-Ac",
            "H_definition": "H=En-Ec",
            "entire_exprel": True,
            "functional_recipe": "DeltaF[%s]=D_Z^%s[En*An*exprel(eps*An)-Ec*Ac*exprel(eps*Ac)]" % (order, order),
            "cocycle_recipe": "DeltaF[%s]=D_Z^%s[H*Ac*exprel(eps*Ac)+En*exp(eps*Ac)*D*exprel(eps*D)]" % (order, order),
            "native_F_function_view": {"operation": "whole_function_ordinary_Z_view", "rows": [refs["Fn_Z%s" % j] for j in OD2_ORDERS], "selected_order": order},
            "C3_F_function_view": {"operation": "whole_function_ordinary_Z_view", "rows": [refs["Fc_Z%s" % j] for j in OD2_ORDERS], "selected_order": order},
            "source_values_evaluated": False,
            "gap_equality_checked": False,
        })
    return rows


def _od2_math_contract(native, c3, creators, clauses, atomrefs):
    cocycle_map = {}
    for order in OD2_ORDERS:
        cocycle_map["En%s" % order] = atomrefs["En_Z%s" % order]
        cocycle_map["Ec%s" % order] = atomrefs["Ec_Z%s" % order]
        cocycle_map["An%s" % order] = _od2_atom("An_Z%s" % order, "native_true_outer_C6_C5", OD2_NATIVE_ARCHIVE_SHA256, native["jets"]["A"][order], order)
        cocycle_map["Ac%s" % order] = _od2_atom("Ac_Z%s" % order, "rebuilt_C3_repair_band_archive", OD2_C3_ARCHIVE_SHA256, c3["jets"]["A"][order], order)
    producer_generic_map = {
        "En%s" % order: cocycle_map["Ec%s" % order]
        for order in OD2_ORDERS
    }
    producer_generic_map.update({
        "Ac%s" % order: cocycle_map["Ac%s" % order]
        for order in OD2_ORDERS
    })
    return {
        "actual_source_bindings": {"passed": True, "native": native, "C3": c3, "creators": creators, "source_values_evaluated": False},
        "shared_symbolic_N_node": OD2_N_NODE,
        "shared_symbolic_N_descriptor": "J" + OD2_J,
        "ordinary_Z_independent_N": True,
        "ordinary_derivative_vs_Taylor_correct": True,
        "parent_normalization_data_only": True,
        "finite_N_producer_C3_exprel_native_expm1_bound": True,
        "actual_C3_exprel_producer_clauses": {"passed": True, "clauses": clauses},
        "actual_C3_exprel_producer_inputs": {"passed": True, "native_E": native["jets"]["E"], "native_A": native["jets"]["A"], "C3_E": c3["jets"]["E"], "C3_A": c3["jets"]["A"], "C3_F": c3["jets"]["F"], "producer_generic_source_mapping": producer_generic_map, "producer_mapping_separate_from_cocycle": True},
        "cocycle_source_atom_mapping": cocycle_map,
        "producer_generic_E_to_C3_Ec_mapping": {"En%s" % order: producer_generic_map["En%s" % order] for order in OD2_ORDERS},
        "producer_generic_Ac_to_C3_A_mapping": {"Ac%s" % order: producer_generic_map["Ac%s" % order] for order in OD2_ORDERS},
    }


def _od2_base(source, blockers, **extra):
    result = {
        "schema": OD2_RESULT_SCHEMA,
        "component": OD2_COMPONENT,
        "math_draft_sha256": OD2_MATH_DRAFT_SHA256,
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": False,
        "source_packet_definition_accepted": False,
        "source_packet_accepted_owned": False,
        "source_packet_accepted_inherited": False,
        "input_digest": source.get("input_digest") if isinstance(source, dict) else None,
        "input_raw_sha256": source.get("input_raw_sha256") if isinstance(source, dict) else None,
        "numeric": False,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "old_generic_identity_cached_data_only": True,
        "old_generic_identity_calls": 0,
        "old_operator_calls": 0,
        "old_source_consumer_calls": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "old_numeric_calls": 0,
        "numeric_calls": 0,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
        "blockers": list(blockers),
        "new_flags": {name: False for name in OD2_FLAGS},
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_generic_identity_calls": 0,
            "old_operator_calls": 0,
            "old_source_consumer_calls": 0,
            "old_provider_calls": 0,
            "old_checker_calls": 0,
            "old_normalizer_calls": 0,
            "old_source_value_reducer_calls": 0,
            "old_numeric_calls": 0,
            "numeric_calls": 0,
            "quadratures": 0,
            "temporal_recurrence_calls": 0,
            "new_exponential_cocycle_order_checks": 0,
            "new_exprel_producer_order_checks": 0,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": 0,
            "actual_producer_order_failures": 0,
            "checked_residual_count": 0,
            "checked_residual_failures": 0,
            "new_exponential_cocycle_order_checks4": 0,
            "new_exprel_producer_order_checks4": 0,
            "new_actual_signed_density_gap_rows": 0,
            "new_actual_signed_density_gap_rows20": 0,
            "new_density_proof_rows": 0,
            "source_only": True,
            "metadata_only": False,
        },
        "open_branch_issues": (
            "primitive source gap unknown",
            "P2 unknown",
            "phase unknown",
            "native/C3 source equality and zero-driver claims withheld",
            "no absolute Rm/Rc or common physical point claim",
            "physical B=2*V_y remains distinct from normalized loop B_over_Pstar",
        ),
    }
    result.update(extra)
    for name in OD2_FLAGS:
        result[name] = False
    return result


def derive(source, journal=None):
    """Bind actual source literals, then invoke the new cocycle exactly once."""

    blockers = []
    spec = _od2_spec(source)
    if spec.get("schema") != OD2_INPUT_SCHEMA:
        blockers.append("input schema changed")
    if spec.get("component") != OD2_COMPONENT:
        blockers.append("input component changed")
    if blockers:
        return _od2_base(source, blockers)
    native_packet, native_payload = _od2_packet(spec, "native_packet", OD2_NATIVE_PACKET_SHA256, blockers)
    c3_packet, c3_payload = _od2_packet(spec, "C3_packet", OD2_C3_PACKET_SHA256, blockers)
    producer_packet, producer_payload = _od2_packet(spec, "C3_producer_packet", OD2_C3_PRODUCER_PACKET_SHA256, blockers)
    scalar_packet, scalar_payload = _od2_packet(spec, "C3_scalar_addendum", OD2_C3_SCALAR_SHA256, blockers)
    cached_packet, cached_defects, cached_result = _od2_cached(spec, blockers)
    slices = _od2_slices(spec, blockers)
    creator_meta, creator_selected = _od2_creator_bindings(spec, blockers)
    if blockers:
        return _od2_base(source, blockers, packet_hashes={"native": native_packet.get("sha256"), "C3": c3_packet.get("sha256"), "C3_producer": producer_packet.get("sha256"), "C3_scalar_addendum": scalar_packet.get("sha256"), "cached": cached_packet.get("sha256")}, creator_source_records=creator_meta, literal_slices=slices)
    native = _od2_native_binding(native_payload, slices, blockers)
    c3 = _od2_c3_binding(c3_payload, producer_payload, scalar_payload, slices, blockers)
    integral_data = _od2_integrals(c3_payload, scalar_payload, slices, blockers)
    if blockers:
        return _od2_base(source, blockers, native_source_binding=native, C3_source_binding=c3, integral_binding=integral_data, creator_source_records=creator_meta, literal_slices=slices)
    atomrefs = _od2_atomrefs(native, c3)
    density_gaps = _od2_density_gaps(native, c3)
    f_gaps = _od2_f_gap_recipes(native, c3)
    instantiated = _od2_instantiate_terms(cached_defects["rows"], atomrefs, blockers)
    h4 = next(row for row in instantiated if row["channel"] == "h" and row["ordinary_Z_order"] == 3)
    h4["correlated_F_gap_recipe"] = f_gaps[3]
    if blockers:
        return _od2_base(
            source,
            blockers,
            native_source_binding=native,
            C3_source_binding=c3,
            integral_binding=integral_data,
            atomrefs=atomrefs,
            actual_correlated_F_gap_functional_recipes=f_gaps,
            signed_density_gap_functions=density_gaps,
            cached_factorized_signed_defects={"formulas": cached_defects.get("formulas"), "rows": instantiated},
        )
    math_contract = _od2_math_contract(native, c3, creator_meta, {"C3": c3["F_literal_rows"], "native": native["higher_jet_semantics"]}, atomrefs)
    _od2_emit(journal, "O3_transition_actual_source_bindings_passed", source_only=True)
    try:
        cocycle = derive_exponential_cocycle(math_contract, journal)
    except Exception as error:
        return _od2_base(source, ["fresh cocycle invocation failed: %s" % error], native_source_binding=native, C3_source_binding=c3, math_contract=math_contract)
    algebra = cocycle.get("algebra", {}) if isinstance(cocycle, dict) else {}
    producer_checks = cocycle.get("c3_exprel_producer4checks", {}) if isinstance(cocycle, dict) else {}
    if cocycle.get("source_gate", {}).get("all_passed") is not True:
        blockers.append("fresh cocycle source gate did not pass")
    if algebra.get("row_count") != 4 or algebra.get("all_rows_zero") is not True:
        blockers.append("fresh cocycle algebra four residuals are not all zero")
    if len(producer_checks.get("checks", ())) != 4 or producer_checks.get("all_checks_zero") is not True:
        blockers.append("fresh C3 exprel producer four checks are not all zero")
    if blockers:
        return _od2_base(source, blockers, native_source_binding=native, C3_source_binding=c3, integral_binding=integral_data, math_contract=math_contract, exponential_driver_cocycle=cocycle)
    result = _od2_base(
        source,
        [],
        passed=True,
        source_packet_defined=True,
        source_packet_definition_accepted=True,
        source_packet_accepted_owned=True,
        native_source_binding=native,
        C3_source_binding=c3,
        integral_binding=integral_data,
        atomrefs=atomrefs,
        actual_correlated_F_gap_functional_recipes=f_gaps,
        dF_atom_whole_function_views={"orders": list(OD2_ORDERS), "refs": [atomrefs["Fn_Z%s" % j] for j in OD2_ORDERS] + [atomrefs["Fc_Z%s" % j] for j in OD2_ORDERS], "source_values_evaluated": False},
        signed_density_gap_functions=density_gaps,
        cached_factorized_signed_defects={
            "formulas": cached_defects["formulas"],
            "unscaled_formulas": cached_defects.get("unscaled_formulas"),
            "rows": instantiated,
            "residual_zero_rows_data_only": True,
            "new_density_polynomial_checks": 0,
            "new_differentiation": 0,
        },
        channel_template_bindings={
            "e2": {"source_field": "e", "cached_formula": cached_defects["formulas"].get("e2"), "coefficient_denominator": 2, "source_physical_channel": True},
            "p2": {"source_field": "p", "cached_formula": cached_defects["formulas"].get("p2"), "coefficient_denominator": 2, "source_physical_channel": True},
            "h4": {"attached_F_gap_recipe": True, "source_values_evaluated": False},
        },
        doubled_internal_channel_denominators={"e2": {"physical_channel": "e", "denominator": 2}, "p2": {"physical_channel": "p", "denominator": 2}},
        exponential_driver_cocycle=cocycle,
        algebra4_residuals=algebra,
        c3_exprel_producer4checks=producer_checks,
        source_inspection_audit={
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_generic_identity_calls": 0,
            "old_operator_calls": 0,
            "old_source_consumer_calls": 0,
            "old_provider_calls": 0,
            "old_checker_calls": 0,
            "old_normalizer_calls": 0,
            "old_source_value_reducer_calls": 0,
            "old_numeric_calls": 0,
            "numeric_calls": 0,
            "quadratures": 0,
            "temporal_recurrence_calls": 0,
            "new_exponential_cocycle_order_checks": 4,
            "new_exprel_producer_order_checks": 4,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": 4,
            "actual_producer_order_failures": 0,
            "checked_residual_count": 4,
            "checked_residual_failures": 0,
            "new_exponential_cocycle_order_checks4": 4,
            "new_exprel_producer_order_checks4": 4,
            "new_actual_signed_density_gap_rows": 20,
            "new_actual_signed_density_gap_rows20": 20,
            "new_density_proof_rows": 0,
            "source_only": True,
            "metadata_only": False,
        },
        new_counts={"cocycle_checks": 4, "new_exprel_producer_order_checks": 4, "new_exprel_producer_order_failures": 0, "actual_producer_order_checks": 4, "actual_producer_order_failures": 0, "actual_driver_gap_rows": 20, "new_density_proof_rows": 0},
        counts={"new_exponential_cocycle_order_checks": 4, "new_exprel_producer_order_checks": 4, "new_actual_signed_density_gap_rows": 20, "ordinary_Z_orders": list(OD2_ORDERS)},
        preservation_constraints={
            "all_212_terms": True,
            "all_errors_retained": True,
            "caps_retained": True,
            "P1_not_m15": True,
            "physical_Q_not_loop_q": True,
            "ordinary_Z_not_temporal_n": True,
            "equal_roots_not_assumed": True,
            "primitive_gap_unknown": True,
            "P2_unknown": True,
            "phase_unknown": True,
            "native_C3_source_gap_zero_claimed": False,
        },
    )
    for name in OD2_FLAGS:
        result[name] = True
    result["new_flags"] = {name: True for name in OD2_FLAGS}
    _od2_emit(journal, "O3_transition_exponential_driver_cocycle_defined", passed=True, source_values_uncomputed=True)
    return result


__all__ = ("derive",)
