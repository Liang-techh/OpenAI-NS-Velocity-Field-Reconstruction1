"""Fresh C3/C5 source bridge and signed target/control/profile differences."""


from fractions import Fraction


cg_input_contract = {
    "source": {
        "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py",
        "sha256": "cb6f11d4a62d59a45f05d0961236d5c6236bd544dca8b18499f68b2d8d98d009",
        "function": "quadratic",
        "metadata_only": True,
    },
    "ordinary_Z_order_max": 3,
    "h_component_order": ["a0", "a2", "e0", "e1", "e2"],
    "h_C": "signed C3 control function jets, ordinary orders 0..3",
    "g": "signed increment h_N-h_C jets, ordinary orders 0..3",
    "target_gap": "signed eta=d_N-d_C jets, ordinary orders 0..3",
    "B": "5x5 signed Z-independent linear matrix with a cached inverse action",
    "W": {
        "cross": ["C0", "C2"],
        "energy": ["E0", "E1", "E2"],
        "pressure": ["P0", "P1", "P2"],
        "mu": "strictly positive and Z-independent",
    },
    "N": "strictly positive and Z-independent common control integer/scalar",
    "required_conditions": [
        "B, W, mu, and N have zero ordinary-Z derivatives",
        "h_N=h_C+g and eta=d_N-d_C use the same five-component order",
        "the residual equation is supplied on the same C1 ball",
        "the cached C1 contraction factor is at most 1/2",
    ],
    "admission": "conditional symbolic identity only; no numeric cap or midpoint is accepted",
}


cg_binomial = ((1,), (1, 1), (1, 2, 1), (1, 3, 3, 1))
cg_order_max = 3


def cg_fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    raise TypeError("cg coefficients require int or Fraction")


def cg_poly_atom(name):
    return {(str(name),): Fraction(1, 1)}


def cg_poly_constant(value):
    value = cg_fraction(value)
    return {(): value} if value else {}


def cg_poly_add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0, 1)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def cg_poly_neg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def cg_poly_sub(left, right):
    return cg_poly_add(left, cg_poly_neg(right))


def cg_poly_scale(poly, coefficient):
    coefficient = cg_fraction(coefficient)
    return {
        monomial: coefficient * value
        for monomial, value in poly.items()
        if coefficient * value
    }


def cg_poly_mul(left, right):
    out = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial))
            out[monomial] = out.get(monomial, Fraction(0, 1)) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def cg_jet_zero():
    return ({}, {}, {}, {})


def cg_jet_atom(name):
    return tuple(cg_poly_atom("%s_Z%d" % (name, order)) for order in range(cg_order_max + 1))


def cg_jet_constant_atom(name):
    return (cg_poly_atom(name), {}, {}, {})


def cg_jet_one():
    return (cg_poly_constant(1), {}, {}, {})


def cg_jet_add(*jets):
    if not jets:
        return cg_jet_zero()
    return tuple(cg_poly_add(*(jet[order] for jet in jets)) for order in range(cg_order_max + 1))


def cg_jet_neg(jet):
    return tuple(cg_poly_neg(value) for value in jet)


def cg_jet_sub(left, right):
    return cg_jet_add(left, cg_jet_neg(right))


def cg_jet_scale(jet, scalar):
    return tuple(cg_poly_mul(scalar, value) for value in jet)


def cg_jet_mul(left, right):
    """Ordinary-Z product rule, including all binomial coefficients through k=3."""
    values = []
    for order in range(cg_order_max + 1):
        values.append(
            cg_poly_add(
                *[
                    cg_poly_scale(
                        cg_poly_mul(left[left_order], right[order - left_order]),
                        cg_binomial[order][left_order],
                    )
                    for left_order in range(order + 1)
                ]
            )
        )
    return tuple(values)


def cg_vector_add(left, right):
    return tuple(cg_jet_add(a, b) for a, b in zip(left, right))


def cg_vector_sub(left, right):
    return tuple(cg_jet_sub(a, b) for a, b in zip(left, right))


def cg_make_weights():
    return {
        "C": [cg_jet_constant_atom("cg_C0")[0], cg_jet_constant_atom("cg_C2")[0]],
        "E": [cg_jet_constant_atom("cg_E%d" % i)[0] for i in range(3)],
        "P": [cg_jet_constant_atom("cg_P%d" % i)[0] for i in range(3)],
        "mu_inverse": cg_jet_constant_atom("cg_mu_inverse")[0],
        "N_inverse": cg_jet_constant_atom("cg_N_inverse")[0],
    }


def cg_make_matrix():
    return [
        [cg_jet_constant_atom("cg_B_%d_%d" % (row, col))[0] for col in range(5)]
        for row in range(5)
    ]


def cg_linear_action(matrix, vector):
    return tuple(
        cg_jet_add(*(cg_jet_scale(vector[col], matrix[row][col]) for col in range(5)))
        for row in range(5)
    )


def cg_quadratic(h, weights):
    """Fresh sparse form of the five-row native quadratic Q."""
    a0, a2, e0, e1, e2 = h
    c0, c2 = weights["C"]
    energy = weights["E"]
    pressure = weights["P"]
    cross = cg_jet_scale(cg_jet_add(cg_jet_scale(cg_jet_mul(a0, e0), c0),
                                    cg_jet_scale(cg_jet_mul(a2, e2), c2)),
                         weights["mu_inverse"])
    energy_square = cg_jet_add(
        cg_jet_scale(cg_jet_mul(a0, a0), energy[0]),
        cg_jet_scale(cg_jet_mul(a2, a2), energy[2]),
    )
    error_square = cg_jet_add(
        cg_jet_scale(cg_jet_mul(e0, e0), energy[0]),
        cg_jet_scale(cg_jet_mul(e1, e1), energy[1]),
        cg_jet_scale(cg_jet_mul(e2, e2), energy[2]),
    )
    pressure_square = cg_jet_add(
        cg_jet_scale(cg_jet_mul(e0, e0), pressure[0]),
        cg_jet_scale(cg_jet_mul(e1, e1), pressure[1]),
        cg_jet_scale(cg_jet_mul(e2, e2), pressure[2]),
    )
    half = cg_poly_constant(Fraction(1, 2))
    return (
        cg_jet_zero(),
        cross,
        cg_jet_zero(),
        cg_jet_sub(energy_square, cg_jet_scale(error_square, half)),
        cg_jet_scale(pressure_square, half),
    )


def cg_quadratic_action(h, direction, weights):
    """The bilinear DQ(h)direction used by the residual equation."""
    a0, a2, e0, e1, e2 = h
    b0, b2, f0, f1, f2 = direction
    c0, c2 = weights["C"]
    energy = weights["E"]
    pressure = weights["P"]
    cross = cg_jet_scale(
        cg_jet_add(
            cg_jet_scale(cg_jet_add(cg_jet_mul(e0, b0), cg_jet_mul(a0, f0)), c0),
            cg_jet_scale(cg_jet_add(cg_jet_mul(e2, b2), cg_jet_mul(a2, f2)), c2),
        ),
        weights["mu_inverse"],
    )
    energy_action = cg_jet_scale(
        cg_jet_add(
            cg_jet_scale(cg_jet_mul(a0, b0), energy[0]),
            cg_jet_scale(cg_jet_mul(a2, b2), energy[2]),
        ),
        cg_poly_constant(2),
    )
    error_action = cg_jet_add(
        cg_jet_scale(cg_jet_mul(e0, f0), energy[0]),
        cg_jet_scale(cg_jet_mul(e1, f1), energy[1]),
        cg_jet_scale(cg_jet_mul(e2, f2), energy[2]),
    )
    pressure_action = cg_jet_add(
        cg_jet_scale(cg_jet_mul(e0, f0), pressure[0]),
        cg_jet_scale(cg_jet_mul(e1, f1), pressure[1]),
        cg_jet_scale(cg_jet_mul(e2, f2), pressure[2]),
    )
    return (
        cg_jet_zero(),
        cross,
        cg_jet_zero(),
        cg_jet_sub(energy_action, error_action),
        pressure_action,
    )


def cg_signed_source_jets():
    names = ("a0", "a2", "e0", "e1", "e2")
    return {
        "h_C": tuple(cg_jet_atom("cg_hC_%s" % name) for name in names),
        "g": tuple(cg_jet_atom("cg_g_%s" % name) for name in names),
        "target_gap": tuple(cg_jet_atom("cg_eta_%s" % name) for name in names),
    }


def cg_zero_residual_report(left, right):
    report = []
    for row, (left_jet, right_jet) in enumerate(zip(left, right)):
        for order in range(cg_order_max + 1):
            difference = cg_poly_sub(left_jet[order], right_jet[order])
            report.append(
                {
                    "row": row,
                    "ordinary_Z_order": order,
                    "lhs_term_count": len(left_jet[order]),
                    "rhs_term_count": len(right_jet[order]),
                    "residual_monomials": [
                        {
                            "monomial": list(monomial),
                            "coefficient": {"numerator": value.numerator, "denominator": value.denominator},
                        }
                        for monomial, value in sorted(difference.items())
                    ],
                    "passed": not difference,
                }
            )
    return report


def cg_encode_poly(poly):
    return [
        {
            "monomial": list(monomial),
            "coefficient": {"numerator": value.numerator, "denominator": value.denominator},
        }
        for monomial, value in sorted(poly.items())
    ]


def cg_encode_jet(jet):
    return [
        {"ordinary_Z_order": order, "terms": cg_encode_poly(jet[order])}
        for order in range(cg_order_max + 1)
    ]


def cg_encode_vector(vector):
    return [cg_encode_jet(jet) for jet in vector]


def cg_product_tensor():
    return [
        {
            "ordinary_Z_order": order,
            "terms": [
                {
                    "left_order": left_order,
                    "right_order": order - left_order,
                    "coefficient": {
                        "numerator": cg_binomial[order][left_order],
                        "denominator": 1,
                    },
                }
                for left_order in range(order + 1)
            ],
        }
        for order in range(cg_order_max + 1)
    ]


def cg_contract_report(contract):
    contract = contract if isinstance(contract, dict) else {}
    return {
        "same_Z_independent_B_W_N": contract.get("same_Z_independent_B_W_N") is True,
        "same_C1_ball": contract.get("same_C1_ball") is True,
        "cached_contraction_at_most_half": contract.get("cached_contraction_at_most_half") is True,
        "residual_equation_supplied": contract.get("residual_equation_supplied") is True,
        "source_binding_admitted": False,
        "numeric_cap_admitted": False,
    }


def cg_derive_control_gap(contract=None, journal=None):
    """Return a JSON-safe signed control-gap expansion and conditional bound."""
    source = cg_signed_source_jets()
    weights = cg_make_weights()
    matrix = cg_make_matrix()
    h_c = source["h_C"]
    increment = source["g"]
    target_gap = source["target_gap"]
    h_n = cg_vector_add(h_c, increment)

    q_c = cg_quadratic(h_c, weights)
    q_n = cg_quadratic(h_n, weights)
    q_g = cg_quadratic(increment, weights)
    dq_g = cg_quadratic_action(h_c, increment, weights)
    expansion_left = cg_vector_sub(q_n, q_c)
    expansion_right = cg_vector_add(dq_g, q_g)
    expansion_checks = cg_zero_residual_report(expansion_left, expansion_right)

    linear_g = cg_linear_action(matrix, increment)
    nonlinear = cg_vector_add(dq_g, q_g)
    nonlinear_over_n = tuple(cg_jet_scale(value, weights["N_inverse"]) for value in nonlinear)
    residual = cg_vector_add(cg_vector_add(linear_g, target_gap), nonlinear_over_n)
    contract_report = cg_contract_report(contract)
    if journal is not None:
        journal(
            "signed_control_gap_tensor",
            {
                "quadratic_expansion_passed": all(row["passed"] for row in expansion_checks),
                "source_binding_admitted": False,
                "numeric_cap_admitted": False,
            },
        )

    return {
        "passed": all(row["passed"] for row in expansion_checks),
        "conditional": True,
        "ordinary_Z_order_max": cg_order_max,
        "input_contract": cg_input_contract,
        "contract_report": contract_report,
        "signed_ports": {
            "h_C": "retained as signed ordinary-Z jets",
            "g": "retained as signed ordinary-Z jets; never forced to zero",
            "target_gap": "retained as signed eta=d_N-d_C jets; never replaced by a bound",
        },
        "product_coefficient_tensor": cg_product_tensor(),
        "quadratic_expansion": {
            "statement": "Q(h_C+g)-Q(h_C)=DQ(h_C)g+Q(g)",
            "per_monomial_checks": expansion_checks,
            "passed": all(row["passed"] for row in expansion_checks),
        },
        "control_residual": {
            "statement": "B*g+target_gap+(DQ(h_C)g+Q(g))/N=0",
            "symbolic_residual": cg_encode_vector(residual),
            "equation_is_input_condition": True,
            "residual_forced_zero": False,
        },
        "conditional_contraction": {
            "assumptions": [
                "same C1 ball for h_C and h_N",
                "cached ||B^{-1} DQ(h)||_C1/N <= 1/2 on that ball",
                "B is invertible with the supplied Z-independent inverse action",
                "the signed control residual equation holds in C1",
            ],
            "bound": "||g||_C1 <= 2 ||B^{-1} target_gap||_C1",
            "zero_gap_consequence": "target_gap=0 implies g=0 and h_N,h_C ordinary Z orders 0..3 agree",
            "numeric_evaluation": False,
        },
        "quadratic_rows": [
            "row0=0",
            "row1=(C0*a0*e0+C2*a2*e2)/mu",
            "row2=0",
            "row3=E0*a0^2+E2*a2^2-(E0*e0^2+E1*e1^2+E2*e2^2)/2",
            "row4=(P0*e0^2+P1*e1^2+P2*e2^2)/2",
        ],
        "source_binding_admitted": False,
        "numeric_cap_admitted": False,
        "old_math_calls": 0,
        "old_source_evaluations": 0,
        "new_numeric_evaluations": 0,
    }


def cg_input_contract_schema():
    return cg_input_contract




from fractions import Fraction


pd_input_contract = {
    "ordinary_Z_order_max": 3,
    "radial_y_orders": [0, 1, 2],
    "amplitude": {
        "A_C": "signed common-profile amplitude jet",
        "a": "signed amplitude difference A_N-A_C",
        "A_N": "defined as A_C+a",
    },
    "profile": {
        "h_C": "five signed common-profile jets [h0,h1,h2,h3,h4]",
        "g": "signed profile difference h_N-h_C",
        "h_N": "defined componentwise as h_C+g",
    },
    "common_physical_bumps": {
        "same_for_N_and_C": True,
        "ordinary_Z_independent": True,
        "y_orders": [0, 1, 2],
        "centers": [0, 1, 2],
    },
    "N": "strictly positive and Z-independent; one inverse-N factor per F/G",
    "alpha_and_power": "common Z-independent source factors for original_E",
    "index_sets": {
        "deltaV": [[0, 0], [2, 1]],
        "deltaE": [[0, 2], [1, 3], [2, 4]],
        "pair_order": "[bump_center, h_index]",
    },
    "source_binding": "conditional tensor identity only; no cap admission",
}


pd_binomial = ((1,), (1, 1), (1, 2, 1), (1, 3, 3, 1))
pd_order_max = 3
pd_deltaV_indices = ((0, 0), (2, 1))
pd_deltaE_indices = ((0, 2), (1, 3), (2, 4))


def pd_fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    raise TypeError("pd coefficients require int or Fraction")


def pd_poly_atom(name):
    return {(str(name),): Fraction(1, 1)}


def pd_poly_constant(value):
    value = pd_fraction(value)
    return {(): value} if value else {}


def pd_poly_add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0, 1)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def pd_poly_neg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def pd_poly_sub(left, right):
    return pd_poly_add(left, pd_poly_neg(right))


def pd_poly_scale(poly, coefficient):
    coefficient = pd_fraction(coefficient)
    return {
        monomial: coefficient * value
        for monomial, value in poly.items()
        if coefficient * value
    }


def pd_poly_mul(left, right):
    out = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial))
            out[monomial] = out.get(monomial, Fraction(0, 1)) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def pd_jet_zero():
    return ({}, {}, {}, {})


def pd_jet_atom(name):
    return tuple(pd_poly_atom("%s_Z%d" % (name, order)) for order in range(pd_order_max + 1))


def pd_jet_constant_atom(name):
    return (pd_poly_atom(name), {}, {}, {})


def pd_jet_one():
    return (pd_poly_constant(1), {}, {}, {})


def pd_jet_add(*jets):
    if not jets:
        return pd_jet_zero()
    return tuple(pd_poly_add(*(jet[order] for jet in jets)) for order in range(pd_order_max + 1))


def pd_jet_neg(jet):
    return tuple(pd_poly_neg(value) for value in jet)


def pd_jet_sub(left, right):
    return pd_jet_add(left, pd_jet_neg(right))


def pd_jet_scale(jet, scalar):
    return tuple(pd_poly_mul(scalar, value) for value in jet)


def pd_jet_mul(left, right):
    """Ordinary-Z product rule through order three."""
    return tuple(
        pd_poly_add(
            *[
                pd_poly_scale(
                    pd_poly_mul(left[left_order], right[order - left_order]),
                    pd_binomial[order][left_order],
                )
                for left_order in range(order + 1)
            ]
        )
        for order in range(pd_order_max + 1)
    )


def pd_jet_pow(jet, exponent):
    if exponent < 0:
        raise ValueError("pd negative powers are not part of this finite tensor")
    result = pd_jet_one()
    for unused in range(exponent):
        result = pd_jet_mul(result, jet)
    return result


def pd_vector_add(left, right):
    return tuple(pd_jet_add(a, b) for a, b in zip(left, right))


def pd_make_symbols():
    h_names = ("h0", "h1", "h2", "h3", "h4")
    return {
        "A_C": pd_jet_atom("pd_A_C"),
        "a": pd_jet_atom("pd_a"),
        "h_C": tuple(pd_jet_atom("pd_hC_%s" % name) for name in h_names),
        "g": tuple(pd_jet_atom("pd_g_%s" % name) for name in h_names),
        "alpha": pd_jet_constant_atom("pd_alpha"),
        "power": pd_jet_constant_atom("pd_power"),
        "N_inverse": pd_jet_constant_atom("pd_N_inverse")[0],
        # The physical bump functions have no ordinary-Z dependence.
        "bump": tuple(
            tuple(pd_jet_constant_atom("pd_bump_y%d_c%d" % (y, center)) for center in range(3))
            for y in range(3)
        ),
    }


def pd_index_sum(profile, bump_row, indices):
    return pd_jet_add(
        *(pd_jet_mul(profile[h_index], bump_row[bump_center]) for bump_center, h_index in indices)
    )


def pd_native_full(amplitude, profile, bump_row, indices, inverse_N):
    # This is the full product for one state, kept separate from the
    # three-term difference construction below.
    return pd_jet_scale(pd_jet_mul(amplitude, pd_index_sum(profile, bump_row, indices)), inverse_N)


def pd_three_term_difference(amplitude, amplitude_difference, profile, profile_difference,
                             bump_row, indices, inverse_N):
    # A fresh expansion of A*g + a*h_C + a*g, independent of pd_native_full.
    common_h = pd_index_sum(profile, bump_row, indices)
    delta_h = pd_index_sum(profile_difference, bump_row, indices)
    return pd_jet_scale(
        pd_jet_add(
            pd_jet_mul(amplitude, delta_h),
            pd_jet_mul(amplitude_difference, common_h),
            pd_jet_mul(amplitude_difference, delta_h),
        ),
        inverse_N,
    )


def pd_original_E_difference(amplitude_difference, alpha, power, y_order):
    fixed = pd_jet_mul(power, pd_jet_pow(alpha, y_order))
    if y_order % 2:
        fixed = pd_jet_neg(fixed)
    return pd_jet_mul(amplitude_difference, fixed)


def pd_zero_report(left, right):
    report = []
    for order in range(pd_order_max + 1):
        difference = pd_poly_sub(left[order], right[order])
        report.append(
            {
                "ordinary_Z_order": order,
                "lhs_term_count": len(left[order]),
                "rhs_term_count": len(right[order]),
                "residual_monomials": [
                    {
                        "monomial": list(monomial),
                        "coefficient": {"numerator": value.numerator, "denominator": value.denominator},
                    }
                    for monomial, value in sorted(difference.items())
                ],
                "passed": not difference,
            }
        )
    return report


def pd_encode_poly(poly):
    return [
        {
            "monomial": list(monomial),
            "coefficient": {"numerator": value.numerator, "denominator": value.denominator},
        }
        for monomial, value in sorted(poly.items())
    ]


def pd_encode_jet(jet):
    return [
        {"ordinary_Z_order": order, "terms": pd_encode_poly(jet[order])}
        for order in range(pd_order_max + 1)
    ]


def pd_encode_rows(rows):
    return [pd_encode_jet(row) for row in rows]


def pd_product_tensor():
    return [
        {
            "ordinary_Z_order": order,
            "terms": [
                {
                    "left_order": left_order,
                    "right_order": order - left_order,
                    "coefficient": {
                        "numerator": pd_binomial[order][left_order],
                        "denominator": 1,
                    },
                }
                for left_order in range(order + 1)
            ],
        }
        for order in range(pd_order_max + 1)
    ]


def pd_derive_profile_difference(contract=None, journal=None):
    """Construct the signed profile difference and its independent expansion."""
    symbols = pd_make_symbols()
    amplitude = symbols["A_C"]
    amplitude_difference = symbols["a"]
    profile = symbols["h_C"]
    profile_difference = symbols["g"]
    amplitude_native = pd_jet_add(amplitude, amplitude_difference)
    profile_native = pd_vector_add(profile, profile_difference)
    inverse_N = symbols["N_inverse"]

    rows = []
    checks = []
    for y_order in range(3):
        bump_row = symbols["bump"][y_order]
        native_v = pd_jet_sub(
            pd_native_full(amplitude_native, profile_native, bump_row, pd_deltaV_indices, inverse_N),
            pd_native_full(amplitude, profile, bump_row, pd_deltaV_indices, inverse_N),
        )
        expanded_v = pd_three_term_difference(
            amplitude, amplitude_difference, profile, profile_difference,
            bump_row, pd_deltaV_indices, inverse_N,
        )
        native_e = pd_jet_sub(
            pd_native_full(amplitude_native, profile_native, bump_row, pd_deltaE_indices, inverse_N),
            pd_native_full(amplitude, profile, bump_row, pd_deltaE_indices, inverse_N),
        )
        expanded_e = pd_three_term_difference(
            amplitude, amplitude_difference, profile, profile_difference,
            bump_row, pd_deltaE_indices, inverse_N,
        )
        original_e = pd_original_E_difference(
            amplitude_difference, symbols["alpha"], symbols["power"], y_order,
        )
        rows.append(
            {
                "radial_y_order": y_order,
                "deltaV_native_difference": native_v,
                "deltaV_three_term": expanded_v,
                "deltaE_native_difference": native_e,
                "deltaE_three_term": expanded_e,
                "original_E_difference": original_e,
            }
        )
        checks.append(
            {
                "radial_y_order": y_order,
                "deltaV": pd_zero_report(native_v, expanded_v),
                "deltaE": pd_zero_report(native_e, expanded_e),
            }
        )

    passed = all(
        item["passed"]
        for row in checks
        for field in ("deltaV", "deltaE")
        for item in row[field]
    )
    if journal is not None:
        journal(
            "signed_profile_difference_tensor",
            {
                "independent_three_term_expansion_passed": passed,
                "source_binding_admitted": False,
                "numeric_cap_admitted": False,
            },
        )
    return {
        "passed": passed,
        "conditional": True,
        "ordinary_Z_order_max": pd_order_max,
        "radial_y_orders": [0, 1, 2],
        "input_contract": pd_input_contract,
        "signed_difference_ports": {
            "a": "A_N-A_C, signed and retained",
            "g": "h_N-h_C, signed and retained componentwise",
            "common_bumps": "shared Z-independent physical bump jets",
        },
        "index_sets": {
            "deltaV": [list(pair) for pair in pd_deltaV_indices],
            "deltaE": [list(pair) for pair in pd_deltaE_indices],
        },
        "product_coefficient_tensor": pd_product_tensor(),
        "rows": [
            {
                "radial_y_order": row["radial_y_order"],
                "deltaV_native_difference": pd_encode_jet(row["deltaV_native_difference"]),
                "deltaV_three_term": pd_encode_jet(row["deltaV_three_term"]),
                "deltaE_native_difference": pd_encode_jet(row["deltaE_native_difference"]),
                "deltaE_three_term": pd_encode_jet(row["deltaE_three_term"]),
                "original_E_difference": pd_encode_jet(row["original_E_difference"]),
            }
            for row in rows
        ],
        "independent_expansion_checks": checks,
        "identity": {
            "statement": "(A_C+a)(h_C+g)-A_C*h_C = A_C*g+a*h_C+a*g, divided by N",
            "ordinary_Z_orders_checked": [0, 1, 2, 3],
            "passed": passed,
        },
        "zero_implication": {
            "condition": "a and every component of g vanish as functions through ordinary Z order 3",
            "result": "deltaV, deltaE, and original_E differences vanish for y=0,1,2",
            "actual_source_equality_admitted": False,
        },
        "source_binding_admitted": False,
        "numeric_cap_admitted": False,
        "old_math_calls": 0,
        "old_source_evaluations": 0,
        "new_numeric_evaluations": 0,
    }


def pd_input_contract_schema():
    return pd_input_contract



def sb_assert(condition,message):
    if not condition:raise ValueError(message)


def sb_literal_pair_proof(left_slice,right_slice,schema,root_pairs,label):
    """New finite DAG induction over typed operands; no algebra normalization."""
    maps=[]
    for side in (left_slice,right_slice):
        sb_assert(side['closure_complete'] is True,'Open typed source slice: '+label)
        mapping={entry['id']:entry for entry in side['nodes']}
        sb_assert(len(mapping)==len(side['nodes']),'Duplicate source node')
        maps.append(mapping)
    pending=[(a,b,[]) for a,b in root_pairs];seen=set();checked=[];failures=[]
    while pending:
        a,b,path=pending.pop()
        if (a,b) in seen:continue
        seen.add((a,b))
        if a not in maps[0] or b not in maps[1]:
            failures.append({'left':a,'right':b,'path':path,'reason':'missing typed operand'});continue
        entries=[maps[0][a],maps[1][b]];rows=[entry['row'] for entry in entries]
        op=rows[0].get('operation')
        if op!=rows[1].get('operation') or op not in schema:
            failures.append({'left':a,'right':b,'path':path,'reason':'operation mismatch or no operand semantics',
                             'left_operation':op,'right_operation':rows[1].get('operation')});continue
        fields=schema[op]
        payload=[{k:v for k,v in row.items() if k not in fields} for row in rows]
        if payload[0]!=payload[1]:
            failures.append({'left':a,'right':b,'path':path,'reason':'literal defining payload differs',
                             'left_payload':payload[0],'right_payload':payload[1]});continue
        refs=[];valid=True
        for field in fields:
            present=[field in row for row in rows]
            if present[0]!=present[1]:valid=False;break
            if not present[0]:continue
            values=[row[field] for row in rows]
            if any(value is None for value in values):
                if values!=[None,None]:valid=False;break
                for entry in entries:
                    sb_assert(field not in entry['typed_reference_fields'],'Null optional field has a typed operand')
                continue
            for entry,value in zip(entries,values):
                sb_assert(entry['typed_reference_fields'].get(field)==value,'Typed operand disagrees with raw row')
            arrays=[v if isinstance(v,list) else [v] for v in values]
            if len(arrays[0])!=len(arrays[1]):valid=False;break
            for index,(x,y) in enumerate(zip(*arrays)):
                if type(x) is not int or type(y) is not int or not 0<=x<a or not 0<=y<b:
                    valid=False;break
                refs.append((x,y,path+[field+'['+str(index)+']']))
            if not valid:break
        if not valid:
            failures.append({'left':a,'right':b,'path':path,'reason':'operand shape or DAG order differs'});continue
        pending.extend(refs);checked.append([a,b])
    return {'schema':'typed_source_DAG_induction_v1','label':label,'proved':not failures,
            'root_pairs':root_pairs,'paired_nodes_checked':len(checked),'visited_pair_count':len(seen),
            'pair_witnesses':checked,'first_mismatches':failures[:12],'mismatch_count':len(failures),
            'argument_reordering_or_rewriting':False,'range_overlap_used':False,
            'same_family_used_as_identity':False,'proof_rule':'Equal defining payload and equal typed operand meanings imply equal node functions; finite backward DAG induction.'}


def sb_signed_graph(ports):
    """Define actual signed source gaps, preserving distinct archive namespaces."""
    graph=[]
    def node(operation,**fields):
        graph.append({'id':len(graph),'operation':operation,**fields});return len(graph)-1
    hashes={'C3':ports['C3_archive_raw_sha256'],'C5':ports['C5_archive_raw_sha256']}
    def source(namespace,source_node,order,quantity):
        return node('authenticated_source_function',namespace=namespace,archive_raw_sha256=hashes[namespace],
                    source_node=source_node,ordinary_Z_order=order,quantity=quantity,
                    numeric_enclosure_not_used=True)
    def difference(n,c):return node('signed_function_difference',native=n,C3=c)
    def jets(namespace,ids,quantity):
        sb_assert(len(ids)==4 and all(type(i) is int for i in ids),'Wrong source jet shape')
        return [source(namespace,i,k,quantity) for k,i in enumerate(ids)]
    target={};incoming={};scaled={};control={}
    for name in ('M','D=(J-M)/mu','I','S','Cp'):
        n=jets('C5',ports['C5_targets'][name],'target_'+name)
        c=jets('C3',ports['C3_targets'][name],'target_'+name)
        target[name]=[difference(a,b) for a,b in zip(n,c)]
    for name in ('m','h','k','e','p'):
        ni=jets('C5',ports['C5_incoming'][name],'incoming_'+name)
        ci=jets('C3',ports['C3_incoming'][name],'incoming_'+name)
        nt=jets('C5',ports['C5_scaled_incoming'][name],'N_scaled_'+name)
        ct=jets('C3',ports['C3_terminal_scaled'][name],'N_scaled_'+name)
        incoming[name]=[difference(a,b) for a,b in zip(ni,ci)]
        scaled[name]=[difference(a,b) for a,b in zip(nt,ct)]
    for i,(nids,cids) in enumerate(zip(ports['C5_controls'],ports['C3_controls'])):
        control[str(i)]=[difference(a,b) for a,b in zip(jets('C5',nids,'h_'+str(i)),jets('C3',cids,'h_'+str(i)))]
    amplitude=[difference(a,b) for a,b in zip(jets('C5',ports['C5_amplitude'],'A'),jets('C3',ports['C3_amplitude'],'A'))]
    return {'schema':'actual_C3_C5_signed_source_gap_graph_v1','nodes':graph,'node_count':len(graph),
            'target_gap_C3':target,'raw_incoming_gap_C3':incoming,'N_scaled_incoming_gap_C3':scaled,
            'control_gap_C3':control,'amplitude_gap_C3':amplitude,
            'ordinary_Z_orders':[0,1,2,3],'gap_values_selected_or_evaluated':False,
            'normalization_formulas':{
                'M':'eta_M=(A_C*T_N_m-A_N*T_C_m)/(A_N*A_C)',
                'I':'eta_I=(A_C*T_N_h-A_N*T_C_h)/(A_N*A_C)',
                'S':'eta_S=(A_C^2*T_N_e-A_N^2*T_C_e)/(A_N^2*A_C^2)',
                'Cp':'eta_Cp=(A_C^2*T_N_p-A_N^2*T_C_p)/(A_N^2*A_C^2)',
                'D=(J-M)/mu':'eta_D=[A_C^2*(T_N_k-A_N*T_N_m)-A_N^2*(T_C_k-A_C*T_C_m)]/(mu*A_N^2*A_C^2)'},
            'formula_premises':['Both cached producer target-normalization definitions','A_N>0, A_C>0 on the common admitted source domain','common exact positive mu and N'],
            'derivative_rule':'The displayed source-node differences are the actual ordinary-Z derivatives by linearity through order 3.',
            'physical_Pstar_added_to_normalized_targets':False}


def sb_normalizations(ports,literal_rows,common_N):
    sb_assert(common_N,'The two integer source definitions have not been identified')
    c3=literal_rows['C3'];c5=literal_rows['C5'];proof=[]
    for namespace,rows in (('C3',c3),('C5',c5)):
        one=rows['1'];integer=rows['18']
        sb_assert(one=={'operation':'exact_rational','numerator':1,'denominator':1},'Integer-one port changed')
        sb_assert(integer['exact_positive_integer'] and integer['Z_independent'] and
                  integer['J']==1358356628656378313,'Selected N source changed')
    for name in ('m','h','k','e','p'):
        for k in range(4):
            ni=ports['C5_incoming'][name][k];nt=ports['C5_scaled_incoming'][name][k]
            ci=ports['C3_incoming'][name][k];ct=ports['C3_terminal_scaled'][name][k]
            nrow=c5[str(nt)];crow=c3[str(ci)]
            sb_assert(nrow['operation']=='product' and nrow['arguments'] in ([ni,18],[18,ni]),'Native N-scaled incoming recipe changed')
            sb_assert(crow['operation']=='product','C3 incoming is not the prescribed invN product')
            args=crow['arguments'];sb_assert(len(args)==2 and ct in args,'Wrong C3 terminal source port')
            inv=args[1] if args[0]==ct else args[0];reciprocal=c3[str(inv)]
            sb_assert(reciprocal['operation']=='positive_quotient' and reciprocal['numerator']==1 and reciprocal['denominator']==18,'Wrong inverse-N source')
            proof.append({'history':name,'ordinary_Z_order':k,'C3_incoming':ci,'C3_N_scaled_history':ct,
                          'C5_incoming':ni,'C5_N_scaled_history':nt,'C3_invN_node':inv,
                          'C3_identity':'N*I_C=T_C','native_identity':'N*I_N=T_N',
                          'signed_gap_identity':'N*(I_N-I_C)=T_N-T_C','extra_N_divisions':0,
                          'incoming_gap_zero_proved':False})
    return {'proved':True,'history_order_identities':proof,'ordinary_Z_orders':[0,1,2,3],
            'count':len(proof),'normalization_is_function_identity':True,
            'native_C3_incoming_equality_claimed':False,'N_materialized':False}


def sb_port_contract(ports,anchors,controls):
    rows=('M','D=(J-M)/mu','I','S','Cp')
    for namespace in ('C3','C5'):
        for name in ('targets','incoming','controls','amplitude'):
            values=ports[namespace+'_'+name]
            arrays=list(values.values()) if isinstance(values,dict) else values if name=='controls' else [values]
            for ids in arrays:
                sb_assert(len(ids)==4 and all(type(i) is int and str(i) in anchors[namespace] for i in ids),'Malformed or missing source jet: '+namespace+' '+name)
    for k,i in enumerate(ports['C3_amplitude']):
        row=anchors['C3'][str(i)]
        sb_assert(row['operation']=='current_original_leading_function_recipe' and row['quantity']=='E' and
            row['native_chart']=='O3_power' and row['coordinate']==48 and row['Z_variable']==20 and
            row['Z_order']==k and row['defining_quantity_not_a_range_value'] and row['terminal_positive_A_rc'],
            'C3 terminal amplitude has the wrong source role')
        if k>=2:sb_assert(row['Taylor_coefficient_factorial']==(2 if k==2 else 6) and row['derivative_of_range_endpoint'] is False,'Wrong ordinary amplitude derivative units')
    for k,i in enumerate(ports['C5_amplitude']):
        row=anchors['C5'][str(i)]
        sb_assert(row['operation']=='function_substitution' and row['variable']==250215 and row['value']==48 and
            row['Z_independent_substitution'] is True,'C5 amplitude endpoint role changed')
    for component,(cids,nids) in enumerate(zip(ports['C3_controls'],ports['C5_controls'])):
        c=anchors['C3'][str(cids[0])];n=anchors['C5'][str(nids[0])];third=anchors['C3'][str(cids[3])]
        sb_assert(c['operation']=='C1_Banach_limit_component' and c['component']==component and c['Z_order']==0,'C3 control value source changed')
        sb_assert(third['operation']=='C3_Banach_limit_third_component' and third['component']==component and
            third['Z_order']==3 and third['genuine_ordinary_third_derivative'] and third['derivative_of_range_endpoint'] is False,'C3 third control jet changed')
        sb_assert(n['operation']=='original_compact_unique_Banach_limit_component' and n['component']==component and
            n['system']==259187 and n['Z_variable']==20 and n['mathematical_function_not_numeric_interval'],
            'C5 control function source changed')
    sb_assert(tuple(ports['C5_targets'])==rows,'Native target row order changed')
    return {'passed':True,'ordinary_Z_orders':[0,1,2,3],'amplitude_history_roles_separated':True,
            'control_values_are_source_functions':True,'source_interval_endpoint_selected':False,
            'provider_ordinary_Z_contract':controls['C5_implicit_J_inverse']['ordinary_Z_contract']}


def derive(source,journal):
    spec=source['input'];sb_assert(spec['schema']=='compact_C3_C5_source_bridge_input_v1','Wrong fresh input schema')
    receipts=spec['receipts'];ports=spec['ports'];parent=source['parent_result']
    authenticated=receipts['weights']['authenticated_records']
    sb_assert(ports['C3_archive_raw_sha256']==authenticated['c3_repair_band_archive']['raw_sha256'] and
              ports['C5_archive_raw_sha256']==authenticated['true_compact_C5_archive']['raw_sha256'],
              'Source gap archive namespaces do not match the authenticated weight sources')
    port_contract=sb_port_contract(ports,spec['source_port_literals'],receipts['controls'])
    sb_assert(parent['passed'] and parent['component']=='compact_prefix_bound_table_certificate','Wrong cached parent')
    weights=receipts['weights'];closure=weights['bounded_typed_closures']
    roots=[]
    for value in weights['weight_root_selectors']['root_ids'].values():
        roots.extend(value if isinstance(value,list) else [value])
    common=sb_literal_pair_proof(closure['c3_repair_band'],closure['true_compact_C5'],
                                weights['typed_reference_schema'],[[i,i] for i in roots],'exact common weights and N')
    sb_assert(common['proved'],'Exact weight-source comparison failed')
    for namespace in ('c3_repair_band','true_compact_C5'):
        for entry in closure[namespace]['nodes']:
            sb_assert(entry['operation'] in ('exact_rational','exact_positive_integer_expression','analytic_unary','negative',
                'sum','product','positive_quotient','bound_variable','compact_raw_beta','definite_integral'),
                'Unexpected Z-dependent source in weight closure')
            if entry['operation']=='bound_variable':
                sb_assert(entry['row']['name']=='repair_normalized_bump_t','Weight depends on an unbound axial variable')
    common['ordinary_Z_independent']=True
    journal({'event':'fresh_common_weight_source_identity','paired_nodes':common['paired_nodes_checked']})
    remap_receipt=receipts['C3_remap'];remap_slices=remap_receipt['typed_closed_nodes']
    all_pairs=spec['C3_remap_root_pairs']
    remap=sb_literal_pair_proof(remap_slices['standalone_C3_source'],remap_slices['repair_C3_band'],
                               remap_receipt['typed_reference_schema'],all_pairs,'C3 standalone to rebuilt source namespace')
    journal({'event':'fresh_C3_namespace_identity_attempt','proved':remap['proved'],
             'paired_nodes':remap['paired_nodes_checked'],'mismatches':remap['mismatch_count']})
    normalization=sb_normalizations(ports,spec['normalization_literal_rows'],common['proved'])
    gap_graph=sb_signed_graph(ports)
    cg=cg_derive_control_gap({'same_Z_independent_B_W_N':True,'same_C1_ball':False,
                              'cached_contraction_at_most_half':True,'residual_equation_supplied':True},
                             lambda event,payload:journal({'event':event,**payload}))
    pd=pd_derive_profile_difference({},lambda event,payload:journal({'event':event,**payload}))
    sb_assert(cg['passed'] and pd['passed'],'Signed difference tensor identity failed')
    sb_assert(receipts['controls']['authenticated_checks']['C3']['all_passed'] is True and
              receipts['controls']['authenticated_checks']['C5']['all_passed'] is True,'Source acceptance receipt missing')
    sb_assert(receipts['controls']['C3_limit_control_slice']['repair_weights']==
              receipts['controls']['C5_unique_system']['original_weights'],'Different weight ports')
    result={'schema':'original_compact_C3_C5_source_bridge','component':'compact_C3_C5_source_bridge',
        'passed':remap['proved'],'metadata_only':False,'source_only':True,'input_digest':source['input_digest'],
        'runtime_model_metadata':None,'original_source_evaluations':0,
        'source_packet_defined':remap['proved'],'source_ftc_transport_defined':True,'source_bridge_defined':remap['proved'],
        'incoming_normalization_identity_proved':True,'control_target_dependency_bridge_defined':True,
        'common_weights_and_selected_N_source_identity_proved':common['proved'],
        'C3_rebuilt_target_C3_namespace_identity_proved':remap['proved'],
        'normalized_target_gap_functions_defined':True,
        'control_gap_C3_equation_transport_proved':cg['passed'],
        'profile_difference_C3_jet_transport_proved':pd['passed'],
        'profile_C3_jet_identity_proved':False,
        'actual_C3_C5_incoming_source_identity_proved':False,
        'actual_C3_C5_control_source_identity_proved':False,
        'actual_C3_cap_to_C5_source_bound_proved':False,
        'compact_prefix_source_binding_proved':False,'compact_volume_cap_same_source_bound_proved':False,
        'actual_common_C1_control_ball_admitted':False,'signed_target_gaps_evaluated':False,
        'c3_c5_source_bridge_certificate':{'common_weights':common,'C3_namespace_remap':remap,
            'incoming_normalization':normalization,'signed_gap_graph':gap_graph,
            'control_gap_transport':cg,'profile_difference_transport':pd},
        'incoming_normalization':normalization,'control_target_dependency':gap_graph,
        'common_weights':common,'C3_namespace_remap':remap,'profile_difference_jet':pd,
        'source_binding':{'passed':remap['proved'],'common_admitted_source_domain_only':True,
            'port_contract':port_contract,
            'actual_gap_functions_defined':True,'gap_functions_used_as_zero':False,
            'fixed_point_residuals_are_cached_source_contracts':True,
            'contraction_stability_requires_same_C1_ball':True,
            'physical_or_global_cap_admitted':False},
        'coordinate_mapping':parent['coordinate_mapping'],
        'source_inspection_audit':{'original_source_evaluations':0,'old_source_consumer_calls':0,
            'old_proof_functions_called':0,'old_operator_calls':0,'old_source_normalizer_calls':0,
            'quadratures':0,'numeric_field_evaluations':0,'temporal_recurrence_calls':0,
            'new_source_namespace_comparisons':2,'new_incoming_history_order_identities':20,
            'new_control_gap_row_order_checks':20,'new_profile_difference_row_order_checks':24},
        'remaining_scope':[
            'Prove the five actual normalized target gap functions vanish or evaluate/bound their signed values, retaining amplitude gap A_N-A_C.',
            'The native all17 incoming histories still require connection to the rebuilt C3 source histories; shared weights/N do not establish that equality.',
            'Discharge the same-C1-ball and common-physical-bump contracts before quantitative control/profile cap transfer.',
            'Signed Evol, independent old reference offset and global physical/temporal/NS stages remain incomplete.']}
    for name in ('retained_cap_provenance','retained_construction_data','retained_integral_term_inventory'):
        result[name]=parent[name]
    for name in ('native_BG_exact_numeric_baseline_C3_identity_proved','relative_mass_C3_source_identity_proved',
                 'compact_prefix_geometry_logC_binding_proved','compact_prefix_fixed_mstar_enclosure_proved',
                 'compact_prefix_12channel_2volume_dominance_proved'):
        result[name]=parent[name]
    for name in source['inherited_false_flags']:result[name]=False
    journal({'event':'fresh_C3_C5_source_bridge_complete','actual_cap_admission':False})
    return result

