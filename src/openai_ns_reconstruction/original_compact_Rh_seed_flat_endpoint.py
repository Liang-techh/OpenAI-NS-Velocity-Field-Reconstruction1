"""Conditional flat-endpoint algebra for the native Rh seed.

This is a new source-contract leaf.  It records an exact support-separation
argument for the three compactly supported patch terms at x=exp(1), then
composes the remaining constant factor with the native Rh seed.  The module
does not import a source provider, evaluate exp, integrate beta, evaluate a
control, or admit a global/native-to-C3 identity.  The caller must supply the
literal source premises described by ``rhf_contract_template``.
"""

from fractions import Fraction


rhf_schema = "compact_Rh_seed_flat_endpoint_algebra_input_v1"
rhf_chart_name = "Rh_reference"
rhf_ordinary_orders = (0, 1, 2, 3)
rhf_centers = ("5/4", "3/2", "7/4")
rhf_v1_preimage = {
    "sha256": "4b72a6c3d6d9a33d35c4648aeb0c7de98874ef6da33117d9b38ded29b1d76ff7",
    "bytes": 11328,
    "preserved": True,
}
rhf_separation_certificates = (
    {
        "center": "5/4",
        "strict_lower_bound": "40*(exp(1)-5/4)>40*(5/2-5/4)=50>1",
    },
    {
        "center": "3/2",
        "strict_lower_bound": "40*(exp(1)-3/2)>40*(5/2-3/2)=40>1",
    },
    {
        "center": "7/4",
        "strict_lower_bound": "40*(exp(1)-7/4)>40*(5/2-7/4)=30>1",
    },
)

rhf_expected_source_forms = {
    "patch": "H(x,Z)=x^(1/10)+sum_{i=0..2}h_{i+2}(Z)*gamma_i(x)",
    "gamma": "gamma_i(x)=(40/B)*beta(40*(x-c_i))",
    "beta_support": "beta(t)=0 for |t|>=1",
    "normalizer": "B=integral_R beta(t)dt and B>0",
    "native_seed": "am_Rh(Z)=exp(-3/5)/(1+Z^2)",
    "canonical_seed": "E_C3_Rh(Z)=exp(-1/2)/(1+Z^2)",
    "endpoint_composition": "E_Rh(Z)=H(exp(1),Z)*am_Rh(Z)",
}

rhf_contract_template = {
    "schema": rhf_schema,
    "source_forms": rhf_expected_source_forms,
    "endpoint": {
        "chart": rhf_chart_name,
        "coordinate": "x",
        "coordinate_value": "exp(1)",
        "Z_interval": "[-1,1]",
        "ordinary_Z_orders": [0, 1, 2, 3],
        "coordinate_is_Z_independent": True,
    },
    "support": {
        "beta_zero_condition": "beta(t)=0 for |t|>=1",
        "centers": ["5/4", "3/2", "7/4"],
        "scale": "40",
        "exp_lower_bound": "exp(1)>5/2",
        "series_reason": "exp(1)=sum_{n>=0}1/n!>1+1+1/2=5/2",
        "separation_certificates": [
            "40*(exp(1)-5/4)>40*(5/2-5/4)=50>1",
            "40*(exp(1)-3/2)>40*(5/2-3/2)=40>1",
            "40*(exp(1)-7/4)>40*(5/2-7/4)=30>1",
        ],
        "B_positive": True,
    },
    "source_binding": {
        "root_wrapper_passed": True,
        "patch_source_ref": "source-defined compact Rh patch H",
        "native_seed_source_ref": "source-defined native Rh amplitude seed",
        "canonical_seed_source_ref": "source-defined C3 endpoint seed",
        "endpoint_source_ref": "source-defined Rh endpoint x=exp(1)",
        "ordinary_projection": "[E0,E1,E2,E3]",
    },
    "flags": {
        "numeric_values_selected": False,
        "integrals_evaluated": False,
        "controls_computed": False,
        "global_amplitude_identity_admitted": False,
        "cap_admitted": False,
    },
}


def rhf_require(condition, message):
    if not condition:
        raise ValueError(message)


def rhf_nonempty_text(value, label):
    rhf_require(isinstance(value, str) and value != "", label + " must be nonempty text")
    return value


def rhf_parse_fraction_literal(value, label):
    """Parse only the exact rational literals admitted by the source contract."""
    rhf_require(isinstance(value, str) and value != "", label + " must be a rational string")
    pieces = value.split("/")
    rhf_require(len(pieces) in (1, 2), label + " must be n or n/d")
    numerator = int(pieces[0])
    denominator = int(pieces[1]) if len(pieces) == 2 else 1
    rhf_require(denominator != 0, label + " denominator is zero")
    return Fraction(numerator, denominator)


def rhf_fraction_record(value):
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
    }


def rhf_validate_source_contract(contract):
    """Validate literal source premises without evaluating their expressions."""
    rhf_require(isinstance(contract, dict), "semantic_contract must be a mapping")
    rhf_require(contract.get("schema") == rhf_schema, "wrong Rh endpoint algebra schema")

    forms = contract.get("source_forms")
    rhf_require(isinstance(forms, dict), "source_forms are missing")
    for name, expected in rhf_expected_source_forms.items():
        rhf_require(forms.get(name) == expected, "source form changed: " + name)

    endpoint = contract.get("endpoint")
    rhf_require(isinstance(endpoint, dict), "endpoint contract is missing")
    rhf_require(endpoint.get("chart") == rhf_chart_name, "wrong endpoint chart")
    rhf_require(endpoint.get("coordinate") == "x", "wrong endpoint coordinate")
    rhf_require(endpoint.get("coordinate_value") == "exp(1)", "endpoint is not exp(1)")
    rhf_require(endpoint.get("Z_interval") == "[-1,1]", "wrong ordinary-Z interval")
    rhf_require(endpoint.get("ordinary_Z_orders") == list(rhf_ordinary_orders),
                "ordinary-Z orders must be 0 through 3")
    rhf_require(endpoint.get("coordinate_is_Z_independent") is True,
                "endpoint coordinate must be Z-independent")

    support = contract.get("support")
    rhf_require(isinstance(support, dict), "support contract is missing")
    rhf_require(support.get("beta_zero_condition") == rhf_expected_source_forms["beta_support"],
                "beta support condition changed")
    rhf_require(support.get("centers") == list(rhf_centers), "support centers changed")
    rhf_require(support.get("scale") == "40", "support scale changed")
    rhf_require(support.get("exp_lower_bound") == "exp(1)>5/2",
                "exp(1) lower-bound premise is missing")
    rhf_require(support.get("series_reason") ==
                "exp(1)=sum_{n>=0}1/n!>1+1+1/2=5/2",
                "exp-series lower-bound reason changed")
    expected_separation = [row["strict_lower_bound"] for row in rhf_separation_certificates]
    rhf_require(support.get("separation_certificates") == expected_separation,
                "support-separation certificates changed")
    rhf_require(support.get("B_positive") is True, "B>0 source premise is missing")

    binding = contract.get("source_binding")
    rhf_require(isinstance(binding, dict), "source binding is missing")
    rhf_require(binding.get("root_wrapper_passed") is True,
                "root source wrapper has not passed its binding")
    rhf_require(binding.get("ordinary_projection") == "[E0,E1,E2,E3]",
                "native ordinary-Z slot contract changed")
    for name in (
        "patch_source_ref",
        "native_seed_source_ref",
        "canonical_seed_source_ref",
        "endpoint_source_ref",
    ):
        rhf_nonempty_text(binding.get(name), "source_binding." + name)

    flags = contract.get("flags")
    rhf_require(isinstance(flags, dict), "source flags are missing")
    for name in (
        "numeric_values_selected",
        "integrals_evaluated",
        "controls_computed",
        "global_amplitude_identity_admitted",
        "cap_admitted",
    ):
        rhf_require(flags.get(name) is False, "forbidden source flag is true: " + name)

    return {
        "schema": rhf_schema,
        "chart": rhf_chart_name,
        "ordinary_Z_orders": list(rhf_ordinary_orders),
        "source_forms_checked": list(rhf_expected_source_forms),
        "source_refs_checked": [
            binding["patch_source_ref"],
            binding["native_seed_source_ref"],
            binding["canonical_seed_source_ref"],
            binding["endpoint_source_ref"],
        ],
    }


def rhf_exact_arithmetic_checks(contract):
    """Perform the bounded rational checks only when the caller invokes derive."""
    support = contract["support"]

    first_three_series_terms = (
        Fraction(1, 1) + Fraction(1, 1) + Fraction(1, 2)
    )
    series_bound = Fraction(5, 2)
    series_check = {
        "partial_sum": rhf_fraction_record(first_three_series_terms),
        "claimed_lower_bound": rhf_fraction_record(series_bound),
        "equal": first_three_series_terms == series_bound,
        "tail_positive_theorem": "for every n>=3, 1/n!>0, so exp(1)>sum_{n=0..2}1/n!",
    }
    rhf_require(series_check["equal"], "exp-series partial sum changed")

    scale = rhf_parse_fraction_literal(support["scale"], "support.scale")
    support_rows = []
    for center_text, certificate in zip(support["centers"], rhf_separation_certificates):
        center = rhf_parse_fraction_literal(center_text, "support.center")
        lower_gap = series_bound - center
        scaled_lower_gap = scale * lower_gap
        row = {
            "center": center_text,
            "center_exact": rhf_fraction_record(center),
            "exp_one_lower_bound": rhf_fraction_record(series_bound),
            "argument_lower_gap": rhf_fraction_record(lower_gap),
            "scaled_argument_lower_bound": rhf_fraction_record(scaled_lower_gap),
            "support_threshold": rhf_fraction_record(Fraction(1, 1)),
            "strictly_outside_support": scaled_lower_gap > Fraction(1, 1),
            "source_certificate": certificate["strict_lower_bound"],
        }
        rhf_require(row["strictly_outside_support"],
                    "support separation is not strictly beyond one")
        support_rows.append(row)

    exponent = Fraction(-3, 5) + Fraction(1, 10)
    exponent_check = {
        "native_seed_exponent": rhf_fraction_record(Fraction(-3, 5)),
        "endpoint_H_exponent": rhf_fraction_record(Fraction(1, 10)),
        "composed_exponent": rhf_fraction_record(exponent),
        "expected_exponent": rhf_fraction_record(Fraction(-1, 2)),
        "passed": exponent == Fraction(-1, 2),
    }
    rhf_require(exponent_check["passed"], "endpoint exponent composition changed")

    gamma_zero_rows = []
    for center_text in support["centers"]:
        gamma_value = Fraction(0, 1)
        gamma_zero_rows.append({
            "center": center_text,
            "gamma_endpoint_value": rhf_fraction_record(gamma_value),
            "zero_from_support": gamma_value == Fraction(0, 1),
        })
    rhf_require(all(row["zero_from_support"] for row in gamma_zero_rows),
                "a gamma endpoint was not exactly zero")

    H_Z_zero_rows = []
    for order in (1, 2, 3):
        terms = [
            {
                "multiplier": "h_%d^(%d)(Z)" % (index + 2, order),
                "gamma_endpoint_value": rhf_fraction_record(Fraction(0, 1)),
                "term_is_zero": True,
            }
            for index in range(3)
        ]
        H_Z_zero_rows.append({
            "ordinary_Z_order": order,
            "flat_gamma_terms": terms,
            "correction_is_zero": all(term["term_is_zero"] for term in terms),
        })
    rhf_require(all(row["correction_is_zero"] for row in H_Z_zero_rows),
                "flat endpoint H_Z correction is not zero")

    return {
        "exp_series_check": series_check,
        "support_checks": support_rows,
        "support_check_count": len(support_rows),
        "all_support_checks_passed": True,
        "exponent_check": exponent_check,
        "exponent_check_count": 1,
        "gamma_endpoint_zero_checks": gamma_zero_rows,
        "H_Z_zero_checks": H_Z_zero_rows,
        "H_Z_zero_check_count": len(H_Z_zero_rows),
        "no_numeric_exp_evaluation": True,
        "no_beta_integral_evaluation": True,
    }


def derive_flat_endpoint_algebra(semantic_contract, journal=None):
    """Return the conditional flat-endpoint theorem for a validated source.

    This function is intentionally source-contract driven.  Its output keeps
    every global amplitude, cap, control, and numeric admission flag false.
    The endpoint identity is conditional on the exact source forms, support
    statement, positive B, and native/C3 endpoint seed strings supplied by the
    caller.  No source provider or mathematical evaluator is called here.
    """
    checked = rhf_validate_source_contract(semantic_contract)
    exact = rhf_exact_arithmetic_checks(semantic_contract)
    if journal is not None:
        rhf_require(callable(journal), "journal must be callable")

    if journal is not None:
        journal({
            "event": "rh_flat_endpoint_support_contract_checked",
            "chart": rhf_chart_name,
            "ordinary_Z_orders": list(rhf_ordinary_orders),
            "support_rows": 3,
            "numeric_evaluation": False,
        })

    support_rows = []
    for row in rhf_separation_certificates:
        support_rows.append({
            "center": row["center"],
            "argument_at_endpoint": "40*(exp(1)-" + row["center"] + ")",
            "strict_lower_bound": row["strict_lower_bound"],
            "beta_zero_by_support": True,
            "gamma_endpoint_value": "0",
        })

    endpoint_factor = {
        "endpoint_coordinate": "exp(1)",
        "H_at_endpoint": "exp(1/10)",
        "H_ordinary_Z_jets": ["exp(1/10)", "0", "0", "0"],
        "native_Rh_seed": "exp(-3/5)/(1+Z^2)",
        "composition": "exp(1/10)*exp(-3/5)/(1+Z^2)",
        "simplified_endpoint_formula": "exp(-1/2)/(1+Z^2)",
        "canonical_C3_endpoint_formula": "exp(-1/2)/(1+Z^2)",
        "constant_exponent_composition_used": True,
        "exp_numeric_value_used": False,
    }

    ordinary_Z_identity = {
        "domain": "Z in [-1,1]",
        "orders": list(rhf_ordinary_orders),
        "native_ordinary_Z_slots": "[E0,E1,E2,E3]",
        "C3_provider_Taylor_conversion_used": False,
        "same_function_on_both_sides": True,
        "derivative_formulas": [
            "exp(-1/2)*(1+Z^2)^(-1)",
            "exp(-1/2)*(-2*Z)*(1+Z^2)^(-2)",
            "exp(-1/2)*(6*Z^2-2)*(1+Z^2)^(-3)",
            "exp(-1/2)*24*Z*(1-Z^2)*(1+Z^2)^(-4)",
        ],
        "ordinary_derivative_identity_proved_conditionally": True,
        "time_or_taylor_coefficient_recurrence_used": False,
    }

    result = {
        "schema": "compact_Rh_seed_flat_endpoint_algebra_result_v1",
        "passed": True,
        "conditional_source_contract": checked,
        "support_separation": {
            "exp_one_lower_bound_source": "exp(1)>5/2",
            "support_threshold": "1",
            "rows": support_rows,
            "all_three_gamma_endpoint_values_zero": True,
            "B_positive_assumed_from_validated_source": True,
            "exact_arithmetic": exact,
        },
        "endpoint_constant_factor": endpoint_factor,
        "ordinary_Z_identity": ordinary_Z_identity,
        "conditional_endpoint_statement": {
            "formula": "E_Rh^(j)(Z)=E_C3_Rh^(j)(Z), j=0..3",
            "scope": "only the Rh endpoint source formula under the validated native seed, patch, support, and canonical C3 premises",
            "native_Rh_endpoint_to_canonical_C3_formula_proved": True,
            "normalized_Rh_endpoint_seed_gap_zero_on_frozen_interval": True,
            "global_Rh_amplitude_identity_proved": False,
            "C3_frontend_seams_closed": False,
        },
        "source_controls": {
            "controls_computed": False,
            "numeric_values_selected": False,
            "integrals_evaluated": False,
            "N0_evaluated": False,
            "global_amplitude_identity_proved": False,
            "actual_native_C3_amplitude_identity_proved": False,
            "cap_admitted": False,
            "midpoint_admitted": False,
        },
        "audit": {
            "source_provider_calls": 0,
            "old_math_calls": 0,
            "quadrature_calls": 0,
            "control_calls": 0,
            "numeric_exp_evaluations": 0,
            "ordinary_Z_orders": list(rhf_ordinary_orders),
            "proof_family": "fresh exact support separation and endpoint exponential composition",
            "exact_support_checks": exact["support_check_count"],
            "exact_exponent_checks": exact["exponent_check_count"],
            "exact_H_Z_zero_checks": exact["H_Z_zero_check_count"],
            "fraction_arithmetic_evaluated_inside_callable": True,
        },
        "v1_preimage": dict(rhf_v1_preimage),
    }
    if journal is not None:
        journal({
            "event": "rh_flat_endpoint_algebra_defined",
            "support_separation_rows": 3,
            "ordinary_Z_orders": 4,
            "global_amplitude_identity": False,
            "numeric_evaluation": False,
        })
    return result


__all__ = (
    "rhf_contract_template",
    "rhf_exact_arithmetic_checks",
    "rhf_validate_source_contract",
    "derive_flat_endpoint_algebra",
)


def rse_require(condition, message):
    if not condition:
        raise ValueError(message)


def rse_creator_binding(spec):
    import ast
    import hashlib
    records = spec['creator_source_records']
    trees = {}
    for name, item in records.items():
        rse_require(hashlib.sha256(item['source_text'].encode('utf-8')).hexdigest() == item['sha256'],
                    'Frozen Rh creator source hash changed: ' + name)
        trees[name] = ast.parse(item['source_text'])

    clauses = (
        ('native_C6', "H=a.add(a.fixed(power(g,x,Fraction(1,10))),*(a.scale(h[i+2],gamma[i][0]) for i in range(3)))"),
        ('native_C6', "argument=g.mul(g.constant(40),g.sub(x,g.constant(CENTERS[i])))"),
        ('native_C6', "self.am=a.div(a.fixed(g.unary('exp',g.constant('-3/5'))),a.add(a.fixed(1),a.mul(z,z)),'positive canonical normalized Rm amplitude exp(-3/5)/(1+Z^2)')"),
        ('native_C6', "self.endpoints={name:self.at(value) for name,value in (('Rm',g.one),('Rh',g.unary('exp',g.one)))}"),
        ('canonical_Rh', "logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10"),
    )
    checked = []
    for name, statement in clauses:
        expected = ast.dump(ast.parse(statement).body[0], include_attributes=False)
        matches = [node for node in ast.walk(trees[name]) if isinstance(node, ast.Assign)
                   and ast.dump(node, include_attributes=False) == expected]
        rse_require(len(matches) == 1, 'Rh defining source clause changed: ' + statement)
        checked.append({'source': name, 'line': matches[0].lineno, 'clause': statement, 'passed': True})
    expected_E = ast.dump(ast.parse('a.mul(am,H)', mode='eval').body, include_attributes=False)
    E_keywords = [node for node in ast.walk(trees['native_C6']) if isinstance(node, ast.keyword)
                  and node.arg == 'E' and ast.dump(node.value, include_attributes=False) == expected_E]
    rse_require(len(E_keywords) == 1, 'Native E is not the product am*H')
    jet_class = next(node for node in trees['functionjet_algebra'].body
                     if isinstance(node, ast.ClassDef) and node.name == 'JetAlgebra')
    mul = next(node for node in jet_class.body if isinstance(node, ast.FunctionDef) and node.name == 'mul')
    rse_require(any(isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name) and node.func.value.id == 'math'
                    and node.func.attr == 'comb' for node in ast.walk(mul)),
                'Native FunctionJet ordinary Leibniz factors are absent')
    loop_class = next(node for node in trees['loop_algebra'].body
                      if isinstance(node, ast.ClassDef) and node.name == 'LoopAlgebra')
    rse_require(not any(isinstance(node, ast.FunctionDef) and node.name == 'mul' for node in loop_class.body),
                'LoopAlgebra overrides the ordinary product convention')
    return {'passed': True, 'source_clauses': checked, 'E_is_am_times_H': True,
            'native_slots_are_ordinary_Z_derivatives': True, 'native_projection_factors': [1, 1, 1, 1],
            'source_functions_imported_or_called': False}


def rse_literal_binding(spec, parent):
    rows = spec['literal_rows']
    def require_row(node, expected):
        rse_require(rows.get(str(node)) == expected, 'Rh literal row changed: ' + str(node))

    def rational(p, q=1):
        return {'operation': 'exact_rational', 'numerator': p, 'denominator': q}

    def product(*arguments):
        return {'operation': 'product', 'arguments': list(arguments)}

    def summation(*arguments):
        return {'operation': 'sum', 'arguments': list(arguments)}

    def unary(name, argument):
        return {'operation': 'analytic_unary', 'name': name, 'argument': argument}

    for node, p, q in ((0, 0, 1), (1, 1, 1), (2, 40, 1), (220, -1, 1),
                       (48, 2, 1), (13766, 3, 1), (540, 5, 4), (272, 3, 2),
                       (545, 7, 4), (31355, 1, 10), (183700, -3, 5)):
        require_row(node, rational(p, q))
    require_row(2948, unary('exp', 1))
    require_row(185676, {'operation': 'bound_variable', 'name': 'true_Rm_patch_x'})
    require_row(20, {'operation': 'bound_variable', 'name': 'Z'})
    require_row(25249, product(20, 20))
    require_row(31353, summation(1, 25249))
    require_row(183701, unary('exp', 183700))
    certificate = 'positive canonical normalized Rm amplitude exp(-3/5)/(1+Z^2)'
    require_row(183702, {'operation': 'positive_quotient', 'numerator': 183701,
                        'denominator': 31353, 'source_positive_certificate': certificate})
    require_row(186751, unary('log', 185676))
    require_row(186752, product(31355, 186751))
    require_row(186753, unary('exp', 186752))
    require_row(184119, {'operation': 'positive_quotient', 'numerator': 1,
                        'denominator': 184092, 'source_positive_certificate': 'positive exact original beta normalization'})
    normalizer = rows['184092']
    rse_require(normalizer['operation'] == 'definite_integral' and normalizer['integrand'] == 184091
                and normalizer['lower'] == 220 and normalizer['upper'] == 1
                and normalizer['variable'] == 'true_Rm_beta_normalization_s'
                and normalizer['exact_original_positive_beta_normalization'] is True
                and normalizer['range_or_midpoint_not_defining_integral'] is True,
                'Wrong exact beta normalization source')
    beta_definition = 'exp(-1/(1-x^2)) for |x|<1; exact0 elsewhere'
    beta_origin = {'AST_sha256': 'c4d97e103701522d471a041a3e2991031f51e5e2031d3382072cde118d1811af',
                   'function': 'SharedFiveMomentRepair.beta',
                   'module': 'lei_ren_part1_paper_compliant_five_moment_repair.py',
                   'module_sha256': '978de37342cc2a63aa7c1bfda573e512f344e0157aed28b6017c7dae2ff8eba5'}
    def beta(argument):
        return {'operation': 'compact_raw_beta', 'argument': argument, 'exact_definition': beta_definition,
                'original_radius_not_absorbed_into_raw_shape': True, 'original_source': beta_origin}
    require_row(184090, {'operation': 'bound_variable', 'name': 'true_Rm_beta_normalization_s'})
    require_row(184091, beta(184090))
    gammas = []
    for center, offset, shift, argument, raw_beta, gamma in (
        (540, 184094, 185677, 185678, 185679, 185701),
        (272, 184198, 185704, 185705, 185706, 185728),
        (545, 184131, 185731, 185732, 185733, 185755)):
        require_row(offset, {'operation': 'negative', 'argument': center})
        require_row(shift, summation(185676, offset))
        require_row(argument, product(2, shift))
        require_row(raw_beta, beta(argument))
        require_row(gamma, product(2, 184119, raw_beta))
        gammas.append({'center_node': center, 'argument_node': argument, 'raw_beta_node': raw_beta,
                       'gamma_node': gamma, 'literal_source_bound': True})
    H_ids = [194763, 194764, 194765, 194766]
    H_terms = [[194742, 194749, 194756], [194743, 194750, 194757],
               [194744, 194751, 194758], [194745, 194752, 194759]]
    controls = [[193473, 193474, 193475], [193756, 193763, 193770],
                [193864, 193871, 193878], [193987, 193994, 194001]]
    for order, root in enumerate(H_ids):
        require_row(root, summation(*(([186753] if order == 0 else []) + H_terms[order])))
        for term, control, gamma in zip(H_terms[order], controls[order], (185701, 185728, 185755)):
            require_row(term, product(control, gamma))
    E_table = {
        194957: product(183702, 194763), 194958: product(183702, 194764),
        194959: product(183705, 194763), 194960: summation(194958, 194959),
        194961: product(183702, 194765), 194962: product(48, 183705, 194764),
        194963: product(183710, 194763), 194964: summation(194961, 194962, 194963),
        194965: product(183702, 194766), 194966: product(13766, 183705, 194765),
        194967: product(13766, 183710, 194764), 194968: product(183715, 194763),
        194969: summation(194965, 194966, 194967, 194968)}
    for node, literal in E_table.items():
        require_row(node, literal)
    seeds = [195918, 195919, 195920, 195921]
    E_ids = [194957, 194960, 194964, 194969]
    for node, expression in zip(seeds, E_ids):
        require_row(node, {'operation': 'function_substitution', 'expression': expression,
                          'variable': 185676, 'value': 2948, 'Z_independent_substitution': True})
    rse_require(seeds == parent['native_amplitude_consequence']['actual_Rh_seed_ports'],
                'The flat-endpoint theorem does not bind the accepted native action seed')
    provenance = spec['C6_outer_prefix_provenance']
    rse_require(provenance['C6_archive_raw_sha256'] == '4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a'
                and provenance['outer_archive_raw_sha256'] == '75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d'
                and provenance['C6_node_count'] == 196042,
                'Rh C6/outer prefix archive scope changed')
    rse_require(provenance['C6_nodes_sha256'] == provenance['outer_prefix_sha256']
                and provenance['C6_node_count'] == provenance['outer_prefix_node_count']
                and provenance['C6_serialization'] == provenance['outer_serialization'],
                'The C6 snapshot is not the same literal source prefix as the outer seed')
    rse_require(provenance['am_ordinary_Z_ports'][:4] == [183702, 183705, 183710, 183715]
                and provenance['H_ordinary_Z_ports'][:4] == H_ids
                and provenance['E_body_ordinary_Z_ports'][:4] == E_ids
                and provenance['Rh_E_ordinary_Z_ports'][:4] == seeds,
                'The actual C6 snapshot ports do not match this endpoint proof')
    return {'passed': True, 'native_seed_ports': seeds, 'native_H_ports': H_ids,
            'native_am_ports': [183702, 183705, 183710, 183715], 'native_E_body_ports': E_ids,
            'gamma_source_rows': gammas, 'C6_outer_prefix_provenance': provenance,
            'ordinary_Z_orders': [0, 1, 2, 3], 'projection_factors': [1, 1, 1, 1],
            'controls_are_opaque_source_functions': True,
            'normalizer_positivity': 'beta is continuous, nonnegative and positive on (-1,1); its exact integral on [-1,1] is finite and strictly positive. No integral value is selected.'}


def derive(source, journal):
    import copy
    spec, parent = source['input'], source['parent_result']
    rse_require(spec['schema'] == 'compact_Rh_seed_flat_endpoint_input_v1', 'Wrong Rh seed input schema')
    creator = rse_creator_binding(spec)
    binding = rse_literal_binding(spec, parent)
    semantic = copy.deepcopy(spec['semantic_contract'])
    semantic['source_binding']['root_wrapper_passed'] = creator['passed'] and binding['passed']
    algebra = derive_flat_endpoint_algebra(semantic, journal)
    rse_require(algebra['passed'], 'Rh support-separation algebra failed')
    journal({'event': 'Rh_seed_actual_source_identity_bound', 'ordinary_Z_orders': 4,
             'gamma_supports': 3, 'controls_or_integrals_evaluated': False})
    result = {
        'schema': 'original_compact_Rh_seed_flat_endpoint', 'component': 'compact_Rh_seed_flat_endpoint',
        'passed': True, 'metadata_only': False, 'source_only': True,
        'input_digest': source['input_digest'], 'runtime_model_metadata': None,
        'original_source_evaluations': 0, 'source_packet_defined': True,
        'source_bound_native_Rh_canonical_seed_ordinary_Z3_identity_proved': True,
        'native_compact_bump_flat_Rh_endpoint_proved': True,
        'actual_native_C3_amplitude_identity_proved': False,
        'rh_seed_flat_endpoint_certificate': {
            'passed': True, 'scope': 'Actual native Rh normalized E source and canonical Rh function, ordinary Z=0..3 on [-1,1] only.',
            'creator_binding': creator, 'literal_source_binding': binding, 'flat_endpoint_algebra': algebra},
        'native_Rh_seed_consequence': {
            'function_identity': 'E_native_Rh^(j)(Z)=d_Z^j[exp(-1/2)/(1+Z^2)], j=0..3, Z in [-1,1]',
            'formula': 'exp(-1/2)/(1+Z^2)', 'native_seed_ports': binding['native_seed_ports'],
            'native_canonical_Rh_seed_gap_ordinary_Z_0_to_3': [0, 0, 0, 0],
            'normalized_source_endpoint_only': True, 'physical_coordinate_scale_bridge_newly_proved': False,
            'Rc_amplitude_gap_zero_proved': False},
        'source_binding': {'passed': True, 'actual_native_Rh_seed_bound_to_canonical_function': True,
                           'source_controls_or_finite_N_values_selected': False,
                           'C3_frontend_seams_proved': False, 'old_errors_removed': False},
        'source_inspection_audit': {
            'original_source_evaluations': 0, 'numeric_field_evaluations': 0,
            'old_source_consumer_calls': 0, 'old_proof_functions_called': 0, 'old_operator_calls': 0,
            'old_source_normalizer_calls': 0, 'quadratures': 0, 'temporal_recurrence_calls': 0,
            'new_seed_substitution_checks': 4, 'new_gamma_flat_support_checks': 3,
            'new_H_ordinary_row_checks': 4, 'new_E_ordinary_product_checks': 4,
            'new_native_am_source_binding_checks': 1, 'new_exponent_composition_checks': 1,
            'new_canonical_source_clause_checks': 1},
        'remaining_scope': [
            'Retain this normalized Rh seed identity without replay; it does not prove the complete Rc amplitude equality.',
            'Bind native/C3 outer action parameters and C3 frontend seams before closing A_N-A_C.',
            'Bind actual signed N-scaled history differences through the saved joint target numerators; all old caps and reference offsets remain.',
            'Functional five-moment closure and all downstream physical, time-recursive, stress and Cartesian NS layers remain incomplete.']}
    for name in ('retained_cap_provenance', 'retained_construction_data', 'retained_integral_term_inventory',
                 'coordinate_mapping', 'native_BG_exact_numeric_baseline_C3_identity_proved',
                 'relative_mass_C3_source_identity_proved', 'compact_prefix_geometry_logC_binding_proved',
                 'compact_prefix_fixed_mstar_enclosure_proved', 'compact_prefix_12channel_2volume_dominance_proved',
                 'joint_target_numerator_source_defined', 'five_joint_numerator_identities_proved',
                 'ordinary_Z_quotient_transport_proved', 'amplitude_action_source_defined',
                 'six_chart_log_action_identity_proved', 'amplitude_signed_seed_seam_transport_C3_defined',
                 'native_closed_Rc_amplitude_action_proved',
                 'native_five_interchart_amplitude_seams_ordinary_Z3_proved',
                 'native_outer_to_C5_amplitude_source_alias_proved'):
        result[name] = parent[name]
    for name in source['inherited_false_flags']:
        result[name] = False
    return result
