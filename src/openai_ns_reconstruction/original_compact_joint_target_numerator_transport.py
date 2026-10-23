"""Jointly normalized C3/C5 target-gap numerators, unentered v2 draft.

The source ports are symbolic ordinary-Z jets.  Their integer node identifiers
are labels, never numerical values.  For C=A_C, A=A_N=C+a and raw signed gaps
delta_T=T_N-T_C, this layer emits the five exact numerator/denominator pairs:

  eta_M  = (C*delta_m - a*S)/(C*A)
  eta_I  = (C*delta_h - a*H)/(C*A)
  eta_S  = (C^2*delta_e - a*(2*C+a)*E)/(C^2*A^2)
  eta_Cp = (C^2*delta_p - a*(2*C+a)*P)/(C^2*A^2)
  eta_D  = (C^2*delta_k - C^2*A*delta_m
            - a*(2*C+a)*K + a*C*A*S)/(mu*C^2*A^2).

The D numerator is kept joint; it is not formed by separately approximating
J and M.  Quotient jets are represented by triangular Leibniz relations, with
positive denominators as a conditional source contract.  A separate direct
ordinary-Z product expansion certifies the triangular reconstruction at every
order without evaluating a quotient or claiming a numerical bound.  No old
consumer, normalizer, proof, cap, midpoint, or source value is imported or
evaluated.
"""

from fractions import Fraction


tn_input_contract = {
    "schema": "compact_C3_C5_target_numerator_ports_v1",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "namespaces": ["C3", "C5"],
    "required_namespace_fields": {
        "archive_raw_sha256": "64-character SHA-256 string",
        "amplitude": "four integer source-node IDs ordered Z0..Z3",
        "histories": {
            "m": "four integer N-scaled raw history node IDs ordered Z0..Z3",
            "k": "four integer N-scaled raw history node IDs ordered Z0..Z3",
            "h": "four integer N-scaled raw history node IDs ordered Z0..Z3",
            "e": "four integer N-scaled raw history node IDs ordered Z0..Z3",
            "p": "four integer N-scaled raw history node IDs ordered Z0..Z3",
        },
    },
    "required_common_fields": {
        "mu": "one symbolic source parameter/node descriptor",
        "mu_Z_independent": True,
        "N": "one positive common source integer descriptor",
        "N_Z_independent": True,
    },
    "normalization": {
        "M": "scaled[m]/A",
        "I": "scaled[h]/A",
        "S": "scaled[e]/A^2",
        "Cp": "scaled[p]/A^2",
        "D": "(scaled[k]-A*scaled[m])/(mu*A^2)",
        "raw_gap": "delta_T=T_N-T_C with signed ordinary-Z jets",
    },
    "admission": "formal conditional source normalization only; no equality or cap admission",
}


tn_binomial = ((1,), (1, 1), (1, 2, 1), (1, 3, 3, 1))
tn_order_max = 3
tn_channels = ("M", "I", "S", "Cp", "D")
tn_history_names = ("m", "k", "h", "e", "p")


def tn_fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    raise TypeError("tn coefficients require int or Fraction")


def tn_poly_atom(name):
    return {(str(name),): Fraction(1, 1)}


def tn_poly_constant(value):
    value = tn_fraction(value)
    return {(): value} if value else {}


def tn_poly_add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0, 1)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def tn_poly_neg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def tn_poly_sub(left, right):
    return tn_poly_add(left, tn_poly_neg(right))


def tn_poly_mul(left, right):
    out = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial))
            out[monomial] = out.get(monomial, Fraction(0, 1)) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def tn_poly_scale(poly, coefficient):
    coefficient = tn_fraction(coefficient)
    return {
        monomial: coefficient * value
        for monomial, value in poly.items()
        if coefficient * value
    }


def tn_jet_zero():
    return ({}, {}, {}, {})


def tn_jet_atom(name):
    return tuple(tn_poly_atom("%s_Z%d" % (name, order)) for order in range(tn_order_max + 1))


def tn_jet_constant_atom(name):
    return (tn_poly_atom(name), {}, {}, {})


def tn_jet_add(*jets):
    if not jets:
        return tn_jet_zero()
    return tuple(tn_poly_add(*(jet[order] for jet in jets)) for order in range(tn_order_max + 1))


def tn_jet_neg(jet):
    return tuple(tn_poly_neg(value) for value in jet)


def tn_jet_sub(left, right):
    return tn_jet_add(left, tn_jet_neg(right))


def tn_jet_scale(jet, scalar):
    return tuple(tn_poly_mul(scalar, value) for value in jet)


def tn_jet_mul(left, right):
    """Ordinary-Z product rule with exact binomial coefficients through order 3."""
    return tuple(
        tn_poly_add(
            *[
                tn_poly_scale(
                    tn_poly_mul(left[left_order], right[order - left_order]),
                    tn_binomial[order][left_order],
                )
                for left_order in range(order + 1)
            ]
        )
        for order in range(tn_order_max + 1)
    )


def tn_jet_square(jet):
    return tn_jet_mul(jet, jet)


def tn_poly_equal(left, right):
    return not tn_poly_sub(left, right)


def tn_jet_equal(left, right):
    return all(tn_poly_equal(left[order], right[order]) for order in range(tn_order_max + 1))


def tn_encode_poly(poly):
    return [
        {
            "monomial": list(monomial),
            "coefficient": {"numerator": value.numerator, "denominator": value.denominator},
        }
        for monomial, value in sorted(poly.items())
    ]


def tn_encode_jet(jet):
    return [
        {"ordinary_Z_order": order, "terms": tn_encode_poly(jet[order])}
        for order in range(tn_order_max + 1)
    ]


def tn_source_node_id(value, label):
    if type(value) is not int or value < 0:
        raise ValueError("%s must be a nonnegative integer source-node ID" % label)
    return value


def tn_archive_sha(value, label):
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError("%s must be a 64-character source archive SHA-256" % label)
    return value


def tn_port_ids(namespace, ports):
    if not isinstance(ports, dict):
        raise TypeError("%s source ports must be a mapping" % namespace)
    archive = tn_archive_sha(ports.get("archive_raw_sha256"), namespace + ".archive_raw_sha256")
    amplitude = ports.get("amplitude")
    if not isinstance(amplitude, list) or len(amplitude) != tn_order_max + 1:
        raise ValueError("%s amplitude must contain four Z-ordered IDs" % namespace)
    amplitude = [tn_source_node_id(value, namespace + ".amplitude") for value in amplitude]
    histories = ports.get("histories")
    if not isinstance(histories, dict):
        raise ValueError("%s histories are missing" % namespace)
    normalized = {}
    for name in tn_history_names:
        values = histories.get(name)
        if not isinstance(values, list) or len(values) != tn_order_max + 1:
            raise ValueError("%s.histories.%s must contain four Z-ordered IDs" % (namespace, name))
        normalized[name] = [tn_source_node_id(value, namespace + ".histories." + name) for value in values]
    return {"archive_raw_sha256": archive, "amplitude": amplitude, "histories": normalized}


def tn_port_jet(namespace, quantity, ids):
    return tuple(
        tn_poly_atom("tn_%s_%s_node_%d" % (namespace, quantity, source_node))
        for source_node in ids
    )


def tn_source_jets(source_ports):
    if not isinstance(source_ports, dict):
        raise TypeError("source_ports must be a mapping")
    c3 = tn_port_ids("C3", source_ports.get("C3"))
    c5 = tn_port_ids("C5", source_ports.get("C5"))
    common = source_ports.get("common")
    if not isinstance(common, dict):
        raise ValueError("source_ports.common is missing")
    if common.get("mu_Z_independent") is not True or common.get("N_Z_independent") is not True:
        raise ValueError("mu and N Z-independence must be explicit")
    if "mu" not in common or "N" not in common:
        raise ValueError("common mu and N descriptors are required")
    source = {}
    for namespace, ports in (("C3", c3), ("C5", c5)):
        source[namespace] = {
            "archive_raw_sha256": ports["archive_raw_sha256"],
            "A": tn_port_jet(namespace, "A", ports["amplitude"]),
            "histories": {
                name: tn_port_jet(namespace, name, ports["histories"][name])
                for name in tn_history_names
            },
            "ids": ports,
        }
    source["common"] = {
        "mu": tn_jet_constant_atom("tn_mu"),
        "N": tn_jet_constant_atom("tn_N"),
        "mu_descriptor": common["mu"],
        "N_descriptor": common["N"],
    }
    return source


def tn_signed_gaps(source):
    c = source["C3"]
    n = source["C5"]
    return {
        "a": tn_jet_sub(n["A"], c["A"]),
        "histories": {
            name: tn_jet_sub(n["histories"][name], c["histories"][name])
            for name in tn_history_names
        },
    }


def tn_target_formulas(source, gaps):
    c = source["C3"]
    C = c["A"]
    A = source["C5"]["A"]
    a = gaps["a"]
    S = c["histories"]["m"]
    K = c["histories"]["k"]
    H = c["histories"]["h"]
    E = c["histories"]["e"]
    P = c["histories"]["p"]
    dm = gaps["histories"]["m"]
    dk = gaps["histories"]["k"]
    dh = gaps["histories"]["h"]
    de = gaps["histories"]["e"]
    dp = gaps["histories"]["p"]
    two_C_plus_a = tn_jet_add(tn_jet_scale(C, tn_poly_constant(2)), a)
    C2 = tn_jet_square(C)
    A2 = tn_jet_square(A)
    CA = tn_jet_mul(C, A)
    C2A2 = tn_jet_mul(C2, A2)
    mu = source["common"]["mu"]
    return {
        "M": {
            "numerator": tn_jet_sub(tn_jet_mul(C, dm), tn_jet_mul(a, S)),
            "denominator": CA,
            "raw_gap": "delta_m",
        },
        "I": {
            "numerator": tn_jet_sub(tn_jet_mul(C, dh), tn_jet_mul(a, H)),
            "denominator": CA,
            "raw_gap": "delta_h",
        },
        "S": {
            "numerator": tn_jet_sub(
                tn_jet_mul(C2, de),
                tn_jet_mul(tn_jet_mul(a, two_C_plus_a), E),
            ),
            "denominator": C2A2,
            "raw_gap": "delta_e",
        },
        "Cp": {
            "numerator": tn_jet_sub(
                tn_jet_mul(C2, dp),
                tn_jet_mul(tn_jet_mul(a, two_C_plus_a), P),
            ),
            "denominator": C2A2,
            "raw_gap": "delta_p",
        },
        "D": {
            "numerator": tn_jet_add(
                tn_jet_mul(C2, dk),
                tn_jet_neg(tn_jet_mul(tn_jet_mul(C2, A), dm)),
                tn_jet_neg(tn_jet_mul(tn_jet_mul(a, two_C_plus_a), K)),
                tn_jet_mul(tn_jet_mul(tn_jet_mul(a, C), A), S),
            ),
            "denominator": tn_jet_mul(mu, C2A2),
            "raw_gap": "joint_delta_k_and_delta_m",
            "joint_numerator": True,
        },
    }


def tn_direct_cross_multipliers(source, gaps):
    """Independent cross-products for checking the five fixed numerators."""
    c = source["C3"]
    n = source["C5"]
    C = c["A"]
    A = n["A"]
    C2 = tn_jet_square(C)
    A2 = tn_jet_square(A)
    mu = source["common"]["mu"]
    direct = {
        "M": tn_jet_sub(tn_jet_mul(C, n["histories"]["m"]), tn_jet_mul(A, c["histories"]["m"])),
        "I": tn_jet_sub(tn_jet_mul(C, n["histories"]["h"]), tn_jet_mul(A, c["histories"]["h"])),
        "S": tn_jet_sub(tn_jet_mul(C2, n["histories"]["e"]), tn_jet_mul(A2, c["histories"]["e"])),
        "Cp": tn_jet_sub(tn_jet_mul(C2, n["histories"]["p"]), tn_jet_mul(A2, c["histories"]["p"])),
        "D": tn_jet_sub(
            tn_jet_mul(C2, tn_jet_sub(n["histories"]["k"], tn_jet_mul(A, n["histories"]["m"]))),
            tn_jet_mul(A2, tn_jet_sub(c["histories"]["k"], tn_jet_mul(C, c["histories"]["m"]))),
        ),
    }
    # The D direct cross-product carries the common mu factor in its denominator.
    direct["D_denominator"] = tn_jet_mul(mu, tn_jet_mul(C2, A2))
    return direct


def tn_quotient_relations(channel, numerator, denominator):
    eta = tn_jet_atom("tn_eta_%s" % channel)
    relations = []
    for order in range(tn_order_max + 1):
        correction = tn_poly_add(
            *[
                tn_poly_scale(
                    tn_poly_mul(denominator[derivative_order], eta[order - derivative_order]),
                    tn_binomial[order][derivative_order],
                )
                for derivative_order in range(1, order + 1)
            ]
        )
        rhs = tn_poly_sub(numerator[order], correction)
        lhs = tn_poly_mul(denominator[0], eta[order])
        relations.append(
            {
                "channel": channel,
                "ordinary_Z_order": order,
                "denominator_zero_order_times_eta": tn_encode_poly(lhs),
                "rhs_numerator_minus_prior_denominator_terms": tn_encode_poly(rhs),
                "triangular": True,
                "conditional_on_source_positive_denominator": True,
            }
        )
    return relations


def tn_direct_product_order(left, right, order):
    """Independently expanded ordinary-Z product coefficient through order 3."""
    if order == 0:
        return tn_poly_mul(left[0], right[0])
    if order == 1:
        return tn_poly_add(
            tn_poly_mul(left[1], right[0]),
            tn_poly_mul(left[0], right[1]),
        )
    if order == 2:
        return tn_poly_add(
            tn_poly_mul(left[2], right[0]),
            tn_poly_scale(tn_poly_mul(left[1], right[1]), 2),
            tn_poly_mul(left[0], right[2]),
        )
    if order == 3:
        return tn_poly_add(
            tn_poly_mul(left[3], right[0]),
            tn_poly_scale(tn_poly_mul(left[2], right[1]), 3),
            tn_poly_scale(tn_poly_mul(left[1], right[2]), 3),
            tn_poly_mul(left[0], right[3]),
        )
    raise ValueError("ordinary-Z quotient certificate order must be 0..3")


def tn_quotient_product_certificate(channel, numerator, denominator):
    """Certify product reconstruction independently of tn_jet_mul.

    The eta jets remain formal atoms.  Each row verifies that the independent
    multinomial product coefficient equals the denominator-zero-order term
    plus the prior denominator terms.  A second link records that the full
    product-minus-numerator residual is exactly the emitted triangular
    relation residual; it is not asserted to vanish until that relation is
    bound by a caller under the positive-denominator contract.
    """
    eta = tn_jet_atom("tn_eta_%s" % channel)
    rows = []
    for order in range(tn_order_max + 1):
        prior = tn_poly_add(
            *[
                tn_poly_scale(
                    tn_poly_mul(denominator[derivative_order], eta[order - derivative_order]),
                    tn_binomial[order][derivative_order],
                )
                for derivative_order in range(1, order + 1)
            ]
        )
        denominator_zero_order_term = tn_poly_mul(denominator[0], eta[order])
        reconstructed = tn_poly_add(denominator_zero_order_term, prior)
        direct = tn_direct_product_order(denominator, eta, order)
        product_residual = tn_poly_sub(direct, reconstructed)
        full_product_minus_numerator = tn_poly_sub(direct, numerator[order])
        triangular_residual = tn_poly_sub(
            denominator_zero_order_term,
            tn_poly_sub(numerator[order], prior),
        )
        link_residual = tn_poly_sub(full_product_minus_numerator, triangular_residual)
        rows.append(
            {
                "channel": channel,
                "ordinary_Z_order": order,
                "independent_direct_product": tn_encode_poly(direct),
                "reconstructed_denominator_zero_plus_prior": tn_encode_poly(reconstructed),
                "product_reconstruction_residual": tn_encode_poly(product_residual),
                "full_product_minus_numerator": tn_encode_poly(full_product_minus_numerator),
                "triangular_relation_residual": tn_encode_poly(triangular_residual),
                "relation_link_residual": tn_encode_poly(link_residual),
                "passed": not product_residual and not link_residual,
                "conditional_on_source_positive_denominator": True,
                "quotient_jet_values_evaluated": False,
                "division_performed": False,
            }
        )
    return rows


def tn_check_report(left, right, channel):
    report = []
    for order in range(tn_order_max + 1):
        difference = tn_poly_sub(left[order], right[order])
        report.append(
            {
                "channel": channel,
                "ordinary_Z_order": order,
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


def tn_source_gap_manifest(source):
    rows = []
    for namespace in ("C3", "C5"):
        rows.append(
            {
                "namespace": namespace,
                "archive_raw_sha256": source[namespace]["archive_raw_sha256"],
                "amplitude_source_nodes_Z0_to_Z3": source[namespace]["ids"]["amplitude"],
                "history_source_nodes_Z0_to_Z3": source[namespace]["ids"]["histories"],
            }
        )
    return rows


def derive_target_numerator(source_ports, journal):
    """Emit signed target-gap rational numerators without evaluating ports."""
    source = tn_source_jets(source_ports)
    gaps = tn_signed_gaps(source)
    formulas = tn_target_formulas(source, gaps)
    direct = tn_direct_cross_multipliers(source, gaps)

    base_checks = []
    order_checks = []
    for channel in tn_channels:
        direct_value = direct[channel]
        formula_value = formulas[channel]["numerator"]
        checks = tn_check_report(direct_value, formula_value, channel)
        base_checks.append(
            {
                "channel": channel,
                "ordinary_Z_orders_checked": [0, 1, 2, 3],
                "passed": all(item["passed"] for item in checks),
            }
        )
        order_checks.extend(checks)

    C = source["C3"]["A"]
    A = source["C5"]["A"]
    a = gaps["a"]
    amplitude_checks = tn_check_report(A, tn_jet_add(C, a), "A=C+a")
    quotient = {
        channel: tn_quotient_relations(channel, formulas[channel]["numerator"], formulas[channel]["denominator"])
        for channel in tn_channels
    }
    quotient_product_certificate = {
        channel: tn_quotient_product_certificate(
            channel,
            formulas[channel]["numerator"],
            formulas[channel]["denominator"],
        )
        for channel in tn_channels
    }
    quotient_certificate_rows = [
        row
        for channel in tn_channels
        for row in quotient_product_certificate[channel]
    ]
    passed = (
        all(item["passed"] for item in order_checks)
        and all(item["passed"] for item in amplitude_checks)
        and all(item["passed"] for item in quotient_certificate_rows)
    )
    if journal is not None:
        journal(
            {
                "event": "joint_target_gap_numerators_constructed",
                "base_checks": len(base_checks),
                "order_checks": len(order_checks),
                "passed": passed,
                "source_equality_admitted": False,
                "numeric_evaluation": False,
            }
        )

    return {
        "schema": "compact_C3_C5_target_numerator_result_v2",
        "passed": passed,
        "conditional": True,
        "ordinary_Z_orders": [0, 1, 2, 3],
        "channels": list(tn_channels),
        "source_gap_manifest": tn_source_gap_manifest(source),
        "signed_gap_definitions": {
            "a": "A_N-A_C",
            "delta_m": "T_N_m-T_C_m",
            "delta_k": "T_N_k-T_C_k",
            "delta_h": "T_N_h-T_C_h",
            "delta_e": "T_N_e-T_C_e",
            "delta_p": "T_N_p-T_C_p",
        },
        "numerators_and_denominators": {
            channel: {
                "raw_gap_role": formulas[channel]["raw_gap"],
                "joint_numerator": formulas[channel].get("joint_numerator", False),
                "numerator": tn_encode_jet(formulas[channel]["numerator"]),
                "denominator": tn_encode_jet(formulas[channel]["denominator"]),
            }
            for channel in tn_channels
        },
        "quotient_jet_relations": quotient,
        "quotient_product_certificate": quotient_product_certificate,
        "quotient_product_certificate_check_count": len(quotient_certificate_rows),
        "quotient_product_certificate_passed": all(item["passed"] for item in quotient_certificate_rows),
        "cross_multiplication_checks": {
            "base_checks": base_checks,
            "order_checks": order_checks,
            "base_check_count": len(base_checks),
            "order_check_count": len(order_checks),
            "amplitude_decomposition_checks": amplitude_checks,
        },
        "conditional_denominator_contract": {
            "C": "positive A_C source amplitude",
            "A": "positive A_N source amplitude",
            "mu": "positive common source parameter",
            "N": "positive common source integer; no extra N division in these quotient numerators",
            "numeric_intervals_used": False,
        },
        "source_equality_admitted": False,
        "gap_zero_assumed": False,
        "cap_admitted": False,
        "old_math_calls": 0,
        "old_source_evaluations": 0,
        "new_numeric_evaluations": 0,
    }


def tn_input_contract_schema():
    return tn_input_contract



def jn_require(condition, message):
    if not condition:
        raise ValueError(message)


def jn_source_recipe_contract(spec):
    """Check ten base target recipes at declared source cuts, without evaluation."""
    literal = spec['normalization_slice']
    jn_require(literal['schema'] == 'compact_joint_target_normalization_literal_slice_v1',
               'Unexpected source normalization slice')
    proofs = []
    for namespace in ('C3', 'C5'):
        ports = spec['source_ports'][namespace]
        roles = literal['ordinary_Z_contract'][namespace]
        history_key = ('terminal_history_ports_Z0_to_Z3' if namespace == 'C3'
                       else 'terminal_scaled_history_ports_Z0_to_Z3')
        jn_require(ports['amplitude'] == roles['amplitude_ports_Z0_to_Z3'] and
                   ports['histories'] == roles[history_key], 'Source jet role changed')
        record_key = 'c3_repair_band' if namespace == 'C3' else 'true_compact_C5'
        raw_hash = literal['authenticated_inputs']['source_archives'][record_key]['raw_sha256']
        jn_require(ports['archive_raw_sha256'] == raw_hash, 'Source archive namespace changed')
        packet = literal[namespace]
        def row(entry):
            return entry.get('row', entry)
        def pair(entry, operation, ids):
            value = row(entry)
            jn_require(value['operation'] == operation and
                       (value['arguments'] == ids or value['arguments'] == list(reversed(ids))),
                       'Wrong two-operand source recipe')
        cuts = packet['cut_boundary_rows']
        for name in ('m', 'k', 'h', 'e', 'p'):
            jn_require(cuts[name]['id'] == ports['histories'][name][0], 'Wrong history source cut')
        amplitude = ports['amplitude'][0]
        jn_require(cuts['A']['id'] == amplitude and cuts['mu']['id'] == 17,
                   'Amplitude and mu source cuts changed')
        intermediate = packet['selected_intermediate_rows']
        squared = intermediate['A_squared']['id']
        mu_squared = intermediate['mu_times_A_squared']['id']
        joint = intermediate['joint_numerator']['id']
        negative = intermediate['joint_negative_A_m']['id']
        product = intermediate['joint_A_times_m']['id']
        pair(intermediate['A_squared'], 'product', [amplitude, amplitude])
        pair(intermediate['mu_times_A_squared'], 'product', [17, squared])
        pair(intermediate['joint_A_times_m'], 'product', [amplitude, ports['histories']['m'][0]])
        neg_row = row(intermediate['joint_negative_A_m'])
        jn_require(neg_row['operation'] == 'negative' and neg_row['argument'] == product,
                   'Wrong signed joint source term')
        pair(intermediate['joint_numerator'], 'sum', [ports['histories']['k'][0], negative])
        for target, history, denominator in (('M', 'm', amplitude), ('I', 'h', amplitude),
                                             ('S', 'e', squared), ('Cp', 'p', squared),
                                             ('D=(J-M)/mu', None, mu_squared)):
            entry = packet['selected_target_Z0_rows'][target]
            value = row(entry)
            numerator = joint if history is None else ports['histories'][history][0]
            jn_require(entry['id'] == roles['target_ports_Z0_to_Z3'][target][0] and
                       value['operation'] == 'positive_quotient' and
                       value['numerator'] == numerator and value['denominator'] == denominator and
                       bool(value['source_positive_certificate']), 'Wrong normalized target source recipe')
            proofs.append({'namespace': namespace, 'target': target, 'source_node': entry['id'],
                           'numerator': numerator, 'denominator': denominator,
                           'base_recipe_checked': True})
    return {'passed': True, 'base_source_recipe_count': len(proofs), 'base_source_recipes': proofs,
            'ordinary_Z_roles': literal['ordinary_Z_contract'],
            'common_mu_is_explicit_premise': True,
            'source_cuts_not_evaluated_or_identified_across_archives': True,
            'dictionary_key_order_has_mathematical_meaning': False,
            'source_positive_certificates_retained': True}


def derive(source, journal):
    spec = source['input']
    jn_require(spec['schema'] == 'compact_joint_target_numerator_transport_input_v1',
               'Unexpected fresh input')
    parent = source['parent_report']
    jn_require(source['parent_result']['passed'] is False and
               source['parent_result']['error'] == 'Native target row order changed',
               'Expected immutable failed source-bridge parent')
    contract = jn_source_recipe_contract(spec)
    journal({'event': 'fresh_base_normalization_source_recipes', 'count': 10})
    certificate = derive_target_numerator(spec['source_ports'], journal)
    jn_require(certificate['passed'] and certificate['cross_multiplication_checks']['order_check_count'] == 20 and
               certificate['quotient_product_certificate_check_count'] == 20 and
               certificate['quotient_product_certificate_passed'],
               'Fresh numerator identities failed')
    result = {
        'schema': 'original_compact_joint_target_numerator_transport',
        'component': 'compact_joint_target_numerator_transport',
        'passed': True, 'metadata_only': False, 'source_only': True,
        'input_digest': source['input_digest'], 'runtime_model_metadata': None,
        'original_source_evaluations': 0, 'source_packet_defined': True,
        'joint_target_numerator_source_defined': True,
        'five_joint_numerator_identities_proved': True,
        'ordinary_Z_quotient_transport_proved': True,
        'source_recipe_contract': contract, 'joint_target_numerator_certificate': certificate,
        'conditional_scope': ['common positive Z-independent mu',
            'positive source amplitudes on their common source domain',
            'the frozen ordinary-Z source roles, not numerical range endpoints'],
        'common_mu_source_identity_newly_proved': False,
        'actual_C3_C5_incoming_source_identity_proved': False,
        'actual_C3_C5_control_source_identity_proved': False,
        'actual_C3_C5_profile_identity_proved': False,
        'actual_C3_cap_to_C5_source_bound_proved': False,
        'compact_prefix_source_binding_proved': False,
        'compact_volume_cap_same_source_bound_proved': False,
        'signed_target_gaps_evaluated': False, 'actual_common_C1_control_ball_admitted': False,
        'source_binding': {'passed': True, 'conditional_target_difference_identity_only': True,
            'actual_target_gap_zero_proved': False, 'signed_gaps_retained': True,
            'finite_N_bounds_selected_as_values': False, 'physical_Pstar_lift_applied': False},
        'source_inspection_audit': {'original_source_evaluations': 0, 'old_source_consumer_calls': 0,
            'old_proof_functions_called': 0, 'old_operator_calls': 0, 'old_source_normalizer_calls': 0,
            'quadratures': 0, 'numeric_field_evaluations': 0, 'temporal_recurrence_calls': 0,
            'new_source_base_recipe_checks': 10, 'new_cross_product_order_checks': 20,
            'new_quotient_relation_checks': 20, 'new_amplitude_decomposition_checks': 4},
        'remaining_scope': ['Bind or bound the five raw N-scaled incoming history differences and the amplitude difference on the common source domain.',
            'Retain their coupled signed joint numerators when computing ordinary-Z target differences.',
            'The failed namespace/control/profile proof family remains immutable and unaccepted; actual C3/C5 equality and cap transfer remain open.',
            'Signed Evol, independent old reference offset, functional terminal closure and all later physical/time-recursion/NS layers remain incomplete.']}
    retained = source.get('retained_parent_data', source['cached_upstream_proof'])
    for name in ('retained_cap_provenance', 'retained_construction_data', 'retained_integral_term_inventory'):
        result[name] = retained[name]
    result['coordinate_mapping'] = spec['inherited_coordinate_mapping']
    for name in ('native_BG_exact_numeric_baseline_C3_identity_proved', 'relative_mass_C3_source_identity_proved',
                 'compact_prefix_geometry_logC_binding_proved', 'compact_prefix_fixed_mstar_enclosure_proved',
                 'compact_prefix_12channel_2volume_dominance_proved'):
        result[name] = parent[name]
    for name in source['inherited_false_flags']:
        result[name] = False
    journal({'event': 'fresh_joint_target_numerator_transport_complete', 'actual_gap_evaluation': False})
    return result
