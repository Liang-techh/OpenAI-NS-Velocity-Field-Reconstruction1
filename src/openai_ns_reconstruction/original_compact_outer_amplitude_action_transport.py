"""Unentered leading-amplitude seam and propagator ledger draft.

This is a fresh symbolic source-contract layer for six ordered outer charts.
Native and C3 attenuation actions are kept separate, so no equality of their
sigma, mu, exp(40), or endpoint data is inferred.  Each of the five
interchart seams is retained as an additive signed amplitude defect and is
propagated by the suffix action on its own side.  The common-factor relation
and the variation-of-constants specialization are conditional source
statements; this module does not admit equality, a cap, a midpoint, or a
numeric source value.
"""

from fractions import Fraction


las_order_max = 3
las_sides = ("native", "C3")
las_action_parameter_names = ("I_slope", "I_trans", "mu", "exp40")
las_chart_names = (
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)


las_chart_specs = (
    {
        "name": "Rh_reference",
        "interval": ("-5", "0"),
        "rate_formula": "1/10",
        "action_formula": "5*(1/10)=1/2",
        "action_kind": "fixed_rate",
    },
    {
        "name": "O2_slope",
        "interval": ("0", "1"),
        "rate_formula": "1/10-(3/5)*sigmaS(x)",
        "action_formula": "1/10-(3/5)*I_slope",
        "action_kind": "sigma_slope",
    },
    {
        "name": "O2_axial",
        "interval": ("0", "1"),
        "rate_formula": "-20*exp(40*x)",
        "action_formula": "-20*(exp(40)-1)/40=(1-exp40)/2",
        "action_kind": "axial_exp40",
    },
    {
        "name": "O2_buffer",
        "interval": ("0", "11"),
        "rate_formula": "-1/2",
        "action_formula": "11*(-1/2)=-11/2",
        "action_kind": "fixed_buffer",
    },
    {
        "name": "O3_transition",
        "interval": ("0", "1"),
        "rate_formula": "-1/2-mu*sigmaT(x)",
        "action_formula": "-1/2-mu*I_trans",
        "action_kind": "sigma_transition",
    },
    {
        "name": "O3_power",
        "interval": ("0", "2"),
        "rate_formula": "-(1/2+mu)",
        "action_formula": "-1-2*mu",
        "action_kind": "fixed_power",
    },
)


las_input_contract = {
    "schema": "compact_leading_amplitude_seam_ports_v2",
    "sides": ["native", "C3"],
    "ordinary_Z_orders": [0, 1, 2, 3],
    "side_fields": {
        "archive_raw_sha256": "64-character raw source archive SHA-256",
        "seed": {
            "ordinary_Z_nodes": "four ordered integer IDs or authenticated query descriptors for E_Rh at Z0..Z3"
        },
        "action_parameters": {
            "I_slope": "source descriptor for integral sigmaS, Z-independent on this side",
            "I_trans": "source descriptor for integral sigmaT, Z-independent on this side",
            "mu": "source descriptor for the O3 transition/power parameter",
            "exp40": "source descriptor for exp(40) in the axial action",
        },
        "charts": {
            "each": {
                "action_ref": "source equation/row reference",
                "inlet": "source_ref plus four ordered integer IDs or query descriptors",
                "outlet": "source_ref plus four ordered integer IDs or query descriptors",
            }
        },
        "adjacent_seam_refs": "five source references, one per next-inlet minus previous-outlet seam",
    },
    "premises": {
        "action_parameters_Z_independent": "must be explicitly supplied separately for each side",
        "common_factor": "unproved conditional specialization only",
        "adjacent_endpoint_equalities": "unproved; signed defects remain in the output",
    },
}


def las_fraction(numerator, denominator=1):
    return Fraction(numerator, denominator)


def las_affine_add(*terms):
    out = {}
    for term in terms:
        for name, coefficient in term.items():
            out[name] = out.get(name, Fraction(0, 1)) + coefficient
            if not out[name]:
                del out[name]
    return out


def las_affine_neg(term):
    return {name: -coefficient for name, coefficient in term.items()}


def las_affine_sub(left, right):
    return las_affine_add(left, las_affine_neg(right))


def las_affine_scale(term, scalar):
    scalar = scalar if isinstance(scalar, Fraction) else Fraction(scalar, 1)
    return {
        name: scalar * coefficient
        for name, coefficient in term.items()
        if scalar * coefficient
    }


def las_affine_constant(numerator, denominator=1):
    value = Fraction(numerator, denominator)
    return {"1": value} if value else {}


def las_affine_atom(name):
    return {str(name): Fraction(1, 1)}


def las_encode_affine(term):
    return [
        {
            "atom": name,
            "coefficient": {
                "numerator": coefficient.numerator,
                "denominator": coefficient.denominator,
            },
        }
        for name, coefficient in sorted(term.items())
    ]


def las_expr_atom(name):
    return {"op": "atom", "name": str(name)}


def las_expr_add(*terms):
    return {"op": "add", "terms": list(terms)}


def las_expr_sub(left, right):
    return {"op": "sub", "left": left, "right": right}


def las_expr_mul(left, right):
    return {"op": "mul", "left": left, "right": right}


def las_expr_exp(argument):
    return {"op": "exp", "argument": argument}


def las_expr_affine(term):
    return {"op": "affine", "terms": las_encode_affine(term)}


def las_source_ref(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError("%s must be a nonempty source reference" % label)
    return value


def las_node_id(value, label):
    if type(value) is not int or value < 0:
        raise ValueError("%s must be a nonnegative integer source-node ID" % label)
    return value


def las_jet_port(value, label, order):
    if type(value) is int:
        return {
            "kind": "node_id",
            "node_id": las_node_id(value, label),
            "order": order,
        }
    if not isinstance(value, dict):
        raise TypeError("%s must be an integer node ID or query descriptor" % label)
    archive = value.get("archivehash", value.get("archive_raw_sha256"))
    archive = las_sha(archive, label + ".archivehash")
    function_node = value.get("functionnode", value.get("function_node"))
    function_node = las_source_ref(function_node, label + ".functionnode")
    if "querycoordinate" in value:
        query_coordinate = value["querycoordinate"]
    else:
        query_coordinate = value.get("query_coordinate")
    if query_coordinate is None:
        raise ValueError(label + ".querycoordinate is missing")
    descriptor_order = value.get("order")
    if type(descriptor_order) is not int or descriptor_order != order:
        raise ValueError("%s.order must equal its ordinary-Z position" % label)
    return {
        "kind": "query_descriptor",
        "archivehash": archive,
        "functionnode": function_node,
        "querycoordinate": query_coordinate,
        "order": order,
    }


def las_port_token(port):
    if port["kind"] == "node_id":
        return "node_%d" % port["node_id"]
    function_node = str(port["functionnode"]).replace(" ", "_")
    query_coordinate = str(port["querycoordinate"]).replace(" ", "_")
    return "query_%s_%s_order_%d" % (function_node, query_coordinate, port["order"])


def las_sha(value, label):
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError("%s must be a 64-character raw archive SHA-256" % label)
    return value


def las_endpoint_port(value, label):
    if not isinstance(value, dict):
        raise TypeError("%s endpoint port must be a mapping" % label)
    source_ref = las_source_ref(value.get("source_ref"), label + ".source_ref")
    raw_ports = value.get("ordinary_Z_nodes", value.get("ordinary_Z_ports"))
    if not isinstance(raw_ports, list) or len(raw_ports) != las_order_max + 1:
        raise ValueError("%s must contain four ordinary-Z ports" % label)
    ports = [
        las_jet_port(port, label + ".ordinary_Z_ports[%d]" % order, order)
        for order, port in enumerate(raw_ports)
    ]
    return {"source_ref": source_ref, "ordinary_Z_ports": ports}


def las_side_ports(side, value):
    if not isinstance(value, dict):
        raise TypeError("%s side contract must be a mapping" % side)
    archive = las_sha(value.get("archive_raw_sha256"), side + ".archive_raw_sha256")
    seed = value.get("seed")
    if not isinstance(seed, dict):
        raise ValueError("%s.seed is missing" % side)
    raw_seed_ports = seed.get("ordinary_Z_nodes", seed.get("ordinary_Z_ports"))
    if not isinstance(raw_seed_ports, list) or len(raw_seed_ports) != las_order_max + 1:
        raise ValueError("%s.seed must contain four ordinary-Z ports" % side)
    seed_ports = [
        las_jet_port(port, side + ".seed.ordinary_Z_ports[%d]" % order, order)
        for order, port in enumerate(raw_seed_ports)
    ]

    parameters = value.get("action_parameters")
    if not isinstance(parameters, dict):
        raise ValueError("%s.action_parameters is missing" % side)
    normalized_parameters = {}
    for name in las_action_parameter_names:
        descriptor = parameters.get(name)
        if not isinstance(descriptor, dict):
            raise ValueError("%s.action_parameters.%s is missing" % (side, name))
        if descriptor.get("ordinary_Z_independent") is not True:
            raise ValueError(
                "%s.action_parameters.%s must explicitly be ordinary-Z independent"
                % (side, name)
            )
        normalized_parameters[name] = {
            "source_ref": las_source_ref(
                descriptor.get("source_ref"),
                side + ".action_parameters." + name + ".source_ref",
            ),
            "ordinary_Z_independent": True,
        }

    charts = value.get("charts")
    if not isinstance(charts, dict):
        raise ValueError("%s.charts is missing" % side)
    normalized_charts = {}
    for chart in las_chart_names:
        row = charts.get(chart)
        if not isinstance(row, dict):
            raise ValueError("%s.charts.%s is missing" % (side, chart))
        normalized_charts[chart] = {
            "action_ref": las_source_ref(
                row.get("action_ref"),
                side + ".charts." + chart + ".action_ref",
            ),
            "inlet": las_endpoint_port(row.get("inlet"), side + ".charts." + chart + ".inlet"),
            "outlet": las_endpoint_port(row.get("outlet"), side + ".charts." + chart + ".outlet"),
        }

    adjacent = value.get("adjacent_seam_refs")
    if not isinstance(adjacent, list) or len(adjacent) != len(las_chart_names) - 1:
        raise ValueError(side + ".adjacent_seam_refs must contain five references")
    adjacent = [
        las_source_ref(ref, side + ".adjacent_seam_refs[%d]" % index)
        for index, ref in enumerate(adjacent)
    ]
    return {
        "archive_raw_sha256": archive,
        "seed": {"ordinary_Z_ports": seed_ports, "source_ref": las_source_ref(seed.get("source_ref"), side + ".seed.source_ref")},
        "action_parameters": normalized_parameters,
        "charts": normalized_charts,
        "adjacent_seam_refs": adjacent,
    }


def las_action_ledger(side, ports):
    parameters = {
        name: las_affine_atom("las_%s_%s" % (side, name))
        for name in las_action_parameter_names
    }
    actions = []
    for spec in las_chart_specs:
        kind = spec["action_kind"]
        if kind == "fixed_rate":
            action = las_affine_constant(1, 2)
        elif kind == "sigma_slope":
            action = las_affine_add(
                las_affine_constant(1, 10),
                las_affine_scale(parameters["I_slope"], Fraction(-3, 5)),
            )
        elif kind == "axial_exp40":
            action = las_affine_add(
                las_affine_constant(1, 2),
                las_affine_scale(parameters["exp40"], Fraction(-1, 2)),
            )
        elif kind == "fixed_buffer":
            action = las_affine_constant(-11, 2)
        elif kind == "sigma_transition":
            action = las_affine_add(
                las_affine_constant(-1, 2),
                las_affine_scale(
                    las_affine_mul_symbol(parameters["mu"], parameters["I_trans"]),
                    Fraction(-1, 1),
                ),
            )
        elif kind == "fixed_power":
            action = las_affine_add(
                las_affine_constant(-1, 1),
                las_affine_scale(parameters["mu"], Fraction(-2, 1)),
            )
        else:
            raise ValueError("unsupported chart action kind")
        actions.append(
            {
                "name": spec["name"],
                "interval": list(spec["interval"]),
                "rate_formula": spec["rate_formula"],
                "action_formula": spec["action_formula"],
                "action_ref": ports["charts"][spec["name"]]["action_ref"],
                "action": las_encode_affine(action),
                "action_affine": action,
            }
        )
    return {"parameters": parameters, "actions": actions}


def las_affine_mul_symbol(left, right):
    """Represent a product of symbolic affine atoms as one exact ledger atom."""
    if len(left) != 1 or len(right) != 1:
        raise ValueError("action parameters must be atomic for the affine ledger")
    left_name = next(iter(left))
    right_name = next(iter(right))
    return {"(%s)*(%s)" % (left_name, right_name): Fraction(1, 1)}


def las_expected_total_action(parameters):
    return las_affine_add(
        las_affine_constant(-59, 10),
        las_affine_scale(parameters["exp40"], Fraction(-1, 2)),
        las_affine_scale(parameters["I_slope"], Fraction(-3, 5)),
        las_affine_scale(
            las_affine_mul_symbol(parameters["mu"], parameters["I_trans"]),
            Fraction(-1, 1),
        ),
        las_affine_scale(parameters["mu"], Fraction(-2, 1)),
    )


def las_endpoint_atom(side, chart, role, order, port):
    return las_expr_atom(
        "las_%s_%s_%s_Z%d_%s" % (side, chart, role, order, las_port_token(port))
    )


def las_seed_jets(side, ports):
    return [
        las_expr_atom("las_%s_seed_Z%d_%s" % (side, order, las_port_token(port)))
        for order, port in enumerate(ports["seed"]["ordinary_Z_ports"])
    ]


def las_eta_jets(side, ports, seam_index):
    previous_chart = las_chart_names[seam_index]
    next_chart = las_chart_names[seam_index + 1]
    previous = ports["charts"][previous_chart]["outlet"]["ordinary_Z_ports"]
    following = ports["charts"][next_chart]["inlet"]["ordinary_Z_ports"]
    return [
        las_expr_sub(
            las_endpoint_atom(side, next_chart, "inlet", order, following[order]),
            las_endpoint_atom(side, previous_chart, "outlet", order, previous[order]),
        )
        for order in range(las_order_max + 1)
    ]


def las_side_ledger(side, ports):
    action_ledger = las_action_ledger(side, ports)
    actions = action_ledger["actions"]
    action_affines = [row["action_affine"] for row in actions]
    total_action = las_affine_add(*action_affines)
    expected_total_action = las_expected_total_action(action_ledger["parameters"])
    suffixes = []
    eta_rows = []
    for seam_index in range(len(las_chart_names) - 1):
        previous_chart = las_chart_names[seam_index]
        next_chart = las_chart_names[seam_index + 1]
        suffix_action = las_affine_add(*action_affines[seam_index + 1 :])
        suffixes.append(
            {
                "seam_index": seam_index,
                "after_chart": las_chart_names[seam_index],
                "before_chart": las_chart_names[seam_index + 1],
                "action": las_encode_affine(suffix_action),
                "action_affine": suffix_action,
                "factor": las_expr_exp(las_expr_affine(suffix_action)),
            }
        )
        eta_rows.append(
            {
                "seam_index": seam_index,
                "source_ref": ports["adjacent_seam_refs"][seam_index],
                "previous_outlet_ref": ports["charts"][previous_chart]["outlet"]["source_ref"],
                "next_inlet_ref": ports["charts"][next_chart]["inlet"]["source_ref"],
                "ordinary_Z_jets": las_eta_jets(side, ports, seam_index),
            }
        )
    seed = las_seed_jets(side, ports)
    factor = las_expr_exp(las_expr_affine(total_action))
    terminal_jets = []
    for order in range(las_order_max + 1):
        seam_terms = [
            las_expr_mul(suffixes[index]["factor"], eta_rows[index]["ordinary_Z_jets"][order])
            for index in range(len(suffixes))
        ]
        terminal_jets.append(
            las_expr_add(las_expr_mul(factor, seed[order]), *seam_terms)
        )
    telescope = []
    cumulative = {}
    for index, row in enumerate(actions):
        cumulative = las_affine_add(cumulative, row["action_affine"])
        telescope.append(
            {
                "chart": row["name"],
                "action_ref": row["action_ref"],
                "cumulative_log_action": las_encode_affine(cumulative),
            }
        )
    return {
        "side": side,
        "archive_raw_sha256": ports["archive_raw_sha256"],
        "seed": seed,
        "action_parameters": {
            name: {
                "source_ref": ports["action_parameters"][name]["source_ref"],
                "ordinary_Z_independent": ports["action_parameters"][name]["ordinary_Z_independent"],
            }
            for name in las_action_parameter_names
        },
        "action_ledger": [
            {
                key: value
                for key, value in row.items()
                if key != "action_affine"
            }
            for row in actions
        ],
        "chart_endpoint_manifest": [
            {
                "chart": chart,
                "inlet_ref": ports["charts"][chart]["inlet"]["source_ref"],
                "inlet_ordinary_Z_ports": ports["charts"][chart]["inlet"]["ordinary_Z_ports"],
                "outlet_ref": ports["charts"][chart]["outlet"]["source_ref"],
                "outlet_ordinary_Z_ports": ports["charts"][chart]["outlet"]["ordinary_Z_ports"],
            }
            for chart in las_chart_names
        ],
        "log_action_total": las_encode_affine(total_action),
        "log_action_total_affine": total_action,
        "closed_total_formula": las_encode_affine(expected_total_action),
        "closed_total_formula_string": "-59/10-exp40/2-(3/5)*I_slope-mu*I_trans-2*mu",
        "closed_total_formula_passed": total_action == expected_total_action,
        "factor": factor,
        "factor_Z_independence": {
            "ordinary_Z_orders": list(range(las_order_max + 1)),
            "F_Z_orders": [factor, 0, 0, 0],
            "interpretation": "[F, F_Z, F_ZZ, F_ZZZ] under the side parameter contract",
            "conditional_on_action_parameter_contract": True,
        },
        "suffix_ledger": [
            {
                key: value
                for key, value in row.items()
                if key not in ("action_affine",)
            }
            for row in suffixes
        ],
        "suffix_internal": suffixes,
        "signed_seam_defects": eta_rows,
        "ordinary_Z_terminal_expression": terminal_jets,
        "normalized_log_telescope": telescope,
    }


def las_book_add(*maps):
    out = {}
    for mapping in maps:
        for basis, coefficient in mapping.items():
            out[basis] = out.get(basis, Fraction(0, 1)) + coefficient
            if not out[basis]:
                del out[basis]
    return out


def las_book_scale(mapping, scalar):
    scalar = scalar if isinstance(scalar, Fraction) else Fraction(scalar, 1)
    return {
        basis: scalar * coefficient
        for basis, coefficient in mapping.items()
        if scalar * coefficient
    }


def las_book_append_factor(mapping, factor):
    return {
        basis + (factor,): coefficient
        for basis, coefficient in mapping.items()
    }


def las_encode_book(mapping):
    return [
        {
            "basis": list(basis),
            "coefficient": {
                "numerator": coefficient.numerator,
                "denominator": coefficient.denominator,
            },
        }
        for basis, coefficient in sorted(mapping.items())
    ]


def las_propagation_bookkeeping(side, order):
    """Independent six-step recurrence over formal factor/seam basis terms."""
    state = {("seed", side, order): Fraction(1, 1)}
    for chart_index in range(len(las_chart_names)):
        state = las_book_append_factor(state, "P_%s_%d" % (side, chart_index))
        if chart_index < len(las_chart_names) - 1:
            state[("seam", side, chart_index, order)] = Fraction(1, 1)
    return state


def las_closed_bookkeeping(side, order):
    """Direct seed-plus-suffix coefficient map, independent of the recurrence."""
    factors = tuple("P_%s_%d" % (side, chart_index) for chart_index in range(len(las_chart_names)))
    out = {("seed", side, order) + factors: Fraction(1, 1)}
    for seam_index in range(len(las_chart_names) - 1):
        suffix = factors[seam_index + 1 :]
        out[("seam", side, seam_index, order) + suffix] = Fraction(1, 1)
    return out


def las_endpoint_certificate():
    propagation_rows = []
    for side in las_sides:
        for order in range(las_order_max + 1):
            recurrence = las_propagation_bookkeeping(side, order)
            closed = las_closed_bookkeeping(side, order)
            residual = las_book_add(recurrence, las_book_scale(closed, -1))
            propagation_rows.append(
                {
                    "side": side,
                    "ordinary_Z_order": order,
                    "recurrence_terms": las_encode_book(recurrence),
                    "closed_seed_suffix_terms": las_encode_book(closed),
                    "residual": las_encode_book(residual),
                    "passed": not residual,
                }
            )

    difference_rows = []
    for order in range(las_order_max + 1):
        native_recurrence = las_propagation_bookkeeping("native", order)
        c3_recurrence = las_propagation_bookkeeping("C3", order)
        direct = las_book_add(native_recurrence, las_book_scale(c3_recurrence, -1))
        native_closed = las_closed_bookkeeping("native", order)
        c3_closed = las_closed_bookkeeping("C3", order)
        native_seed = {basis: coefficient for basis, coefficient in native_closed.items() if basis[0] == "seed"}
        c3_seed = {basis: coefficient for basis, coefficient in c3_closed.items() if basis[0] == "seed"}
        native_seams = {basis: coefficient for basis, coefficient in native_closed.items() if basis[0] == "seam"}
        c3_seams = {basis: coefficient for basis, coefficient in c3_closed.items() if basis[0] == "seam"}
        mixed_basis = ("mixed", "native_factor_times_C3_seed", order)
        decomposition = las_book_add(
            native_seed,
            las_book_scale(c3_seed, -1),
            {mixed_basis: Fraction(-1, 1)},
            {mixed_basis: Fraction(1, 1)},
            native_seams,
            las_book_scale(c3_seams, -1),
        )
        residual = las_book_add(direct, las_book_scale(decomposition, -1))
        difference_rows.append(
            {
                "ordinary_Z_order": order,
                "decomposition_identity": "F_native*(seed_native-seed_C3)+(F_native-F_C3)*seed_C3+signed_seam_terms",
                "direct_terminal_difference": las_encode_book(direct),
                "decomposed_seed_factor_and_seam_terms": las_encode_book(decomposition),
                "residual": las_encode_book(residual),
                "mixed_seed_cancellation_basis": list(mixed_basis),
                "passed": not residual,
            }
        )
    all_rows = propagation_rows + difference_rows
    return {
        "propagation_rows": propagation_rows,
        "difference_rows": difference_rows,
        "check_count": len(all_rows),
        "passed": all(row["passed"] for row in all_rows),
        "basis_is_formal": True,
        "numeric_evaluation": False,
    }


def las_difference(native, c3):
    seed_gap = [las_expr_sub(native["seed"][order], c3["seed"][order]) for order in range(las_order_max + 1)]
    factor_difference = las_expr_sub(native["factor"], c3["factor"])
    ordinary_rows = []
    seam_rows = []
    for seam_index in range(len(las_chart_names) - 1):
        native_eta = native["signed_seam_defects"][seam_index]["ordinary_Z_jets"]
        c3_eta = c3["signed_seam_defects"][seam_index]["ordinary_Z_jets"]
        native_suffix = native["suffix_internal"][seam_index]["factor"]
        c3_suffix = c3["suffix_internal"][seam_index]["factor"]
        seam_rows.append(
            {
                "seam_index": seam_index,
                "native_source_ref": native["signed_seam_defects"][seam_index]["source_ref"],
                "C3_source_ref": c3["signed_seam_defects"][seam_index]["source_ref"],
                "ordinary_Z_signed_terms": [
                    las_expr_sub(
                        las_expr_mul(native_suffix, native_eta[order]),
                        las_expr_mul(c3_suffix, c3_eta[order]),
                    )
                    for order in range(las_order_max + 1)
                ],
            }
        )
    for order in range(las_order_max + 1):
        seed_part = las_expr_mul(native["factor"], seed_gap[order])
        factor_mismatch_part = las_expr_mul(factor_difference, c3["seed"][order])
        signed_seams = [row["ordinary_Z_signed_terms"][order] for row in seam_rows]
        ordinary_rows.append(
            {
                "ordinary_Z_order": order,
                "seed_gap": seed_gap[order],
                "F_N_times_seed_gap": seed_part,
                "F_N_minus_F_C_times_C3_seed": factor_mismatch_part,
                "signed_seam_defect_sum": las_expr_add(*signed_seams),
                "general_endpoint_difference": las_expr_add(
                    seed_part,
                    factor_mismatch_part,
                    *signed_seams,
                ),
                "direct_terminal_difference": las_expr_sub(
                    native["ordinary_Z_terminal_expression"][order],
                    c3["ordinary_Z_terminal_expression"][order],
                ),
            }
        )
    return {
        "ordinary_Z_orders": list(range(las_order_max + 1)),
        "factor_difference": factor_difference,
        "seed_gap_jets": seed_gap,
        "signed_seam_defect_terms": seam_rows,
        "general_endpoint_difference_jets": ordinary_rows,
        "common_factor_specialization": {
            "unproved": True,
            "conditions": [
                "native and C3 action parameters are source-aligned",
                "native and C3 total log actions agree",
                "all five native and C3 signed interchart seam defects vanish",
                "the common action parameters are ordinary-Z independent",
            ],
            "F_common": las_expr_atom("F_common_conditional"),
            "F_Z_orders": [
                {"ordinary_Z_order": 0, "derivative": las_expr_atom("F_common_conditional")},
                {"ordinary_Z_order": 1, "derivative": 0},
                {"ordinary_Z_order": 2, "derivative": 0},
                {"ordinary_Z_order": 3, "derivative": 0},
            ],
            "a_j": [
                las_expr_mul(las_expr_atom("F_common_conditional"), seed_gap[order])
                for order in range(las_order_max + 1)
            ],
            "admission": False,
        },
    }


def derive_leading_amplitude_seam(source_contract, journal=None):
    """Construct the conditional six-chart seam ledger without source evaluation."""
    if not isinstance(source_contract, dict):
        raise TypeError("source_contract must be a mapping")
    native = las_side_ports("native", source_contract.get("native"))
    c3 = las_side_ports("C3", source_contract.get("C3"))
    premise = source_contract.get("premises")
    if not isinstance(premise, dict):
        raise ValueError("source_contract.premises is missing")
    common_factor_status = premise.get("common_factor")
    seam_status = premise.get("adjacent_endpoint_equalities")
    native_ledger = las_side_ledger("native", native)
    c3_ledger = las_side_ledger("C3", c3)
    difference = las_difference(native_ledger, c3_ledger)
    endpoint_certificate = las_endpoint_certificate()
    formal_certificate_passed = (
        endpoint_certificate["passed"]
        and native_ledger["closed_total_formula_passed"]
        and c3_ledger["closed_total_formula_passed"]
    )
    result = {
        "schema": "compact_leading_amplitude_seam_ledger_result_v2",
        "input_contract": las_input_contract,
        "formal_ledger_defined": True,
        "passed": formal_certificate_passed,
        "passed_scope": "symbolic six-step ledger and exact coefficient bookkeeping only",
        "formal_certificate_passed": formal_certificate_passed,
        "ordinary_Z_orders": [0, 1, 2, 3],
        "chart_order": list(las_chart_names),
        "action_parameter_contract": {
            "native": native_ledger["action_parameters"],
            "C3": c3_ledger["action_parameters"],
            "native_Z_independence_supplied": all(
                row["ordinary_Z_independent"] for row in native_ledger["action_parameters"].values()
            ),
            "C3_Z_independence_supplied": all(
                row["ordinary_Z_independent"] for row in c3_ledger["action_parameters"].values()
            ),
        },
        "native": native_ledger,
        "C3": c3_ledger,
        "symbolic_closed_sum": {
            "native_log_action": native_ledger["log_action_total"],
            "C3_log_action": c3_ledger["log_action_total"],
            "native_factor": native_ledger["factor"],
            "C3_factor": c3_ledger["factor"],
        },
        "signed_target_difference": difference,
        "endpoint_propagation_certificate": endpoint_certificate,
        "source_premise_status": {
            "common_factor_input_status": common_factor_status,
            "adjacent_endpoint_input_status": seam_status,
            "common_factor_admitted": False,
            "adjacent_endpoint_equalities_admitted": False,
            "seed_gap_forced_zero": False,
            "numeric_source_values_evaluated": False,
            "cap_or_midpoint_admitted": False,
        },
        "source_evaluation_counts": {
            "old_modules": 0,
            "old_proofs": 0,
            "old_normalizers": 0,
            "numeric_evaluations": 0,
        },
    }
    if journal is not None:
        journal(
            {
                "event": "leading_amplitude_seam_ledger_constructed",
                "ordinary_Z_orders": 4,
                "charts": 6,
                "interchart_seams_per_side": 5,
                "common_factor_admitted": False,
                "numeric_evaluation": False,
            }
        )
    return result



def oa_require(condition, message):
    if not condition:
        raise ValueError(message)


def oa_exact_json(value):
    """Encode internal rational coefficients without rounding or evaluation."""
    if isinstance(value, Fraction):
        return {'numerator': value.numerator, 'denominator': value.denominator}
    if isinstance(value, dict):
        return {key: oa_exact_json(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [oa_exact_json(child) for child in value]
    return value


def oa_creator_contract(records):
    """Bind the six amplitude formulas to the frozen native creator source."""
    import ast
    import hashlib
    native = records['native_outer']
    compact = records['native_compact']
    for item in (native, compact):
        oa_require(hashlib.sha256(item['source_text'].encode('utf-8')).hexdigest() == item['sha256'],
                   'Frozen creator source text/hash mismatch')
    formulas = (
        "E=a.scale(self.leading_Rh['fields']['E'],g.unary('exp',g.mul(g.constant('1/10'),g.add(x,g.constant(5)))))",
        "E=a.scale(self.background_endpoints['Rh_reference']['fields']['E'],g.unary('exp',g.sub(g.mul(g.constant('1/10'),x),g.mul(g.constant('3/5'),J))))",
        "E=a.scale(self.background_endpoints['O2_slope']['fields']['E'],g.unary('exp',g.mul(g.constant('-1/2'),distance)))",
        "E=a.scale(self.background_endpoints['O2_axial']['fields']['E'],g.unary('exp',g.mul(g.constant('-1/2'),x)))",
        "E=a.scale(self.background_endpoints['O2_buffer']['fields']['E'],g.unary('exp',g.sub(g.mul(g.constant('-1/2'),x),g.mul(self.mu,J))))",
        "E=a.scale(self.background_endpoints['O3_transition']['fields']['E'],g.unary('exp',g.neg(g.mul(g.add(g.constant('1/2'),self.mu),x))))")
    trees = [ast.parse(native['source_text']), ast.parse(compact['source_text'])]
    checks = []
    for chart, statement in zip(las_chart_names, formulas):
        target = ast.dump(ast.parse(statement).body[0], include_attributes=False)
        matches = [node for node in ast.walk(trees[0])
                   if isinstance(node, ast.Assign) and ast.dump(node, include_attributes=False) == target]
        oa_require(len(matches) == 1, 'Native amplitude assignment changed: ' + chart)
        checks.append({'chart': chart, 'creator_line': matches[0].lineno, 'matched': True})
    for statement, expected_count in (
        ("distance=g.sub(g.unary('exp',g.mul(g.constant(40),x)),g.one)", 1),
        ("J=controls.integral(g,sigma(g,t),g.nodes[t.node]['name'],g.zero,x,exact_original_sigma_integral=True)", 2)):
        target = ast.dump(ast.parse(statement).body[0], include_attributes=False)
        oa_require(sum(isinstance(node, ast.Assign) and ast.dump(node, include_attributes=False) == target
                       for node in ast.walk(trees[0])) == expected_count,
                   'Native axial coordinate or exact sigma action changed')
    for statement in ("self.g=g=outer.g", "self.A=outer.Rc_background['fields']['E']"):
        target = ast.dump(ast.parse(statement).body[0], include_attributes=False)
        oa_require(sum(isinstance(node, ast.Assign) and ast.dump(node, include_attributes=False) == target
                       for node in ast.walk(trees[1])) == 1, 'Compact amplitude owner binding changed')
    return {'passed': True, 'six_creator_assignments': checks,
            'shared_outer_graph_and_compact_amplitude_owner_matched': True,
            'source_functions_imported_or_called': False}


def oa_native_source_transport(spec):
    contract = spec['source_contract']['native']
    rows = spec['native_literal_rows']
    source_rows = []
    seed = spec['native_Rh_seed_nodes']
    oa_require(len(seed) == 4, 'Wrong native seed jet')
    oa_require(contract['seed']['ordinary_Z_ports'] == seed,
               'The action ledger does not use the actual native Rh seed')
    expected_domains = ((-5, 0), (0, 1), (0, 1), (0, 11), (0, 1), (0, 2))
    for chart, domain in zip(las_chart_names, expected_domains):
        port = spec['native_chart_ports'][chart]
        oa_require(tuple(port['domain_exact']) == domain, 'Native amplitude domain changed')
        body = port['body_nodes']
        outlet = port['outlet_nodes']
        variable, lower, upper = port['phase_node'], port['lower_node'], port['upper_node']
        for node, expected in ((lower, domain[0]), (upper, domain[1])):
            literal = rows[str(node)]
            oa_require(literal['operation'] == 'exact_rational' and
                       Fraction(literal['numerator'], literal['denominator']) == expected,
                       'Native source endpoint literal changed')
        oa_require(len(body) == len(outlet) == 4 and max(seed) < variable,
                   'Native scalar source variable is not fresh')
        factors = []
        for order in range(4):
            value = rows[str(body[order])]
            oa_require(value['operation'] == 'product' and len(value['arguments']) == 2 and
                       seed[order] in value['arguments'], 'Native amplitude is not the prescribed scalar scale')
            arguments = value['arguments']
            factor = arguments[1] if arguments[0] == seed[order] else arguments[0]
            factors.append(factor)
            end = rows[str(outlet[order])]
            oa_require(end['operation'] == 'function_substitution' and end['expression'] == body[order] and
                       end['variable'] == variable and end['value'] == upper and
                       end['Z_independent_substitution'] is True, 'Wrong native amplitude outlet query')
            query = contract['charts'][chart]['inlet']['ordinary_Z_ports'][order]
            oa_require(query['archivehash'] == contract['archive_raw_sha256'] and
                       query['functionnode'] == 'graph:' + str(body[order]) and
                       query['querycoordinate'] == {'variable_node': variable, 'value_node': lower} and
                       query['order'] == order, 'A free radial body was used as an inlet value')
            oa_require(contract['charts'][chart]['outlet']['ordinary_Z_ports'][order] == outlet[order],
                       'Native terminal amplitude handle changed')
            source_rows.append({'chart': chart, 'ordinary_Z_order': order,
                                'predecessor': seed[order], 'body': body[order],
                                'outlet': outlet[order], 'scale_row_passed': True,
                                'outlet_substitution_passed': True})
        oa_require(len(set(factors)) == 1, 'The ordinary-Z rows have different attenuation factors')
        oa_require(factors[0] == spec['native_scalar_action_ports'][chart]['factor_node'],
                   'The body attenuation factor is not the action source factor')
        factor_row = rows[str(factors[0])]
        oa_require(factor_row['operation'] == 'analytic_unary' and factor_row['name'] == 'exp',
                   'Native attenuation is not the prescribed exponential factor')
        seed = outlet
    oa_require(seed == spec['target_native_amplitude_nodes'], 'The action does not end at the actual target amplitude')
    prefix = spec['native_outer_C5_prefix_authentication']
    serialization = 'JSON sort_keys=True, compact separators, UTF-8; no normalization or comparison'
    oa_require(prefix['outer_archive_raw_sha256'] == contract['archive_raw_sha256'] and
               prefix['C5_archive_raw_sha256'] == '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269' and
               prefix['outer_serialization'] == prefix['C5_serialization'] == serialization,
               'Native outer/compact prefix archive or serialization scope changed')
    oa_require(prefix['outer_node_count'] == prefix['C5_prefix_node_count'] == 258857 and
               prefix['outer_prefix_sha256'] == prefix['C5_prefix_sha256'],
               'Native outer source prefix does not equal the compact source prefix')
    return {'passed': True, 'body_scale_checks': 24, 'outlet_substitution_checks': 24,
            'rows': source_rows, 'native_outer_to_C5_amplitude_source_alias_proved': True,
            'literal_prefix_authentication': prefix,
            'seed_not_assumed_canonical': True}


def oa_literal_scalar_source_certificate(spec):
    """Bind fixed source IDs by a finite literal table, without a tree matcher.

    This contract uses ordered source edges and exact stored operations.  It
    does not normalize expressions, compare rational values, traverse a DAG,
    or call the withdrawn recursive matching family.  The six formulas below
    are read directly in the authenticated creator AST and these literal rows.
    """
    def rational(p, q=1):
        return {'operation': 'exact_rational', 'numerator': p, 'denominator': q}

    def binary(operation, left, right):
        return {'operation': operation, 'arguments': [left, right]}

    def unary(name, argument):
        return {'operation': 'analytic_unary', 'name': name, 'argument': argument}

    def negative(argument):
        return {'operation': 'negative', 'argument': argument}

    def bound(name):
        return {'operation': 'bound_variable', 'name': name}

    expected = {
        0: rational(0), 1: rational(1), 2: rational(40),
        3: unary('exp', 2), 4: rational(11), 5: binary('sum', 3, 4),
        11: rational(1, 1000), 12: unary('log', 11), 13: rational(4),
        14: binary('product', 13, 5), 15: negative(14),
        16: binary('sum', 12, 15), 17: unary('exp', 16),
        24: rational(-1, 2), 28: rational(1, 2), 2302: rational(5),
        25085: negative(1), 31355: rational(1, 10), 31357: rational(3, 5),
        4993: binary('sum', 28, 17),
        205148: bound('true_outer_Rh_reference'),
        205158: binary('sum', 205148, 2302),
        205159: binary('product', 31355, 205158), 205160: unary('exp', 205159),
        214116: bound('true_outer_O2_slope'),
        214124: bound('true_outer_slope_sigma_t'),
        214127: binary('product', 31355, 214116),
        214128: binary('product', 31357, 214126), 214129: negative(214128),
        214130: binary('sum', 214127, 214129), 214131: unary('exp', 214130),
        223058: bound('true_outer_O2_axial'),
        223060: binary('product', 2, 223058), 223061: unary('exp', 223060),
        223067: binary('sum', 223061, 25085),
        223068: binary('product', 24, 223067), 223069: unary('exp', 223068),
        232895: bound('true_outer_O2_buffer'),
        232903: binary('product', 24, 232895), 232904: unary('exp', 232903),
        241536: bound('true_outer_O3_transition'),
        241544: bound('true_outer_transition_sigma_t'),
        241547: binary('product', 24, 241536),
        241548: binary('product', 17, 241546), 241549: negative(241548),
        241550: binary('sum', 241547, 241549), 241551: unary('exp', 241550),
        250215: bound('true_outer_O3_power'),
        250223: binary('product', 4993, 250215), 250224: negative(250223),
        250225: unary('exp', 250224),
    }
    sigma_binding = {
        'AST_sha256': '4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c',
        'function': 'stable_sigma', 'module': 'lei_ren_part1_paper_compliant_outer_initial.py',
        'module_sha256': 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657'}
    for j, sigma, variable, phase, name in (
        (214126, 214125, 214124, 214116, 'true_outer_slope_sigma_t'),
        (241546, 241545, 241544, 241536, 'true_outer_transition_sigma_t')):
        expected[sigma] = {
            'operation': 'exact_original_source_flat_sigma', 'argument': variable,
            'definition': '0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside',
            'exact_function_not_its_interval_cap': True, 'source_binding': sigma_binding}
        expected[j] = {
            'operation': 'definite_integral', 'exact_function_integral': True,
            'exact_original_sigma_integral': True, 'numerical_value_not_installed': True,
            'integrand': sigma, 'lower': 0, 'upper': phase, 'variable': name}
    rows = spec['native_literal_rows']
    literal_checks = []
    for node, literal in expected.items():
        oa_require(rows.get(str(node)) == literal, 'Scalar literal source-row contract changed: ' + str(node))
        literal_checks.append({'node': node, 'literal_row': literal, 'passed': True})
    chart_table = (
        ('Rh_reference', 205148, 205159, 205160, [], '(x+5)/10'),
        ('O2_slope', 214116, 214130, 214131, [214126], 'x/10-(3/5)*J_slope(x)'),
        ('O2_axial', 223058, 223068, 223069, [], '-(exp(40*x)-1)/2'),
        ('O2_buffer', 232895, 232903, 232904, [], '-x/2'),
        ('O3_transition', 241536, 241550, 241551, [241546], '-x/2-mu*J_trans(x)'),
        ('O3_power', 250215, 250224, 250225, [], '-(1/2+mu)*x'))
    proofs = []
    for chart, phase, exponent, factor, integrals, formula in chart_table:
        ports = spec['native_scalar_action_ports'][chart]
        oa_require(ports['phase_node'] == phase == spec['native_chart_ports'][chart]['phase_node'] and
                   ports['exponent_node'] == exponent and ports['factor_node'] == factor and
                   ports['J_nodes'] == integrals, 'Scalar literal source-port contract changed: ' + chart)
        proofs.append({'chart': chart, 'factor_node': factor, 'exponent_node': exponent,
                       'source_integral_nodes': integrals, 'source_formula': formula, 'passed': True})
    parameters = spec['source_contract']['native']['action_parameters']
    parameter_refs = {'I_slope': 'graph:214126@214116=1', 'I_trans': 'graph:241546@241536=1',
                      'mu': 'graph:17', 'exp40': 'exact exp(40)'}
    for name, reference in parameter_refs.items():
        oa_require(parameters[name]['source_ref'] == reference, 'Wrong literal action parameter source: ' + name)
    return {'passed': True, 'scalar_exponent_source_checks': 6, 'rows': proofs,
            'literal_row_checks': literal_checks, 'literal_row_check_count': len(literal_checks),
            'source_intervals_or_midpoints_used': False,
            'recursive_source_matcher_calls': 0, 'rational_value_comparison_calls': 0,
            'native_mu_and_exp40_Z_independence_from_literal_constant_rows': True}


def oa_inlet_exponent_certificate(spec):
    rows = []
    for index, chart in enumerate(las_chart_names):
        lower = Fraction(spec['native_chart_ports'][chart]['domain_exact'][0])
        if index == 0:
            constant, mu_coefficient = (lower + 5) / 10, Fraction(0)
            rule = 'exp((x+5)/10) at x=-5'
        elif index == 1:
            constant, mu_coefficient = lower / 10, Fraction(0)
            rule = 'integral_0^0 sigma = 0; exp(0/10-(3/5)*0) = 1'
        elif index == 2:
            oa_require(lower == 0, 'Exponential axial coordinate must start at zero')
            constant, mu_coefficient = Fraction(-1, 2) * (1 - 1), Fraction(0)
            rule = 'exp(40*0)=1; exp(-(1-1)/2)=1'
        elif index == 3:
            constant, mu_coefficient = -lower / 2, Fraction(0)
            rule = 'exp(-x/2) at x=0'
        elif index == 4:
            constant, mu_coefficient = -lower / 2, Fraction(0)
            rule = 'integral_0^0 sigma = 0; exp(-0/2-mu*0)=1'
        else:
            constant, mu_coefficient = -lower / 2, -lower
            rule = 'exp(-(1/2+mu)*0)=1'
        oa_require(constant == mu_coefficient == 0, 'Nonzero inlet attenuation action')
        rows.append({'chart': chart, 'rule': rule, 'log_factor_constant': 0,
                     'log_factor_mu_coefficient': 0, 'factor': 1,
                     'ordinary_Z_orders': [0, 1, 2, 3], 'passed': True})
    return {'passed': True, 'inlet_factor_identities': rows, 'count': 6,
            'native_interchart_seams_proved': 5,
            'ordinary_Z_function_identity_count': 24,
            'proof_rule': 'The exact scalar inlet factor is one; each source body scales the actual previous outlet jet. Fresh Z-independent radial substitution preserves all four ordinary-Z rows.',
            'numerical_intervals_or_source_values_used': False,
            'C3_interchart_seams_admitted': False}


def derive(source, journal):
    spec = source['input']
    oa_require(spec['schema'] == 'compact_outer_amplitude_action_transport_input_v1', 'Wrong new action input')
    creator = oa_creator_contract(spec['creator_source_records'])
    native = oa_native_source_transport(spec)
    scalar = oa_literal_scalar_source_certificate(spec)
    inlet = oa_inlet_exponent_certificate(spec)
    journal({'event': 'new_native_amplitude_source_transport', 'scale_checks': 24,
             'outlet_checks': 24, 'inlet_factor_identities': 6})
    ledger = derive_leading_amplitude_seam(spec['source_contract'], journal)
    oa_require(ledger['passed'] and ledger['endpoint_propagation_certificate']['check_count'] == 12,
               'New amplitude action ledger failed')
    parent = source['parent_result']
    result = {'schema': 'original_compact_outer_amplitude_action_transport',
        'component': 'compact_outer_amplitude_action_transport', 'passed': True,
        'metadata_only': False, 'source_only': True, 'input_digest': source['input_digest'],
        'runtime_model_metadata': None, 'original_source_evaluations': 0,
        'source_packet_defined': True, 'amplitude_action_source_defined': True,
        'unaccepted_preentry_worker_consumer': spec['unaccepted_preentry_worker_consumer'],
        'six_chart_log_action_identity_proved': True,
        'amplitude_signed_seed_seam_transport_C3_defined': True,
        'native_closed_Rc_amplitude_action_proved': True,
        'native_five_interchart_amplitude_seams_ordinary_Z3_proved': True,
        'native_outer_to_C5_amplitude_source_alias_proved': True,
        'actual_native_C3_amplitude_identity_proved': False,
        'actual_C3_C5_incoming_source_identity_proved': False,
        'actual_C3_C5_control_source_identity_proved': False,
        'actual_C3_C5_profile_identity_proved': False,
        'actual_C3_cap_to_C5_source_bound_proved': False,
        'compact_prefix_source_binding_proved': False,
        'compact_volume_cap_same_source_bound_proved': False,
        'actual_common_C1_control_ball_admitted': False, 'signed_target_gaps_evaluated': False,
        'native_creator_contract': creator, 'native_source_transport_certificate': native,
        'native_scalar_action_certificate': scalar,
        'native_inlet_seam_certificate': inlet, 'amplitude_action_certificate': ledger,
        'native_amplitude_consequence': {
            'source_ports': spec['target_native_amplitude_nodes'],
            'actual_Rh_seed_ports': spec['native_Rh_seed_nodes'],
            'ordinary_Z_orders': [0, 1, 2, 3],
            'function_identity': 'A_N^(j)(Z)=exp(L_N)*E_native_Rh^(j)(Z), j=0..3',
            'L_N': ledger['native']['log_action_total'],
            'native_seed_canonical_identity_proved': False,
            'C3_seam_terms_retained': True,
            'signed_seed_and_factor_differences_retained': True},
        'source_binding': {'passed': True, 'native_action_source_identity_only': True,
            'actual_native_C3_amplitude_identity_proved': False,
            'source_domains_and_action_parameters_Z_independent': True,
            'C3_action_transport_conditional_on_frozen_frontend_contract': True,
            'archived_caps_used_as_values': False, 'old_errors_removed': False},
        'source_inspection_audit': {'original_source_evaluations': 0, 'old_source_consumer_calls': 0,
            'old_proof_functions_called': 0, 'old_operator_calls': 0, 'old_source_normalizer_calls': 0,
            'quadratures': 0, 'numeric_field_evaluations': 0, 'temporal_recurrence_calls': 0,
            'new_native_body_scale_checks': 24, 'new_native_outlet_substitution_checks': 24,
            'new_native_scalar_exponent_source_checks': 6,
            'new_native_inlet_factor_identities': 6, 'new_action_bookkeeping_checks': 12},
        'remaining_scope': ['Prove or bound the native Rh seed minus the C3 canonical seed as an actual signed source function.',
            'Bind the native/C3 action parameters and discharge C3 frontend seams before setting the amplitude gap to zero.',
            'Use the retained joint target numerators with the actual signed raw incoming history differences; do not replace them with caps or midpoints.',
            'Functional terminal closure, signed Evol and independent reference offset, and all later physical/time/NS layers remain incomplete.']}
    for name in ('retained_cap_provenance', 'retained_construction_data', 'retained_integral_term_inventory',
                 'coordinate_mapping', 'native_BG_exact_numeric_baseline_C3_identity_proved',
                 'relative_mass_C3_source_identity_proved', 'compact_prefix_geometry_logC_binding_proved',
                 'compact_prefix_fixed_mstar_enclosure_proved', 'compact_prefix_12channel_2volume_dominance_proved',
                 'joint_target_numerator_source_defined', 'five_joint_numerator_identities_proved',
                 'ordinary_Z_quotient_transport_proved'):
        result[name] = parent[name]
    for name in source['inherited_false_flags']:
        result[name] = False
    journal({'event': 'new_outer_amplitude_action_complete', 'actual_native_C3_amplitude_equality': False})
    return oa_exact_json(result)
