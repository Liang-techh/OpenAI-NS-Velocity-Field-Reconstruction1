"""Fresh signed history-density defect algebra.

The callable in this draft uses integer sparse exponent-vector polynomials and
ordinary-Z Leibniz jets through order three.  It compares two independent
source packets, native and C3, while retaining every signed field difference.
The e and p channels are represented in doubled form so every coefficient is
an integer.  A separate paired Duhamel identity records two-chain history
transport without identifying memories, densities, caps, or source values.

This module is a source-contract draft only.  It does not import a provider,
replay an entered proof, evaluate a finite-N value, or admit a physical/global
claim.  The root wrapper must provide the source ports, units, chart names,
and flags before invoking the callable.
"""


sdh_schema = "compact_signed_history_density_algebra_input_v1"
sdh_ordinary_orders = (0, 1, 2, 3)
sdh_field_names = (
    "En",
    "Ec",
    "Vn",
    "Vc",
    "Fn",
    "Fc",
    "Bn",
    "Bc",
)
sdh_channel_names = ("m", "h", "k", "e2", "p2")
sdh_basis_note = "sparse monomial exponent maps with integer coefficients"

sdh_factorized_formulas = {
    "m": "dB",
    "h": "dF",
    "k": "Vn*dF+dV*Fc+En*dB+dE*Bc+eps*(Fn*dB+dF*Bc)",
    "e2": "4*(Vn*dB+dV*Bc)-2*(En*dF+dE*Fc)+eps*(2*(Bn+Bc)*dB-(Fn+Fc)*dF)",
    "p2": "2*(En*dF+dE*Fc)+eps*(Fn+Fc)*dF",
}
sdh_unscaled_factorized_formulas = {
    "m": "dB",
    "h": "dF",
    "k": "Vn*dF+dV*Fc+En*dB+dE*Bc+eps*(Fn*dB+dF*Bc)",
    "e": "2*(Vn*dB+dV*Bc)-(En*dF+dE*Fc)+eps*((Bn+Bc)*dB-(Fn+Fc)*dF/2)",
    "p": "En*dF+dE*Fc+eps*(Fn+Fc)*dF/2",
}

sdh_contract_template = {
    "schema": sdh_schema,
    "required_fields": [
        "schema",
        "source_binding",
        "field_source_ports",
        "history_transport",
        "flags",
    ],
    "field_names": [
        "En",
        "Ec",
        "Vn",
        "Vc",
        "Fn",
        "Fc",
        "Bn",
        "Bc",
    ],
    "ordinary_Z_contract": {
        "orders": [0, 1, 2, 3],
        "jet_kind": "ordinary Z derivatives",
        "source_row_count_per_field": 4,
        "source_calls_at_this_layer": 0,
    },
    "density_formula_contract": {
        "native_and_C3_form": "(B,F,V*F+E*B+eps*F*B,2*V*B-E*F+eps*(B^2-F^2/2),E*F+eps*F^2/2)",
        "integer_channel_form": {
            "m": "B",
            "h": "F",
            "k": "V*F+E*B+eps*F*B",
            "e2": "4*V*B-2*E*F+eps*(2*B^2-F^2)",
            "p2": "2*E*F+eps*F^2",
        },
        "scaled_channels": "e2=2*e and p2=2*p",
        "ordinary_Z_rule": "D^k(f*g)=sum(binomial(k,j)*D^j(f)*D^(k-j)(g), k=0..3)",
    },
    "source_binding": {
        "root_wrapper_passed": True,
        "eps_definition": "eps=1/N",
        "eps_source_ref": "root-authenticated inverse-N source contract",
        "units_source_ref": "root-authenticated native/C3 physical units",
        "numeric_source_values_selected": False,
        "source_calls_allowed": False,
        "caps_admitted": False,
    },
    "field_source_ports": {
        "En": {
            "source_ref": "root supplies four ordinary-Z source ports for En",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-En-Z0", "root-port-En-Z1", "root-port-En-Z2", "root-port-En-Z3"],
        },
        "Ec": {
            "source_ref": "root supplies four ordinary-Z source ports for Ec",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Ec-Z0", "root-port-Ec-Z1", "root-port-Ec-Z2", "root-port-Ec-Z3"],
        },
        "Vn": {
            "source_ref": "root supplies four ordinary-Z source ports for Vn",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Vn-Z0", "root-port-Vn-Z1", "root-port-Vn-Z2", "root-port-Vn-Z3"],
        },
        "Vc": {
            "source_ref": "root supplies four ordinary-Z source ports for Vc",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Vc-Z0", "root-port-Vc-Z1", "root-port-Vc-Z2", "root-port-Vc-Z3"],
        },
        "Fn": {
            "source_ref": "root supplies four ordinary-Z source ports for Fn",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Fn-Z0", "root-port-Fn-Z1", "root-port-Fn-Z2", "root-port-Fn-Z3"],
        },
        "Fc": {
            "source_ref": "root supplies four ordinary-Z source ports for Fc",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Fc-Z0", "root-port-Fc-Z1", "root-port-Fc-Z2", "root-port-Fc-Z3"],
        },
        "Bn": {
            "source_ref": "root supplies four ordinary-Z source ports for Bn",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Bn-Z0", "root-port-Bn-Z1", "root-port-Bn-Z2", "root-port-Bn-Z3"],
        },
        "Bc": {
            "source_ref": "root supplies four ordinary-Z source ports for Bc",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "ordinary_Z_ports": ["root-port-Bc-Z0", "root-port-Bc-Z1", "root-port-Bc-Z2", "root-port-Bc-Z3"],
        },
    },
    "history_transport": {
        "chart_names": "root supplies arbitrary chart-name list; no fixed six-chart assumption",
        "history_channels": ["channel_0", "channel_1", "channel_2", "channel_3", "channel_4"],
        "identity_scope": "one generic five-channel identity applies independently on every supplied chart",
        "symbolic_check_rows": 20,
        "memory_equality_assumed": False,
        "density_equality_assumed": False,
    },
    "flags": {
        "global_gap_zero": False,
        "numeric_values_selected": False,
        "cap_reduction_admitted": False,
        "physical_completion_proved": False,
    },
}


def sdh_require(condition, message):
    if not condition:
        raise ValueError(message)


def sdh_text(value, label):
    sdh_require(isinstance(value, str) and value != "", label + " must be nonempty text")
    return value


def sdh_zero():
    return {}


def sdh_constant(value):
    sdh_require(type(value) is int, "polynomial constants must be integers")
    return {} if value == 0 else {(): value}


def sdh_atom(name):
    return {((str(name), 1),): 1}


def sdh_normalize_monomial(pairs):
    powers = {}
    for name, power in pairs:
        if power:
            powers[name] = powers.get(name, 0) + power
    return tuple(sorted((name, power) for name, power in powers.items() if power))


def sdh_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] = result.get(monomial, 0) + coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def sdh_scale(polynomial, scalar):
    sdh_require(type(scalar) is int, "polynomial scale must be an integer")
    if scalar == 0:
        return {}
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in polynomial.items()
        if coefficient * scalar
    }


def sdh_mul(left, right):
    result = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = sdh_normalize_monomial(left_monomial + right_monomial)
            coefficient = left_coefficient * right_coefficient
            result[monomial] = result.get(monomial, 0) + coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def sdh_sub(left, right):
    return sdh_add(left, sdh_scale(right, -1))


def sdh_encode(polynomial):
    terms = []
    for monomial, coefficient in sorted(polynomial.items()):
        terms.append({
            "exponents": [
                {"atom": name, "power": power}
                for name, power in monomial
            ],
            "coefficient": coefficient,
        })
    return terms


def sdh_is_zero(polynomial):
    return not polynomial


def sdh_binomial(order, position):
    table = (
        (1,),
        (1, 1),
        (1, 2, 1),
        (1, 3, 3, 1),
    )
    return table[order][position]


def sdh_jet_atom(prefix):
    return [sdh_atom(prefix + "_Z" + str(order)) for order in sdh_ordinary_orders]


def sdh_jet_constant_atom(name):
    """Represent a Z-independent scalar as a value atom and zero derivatives."""
    return [sdh_atom(name), sdh_zero(), sdh_zero(), sdh_zero()]


def sdh_jet_add(*jets):
    return [
        sdh_add(*(jet[order] for jet in jets))
        for order in range(len(sdh_ordinary_orders))
    ]


def sdh_jet_scale(jet, scalar):
    return [sdh_scale(value, scalar) for value in jet]


def sdh_jet_sub(left, right):
    return sdh_jet_add(left, sdh_jet_scale(right, -1))


def sdh_jet_mul(left, right):
    return [
        sdh_add(*(
            sdh_scale(
                sdh_mul(left[position], right[order - position]),
                sdh_binomial(order, position),
            )
            for position in range(order + 1)
        ))
        for order in sdh_ordinary_orders
    ]


def sdh_validate_contract(contract):
    """Validate ports and flags without reading or evaluating source rows."""
    sdh_require(isinstance(contract, dict), "contract must be a mapping")
    sdh_require(contract.get("schema") == sdh_schema, "wrong signed density schema")

    binding = contract.get("source_binding")
    sdh_require(isinstance(binding, dict), "source_binding is missing")
    sdh_require(binding.get("root_wrapper_passed") is True,
                "root wrapper has not passed source binding")
    sdh_require(binding.get("eps_definition") == "eps=1/N", "eps definition changed")
    sdh_text(binding.get("eps_source_ref"), "source_binding.eps_source_ref")
    sdh_text(binding.get("units_source_ref"), "source_binding.units_source_ref")
    for name in ("numeric_source_values_selected", "source_calls_allowed", "caps_admitted"):
        sdh_require(binding.get(name) is False, "forbidden source binding flag is true: " + name)

    ports = contract.get("field_source_ports")
    sdh_require(isinstance(ports, dict), "field_source_ports are missing")
    for name in sdh_field_names:
        row = ports.get(name)
        sdh_require(isinstance(row, dict), "missing field source ports: " + name)
        sdh_text(row.get("source_ref"), name + ".source_ref")
        sdh_require(row.get("ordinary_Z_orders") == list(sdh_ordinary_orders),
                    "wrong ordinary-Z orders for " + name)
        raw_ports = row.get("ordinary_Z_ports")
        sdh_require(isinstance(raw_ports, list) and len(raw_ports) == 4,
                    "four ordinary-Z ports required for " + name)
        for order, port in enumerate(raw_ports):
            sdh_text(port, "%s.ordinary_Z_ports[%d]" % (name, order))

    history = contract.get("history_transport")
    sdh_require(isinstance(history, dict), "history_transport is missing")
    chart_names = history.get("chart_names")
    sdh_require(isinstance(chart_names, list) and chart_names,
                "root must supply a nonempty chart-name list")
    for index, chart in enumerate(chart_names):
        sdh_text(chart, "history_transport.chart_names[%d]" % index)
    history_channels = history.get("history_channels")
    sdh_require(isinstance(history_channels, list) and len(history_channels) == 5,
                "five generic history channels are required")
    for index, channel in enumerate(history_channels):
        sdh_text(channel, "history_transport.history_channels[%d]" % index)
    sdh_require(history.get("memory_equality_assumed") is False,
                "memory equality must remain unassumed")
    sdh_require(history.get("density_equality_assumed") is False,
                "density equality must remain unassumed")

    flags = contract.get("flags")
    sdh_require(isinstance(flags, dict), "flags are missing")
    for name in (
        "global_gap_zero",
        "numeric_values_selected",
        "cap_reduction_admitted",
        "physical_completion_proved",
    ):
        sdh_require(flags.get(name) is False, "forbidden flag is true: " + name)

    return {
        "schema": sdh_schema,
        "field_source_port_count": len(sdh_field_names) * 4,
        "chart_count": len(chart_names),
        "chart_names": list(chart_names),
        "history_channels": list(history_channels),
        "ordinary_Z_orders": list(sdh_ordinary_orders),
        "source_calls_allowed": False,
    }


def sdh_density_channels(E, V, F, B, eps):
    """Return doubled integer-coefficient source density channels."""
    return {
        "m": B,
        "h": F,
        "k": sdh_jet_add(
            sdh_jet_mul(V, F),
            sdh_jet_mul(E, B),
            sdh_jet_mul(eps, sdh_jet_mul(F, B)),
        ),
        "e2": sdh_jet_add(
            sdh_jet_scale(sdh_jet_mul(V, B), 4),
            sdh_jet_scale(sdh_jet_mul(E, F), -2),
            sdh_jet_mul(eps, sdh_jet_add(
                sdh_jet_scale(sdh_jet_mul(B, B), 2),
                sdh_jet_scale(sdh_jet_mul(F, F), -1),
            )),
        ),
        "p2": sdh_jet_add(
            sdh_jet_scale(sdh_jet_mul(E, F), 2),
            sdh_jet_mul(eps, sdh_jet_mul(F, F)),
        ),
    }


def sdh_factorized_gaps(En, Ec, Vn, Vc, Fn, Fc, Bn, Bc, eps):
    dE = sdh_jet_sub(En, Ec)
    dV = sdh_jet_sub(Vn, Vc)
    dF = sdh_jet_sub(Fn, Fc)
    dB = sdh_jet_sub(Bn, Bc)
    return {
        "m": dB,
        "h": dF,
        "k": sdh_jet_add(
            sdh_jet_mul(Vn, dF),
            sdh_jet_mul(dV, Fc),
            sdh_jet_mul(En, dB),
            sdh_jet_mul(dE, Bc),
            sdh_jet_mul(eps, sdh_jet_add(
                sdh_jet_mul(Fn, dB),
                sdh_jet_mul(dF, Bc),
            )),
        ),
        "e2": sdh_jet_add(
            sdh_jet_scale(sdh_jet_add(
                sdh_jet_mul(Vn, dB),
                sdh_jet_mul(dV, Bc),
            ), 4),
            sdh_jet_scale(sdh_jet_add(
                sdh_jet_mul(En, dF),
                sdh_jet_mul(dE, Fc),
            ), -2),
            sdh_jet_mul(eps, sdh_jet_add(
                sdh_jet_scale(sdh_jet_mul(sdh_jet_add(Bn, Bc), dB), 2),
                sdh_jet_scale(sdh_jet_mul(sdh_jet_add(Fn, Fc), dF), -1),
            )),
        ),
        "p2": sdh_jet_add(
            sdh_jet_scale(sdh_jet_add(
                sdh_jet_mul(En, dF),
                sdh_jet_mul(dE, Fc),
            ), 2),
            sdh_jet_mul(eps, sdh_jet_mul(sdh_jet_add(Fn, Fc), dF)),
        ),
    }


def sdh_history_rows(history_channels):
    rows = []
    for channel in history_channels:
        MN = sdh_jet_atom("history_" + channel + "_MN")
        MC = sdh_jet_atom("history_" + channel + "_MC")
        DN = sdh_jet_atom("history_" + channel + "_DN_in")
        DC = sdh_jet_atom("history_" + channel + "_DC_in")
        IN = sdh_jet_atom("history_" + channel + "_IN")
        IC = sdh_jet_atom("history_" + channel + "_IC")
        delta_in = sdh_jet_sub(DN, DC)
        left = sdh_jet_sub(
            sdh_jet_add(sdh_jet_mul(MN, DN), IN),
            sdh_jet_add(sdh_jet_mul(MC, DC), IC),
        )
        right = sdh_jet_add(
            sdh_jet_mul(MN, delta_in),
            sdh_jet_mul(sdh_jet_sub(MN, MC), DC),
            sdh_jet_sub(IN, IC),
        )
        for order in sdh_ordinary_orders:
            residual = sdh_sub(left[order], right[order])
            rows.append({
                "history_channel": channel,
                "ordinary_Z_order": order,
                "residual": sdh_encode(residual),
                "passed": sdh_is_zero(residual),
                "formula": "Delta_out=MN*Delta_in+(MN-MC)*DC_in+(IN-IC)",
            })
    return rows


def derive_signed_density_defect_algebra(contract, journal=None):
    """Build signed density defects and paired history transport symbolically."""
    checked = sdh_validate_contract(contract)
    if journal is not None:
        sdh_require(callable(journal), "journal must be callable")

    En = sdh_jet_atom("En")
    Ec = sdh_jet_atom("Ec")
    Vn = sdh_jet_atom("Vn")
    Vc = sdh_jet_atom("Vc")
    Fn = sdh_jet_atom("Fn")
    Fc = sdh_jet_atom("Fc")
    Bn = sdh_jet_atom("Bn")
    Bc = sdh_jet_atom("Bc")
    eps = sdh_jet_constant_atom("eps")

    native = sdh_density_channels(En, Vn, Fn, Bn, eps)
    c3 = sdh_density_channels(Ec, Vc, Fc, Bc, eps)
    direct_gap = {
        name: sdh_jet_sub(native[name], c3[name])
        for name in sdh_channel_names
    }
    factorized_gap = sdh_factorized_gaps(En, Ec, Vn, Vc, Fn, Fc, Bn, Bc, eps)

    residual_checks = []
    factorized_rows = []
    for channel in sdh_channel_names:
        for order in sdh_ordinary_orders:
            direct = direct_gap[channel][order]
            factorized = factorized_gap[channel][order]
            residual = sdh_sub(direct, factorized)
            check = {
                "channel": channel,
                "ordinary_Z_order": order,
                "residual": sdh_encode(residual),
                "passed": sdh_is_zero(residual),
            }
            residual_checks.append(check)
            factorized_rows.append({
                "channel": channel,
                "ordinary_Z_order": order,
                "factorized_gap": sdh_encode(factorized),
                "direct_signed_difference": sdh_encode(direct),
                "residual_zero": sdh_is_zero(residual),
            })
    sdh_require(len(residual_checks) == 20, "wrong density residual count")
    sdh_require(len(factorized_rows) == 20, "wrong factorized gap row count")
    sdh_require(all(row["passed"] for row in residual_checks),
                "signed density factorization residual is nonzero")

    shared_sum = sdh_jet_add(factorized_gap["e2"], factorized_gap["p2"])
    dE = sdh_jet_sub(En, Ec)
    dV = sdh_jet_sub(Vn, Vc)
    dB = sdh_jet_sub(Bn, Bc)
    shared_expected = sdh_jet_add(
        sdh_jet_scale(sdh_jet_add(
            sdh_jet_mul(Vn, dB),
            sdh_jet_mul(dV, Bc),
        ), 4),
        sdh_jet_mul(eps, sdh_jet_scale(sdh_jet_mul(sdh_jet_add(Bn, Bc), dB), 2)),
    )
    shared_rows = []
    for order in sdh_ordinary_orders:
        shared_residual = sdh_sub(shared_sum[order], shared_expected[order])
        shared_rows.append({
            "ordinary_Z_order": order,
            "residual": sdh_encode(shared_residual),
            "passed": sdh_is_zero(shared_residual),
        })
    sdh_require(all(row["passed"] for row in shared_rows),
                "shared e2+p2 correlation residual is nonzero")

    history_channels = contract["history_transport"]["history_channels"]
    history_rows = sdh_history_rows(history_channels)
    sdh_require(len(history_rows) == 20, "wrong paired history check count")
    sdh_require(all(row["passed"] for row in history_rows),
                "paired history Duhamel residual is nonzero")

    if journal is not None:
        journal({
            "event": "signed_history_density_algebra_defined",
            "density_residual_checks": 20,
            "factorized_gap_rows": 20,
            "paired_history_checks": 20,
            "ordinary_Z_orders": 4,
            "numeric_evaluation": False,
        })

    return {
        "schema": "compact_signed_history_density_algebra_result_v1",
        "passed": True,
        "conditional_source_contract": checked,
        "basis": sdh_basis_note,
        "ordinary_Z_orders": list(sdh_ordinary_orders),
        "signed_field_differences": {
            "dE": "En-Ec",
            "dV": "Vn-Vc",
            "dF": "Fn-Fc",
            "dB": "Bn-Bc",
            "eps": "1/N",
        },
        "scaled_density_channels": {
            "m": "Bn",
            "h": "Fn",
            "k": "Vn*Fn+En*Bn+eps*Fn*Bn",
            "e2": "4*Vn*Bn-2*En*Fn+eps*(2*Bn^2-Fn^2)",
            "p2": "2*En*Fn+eps*Fn^2",
            "scale_note": "e=e2/2 and p=p2/2; integer coefficient maps are used",
        },
        "factorized_signed_defects": {
            "formulas": dict(sdh_factorized_formulas),
            "unscaled_formulas": dict(sdh_unscaled_factorized_formulas),
            "physical_channel_names": ["m", "h", "k", "e", "p"],
            "doubled_internal_channels": ["m", "h", "k", "e2", "p2"],
            "rows": factorized_rows,
            "row_count": len(factorized_rows),
            "all_rows_passed": True,
        },
        "residual_proof": {
            "checks": residual_checks,
            "check_count": len(residual_checks),
            "all_zero": True,
        },
        "shared_e_p_correlation": {
            "sum": "e2+p2",
            "F_terms_removed_algebraically": True,
            "rows": shared_rows,
            "rate_weighted_direct_cancellation_claimed": False,
            "rate_note": "e and p retain distinct rate weights; no T_e+T_p cancellation is asserted",
        },
        "paired_history_transport": {
            "formula": "DN_out=MN*DN_in+IN; DC_out=MC*DC_in+IC; Delta_out=MN*Delta_in+(MN-MC)*DC_in+(IN-IC)",
            "checks": history_rows,
            "check_count": len(history_rows),
            "all_zero": True,
            "chart_names_are_contract_supplied": True,
            "chart_names": list(checked["chart_names"]),
            "identity_scope": "the 20 symbolic channel/order rows apply independently on each supplied chart",
            "memory_equality_assumed": False,
            "density_equality_assumed": False,
            "quiet_local_zero_does_not_reset_incoming": True,
        },
        "scope": {
            "global_gap_zero": False,
            "numeric_values_selected": False,
            "cap_reduction_admitted": False,
            "physical_completion_proved": False,
            "source_provider_calls": 0,
            "old_kernel_or_helper_replayed": False,
        },
        "audit": {
            "new_integer_sparse_polynomial_algorithm": True,
            "fraction_arithmetic_used": False,
            "source_calls": 0,
            "quadrature_calls": 0,
            "numeric_field_evaluations": 0,
        },
    }


__all__ = (
    "sdh_contract_template",
    "sdh_validate_contract",
    "derive_signed_density_defect_algebra",
)


"""Unexecuted literal source binder for signed compact-history transport.

This wrapper consumes embedded authenticated receipts only.  It does not read
files, import source producers, replay old proofs, evaluate source values, or
perform numerical integration.  The pure algebra helper is concatenated before
this file and is called once after the finite source guards pass.
"""

import ast
import hashlib
import json


SHD_INPUT_SCHEMA = "compact_signed_history_defect_transport_input_v1"
SHD_RESULT_SCHEMA = "original_compact_signed_history_defect_transport"
SHD_COMPONENT = "compact_signed_history_defect_transport"
SHD_ORDERS = (0, 1, 2, 3)
SHD_FIELDS = ("m", "h", "k", "e", "p")
SHD_N = 18
SHD_J = 1358356628656378313
SHD_ENDPOINT_NODE = 48
SHD_NATIVE_ARCHIVE_SHA = (
    "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
)
SHD_C3_ARCHIVE_SHA = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
SHD_NATIVE_RECEIPT_SHA = (
    "7fe08adc82d74a269ef0fdb4d5062f4805de6e4ed5fa1a0e0f8d71ed5cb2f520"
)
SHD_C3_RECEIPT_SHA = (
    "e42566ad2a13528b6a302c1cb1130f0cc4c1801b51eefaecd01a785267334f88"
)
SHD_CACHED_JOINT_SHA = (
    "13cad9f7578258f209d15b84a74efc0df3594049da6c151f4859453add08e3c1"
)
SHD_NATIVE_RAW = {
    "m": [258646, 258647, 258648, 258649],
    "h": [258652, 258653, 258654, 258655],
    "k": [258658, 258659, 258660, 258661],
    "e": [258664, 258665, 258666, 258667],
    "p": [258670, 258671, 258672, 258673],
}
SHD_NATIVE_SCALED = {
    "m": [258883, 258884, 258885, 258886],
    "h": [258889, 258890, 258891, 258892],
    "k": [258895, 258896, 258897, 258898],
    "e": [258901, 258902, 258903, 258904],
    "p": [258907, 258908, 258909, 258910],
}
SHD_C3_TERMINAL = {
    "m": [4920, 4921, 12042, 24177],
    "h": [4928, 4929, 12043, 24178],
    "k": [4930, 4931, 12044, 24179],
    "e": [4932, 4933, 12045, 24180],
    "p": [4888, 4889, 12041, 24176],
}
SHD_NATIVE_EXPRESSION = {
    "m": [253257, 253258, 253259, 253260],
    "h": [253290, 253291, 253292, 253293],
    "k": [253320, 253321, 253322, 253323],
    "e": [253350, 253351, 253352, 253353],
    "p": [253374, 253375, 253376, 253377],
}


def shd_require(condition, message):
    if not condition:
        raise ValueError(message)


def shd_clone(value):
    return json.loads(json.dumps(value, sort_keys=True, separators=(",", ":")))


def shd_text_sha(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def shd_journal(journal, event, **payload):
    if journal is not None:
        shd_require(callable(journal), "journal must be callable")
        journal({"event": event, **payload})


def shd_unwrap(value):
    """Unwrap an embedded receipt without path access or decoding."""
    current = value
    for _ in range(4):
        if not isinstance(current, dict):
            return current
        for key in ("payload", "receipt", "data", "record", "value"):
            child = current.get(key)
            if isinstance(child, dict):
                current = child
                break
        else:
            return current
    return current


def shd_find(spec, name, aliases=()):
    names = (name,) + tuple(aliases)
    containers = [spec]
    for key in ("receipts", "source_receipts", "source_packets"):
        child = spec.get(key)
        if isinstance(child, dict):
            containers.append(child)
    for container in containers:
        for key in names:
            value = container.get(key)
            if value is None:
                continue
            if isinstance(value, dict) and "path" in value and not any(
                k in value for k in ("payload", "receipt", "data", "record", "value")
            ):
                continue
            outer_sha = next((value.get(k) for k in
                              ("sha256", "receipt_sha256", "artifact_sha256")
                              if isinstance(value.get(k), str)), None)
            unwrapped = shd_unwrap(value)
            if outer_sha is not None and isinstance(unwrapped, dict):
                unwrapped = dict(unwrapped)
                unwrapped["_embedded_receipt_sha256"] = outer_sha
            return unwrapped
    raise ValueError("missing embedded receipt: " + name)


def shd_embedded_sha(value, expected, label):
    if not isinstance(value, dict):
        return
    for key in ("sha256", "receipt_sha256", "artifact_sha256",
                "_embedded_receipt_sha256"):
        if isinstance(value.get(key), str):
            shd_require(value[key] == expected, label + " SHA changed")
            return


def shd_row(entry):
    shd_require(isinstance(entry, dict), "source row is not a mapping")
    row = entry.get("row", entry)
    shd_require(isinstance(row, dict), "source row body is not a mapping")
    return row


def shd_port_map(value, expected, label):
    shd_require(isinstance(value, dict), label + " missing")
    for field in SHD_FIELDS:
        shd_require(value.get(field) == expected[field],
                    label + "." + field + " changed")


def shd_native_binding(receipt):
    receipt = shd_unwrap(receipt)
    shd_require(isinstance(receipt, dict), "native receipt is not a mapping")
    shd_embedded_sha(receipt, SHD_NATIVE_RECEIPT_SHA, "native receipt")
    shd_require(receipt.get("schema") == "compact_native_history_source_ports_v1",
                "native receipt schema changed")
    for key in ("static_only", "source_execution_performed",
                "source_imports_performed", "math_or_proof_replay_performed"):
        shd_require(receipt.get(key) is (True if key == "static_only" else False),
                    "native receipt scope changed: " + key)
    archive = receipt.get("authenticated_native_archive")
    shd_require(isinstance(archive, dict) and
                archive.get("raw_sha256_before_gzip_decode") == SHD_NATIVE_ARCHIVE_SHA,
                "native archive raw SHA changed")
    params = receipt.get("native_parameter_and_namespace_slice")
    shd_require(isinstance(params, dict) and params.get("N") == SHD_N and
                params.get("Pstar") == 33572 and params.get("mu") == 17,
                "native parameter slice changed")
    ports = receipt.get("native_port_maps_Z0_to_Z3")
    shd_port_map(ports.get("raw_selected_N_correction_ports"), SHD_NATIVE_RAW,
                 "native raw ports")
    shd_port_map(ports.get("N_scaled_terminal_history_ports"), SHD_NATIVE_SCALED,
                 "native N-scaled ports")
    raw_rows = receipt.get("native_raw_correction_rows_Z0_to_Z3")
    scaled_rows = receipt.get("native_N_scaled_rows_Z0_to_Z3")
    history_rows = receipt.get("native_history_rows_Z0_to_Z3")
    shd_require(all(isinstance(rows, list) and len(rows) == 20
                    for rows in (raw_rows, scaled_rows, history_rows)),
                "native source row count changed")
    raw, scaled = {}, {}
    for entry in raw_rows:
        field, order, node = entry.get("field"), entry.get("ordinary_Z_order"), entry.get("node_id")
        shd_require(field in SHD_FIELDS and order in SHD_ORDERS and
                    node == SHD_NATIVE_RAW[field][order],
                    "native raw source port changed")
        row = shd_row(entry)
        shd_require(row.get("operation") == "function_substitution" and
                    row.get("variable") == 250215 and row.get("value") == SHD_ENDPOINT_NODE and
                    row.get("Z_independent_substitution") is True and
                    row.get("expression") == SHD_NATIVE_EXPRESSION[field][order],
                    "native raw substitution guard failed")
        raw[field, order] = shd_clone(entry)
    for entry in scaled_rows:
        field, order, node = entry.get("field"), entry.get("ordinary_Z_order"), entry.get("node_id")
        shd_require(field in SHD_FIELDS and order in SHD_ORDERS and
                    node == SHD_NATIVE_SCALED[field][order],
                    "native scaled source port changed")
        row = shd_row(entry)
        shd_require(row.get("operation") == "product" and
                    row.get("arguments") == [SHD_NATIVE_RAW[field][order], SHD_N] and
                    entry.get("scaled_from") == SHD_NATIVE_RAW[field][order] and
                    entry.get("scale_factor") == SHD_N,
                    "native N scaling is not raw correction times N")
        scaled[field, order] = shd_clone(entry)
    for entry in history_rows:
        field, order = entry.get("field"), entry.get("ordinary_Z_order")
        shd_require(entry.get("node_id") == SHD_NATIVE_SCALED[field][order] and
                    shd_row(entry) == shd_row(scaled[field, order]),
                    "native history alias changed")
    return {"archive_sha256": SHD_NATIVE_ARCHIVE_SHA, "raw": raw, "scaled": scaled,
            "ports": shd_clone(ports), "params": shd_clone(params)}


def shd_c3_binding(receipt, cached_contract):
    receipt = shd_unwrap(receipt)
    shd_require(isinstance(receipt, dict), "C3 receipt is not a mapping")
    shd_embedded_sha(receipt, SHD_C3_RECEIPT_SHA, "C3 receipt")
    shd_require(receipt.get("schema") == "compact_C3_rebuilt_history_literal_addendum_v1",
                "C3 receipt schema changed")
    archives = receipt.get("archive_records", {})
    archive = archives.get("C3_repair_band_rebuilt", {}) if isinstance(archives, dict) else {}
    shd_require(archive.get("raw_sha256") == SHD_C3_ARCHIVE_SHA and
                archive.get("raw_hash_authenticated_before_decode") is True,
                "C3 archive raw authentication changed")
    distinction = receipt.get("distinction_and_nonclaims", {})
    for key in ("primary_vs_rebuilt_history_equivalence_claimed",
                "C3_to_C5_history_equivalence_claimed", "cap_or_range_value_used",
                "N_numeric_value_materialized"):
        shd_require(distinction.get(key) is False, "C3 nonclaim changed: " + key)
    ports = receipt.get("accepted_rebuilt_terminal_history_ports", {})
    for field in SHD_FIELDS:
        shd_require(ports.get(field) == SHD_C3_TERMINAL[field],
                    "C3 rebuilt terminal port changed: " + field)
    packet = receipt.get("literal_terminal_history_rows", {})
    rows = packet.get("rows") if isinstance(packet, dict) else None
    shd_require(isinstance(rows, list) and len(rows) == 20, "C3 terminal rows changed")
    by_key = {}
    for entry in rows:
        field, order, node = entry.get("field"), entry.get("ordinary_Z_order"), entry.get("node")
        shd_require(field in SHD_FIELDS and order in SHD_ORDERS and
                    node == SHD_C3_TERMINAL[field][order],
                    "C3 terminal source port changed")
        shd_require(shd_row(entry).get("operation") in ("product", "sum"),
                    "C3 terminal operation changed")
        by_key[field, order] = shd_clone(entry)
    role = receipt.get("source_recipe_C3_role_metadata", {})
    shd_require(role.get("namespace") == "repair_C3_band",
                "C3 role namespace changed")
    ordinary_roles = cached_contract.get("ordinary_Z_roles", {})
    c3_roles = ordinary_roles.get("C3", {}) if isinstance(ordinary_roles, dict) else {}
    shd_require(c3_roles.get("namespace") == "repair_C3_band" and
                c3_roles.get("terminal_history_ports_Z0_to_Z3") == SHD_C3_TERMINAL and
                ordinary_roles.get("range_values_are_not_function_values") is True,
                "cached rebuilt C3 roles changed")
    nmeta = receipt.get("N_and_endpoint_literal_metadata", {})
    shd_require(nmeta.get("C3_J_literal") == SHD_J and
                nmeta.get("native_C5_J_literal") == SHD_J,
                "native/C3 J literals differ")
    for key in ("C3_repair_node_18", "native_C5_node_18"):
        row = nmeta.get(key, {})
        shd_require(row.get("J") == SHD_J and row.get("Z_independent") is True and
                    row.get("exact_positive_integer") is True and
                    row.get("operation") == "exact_positive_integer_expression" and
                    row.get("definition") == "N=2^(2^J)" and
                    row.get("type") == "DoubleDyadicInteger" and
                    row.get("outer_integer_and_inner_power_not_materialized") is True,
                    key + " changed")
    for key in ("C3_repair_invN_node_19", "native_C5_invN_node_19"):
        row = nmeta.get(key, {})
        shd_require(row.get("operation") == "positive_quotient" and
                    row.get("numerator") == 1 and row.get("denominator") == SHD_N,
                    key + " changed")
    shd_require(nmeta.get("endpoint_value_node_from_joint_contract") == SHD_ENDPOINT_NODE and
                nmeta.get("endpoint_value_not_evaluated") is True,
                "endpoint metadata changed")
    return {"archive_sha256": SHD_C3_ARCHIVE_SHA, "rows": by_key,
            "ports": shd_clone(ports), "N_metadata": shd_clone(nmeta)}


def shd_cached_contract(cached):
    cached = shd_unwrap(cached)
    shd_require(isinstance(cached, dict), "cached joint result is not a mapping")
    cached_sha = cached.get("sha256", cached.get("_embedded_receipt_sha256"))
    if isinstance(cached_sha, str):
        shd_require(cached_sha == SHD_CACHED_JOINT_SHA,
                    "cached joint SHA changed")
    contract = cached.get("source_recipe_contract")
    if not isinstance(contract, dict) and isinstance(cached.get("result"), dict):
        contract = cached["result"].get("source_recipe_contract")
    shd_require(isinstance(contract, dict), "cached source recipe contract missing")
    ordinary = contract.get("ordinary_Z_roles", {})
    c3 = ordinary.get("C3", {}) if isinstance(ordinary, dict) else {}
    shd_require(c3.get("namespace") == "repair_C3_band" and
                c3.get("terminal_history_ports_Z0_to_Z3") == SHD_C3_TERMINAL and
                ordinary.get("range_values_are_not_function_values") is True and
                contract.get("dictionary_key_order_has_mathematical_meaning") is False,
                "cached C3 role contract changed")
    return contract


def shd_normalization_binding(addendum):
    addendum = shd_unwrap(addendum)
    shd_require(isinstance(addendum, dict), "normalization addendum missing")
    schema = addendum.get("schema")
    if schema == "compact_history_normalization_literal_addendum_v1":
        shd_require(addendum.get("static_only") is True and
                    addendum.get("source_execution") is False and
                    addendum.get("proof_or_normalization_execution") is False and
                    addendum.get("old_or_new_checker_imported") is False,
                    "history normalization addendum scope changed")
        selectors = addendum.get("archive_literal_selectors", {})
        c3 = selectors.get("C3_repair_band", {}) if isinstance(selectors, dict) else {}
        native = selectors.get("native_true_compact_C5", {}) if isinstance(selectors, dict) else {}
        c3n = c3.get("N_graph_node", {}).get("literal_row", {})
        nsem = native.get("N_semantics", {}).get("literal_row", {})
        shd_require(c3n.get("J") == SHD_J and nsem.get("J") == SHD_J and
                    c3n.get("Z_independent") is True and nsem.get("Z_independent") is True and
                    c3n.get("operation") == "exact_positive_integer_expression" and
                    nsem.get("operation") == "exact_positive_integer_expression" and
                    c3n.get("definition") == "N=2^(2^J)" and
                    nsem.get("definition") == "N=2^(2^J)" and
                    c3n.get("type") == "DoubleDyadicInteger" and
                    nsem.get("type") == "DoubleDyadicInteger" and
                    c3n.get("outer_integer_and_inner_power_not_materialized") is True and
                    nsem.get("outer_integer_and_inner_power_not_materialized") is True,
                    "history normalization N literal changed")
        contracts = addendum.get("unit_and_normalization_literal_contracts", {})
        shd_require(isinstance(contracts, dict) and
                    "native_compact_branch" in contracts and
                    "physical_pressure_lift" in contracts and
                    "Rm_generic_inputs_recover_inputs" in contracts,
                    "history normalization unit contracts missing")
        return {
            "schema": schema,
            "ordinary_Z_contract": {
                "orders": list(SHD_ORDERS),
                "orders_are_source_roles": True,
                "range_values_are_not_function_values": True,
            },
            "formula_contract": {
                "C3_scaling": "C3 terminal rows are already N-scaled by the selected source construction.",
                "C5_scaling": "C5 raw correction rows are multiplied by the graph N once.",
                "no_value_evaluation": True,
            },
            "unit_and_normalization_literal_contracts": shd_clone(contracts),
            "no_value_evaluation": True,
        }
    shd_require(schema in (
        "compact_joint_target_normalization_literal_slice_v1",
        "compact_history_normalization_literal_addendum_v1",
        "compact_signed_history_normalization_literal_addendum_v1",
    ), "normalization addendum schema changed")
    ordinary = addendum.get("ordinary_Z_contract", {})
    shd_require(ordinary.get("orders") == list(SHD_ORDERS) and
                ordinary.get("orders_are_source_roles") is True and
                ordinary.get("range_values_are_not_function_values") is True,
                "normalization ordinary-Z contract changed")
    formula = addendum.get("formula_contract", {})
    shd_require(isinstance(formula, dict) and formula.get("no_value_evaluation") is True,
                "normalization formula contract changed")
    c3_text, c5_text = str(formula.get("C3_scaling", "")), str(formula.get("C5_scaling", ""))
    shd_require("already" in c3_text.lower() and "N" in c3_text,
                "C3 terminal scaling is not declared already scaled")
    shd_require("raw correction" in c5_text.lower() and "N" in c5_text,
                "C5 raw correction N scaling is not declared")
    return {"schema": addendum.get("schema"), "ordinary_Z_contract": shd_clone(ordinary),
            "formula_contract": shd_clone(formula), "no_value_evaluation": True}


def shd_source_text_guard(spec):
    records = spec.get("creator_source_records")
    shd_require(isinstance(records, dict) and records, "creator source records missing")
    explicit = spec.get("source_text_clauses", spec.get("source_clause_contract", {}))
    shd_require(isinstance(explicit, dict), "source text clause contract malformed")
    defaults = {
        "prescribed_loop_C5_functions.py": (
            "deffinite_N_functions", "invN=g.quotient(g.one,N",
            "r=a.mul(roots['E'],exm1)", "s=a.scale(loop['B_over_Pstar'],invN)",
            "inc=dict(m=s,h=r", "k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s))",
            "e=a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s))",
            "p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))",
        ),
        "true_compact_C5_functions.py": (
            "defnormalized_targets", "scaled={k:a.scale(v,N)for k,v in correction.items()}",
            "scaled['k']", "actual_C3_N_scaled_incoming",
        ),
        "C3_source_targets.py": (
            "defcoefficient_C3(g,E,V,A,B,N,old)",
            "first=dict(m=B,h=F,k=alg.add(VF,EB),e=alg.add(alg.scale(VB,2),alg.neg(EF)),p=EF)",
            "second=dict(m=zero,h=zero,k=alg.mul(F,B),e=alg.add(B2,alg.neg(alg.scale(F2,'1/2'))),p=alg.scale(F2,'1/2'))",
            "epsilon=alg.fixed(g.quotient(g.one,N,'same exact positive integer N'))",
            "density={key:preserve(old['density'][key],alg.add(coefficients[-1][key],alg.mul(epsilon,coefficients[-2][key])).ZZZ)",
            "ifbase['quiet_local_source_zero']",
            "predecessor_ZZZ_memory_retained=True",
            "terminal_N_scaled_histories=history", "actual_N_scaled_five_targets=targets",
            "history[key]", "joint=alg.add(history['k']",
        ),
        "C3_repair_band.py": (
            "defleading_C3(field,prior,x)", "defbuild(field)",
            "signed_density_C3", "invN=g.quotient(g.one,built['N']",
            "actual_incoming_correction_C2", "ordinary_slow_Z_derivative_order=3",
            "existing_C2_profiles_histories_and_P0_retained=True",
        ),
        "modified_Q4_pressure_C3_functions.py": (
            "defphysical_pressure_fields(g,fields,radial,delta,Pstar)",
            "U=a.scale(fields['V'],Pstar)", "Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar)", "my=a.sub(U,m)",
            "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))",
        ),
        "Rm_generic_inputs.py": (
            "invS=f.factor((0,-.5,0,0,0))", "E=raw['velocity']['theta'][0]",
            "V=scale(raw['velocity']['axial'][0],invS)",
            "histories={key:scale(rows[0],invS)",
        ),
        "core_inlet_axial6_generic_source.py": (
            "Pstar=f.factor((0,.5,0,0,0))", "generic.recover_inputs(op,actual)",
        ),
    }
    checked = []
    physical_lift_guard = False
    recover_guard = False
    for name, entry in records.items():
        shd_require(isinstance(entry, dict) and isinstance(entry.get("source_text"), str),
                    "source text missing: " + str(name))
        text = entry["source_text"]
        expected = entry.get("sha256", entry.get("raw_sha256"))
        if isinstance(expected, str):
            shd_require(shd_text_sha(text) == expected, "source text SHA changed: " + str(name))
        try:
            ast.parse(text, filename=str(name), mode="exec")
        except SyntaxError as error:
            raise ValueError("source AST changed: " + str(name)) from error
        basename = str(name).replace("\\", "/").split("/")[-1]
        clauses = explicit.get(name, explicit.get(basename))
        if clauses is None:
            clauses = ()
            for suffix, candidate in defaults.items():
                if basename == suffix or basename.endswith(suffix):
                    clauses = candidate
                    break
        shd_require(isinstance(clauses, (list, tuple)),
                    "source clause list malformed: " + basename)
        compact = "".join(text.split())
        missing = [str(c) for c in clauses if str(c).replace(" ", "") not in compact]
        shd_require(not missing, "source formula clause missing: " + basename)
        if basename.endswith("modified_Q4_pressure_C3_functions.py"):
            physical_lift_guard = True
        if basename.endswith("Rm_generic_inputs.py") or basename.endswith(
            "core_inlet_axial6_generic_source.py"
        ):
            recover_guard = True
        checked.append({"basename": basename, "sha256": expected, "AST_parsed": True,
                        "clause_count": len(clauses)})
    for suffix in defaults:
        shd_require(any(row["basename"].endswith(suffix) and row["clause_count"] > 0
                        for row in checked), "required source formula family missing: " + suffix)
    shd_require(physical_lift_guard and recover_guard,
                "normalized source-unit convention clauses are missing")
    return {
        "records": checked,
        "unit_conventions_bound": bool(physical_lift_guard and recover_guard),
        "physical_lift_source_clause_bound": physical_lift_guard,
        "recover_source_clause_bound": recover_guard,
    }


def shd_gap_nodes(native, c3):
    gaps = []
    for field in SHD_FIELDS:
        for order in SHD_ORDERS:
            left = native["scaled"][field, order]
            raw = native["raw"][field, order]
            right = c3["rows"][field, order]
            gaps.append({
                "gap_index": len(gaps), "field": field, "ordinary_Z_order": order,
                "operation": "signed_difference",
                "left": {"archive_sha256": SHD_NATIVE_ARCHIVE_SHA,
                          "namespace": "true_compact_C5",
                          "role": "native_N_scaled_terminal_history",
                          "source_node": left["node_id"],
                          "raw_correction_node": raw["node_id"],
                          "raw_correction_role": "native_function_substitution_at_coordinate_2",
                          "N_factor": SHD_N, "N_scaling_applied_once": True},
                "right": {"archive_sha256": SHD_C3_ARCHIVE_SHA,
                          "namespace": "repair_C3_band",
                          "role": "C3_rebuilt_terminal_N_scaled_history",
                          "source_node": right["node"],
                          "source_role": right.get("mapped_role"),
                          "already_N_scaled": True},
                "common_N_literal_binding": {
                    "J": SHD_J, "native_node_18": 18, "C3_node_18": 18,
                    "native_invN_node_19": 19, "C3_invN_node_19": 19,
                    "N_integer_materialized": False,
                },
                "values_uncomputed": True, "gap_zero_claimed": False,
                "same_units_claimed": False, "inlet_cap_equality_claimed": False,
                "normalization": (
                    "generic common E=Utheta/Pstar, V=Uz/Pstar, "
                    "m/k scaled1/Pstar, h/e/p retained per original moment definitions"
                ),
                "conditional_common_point_pullback_unbound": True,
            })
    shd_require(len(gaps) == 20, "signed gap count changed")
    return gaps


def shd_endpoint_literal_binding(spec):
    packet = spec.get("endpoint_literal_rows")
    shd_require(isinstance(packet, dict), "endpoint literal rows are missing")
    checked = {}
    for side in ("native", "C3"):
        entry = packet.get(side)
        shd_require(isinstance(entry, dict) and entry.get("node_id") == SHD_ENDPOINT_NODE,
                    "endpoint node changed: " + side)
        row = entry.get("row", entry)
        shd_require(isinstance(row, dict) and row.get("operation") == "exact_rational" and
                    row.get("numerator") == 2 and row.get("denominator") == 1,
                    "endpoint literal row changed: " + side)
        checked[side] = shd_clone(entry)
    return checked


def shd_contract(spec, gaps, source_ready, unit_conventions_bound):
    provided = spec.get("algebra_contract")
    shd_require(isinstance(provided, dict), "algebra_contract missing")
    contract = shd_clone(provided)
    contract["schema"] = "compact_signed_history_density_algebra_input_v1"
    binding = contract.setdefault("source_binding", {})
    shd_require(isinstance(binding, dict), "algebra source binding malformed")
    binding.update({
        "root_wrapper_passed": bool(source_ready),
        "eps_definition": "eps=1/N",
        "eps_source_ref": "literal node 18 and positive-quotient node 19 in native and C3 receipts",
        "units_source_ref": (
            "source-unit conventions guarded from generic recover and physical-lift clauses; "
            "native/C3 physical Pstar equality remains unclaimed"
        ),
        "normalized_source_unit_conventions_bound": bool(unit_conventions_bound),
        "numeric_source_values_selected": False, "source_calls_allowed": False,
        "caps_admitted": False,
    })
    ports = contract.get("field_source_ports")
    shd_require(isinstance(ports, dict), "algebra field source ports missing")
    for field in ("En", "Ec", "Vn", "Vc", "Fn", "Fc", "Bn", "Bc"):
        row = ports.get(field)
        shd_require(isinstance(row, dict) and row.get("ordinary_Z_orders") == list(SHD_ORDERS),
                    "algebra field ports changed: " + field)
        shd_require(isinstance(row.get("ordinary_Z_ports"), list) and
                    len(row["ordinary_Z_ports"]) == 4,
                    "algebra field port count changed: " + field)
        recipes = {
            "En": "native.prescribed_loop_C5_functions.finite_N_functions.roots.E",
            "Vn": "native.prescribed_loop_C5_functions.finite_N_functions.roots.V",
            "Fn": "native.prescribed_loop_C5_functions.N_times_delta_E",
            "Bn": "native.prescribed_loop_C5_functions.N_times_delta_V",
            "Ec": "C3.rebuilt.coefficient_C3.E",
            "Vc": "C3.rebuilt.coefficient_C3.V",
            "Fc": "C3.rebuilt.F_N",
            "Bc": "C3.rebuilt.B_source_coefficient",
        }
        row["source_ref"] = "FUNCTION_RECIPE:" + recipes[field]
        if all(isinstance(item, str) and item.startswith("root-port-")
               for item in row["ordinary_Z_ports"]):
            row["ordinary_Z_ports"] = [
                "FUNCTION_RECIPE:" + recipes[field] + ":ordinary_Z_" + str(order) +
                ":conditional_common_point_pullback_unbound"
                for order in SHD_ORDERS
            ]
    history = contract.setdefault("history_transport", {})
    shd_require(isinstance(history, dict), "algebra history transport malformed")
    history.setdefault("chart_names", list(spec.get("chart_names", ["source_signed_history"])))
    history.setdefault("history_channels", list(SHD_FIELDS))
    history["memory_equality_assumed"] = False
    history["density_equality_assumed"] = False
    flags = contract.setdefault("flags", {})
    shd_require(isinstance(flags, dict), "algebra flags malformed")
    flags.update({"global_gap_zero": False, "numeric_values_selected": False,
                  "cap_reduction_admitted": False, "physical_completion_proved": False})
    contract["source_gap_function_nodes"] = shd_clone(gaps)
    contract["source_gap_function_node_count"] = 20
    contract["typed_function_recipe_contract"] = {
        "native": {
            "En": "background E",
            "Vn": "background V",
            "Fn": "N*delta_E",
            "Bn": "N*delta_V",
        },
        "C3": {
            "Ec": "source E",
            "Vc": "source V",
            "Fc": "F_N",
            "Bc": "B source coefficient",
        },
        "conditional_common_point_pullback_unbound": True,
        "actual_32_node_bindings_claimed": False,
    }
    contract["terminal_history_common_N_binding"] = {
        "native_archive_sha256": SHD_NATIVE_ARCHIVE_SHA,
        "C3_archive_sha256": SHD_C3_ARCHIVE_SHA, "J": SHD_J,
        "N_node": 18, "invN_node": 19, "N_integer_materialized": False,
        "same_units_claimed": False, "inlet_cap_equality_claimed": False,
    }
    return contract


def shd_parent_data(source, parent):
    retained = source.get("retained_parent_data", {})
    if not isinstance(retained, dict):
        retained = {}
    result = {}
    for key in ("retained_cap_provenance", "retained_construction_data",
                "retained_integral_term_inventory"):
        if key in retained:
            result[key] = shd_clone(retained[key])
        elif isinstance(parent, dict) and key in parent:
            result[key] = shd_clone(parent[key])
    if "coordinate_mapping" in source:
        result["coordinate_mapping"] = shd_clone(source["coordinate_mapping"])
    elif isinstance(parent, dict) and "coordinate_mapping" in parent:
        result["coordinate_mapping"] = shd_clone(parent["coordinate_mapping"])
    return result


def shd_flags(source, parent):
    result = {name: False for name in source.get("inherited_false_flags", [])
              if isinstance(name, str)}
    inherited = source.get("inherited_parent_true_flags", {})
    if isinstance(inherited, dict):
        for name, value in inherited.items():
            if value is True and name not in result and name not in {
                "passed", "source_packet_defined", "kernel_imported", "metadata_only",
            }:
                result[name] = True
    if isinstance(parent, dict) and parent.get("actual_native_C3_amplitude_identity_proved") is True:
        result["actual_native_C3_amplitude_identity_proved"] = True
    return result


def shd_failure(source, parent, error=None):
    result = {
        "schema": SHD_RESULT_SCHEMA, "component": SHD_COMPONENT, "passed": False,
        "metadata_only": False, "source_only": True, "input_digest": source.get("input_digest"),
        "runtime_model_metadata": None, "original_source_evaluations": 0,
        "source_packet_defined": False,
        "signed_history_source_gap_functions_defined": False,
        "native_N_scaling_applied_once_bound": False,
        "five_signed_density_defect_identities_proved": False,
        "paired_history_defect_transport_proved": False,
        "normalized_source_unit_conventions_bound": False,
        "actual_native_C3_amplitude_identity_proved": bool(
            isinstance(parent, dict) and parent.get("actual_native_C3_amplitude_identity_proved") is True
        ),
        "source_inspection_audit": {
            "original_source_evaluations": 0, "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0, "old_proof_functions_called": 0,
            "old_operator_calls": 0, "old_source_normalizer_calls": 0,
            "quadratures": 0, "temporal_recurrence_calls": 0,
            "whole_archive_DAG_walks": 0, "generic_archive_wide_walks": 0,
            "new_signed_history_source_rows": 0,
            "new_signed_density_order_checks": 0,
            "new_paired_transport_order_checks": 0,
            "new_energy_pressure_local_correlation_checks": 0,
        },
        "source_binding": {"passed": False, "values_uncomputed": True,
                           "gap_zero_claimed": False, "inlet_cap_equality_claimed": False},
    }
    result.update(shd_flags(source, parent))
    result.update(shd_parent_data(source, parent))
    if error is not None:
        result["source_binding_error"] = str(error)
    return result


def derive(source, journal=None):
    """Bind 20 source differences and call the pure algebra helper once."""
    parent = source.get("parent_result", {}) if isinstance(source, dict) else {}
    try:
        shd_require(isinstance(source, dict), "source must be a mapping")
        spec = source.get("input")
        shd_require(isinstance(spec, dict) and spec.get("schema") == SHD_INPUT_SCHEMA and
                    spec.get("component") == SHD_COMPONENT,
                    "signed-history input schema changed")
        shd_require(isinstance(parent, dict) and parent.get("passed") is True and
                    parent.get("actual_native_C3_amplitude_identity_proved") is True,
                    "authenticated parent amplitude result not retained")
        endpoint_rows = shd_endpoint_literal_binding(spec)
        native = shd_native_binding(shd_find(spec, "native_history_receipt"))
        cached = shd_cached_contract(shd_find(spec, "cached_joint_target_result"))
        c3 = shd_c3_binding(shd_find(spec, "C3_rebuilt_history_receipt",
                                     ("c3_rebuilt_history_receipt",)), cached)
        normalization = shd_normalization_binding(
            shd_find(spec, "normalization_literal_addendum", ("normalization_slice",)))
        text_packet = shd_source_text_guard(spec)
        text_checks = text_packet["records"]
        unit_conventions_bound = text_packet["unit_conventions_bound"]
        gaps = shd_gap_nodes(native, c3)
        contract = shd_contract(spec, gaps, unit_conventions_bound,
                                unit_conventions_bound)
        shd_journal(journal, "signed_history_source_rows_bound",
                    native_rows=20, c3_rows=20, gap_nodes=20,
                    N_scaling_once=True, values_uncomputed=True)
        algebra_result = None
        algebra_error = None
        try:
            algebra_result = derive_signed_density_defect_algebra(contract, journal)
        except Exception as error:
            algebra_error = error
        algebra_ok = isinstance(algebra_result, dict) and algebra_result.get("passed") is True
        density_ok = bool(algebra_ok and
                          algebra_result.get("residual_proof", {}).get("all_zero") is True)
        paired_ok = bool(algebra_ok and
                         algebra_result.get("paired_history_transport", {}).get("all_zero") is True)
        passed = bool(density_ok and paired_ok)
        audit = {
            "original_source_evaluations": 0, "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0, "old_proof_functions_called": 0,
            "old_operator_calls": 0, "old_source_normalizer_calls": 0,
            "quadratures": 0, "temporal_recurrence_calls": 0,
            "whole_archive_DAG_walks": 0, "generic_archive_wide_walks": 0,
            "new_signed_history_source_rows": 20,
            "new_signed_density_order_checks": 20 if algebra_ok else 0,
            "new_paired_transport_order_checks": 20 if algebra_ok else 0,
            "new_energy_pressure_local_correlation_checks": 4 if algebra_ok else 0,
            "native_raw_correction_rows": 20, "native_N_scaled_rows": 20,
            "C3_rebuilt_terminal_rows": 20, "N_literal_descriptor_checks": 4,
            "creator_source_AST_checks": len(text_checks),
            "source_formula_clause_checks": sum(x["clause_count"] for x in text_checks),
        }
        result = {
            "schema": SHD_RESULT_SCHEMA, "component": SHD_COMPONENT, "passed": passed,
            "metadata_only": False, "source_only": True, "input_digest": source.get("input_digest"),
            "runtime_model_metadata": None, "original_source_evaluations": 0,
            "source_packet_defined": unit_conventions_bound,
            "signed_history_source_gap_functions_defined": passed,
            "native_N_scaling_applied_once_bound": passed,
            "five_signed_density_defect_identities_proved": density_ok,
            "paired_history_defect_transport_proved": paired_ok,
            "actual_native_C3_amplitude_identity_proved": True,
            "normalized_source_unit_conventions_bound": unit_conventions_bound,
            "endpoint_literal_rows_bound": endpoint_rows,
            "source_binding": {
                "passed": passed, "source_gap_functions_defined": True,
                "signed_gap_node_count": 20,
                "native_archive_sha256": SHD_NATIVE_ARCHIVE_SHA,
                "C3_archive_sha256": SHD_C3_ARCHIVE_SHA,
                "native_raw_correction_row_count": 20,
                "native_N_scaled_row_count": 20,
                "C3_rebuilt_terminal_row_count": 20,
                "N_scaling_applied_once": True,
                "C3_terminal_rows_already_N_scaled": True,
                "normalized_source_unit_conventions_bound": unit_conventions_bound,
                "values_uncomputed": True, "numeric_field_values_selected": False,
                "cap_or_range_value_used": False, "gap_zero_claimed": False,
                "inlet_cap_equality_claimed": False, "physical_Pstar_lift_claimed": False,
                "same_units_claimed": False, "primary_C3_namespace_compared": False,
                "common_N_literal_binding": {
                    "J": SHD_J, "native_node_18": 18, "C3_node_18": 18,
                    "native_invN_node_19": 19, "C3_invN_node_19": 19,
                    "N_integer_materialized": False,
                },
                "endpoint_literal_rows_bound": endpoint_rows,
            },
            "native_history_source_ports": {
                "archive_sha256": SHD_NATIVE_ARCHIVE_SHA,
                "raw_correction_ports": shd_clone(SHD_NATIVE_RAW),
                "N_scaled_ports": shd_clone(SHD_NATIVE_SCALED),
                "ordinary_Z_orders": list(SHD_ORDERS),
                "function_substitution_coordinate": {
                    "variable_node": 250215, "value_node": SHD_ENDPOINT_NODE,
                    "coordinate_value": 2, "Z_independent_substitution": True,
                },
            },
            "C3_rebuilt_history_source_ports": {
                "archive_sha256": SHD_C3_ARCHIVE_SHA,
                "terminal_ports": shd_clone(SHD_C3_TERMINAL),
                "ordinary_Z_orders": list(SHD_ORDERS), "namespace": "repair_C3_band",
                "range_values_not_function_values": True,
            },
            "signed_history_defect_transport_certificate": {
                "passed": passed,
                "scope": "20 archive-keyed terminal history gap definitions, normalized source conventions, and conditional local density/paired-history algebra through ordinary Z3",
                "source_history_gap_rows": 20,
                "density_identity_checks": 20 if density_ok else 0,
                "paired_transport_checks": 20 if paired_ok else 0,
                "local_energy_pressure_correlation_checks": 4 if algebra_ok else 0,
                "algebra_result_selector": "$.signed_density_algebra_result",
                "actual_common_point_pullback_bound": False,
                "actual_17_kernel_pairs_bound": False,
                "terminal_gap_zero_proved": False,
            },
            "signed_history_defect_transport_consequence": {
                "function_identity": "deltaT_i^(j)(Z)=T_N,i^(j)(Z)-T_C,i^(j)(Z), i=m,h,k,e,p; j=0..3",
                "formula": "Delta_out=M_N*Delta_in+(M_N-M_C)*T_C_in+(I_N-I_C)",
                "local_density_formula_selector": "$.signed_density_algebra_result.factorized_signed_defects.unscaled_formulas",
                "native_N_scaling_once": True,
                "C3_terminal_already_N_scaled": True,
                "unequal_memory_term_retained": True,
                "predecessor_incoming_retained": True,
                "local_e_plus_p_correlation_only": True,
                "integrated_energy_pressure_cancellation_claimed": False,
                "numeric_gap_values_computed": False,
                "terminal_gap_zero_proved": False,
            },
            "signed_history_gap_function_nodes": gaps,
            "terminal_history_common_N_binding": {
                "native_archive_sha256": SHD_NATIVE_ARCHIVE_SHA,
                "C3_archive_sha256": SHD_C3_ARCHIVE_SHA, "J": SHD_J,
                "N_node": 18, "invN_node": 19, "N_integer_materialized": False,
                "same_units_claimed": False, "inlet_cap_equality_claimed": False,
                "endpoint_literal_rows_bound": endpoint_rows,
            },
            "source_recipe_contract": cached,
            "normalization_literal_binding": normalization,
            "creator_source_guards": text_checks,
            "signed_density_algebra_contract": contract,
            "signed_density_algebra_result": algebra_result,
            "source_inspection_audit": audit,
            "paired_history_transport_scope": {
                "source_independent_algebra": True,
                "actual_native_C3_kernel_pair_bound": False,
                "whole_17_actual_paired_kernels_matched": False,
                "zero_claimed": False, "common_point_domain_bound": False,
                "root_difference_transport_next": True,
            },
            "cached_joint_target_reference": {
                "sha256": SHD_CACHED_JOINT_SHA, "data_only": True,
                "old_polynomial_results_mapped_to_new_source_ports": False,
                "source_recipe_roles_match": True,
            },
            "conditional_scope": [
                "20 archive-keyed signed source gap functions are defined; values remain uncomputed.",
                "Node 18 and positive-quotient node 19 bind one common N literal; terminal decay rows are not N factors.",
                "The pure density and generic paired-history identities are source-independent.",
                "Physical units, inlet-cap transfer, whole-17 matching, and zero gaps remain unclaimed.",
            ],
            "remaining_scope": [
                "Bind the common point domain before quantitative or physical claims.",
                "Do not map primary C3 source-target rows to rebuilt C3 rows.",
                "No physical, temporal, global, cap, numeric, or whole-17 completion is asserted.",
            ],
        }
        result.update(shd_flags(source, parent))
        result.update(shd_parent_data(source, parent))
        if algebra_error is not None:
            result["algebra_error"] = str(algebra_error)
        shd_journal(journal, "signed_history_defect_transport_defined",
                    passed=passed, algebra_passed=algebra_ok,
                    density_checks=audit["new_signed_density_order_checks"],
                    paired_checks=audit["new_paired_transport_order_checks"],
                    source_values_uncomputed=True)
        return result
    except Exception as error:
        return shd_failure(source if isinstance(source, dict) else {},
                           parent if isinstance(parent, dict) else {}, error)


__all__ = ("derive",)
