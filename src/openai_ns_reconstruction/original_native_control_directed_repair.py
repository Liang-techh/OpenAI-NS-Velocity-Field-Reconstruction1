"""Directed local repairs for the affected native-control logarithmic caps.

This module is a narrow corrective kernel. It reads the accepted integrity
packet and its preserved signed interval rows, then rebuilds three local
families of inputs:

* direct abs_log conversion checkpoints for signed U/M leaves;
* the signed L=1-delta*Z**2 five-row jet and its directed reciprocal.
* the exact polynomial reciprocal Utheta on restoration/post_restore.

No accepted native-control builder, producer, checker, owner, recurrence,
quadrature, or self-check is imported or replayed. The returned packet keeps
the inherited historical caps as evidence, but does not use those caps as
signed values and does not certify downstream native, mixed-loop, Rm, or
finite-N results.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_directed_abs_log import abs_log_upper
from .original_directed_log_reciprocal import (
    directed_log_reciprocal, fixture_bound_check as reciprocal_fixture_check,
)


BASE = "experiments/root_st073/"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
INTEGRITY_REPORT = (
    BASE + "lei_ren_part1_paper_compliant_current_original_abs_log_integrity.json"
)
INTEGRITY_REPORT_SHA256 = (
    "dd4cb5741cd7e04c112d2968bf9a81aaceecfac3351a7ca93475c5c9aae646dd"
)
CONSTRUCTION_PARENT = (
    BASE + "lei_ren_part1_paper_compliant_current_original_restoration_finite_N_I1.json"
)
CONSTRUCTION_PARENT_SHA256 = (
    "fbede8e8c948a2025c65671bce6265af84e57fcb25d4a7e497bf456769d647e4"
)

MODULE = "src/openai_ns_reconstruction/original_native_control_directed_repair.py"
DIRECTED_PRIMITIVE = "src/openai_ns_reconstruction/original_directed_abs_log.py"
DIRECTED_PRIMITIVE_SHA256 = (
    "96126cf50406cfc87c8d9af854a42ac9b4ed9a0222b55587d5c913bfcd9ec767"
)
RECIPROCAL_HELPER = (
    "src/openai_ns_reconstruction/original_directed_log_reciprocal.py"
)
RECIPROCAL_HELPER_SHA256 = (
    "5f948b11559c0c7533f3b6f98a210a18b06f3a959ae96be31d3ef139fcad72df"
)
LONG_NATIVE = "src/openai_ns_reconstruction/original_long_reference_native_controls.py"
LONG_NATIVE_SHA256 = (
    "b15f020def2316b4e91a38dea0faa6b6a04cc35fa50f1cf6f0b78fa4d2235094"
)
RESTORATION_NATIVE = (
    "src/openai_ns_reconstruction/original_restoration_native_controls.py"
)
RESTORATION_NATIVE_SHA256 = (
    "9682dd93ae1637467cdc068542d4accdd03aab3e413cdc75f6fd14019b8ad69d"
)

ORDERS = tuple(range(5))
CHARTS = (
    ("long_reference", "long_reshape"),
    ("long_reference", "reference"),
    ("restoration", "restoration"),
    ("restoration", "post_restore"),
)
Z_SOURCE_FORMULA = 'z = c.mpf(["37/100", "93/250"])'
L_FORMULA = (
    "L(Z)=1-delta*Z^2; L_Z=-2*delta*Z; "
    "L_ZZ=-2*delta; L_Z^k=0 for k>=3"
)


def _lower(value):
    return make_mpf(value._mpi_[0])


def _upper(value):
    return make_mpf(value._mpi_[1])


def _point(context, value):
    return context.mpf([value, value])


def _record_or_none(value):
    return None if value is None else record(value)


def _path_get(root, dotted):
    value = root
    for key in dotted.split("."):
        require(isinstance(value, Mapping) and key in value,
                "Missing preserved source handle: " + dotted)
        value = value[key]
    return value


def _read_signed_row(context, row, label):
    if isinstance(row, Mapping) and "value" in row:
        row = row["value"]
    require(isinstance(row, Mapping), label + " must be a directed interval record")
    require(
        "lower_exact_mpf_tuple" in row
        and "upper_exact_mpf_tuple" in row
        and row.get("directed") is True,
        label + " is not a directed signed interval record",
    )
    value = read(context, row)
    require(_lower(value) <= _upper(value), label + " has reversed endpoints")
    return value


def _read_rows(context, rows, label):
    require(isinstance(rows, list), label + " must be a list")
    found = {}
    for position, row in enumerate(rows):
        require(isinstance(row, Mapping), label + " row must be a mapping")
        order = row.get("ordinary_Z_order", position)
        require(type(order) is int and order in ORDERS,
                label + " has an invalid ordinary Z order")
        require(order not in found, label + " repeats an ordinary Z order")
        found[order] = _read_signed_row(context, row, label + "[" + str(order) + "]")
    require(set(found) == set(ORDERS), label + " must contain five ordered rows")
    return [found[order] for order in ORDERS]


def _safe_logs(context, values):
    """Convert signed intervals only through the tuple-safe primitive."""
    return [abs_log_upper(context, value) for value in values]


def _directed_log_reciprocal(context, signed_rows, log_rows):
    """Apply the new tuple-safe five-row reciprocal helper to signed L rows."""
    require(len(signed_rows) == len(log_rows) == len(ORDERS),
            "The reciprocal requires five signed and logarithmic rows")
    require(_lower(signed_rows[0]) > 0,
            "The reciprocal source L must have a positive lower endpoint")
    require(log_rows[0] is not None, "The reciprocal value row cannot be zero")
    lower_value = _point(context, _lower(signed_rows[0]))
    lower_log = context.ln(lower_value)
    inverse = directed_log_reciprocal(context, log_rows, lower_log)
    return inverse, dict(
        passed=True,
        formula="directed_log_reciprocal(logjet5, certifiedlowerlog)",
        ordinary_Z_orders=list(ORDERS),
        lower_L_value=record(lower_value),
        lower_L_log=record(lower_log),
        helper_module=(
            "src/openai_ns_reconstruction/"
            "original_directed_log_reciprocal.py"
        ),
        helper_function="directed_log_reciprocal",
        interval_negation="helper reboxes and flips the lower endpoint tuple",
        scalar_endpoint_negation=False,
        signed_value_selected=False,
        source_executed=False,
    )


def _source_line_guard(root, path, snippets):
    source_path = root / path
    text = source_path.read_text(encoding="utf-8")
    compact_lines = ["".join(line.split()) for line in text.splitlines()]
    found = []
    for name, snippet in snippets:
        needle = "".join(snippet.split())
        matches = [
            index + 1
            for index, line in enumerate(compact_lines)
            if needle in line
        ]
        require(len(matches) == 1, "Source guard changed: " + path + ":" + name)
        found.append(dict(name=name, line=matches[0], required_text=snippet))
    return dict(
        path=path,
        lines=found,
        source_sha256=sha(source_path),
        source_executed=False,
    )


def _source_guards(root, parent):
    primitive_path = root / DIRECTED_PRIMITIVE
    primitive_text = primitive_path.read_text(encoding="utf-8")
    primitive_tree = ast.parse(primitive_text, filename=DIRECTED_PRIMITIVE)
    scalar_abs_calls = [
        node for node in ast.walk(primitive_tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "abs"
    ]
    require(not scalar_abs_calls,
            "The directed primitive contains a scalar abs call")

    long_section = parent["actual_long_reference_native_controls"]
    restoration_section = parent["actual_restoration_native_controls"]
    require(long_section.get("source_executed") is not True
            and restoration_section.get("source_executed") is not True,
            "An inherited native packet claims source execution")
    require(long_section.get("accepted_ancestors_replayed") is not True
            and restoration_section.get("accepted_ancestors_replayed") is not True,
            "An inherited native packet claims ancestor replay")
    for section, path, expected in (
        (long_section, LONG_NATIVE, LONG_NATIVE_SHA256),
        (restoration_section, RESTORATION_NATIVE, RESTORATION_NATIVE_SHA256),
    ):
        pin = section.get("source_pins", {}).get(path)
        require(isinstance(pin, Mapping)
                and pin.get("expected_sha256") == expected,
                "The accepted source pin is missing: " + path)
        require(sha(root / path) == expected, "Pinned native source changed: " + path)

    require(sha(primitive_path) == DIRECTED_PRIMITIVE_SHA256,
            "Pinned directed primitive changed")
    require(sha(root / RECIPROCAL_HELPER) == RECIPROCAL_HELPER_SHA256,
            "Pinned directed reciprocal helper changed")
    return dict(
        primitive=dict(
            path=DIRECTED_PRIMITIVE,
            expected_sha256=DIRECTED_PRIMITIVE_SHA256,
            actual_sha256=sha(primitive_path),
            abs_upper_definition=True,
            abs_log_upper_definition=True,
            scalar_abs_call_count=0,
            source_executed=False,
        ),
        reciprocal_helper=dict(
            path=RECIPROCAL_HELPER,
            expected_sha256=RECIPROCAL_HELPER_SHA256,
            actual_sha256=sha(root / RECIPROCAL_HELPER),
            function="directed_log_reciprocal",
            five_rows=True,
            source_executed=False,
        ),
        long_reference=_source_line_guard(
            root,
            LONG_NATIVE,
            (
                ("z_interval", Z_SOURCE_FORMULA),
                ("L_formula",
                 "Lv = [1-delta*z*z, -2*delta*z, -2*delta, c.mpf(0)]"),
            ),
        ),
        restoration=_source_line_guard(
            root,
            RESTORATION_NATIVE,
            (
                ("z_interval", 'z = context.mpf(["37/100", "93/250"])'),
                ("L_formula",
                 "L_value = [1 - delta * z * z, -2 * delta * z, -2 * delta,"),
            ),
        ),
        reciprocal_Utheta=source_guard(root, C6_SOURCE, "__init__", (
            "self.logq=logq=a.log(a.add(a.fixed(1),a.mul(z,z))",
            "t=g.symbol('true_prescribed_restoration_t');offset=g.add(c(-8),t)",
            "p=g.symbol('true_prescribed_post_restore_fraction');offset=g.add(c(-7),p)",
            "logUtheta=a.add(a.neg(logq),a.fixed(g.add(self.core.parameters['logP'],g.mul(c('1/10'),offset))))",
            "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
            "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
        )),
        inherited_packets=dict(
            long_reference_native=True,
            restoration_native=True,
            accepted_builders_or_checkers_called=False,
            source_executed=False,
        ),
        source_executed=False,
    )


def _pin(root, path, expected, *, self_pin=False):
    actual = sha(root / path)
    require(actual == expected, "Pinned source changed: " + path)
    result = dict(
        path=path,
        expected_sha256=expected,
        actual_sha256=actual,
        passed=True,
        source_executed=False,
    )
    if self_pin:
        result["self_pin"] = True
    return result


def _validate_parent(root, parent):
    report_path = root / INTEGRITY_REPORT
    require(report_path.is_file(), "The accepted integrity report is missing")
    report_bytes = report_path.read_bytes()
    require(hashlib.sha256(report_bytes).hexdigest() == INTEGRITY_REPORT_SHA256,
            "The accepted integrity report hash changed")
    disk_report = json.loads(report_bytes)
    require(dict(parent) == disk_report,
            "build requires the exact accepted integrity report mapping")
    require(parent.get("parent_report") in
            (CONSTRUCTION_PARENT, Path(CONSTRUCTION_PARENT).name)
            and parent.get("parent_report_sha256") == CONSTRUCTION_PARENT_SHA256,
            "The accepted construction parent identity changed")
    construction_path = root / CONSTRUCTION_PARENT
    require(construction_path.is_file(), "The accepted construction parent is missing")
    require(sha(construction_path) == CONSTRUCTION_PARENT_SHA256,
            "The accepted construction parent hash changed")
    require(isinstance(parent.get("source_family"), Mapping),
            "The accepted source family is missing")
    for name in ("actual_long_reference_native_controls",
                 "actual_restoration_native_controls"):
        section = parent.get(name)
        require(isinstance(section, Mapping) and section.get("passed") is True,
                "The accepted native section is missing: " + name)
        require(section.get("source_executed") is not True,
                "The accepted native section claims source execution: " + name)
    return dict(
        passed=True,
        report=INTEGRITY_REPORT,
        report_sha256=INTEGRITY_REPORT_SHA256,
        construction_parent=CONSTRUCTION_PARENT,
        construction_parent_sha256=CONSTRUCTION_PARENT_SHA256,
        source_family=copy.deepcopy(parent["source_family"]),
        exact_parent_mapping_checked=True,
        historical_native_packets_retained=True,
        source_executed=False,
    )


def _source_z(context):
    return context.mpf(["37/100", "93/250"])


def _l_rows(context, parent):
    parameter_path = "actual_core_signed_direction.parameters.delta"
    delta = read(context, _path_get(parent, parameter_path))
    z = _source_z(context)
    require(_lower(delta) > 0, "The preserved delta interval is not positive")
    require(_upper(delta) < 1, "The preserved delta interval left the source domain")
    one = context.mpf(1)
    two = context.mpf(2)
    zero = context.mpf(0)
    signed = [
        one - delta * z * z,
        -two * delta * z,
        -two * delta,
        zero,
        zero,
    ]
    require(_lower(signed[0]) > 0,
            "The directed L denominator lost source positivity")
    logs = _safe_logs(context, signed)
    inverse, reciprocal_receipt = _directed_log_reciprocal(
        context, signed, logs
    )
    return dict(
        formula=L_FORMULA,
        ordinary_Z_orders=list(ORDERS),
        z_source_handle=(
            LONG_NATIVE + ":build.z / " + RESTORATION_NATIVE + ":build.z"
        ),
        z_signed_interval=record(z),
        delta_source_handle=parameter_path,
        delta_signed_interval=record(delta),
        signed_L_rows=[
            dict(
                ordinary_Z_order=order,
                signed_interval=record(value),
                exact_zero=_lower(value) == 0 and _upper(value) == 0,
            )
            for order, value in zip(ORDERS, signed)
        ],
        directed_abs_log_upper_rows=[
            dict(
                ordinary_Z_order=order,
                exact_zero=log_value is None,
                log_absolute_upper_bound=_record_or_none(log_value),
            )
            for order, log_value in zip(ORDERS, logs)
        ],
        directed_inverse_L_rows=[
            dict(
                ordinary_Z_order=order,
                exact_zero=log_value is None,
                log_absolute_upper_bound=_record_or_none(log_value),
            )
            for order, log_value in zip(ORDERS, inverse)
        ],
        reciprocal_receipt=reciprocal_receipt,
        safe_primitive=DIRECTED_PRIMITIVE + ":abs_log_upper",
        signed_value_selected=False,
        inherited_log_caps_used_as_signed_inputs=False,
        source_executed=False,
    )


def _source_rows(parent, path):
    return _path_get(parent, path)


def _leaf_checkpoints(context, parent, path, label):
    values = _read_rows(context, _source_rows(parent, path), label)
    logs = _safe_logs(context, values)
    return dict(
        source_handle=path,
        label=label,
        ordinary_Z_orders=list(ORDERS),
        rows=[
            dict(
                ordinary_Z_order=order,
                signed_interval=record(value),
                directed_abs_log_upper=_record_or_none(log_value),
                exact_zero=log_value is None,
            )
            for order, (value, log_value) in enumerate(zip(values, logs))
        ],
        source_kind="preserved_signed_interval_data",
        signed_value_selected=False,
        source_executed=False,
    )


def _derived_signed_checkpoints(context, values, handle_prefix):
    logs = _safe_logs(context, values)
    return dict(
        source_handle=handle_prefix,
        ordinary_Z_orders=list(ORDERS),
        rows=[
            dict(
                ordinary_Z_order=order,
                signed_interval=record(value),
                directed_abs_log_upper=_record_or_none(log_value),
                exact_zero=log_value is None,
            )
            for order, (value, log_value) in enumerate(zip(values, logs))
        ],
        source_kind="new_signed_interval_algebra",
        signed_value_selected=False,
        source_executed=False,
    )


def _historical_inverse_cap(chart, section_name, chart_name):
    source = chart.get("p2_source", {})
    cap = source.get("inverse_L")
    require(isinstance(cap, list),
            "The inherited inverse-L cap is missing: " + chart_name)
    return dict(
        source_handle=section_name + ".charts." + chart_name
                        + ".p2_source.inverse_L",
        rows=copy.deepcopy(cap),
        retained_as_historical_evidence=True,
        used_as_signed_input=False,
        directed_revalidated=False,
        downstream_certification=False,
    )


def _long_chart(context, parent, chart_name, l_block):
    section_name = "actual_long_reference_native_controls"
    chart = parent[section_name]["charts"][chart_name]
    require(chart.get("passed") is True, "The inherited long chart is not passed")
    u_path = (
        "actual_switch_post_background_I1_numeric.charts.post_switch."
        "native_fields.U"
    )
    m_path = (
        "actual_long_reference_background_I1_numeric.charts."
        + chart_name + ".whole_segment_M_hull"
    )
    return dict(
        chart=chart_name,
        family="long_reference",
        ordinary_Z_orders=list(ORDERS),
        directed_L_block=copy.deepcopy(l_block),
        signed_leaf_checkpoints=dict(
            U=_leaf_checkpoints(context, parent, u_path, "long U_R110"),
            M=_leaf_checkpoints(context, parent, m_path,
                                "long whole-segment M hull"),
        ),
        historical_inverse_L_cap=_historical_inverse_cap(
            chart, section_name, chart_name
        ),
        formulas=dict(
            source_L="native_long_reference._native_p2",
            signed_U_source=u_path,
            signed_M_source=m_path,
            direct_abs_replacement="abs_log_upper(context, signed_interval)",
        ),
        status="local_directed_input_and_inverse_L_repair_only",
        downstream_recertification=False,
        source_executed=False,
    )


def _restoration_chart(context, parent, chart_name, l_block):
    section_name = "actual_restoration_native_controls"
    chart = parent[section_name]["charts"][chart_name]
    require(chart.get("passed") is True, "The inherited restoration chart is not passed")
    u_path = (
        "actual_switch_post_background_I1_numeric.charts.post_switch."
        "native_fields.U"
    )
    m_path = (
        "actual_long_reference_background_I1_numeric.charts.reference.endpoint_M"
    )
    u_rows = _read_rows(context, _source_rows(parent, u_path), "restoration U_R110")
    m_rows = _read_rows(context, _source_rows(parent, m_path), "restoration M_Rz")
    z = _source_z(context)
    four_z = [4 * z, context.mpf(4), context.mpf(0),
              context.mpf(0), context.mpf(0)]
    E = [u - base for u, base in zip(u_rows, four_z)]
    K = [m - u for m, u in zip(m_rows, u_rows)]
    return dict(
        chart=chart_name,
        family="restoration",
        ordinary_Z_orders=list(ORDERS),
        directed_L_block=copy.deepcopy(l_block),
        directed_inverse_Utheta_block=_inverse_utheta(context, parent, chart_name),
        signed_leaf_checkpoints=dict(
            U_R110=_leaf_checkpoints(context, parent, u_path, "restoration U_R110"),
            M_Rz=_leaf_checkpoints(context, parent, m_path, "restoration M_Rz"),
            four_Z=_derived_signed_checkpoints(context, four_z,
                                                "new restoration 4Z"),
            E=_derived_signed_checkpoints(context, E, "new restoration E=U_R110-4Z"),
            K=_derived_signed_checkpoints(context, K, "new restoration K=M_Rz-U_R110"),
        ),
        inherited_cap_inputs=dict(
            post_restore_M_restore_source=(
                "actual_restoration_native_controls.charts.restoration."
                "source_detail.M_memory"
            ),
            post_restore_M_restore_is_inherited_cap=True,
            inherited_cap_used_as_signed_input=False,
        ),
        historical_inverse_L_cap=_historical_inverse_cap(
            chart, section_name, chart_name
        ),
        formulas=dict(
            source_L="restoration_native._native_p2",
            signed_U_source=u_path,
            signed_M_source=m_path,
            direct_abs_replacement="abs_log_upper(context, signed_interval)",
        ),
        status="local_directed_input_and_inverse_L_repair_only",
        downstream_recertification=False,
        source_executed=False,
    )


def _inverse_utheta(context, parent, chart_name):
    """Use 1/Utheta=(1+Z²)*exp(-logP-offset/10) before differentiating.

    The factor independent of Z is bounded logarithmically on the full chart.
    Consequently Z3/Z4 vanish exactly. No inherited relative derivative cap
    or rounded scalar endpoint negation enters this replacement.
    """
    require(chart_name in ("restoration", "post_restore"), "Wrong Utheta chart")
    require(parent["actual_restoration_native_controls"]["source_function_identities"]
            ["Utheta_definition"] == "Utheta=exp(logP+offset/10-log(1+Z^2))",
            "Preserved native Utheta identity changed")
    z = _source_z(context)
    log_p_path = "actual_long_reference_background_I1_numeric.scalar_graph.logP"
    log_p = _read_signed_row(context, _path_get(parent, log_p_path), "source logP")
    offset = context.mpf([-8, -7] if chart_name == "restoration" else [-7, -6])
    log_amplitude = log_p + offset / 10
    # This subtraction acts on an interval; extracting the endpoint comes last.
    inverse_log_amplitude = -_point(context, _lower(log_amplitude))
    q_rows = [1 + z*z, 2*z, context.mpf(2), context.mpf(0), context.mpf(0)]
    q_logs = _safe_logs(context, q_rows)
    inverse = [None if row is None else row + inverse_log_amplitude for row in q_logs]
    old = parent["actual_restoration_native_controls"]["charts"][chart_name]
    return dict(
        passed=True,
        formula="1/Utheta=(1+Z^2)*exp(-logP-offset/10)",
        offset_domain=["-8", "-7"] if chart_name == "restoration" else ["-7", "-6"],
        ordinary_Z_orders=list(ORDERS),
        logP_source_handle=log_p_path,
        logP_signed_interval=record(log_p),
        offset_signed_interval=record(offset),
        log_amplitude_interval=record(log_amplitude),
        inverse_log_amplitude_upper=record(inverse_log_amplitude),
        polynomial_signed_rows=[record(row) for row in q_rows],
        polynomial_abs_log_rows=[_record_or_none(row) for row in q_logs],
        directed_inverse_Utheta_rows=[dict(
            ordinary_Z_order=k, exact_zero=row is None,
            log_absolute_upper_bound=_record_or_none(row)) for k, row in enumerate(inverse)],
        historical_inverse_Utheta_caps=copy.deepcopy(old["p2_source"]["inverse_Utheta"]),
        exact_zero_orders=[3, 4],
        whole_chart_ordinary_Z_only=True,
        physical_y_derivatives_claimed=False,
        inherited_relative_caps_used=False,
        stored_logarithms_exponentiated=False,
        scalar_endpoint_negation=False,
        signed_function_value_selected=False,
        complete_p2_recertified=False,
    )


def _fixture_bound_check():
    """Exercise a nonzero safe log and reciprocal block independently."""
    from fractions import Fraction
    from mpmath import mp

    previous = mp.prec
    tested = 0
    nonzero = 0
    try:
        for scalar_precision in (37, 53, 113):
            mp.prec = scalar_precision
            context = MPIntervalContext()
            context.prec = 1600
            for fraction in (
                Fraction(1, 3),
                Fraction(7, 10),
                Fraction(-7, 10),
                Fraction(4, 5),
                Fraction(10**100, 3),
                Fraction(1, 10**100),
            ):
                value = (
                    context.mpf(fraction.numerator)
                    / fraction.denominator
                )
                cap = abs_log_upper(context, value)
                require(cap is not None, "A nonzero fixture log became exact zero")
                from .original_directed_abs_log import abs_upper

                magnitude = abs_upper(context, value)
                reference = context.ln(magnitude)
                require(_upper(cap) >= _upper(reference),
                        "The safe log fixture underbounded its reference")
                tested += 1
                nonzero += 1

            crossing = context.mpf(
                [-context.mpf(7) / 10, context.mpf(1) / 3]
            )
            crossing_cap = abs_log_upper(context, crossing)
            require(crossing_cap is not None,
                    "The crossing interval log became exact zero")
            require(
                _upper(crossing_cap)
                >= _upper(abs_log_upper(context, context.mpf("-7/10"))),
                "The signed crossing interval was underbounded",
            )

            z = context.mpf("37/100")
            delta = context.mpf("1/1000")
            signed = [
                1 - delta * z * z,
                -2 * delta * z,
                -2 * delta,
                context.mpf(0),
                context.mpf(0),
            ]
            logs = _safe_logs(context, signed)
            inverse, receipt = _directed_log_reciprocal(
                context, signed, logs
            )
            require(all(value is not None for value in inverse),
                    "The nonzero reciprocal fixture lost a derivative row")
            require(receipt["scalar_endpoint_negation"] is False,
                    "The reciprocal fixture used scalar endpoint negation")
            require(_lower(signed[0]) > 0 and inverse[0] is not None,
                    "The reciprocal fixture lost positivity")
            nonzero += sum(value is not None for value in logs)
    finally:
        mp.prec = previous
    require(mp.prec == previous, "Global scalar precision was not restored")
    return dict(
        passed=True,
        rational_cases=tested,
        scalar_precisions=[37, 53, 113],
        interval_precision_bits=1600,
        signed_crossing_checked=True,
        nonzero_reciprocal_rows=5,
        nonzero_check_count=nonzero,
        scalar_abs_called=False,
        scalar_endpoint_negation=False,
        signed_function_value_selected=False,
        accepted_self_checks_replayed=False,
        source_executed=False,
    )


def fixture_bound_check():
    """Public local fixture; it does not read or execute an accepted packet."""
    return _fixture_bound_check()


def build(root, parent):
    """Build local directed checkpoints over the accepted integrity packet."""
    root = Path(root)
    parent_scope = _validate_parent(root, parent)
    guards = _source_guards(root, parent)
    context = MPIntervalContext()
    context.prec = 1600
    l_block = _l_rows(context, parent)
    charts = {}
    for family, chart_name in CHARTS:
        if family == "long_reference":
            chart = _long_chart(context, parent, chart_name, l_block)
        else:
            chart = _restoration_chart(context, parent, chart_name, l_block)
        charts[family + ":" + chart_name] = chart
    fixture = _fixture_bound_check()
    reciprocal_fixture = reciprocal_fixture_check()

    pins = {
        INTEGRITY_REPORT: _pin(
            root, INTEGRITY_REPORT, INTEGRITY_REPORT_SHA256
        ),
        CONSTRUCTION_PARENT: _pin(
            root, CONSTRUCTION_PARENT, CONSTRUCTION_PARENT_SHA256
        ),
        DIRECTED_PRIMITIVE: _pin(
            root, DIRECTED_PRIMITIVE, DIRECTED_PRIMITIVE_SHA256
        ),
        RECIPROCAL_HELPER: _pin(
            root, RECIPROCAL_HELPER, RECIPROCAL_HELPER_SHA256
        ),
        C6_SOURCE: _pin(root, C6_SOURCE, parent["input_hashes"][Path(C6_SOURCE).name]),
        LONG_NATIVE: _pin(root, LONG_NATIVE, LONG_NATIVE_SHA256),
        RESTORATION_NATIVE: _pin(
            root, RESTORATION_NATIVE, RESTORATION_NATIVE_SHA256
        ),
        MODULE: _pin(root, MODULE, sha(root / MODULE), self_pin=True),
    }
    result = dict(
        passed=True,
        stage="original_native_control_directed_local_repair",
        schema="original_native_control_directed_repair_v1",
        passed_semantics=(
            "local_directed_abs_input_inverse_L_and_restoration_inverse_Utheta_valid; "
            "inherited_and_downstream_caps_not_recertified"
        ),
        parent_scope=parent_scope,
        source_family=copy.deepcopy(parent["source_family"]),
        source_pins=pins,
        source_guards=guards,
        primitive_api=dict(
            module=DIRECTED_PRIMITIVE,
            function="abs_log_upper",
            operation="exact endpoint tuple magnitude then interval-context ln",
            exact_zero_sentinel=True,
            scalar_abs_endpoint_rounding=False,
            signed_function_value_selected=False,
            source_executed=False,
        ),
        reciprocal_helper_api=dict(
            module=RECIPROCAL_HELPER,
            function="directed_log_reciprocal",
            ordinary_Z_orders=list(ORDERS),
            lower_log_endpoint_tuple_flipped_in_new_helper=True,
            stored_logarithms_exponentiated=False,
            source_executed=False,
        ),
        parameter_logs=dict(
            delta_source_handle="actual_core_signed_direction.parameters.delta",
            delta_is_signed_interval=True,
            z_source_handles=[
                LONG_NATIVE + ":build.z",
                RESTORATION_NATIVE + ":build.z",
            ],
            z_interval_formula=Z_SOURCE_FORMULA,
            source_handle_is_not_numeric_value=True,
        ),
        charts=charts,
        fixture_bound_check=fixture,
        reciprocal_fixture=reciprocal_fixture,
        preservation=dict(
            inherited_historical_inverse_L_caps_retained=True,
            inherited_log_caps_not_treated_as_signed_values=True,
            signed_interval_data_not_deleted_or_midpointed=True,
            old_native_modules_edited=False,
            accepted_builders_or_producers_called=False,
            accepted_checkers_or_self_checks_called=False,
            accepted_recurrences_or_quadratures_called=False,
            generic_helpers_used_on_changed_inputs=[
                RECIPROCAL_HELPER,
                DIRECTED_PRIMITIVE,
            ],
        ),
        checks=dict(
            passed=True,
            parent_integrity_hash_checked=True,
            construction_parent_hash_checked=True,
            source_pins_checked=True,
            source_guards_checked=True,
            safe_abs_log_fixture_passed=True,
            signed_U_M_leaf_checkpoints_built=True,
            directed_L_rows_built=True,
            directed_inverse_L_rows_built=True,
            directed_restoration_post_inverse_Utheta_rows_built=True,
            inverse_Utheta_formula_guarded=True,
            inverse_L_uses_context_interval_negation=True,
            inherited_caps_recertified=False,
            downstream_finite_N_recertified=False,
            downstream_Rm_recertified=False,
            source_executed=False,
        ),
        limits=dict(
            local_blocks_only=True,
            restoration_post_inverse_Utheta_repaired=True,
            long_reference_inverse_Utheta_not_repaired=True,
            gradient_and_value_log_recurrences_not_repaired=True,
            native_full_p2_not_rebuilt=True,
            long_reference_full_caps_not_recertified=True,
            restoration_full_caps_not_recertified=True,
            mixed_loop_and_Rm_not_recertified=True,
            finite_N_not_recertified=True,
            full_goal_complete=False,
        ),
        execution_audit=dict(
            source_executed=False,
            accepted_ancestors_replayed=False,
            accepted_builders_or_owners_replayed=False,
            accepted_checkers_replayed=False,
            accepted_recurrences_replayed=False,
            accepted_quadratures_replayed=False,
            legacy_generic_reciprocal_called_only_on_six_fresh_inputs=True,
            new_signed_interval_algebra_only=True,
            no_legacy_native_builder_imports=True,
            generic_new_helper_imports_are_pinned=True,
            generic_helpers_used_on_changed_inputs=[
                RECIPROCAL_HELPER,
                DIRECTED_PRIMITIVE,
            ],
            no_overclaim=True,
        ),
        source_executed=False,
        accepted_ancestors_replayed=False,
        full_goal_complete=False,
    )
    json.dumps(result, allow_nan=False)
    return result


__all__ = (
    "INTEGRITY_REPORT",
    "INTEGRITY_REPORT_SHA256",
    "CONSTRUCTION_PARENT",
    "CONSTRUCTION_PARENT_SHA256",
    "DIRECTED_PRIMITIVE",
    "DIRECTED_PRIMITIVE_SHA256",
    "LONG_NATIVE",
    "RESTORATION_NATIVE",
    "MODULE",
    "fixture_bound_check",
    "build",
)
