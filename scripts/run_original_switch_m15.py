"""Write or inspect the first-switch to second-switch m15 packet.

Write mode constructs the new switch owner once and runs the independent local
germ fixture once.  Inspection reads only the serialized report and verifies
its hashes, source certificates, signed contraction, and packet invariants;
it never constructs the owner or imports the fixture.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments" / "root_st073"
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
import sys

sys.path[:0] = [str(SRC), str(SCRIPTS), str(ROOT / "experiments" / "root_st073")]

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402


NAME = PREFIX + "compliant_current_original_switch_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_switch_m15.py"
RUNNER_NAME = "scripts/run_original_switch_m15.py"
PARENT = PREFIX + "compliant_current_original_upstream_m15_numeric.json"
PARENT_SHA = "6028d49ac9ec1b94ae967b15d698366061b4a15b9e0ddef9b45be2d3971e355e"

SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PARENT_REDUCED_COUNT = 187
NEW_REDUCED_COUNT = 186
UNCHANGED_COUNT = 185
FIRST_SWITCH_RAW_INDEX = 50
OLD_FIRST_TERM = "first_switch_finite_N_viscous_error"
OLD_EXTERNAL_TERM = "second_switch_long_finite_N_viscous_external_endpoints"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def require(condition: bool, message: Any) -> None:
    if not condition:
        raise AssertionError(message)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve(name: str, base: Path) -> Path:
    candidate = Path(name)
    choices = []
    if candidate.is_absolute():
        choices.append(candidate)
    else:
        choices.extend((base / candidate, ROOT / candidate, DIRECTORY / candidate,
                        ROOT / "experiments" / "root_st073" / candidate))
    for path in choices:
        if path.is_file():
            return path
    return choices[0]


def check_hash_group(group: Any, label: str, base: Path) -> None:
    require(isinstance(group, Mapping) and group, f"{label} hash map missing")
    for name, expected in group.items():
        require(isinstance(name, str) and isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
                (label, name, "invalid SHA 256"))
        path = _resolve(name, base)
        require(path.is_file() and sha(path) == expected, (label, name, str(path)))


def _find_flag(value: Any, names: set[str]) -> bool:
    if not isinstance(value, Mapping):
        return False
    return any(value.get(name) is True for name in names)


def _quantity(row: Mapping[str, Any]) -> str:
    value = row.get("quantity", row.get("field", ""))
    return str(value).replace("_y", "y").replace("Y", "y")


def _identity_quantity(row: Mapping[str, Any]) -> str:
    """Use the canonical nested source identity quantity when available."""
    identity = _mapping(row.get("actual_function_identity"))
    value = identity.get("quantity")
    if value is None:
        value = row.get("quantity", row.get("field", ""))
    token = str(value).strip().replace("_", "").lower()
    if token in {"uy", "uyrow", "uyrows", "uyphysical", "uyphysicalrow", "uyphysicalrows", "modifieduy", "modifieduyrows"}:
        return "Uy"
    if token in {"vy", "vyrow", "vyrows", "vyphysical", "vyphysicalrow", "vyphysicalrows"}:
        return "Vy"
    if token in {"m", "mrow", "mrows", "modifiedm", "modifiedmrow", "modifiedmrows", "leadingm", "leadingmrow", "leadingmrows"}:
        return "m"
    return _quantity({"quantity": value})


def _check_record(record: Any, label: str) -> None:
    row = _mapping(record)
    require(row.get("positive_error_bound_atom_only") is True, (label, "not a positive atom"))
    require(row.get("source_value_selected") is False, (label, "source value selected"))
    tuples = row.get("source_enclosure_exact_mpf_tuples")
    require(isinstance(tuples, list) and len(tuples) >= 1, (label, "source enclosure missing"))
    require(not row.get("cap_as_value", False) and not row.get("signed_function_value", False),
            (label, "cap used as value"))


def _check_positive_union(inputs: Any, union: Any, label: str) -> dict[str, Any]:
    require(isinstance(inputs, list) and len(inputs) == 2, (label, "two active atoms required"))
    for index, item in enumerate(inputs):
        _check_record(item, f"{label} input {index}")
    result = _mapping(union)
    require(result.get("positive_error_bound_atom_only") is True or result.get("source_value_selected") is False,
            (label, "union is not an enclosure"))
    require(result.get("source_value_selected", False) is False, (label, "union source value selected"))
    require(not result.get("cap_as_value", False) and not result.get("overlap_subtraction", False),
            (label, "union overlap or cap subtraction"))
    return {"passed": True, "components": len(inputs), "cap_as_value": False}


def _rows(value: Mapping[str, Any], names: tuple[str, ...]) -> list[Any]:
    for name in names:
        candidate = value.get(name)
        if isinstance(candidate, list):
            return candidate
    return []


def _check_identity(row: Mapping[str, Any], seam: Mapping[str, Any], label: str,
                    *, loop_digest: str | None = None,
                    derivative_digest: str | None = None) -> None:
    identity = row.get("actual_function_identity")
    # The source theorem stores modified identities in a compact sibling
    # array, while leading rows wrap theirs under this field.
    if not isinstance(identity, Mapping) and "left_endpoint_sha256" in row:
        identity = row
    require(isinstance(identity, Mapping) and identity.get("passed") is True, (label, "identity missing"))
    for side in ("left", "right"):
        if f"{side}_node" in identity and f"{side}_node" in row:
            require(identity[f"{side}_node"] == row[f"{side}_node"], (label, side, "foreign node"))
        scope_key = f"{side}_scope"
        if scope_key in identity and scope_key in seam:
            require(identity[scope_key] == seam[scope_key], (label, side, "foreign scope"))
    if identity.get("exact_endpoint_fingerprints") is True:
        require(identity.get("left_endpoint_sha256") == identity.get("right_endpoint_sha256") and
                re.fullmatch(r"[0-9a-f]{64}", str(identity.get("left_endpoint_sha256"))),
                (label, "unequal direct fingerprints"))
    else:
        require(identity.get("actual_loop_operator_identity_proved") is True and
                re.fullmatch(r"[0-9a-f]{64}", str(identity.get("actual_loop_identity_sha256"))),
                (label, "unproved loop identity"))
        if loop_digest is not None:
            require(identity.get("actual_loop_identity_sha256") == loop_digest,
                    (label, "loop certificate digest mismatch"))
    if _identity_quantity(row) == "Uy" and identity.get("exact_endpoint_fingerprints") is not True:
        value = identity.get("issued_derivative_wiring_sha256")
        require(re.fullmatch(r"[0-9a-f]{64}", str(value)), (label, "issued derivative certificate missing"))
        if derivative_digest is not None:
            require(value == derivative_digest, (label, "derivative chart digest mismatch"))


def _mean_quantity(value: Any) -> bool:
    """Recognize only serialized mean-row labels, without aliasing fields."""
    token = str(value).strip().lower().replace("_", "")
    return token in {
        "m", "mrow", "mrows", "modifiedm", "modifiedmrow", "modifiedmrows",
        "leadingm", "leadingmrow", "leadingmrows",
    }


def _is_mean_row(row: Any) -> bool:
    if not isinstance(row, Mapping):
        return False
    if _mean_quantity(row.get("quantity", row.get("field", ""))):
        return True
    identity = _mapping(row.get("actual_function_identity"))
    return _mean_quantity(identity.get("quantity"))


def _mean_nodes(rows: Any, *, explicit_group: bool = False) -> set[Any]:
    if not isinstance(rows, list):
        return set()
    output = set()
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        if not explicit_group and not _is_mean_row(row):
            continue
        for key in ("left_node", "right_node"):
            node = row.get(key)
            if node is not None:
                output.add(node)
    return output


def _direct_constant_nodes(source_row: Mapping[str, Any], groups: tuple[str, ...]) -> set[Any]:
    """Return nodes supported as exact constants by both endpoint fingerprints."""
    direct = _mapping(source_row.get("direct_endpoint_fingerprints"))
    group_map = _mapping(direct.get("groups"))
    output = set()
    for name in groups:
        rows = _mapping(group_map.get(name)).get("rows")
        if not isinstance(rows, list):
            continue
        for row in rows:
            if not isinstance(row, Mapping):
                continue
            if (row.get("left_rational") is None or row.get("right_rational") is None or
                    row.get("left_endpoint_sha256") != row.get("right_endpoint_sha256")):
                continue
            for key in ("left_node", "right_node"):
                node = row.get(key)
                if node is not None:
                    output.add(node)
    return output


def _check_source_rows(seam: Mapping[str, Any], expected: Mapping[str, int], label: str) -> dict[str, Any]:
    rows = _rows(seam, ("modified_source_rows", "modified_trace_rows", "actual_modified_rows"))
    require(len(rows) == 9, (label, "modified source row count", len(rows)))
    counts = {key: 0 for key in expected}
    for row in rows:
        require(isinstance(row, Mapping) and row.get("passed") is True, (label, "modified row failed"))
        quantity = _quantity(row)
        if quantity in counts:
            counts[quantity] += 1
        require(isinstance(row.get("source_row_sha256"), str) or isinstance(row.get("actual_function_identity"), Mapping),
                (label, "source row digest missing"))
    require(counts == dict(expected), (label, "modified quantities", counts))
    leading = _rows(seam, ("leading_source_rows", "leading_trace_rows", "actual_leading_rows"))
    require(len(leading) == 8, (label, "leading source row count", len(leading)))
    lead_counts = {"V": 0, "m": 0, "Vy": 0}
    for row in leading:
        require(isinstance(row, Mapping) and row.get("passed") is True, (label, "leading row failed"))
        quantity = str(row.get("quantity", row.get("field", "")))
        quantity = {"V_y": "Vy", "VY": "Vy", "Uy": "Vy", "Vy": "Vy"}.get(quantity, quantity)
        if quantity in lead_counts:
            lead_counts[quantity] += 1
    require(lead_counts == {"V": 3, "m": 3, "Vy": 2}, (label, "leading quantities", lead_counts))
    return {"modified": counts, "leading": lead_counts}


def _check_pressure_witness(seam: Mapping[str, Any], label: str) -> dict[str, Any]:
    witness = _mapping(seam.get("pressure_algebra_witness", seam.get("pressure_algebra_certificate")))
    require(witness.get("passed") is True, (label, "pressure algebra witness"))
    rows = witness.get("rows")
    require(isinstance(rows, list) and len(rows) == 5, (label, "five pressure rows"))
    quantities = [_quantity(row) for row in rows if isinstance(row, Mapping)]
    require(quantities.count("Q") == 3 and quantities.count("Qy") == 2,
            (label, "pressure quantity coverage", quantities))
    for row in rows:
        require(isinstance(row, Mapping) and row.get("passed") is True and
                row.get("reduced_difference_monomials", 0) == 0 and
                row.get("ordinary_derivative_product_rule") is True,
                (label, "pressure row"))
        identities = row.get("actual_input_identity_sha256", {})
        require(isinstance(identities, Mapping) and identities and
                all(re.fullmatch(r"[0-9a-f]{64}", str(v)) for v in identities.values()),
                (label, "pressure input identities"))
    leading = witness.get("independent_leading_rows")
    require(isinstance(leading, list) and len(leading) == 2 and all(row.get("passed") is True for row in leading),
            (label, "independent leading rows"))
    require(re.fullmatch(r"[0-9a-f]{64}", str(witness.get("formula_AST_sha256"))),
            (label, "formula AST pin"))
    require(witness.get("source_values_selected") is False and witness.get("numeric_caps_used") is False,
            (label, "pressure cap/value contract"))
    return {"passed": True, "rows": 5, "leading_rows": 2}


def _seam_map(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    seams = _mapping(packet.get("source_seam_proofs"))
    value = seams.get("joins", seams.get("seams", seams))
    return value if isinstance(value, Mapping) else {}


def _check_function_certificate(source: Mapping[str, Any], native: Mapping[str, Any], seam: Mapping[str, Any]) -> dict[str, Any]:
    loop = _mapping(source.get("actual_loop_function_identity"))
    require(loop.get("passed") is True, "actual loop function identity failed")
    # The current source adapter emits one first_second loop certificate at
    # the packet root.  Keep support for the older nested form only when it
    # is explicitly present; never manufacture a seam wrapper in the check.
    if "seams" in loop or "source_seam_proofs" in loop:
        loop_seams = _mapping(loop.get("seams", loop.get("source_seam_proofs")))
        loop_row = _mapping(loop_seams.get("first_second"))
    else:
        loop_row = loop
    require(loop_row.get("passed") is True, "first_second loop identity missing")
    roots = loop_row.get("root_rows")
    loops = loop_row.get("loop_rows")
    require(isinstance(roots, list) and len(roots) == 27 and all(row.get("passed") is True for row in roots),
            "first_second root rows must be 27")
    require(isinstance(loops, list) and len(loops) == 9 and all(row.get("passed") is True for row in loops),
            "first_second loop rows must be 9")
    for identity_key in ("scoped_fingerprint", "physical_y_derivative_congruence"):
        identity_value = loop_row.get(identity_key, True)
        require(identity_value is not False and
                (not isinstance(identity_value, Mapping) or identity_value.get("passed") is True),
                ("first_second", identity_key, "scoped source identity incomplete"))

    # Validate the certificate hash in the same pre-field form used by the
    # source adapter.  Its certificate_sha256 is a receipt for the object
    # before that field is attached, while operator rows use the digest of the
    # complete object that pressure_congruence consumes.
    loop_certificate = loop.get("certificate_sha256")
    if loop_certificate is not None:
        loop_payload = dict(loop)
        loop_payload.pop("certificate_sha256", None)
        require(re.fullmatch(r"[0-9a-f]{64}", str(loop_certificate)) and
                loop_certificate == _digest(loop_payload),
                "loop certificate digest mismatch")
    loop_digest = _digest(loop)

    # The source theorem carries the actual field identities separately from
    # the compact pressure trace.  Keep this check distinct so a pressure
    # trace row cannot stand in for a source-function identity.  The current
    # source packet stores identities on each modified row; an older packet
    # may additionally supply an aggregate rows array.
    source_seams = _mapping(source.get("seams", source.get("source_seam_proofs")))
    source_row = _mapping(source_seams.get("first_second"))
    require(source_row.get("passed") is True, "first_second source theorem missing")
    identity = _mapping(source_row.get("actual_function_identity"))
    identity_rows = identity.get("rows")
    if not (identity.get("passed") is True and isinstance(identity_rows, list) and len(identity_rows) == 9):
        identity_rows = []
        for key in ("modified_U_rows", "modified_Uy_rows", "modified_m_rows"):
            rows = source_row.get(key)
            require(isinstance(rows, list), ("first_second", key, "source rows missing"))
            identity_rows.extend(rows)
        require(len(identity_rows) == 9, "first_second modified source identities missing")
    for row in identity_rows:
        _check_identity(row, source_row, "first_second actual modified",
                        loop_digest=loop_digest)
    # Use the explicit modified_m_rows group when present.  The serialized
    # row quantity is allowed to be ``m_rows`` in the source adapter and must
    # not be confused with the nested identity quantity ``m``.
    modified_mean_rows = source_row.get("modified_m_rows")
    if not isinstance(modified_mean_rows, list):
        modified_mean_rows = [row for row in identity_rows
                              if _is_mean_row(row)]
    require(len(modified_mean_rows) == 4, "modified m source namespace order coverage changed")
    modified_m_nodes = _mean_nodes(modified_mean_rows, explicit_group=True)
    leading_inputs = _mapping(source_row.get("leading_Qy_inputs"))
    leading_rows = leading_inputs.get("rows")
    require(leading_inputs.get("passed") is True and isinstance(leading_rows, list) and len(leading_rows) == 8,
            "first_second leading source identities missing")
    for row in leading_rows:
        _check_identity(row, source_row, "first_second actual leading",
                        loop_digest=loop_digest)
    leading_mean_rows = [row for row in leading_rows
                         if _is_mean_row(row)]
    require(len(leading_mean_rows) == 3, "leading m source namespace order coverage changed")
    leading_m_nodes = _mean_nodes(leading_mean_rows, explicit_group=True)
    require(modified_m_nodes, "modified m source namespace is empty")
    require(leading_m_nodes, "leading m source namespace is empty")
    overlap = modified_m_nodes & leading_m_nodes
    shared_exact_constants: set[Any] = set()
    nonconstant_overlap: set[Any] = set()
    if overlap:
        modified_constants = _direct_constant_nodes(source_row, ("modified_m",))
        leading_constants = _direct_constant_nodes(source_row, ("leading_m",))
        shared_exact_constants = overlap & modified_constants & leading_constants
        nonconstant_overlap = overlap - shared_exact_constants
        require(not nonconstant_overlap,
                ("leading m namespace overlaps modified m namespace at nonconstant nodes",
                 sorted(nonconstant_overlap)))
    premises = _mapping(source_row.get("source_proof_premises"))
    require(premises.get("nonzero_histories_retained") is True,
            "nonzero leading and incoming histories are not retained")

    derivative = _mapping(source.get("native_derivative_congruence"))
    if not derivative:
        derivative = _mapping(native.get("native_derivative_congruence"))
    require(derivative.get("passed") is True, "issued native derivative certificate missing")
    parent_pin = _mapping(source.get("parent_report_pin"))
    require(source.get("parent_report_sha256") == PARENT_SHA and
            parent_pin.get("passed") is True and
            parent_pin.get("selected_native_guards_match_canonical") is True,
            "canonical parent native evidence is not hash-bound")
    charts = _mapping(derivative.get("charts", derivative.get("chart_guards")))
    require(set(charts) == {"first_switch", "second_switch"},
            "first/second native derivative chart scope changed")
    charts_digest = derivative.get("charts_sha256")
    require(re.fullmatch(r"[0-9a-f]{64}", str(charts_digest)) and
            charts_digest == _digest(charts),
            "native derivative chart digest mismatch")
    derivative_certificate = derivative.get("certificate_sha256")
    if derivative_certificate is not None:
        derivative_payload = dict(derivative)
        derivative_payload.pop("certificate_sha256", None)
        require(re.fullmatch(r"[0-9a-f]{64}", str(derivative_certificate)) and
                derivative_certificate == _digest(derivative_payload),
                "native derivative certificate digest mismatch")
    chart_rows = {}
    for chart_name in ("first_switch", "second_switch"):
        chart = _mapping(charts.get(chart_name))
        require(chart.get("passed") is True, (chart_name, "native derivative chart"))
        radius = _mapping(chart.get("radius_phase_J"))
        slow = _mapping(chart.get("slow_B_native_derivative"))
        pressure = _mapping(chart.get("physical_pressure_ports"))
        radial = chart.get("radial_Bphi_source_rows")
        finite = chart.get("finite_N_source_rows")
        positive = _mapping(radius.get("positive_native_J_source_premise"))
        require(positive.get("passed") is True and
                positive.get("parent_report_sha256") == PARENT_SHA and
                positive.get("actual_source_node_bound") is True and
                positive.get("numeric_value_selected") is False and
                positive.get("coordinate") == radius.get("coordinate") and
                positive.get("radius") == radius.get("radius") and
                positive.get("native_J") == radius.get("native_J"),
                (chart_name, "canonical native Jacobian positivity receipt"))
        require(radius.get("passed") is True and
                radius.get("d_native_equals_J_d_y") is True and
                radius.get("d_y_equals_d_native_over_J") is True and
                radius.get("native_J_positive") is True,
                (chart_name, "native Jacobian source guard"))
        slow_rows = slow.get("rows")
        require(slow.get("passed") is True and isinstance(slow_rows, list) and len(slow_rows) == 2,
                (chart_name, "slow B native derivative rows"))
        for row in slow_rows:
            require(row.get("inverse_J_binding") is True and
                    row.get("lazy_native_derivative_binding") is True and
                    row.get("loop_B_source_binding") is True,
                    (chart_name, "native J/B source row"))
        require(isinstance(radial, list) and len(radial) == 2 and all(row.get("passed") is True for row in radial),
                (chart_name, "B phi source rows"))
        for row in radial:
            require(row.get("Bphi_in_total_radial_sum") is True and
                    row.get("slow_B_y_product") is True and
                    row.get("native_inverse_J_once") is True and
                    row.get("common_inverse_N_once") is True and
                    row.get("total_fast_phase_terms_included") is True and
                    row.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
                    (chart_name, "B phi/J source wiring"))
        require(isinstance(finite, list) and len(finite) == 3 and
                all(row.get("root_plus_B_over_N") is True for row in finite),
                (chart_name, "finite N source rows"))
        require(pressure.get("passed") is True and pressure.get("Pstar_once") is True and
                pressure.get("R_Q_y_radius_once") is True and pressure.get("ordinary_Z1") is True and
                pressure.get("Q_formula_source") is True and pressure.get("Q_plus_Q_y_rows") is True,
                (chart_name, "physical pressure source ports"))
        orders = chart.get("ordinary_Z_orders", [0, 1])
        require(orders == [0, 1], (chart_name, "ordinary Z 0/1"))
        chart_rows[chart_name] = {"passed": True, "ordinary_Z_orders": [0, 1]}
    uy_rows = source_row.get("modified_Uy_rows")
    if not isinstance(uy_rows, list):
        uy_rows = [row for row in identity_rows if _identity_quantity(row) == "Uy"]
    require(len(uy_rows) == 2, "modified Uy source namespace order coverage changed")
    derivative_hashes = []
    for row in uy_rows:
        identity_value = _mapping(row.get("actual_function_identity"))
        if _identity_quantity(row) == "Uy" and identity_value.get("exact_endpoint_fingerprints") is not True:
            derivative_hashes.append(identity_value.get("issued_derivative_wiring_sha256"))
    require(derivative_hashes and all(value == charts_digest for value in derivative_hashes),
            "native derivative chart wiring SHA 256 missing or mismatched")
    return {"passed": True, "root_rows": 27, "loop_rows": 9,
            "charts": chart_rows, "issued_certificate_sha256": charts_digest,
            "loop_certificate_sha256": loop_digest,
            "native_charts_sha256": charts_digest,
            "modified_m_node_count": len(modified_m_nodes),
            "leading_m_node_count": len(leading_m_nodes),
            "shared_exact_constant_m_nodes": sorted(shared_exact_constants) if overlap else [],
            "nonconstant_m_namespace_overlap": sorted(nonconstant_overlap) if overlap else []}


def _check_source_certificates(report: Mapping[str, Any]) -> dict[str, Any]:
    source = _mapping(report.get("source_derivative_theorem"))
    native = _mapping(report.get("native_source_guards"))
    seams = _seam_map(report)
    require(source.get("passed") is True and native.get("passed") is True, "source certificates failed")
    for label, cert in (("source", source), ("native", native)):
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if key in cert:
                require(cert.get(key) is False, (label, key))
    require(set(seams) == {"first_second"}, "exactly one first_second source seam required")
    seam = _mapping(seams["first_second"])
    require(seam.get("passed") is True, "first_second source seam failed")
    trace = seam.get("trace_rows")
    require(isinstance(trace, list) and len(trace) == 15 and all(row.get("passed") is True for row in trace),
            "first_second physical trace must have 15 rows")
    quantities = [_quantity(row) for row in trace if isinstance(row, Mapping)]
    require(quantities.count("U") == 3 and quantities.count("Uy") == 2 and
            quantities.count("m") == 3 and quantities.count("Q") == 3 and
            quantities.count("Qy") == 2 and quantities.count("RQy") == 2,
            ("first_second", "trace quantities", quantities))
    rows = _check_source_rows(seam, {"U": 3, "Uy": 2, "m": 4}, "first_second")
    pressure = _check_pressure_witness(seam, "first_second")
    d_rows = seam.get("D_equality_rows")
    require(isinstance(d_rows, list) and len(d_rows) == 2 and
            all(row.get("passed") is True and row.get("physical_R_once") is True for row in d_rows),
            "first_second D rows")
    require(seam.get("endpoint_wrapper_binding_passed") is True and
            seam.get("missing_physical_wrappers_are_local_recipes") is True,
            "first_second endpoint wrappers")
    native_flags = set(native) | set(_mapping(native.get("premises")))
    require(any(native.get(key) is True for key in ("physical_R_once", "physical_radius_once", "R_factor_once")),
            "physical R source premise missing")
    require(any(native.get(key) is True for key in ("ordinary_Z1", "ordinary_Z_orders_0_1")) or
            _mapping(native.get("scope")).get("ordinary_Z_orders") == [0, 1],
            "ordinary Z 1 source premise missing")
    require(any(native.get(key) is True for key in
                ("first_native_J_positive", "first_switch_native_J_positive", "positive_native_J")) and
            any(native.get(key) is True for key in ("second_native_J_positive", "second_switch_native_J_positive", "positive_native_J")),
            "first/second positive native Jacobian premise missing")
    function = _check_function_certificate(source, native, seam)
    return {"passed": True, "interfaces": ["first_second"], "native_source_only": True,
            "trace_rows_each": 15, "modified_rows": rows, "pressure_algebra_witness": pressure,
            "actual_function_identity": function}


def _check_fixture_shape(fixture: Any) -> dict[str, Any]:
    value = _mapping(fixture)
    require(value.get("passed") is True and value.get("synthetic_local_germs_only") is True,
            "switch fixture failed")
    require(value.get("finite_N") == 17 and value.get("interface_count") == 1,
            "switch fixture N/interface count changed")
    for key in ("original_source_values_evaluated", "original_Poisson_inverse_evaluated", "NS_residual"):
        require(value.get(key) is False, ("fixture", key))
    require(value.get("physical_R_once") is True and value.get("ordinary_Z1_nonzero") is True and
            value.get("leading_subtraction_nonzero") is True,
            "fixture physical source scope changed")
    rows = value.get("finite_N_error_rows")
    require(isinstance(rows, list) and len(rows) == 1 and rows[0].get("interface") == "first_second" and
            rows[0].get("ordinary_Z_orders") == [0, 1] and
            rows[0].get("positive_error_bound_atom_only") is True and
            rows[0].get("source_value_selected") is False,
            "fixture finite N atom row changed")
    controls = _mapping(value.get("controls"))
    required = {
        "shared_finite_N_once", "phase_reset_rejected", "first_positive_J",
        "second_positive_varying_J", "first_physical_y_derivative_required",
        "second_physical_y_derivative_required", "nonzero_incoming_mean",
        "leading_subtraction_retained", "physical_R_required", "ordinary_Z1_nonzero",
        "finite_N_error_atoms_are_not_values", "no_cap_overlap",
    }
    require(required <= set(controls) and all(controls[key] is True for key in required),
            "switch fixture controls changed")
    chain = _mapping(value.get("signed_chain"))
    require(chain.get("ordinary_Z_orders") == [0, 1] and chain.get("physical_R_once") is True and
            "first_switch native0" in str(chain.get("function")), "fixture signed chain changed")
    differences = value.get("signed_chain_C0_Z1_differences")
    require(isinstance(differences, list) and len(differences) == 2 and
            all(float(str(item)) < 1e-70 for item in differences), "fixture signed chain is not exact")
    require(float(str(value.get("maximum_supporting_difference"))) < 1e-70,
            "fixture supporting difference changed")
    return {"passed": True, "finite_N": 17, "interface_count": 1,
            "control_count": value.get("control_count"),
            "maximum_supporting_difference": value.get("maximum_supporting_difference")}


def _check_budget(sample: Mapping[str, Any]) -> dict[str, Any]:
    budget = _mapping(sample.get("external_endpoint_budget_witness"))
    require(budget, "switch external endpoint budget missing")
    own = None
    for key, value in budget.items():
        if isinstance(value, Mapping) and "own" in key.lower() and "first" in key.lower():
            own = value
            break
    require(isinstance(own, Mapping), "own first-switch native 0 budget missing")
    atom_rows = None
    for key in ("radius_weighted_C0_C1_log_bounds", "radius_weighted_C0_C1_bounds", "C0_C1_log_bounds"):
        if isinstance(own.get(key), list):
            atom_rows = own[key]
            break
    require(isinstance(atom_rows, list) and len(atom_rows) == 2, "own first-switch C0/Z1 rows changed")
    for row in atom_rows:
        _check_record(row, "own first-switch")
    inputs = budget.get("positive_union_inputs")
    union = _check_positive_union(inputs, budget.get("positive_union"), "switch positive union")
    require(budget.get("no_cap_subtraction_value_midpoint_or_overlap") is True or
            budget.get("no_cap_as_value_or_overlap") is True,
            "switch union subtraction contract changed")
    require(budget.get("raw50_and_previous_external_provenance_only") is True or
            budget.get("raw50_full_aggregate_used") is False or
            budget.get("raw50_and_old_external_retained_as_provenance_only") is True,
            "switch raw provenance contract missing")
    source_binding = _mapping(budget.get("source_binding"))
    require(source_binding.get("physical_R_once", True) is True and
            source_binding.get("ordinary_Z_orders", [0, 1]) == [0, 1] and
            source_binding.get("cap_is_not_signed_function_value", True) is True,
            "switch budget source scope changed")
    return {"passed": True, "components": union["components"], "cap_as_value": False}


def _check_samples(report: Mapping[str, Any], parent: Mapping[str, Any]) -> list[dict[str, Any]]:
    samples = report.get("samples")
    parent_samples = {row.get("Z"): row for row in parent.get("samples", []) if isinstance(row, Mapping)}
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "switch saved Z set changed")
    output = []
    replacement_names = set()
    removed_set = {OLD_FIRST_TERM, OLD_EXTERNAL_TERM}
    for sample in samples:
        require(isinstance(sample, Mapping) and sample.get("Z") in parent_samples, ("sample", sample.get("Z")))
        base = parent_samples[sample["Z"]]
        raw = sample.get("raw_m15_C0_C1_terms")
        previous = base.get("correlated_reduced_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == base.get("raw_m15_C0_C1_terms"),
                (sample.get("Z"), "raw 213 changed"))
        require(raw[FIRST_SWITCH_RAW_INDEX].get("name") == OLD_FIRST_TERM,
                (sample.get("Z"), "raw first-switch term changed"))
        require(isinstance(previous, list) and len(previous) == PARENT_REDUCED_COUNT and
                isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
                (sample.get("Z"), "term counts changed"))
        removed = [row.get("name") for row in previous if isinstance(row, Mapping) and row.get("name") in removed_set]
        require(sorted(removed) == sorted(removed_set), (sample.get("Z"), "removed terms", removed))
        keep = [row for row in previous if row.get("name") not in removed_set]
        require(len(keep) == UNCHANGED_COUNT, (sample.get("Z"), "unchanged count"))
        keep_names = {row.get("name") for row in keep}
        replacement = [row for row in reduced if row.get("name") not in keep_names]
        require(len(replacement) == 1, (sample.get("Z"), "single replacement"))
        replacement = replacement[0]
        replacement_name = replacement.get("name")
        require(isinstance(replacement_name, str) and all(token in replacement_name.lower() for token in ("first", "second", "external", "viscous")) and
                replacement_name not in removed_set, (sample.get("Z"), "replacement identity", replacement_name))
        replacement_names.add(replacement_name)
        require([row for row in reduced if row.get("name") in keep_names] == keep,
                (sample.get("Z"), "other 185 sectors changed"))
        _binding = _mapping(replacement.get("actual_function_binding"))
        require(_binding and (_find_flag(_binding, {"source_functions_not_cap_values", "error_enclosure_only", "source_bound", "signed_source_function"}) or
                              isinstance(replacement.get("log_scale_terms"), list)), (sample.get("Z"), "replacement source binding"))
        require(sum(row.get("name") == "reserved_n1_own_swirl_once" for row in reduced) == 1,
                (sample.get("Z"), "reserved term count"))
        require(sample.get("independent_target") == base.get("independent_target"), (sample.get("Z"), "target changed"))
        chain = _mapping(sample.get("signed_trace_chain"))
        aliases = _mapping(chain.get("aliases"))
        require(aliases.get("first_switch_native1") == "second_switch_native0" and
                chain.get("ordinary_Z_orders") == [0, 1] and chain.get("physical_R_once") is True,
                (sample.get("Z"), "signed first-second chain"))
        contract = _mapping(sample.get("contract"))
        for key in ("raw_213_verbatim", "other_185_verbatim", "independent_target_unchanged",
                    "reserved_retained_once", "ordinary_Z1_and_physical_R_once",
                    "leading_subtraction_retained", "source_functions_not_cap_values",
                    "no_error_cap_as_value", "all_other_volume_width_numerical_errors_retained",
                    "ancestors_not_reexecuted"):
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        require(any(contract.get(key) is True for key in ("first_second_source_join_proved", "first_switch_second_switch_source_join_proved")),
                (sample.get("Z"), "first-second source join contract"))
        limits = sample.get("limits")
        if isinstance(limits, Mapping):
            require(all(value is False or ("internal" in str(key).lower() and "trace" in str(key).lower())
                        for key, value in limits.items()), (sample.get("Z"), "sample limits"))
        output.append({"Z": sample.get("Z"), "raw_213": True, "reduced_186": True,
                       "unchanged_185": True, "budget": _check_budget(sample),
                       "replacement": replacement_name, "target_unchanged": True})
    require(len(replacement_names) == 1, "replacement name changed across Z")
    return output


def _load_parent() -> Mapping[str, Any]:
    path = DIRECTORY / PARENT
    require(path.is_file() and sha(path) == PARENT_SHA, "accepted upstream parent changed")
    parent = read_json(path)
    require(isinstance(parent, Mapping) and _mapping(parent.get("checks")).get("passed") is True,
            "accepted upstream parent is not passed")
    return parent


def _check_packet(report: Mapping[str, Any], parent: Mapping[str, Any], fixture_summary: Mapping[str, Any] | None = None) -> dict[str, Any]:
    require(report.get("source_family") == parent.get("source_family"), "switch source family changed")
    controls = report.get("pressure_input_identity_guard_controls")
    require(isinstance(controls, Mapping) and controls.get("passed") is True and controls.get("control_count") == 3 and
            controls.get("controls") == {"missing_identity": True, "foreign_node_pair": True, "unequal_fingerprints": True} and
            controls.get("source_values_evaluated") is False, "pressure input identity controls changed")
    source_check = _check_source_certificates(report)
    samples = _check_samples(report, parent)
    binding = _mapping(report.get("numerical_source_binding"))
    require(binding.get("passed") is True and binding.get("physical_radius_power", 1) in (1, True) and
            binding.get("physical_R_once", True) in (1, True) and
            binding.get("ordinary_Z_orders", [0, 1]) == [0, 1] and
            binding.get("leading_subtraction_retained", True) is True and
            binding.get("cap_is_not_signed_function_value", True) is True,
            "switch numerical source binding changed")
    require(report.get("contract") == {
        "accepted_ancestor_constructors_not_run": True,
        "source_functions_not_cap_values": True,
        "raw_213_and_independent_target_retained": True,
    }, "top contract changed")
    limits = report.get("limits")
    require(isinstance(limits, Mapping) and all(value is False for value in limits.values()),
            "global completion limits changed")
    if fixture_summary is not None:
        require(fixture_summary.get("passed") is True, "fixture summary failed")
    return {"source": source_check, "samples": samples, "global_limits_all_false": True, "cap_as_value": False}


def run() -> dict[str, Any]:
    started = time.perf_counter()
    from check_original_switch_m15 import run_fixture

    fixture = run_fixture()
    fixture_summary = _check_fixture_shape(fixture)
    from openai_ns_reconstruction.original_switch_m15 import OriginalSwitchM15

    owner = OriginalSwitchM15(precision=150)
    packet = owner.build()
    code = dict(getattr(owner, "code_hashes", {}))
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(getattr(owner, "hashes", {})), implementation_hashes=code)
    parent = _load_parent()
    packet_check = _check_packet(report, parent, fixture_summary)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_summary,
        "packet_invariants": packet_check,
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_187_to_new_186": True,
        "unchanged_other_sectors_185": True,
        "first_second_source_join_proved": True,
        "pressure_algebra_witness_checked": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    data = read_json(path)
    require(isinstance(data, Mapping) and _mapping(data.get("checks")).get("passed") is True,
            "saved switch report is not passed")
    check_hash_group(data.get("input_hashes"), "input", DIRECTORY)
    check_hash_group(data.get("implementation_hashes"), "implementation", ROOT)
    parent = _load_parent()
    fixture_summary = _check_fixture_shape(_mapping(data.get("checks")).get("independent_fixture"))
    packet_check = _check_packet(data, parent, fixture_summary)
    checks = _mapping(data.get("checks"))
    for key, value in {
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_187_to_new_186": True,
        "unchanged_other_sectors_185": True,
        "first_second_source_join_proved": True,
        "pressure_algebra_witness_checked": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
    }.items():
        require(checks.get(key) == value, ("saved check", key))
    return {"checks": checks, "fixture_invariants": fixture_summary,
            "packet_invariants": packet_check, "execution_seconds": data.get("execution_seconds")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-report", action="store_true")
    modes.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        data = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
    else:
        require(path.is_file(), ("missing report", str(path)))
        data = inspect(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")}, indent=2), flush=True)
