"""Write or inspect the source-bound macro endpoint packet.

Write mode constructs ``OriginalMacroEndpointM15`` once and executes the
independent synthetic fixture once.  Inspection is serialized only: it reads
the report, verifies file pins and packet shape, and never imports the owner,
fixture, producer, checker, or an accepted ancestor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path
import sys
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DIRECTORY = ROOT / "experiments" / "root_st073"
for path in (str(ROOT), str(SRC), str(SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402


NAME = PREFIX + "compliant_current_original_macro_endpoint_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_macro_endpoint_m15.py"
RUNNER_NAME = "scripts/run_original_macro_endpoint_m15.py"
PARENT = PREFIX + "compliant_current_original_switch_m15_numeric.json"
PARENT_SHA = "b4575dd765e43f9f02aac18b5fa11e3825380bd42f6a30aefd5b420bc86ba0e5"
SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PARENT_REDUCED_COUNT = 186
NEW_REDUCED_COUNT = 183
UNCHANGED_COUNT = 183
REMOVED_TERMS = (
    "micro_shared_finite_N_viscous_outer_endpoint_error",
    "macro_finite_N_viscous_error",
    "first_switch_second_finite_N_viscous_external_endpoints",
)
CHARTS = ("second_micro", "frozen_macro", "first_switch")
JOINS = ("micro_macro", "macro_first")
EXPECTED_CHART_HANDLES = {
    "second_micro": {"coordinate": 66988, "radius": 68644, "native_J": 8,
                      "physical_measure": 8, "phase": 97301},
    "frozen_macro": {"coordinate": 70265, "radius": 71448, "native_J": 66987,
                      "physical_measure": 66987, "phase": 101909},
    "first_switch": {"coordinate": 73190, "radius": 85905, "native_J": 9,
                      "physical_measure": 9, "phase": 104537},
}
EXPECTED_JOIN_SCOPES = {
    "micro_macro": (
        ("second_micro", 66988, 48, 68644, 8),
        ("frozen_macro", 70265, 0, 71448, 66987),
    ),
    "macro_first": (
        ("frozen_macro", 70265, 1, 71448, 66987),
        ("first_switch", 73190, 0, 85905, 9),
    ),
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def require(condition: bool, message: Any) -> None:
    if not condition:
        raise AssertionError(message)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve(name: str, base: Path | None = None) -> Path:
    candidate = Path(name)
    if candidate.is_absolute():
        return candidate
    choices = (ROOT / candidate, (base or DIRECTORY) / candidate,
               DIRECTORY / candidate, ROOT / "experiments" / "root_st073" / candidate)
    for path in choices:
        if path.is_file():
            return path
    return choices[0]


def check_hash_group(group: Any, label: str, base: Path | None = None) -> None:
    require(isinstance(group, Mapping) and group, (label, "hash map missing"))
    for name, expected in group.items():
        require(isinstance(name, str) and isinstance(expected, str) and
                re.fullmatch(r"[0-9a-f]{64}", expected), (label, name, "malformed SHA 256"))
        path = _resolve(name, base)
        require(path.is_file() and sha(path) == expected, (label, name, str(path)))


def _rows(value: Any, names: tuple[str, ...]) -> list[Any]:
    if isinstance(value, Mapping):
        for name in names:
            candidate = value.get(name)
            if isinstance(candidate, list):
                return candidate
    return []


def _quantity(row: Mapping[str, Any]) -> str:
    identity = _mapping(row.get("actual_function_identity"))
    raw = identity.get("quantity", row.get("quantity", row.get("field", row.get("name", ""))))
    token = str(raw).strip().replace("_", "").lower()
    aliases = {
        "u": "U", "urow": "U", "urows": "U", "modifiedu": "U", "modifiedurows": "U",
        "uy": "Uy", "uyrow": "Uy", "uyrows": "Uy", "uyphysical": "Uy", "uyphysicalrows": "Uy",
        "v": "V", "vrow": "V", "vrows": "V", "leadingv": "V", "leadingvrows": "V",
        "vy": "Vy", "vyrow": "Vy", "vyrows": "Vy", "vyphysical": "Vy", "vyphysicalrows": "Vy",
        "m": "m", "mrow": "m", "mrows": "m", "modifiedm": "m", "modifiedmrows": "m",
        "leadingm": "m", "leadingmrows": "m",
        "q": "Q", "qrow": "Q", "qrows": "Q", "qy": "Qy", "qyrow": "Qy", "qyrows": "Qy",
        "rqy": "RQy", "rqyrow": "RQy", "rqyrows": "RQy", "rq_y": "RQy",
    }
    return aliases.get(token, {"u_y": "Uy", "v_y": "Vy", "q_y": "Qy", "r_q_y": "RQy"}.get(str(raw), str(raw)))


def _identity_quantity(row: Mapping[str, Any]) -> str:
    return _quantity(row)


def _is_mean(row: Any) -> bool:
    return isinstance(row, Mapping) and _identity_quantity(row) == "m"


def _check_identity(row: Mapping[str, Any], seam: Mapping[str, Any], label: str,
                    *, loop_digest: str | None = None, derivative_digest: str | None = None) -> None:
    identity = _mapping(row.get("actual_function_identity"))
    if not identity and "left_endpoint_sha256" in row:
        identity = row
    require(identity.get("passed") is True, (label, "actual identity missing"))
    for side in ("left", "right"):
        scope_key = side + "_scope"
        if scope_key in identity and scope_key in seam:
            require(identity[scope_key] == seam[scope_key], (label, side, "foreign scope"))
        node_key = side + "_node"
        if node_key in row and node_key in identity:
            require(row[node_key] == identity[node_key], (label, side, "foreign node"))
    if identity.get("exact_endpoint_fingerprints") is True:
        require(identity.get("left_endpoint_sha256") == identity.get("right_endpoint_sha256") and
                re.fullmatch(r"[0-9a-f]{64}", str(identity.get("left_endpoint_sha256"))),
                (label, "direct endpoint fingerprints"))
    else:
        require(identity.get("actual_loop_operator_identity_proved") is True and
                re.fullmatch(r"[0-9a-f]{64}", str(identity.get("actual_loop_identity_sha256"))),
                (label, "loop operator identity"))
        if loop_digest is not None:
            require(identity.get("actual_loop_identity_sha256") == loop_digest,
                    (label, "loop identity digest"))
    if _identity_quantity(row) in ("Uy", "Vy"):
        value = identity.get("issued_derivative_wiring_sha256")
        require(re.fullmatch(r"[0-9a-f]{64}", str(value)), (label, "derivative wiring digest"))
        if derivative_digest is not None:
            require(value == derivative_digest, (label, "derivative wiring digest mismatch"))


def _check_source_rows(join: Mapping[str, Any], label: str) -> dict[str, Any]:
    modified = _rows(join, ("modified_source_rows", "modified_trace_rows", "actual_modified_rows"))
    leading = _rows(join, ("leading_source_rows", "leading_trace_rows", "actual_leading_rows"))
    require(len(modified) == 9 and all(_mapping(row).get("passed") is True for row in modified),
            (label, "modified rows"))
    require(len(leading) == 8 and all(_mapping(row).get("passed") is True for row in leading),
            (label, "leading rows"))
    modified_counts = {name: sum(_quantity(row) == name for row in modified) for name in ("U", "Uy", "m")}
    leading_counts = {name: sum(_quantity(row) == name for row in leading) for name in ("V", "m", "Vy")}
    require(modified_counts == {"U": 3, "Uy": 2, "m": 4}, (label, "modified quantities", modified_counts))
    require(leading_counts == {"V": 3, "m": 3, "Vy": 2}, (label, "leading quantities", leading_counts))
    for row in modified + leading:
        require(re.fullmatch(r"[0-9a-f]{64}", str(row.get("source_row_sha256"))) or
                isinstance(row.get("actual_function_identity"), Mapping),
                (label, "source row digest"))
    return {"modified": modified_counts, "leading": leading_counts}


def _check_pressure_witness(join: Mapping[str, Any], label: str) -> dict[str, Any]:
    witness = _mapping(join.get("pressure_algebra_witness", join.get("pressure_algebra_certificate")))
    require(witness.get("passed") is True, (label, "pressure witness"))
    rows = witness.get("rows")
    require(isinstance(rows, list) and len(rows) == 5, (label, "pressure rows"))
    counts = {name: sum(_quantity(row) == name for row in rows) for name in ("Q", "Qy")}
    require(counts == {"Q": 3, "Qy": 2}, (label, "pressure quantities", counts))
    for row in rows:
        require(_mapping(row).get("passed") is True and
                row.get("reduced_difference_monomials", 0) == 0 and
                row.get("ordinary_derivative_product_rule") is True,
                (label, "pressure identity row"))
        ids = row.get("actual_input_identity_sha256")
        require(isinstance(ids, Mapping) and ids and all(re.fullmatch(r"[0-9a-f]{64}", str(v)) for v in ids.values()),
                (label, "pressure input identities"))
    leading = witness.get("independent_leading_rows")
    require(isinstance(leading, list) and len(leading) == 2 and
            all(_mapping(row).get("passed") is True for row in leading), (label, "leading pressure rows"))
    require(re.fullmatch(r"[0-9a-f]{64}", str(witness.get("formula_AST_sha256"))),
            (label, "pressure formula hash"))
    require(witness.get("source_values_selected") is False and witness.get("numeric_caps_used") is False,
            (label, "pressure value provenance"))
    return {"passed": True, "rows": 5, "leading_rows": 2}


def _join_map(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    value = packet.get("joins", packet.get("seams", packet))
    return value if isinstance(value, Mapping) else {}


def _loop_join(loop: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    for key in ("seams", "joins", "source_seam_proofs"):
        nested = _mapping(loop.get(key))
        if isinstance(nested.get(name), Mapping):
            return nested[name]
    if isinstance(loop.get(name), Mapping):
        return loop[name]
    return loop


def _check_loop_rows(source: Mapping[str, Any], name: str) -> dict[str, Any]:
    loop = _mapping(source.get("actual_loop_function_identity"))
    require(loop.get("passed") is True, (name, "loop identity"))
    row = _loop_join(loop, name)
    roots = row.get("root_rows")
    loops = row.get("loop_rows")
    require(isinstance(roots, list) and len(roots) == 27 and
            all(_mapping(item).get("passed") is True for item in roots), (name, "27 root rows"))
    require(isinstance(loops, list) and len(loops) == 9 and
            all(_mapping(item).get("passed") is True for item in loops), (name, "9 loop rows"))
    require(_mapping(row.get("scoped_fingerprint", {"passed": True})).get("passed", True) is True,
            (name, "scoped fingerprint"))
    require(_mapping(row.get("physical_y_derivative_congruence", {"passed": True})).get("passed", True) is True,
            (name, "physical derivative congruence"))
    # Root/loop source names are case-sensitive; loop q is distinct from pressure Q.
    quantities = {"root": [item.get("quantity") for item in roots],
                  "loop": [item.get("quantity") for item in loops]}
    require(set(quantities["root"]) >= {"E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa"},
            (name, "root quantity coverage"))
    require(set(quantities["loop"]) >= {"q", "B_over_Pstar", "B_over_Pstar_phi"},
            (name, "loop quantity coverage"))
    # Source identities bind to the complete combined loop certificate.  Keep
    # the seam-local ring digest as separate audit evidence rather than using
    # it as the row identity digest.
    return {"passed": True, "root_rows": 27, "loop_rows": 9,
            "loop_digest": _digest(loop), "seam_ring_digest": _digest(row),
            "loop": loop, "row": row}


def _check_actual_function_identities(source_join: Mapping[str, Any], loop_check: Mapping[str, Any],
                                      derivative: Mapping[str, Any], name: str) -> dict[str, Any]:
    """Check the actual source ports, keeping loop/chart receipts separate."""

    loop = _mapping(loop_check.get("loop"))
    loop_digest = str(loop_check.get("loop_digest"))
    charts = derivative.get("charts", derivative.get("chart_guards"))
    charts_digest = derivative.get("charts_sha256")
    if isinstance(charts, Mapping):
        require(re.fullmatch(r"[0-9a-f]{64}", str(charts_digest)) and
                charts_digest == _digest(charts), (name, "native charts digest"))
    identity = _mapping(source_join.get("actual_function_identity"))
    identity_rows = identity.get("rows")
    if not (identity.get("passed") is True and isinstance(identity_rows, list) and len(identity_rows) == 9):
        identity_rows = []
        for key in ("modified_U_rows", "modified_Uy_rows", "modified_m_rows"):
            rows = source_join.get(key)
            require(isinstance(rows, list), (name, key, "actual input rows"))
            identity_rows.extend(rows)
        require(len(identity_rows) == 9, (name, "actual modified identity rows"))
    for row in identity_rows:
        require(isinstance(row, Mapping), (name, "identity row"))
        _check_identity(row, source_join, name + " modified", loop_digest=loop_digest,
                        derivative_digest=str(charts_digest) if charts_digest else None)
    modified_m = source_join.get("modified_m_rows")
    if not isinstance(modified_m, list):
        modified_m = [row for row in identity_rows if _identity_quantity(row) == "m"]
    require(len(modified_m) == 4 and any(_identity_quantity(row) == "m" for row in modified_m),
            (name, "modified m namespace"))
    leading_inputs = _mapping(source_join.get("leading_Qy_inputs"))
    leading_rows = leading_inputs.get("rows")
    require(leading_inputs.get("passed") is True and isinstance(leading_rows, list) and len(leading_rows) == 8,
            (name, "leading source identity rows"))
    for row in leading_rows:
        _check_identity(row, source_join, name + " leading", loop_digest=loop_digest,
                        derivative_digest=str(charts_digest) if charts_digest else None)
    return {"passed": True, "modified_rows": 9, "leading_rows": 8,
            "modified_m_rows": 4, "leading_m_rows": sum(_identity_quantity(row) == "m" for row in leading_rows),
            "loop_digest": loop_digest, "charts_digest": charts_digest}


def _check_native(native: Mapping[str, Any]) -> dict[str, Any]:
    charts = native.get("chart_guards", native.get("charts"))
    require(isinstance(charts, Mapping) and set(CHARTS) <= set(charts), "macro native chart scope changed")
    pressure_ports = _mapping(native.get("pressure_ports"))
    out = {}
    for name in CHARTS:
        chart = _mapping(charts.get(name))
        require(chart.get("passed") is True, (name, "native chart"))
        # NativeDerivative packets may expose the radius guard directly or
        # under the per-chart ``radius_phase_J`` wrapper.  Both are source
        # certificates; do not infer a handle from a numerical field.
        radius_guard = _mapping(chart.get("radius_phase_J"))
        guard = dict(radius_guard)
        guard.update({key: value for key, value in chart.items() if key not in guard})
        expected = EXPECTED_CHART_HANDLES[name]
        for key in ("coordinate", "radius", "native_J", "physical_measure"):
            require(guard.get(key) == expected[key], (name, key, guard.get(key), expected[key]))
        phase = _mapping(guard.get("phase"))
        require(phase.get("actual_phase") == expected["phase"], (name, "phase handle", phase))
        require(guard.get("native_J_positive") is True and
                guard.get("d_native_equals_J_d_y") is True and
                guard.get("d_y_equals_d_native_over_J") is True,
                (name, "positive native J calculus"))
        premise = guard.get("positive_native_J_source_premise")
        if isinstance(premise, Mapping):
            require(premise.get("passed") is True or premise.get("native_J_positive") is True,
                    (name, "native J receipt"))
            require(premise.get("numeric_value_selected", False) is False,
                    (name, "native J numeric value"))
        ports = _mapping(chart.get("physical_pressure_ports"))
        if not ports:
            ports = _mapping(pressure_ports.get(name))
        require(ports.get("passed") is True and ports.get("Pstar_once") is True and
                ports.get("R_Q_y_radius_once") is True and ports.get("ordinary_Z1") is True,
                (name, "physical pressure ports"))
        out[name] = {"passed": True, "coordinate": guard.get("coordinate"),
                     "native_J": guard.get("native_J"), "phase": expected["phase"], "ordinary_Z1": True}
    require(native.get("physical_R_once") is True or _mapping(native.get("premises")).get("physical_R_once") is True,
            "native physical R premise")
    require(native.get("ordinary_Z1") is True or _mapping(native.get("premises")).get("ordinary_Z1") is True,
            "native ordinary Z1 premise")
    return {"passed": True, "charts": out}


def _check_source_pins(report: Mapping[str, Any], packets: tuple[Mapping[str, Any], ...]) -> dict[str, Any]:
    """Require every serialized source pin to be a registered, live hash."""

    input_hashes = _mapping(report.get("input_hashes"))
    implementation_hashes = _mapping(report.get("implementation_hashes"))
    checked = 0
    for packet in packets:
        pins = packet.get("source_pins", packet.get("source_hashes", {}))
        if not isinstance(pins, Mapping):
            continue
        for key, value in pins.items():
            pin = _mapping(value)
            relative = pin.get("path", pin.get("relative", key))
            expected = pin.get("sha256", pin.get("actual_sha256"))
            require(isinstance(relative, str) and re.fullmatch(r"[0-9a-f]{64}", str(expected)),
                    ("source pin", key, "malformed"))
            path = _resolve(relative, DIRECTORY)
            require(path.is_file() and sha(path) == expected, ("source pin", relative, "bytes changed"))
            registered = (input_hashes.get(relative) or input_hashes.get(Path(relative).name) or
                          implementation_hashes.get(relative) or implementation_hashes.get(Path(relative).name))
            require(registered == expected, ("source pin", relative, "not registered in report hashes"))
            checked += 1
    require(checked > 0, "macro source pin map is empty")
    return {"passed": True, "pins": checked}


def _check_source_certificates(report: Mapping[str, Any]) -> dict[str, Any]:
    source = _mapping(report.get("source_derivative_theorem"))
    native = _mapping(report.get("native_source_guards"))
    seams = _mapping(report.get("source_seam_proofs"))
    require(source.get("passed") is True and native.get("passed") is True and seams.get("passed") is True,
            "macro source certificates failed")
    for label, packet in (("source", source), ("native", native), ("seams", seams)):
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if key in packet:
                require(packet.get(key) is False, (label, key, "source-only flag"))
    joins = _join_map(seams)
    require(set(joins) == set(JOINS), ("macro seam names", sorted(joins)))
    source_joins = _join_map(source)
    traces = {}
    loop_checks = {}
    source_join_rows = {}
    for name in JOINS:
        join = _mapping(joins[name])
        require(join.get("passed") is True, (name, "seam failed"))
        for side, expected in zip(("left_scope", "right_scope"), EXPECTED_JOIN_SCOPES[name]):
            scope = _mapping(join.get(side))
            require(scope.get("chart") == expected[0] and
                    scope.get("coordinate_node") == expected[1] and
                    scope.get("native_query_node") == expected[2] and
                    scope.get("radius_node") == expected[3] and
                    scope.get("physical_measure_node") == expected[4],
                    (name, side, "actual endpoint scope", scope))
            require(scope.get("native_query_node") != 40 if name == "micro_macro" and side == "left" else True,
                    (name, side, "incorrect micro query node 40"))
        trace = _rows(join, ("trace_rows", "physical_trace_rows"))
        require(len(trace) == 15 and all(_mapping(row).get("passed") is True for row in trace),
                (name, "15 physical trace rows"))
        quantities = [_quantity(row) for row in trace]
        require({q: quantities.count(q) for q in ("U", "Uy", "m", "Q", "Qy", "RQy")} ==
                {"U": 3, "Uy": 2, "m": 3, "Q": 3, "Qy": 2, "RQy": 2},
                (name, "trace quantity coverage", quantities))
        row_check = _check_source_rows(join, name)
        pressure = _check_pressure_witness(join, name)
        d_rows = join.get("D_equality_rows")
        require(isinstance(d_rows, list) and len(d_rows) == 2 and
                all(_mapping(row).get("passed") is True and row.get("physical_R_once") is True for row in d_rows),
                (name, "D rows"))
        require(join.get("endpoint_wrapper_binding_passed") is True and
                join.get("missing_physical_wrappers_are_local_recipes") is True,
                (name, "endpoint wrappers"))
        radius_phase = _mapping(join.get("radius_phase_identity"))
        if radius_phase:
            require(radius_phase.get("passed") is True and
                    radius_phase.get("global_phase_transfer") is True and
                    radius_phase.get("source_values_selected") is False,
                    (name, "radius/global-phase source identity"))
        # The pressure kernel intentionally carries the derived trace rows but
        # does not duplicate the source theorem receipt.  Bind the premise to
        # the matching source seam instead of accepting a pressure-only flag.
        source_join = _mapping(source_joins.get(name))
        require(source_join.get("passed") is True, (name, "source theorem join"))
        premises = _mapping(join.get("source_proof_premises"))
        if not premises:
            premises = _mapping(source_join.get("source_proof_premises"))
        require(premises.get("passed") is True, (name, "source proof premises"))
        require(premises.get("source_values_selected", False) is False and
                premises.get("numeric_caps_used", False) is False,
                (name, "source proof value provenance"))
        history_flag = premises.get("nonzero_histories_retained")
        if history_flag is None:
            receipt = _mapping(premises.get("receipt"))
            history_flag = _mapping(receipt.get("source_flags")).get("nonzero_incoming_histories")
        require(history_flag is True, (name, "nonzero incoming histories"))
        for key in ("modified_U_rows", "modified_Uy_rows", "modified_m_rows"):
            rows = source_join.get(key)
            require(isinstance(rows, list), (name, key, "missing source rows"))
        ids = source_join.get("actual_function_identity")
        if isinstance(ids, Mapping):
            require(ids.get("passed") is True, (name, "actual function identity"))
        loop_checks[name] = _check_loop_rows(source, name)
        require(source_join.get("global_loop_identity_sha256") == loop_checks[name]["loop_digest"],
                (name, "combined loop identity digest"))
        source_join_rows[name] = source_join
        traces[name] = {"passed": True, "trace_rows": 15, "modified": row_check,
                        "pressure": pressure, "loop": loop_checks[name]}
    derivative = _mapping(source.get("native_derivative_congruence", native.get("native_derivative_congruence")))
    require(derivative.get("passed") is True, "macro native derivative certificate")
    chart_map = derivative.get("charts", derivative.get("chart_guards"))
    if isinstance(chart_map, Mapping):
        require(set(CHARTS) <= set(chart_map), "macro derivative chart scope")
        charts_digest = derivative.get("charts_sha256")
        if charts_digest is not None:
            require(charts_digest == _digest(chart_map), "macro derivative chart digest")
    identities = {}
    for name in JOINS:
        identities[name] = _check_actual_function_identities(
            source_join_rows[name], loop_checks[name], derivative, name)
    _check_native(native)
    pin_check = _check_source_pins(report, (source, native, seams))
    return {"passed": True, "interfaces": list(JOINS), "trace_rows_each": 15,
            "root_rows_each": 27, "loop_rows_each": 9, "traces": traces,
            "actual_function_identities": identities, "source_pins": pin_check,
            "native_source_only": True}


def _check_fixture_shape(fixture: Any) -> dict[str, Any]:
    value = _mapping(fixture)
    require(value.get("passed") is True and value.get("synthetic_local_germs_only") is True,
            "macro endpoint fixture failed")
    require(value.get("finite_N") == 17 and value.get("interface_count") == 2,
            "macro endpoint fixture N/interfaces")
    for key in ("original_source_values_evaluated", "original_Poisson_inverse_evaluated",
                "original_inverse_or_NS_evaluated", "NS_residual"):
        require(value.get(key) is False, ("fixture", key))
    require(value.get("physical_R_once") is True and value.get("ordinary_Z1_nonzero") is True and
            value.get("nonzero_incoming_mean") is True and value.get("genuine_physical_y_derivatives") is True,
            "fixture physical scope")
    rows = value.get("finite_N_error_rows")
    require(isinstance(rows, list) and len(rows) == 2 and
            all(_mapping(row).get("ordinary_Z_orders") == [0, 1] and
                row.get("positive_error_bound_atom_only") is True and row.get("source_value_selected") is False
                for row in rows), "fixture finite-N rows")
    controls = _mapping(value.get("controls"))
    required = {
        "shared_finite_N_once", "phase_reset_rejected", "micro_positive_varying_J",
        "macro0_positive_varying_J", "macro1_positive_varying_J", "first_positive_varying_J",
        "micro_physical_y_derivative_required", "macro0_physical_y_derivative_required",
        "macro1_physical_y_derivative_required", "first_physical_y_derivative_required",
        "nonzero_incoming_mean", "leading_subtraction_retained", "physical_R_required",
        "ordinary_Z1_nonzero", "finite_N_error_atoms_are_not_values", "no_cap_overlap",
        "original_source_values_not_selected", "not_original_inverse_or_NS",
    }
    require(required <= set(controls) and all(controls[key] is True for key in required),
            ("fixture controls", controls))
    chain = _mapping(value.get("signed_chain"))
    require(chain.get("ordinary_Z_orders") == [0, 1] and chain.get("physical_R_once") is True and
            all(token in str(chain.get("function")) for token in ("second_micro", "frozen_macro", "first_switch")),
            "fixture signed macro chain")
    differences = value.get("signed_chain_C0_Z1_differences")
    require(isinstance(differences, list) and len(differences) == 2 and differences == ["0", "0"],
            "fixture exact signed chain")
    require(str(value.get("maximum_supporting_difference")) == "0", "fixture exact supporting difference")
    return {"passed": True, "finite_N": 17, "interface_count": 2,
            "control_count": value.get("control_count"), "maximum_supporting_difference": "0"}


def _check_budget(sample: Mapping[str, Any]) -> dict[str, Any]:
    budget = _mapping(sample.get("external_endpoint_budget_witness"))
    if not budget:
        # The endpoint adapter may expose its two positive atoms under the
        # numerical binding instead of repeating the old switch key.  In that
        # case the packet-level source binding still has to be explicit.
        return {"passed": True, "not_serialized": True, "cap_as_value": False}
    for key in ("positive_union_inputs", "positive_union"):
        if key in budget:
            if key == "positive_union_inputs":
                require(isinstance(budget[key], list) and len(budget[key]) >= 1,
                        (sample.get("Z"), "positive budget inputs"))
            else:
                require(_mapping(budget[key]).get("source_value_selected", False) is False and
                        _mapping(budget[key]).get("cap_as_value", False) is False,
                        (sample.get("Z"), "positive budget value provenance"))
    require(not budget.get("cap_as_value", False) and not budget.get("overlap_subtraction", False),
            (sample.get("Z"), "budget cap/overlap"))
    return {"passed": True, "cap_as_value": False}


def _check_samples(report: Mapping[str, Any], parent: Mapping[str, Any]) -> list[dict[str, Any]]:
    samples = report.get("samples")
    parent_samples = {row.get("Z"): row for row in parent.get("samples", []) if isinstance(row, Mapping)}
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "macro endpoint saved Z set")
    output = []
    for sample in samples:
        require(isinstance(sample, Mapping) and sample.get("Z") in parent_samples, ("sample", sample.get("Z")))
        base = parent_samples[sample["Z"]]
        raw = sample.get("raw_m15_C0_C1_terms")
        previous = base.get("correlated_reduced_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == base.get("raw_m15_C0_C1_terms"),
                (sample.get("Z"), "raw 213 changed"))
        require(isinstance(previous, list) and len(previous) == PARENT_REDUCED_COUNT and
                isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
                (sample.get("Z"), "reduced counts changed"))
        removed = [row.get("name") for row in previous if _mapping(row).get("name") in REMOVED_TERMS]
        require(sorted(removed) == sorted(REMOVED_TERMS), (sample.get("Z"), "removed terms", removed))
        keep = [row for row in previous if row.get("name") not in REMOVED_TERMS]
        require(len(keep) == UNCHANGED_COUNT and reduced == keep,
                (sample.get("Z"), "other 183 terms changed or replacement added"))
        require(sum(row.get("name") == "reserved_n1_own_swirl_once" for row in reduced) == 1,
                (sample.get("Z"), "reserved term count"))
        require(sample.get("independent_target") == base.get("independent_target"),
                (sample.get("Z"), "independent target changed"))
        provenance = sample.get("removed_endpoint_provenance")
        require(isinstance(provenance, list) and len(provenance) == len(REMOVED_TERMS) and
                {row.get("name") for row in provenance} == set(REMOVED_TERMS) and
                all(row.get("active") is False and row.get("provenance_only") is True and
                    row.get("no_cap_value_or_midpoint") is True for row in provenance),
                (sample.get("Z"), "removed endpoint provenance"))
        for row in reduced:
            if isinstance(row, Mapping):
                binding = _mapping(row.get("actual_function_binding"))
                require(binding.get("cap_as_value", False) is False and
                        binding.get("source_value_selected", False) is False,
                        (sample.get("Z"), row.get("name"), "cap/value provenance"))
        chain = _mapping(sample.get("signed_trace_chain", sample.get("signed_macro_chain")))
        require(chain, (sample.get("Z"), "signed macro chain missing"))
        require(chain.get("ordinary_Z_orders", [0, 1]) == [0, 1] and
                chain.get("physical_R_once") is True and
                all(token in str(chain.get("function", "")) for token in ("micro", "macro", "first")),
                (sample.get("Z"), "signed macro chain"))
        require(chain.get("exact_function_identity_zero") is True and
                chain.get("individual_endpoint_values_zero", False) is False,
                (sample.get("Z"), "category identity versus individual endpoint values"))
        contract = _mapping(sample.get("contract"))
        for key in ("raw_213_verbatim", "other_183_verbatim", "independent_target_unchanged",
                    "reserved_retained_once", "ordinary_Z1_and_physical_R_once",
                    "source_functions_not_cap_values", "ancestors_not_reexecuted"):
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        for key in ("six_signed_inner_lower_retained", "three_endpoint_sectors_removed_once",
                    "deleted_endpoint_atoms_provenance_only", "no_active_endpoint_budget",
                    "source_bound_signed_trace_contraction", "previous_186_to_new_183"):
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        require(any(contract.get(key) is True for key in
                    ("no_error_cap_as_value", "no_error_cap_subtraction_or_midpoint")),
                (sample.get("Z"), "sample cap contract"))
        require(any(contract.get(key) is True for key in
                    ("micro_macro_source_join_proved", "micro_macro_source_seam_proved")),
                (sample.get("Z"), "micro/macro source contract"))
        require(any(contract.get(key) is True for key in
                    ("macro_first_source_join_proved", "macro_first_source_seam_proved")),
                (sample.get("Z"), "macro/first source contract"))
        limits = sample.get("limits")
        if isinstance(limits, Mapping):
            require(all(value is False or ("internal" in str(key).lower() and "trace" in str(key).lower())
                        for key, value in limits.items()), (sample.get("Z"), "sample limits"))
        output.append({"Z": sample.get("Z"), "raw_213": True, "reduced_183": True,
                       "unchanged_183": True, "budget": _check_budget(sample), "target_unchanged": True})
    return output


def _load_parent() -> Mapping[str, Any]:
    path = DIRECTORY / PARENT
    require(path.is_file() and sha(path) == PARENT_SHA, "accepted switch parent changed")
    parent = read_json(path)
    require(isinstance(parent, Mapping) and _mapping(parent.get("checks")).get("passed") is True,
            "accepted switch parent is not passed")
    return parent


def _check_packet(report: Mapping[str, Any], parent: Mapping[str, Any],
                  fixture_summary: Mapping[str, Any] | None = None) -> dict[str, Any]:
    require(report.get("source_family") == parent.get("source_family"), "macro source family changed")
    source = _check_source_certificates(report)
    samples = _check_samples(report, parent)
    binding = _mapping(report.get("numerical_source_binding"))
    require(binding.get("passed") is True, "macro numerical source binding")
    require(binding.get("physical_radius_power", 1) in (1, True) and
            binding.get("physical_R_once", True) is True and
            binding.get("ordinary_Z_orders", [0, 1]) == [0, 1] and
            binding.get("leading_subtraction_retained", True) is True and
            binding.get("cap_is_not_signed_function_value", True) is True,
            "macro numerical source scope")
    require(binding.get("deleted_endpoint_atoms_provenance_only") is True and
            binding.get("no_active_endpoint_budget") is True and
            binding.get("no_cap_subtraction_or_midpoint_closure") is True,
            "macro endpoint category provenance")
    for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
        if key in binding:
            require(binding[key] is False, ("numeric binding", key))
    contract = _mapping(report.get("contract"))
    require(contract.get("accepted_ancestor_constructors_not_run") is True and
            contract.get("source_functions_not_cap_values") is True and
            contract.get("raw_213_and_independent_target_retained") is True,
            "macro top contract")
    limits = report.get("limits")
    require(isinstance(limits, Mapping) and all(value is False for value in limits.values()),
            "macro global limits changed")
    if fixture_summary is not None:
        require(fixture_summary.get("passed") is True, "fixture summary failed")
    return {"source": source, "samples": samples, "global_limits_all_false": True,
            "cap_as_value": False}


def run() -> dict[str, Any]:
    started = time.perf_counter()
    from check_original_macro_endpoint_m15 import run_fixture

    fixture = run_fixture()
    fixture_summary = _check_fixture_shape(fixture)
    from openai_ns_reconstruction.original_macro_endpoint_m15 import OriginalMacroEndpointM15

    owner = OriginalMacroEndpointM15(precision=150)
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
        "parent_reduced_186_to_new_183": True,
        "unchanged_other_sectors_183": True,
        "removed_terms_exactly_once": list(REMOVED_TERMS),
        "two_actual_source_interfaces": True,
        "micro_macro_and_macro_first_source_proved": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    data = read_json(path)
    require(isinstance(data, Mapping) and _mapping(data.get("checks")).get("passed") is True,
            "saved macro endpoint report is not passed")
    check_hash_group(data.get("input_hashes"), "input", DIRECTORY)
    check_hash_group(data.get("implementation_hashes"), "implementation", ROOT)
    parent = _load_parent()
    fixture_summary = _check_fixture_shape(_mapping(data.get("checks")).get("independent_fixture"))
    packet_check = _check_packet(data, parent, fixture_summary)
    checks = _mapping(data.get("checks"))
    for key, expected in {
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_186_to_new_183": True,
        "unchanged_other_sectors_183": True,
        "removed_terms_exactly_once": list(REMOVED_TERMS),
        "two_actual_source_interfaces": True,
        "micro_macro_and_macro_first_source_proved": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
    }.items():
        require(checks.get(key) == expected, ("saved check", key))
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
        inspect(path)
