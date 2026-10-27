"""Unentered full-W signed-sector symbolic kernel draft.

The module is inert at import time.  It does not import the project, execute
source providers, read files, or construct symbolic objects.  ``derive`` first
validates the literal source packet and emits a durable algebra-entry marker;
only then does it lazily import the generic symbolic engine and build fresh
formal expressions.

The direct full-source numerator graph and the eight-sector graph are built
independently.  This is source-conditional algebra and metadata only: it does
not claim source values, a physical gap, a P2 gap, or any radius, phase,
temporal, or global closure identity.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping
from copy import deepcopy


INPUT_SCHEMA = "compact_O3_transition_W_sector_correlation_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_W_sector_correlation"
API_SCHEMA = "compact_O3_transition_W_sector_correlation_api_v1"
ROOT_ENTRY_SCHEMA = "compact_O3_transition_W_sector_correlation_algebra_entry_v1"
SOURCE_FACTS_SHA256 = "eb0af8bce8ab9af4bb72f8bc8a3eb5e03f5a2a6db9688ab30f37104084003dcf"
SYMPY_VERSION = "1.14.0"

ORDERS = (0, 1, 2, 3)
SECTOR_NAMES = ("R", "S", "delta", "V", "m", "e", "P0", "p")
QUIET_POWER_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
)

# These are source AST anchors, not executable source code.  They are checked
# as machine provenance pins; no display-text guard is treated as evidence.
KNOWN_AST_DIGESTS = frozenset(
    {
        "e63ae3283447c56ae542f03f2ab1480bb7e96d1624cc1b83a3cb53b88354f11c",
        "fbfe367ab432bd721d16129adf80595e0f512fe7b8a8bff62cb27378a29bb005",
        "ab33333c0bff2bd25ce8b8a0aa0472f656d43d1628b350103c620311d8205aaa",
        "3ebe8c6fa524fb0c19e0982313fa02c1957b6f919cd9b4ea793988c51560d1fe",
        "ae7ea5b8440f9f654dc556acfa864a7df15cacb594642470eae6efb4b7e83da3",
        "58101254f7a8609a45ed3ee832c47e7bd185338c326795f334e49b3635422648",
        "984da1fa05f96f9d1561df51e67b26ef77f0cf223806d96c8ededadb9d08b87b",
        "f68ccb8f07e6244113f92404aa30cf904176df03bda8f6508600182c4efe0dfc",
    }
)

NORMALIZED_FORMULAS = {
    "A": "-V+(1-delta)*(m-Z*m_Z)/2",
    "B": "V*((1-delta)*Z*m+(1-Z**2)*m_Z)+2*delta*Z*e-(1-Z**2)*e_Z+2*(1+delta)*Z*(P0+p)-(1-Z**2)*(P0_Z+p_Z)",
    "L": "1-delta*Z**2",
    "W": "R*(A+S*B)/L",
}

EXPECTED_FORMULA_FACTS = {
    "native": {
        "A_numerator": {
            "original_ast_dump_sha256": "e63ae3283447c56ae542f03f2ab1480bb7e96d1624cc1b83a3cb53b88354f11c",
            "normalized_text": "A=-V+(1-delta)*(m-Z*m_Z)/2",
            "normalized_text_sha256": "b7dab91730dfb2b76348142a9704f88e414ac939b51dafd33bf3a716c92f3635",
            "quantity": "A",
            "side": "native",
        },
        "B_numerator": {
            "original_ast_dump_sha256": "fbfe367ab432bd721d16129adf80595e0f512fe7b8a8bff62cb27378a29bb005",
            "normalized_text": "B=V*((1-delta)*Z*m+(1-Z**2)*m_Z)+2*delta*Z*e-(1-Z**2)*e_Z+2*(1+delta)*Z*(P0+p)-(1-Z**2)*(P0_Z+p_Z)",
            "normalized_text_sha256": "0849d9ab6400b2a063a535f949eae2311bed4dea7a442b4d8583dd90d93bacc5",
            "quantity": "B",
            "side": "native",
        },
        "L_denominator": {
            "original_ast_dump_sha256": "ab33333c0bff2bd25ce8b8a0aa0472f656d43d1628b350103c620311d8205aaa",
            "normalized_text": "L=1-delta*Z**2",
            "normalized_text_sha256": "74028a72ffb53debaad11af3909ce4173c0db5a0cc4d69279efdcf5b6c4fd1e4",
            "quantity": "L",
            "side": "native",
        },
        "W": {
            "original_ast_dump_sha256": "3ebe8c6fa524fb0c19e0982313fa02c1957b6f919cd9b4ea793988c51560d1fe",
            "normalized_text": "W=R*(A+S*B)/L",
            "normalized_text_sha256": "04b5147a41eeec4dd92c81e254ef407f22f5c9fb94d3b940700770b334902d5b",
            "quantity": "W",
            "side": "native",
        },
    },
    "rebuilt_C3": {
        "A_numerator": {
            "original_ast_dump_sha256": "ae7ea5b8440f9f654dc556acfa864a7df15cacb594642470eae6efb4b7e83da3",
            "normalized_text": "A=-V+(1-delta)*(m-Z*m_Z)/2",
            "normalized_text_sha256": "b7dab91730dfb2b76348142a9704f88e414ac939b51dafd33bf3a716c92f3635",
            "quantity": "A",
            "side": "rebuilt_C3",
        },
        "B_numerator": {
            "original_ast_dump_sha256": "58101254f7a8609a45ed3ee832c47e7bd185338c326795f334e49b3635422648",
            "normalized_text": "B=V*((1-delta)*Z*m+(1-Z**2)*m_Z)+2*delta*Z*e-(1-Z**2)*e_Z+2*(1+delta)*Z*(P0+p)-(1-Z**2)*(P0_Z+p_Z)",
            "normalized_text_sha256": "0849d9ab6400b2a063a535f949eae2311bed4dea7a442b4d8583dd90d93bacc5",
            "quantity": "B",
            "side": "rebuilt_C3",
        },
        "L_denominator": {
            "original_ast_dump_sha256": "984da1fa05f96f9d1561df51e67b26ef77f0cf223806d96c8ededadb9d08b87b",
            "normalized_text": "L=1-delta*Z**2",
            "normalized_text_sha256": "74028a72ffb53debaad11af3909ce4173c0db5a0cc4d69279efdcf5b6c4fd1e4",
            "quantity": "L",
            "side": "rebuilt_C3",
        },
        "W": {
            "original_ast_dump_sha256": "f68ccb8f07e6244113f92404aa30cf904176df03bda8f6508600182c4efe0dfc",
            "normalized_text": "W=R*(A+S*B)/L",
            "normalized_text_sha256": "04b5147a41eeec4dd92c81e254ef407f22f5c9fb94d3b940700770b334902d5b",
            "quantity": "W",
            "side": "rebuilt_C3",
        },
    },
}

EXPECTED_SOURCE_FUNCTIONS = {
    "native": {
        "record": "native_generic_functions",
        "function": "generic_functions",
        "function_ast_dump_sha256": "cfc50a5e623c1d71b1c2551b6d8ed72c83565ff3c74329c58fa0e92a4c7be52f",
    },
    "rebuilt_C3": {
        "record": "rebuilt_C3_recover_inputs",
        "function": "recover_inputs",
        "function_ast_dump_sha256": "b4033a6990b784a673874439d90a6f772737f584b20b27cee823fb303dcf9cac",
    },
}

EXPECTED_PARENT_SHAS = frozenset(
    {
        "2956baf8c82250ffd789714ceb1ee80170bd69e57f3c84d1d6b3b8889a1d1466",
        "958ba5154ea7a532f33d2e48cafeb12209bf94e376b9a0cc5dafae83bfb676a7",
    }
)


class SourceContractError(ValueError):
    """Raised internally for a source packet that is not machine-anchored."""


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(
        ch in "0123456789abcdef" for ch in value.lower()
    )


def _canonical_sha(value: object) -> str:
    """Return a deterministic digest for a loaded literal packet."""

    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _walk_strings(value: object, path: str = ""):
    """Yield machine string leaves and their paths without interpreting text."""

    if isinstance(value, Mapping):
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            yield from _walk_strings(child, child_path)
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            yield from _walk_strings(child, f"{path}[{index}]")
    elif isinstance(value, str):
        yield path, value


def _find_values(value: object, wanted: set[str]) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {item: [] for item in wanted}
    for path, text in _walk_strings(value):
        if text in wanted:
            found[text].append(path)
    return found


def _find_key_values(value: object, wanted_keys: set[str]) -> dict[str, list[tuple[str, str]]]:
    """Collect exact key/value pairs for typed provenance fields."""

    found: dict[str, list[tuple[str, str]]] = {key: [] for key in wanted_keys}

    def visit(node: object, path: str = "") -> None:
        if isinstance(node, Mapping):
            for key, child in node.items():
                child_path = f"{path}.{key}" if path else str(key)
                if key in wanted_keys and isinstance(child, str):
                    found[key].append((child_path, child))
                visit(child, child_path)
        elif isinstance(node, (list, tuple)):
            for index, child in enumerate(node):
                visit(child, f"{path}[{index}]")

    visit(value)
    return found


def _formula_paths(facts: Mapping[str, object]) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {name: [] for name in NORMALIZED_FORMULAS}

    def visit(node: object, path: str = "") -> None:
        if isinstance(node, Mapping):
            for key, child in node.items():
                child_path = f"{path}.{key}" if path else str(key)
                if key in {"normalized_text", "formula", "normalized_formula", "text"}:
                    if isinstance(child, str):
                        for name, expected in NORMALIZED_FORMULAS.items():
                            if child == expected or child == f"{name}={expected}":
                                found[name].append(child_path)
                visit(child, child_path)
        elif isinstance(node, (list, tuple)):
            for index, child in enumerate(node):
                visit(child, f"{path}[{index}]")

    visit(facts)
    return found


def _parent_sha_values(source: Mapping[str, object]) -> list[str]:
    values: list[str] = []
    for key, value in source.items():
        lowered = str(key).lower()
        if "parent" not in lowered or "sha" not in lowered:
            continue
        if isinstance(value, str) and _is_sha256(value):
            values.append(value)
        elif isinstance(value, (list, tuple)):
            values.extend(item for item in value if _is_sha256(item))
        elif isinstance(value, Mapping):
            values.extend(text for _path, text in _walk_strings(value) if _is_sha256(text))
    return sorted(set(values))


def _source_refs(source: Mapping[str, object], facts: Mapping[str, object]) -> dict[str, object]:
    """Extract compact, machine-oriented references for every emitted row."""

    ast_values = _find_values(facts, set(KNOWN_AST_DIGESTS))
    formula_paths = _formula_paths(facts)
    return {
        "input_schema": source.get("schema"),
        "sector_source_facts_raw_sha256": source.get("sector_source_facts_raw_sha256"),
        "parent_result_raw_sha256s": _parent_sha_values(source),
        "ast_digest_paths": {digest: paths for digest, paths in ast_values.items() if paths},
        "normalized_formula_paths": formula_paths,
        "exact_formula_fact_refs": {
            side: sorted(records) for side, records in EXPECTED_FORMULA_FACTS.items()
        },
        "exact_source_function_refs": {
            side: expected["record"] for side, expected in EXPECTED_SOURCE_FUNCTIONS.items()
        },
        "source_fact_packet_canonical_sha256": _canonical_sha(facts),
    }


def _validate_exact_formula_facts(facts: Mapping[str, object]) -> None:
    """Bind each side's four v2 records, rather than recursive presence scans."""

    actual = facts.get("v2_normalized_formula_facts")
    if not isinstance(actual, Mapping):
        raise SourceContractError("v2_normalized_formula_facts mapping is missing")
    for side, expected_records in EXPECTED_FORMULA_FACTS.items():
        side_records = actual.get(side)
        if not isinstance(side_records, Mapping):
            raise SourceContractError(f"formula facts for side {side} are missing")
        for name, expected in expected_records.items():
            record = side_records.get(name)
            if not isinstance(record, Mapping):
                raise SourceContractError(f"formula fact {side}.{name} is missing")
            for field in (
                "original_ast_dump_sha256",
                "normalized_text",
                "normalized_text_sha256",
                "quantity",
                "side",
            ):
                if record.get(field) != expected[field]:
                    raise SourceContractError(f"formula fact {side}.{name}.{field} changed")


def _validate_exact_source_functions(facts: Mapping[str, object]) -> None:
    """Bind exact producer metadata for native and rebuilt C3 source sides."""

    actual = facts.get("source_functions")
    if not isinstance(actual, Mapping):
        raise SourceContractError("source_functions mapping is missing")
    for side, expected in EXPECTED_SOURCE_FUNCTIONS.items():
        record = actual.get(expected["record"])
        if not isinstance(record, Mapping):
            raise SourceContractError(f"source function record for side {side} is missing")
        for field in ("function", "function_ast_dump_sha256"):
            if record.get(field) != expected[field]:
                raise SourceContractError(f"source function {side}.{field} changed")


def _validate_source(source: object) -> tuple[dict[str, object], dict[str, object]]:
    if not isinstance(source, Mapping):
        raise SourceContractError("source must be a mapping")
    if source.get("schema") != INPUT_SCHEMA:
        raise SourceContractError("source schema is not the W-sector input schema")
    if source.get("component") != "compact_O3_transition_W_sector_correlation":
        raise SourceContractError("source component is not the W-sector component")
    facts = source.get("sector_source_facts")
    if not isinstance(facts, Mapping):
        raise SourceContractError("sector_source_facts must be a literal mapping")
    packet_sha = source.get("sector_source_facts_raw_sha256")
    if packet_sha != SOURCE_FACTS_SHA256:
        raise SourceContractError("sector_source_facts_raw_sha256 is not the authenticated v2 trace SHA")
    parent_result = source.get("parent_result")
    if not isinstance(parent_result, Mapping):
        raise SourceContractError("parent_result saved data is missing")

    # SymPy metadata is optional at this leaf.  When the frozen packet carries
    # a version field, it must identify the reviewed generic engine; no missing
    # or alternate metadata key is invented here.
    sympy_metadata = source.get("sympy_raw_pins")
    if isinstance(sympy_metadata, Mapping):
        provided_version = sympy_metadata.get("version")
        if provided_version is not None and provided_version != SYMPY_VERSION:
            raise SourceContractError("provided SymPy version metadata is not 1.14.0")

    # Retained parent anchors are required as raw SHA fields.  Parent text and
    # provider output are never interpreted.
    parent_values = set(_parent_sha_values(source))
    parent_values.update(text for _path, text in _walk_strings(parent_result) if _is_sha256(text))
    missing_parent = EXPECTED_PARENT_SHAS - parent_values
    if missing_parent:
        raise SourceContractError("retained parent raw SHA anchors are incomplete")

    # Exact side-specific facts are the binding guard.  Recursive scans remain
    # only as compact provenance references in emitted rows.
    _validate_exact_formula_facts(facts)
    _validate_exact_source_functions(facts)
    found_ast = _find_values(facts, set(KNOWN_AST_DIGESTS))
    formula_paths = _formula_paths(facts)

    refs = _source_refs(source, facts)
    refs["validated_ast_digest_count"] = len(found_ast)
    refs["validated_formula_count"] = len(formula_paths)
    refs["parent_sha_anchor_count"] = len(parent_values & EXPECTED_PARENT_SHAS)
    return dict(source), refs


def _expr_record(expr: object, sp: object, *, label: str) -> dict[str, object]:
    """Serialize a complete symbolic expression without evaluating source text."""

    return {"label": label, "text": str(expr), "srepr": sp.srepr(expr)}


def _jet_row(
    *,
    order: int,
    name: str,
    expr: object,
    sp: object,
    source_refs: Mapping[str, object],
    counts: Mapping[str, int],
) -> dict[str, object]:
    return {
        "order": order,
        "sector": name,
        "expression": _expr_record(expr, sp, label=f"D_Z^{order}[D_{name}]"),
        "source_references": deepcopy(dict(source_refs)),
        "actual_api_counts": dict(counts),
        "actual_counts": dict(counts),
    }


def _set_counts(result: dict[str, object], counts: Mapping[str, int]) -> None:
    """Keep the workflow-facing and internal count maps byte-for-byte equal."""

    result["actual_api_counts"] = dict(counts)
    result["actual_counts"] = dict(counts)
    for name, value in counts.items():
        result[name] = value


def _set_flags(result: dict[str, object], values: Mapping[str, bool]) -> None:
    """Keep the public flag map and the internal flag mirror synchronized."""

    flags = dict(values)
    result["flags"] = flags
    result["new_flags"] = dict(flags)


def _row_ids(result: Mapping[str, object]) -> dict[str, list[str]]:
    """Return compact durable-prefix identifiers, never duplicated row bodies."""

    ids: dict[str, list[str]] = {}
    for key in ("numerator_jet_rows", "direct_source_jet_rows", "proof_rows"):
        rows = result.get(key)
        if isinstance(rows, list):
            values: list[str] = []
            for row in rows:
                if not isinstance(row, Mapping):
                    continue
                if key == "numerator_jet_rows":
                    values.append(f"{row.get('order')}:{row.get('sector')}")
                elif key == "direct_source_jet_rows":
                    values.append(f"{row.get('order')}:Dsrc")
                else:
                    values.append(f"{row.get('order')}:identity")
            ids[key] = values
    return ids


def _journal_prefix(result: Mapping[str, object], counts: Mapping[str, int]) -> dict[str, object]:
    return {
        "row_ids": _row_ids(result),
        "actual_api_counts": dict(counts),
        "actual_counts": dict(counts),
    }


def _base_result(source_refs: Mapping[str, object], counts: Mapping[str, int]) -> dict[str, object]:
    result = {
        "schema": RESULT_SCHEMA,
        "component": "compact_O3_transition_W_sector_correlation",
        "status": "failed",
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "conditional_source_functions": True,
        "conditional_on_source_port_semantics": True,
        "physical_layer_completed": False,
        "derive_entry_calls": 1,
        "derive_calls": 1,
        "source_references": deepcopy(dict(source_refs)),
        "numerator_jet_rows": [],
        "direct_source_jet_rows": [],
        "proof_rows": [],
        "whole_W_gap_derivative_views": [],
        "flags": {
            "O3_transition_full_W_signed_sector_functions_defined": False,
            "O3_transition_full_W_cross_numerator_correlation_ordinary_Z3_proved": False,
        },
        "claims": {
            "actual_numeric_jets": False,
            "source_W_gap_zero": False,
            "P2_gap_zero": False,
            "phase_equality": False,
            "radius_equality": False,
            "Pstar_equality": False,
            "delta_equality": False,
            "cap_reduction": False,
            "error_reduction": False,
            "physical_closure": False,
            "global_closure": False,
            "temporal_n_recursion": False,
        },
        "preserved_metadata": {
            "old_errors_removed": False,
            "old_error_terms_removed": False,
            "old_errors_removed_global": False,
            "global": False,
            "quiet_power4_retained": True,
        },
        "runtime_model_metadata": None,
        "whole_W_gap_views_unevaluated": True,
        "whole_W_gap_definition": "D_Z^j[(sum_of_8_cross_numerator_sectors)/(L_n*L_c)]",
        "kernel_workflow_coordinate": {
            "ordinary_Z_variable": "Z",
            "held_fixed_scalars": ["Rn", "Rc", "Sn", "Sc", "dn", "dc"],
        },
    }
    for name in QUIET_POWER_FLAGS:
        result[name] = False
    result.update(
        {
            "old_errors_removed": False,
            "old_error_terms_removed": False,
            "old_errors_removed_global": False,
            "global": False,
            "global_five_moment_closure": False,
            "global_closure_proved": False,
            "numeric_values_computed": False,
            "numeric_field_values_present": False,
            "time_order_recursion_claimed": False,
            "full_goal_complete": False,
        }
    )
    _set_counts(result, counts)
    _set_flags(result, result["flags"])
    return result


def _emit(diagnostic: Callable[[dict[str, object]], object] | None, payload: Mapping[str, object]) -> None:
    if diagnostic is not None:
        diagnostic(deepcopy(dict(payload)))


def _make_formal_ports(sp: object) -> dict[str, object]:
    """Construct independent side-local ports and one common formal Z."""

    Z = sp.Symbol("Z")
    d = 1 - Z**2
    ports: dict[str, object] = {"Z": Z, "d": d}
    for name in ("Rn", "Rc", "Sn", "Sc", "dn", "dc"):
        ports[name] = sp.Symbol(name)
    for name in ("Vn", "Vc", "mn", "mc", "en", "ec", "P0n", "P0c", "pn", "pc"):
        ports[name] = sp.Function(name)(Z)
    return ports


def _derivative(expr: object, Z: object, sp: object) -> object:
    """Make an unevaluated first-derivative node without a diff API call."""

    return sp.Derivative(expr, Z, evaluate=False)


def _build_direct_source(ports: Mapping[str, object], sp: object) -> tuple[object, object, object]:
    """Build independent inline A/B/L source facts for n and c."""

    Z, d = ports["Z"], ports["d"]
    Vn, Vc = ports["Vn"], ports["Vc"]
    mn, mc = ports["mn"], ports["mc"]
    en, ec = ports["en"], ports["ec"]
    P0n, P0c = ports["P0n"], ports["P0c"]
    pn, pc = ports["pn"], ports["pc"]
    dn, dc = ports["dn"], ports["dc"]
    Rn, Rc, Sn, Sc = ports["Rn"], ports["Rc"], ports["Sn"], ports["Sc"]

    # Deliberately inline both sides.  This branch does not call the sector
    # graph's h/T/F helpers and does not share those helper aliases.
    Ln_inline = 1 - dn * Z**2
    Lc_inline = 1 - dc * Z**2
    An_inline = -Vn + (1 - dn) * (mn - Z * _derivative(mn, Z, sp)) / 2
    Ac_inline = -Vc + (1 - dc) * (mc - Z * _derivative(mc, Z, sp)) / 2
    Bn_inline = (
        Vn * ((1 - dn) * Z * mn + d * _derivative(mn, Z, sp))
        + 2 * dn * Z * en
        - d * _derivative(en, Z, sp)
        + 2 * (1 + dn) * Z * (P0n + pn)
        - d * (_derivative(P0n, Z, sp) + _derivative(pn, Z, sp))
    )
    Bc_inline = (
        Vc * ((1 - dc) * Z * mc + d * _derivative(mc, Z, sp))
        + 2 * dc * Z * ec
        - d * _derivative(ec, Z, sp)
        + 2 * (1 + dc) * Z * (P0c + pc)
        - d * (_derivative(P0c, Z, sp) + _derivative(pc, Z, sp))
    )
    Fn_inline = An_inline + Sn * Bn_inline
    Fc_inline = Ac_inline + Sc * Bc_inline
    Dsrc = Rn * Lc_inline * Fn_inline - Rc * Ln_inline * Fc_inline
    return Dsrc, Ln_inline, Lc_inline


def _build_sector_graph(ports: Mapping[str, object], sp: object) -> tuple[dict[str, object], object, object]:
    """Construct the frozen eight-sector graph with separate helper aliases."""

    Z, d = ports["Z"], ports["d"]
    Rn, Rc, Sn, Sc = ports["Rn"], ports["Rc"], ports["Sn"], ports["Sc"]
    dn, dc = ports["dn"], ports["dc"]
    Vn, Vc = ports["Vn"], ports["Vc"]
    mn, mc = ports["mn"], ports["mc"]
    en, ec = ports["en"], ports["ec"]
    P0n, P0c = ports["P0n"], ports["P0c"]
    pn, pc = ports["pn"], ports["pc"]

    Ln = 1 - dn * Z**2
    Lc = 1 - dc * Z**2
    An = -Vn + (1 - dn) * (mn - Z * _derivative(mn, Z, sp)) / 2
    Ac = -Vc + (1 - dc) * (mc - Z * _derivative(mc, Z, sp)) / 2
    Bn = (
        Vn * ((1 - dn) * Z * mn + d * _derivative(mn, Z, sp))
        + 2 * dn * Z * en
        - d * _derivative(en, Z, sp)
        + 2 * (1 + dn) * Z * (P0n + pn)
        - d * (_derivative(P0n, Z, sp) + _derivative(pn, Z, sp))
    )
    Bc = (
        Vc * ((1 - dc) * Z * mc + d * _derivative(mc, Z, sp))
        + 2 * dc * Z * ec
        - d * _derivative(ec, Z, sp)
        + 2 * (1 + dc) * Z * (P0c + pc)
        - d * (_derivative(P0c, Z, sp) + _derivative(pc, Z, sp))
    )
    Fn = An + Sn * Bn
    Fc = Ac + Sc * Bc

    h_c = mc - Z * _derivative(mc, Z, sp)
    T_c = (1 - dc) * Z * mc + d * _derivative(mc, Z, sp)
    q_c = P0c + pc
    delta_R = Rn - Rc
    delta_S = Sn - Sc
    delta_delta = dn - dc
    delta_V = Vn - Vc
    delta_m = mn - mc
    delta_m_prime = _derivative(mn, Z, sp) - _derivative(mc, Z, sp)
    delta_e = en - ec
    delta_e_prime = _derivative(en, Z, sp) - _derivative(ec, Z, sp)
    delta_P0 = P0n - P0c
    delta_P0_prime = _derivative(P0n, Z, sp) - _derivative(P0c, Z, sp)
    delta_p = pn - pc
    delta_p_prime = _derivative(pn, Z, sp) - _derivative(pc, Z, sp)

    sectors = {
        "R": delta_R * Lc * Fc,
        "S": Rn * Lc * delta_S * Bc,
        "delta": (
            Rc * Z**2 * delta_delta * Fc
            + Rn * Lc * delta_delta * (
                -h_c / 2 + Sn * (-Vn * Z * mc + 2 * Z * ec + 2 * Z * q_c)
            )
        ),
        "V": Rn * Lc * delta_V * (-1 + Sn * T_c),
        "m": Rn
        * Lc
        * (
            (1 - dn) * (delta_m - Z * delta_m_prime) / 2
            + Sn * Vn * ((1 - dn) * Z * delta_m + d * delta_m_prime)
        ),
        "e": Rn * Lc * Sn * (2 * dn * Z * delta_e - d * delta_e_prime),
        "P0": Rn * Lc * Sn * (2 * (1 + dn) * Z * delta_P0 - d * delta_P0_prime),
        "p": Rn * Lc * Sn * (2 * (1 + dn) * Z * delta_p - d * delta_p_prime),
    }
    return sectors, Ln, Lc, Rn * Fn / Ln, Rc * Fc / Lc


def _build_derivative_view(*, order: int, Z: object) -> dict[str, object]:
    """Return an unevaluated data tree, without another symbolic diff call."""

    quotient = "(D_R + D_S + D_delta + D_V + D_m + D_e + D_P0 + D_p)/(Ln*Lc)"
    if order == 0:
        return {
            "order": 0,
            "operator": "Identity",
            "evaluated": False,
            "argument": {
                "operator": "quotient",
                "numerator": "sum_of_eight_cross_numerator_sectors",
                "denominator": "Ln*Lc",
            },
            "formal_expression": {
                "head": "Identity",
                "variables": [],
                "argument": quotient,
            },
            "domain_premises": ["Ln != 0", "Lc != 0"],
            "source_conditional": True,
            "actual_symbolic_diff_calls": 0,
        }
    return {
        "order": order,
        "operator": "Derivative_Z",
        "evaluated": False,
        "variable": str(Z),
        "argument": {
            "operator": "quotient",
            "numerator": "sum_of_eight_cross_numerator_sectors",
            "denominator": "Ln*Lc",
        },
        "formal_expression": {
            "head": "Derivative",
            "variables": [str(Z)] * order,
            "argument": quotient,
        },
        "domain_premises": ["Ln != 0", "Lc != 0"],
        "source_conditional": True,
        "actual_symbolic_diff_calls": 0,
    }


def derive(source: Mapping[str, object], diagnostic: Callable[[dict[str, object]], object] | None = None) -> dict[str, object]:
    """Derive the fresh source-conditional signed W-sector algebra.

    Validation is AST/literal-only.  The algebra-entry event is emitted before
    the lazy generic symbolic import and before symbolic construction.  Every
    counted API call is incremented immediately before that call; a failure
    returns the durable prefix and actual counts.
    """

    counts: dict[str, int] = {
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "source_bound_sides": 0,
        "sector_families": 0,
        "numerator_jet_rows": 0,
        "cross_numerator_identity_orders": 0,
        "whole_W_gap_derivative_views": 0,
        "symbolic_diff_api_calls": 0,
        "symbolic_expand_api_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
    }
    source_refs: dict[str, object] = {}
    try:
        source_value, source_refs = _validate_source(source)
    except Exception as exc:
        result = _base_result(source_refs, counts)
        result.update(
            {
                "status": "source_validation_failed",
                "algebra_entry": False,
                "failure": {
                    "stage": "source_validation",
                    "type": type(exc).__name__,
                    "message": str(exc),
                    "saved_prefix": True,
                    "mathwork_counters_zero": True,
                    "actual_api_counts": dict(counts),
                    "actual_counts": dict(counts),
                },
            }
        )
        return result

    facts = source_value["sector_source_facts"]
    counts["source_bound_sides"] = 2
    result = _base_result(source_refs, counts)
    result["source_formula_premises"] = {
        "conditional": True,
        "normalized_formulas": dict(NORMALIZED_FORMULAS),
        "anchored_to_source_asts": True,
        "provider_value_equality": False,
        "source_fact_packet_is_literal": True,
        "source_fact_packet_type": type(facts).__name__,
    }
    result["domain_premises"] = [
        "common formal coordinate Z",
        "Ln != 0 wherever Wn is interpreted as a quotient",
        "Lc != 0 wherever Wc is interpreted as a quotient",
        "R, S, delta and all source fields remain side-local; no parameter equality",
    ]

    # Durable root entry: this event precedes the only symbolic-engine import
    # and all expression construction.
    counts["new_algebra_family_calls"] += 1
    result["algebra_entry"] = {
        "schema": ROOT_ENTRY_SCHEMA,
        "entered": True,
        "source_validated": True,
        "symbolic_engine_imported": False,
        "symbolic_construction_started": False,
        "callback_completed": False,
        "source_provider_calls": counts["source_provider_calls"],
        "old_math_calls": counts["old_math_calls"],
    }
    _set_counts(result, counts)

    # The first diagnostic write is itself guarded.  If its fsync-backed
    # manager callback fails, return a durable failed prefix without importing
    # Sympy or constructing any symbolic object.
    try:
        _emit(
            diagnostic,
            {
                "schema": ROOT_ENTRY_SCHEMA,
                "event": "algebra_entry",
                "payload": deepcopy(result),
                "prefix": _journal_prefix(result, counts),
            },
        )
        result["algebra_entry"]["callback_completed"] = True
    except Exception as exc:
        result["status"] = "failed_initial_algebra_entry_callback"
        result["passed"] = False
        _set_flags(
            result,
            {
                "O3_transition_full_W_signed_sector_functions_defined": False,
                "O3_transition_full_W_cross_numerator_correlation_ordinary_Z3_proved": False,
            },
        )
        _set_counts(result, counts)
        result["failure"] = {
            "stage": "algebra_entry_callback",
            "type": type(exc).__name__,
            "message": str(exc),
            "saved_prefix": True,
            "mathwork_counters_zero": (
                counts["symbolic_diff_api_calls"] == 0
                and counts["symbolic_expand_api_calls"] == 0
                and counts["new_mathematical_proofs"] == 0
            ),
            "actual_api_counts": dict(counts),
            "actual_counts": dict(counts),
        }
        return result

    try:
        # Lazy import is intentionally after the root-entry marker.  No project
        # module, ancestor helper, provider, or source code is imported.
        import sympy as sp  # noqa: PLC0415  (root-entry gated generic engine)

        result["algebra_entry"]["symbolic_engine_imported"] = True
        result["algebra_entry"]["symbolic_engine"] = "sympy"
        result["algebra_entry"]["symbolic_engine_version"] = getattr(sp, "__version__", None)
        result["algebra_entry"]["symbolic_construction_started"] = True

        ports = _make_formal_ports(sp)
        direct_target, direct_Ln, direct_Lc = _build_direct_source(ports, sp)
        sectors, sector_Ln, sector_Lc, actual_Wn, actual_Wc = _build_sector_graph(ports, sp)
        # This is only a structural denominator-port check; the direct branch
        # still owns its independently constructed expressions.
        if str(direct_Ln) != str(sector_Ln) or str(direct_Lc) != str(sector_Lc):
            raise RuntimeError("direct and sector denominator ports diverged")
        counts["sector_families"] = len(sectors)
        result["formal_ports"] = {
            "common": {"Z": str(ports["Z"]), "d": str(ports["d"])},
            "side_local_symbols": ["Rn", "Rc", "Sn", "Sc", "dn", "dc"],
            "side_local_functions": [
                "Vn(Z)", "Vc(Z)", "mn(Z)", "mc(Z)", "en(Z)", "ec(Z)",
                "P0n(Z)", "P0c(Z)", "pn(Z)", "pc(Z)",
            ],
            "parameters_held_fixed": ["Rn", "Rc", "Sn", "Sc", "dn", "dc"],
            "differentiation": "ordinary raw Sympy.diff in Z; no Taylor factorial normalization",
        }
        result["sector_formulas"] = {
            name: {"text": str(expr), "srepr": sp.srepr(expr)} for name, expr in sectors.items()
        }
        result["direct_source_target"] = {
            "Dsrc": _expr_record(direct_target, sp, label="Dsrc"),
            "formula": "Rn*Lc*Fn_inline - Rc*Ln*Fc_inline",
            "graph": "independent_inline_A_B_L",
            "uses_sector_helpers": False,
        }
        result["actual_W_forms"] = {
            "Wn": "Rn*Fn/Ln",
            "Wc": "Rc*Fc/Lc",
            "cross_numerator_is_not_W_gap_jet": True,
        }
        result["actual_W_expressions"] = {
            "Wn": _expr_record(actual_Wn, sp, label="ActualW_n"),
            "Wc": _expr_record(actual_Wc, sp, label="ActualW_c"),
        }

        sector_jets: dict[int, dict[str, object]] = {}
        direct_jets: dict[int, object] = {}
        for order in ORDERS:
            sector_jets[order] = {}
            for name in SECTOR_NAMES:
                expr = sectors[name]
                counts["symbolic_diff_api_calls"] += 1
                jet = sp.diff(expr, ports["Z"], order)
                row = _jet_row(
                    order=order,
                    name=name,
                    expr=jet,
                    sp=sp,
                    source_refs=source_refs,
                    counts=counts,
                )
                sector_jets[order][name] = jet
                result["numerator_jet_rows"].append(row)
                counts["numerator_jet_rows"] += 1
                row["actual_api_counts"] = dict(counts)
                row["actual_counts"] = dict(counts)
                # Emit the complete row immediately at the row boundary.
                _emit(
                    diagnostic,
                    {
                        "schema": RESULT_SCHEMA,
                        "event": "numerator_jet_row",
                        "payload": deepcopy(row),
                        "prefix": _journal_prefix(result, counts),
                    },
                )

            counts["symbolic_diff_api_calls"] += 1
            direct_jet = sp.diff(direct_target, ports["Z"], order)
            direct_jets[order] = direct_jet
            direct_row = {
                "order": order,
                "expression": _expr_record(direct_jet, sp, label=f"D_Z^{order}[Dsrc]"),
                "source_references": deepcopy(source_refs),
                "actual_api_counts": dict(counts),
                "actual_counts": dict(counts),
            }
            result["direct_source_jet_rows"].append(direct_row)

            # Persist the complete direct-source row before attempting the
            # order's expansion.  If expand fails, this source expression stays
            # in the durable result prefix and is not silently discarded.
            _emit(
                diagnostic,
                {
                    "schema": RESULT_SCHEMA,
                    "event": "direct_source_jet_row",
                    "payload": deepcopy(direct_row),
                    "prefix": _journal_prefix(result, counts),
                },
            )

            counts["symbolic_expand_api_calls"] += 1
            residual = sp.expand(direct_jet - sum(sector_jets[order].values()))
            proof = {
                "order": order,
                "direct_source_jet": _expr_record(direct_jet, sp, label=f"D_Z^{order}[Dsrc]"),
                "sector_sum_jet": _expr_record(sum(sector_jets[order].values()), sp, label=f"sum(D_Z^{order}[D_sector])"),
                "expanded_residual": _expr_record(residual, sp, label=f"expand(Dsrc_{order}-sum_sector_{order})"),
                "identity_zero": bool(residual == 0),
                "ordinary_derivative": True,
                "taylor_coefficient": False,
                "factorial_reapplied": False,
                "source_references": deepcopy(source_refs),
                "actual_api_counts": dict(counts),
            }
            result["proof_rows"].append(proof)
            counts["cross_numerator_identity_orders"] += 1
            proof["actual_api_counts"] = dict(counts)
            proof["actual_counts"] = dict(counts)
            _emit(
                diagnostic,
                {
                    "schema": RESULT_SCHEMA,
                    "event": "cross_numerator_identity_row",
                    "payload": deepcopy(proof),
                    "prefix": _journal_prefix(result, counts),
                },
            )

        # Formal unevaluated quotient views.  No whole-quotient diff calls are
        # made, and no previous quotient recursion is replayed.
        result["whole_W_gap_derivative_views"] = [
            _build_derivative_view(order=order, Z=ports["Z"]) for order in ORDERS
        ]
        counts["whole_W_gap_derivative_views"] = len(result["whole_W_gap_derivative_views"])
        complete = (
            len(result["numerator_jet_rows"]) == 32
            and len(result["proof_rows"]) == 4
            and len(result["whole_W_gap_derivative_views"]) == 4
            and all(row.get("identity_zero") is True for row in result["proof_rows"])
        )
        if not complete:
            raise RuntimeError("ordinary-Z cross-numerator proof prefix is incomplete or nonzero")
        # A mathematical-proof count is earned only after all four explicit
        # zero checks and all required row/view cardinalities pass.
        counts["new_mathematical_proofs"] += 1
        flags = {
            "O3_transition_full_W_signed_sector_functions_defined": True,
            "O3_transition_full_W_cross_numerator_correlation_ordinary_Z3_proved": True,
        }
        _set_flags(result, flags)
        result["status"] = "passed"
        result["passed"] = True
        _set_counts(result, counts)
        return result
    except Exception as exc:
        result["status"] = "failed_after_algebra_entry"
        result["passed"] = False
        _set_flags(
            result,
            {
                "O3_transition_full_W_signed_sector_functions_defined": False,
                "O3_transition_full_W_cross_numerator_correlation_ordinary_Z3_proved": False,
            },
        )
        _set_counts(result, counts)
        result["failure"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "prefix_preserved": True,
            "actual_api_counts": dict(counts),
            "actual_counts": dict(counts),
        }
        return result


__all__ = ["derive", "INPUT_SCHEMA", "RESULT_SCHEMA", "API_SCHEMA"]
