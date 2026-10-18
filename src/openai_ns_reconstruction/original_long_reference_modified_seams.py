"""Source-bound composition of the two modified long/reference seams.

The packet made by :func:`build` is deliberately a conditional theorem.  It
reads the issued C6, C5, and Q4 archives as data and checks the actual field,
history, native-derivative, and pressure ports.  It never imports an owner or
replays an accepted builder/checker/producer.  The two worker packets are
optional because this module must remain honest while either source proof is
still open; a missing or weak packet leaves the full seam claim false.

The source theorem is written in terms of the actual maps in the archives.
In particular, a generic ``C1`` boolean is not accepted in place of the C6
field-to-Q4-radial seed table, the six A/B derivative rows on every chart, or
the actual C5 incoming-memory rows.  The pressure identity is the issued
rational operator with ``L = 1 - delta*Z**2`` and is checked independently by
exact SymPy algebra.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
MODULE = "src/openai_ns_reconstruction/original_long_reference_modified_seams.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
INLET_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"

# These are source data, not executable dependencies.  Keeping the exact
# bytes pinned is important because the worker packets are allowed to refer to
# source-map rows only from this family.
PINNED_HASHES = {
    PARENT: PARENT_SHA256,
    C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    C6_CACHE: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    C5_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    C5_CACHE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    Q4_SOURCE: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    Q4_CACHE: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
    BASE + "lei_ren_part1_paper_compliant_outer_initial.py": "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
    BASE + "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py": "9b41c9a587c52f3af5de35a39695c16e07de509155ce466680d984cabd54b231",
    BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz": "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92",
    BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json": "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
    BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py": "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
    BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py": "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3",
    BASE + "lei_ren_part1_paper_compliant_anchored_axis_amplitude.py": "f3ad62f07b7631ccbfe67c787b4989c0d70db2ff49c96325ca3ab2be284ad87c",
}

SOURCE_FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)
CHARTS = ("post_power", "long_reshape", "reference")
C6_CHART = {"post_power": "post_switch", "long_reshape": "long_reshape", "reference": "reference"}
FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
CORRECTION_FIELDS = ("m", "h", "k", "e", "p")
ORDERS = tuple(range(6))
AXIAL_SCOPE = {
    "Z_lower": "37/100",
    "Z_upper": "93/250",
    "Z_interval": ["37/100", "93/250"],
    "ordinary_Z_orders": [0, 1],
    "source": "accepted inherited source domain; no extension to all axial Z",
    "inherited": True,
    "global_five_moment_closure": False,
}

# The endpoint names and roots are part of the contract with the final graph
# producer.  The chart's finite_boundary rows are also retained below; the
# endpoint rows are the actual upper/lower traces used by that producer.
SEAM_DATA = {
    "post_power_upper_to_long_lower": {
        "upper_chart": "post_power",
        "lower_chart": "long_reshape",
        "upper_endpoint_R_Q_y": [160287, 160288, 160289, 160290],
        "lower_chart_R_Q_y": [166062, 166063, 166064, 166065],
        "graph_lower_wrappers": [185489, 185490],
    },
    "long_upper_to_reference_lower": {
        "upper_chart": "long_reshape",
        "lower_chart": "reference",
        "upper_endpoint_R_Q_y": [166149, 166150, 166151, 166152],
        "lower_chart_R_Q_y": [171782, 171783, 171784, 171785],
        "graph_lower_wrappers": [185672, 185673],
    },
}

SOURCE_GUARDS = (
    (C6_SOURCE, "__init__", "CurrentOriginalPrescribedR100RmC6Functions", (
        "self.add_chart('long_reshape'",
        "self.reference_length=length=g.sub(g.sub(logref,T),c(8))",
        "self.add_chart('reference'",
    )),
    (C5_SOURCE, "loop_functions", None, (
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
        "Aphi=lazy(a.scale(a.mul(roots['a'],a.sub(one,a.scale(invlambda,g.quotient(g.one,twopi,'positive 2*pi')))),'1/2'))",
        "Bphi=lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),g.quotient(g.one,twopi,'positive 2*pi'))),roots['b'])),'1/2'))",
        "partials={(k,0):g.mul(c(Fraction(1,math.factorial(k))),rows(Phi)[k]) for k in range(6)}",
        "partials[k,j]=g.mul(c(Fraction(1,j)),lam[k,j-1])",
    )),
    (C5_SOURCE, "finite_N_functions", None, (
        "x=a.scale(loop['A'],invN);ex=a.exp(x)",
        "s=a.scale(loop['B_over_Pstar'],invN)",
        "inc=dict(m=s",
        "return dict(E_N=a.mul(E,ex),V_N=a.add(V,s)",
        "five_signed_density_increments=inc",
    )),
    (C5_SOURCE, "transport", None, (
        "predecessor_memory_never_reset=True",
        "zero_correction_inlet_only_at_original_rminus=True",
        "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
    )),
    (Q4_SOURCE, "__init__", "NativeDerivative", (
        "self.cache={phase.node:g.zero};self.core_extensions=[]",
        "self.shared['factories'][variable.node]=self",
    )),
    (Q4_SOURCE, "seed", "NativeDerivative", (
        "self.cache[left.node]=right",
        "Aconstantcannot have a nonzero radial derivative",
    )),
    (Q4_SOURCE, "sub", "NativeDerivative", (
        "function_substitution",
        "Z_independent_substitution=True",
        "exact_native_chain_rule_substitution=True",
    )),
    (Q4_SOURCE, "derivative", "NativeDerivative", (
        "exact_native_partial_of_same_integrand",
        "original positive phase lambda; exact implicit native derivative",
        "same_original_lazy_native_derivative",
    )),
    (Q4_SOURCE, "total_radial_functions", None, (
        "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))",
        "xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
        "physical_y_equals_native_partial_divided_by_Jacobian=True",
    )),
    (Q4_SOURCE, "physical_pressure_fields", None, (
        "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
        "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
        "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
    )),
)


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _locate(root: Path, relative: str) -> Path:
    candidate = root / relative
    if candidate.is_file():
        return candidate
    # A focused review can pass the experiments directory itself.  Preserve
    # the same repo-relative names in the returned pin records.
    candidate = root / Path(relative).name
    if candidate.is_file():
        return candidate
    for found in root.rglob(Path(relative).name):
        if found.is_file():
            return found
    raise FileNotFoundError(relative)


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def _gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path}")
    return value


def _compact(value: str) -> str:
    return "".join(value.split())


def _ast_target(tree: ast.AST, name: str, owner: str | None) -> ast.AST | None:
    if owner is None:
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
                return node
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == owner:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                    return child
    return None


def _source_guard(root: Path, path_name: str, name: str, owner: str | None,
                  snippets: tuple[str, ...]) -> dict[str, Any]:
    path = _locate(root, path_name)
    raw = path.read_bytes()
    _require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
             f"strict LF single EOF guard failed: {path_name}")
    source = raw.decode("utf-8")
    tree = ast.parse(source, filename=str(path))
    node = _ast_target(tree, name, owner)
    _require(node is not None, f"AST target missing: {path_name}:{owner or ''}{name}")
    segment = ast.get_source_segment(source, node) or ""
    compact = _compact(segment)
    # A small typo in an error-only snippet must fail closed, rather than turn
    # an absent source rule into a passing packet.
    matched = all(_compact(token) in compact for token in snippets)
    _require(matched, f"source guard failed: {path_name}:{owner or ''}{name}")
    return {
        "path": path_name,
        "target": name,
        "owner": owner,
        "tokens": snippets,
        "token_count": len(snippets),
        "AST_sha256": hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        "passed": True,
        "source_executed": False,
    }


def _source_guards(root: Path) -> tuple[dict[str, Any], ...]:
    guards = tuple(_source_guard(root, *spec) for spec in SOURCE_GUARDS)
    _require(isinstance(guards, tuple), "source guards must remain a tuple")
    return guards


def _pin(root: Path, relative: str, expected: str | None) -> dict[str, Any]:
    path = _locate(root, relative)
    actual = _sha(path)
    expected_value = actual if expected is None else expected
    _require(actual == expected_value, f"source pin changed: {relative}")
    return {
        "path": relative,
        "expected_sha256": expected_value,
        "actual_sha256": actual,
        "passed": actual == expected_value,
        "source_executed": False,
    }


def _source_pins(root: Path) -> dict[str, Any]:
    pins = {relative: _pin(root, relative, expected) for relative, expected in PINNED_HASHES.items()}
    module_pin = _pin(root, MODULE, None)
    pins[MODULE] = module_pin
    _require(all(pin["passed"] and not pin["source_executed"] for pin in pins.values()),
             "source pins must all pass without execution")
    return pins


def _load_parent(root: Path, parent: Mapping[str, Any] | str | Path | None) -> tuple[Path, Mapping[str, Any]]:
    path = _locate(root, PARENT)
    disk = _json(path)
    if parent is None:
        supplied: Mapping[str, Any] = disk
    elif isinstance(parent, (str, Path)):
        supplied = _json(_locate(root, str(parent)))
    else:
        supplied = parent
    # The accepted report contains a few legacy non-finite diagnostic values;
    # ordinary mapping equality treats two independently decoded NaNs as
    # unequal even when the bytes and every serialised field agree.
    supplied_dump = json.dumps(supplied, sort_keys=True, separators=(",", ":"), allow_nan=True)
    disk_dump = json.dumps(disk, sort_keys=True, separators=(",", ":"), allow_nan=True)
    if supplied_dump != disk_dump:
        # Some Windows Python callers decode the report through the console
        # code page while this module reads UTF-8.  The only observed
        # discrepancy is a redundant U+00C2 lead marker before a UTF-8
        # superscript; normalize that transport artefact while retaining the
        # full key/value comparison.
        def normalize(value: Any) -> Any:
            if isinstance(value, str):
                return value.replace("\u00c2", "")
            if isinstance(value, list):
                return [normalize(item) for item in value]
            if isinstance(value, Mapping):
                return {key: normalize(item) for key, item in value.items()}
            return value
        supplied_dump = json.dumps(normalize(supplied), sort_keys=True, separators=(",", ":"), allow_nan=True)
        disk_dump = json.dumps(normalize(disk), sort_keys=True, separators=(",", ":"), allow_nan=True)
    _require(supplied_dump == disk_dump, "the exact accepted current seam-prerequisite parent is required")
    _require(_sha(path) == PARENT_SHA256, "accepted current parent bytes changed")
    _require(isinstance(disk.get("source_family"), Mapping), "accepted source family missing")
    for section in ("actual_long_reference_native_seams", "actual_long_reference_history_seams",
                    "actual_long_reference_loop_seams"):
        _require(isinstance(disk.get(section), Mapping), f"accepted parent section missing: {section}")
    return path, disk


def _copy(value: Any) -> Any:
    # JSON round-tripping keeps receipts detached from mutable archive maps and
    # also rejects accidental non-serializable source objects.
    return json.loads(json.dumps(value, allow_nan=False))


def _load_archives(root: Path, family: Mapping[str, Any]) -> dict[str, Any]:
    c6 = _gzip_json(_locate(root, C6_CACHE))
    c5 = _gzip_json(_locate(root, C5_CACHE))
    q4 = _gzip_json(_locate(root, Q4_CACHE))
    for packet, key, length, label in (
        (c6, "prescribed_R100_Rm_C6_bindings", 97295, "C6"),
        (c5, "prescribed_loop_C5_bindings", 121834, "C5"),
        (q4, "modified_Q4_pressure_C3_bindings", 183700, "Q4"),
    ):
        _require(packet.get("source_family") == family, f"{label} source family changed")
        nodes = packet.get("exact_graph_nodes")
        _require(isinstance(nodes, list) and len(nodes) == length, f"{label} graph node length changed")
        _require(isinstance(packet.get(key), Mapping), f"{label} source binding missing")
    _require(q4["exact_graph_nodes"][:len(c5["exact_graph_nodes"])] == c5["exact_graph_nodes"],
             "Q4/C5 graph prefix changed")
    _require(c5["exact_graph_nodes"][:len(c6["exact_graph_nodes"])] == c6["exact_graph_nodes"],
             "C5/C6 graph prefix changed")
    return {"c6": c6, "c5": c5, "q4": q4}


def _node(nodes: Sequence[Any], index: int) -> Mapping[str, Any]:
    _require(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(nodes),
             f"invalid graph node {index!r}")
    value = nodes[index]
    _require(isinstance(value, Mapping), f"graph node {index} is not an object")
    return value


def _six_rows(value: Any, label: str) -> list[int]:
    _require(isinstance(value, list) and len(value) == 6 and all(isinstance(x, int) and not isinstance(x, bool) for x in value),
             f"six ordinary rows required: {label}")
    return [int(x) for x in value]


def _field_rows(value: Any, label: str) -> list[int]:
    # C6 source fields and Q4 radial backgrounds retain the seven rows
    # through ordinary order six.  The A/B loop and physical U rows below are
    # six-row jets because the physical-y derivative consumes one order.
    _require(isinstance(value, list) and len(value) == 7 and all(isinstance(x, int) and not isinstance(x, bool) for x in value),
             f"seven source field rows required: {label}")
    return [int(x) for x in value]


def _endpoint_rows(value: Any, label: str) -> list[int]:
    _require(isinstance(value, list) and len(value) == 4 and all(isinstance(x, int) and not isinstance(x, bool) for x in value),
             f"four endpoint rows required: {label}")
    return [int(x) for x in value]


def _five_rows(value: Any, label: str) -> list[int]:
    _require(isinstance(value, list) and len(value) == 5 and all(isinstance(x, int) and not isinstance(x, bool) for x in value),
             f"five pressure rows required: {label}")
    return [int(x) for x in value]


def _same_rows(left: Mapping[str, Any], right: Mapping[str, Any], fields: Sequence[str]) -> bool:
    return all(left.get(field) == right.get(field) for field in fields)


def _graph_port_receipt(archives: Mapping[str, Any]) -> dict[str, Any]:
    c6 = archives["c6"]["prescribed_R100_Rm_C6_bindings"]
    c5 = archives["c5"]["prescribed_loop_C5_bindings"]
    q4 = archives["q4"]["modified_Q4_pressure_C3_bindings"]
    c6charts, c5charts, q4charts = c6["charts"], c5["charts"], q4["charts"]
    c6nodes = archives["c6"]["exact_graph_nodes"]
    c5nodes = archives["c5"]["exact_graph_nodes"]
    q4nodes = archives["q4"]["exact_graph_nodes"]
    pstar = int(c6["Pstar"])
    chart_receipts: dict[str, Any] = {}
    all_ports = True
    all_histories = True
    all_phase_rows = True
    seed_rows: dict[str, Any] = {}
    generic_rows: dict[str, Any] = {}
    loop_source_receipts: dict[str, Any] = {}
    for chart in CHARTS:
        c6chart = c6charts[C6_CHART[chart]]
        c5chart = c5charts[chart]
        q4chart = q4charts[chart]
        _require(c5chart["physical_measure"] == c6chart["physical_measure"], f"J mismatch: {chart}")
        chart_N = int(c5chart["N"]) if "N" in c5chart else int(c5["N"])
        _require(chart_N == int(c5["N"]) == 18, f"N mismatch: {chart}")
        J = int(c6chart["physical_measure"])
        c6fields = {name: _field_rows(c6chart["fields"][name], f"{chart}.C6.{name}") for name in FIELDS}
        q4bg = {name: _field_rows(q4chart["radial_background"][name], f"{chart}.Q4.{name}") for name in FIELDS}
        seed_rows[chart] = {
            field: [{"source_field_node": source, "physical_radial_row_node": physical,
                     "ordinary_order": order, "physical_J_node": J}
                    for order, (source, physical) in enumerate(zip(c6fields[field], q4bg[field]))]
            for field in FIELDS
        }
        generic_rows[chart] = {
            "C6": {key: _six_rows(value, f"{chart}.C6.generic.{key}")
                   for key, value in c6chart["generic"].items() if isinstance(value, list)},
            "Q4": {key: _five_rows(value, f"{chart}.Q4.generic.{key}")
                   for key, value in q4chart["generic"].items() if isinstance(value, list)},
        }
        loop = c5chart["loop"]
        for name in ("A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi"):
            _six_rows(loop[name], f"{chart}.loop.{name}")
        phase_rows = loop.get("phase_partial_coefficients")
        if not isinstance(phase_rows, list) or len(phase_rows) != 21:
            all_phase_rows = False
        inverse_rows = _six_rows(loop["inverse"], f"{chart}.loop.inverse")
        phase_lambda_rows = _six_rows(loop["phase_lambda_at_inverse"], f"{chart}.loop.phase_lambda_at_inverse")
        aphi_rows = _six_rows(loop["A_phi"], f"{chart}.loop.A_phi")
        bphi_rows = _six_rows(loop["B_over_Pstar_phi"], f"{chart}.loop.B_over_Pstar_phi")
        inverse_checks = []
        for inverse_order, node_id in enumerate(inverse_rows):
            row = _node(c5nodes, node_id)
            inverse_checks.append({
                "node": node_id,
                "operation": row.get("operation"),
                "phase_function": row.get("phase_function"),
                "target_phase": row.get("target_phase"),
                "angle_endpoints": row.get("phase_endpoints"),
                "positive_certificate": row.get("derivative_positive_certificate"),
                "actual_original_inverse": row.get("branch_inverse_is_original_not_support_cap") is True,
                # The first row is the actual inverse root.  The remaining
                # rows are its ordinary-Z jet derivatives and intentionally
                # have product/negative operations.
                "passed": (row.get("operation") == "original_monotone_phase_inverse"
                           and row.get("phase_endpoints") == [0, 1]
                           and row.get("branch_inverse_is_original_not_support_cap") is True)
                if inverse_order == 0 else True,
            })
        lazy_checks = []
        for kind, rows in (("A_phi", aphi_rows), ("B_over_Pstar_phi", bphi_rows)):
            for node_id in rows:
                row = _node(c5nodes, node_id)
                lazy_checks.append({
                    "kind": kind,
                    "node": node_id,
                    "operation": row.get("operation"),
                    "Delta": row.get("Delta"),
                    "eta": row.get("eta"),
                    "flat_value": row.get("flat_value"),
                    "flat_predicate": row.get("flat_predicate"),
                    "passed": row.get("operation") == "original_lazy_flat_branch"
                    and row.get("flat_value") == 0
                    and row.get("active_body_not_evaluated_on_flat_branch") is True,
                })
        loop_source_receipts[chart] = {
            "actual_inverse_roots": inverse_rows,
            "actual_phase_lambda_roots": phase_lambda_rows,
            "actual_A_phi_roots": aphi_rows,
            "actual_B_over_Pstar_phi_roots": bphi_rows,
            "inverse_checks": inverse_checks,
            "lazy_flat_checks": lazy_checks,
            "actual_inverse_target_and_positive_lambda_bound": all(row["passed"] for row in inverse_checks),
            "actual_lazy_flat_metadata_bound": all(row["passed"] for row in lazy_checks),
            "source_formula_status": "conditional on positive phase controls and IFT source theorem",
        }
        finite = q4chart["finite_boundary"]
        endpoint = q4["endpoints"][chart]["boundary"]
        finite_rqy = _endpoint_rows(finite["R_Q_y"], f"{chart}.finite_boundary.R_Q_y")
        endpoint_rqy = _endpoint_rows(endpoint["R_Q_y"], f"{chart}.endpoint.R_Q_y")
        physical = q4chart["physical"]
        radial = q4chart["radial"]
        physical_rows = {key: _six_rows(physical[key], f"{chart}.physical.{key}")
                         for key in ("U", "U_y", "m", "m_y")}
        physical_rows["Q_y"] = _five_rows(physical["Q_y"], f"{chart}.physical.Q_y")[:]
        physical_checks: list[dict[str, Any]] = []
        for order, (u, v, m, my) in enumerate(zip(physical_rows["U"], physical_rows["U_y"], physical_rows["m"], physical_rows["m_y"])):
            urow, vrow, mrow, myrow = _node(q4nodes, u), _node(q4nodes, v), _node(q4nodes, m), _node(q4nodes, my)
            expected = {
                "U": urow.get("operation") == "product" and urow.get("arguments") == [c5chart["modified_fields"]["V"][order], pstar],
                "U_y": vrow.get("operation") == "product" and vrow.get("arguments") == [radial["V_y"][order], pstar],
                "m": mrow.get("operation") == "product" and mrow.get("arguments") == [c5chart["modified_fields"]["m"][order], pstar],
                "m_y": myrow.get("operation") == "sum" and myrow.get("arguments") == [u, myrow.get("arguments", [None, None])[1]],
            }
            # The second summand is the issued negative m root.  Check its
            # operation and source explicitly rather than accepting any sum.
            if expected["m_y"]:
                neg = _node(q4nodes, int(myrow["arguments"][1]))
                expected["m_y"] = neg.get("operation") == "negative" and neg.get("argument") == m
            physical_checks.append({"order": order, "passed": all(expected.values()), "checks": expected})
        qy_checks: list[dict[str, Any]] = []
        for order, qy in enumerate(physical_rows["Q_y"]):
            row = _node(q4nodes, qy)
            qy_checks.append({
                "order": order,
                "node": qy,
                "operation": row.get("operation"),
                "denominator": row.get("denominator"),
                "passed": row.get("operation") == "positive_quotient" and row.get("denominator") == 25272
                and row.get("source_positive_certificate") == "same original L=1-delta*Z^2",
            })
        # Bind the finite physical radius times Q_y rows and the endpoint
        # substitution rows separately.  The latter are the source endpoint
        # theorem used by the parent graph patch; a bare Q_y root is not an
        # R_Q_y trace.
        finite_rqy_checks: list[dict[str, Any]] = []
        endpoint_rqy_checks: list[dict[str, Any]] = []
        coordinate = int(c5chart["coordinate"])
        for order, (rqy, qy) in enumerate(zip(finite_rqy, physical_rows["Q_y"])):
            row = _node(q4nodes, rqy)
            check = (row.get("operation") == "product"
                     and row.get("arguments") == [qy, int(c5chart["radius"])])
            finite_rqy_checks.append({"order": order, "node": rqy,
                                      "Q_y": qy, "radius": int(c5chart["radius"]),
                                      "passed": check})
        for order, (endpoint_row, finite_row) in enumerate(zip(endpoint_rqy, finite_rqy)):
            row = _node(q4nodes, endpoint_row)
            check = (row.get("operation") == "function_substitution"
                     and row.get("expression") == finite_row
                     and row.get("variable") == coordinate
                     and row.get("value") == 1
                     and row.get("Z_independent_substitution") is True)
            endpoint_rqy_checks.append({"order": order, "node": endpoint_row,
                                        "expression": finite_row, "coordinate": coordinate,
                                        "value": 1, "passed": check})
        port_passed = (all(row["passed"] for row in physical_checks)
                       and all(row["passed"] for row in qy_checks)
                       and all(row["passed"] for row in finite_rqy_checks)
                       and all(row["passed"] for row in endpoint_rqy_checks))
        all_ports = all_ports and port_passed and radial.get("total_fast_phase_terms_included") is True and radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True
        # Every correction history is non-reset and the next chart's incoming
        # rows are the preceding endpoint rows.  This is pointer identity, not
        # a numerical value selection.
        if chart != "post_power":
            previous = "long_reshape" if chart == "reference" else "post_power"
            incoming = c5chart["incoming"]
            prior = c5["endpoints"][previous]
            same = all(incoming[name] == prior[name] for name in CORRECTION_FIELDS)
            all_histories = all_histories and same
        chart_receipts[chart] = {
            "C6_chart": C6_CHART[chart],
            "physical_J_node": J,
            "C5_coordinate_node": int(c5chart["coordinate"]),
            "C5_actual_phase_node": int(c5chart["actual_phase"]),
            "C5_global_offset_node": int(c5chart["global_offset"]),
            "common_N_graph_handle": int(c5["N"]),
            "seed_rows": seed_rows[chart],
            "generic_rows": generic_rows[chart],
            "phase_partial_coefficients": _copy(phase_rows) if isinstance(phase_rows, list) else phase_rows,
            "loop_source": loop_source_receipts[chart],
            "physical_ports": physical_rows,
            "physical_port_checks": physical_checks,
            "Q_y_port_checks": qy_checks,
            "finite_R_Q_y_checks": finite_rqy_checks,
            "endpoint_R_Q_y_checks": endpoint_rqy_checks,
            "finite_boundary_R_Q_y": finite_rqy,
            "endpoint_R_Q_y": endpoint_rqy,
            "history": {
                "incoming": _copy(c5chart["incoming"]),
                "endpoints": _copy(c5["endpoints"][chart]),
                "memory": _copy(c5chart["memory"]),
                "predecessor_memory_never_reset": True,
                "zero_rate_p_root": c5chart["memory"].get("p") == 1,
            },
            "passed": bool(port_passed and all_phase_rows),
        }
    endpoint_map = {chart: chart_receipts[chart]["endpoint_R_Q_y"] for chart in CHARTS}
    rqy_traces = {
        name: {
            "upper_chart": data["upper_chart"],
            "lower_chart": data["lower_chart"],
            "upper_endpoint_R_Q_y": _copy(endpoint_map[data["upper_chart"]]),
            "lower_chart_R_Q_y": _copy(chart_receipts[data["lower_chart"]]["finite_boundary_R_Q_y"]),
            "expected_upper_endpoint_R_Q_y": list(data["upper_endpoint_R_Q_y"]),
            "expected_lower_chart_R_Q_y": list(data["lower_chart_R_Q_y"]),
            "graph_lower_wrappers": list(data["graph_lower_wrappers"]),
            "cache_endpoint_trace_matches_contract": endpoint_map[data["upper_chart"]] == data["upper_endpoint_R_Q_y"]
            and chart_receipts[data["lower_chart"]]["finite_boundary_R_Q_y"] == data["lower_chart_R_Q_y"],
        }
        for name, data in SEAM_DATA.items()
    }
    return {
        "charts": chart_receipts,
        "seed_rows": seed_rows,
        "generic_rows": generic_rows,
        "loop_source": loop_source_receipts,
        "source_formula_identities": {
            "A_phi": "a/2*(1-1/(2*pi*lambda))",
            "B_phi": "E/2*(-a*t(psi)/(2*pi*lambda)-b)",
            "phase_inverse": "Phi(psi)=target_phase with lambda=partial_psi Phi and Phi_psi>0",
            "implicit_function_theorem": "actual inverse rows use the original monotone phase inverse and endpoint phases [0,1]",
            "lazy_flat_extension": "A_phi and B_phi retain Delta/eta/flat_value/predicate and defer active body off the flat branch",
            "universal_functional_map_binding_required": True,
            "physical_chain_rule": "D_y = J^{-1} D_native on actual Q4 slow_A_y/slow_B_y roots; fast phase terms are added explicitly",
            "source_executed": False,
        },
        "physical_port_roots": {chart: _copy(chart_receipts[chart]["physical_ports"]) for chart in CHARTS},
        "rqy_endpoint_traces": rqy_traces,
        "graph_prefix_lengths": {"C6": 97295, "C5": 121834, "Q4": 183700},
        "common_N_graph_handle": int(c5["N"]),
        "all_physical_ports_bound": all_ports,
        "all_R_Q_y_source_rows_bound": all(
            all(row["passed"] for row in chart_receipts[chart]["finite_R_Q_y_checks"])
            and all(row["passed"] for row in chart_receipts[chart]["endpoint_R_Q_y_checks"])
            for chart in CHARTS),
        "all_nonreset_history_pointers_bound": all_histories,
        "all_actual_phase_partial_rows_present": all_phase_rows,
        "all_endpoint_traces_match_contract": all(row["cache_endpoint_trace_matches_contract"] for row in rqy_traces.values()),
        "passed": bool(all_ports and all_histories and all_phase_rows
                       and all(row["cache_endpoint_trace_matches_contract"] for row in rqy_traces.values())
                       and all(all(item["passed"] for item in chart_receipts[chart]["finite_R_Q_y_checks"])
                               and all(item["passed"] for item in chart_receipts[chart]["endpoint_R_Q_y_checks"])
                               for chart in CHARTS)),
    }


def _qy_denominator_receipt(archives: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    """Consume the saved positive ``L=1-delta*Z^2`` source theorem.

    The Q4 graph carries one canonical denominator row.  This receipt binds
    every physical Q_y chart to that row and consumes the inherited
    long/reference tail inequalities for the same source window.  It only
    checks saved booleans and graph topology; it does not recompute the
    accepted interval arithmetic.
    """

    q4 = archives["q4"]["modified_Q4_pressure_C3_bindings"]
    nodes = archives["q4"]["exact_graph_nodes"]
    denominator_id = 25272
    canonical_delta_id = 25220
    canonical_z_squared_id = 25249
    denominator = _node(nodes, denominator_id)
    negative = _node(nodes, 25271)
    product = _node(nodes, 25270)
    z_squared = _node(nodes, canonical_z_squared_id)
    delta = _node(nodes, canonical_delta_id)
    topology_checks = {
        "denominator_row": denominator.get("operation") == "sum"
        and denominator.get("arguments") == [1, 25271],
        "negative_delta_z_squared": negative.get("operation") == "negative"
        and negative.get("argument") == 25270,
        "canonical_delta_product": product.get("operation") == "product"
        and product.get("arguments") == [canonical_delta_id, canonical_z_squared_id],
        "canonical_Z_squared": z_squared.get("operation") == "product"
        and z_squared.get("arguments") == [20, 20],
        "canonical_delta_source_row": delta.get("operation") == "analytic_unary"
        and delta.get("name") == "exp",
    }
    chart_rows: dict[str, Any] = {}
    graph_rows_ok = True
    for chart in CHARTS:
        qy_rows = q4["charts"][chart]["physical"]["Q_y"]
        checks = []
        for order, node_id in enumerate(qy_rows):
            row = _node(nodes, node_id)
            check = {
                "order": order,
                "node": node_id,
                "denominator": row.get("denominator"),
                "same_canonical_denominator": row.get("denominator") == denominator_id,
                "operation": row.get("operation") == "positive_quotient",
                "source_positive_certificate": row.get("source_positive_certificate")
                == "same original L=1-delta*Z^2",
            }
            check["passed"] = all(value for key, value in check.items()
                                   if key not in ("order", "node", "denominator"))
            checks.append(check)
        chart_rows[chart] = {"Q_y_rows": checks,
                             "passed": all(item["passed"] for item in checks)}
        graph_rows_ok = graph_rows_ok and chart_rows[chart]["passed"]

    background = parent.get("actual_long_reference_background_I1_numeric")
    background_charts = background.get("charts") if isinstance(background, Mapping) else None
    saved_rows: dict[str, Any] = {}
    saved_positive = True
    for chart in ("long_reshape", "reference"):
        section = background_charts.get(chart) if isinstance(background_charts, Mapping) else None
        tail = section.get("tail_proof") if isinstance(section, Mapping) else None
        inequalities = tail.get("inequalities") if isinstance(tail, Mapping) else None
        checks = {
            "chart_present": isinstance(section, Mapping),
            "chart_passed": isinstance(section, Mapping) and section.get("passed") is True,
            "source_executed_false": isinstance(section, Mapping)
            and section.get("source_executed") is False,
            "tail_source_executed_false": isinstance(tail, Mapping)
            and tail.get("source_executed") is False,
            "L_positive": isinstance(inequalities, Mapping)
            and inequalities.get("L_positive") is True,
        }
        saved_rows[chart] = {"checks": checks, "passed": all(checks.values())}
        saved_positive = saved_positive and saved_rows[chart]["passed"]
    directions = background.get("directions") if isinstance(background, Mapping) else None
    inherited_window = (isinstance(directions, Mapping)
                         and directions.get("Z_domain") == "[37/100,93/250]")
    passed = bool(all(topology_checks.values()) and graph_rows_ok and saved_positive
                  and inherited_window)
    return {
        "passed": passed,
        "denominator_node": denominator_id,
        "canonical_delta_node": canonical_delta_id,
        "canonical_Z_squared_node": canonical_z_squared_id,
        "topology_checks": topology_checks,
        "graph_rows": chart_rows,
        "saved_background_rows": saved_rows,
        "inherited_Z_window": "[37/100,93/250]",
        "inherited_window": inherited_window,
        "source_formula": "L=1-delta*Z^2",
        "source_executed": False,
    }


def _exact_rational_row(nodes: Sequence[Any], node_id: Any, numerator: int,
                        denominator: int = 1) -> bool:
    if not isinstance(node_id, int) or not (0 <= node_id < len(nodes)):
        return False
    row = nodes[node_id]
    return (isinstance(row, Mapping) and row.get("operation") == "exact_rational"
            and row.get("numerator") == numerator and row.get("denominator") == denominator)


def _product_factors(nodes: Sequence[Any], node_id: Any) -> list[int]:
    """Flatten only issued product nodes, preserving every actual handle."""

    if not isinstance(node_id, int) or not (0 <= node_id < len(nodes)):
        return []
    row = nodes[node_id]
    if not isinstance(row, Mapping):
        return []
    if row.get("operation") != "product":
        return [node_id]
    factors: list[int] = []
    for child in row.get("arguments", []):
        factors.extend(_product_factors(nodes, child))
    return factors


def _descendants(nodes: Sequence[Any], node_id: Any) -> list[int]:
    """Return actual integer graph handles reachable through operation edges."""

    out: list[int] = []
    seen: set[int] = set()
    pending = [node_id]
    while pending:
        current = pending.pop()
        if not isinstance(current, int) or current in seen or not (0 <= current < len(nodes)):
            continue
        seen.add(current)
        out.append(current)
        row = nodes[current]
        if not isinstance(row, Mapping):
            continue
        for key in ("active_body", "argument", "expression", "numerator", "denominator",
                    "body", "inverse_angle"):
            child = row.get(key)
            if isinstance(child, int):
                pending.append(child)
        for child in row.get("arguments", []):
            if isinstance(child, int):
                pending.append(child)
    return out


def _positive_inverse_row(nodes: Sequence[Any], node_id: Any, denominator: int,
                          certificate_fragment: str | None = None) -> bool:
    if not isinstance(node_id, int) or not (0 <= node_id < len(nodes)):
        return False
    row = nodes[node_id]
    if not isinstance(row, Mapping):
        return False
    certificate = str(row.get("source_positive_certificate", ""))
    return (row.get("operation") == "positive_quotient"
            and row.get("numerator") == 1
            and row.get("denominator") == denominator
            and (certificate_fragment is None or certificate_fragment in certificate))


def _formula_core(nodes: Sequence[Any], root_id: int, lead_id: int,
                  inverse_lambda_id: int, inverse_twopi_id: int,
                  kind: str) -> dict[str, Any]:
    """Check an issued A_phi/B_phi graph against its actual source formula.

    This is a structural check over the saved C5 DAG.  It does not evaluate
    any source function, and it deliberately keeps the actual node handles in
    the receipt so a later reviewer can distinguish the universal formula
    from a free-standing equality flag.
    """

    root = nodes[root_id] if isinstance(root_id, int) and root_id < len(nodes) else {}
    passed = isinstance(root, Mapping)
    active_id = root.get("active_body") if isinstance(root, Mapping) else None
    active = nodes[active_id] if isinstance(active_id, int) and active_id < len(nodes) else {}
    args = list(active.get("arguments", [])) if isinstance(active, Mapping) else []
    half_ids = [value for value in args if _exact_rational_row(nodes, value, 1, 2)]
    core_ids = [value for value in args if value not in half_ids]
    passed = passed and root.get("operation") == "original_lazy_flat_branch" and len(half_ids) == 1 and len(core_ids) == 1
    core_id = core_ids[0] if core_ids else None
    core = nodes[core_id] if isinstance(core_id, int) and core_id < len(nodes) else {}
    core_args = list(core.get("arguments", [])) if isinstance(core, Mapping) else []
    if kind == "A_phi":
        lead_ok = lead_id in core_args
        residual_ids = [value for value in core_args if value != lead_id]
        residual = nodes[residual_ids[0]] if len(residual_ids) == 1 and isinstance(residual_ids[0], int) else {}
        residual_args = list(residual.get("arguments", [])) if isinstance(residual, Mapping) else []
        one_ids = [value for value in residual_args if _exact_rational_row(nodes, value, 1)]
        neg_ids = [value for value in residual_args if isinstance(value, int) and isinstance(nodes[value], Mapping)
                   and nodes[value].get("operation") == "negative"]
        inv_product_id = nodes[neg_ids[0]].get("argument") if len(neg_ids) == 1 else None
        factor_ids = _product_factors(nodes, inv_product_id)
        inv_lambda_ok = inverse_lambda_id in factor_ids and _positive_inverse_row(nodes, inverse_lambda_id,
                                                                                    nodes[inverse_lambda_id].get("denominator") if isinstance(nodes[inverse_lambda_id], Mapping) else -1,
                                                                                    "lambda at true inverse")
        inv_twopi_ok = any(_positive_inverse_row(nodes, value,
                                                  nodes[value].get("denominator") if isinstance(nodes[value], Mapping) else -1,
                                                  "positive 2*pi") for value in factor_ids)
        passed = (passed and core.get("operation") == "product" and lead_ok
                  and len(one_ids) == 1 and len(neg_ids) == 1
                  and isinstance(inv_product_id, int) and inv_lambda_ok and inv_twopi_ok)
        formula = "a/2*(1-1/(2*pi*lambda))"
        details = {"lead_handle": lead_id, "one_handle": one_ids[0] if one_ids else None,
                   "negative_inverse_product": inv_product_id, "factors": factor_ids,
                   "inverse_lambda_factor": inverse_lambda_id, "passed": passed}
    else:
        lead_ok = lead_id in core_args
        bracket_ids = [value for value in core_args if value != lead_id]
        bracket_id = bracket_ids[0] if len(bracket_ids) == 1 else None
        bracket = nodes[bracket_id] if isinstance(bracket_id, int) and bracket_id < len(nodes) else {}
        bracket_terms = list(bracket.get("arguments", [])) if isinstance(bracket, Mapping) else []
        if isinstance(bracket, Mapping) and bracket.get("operation") == "negative":
            bracket_terms = [bracket_id]
        negative_term_ids = []
        b_term_ids = []
        for value in bracket_terms:
            row = nodes[value] if isinstance(value, int) and value < len(nodes) else {}
            if isinstance(row, Mapping) and row.get("operation") == "negative":
                negative_term_ids.append(row.get("argument"))
            else:
                b_term_ids.append(value)
        term_factors = _product_factors(nodes, negative_term_ids[0]) if negative_term_ids else []
        inv_lambda_ok = inverse_lambda_id in term_factors
        inv_twopi_ok = any(_positive_inverse_row(nodes, value,
                                                  nodes[value].get("denominator") if isinstance(nodes[value], Mapping) else -1,
                                                  "positive 2*pi") for value in term_factors)
        passed = (passed and core.get("operation") == "product" and lead_ok
                  and len(bracket_ids) == 1 and len(negative_term_ids) == 1
                  and inv_lambda_ok and inv_twopi_ok and not b_term_ids)
        formula = "E/2*(-a*t(psi)/(2*pi*lambda)-b)"
        details = {"lead_handle": lead_id, "bracket_handle": bracket_id,
                   "negative_term": negative_term_ids[0] if negative_term_ids else None,
                   "term_factors": term_factors, "inverse_lambda_factor": inverse_lambda_id,
                   "b_zero_omitted": not b_term_ids, "passed": passed}
    return {"formula": formula, "root": root_id, "active_body": active_id,
            "active_half_factor": half_ids[0] if half_ids else None,
            "details": details, "passed": bool(passed)}


def _mixed_lambda_source_receipt(nodes: Sequence[Any], loop: Mapping[str, Any],
                                 roots: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the issued normalized MixedTaylor phase source graph.

    The phase partial table is emitted by the C5 source before the inverse
    substitution.  A row count and a positive ``K`` descendant do not show
    that the row is the coefficient generated by
    ``MixedTaylor(K)*(1+t**2)``.  This receipt checks the actual graph
    convolutions, the source ``t = alpha*(cos-r)/(1-2*r*cos+r**2)`` rows,
    and the common denominator handles.  It only reads graph data.
    """

    def valid_id(value: Any) -> bool:
        return isinstance(value, int) and 0 <= value < len(nodes)

    def row(value: Any) -> Mapping[str, Any]:
        return nodes[value] if valid_id(value) and isinstance(nodes[value], Mapping) else {}

    def product(value: Any) -> bool:
        return row(value).get("operation") == "product"

    def terms(value: Any) -> list[int]:
        current = row(value)
        if current.get("operation") == "sum":
            return [item for item in current.get("arguments", []) if isinstance(item, int)]
        return [value] if valid_id(value) else []

    def exact(value: Any, numerator: int, denominator: int = 1) -> bool:
        return _exact_rational_row(nodes, value, numerator, denominator)

    def positive_quotient(value: Any, source_token: str) -> bool:
        current = row(value)
        return (current.get("operation") == "positive_quotient"
                and source_token in str(current.get("source_positive_certificate", "")))

    expected = [(k, j) for k in range(5) for j in range(5 - k)]
    partial_rows = loop.get("phase_partial_coefficients")
    indexed = {(item.get("Z"), item.get("angle")): item.get("node")
               for item in partial_rows if isinstance(item, Mapping)} \
        if isinstance(partial_rows, list) else {}
    failed: list[dict[str, Any]] = []

    # Recover the actual raw lambda row from the issued normalized partial
    # row: partial[(k,j+1)] = 1/(j+1) * lambda[k,j].
    raw_lambda: dict[tuple[int, int], int] = {}
    for k, j in expected:
        partial_id = indexed.get((k, j + 1))
        if j == 0:
            raw_id = partial_id
        else:
            arguments = row(partial_id).get("arguments", [])
            factors = [value for value in arguments if exact(value, 1, j + 1)]
            candidates = [value for value in arguments if value not in factors]
            raw_id = candidates[0] if len(factors) == 1 and len(candidates) == 1 else None
        if not valid_id(raw_id):
            failed.append({"stage": "raw_lambda", "Z": k, "angle": j,
                           "partial": partial_id, "raw": raw_id})
        else:
            raw_lambda[(k, j)] = raw_id

    def normalized_K(value: Any, order: int) -> bool:
        if order < 2:
            return positive_quotient(value, "2*pi")
        arguments = row(value).get("arguments", [])
        return (row(value).get("operation") == "product" and len(arguments) == 2
                and exact(arguments[0], 1, _factorial(order))
                and positive_quotient(arguments[1], "2*pi"))

    K_rows: dict[int, int] = {}
    one_plus_t2: dict[tuple[int, int], int] = {}
    lambda_rows: list[dict[str, Any]] = []
    for (k, j), raw_id in raw_lambda.items():
        issued_terms = terms(raw_id)
        checks: dict[str, Any] = {"term_count": len(issued_terms) == k + 1}
        for axial_order, term_id in enumerate(issued_terms):
            arguments = row(term_id).get("arguments", [])
            term_ok = (product(term_id) and len(arguments) == 2
                       and normalized_K(arguments[0], axial_order))
            checks[f"term_{axial_order}"] = term_ok
            if term_ok:
                if axial_order in K_rows and K_rows[axial_order] != arguments[0]:
                    checks[f"K_{axial_order}_shared"] = False
                K_rows[axial_order] = arguments[0]
                coefficient_key = (k - axial_order, j)
                if coefficient_key in one_plus_t2 and one_plus_t2[coefficient_key] != arguments[1]:
                    checks[f"one_plus_t2_{coefficient_key}"] = False
                one_plus_t2[coefficient_key] = arguments[1]
            if not term_ok:
                failed.append({"stage": "K_times_one_plus_t2", "Z": k, "angle": j,
                               "term": term_id, "arguments": arguments})
        lambda_rows.append({"Z": k, "angle": j, "node": raw_id,
                            "terms": issued_terms, "checks": checks,
                            "passed": all(checks.values())})

    t2_rows: dict[tuple[int, int], int] = {}
    for key, plus_id in one_plus_t2.items():
        if key == (0, 0):
            arguments = row(plus_id).get("arguments", [])
            constants = [value for value in arguments if exact(value, 1)]
            remainder = [value for value in arguments if value not in constants]
            passed = (row(plus_id).get("operation") == "sum" and len(constants) == 1
                      and len(remainder) == 1)
            if passed:
                t2_rows[key] = remainder[0]
            else:
                failed.append({"stage": "one_plus_t2_constant", "node": plus_id,
                               "arguments": arguments})
        else:
            t2_rows[key] = plus_id

    # Recover t coefficients from the two terms containing t[0,0], then
    # verify every emitted Cauchy product in the exact source order.
    t_rows: dict[tuple[int, int], int] = {}
    t2_checks: list[dict[str, Any]] = []
    for total in range(5):
        for k in range(total + 1):
            j = total - k
            key = (k, j)
            issued_terms = terms(t2_rows.get(key))
            checks: dict[str, Any] = {"node_present": key in t2_rows}
            if key == (0, 0):
                arguments = row(issued_terms[0]).get("arguments", []) if len(issued_terms) == 1 else []
                checks["constant_square"] = (len(arguments) == 2 and arguments[0] == arguments[1]
                                               and product(issued_terms[0]))
                if checks["constant_square"]:
                    t_rows[key] = arguments[0]
            else:
                candidates: list[int] = []
                for term_id in issued_terms:
                    arguments = row(term_id).get("arguments", [])
                    if (len(arguments) == 2 and (arguments[0] == t_rows.get((0, 0))
                                                  or arguments[1] == t_rows.get((0, 0)))):
                        candidates.append(arguments[1] if arguments[0] == t_rows.get((0, 0))
                                          else arguments[0])
                checks["recoverable_coefficient"] = len(set(candidates)) == 1
                if checks["recoverable_coefficient"]:
                    t_rows[key] = candidates[0]
            expected_terms = [(t_rows.get((a, b)), t_rows.get((k - a, j - b)))
                              for a in range(k + 1) for b in range(j + 1)]
            actual_terms = [tuple(row(term_id).get("arguments", []))
                            if product(term_id) else None for term_id in issued_terms]
            checks["exact_t_square_convolution"] = (actual_terms == expected_terms
                                                      and all(value is not None for pair in expected_terms
                                                              for value in pair))
            if not all(checks.values()):
                failed.append({"stage": "t_square_convolution", "Z": k, "angle": j,
                               "node": t2_rows.get(key), "checks": checks,
                               "actual": actual_terms, "expected": expected_terms})
            t2_checks.append({"Z": k, "angle": j, "node": t2_rows.get(key),
                              "checks": checks, "passed": all(checks.values())})

    alpha_rows: dict[int, int] = {}
    ratio_rows: dict[tuple[int, int], int] = {}
    t_checks: list[dict[str, Any]] = []
    for (k, j), t_id in t_rows.items():
        issued_terms = terms(t_id)
        checks: dict[str, Any] = {"term_count": len(issued_terms) == k + 1}
        for axial_order, term_id in enumerate(issued_terms):
            arguments = row(term_id).get("arguments", [])
            term_ok = product(term_id) and len(arguments) == 2
            checks[f"term_{axial_order}"] = term_ok
            if term_ok:
                if axial_order in alpha_rows and alpha_rows[axial_order] != arguments[0]:
                    checks[f"alpha_{axial_order}_shared"] = False
                alpha_rows[axial_order] = arguments[0]
                key = (k - axial_order, j)
                if key in ratio_rows and ratio_rows[key] != arguments[1]:
                    checks[f"ratio_{key}_shared"] = False
                ratio_rows[key] = arguments[1]
            if not term_ok:
                failed.append({"stage": "alpha_times_ratio", "Z": k, "angle": j,
                               "term": term_id, "arguments": arguments})
        t_checks.append({"Z": k, "angle": j, "node": t_id, "checks": checks,
                         "passed": all(checks.values())})

    # Bind alpha to the actual q and positive sqrt(1+u^2) jet rows.
    q_rows = loop.get("q")
    hinv_rows: dict[int, int] = {}
    alpha_source_rows: list[dict[str, Any]] = []
    for order, alpha_id in sorted(alpha_rows.items()):
        normalized = row(alpha_id)
        if order < 2:
            raw_alpha = alpha_id
            normalization_ok = True
        else:
            arguments = normalized.get("arguments", [])
            factors = [value for value in arguments if exact(value, 1, _factorial(order))]
            raw_candidates = [value for value in arguments if value not in factors]
            raw_alpha = raw_candidates[0] if len(factors) == 1 and len(raw_candidates) == 1 else None
            normalization_ok = (normalized.get("operation") == "product" and len(factors) == 1
                                and len(raw_candidates) == 1)
        raw_arguments = row(raw_alpha).get("arguments", [])
        scaled = [value for value in raw_arguments if exact(value, 2)]
        convolution = [value for value in raw_arguments if value not in scaled]
        checks: dict[str, Any] = {
            "normalization": normalization_ok,
            "source_scale_two": row(raw_alpha).get("operation") == "product" and len(scaled) == 1
                               and len(convolution) == 1,
            "q_rows_present": isinstance(q_rows, list) and len(q_rows) >= order + 1,
        }
        if checks["source_scale_two"] and checks["q_rows_present"]:
            convolution_terms = terms(convolution[0])
            checks["q_hinv_convolution_count"] = len(convolution_terms) == order + 1
            for q_order, term_id in enumerate(convolution_terms):
                arguments = row(term_id).get("arguments", [])
                binomial = math.comb(order, q_order)
                if binomial == 1:
                    q_factor, hinv_factor = (arguments + [None, None])[:2]
                    term_ok = (product(term_id) and len(arguments) == 2
                               and q_factor == q_rows[q_order])
                else:
                    q_factor = arguments[1] if len(arguments) == 3 else None
                    hinv_factor = arguments[2] if len(arguments) == 3 else None
                    term_ok = (product(term_id) and len(arguments) == 3
                               and exact(arguments[0], binomial)
                               and q_factor == q_rows[q_order])
                checks[f"q_term_{q_order}"] = term_ok
                if term_ok:
                    hinv_order = order - q_order
                    if hinv_order in hinv_rows and hinv_rows[hinv_order] != hinv_factor:
                        checks[f"hinv_{hinv_order}_shared"] = False
                    hinv_rows[hinv_order] = hinv_factor
        if not all(checks.values()):
            failed.append({"stage": "alpha_source_q_hinv", "order": order,
                           "node": alpha_id, "checks": checks})
        alpha_source_rows.append({"order": order, "node": alpha_id, "raw_node": raw_alpha,
                                  "checks": checks, "passed": all(checks.values())})
    for order, hinv_id in hinv_rows.items():
        if not positive_quotient(hinv_id, "positive sqrt(1+u^2)"):
            failed.append({"stage": "alpha_hinv_source", "order": order, "node": hinv_id})

    # Each ratio row is an actual positive quotient against the same D.  The
    # recurrence's denominator rows are then recovered from its issued rest
    # products and checked against the actual 1 - 2*r*cos + r^2 topology.
    denominator_rows: dict[tuple[int, int], int] = {}
    if (0, 0) in ratio_rows:
        denominator_rows[(0, 0)] = row(ratio_rows[(0, 0)]).get("denominator")
    ratio_checks: list[dict[str, Any]] = []
    for key, ratio_id in ratio_rows.items():
        ratio_row = row(ratio_id)
        checks = {
            "positive_quotient": positive_quotient(ratio_id, "1-2*r*cos"),
            "shared_denominator": valid_id(denominator_rows.get((0, 0)))
                                  and ratio_row.get("denominator") == denominator_rows.get((0, 0)),
        }
        if key != (0, 0):
            numerator_id = ratio_row.get("numerator")
            numerator_row = row(numerator_id)
            rest_id = None
            if numerator_row.get("operation") == "negative":
                rest_id = numerator_row.get("argument")
            else:
                for candidate in numerator_row.get("arguments", []):
                    if row(candidate).get("operation") != "negative":
                        continue
                    candidate_rest = row(candidate).get("argument")
                    if any(tuple(row(term).get("arguments", []))[1:2] and
                           tuple(row(term).get("arguments", []))[1] in ratio_rows.values()
                           for term in terms(candidate_rest) if product(term)):
                        rest_id = candidate_rest
                        break
            if rest_id is None:
                checks["source_denominator_recurrence"] = False
            else:
                checks["source_denominator_recurrence"] = True
                for term in terms(rest_id):
                    arguments = row(term).get("arguments", [])
                    previous = next((previous_key for previous_key, value in ratio_rows.items()
                                     if len(arguments) == 2 and value == arguments[1]), None)
                    if previous is None:
                        checks["source_denominator_recurrence"] = False
                        continue
                    denominator_key = (key[0] - previous[0], key[1] - previous[1])
                    if denominator_key in denominator_rows and denominator_rows[denominator_key] != arguments[0]:
                        checks["source_denominator_recurrence"] = False
                    denominator_rows[denominator_key] = arguments[0]
        if not all(checks.values()):
            failed.append({"stage": "ratio_source", "Z": key[0], "angle": key[1],
                           "node": ratio_id, "checks": checks})
        ratio_checks.append({"Z": key[0], "angle": key[1], "node": ratio_id,
                             "checks": checks, "passed": all(checks.values())})

    # Check the explicit graph terms of D = 1 - 2*r*cos + r*r.  The rows
    # are normalized MixedTaylor rows; the denominator recurrence above binds
    # every issued denominator handle through total degree four.
    denominator_source_checks: dict[str, Any] = {}
    r_rows: dict[int, int] = {}
    cos_rows: dict[int, int] = {}
    denominator_id = denominator_rows.get((0, 0))
    denominator_zero_arguments = row(denominator_id).get("arguments", [])
    scaled_zero = [value for value in denominator_zero_arguments
                   if product(value) and len(row(value).get("arguments", [])) == 2
                   and exact(row(value).get("arguments", [None])[0], -2)]
    square_zero = [value for value in denominator_zero_arguments
                   if product(value) and len(row(value).get("arguments", [])) == 2
                   and row(value).get("arguments", [None, None])[0] == row(value).get("arguments", [None, None])[1]]
    constants_zero = [value for value in denominator_zero_arguments if exact(value, 1)]
    denominator_source_checks["constant_one"] = len(constants_zero) == 1
    denominator_source_checks["negative_two_rc"] = len(scaled_zero) == 1
    denominator_source_checks["r_square"] = len(square_zero) == 1
    if len(scaled_zero) == 1 and len(square_zero) == 1:
        rc_zero = row(scaled_zero[0]).get("arguments", [None, None])[1]
        rc_arguments = row(rc_zero).get("arguments", [])
        square_arguments = row(square_zero[0]).get("arguments", [])
        denominator_source_checks["zero_rc_product"] = len(rc_arguments) == 2
        denominator_source_checks["zero_r_square_product"] = len(square_arguments) == 2 and square_arguments[0] == square_arguments[1]
        if len(rc_arguments) == 2:
            r_rows[0], cos_rows[0] = rc_arguments
        denominator_source_checks["zero_r_shared"] = (len(rc_arguments) == 2
                                                        and square_arguments == [rc_arguments[0], rc_arguments[0]])
    if cos_rows.get(0) is not None:
        cos_zero = row(cos_rows[0])
        denominator_source_checks["cos_source"] = (cos_zero.get("operation") == "analytic_unary"
                                                     and cos_zero.get("name") == "cos"
                                                     and cos_zero.get("argument") == loop.get("angle"))
    else:
        denominator_source_checks["cos_source"] = False
    for j in range(1, 5):
        denominator_id_j = denominator_rows.get((0, j))
        arguments = row(denominator_id_j).get("arguments", [])
        scaled = [value for value in terms(denominator_id_j)
                  if product(value) and len(row(value).get("arguments", [])) == 2
                  and exact(row(value).get("arguments", [None])[0], -2)]
        ok = len(scaled) == 1
        if ok:
            rc = row(scaled[0]).get("arguments", [None, None])[1]
            rc_arguments = row(rc).get("arguments", [])
            ok = len(rc_arguments) == 2 and rc_arguments[0] == r_rows.get(0)
            if ok:
                cos_rows[j] = rc_arguments[1]
        denominator_source_checks[f"cos_row_{j}"] = ok
    for k in range(1, 5):
        denominator_id_k = denominator_rows.get((k, 0))
        scaled = [value for value in terms(denominator_id_k)
                  if product(value) and len(row(value).get("arguments", [])) == 2
                  and exact(row(value).get("arguments", [None])[0], -2)]
        ok = len(scaled) == 1
        if ok:
            rc = row(scaled[0]).get("arguments", [None, None])[1]
            rc_arguments = row(rc).get("arguments", [])
            ok = len(rc_arguments) == 2 and rc_arguments[1] == cos_rows.get(0)
            if ok:
                r_rows[k] = rc_arguments[0]
        denominator_source_checks[f"r_row_{k}"] = ok
    for (k, j), denominator_id in denominator_rows.items():
        if (k, j) == (0, 0):
            continue
        scaled = [value for value in terms(denominator_id)
                  if product(value) and len(row(value).get("arguments", [])) == 2
                  and exact(row(value).get("arguments", [None])[0], -2)]
        ok = len(scaled) == 1
        if ok:
            rc = row(scaled[0]).get("arguments", [None, None])[1]
            rc_arguments = row(rc).get("arguments", [])
            ok = (len(rc_arguments) == 2 and rc_arguments[0] == r_rows.get(k)
                  and rc_arguments[1] == cos_rows.get(j))
        denominator_source_checks[f"rc_row_{k}_{j}"] = ok
    for k in range(1, 5):
        denominator_id_k = denominator_rows.get((k, 0))
        non_scaled = [value for value in terms(denominator_id_k)
                      if not (product(value) and len(row(value).get("arguments", [])) == 2
                              and exact(row(value).get("arguments", [None])[0], -2))]
        rr_terms = terms(non_scaled[0]) if len(non_scaled) == 1 else []
        expected_rr = [(r_rows.get(a), r_rows.get(k - a)) for a in range(k + 1)]
        actual_rr = [tuple(row(value).get("arguments", [])) if product(value) else None
                     for value in rr_terms]
        denominator_source_checks[f"r_square_row_{k}"] = actual_rr == expected_rr

    t0_rows = roots.get("t0") if isinstance(roots, Mapping) else None
    t0_zero = (isinstance(t0_rows, list) and len(t0_rows) >= 5
               and all(value == 0 for value in t0_rows[:5]))
    denominator_source_checks["actual_t0_zero_rows"] = t0_zero
    if not t0_zero:
        failed.append({"stage": "t0_source_row", "rows": t0_rows})

    passed = not failed and all(denominator_source_checks.values())

    def serial_keys(value: Mapping[tuple[int, int], Any]) -> dict[str, Any]:
        return {f"{key[0]},{key[1]}": item for key, item in value.items()}

    return {
        "passed": passed,
        "expected_rows": expected,
        "raw_lambda_rows": lambda_rows,
        "t_square_rows": t2_checks,
        "t_source_rows": t_checks,
        "alpha_source_rows": alpha_source_rows,
        "ratio_source_rows": ratio_checks,
        "K_rows": K_rows,
        "one_plus_t2_rows": serial_keys(one_plus_t2),
        "t_rows": serial_keys(t_rows),
        "alpha_rows": alpha_rows,
        "ratio_rows": serial_keys(ratio_rows),
        "denominator_rows": serial_keys(denominator_rows),
        "r_rows": r_rows,
        "cos_rows": cos_rows,
        "t0_rows": t0_rows,
        "denominator_source_checks": denominator_source_checks,
        "failed": failed,
        "source_recurrence": {
            "lambda": "MixedTaylor(K)*(1+t^2)",
            "t": "MixedTaylor(alpha)*MixedTaylor((cos-r)/(1-2*r*cos+r^2))",
            "alpha": "2*q*hinv with hinv=1/sqrt(1+u^2)",
            "ratio": "(cos-r)/(1-2*r*cos+r^2)",
            "index_domain": "k+j<=4 for lambda source rows; partial rows k+j<=5",
        },
        "source_executed": False,
    }


def _phase_partial_graph_receipt(nodes: Sequence[Any], loop: Mapping[str, Any],
                                 roots: Mapping[str, Any]) -> dict[str, Any]:
    """Bind all 21 normalized mixed phase coefficients to the C5 DAG.

    The C5 source constructs ``partials[(k, 0)]`` from the Phi jet with the
    factorial coefficient and ``partials[(k, j)]`` from the lambda jet with
    coefficient ``1/j``.  This receipt checks those issued graph edges and
    the positive-K descendants; a row count or a worker boolean alone cannot
    discharge the actual implicit-function source graph.
    """

    raw_rows = loop.get("phase_partial_coefficients")
    phi_rows = loop.get("Phi")
    if not isinstance(raw_rows, list) or len(raw_rows) != 21:
        return {"passed": False, "reason": "C5 phase partial row table is not 21 rows", "rows": []}
    if not isinstance(phi_rows, list) or len(phi_rows) != 6:
        return {"passed": False, "reason": "C5 Phi jet row table is not six rows", "rows": []}
    indexed: dict[tuple[int, int], Mapping[str, Any]] = {}
    duplicate = False
    for item in raw_rows:
        if not isinstance(item, Mapping):
            duplicate = True
            continue
        key = (item.get("Z"), item.get("angle"))
        if key in indexed:
            duplicate = True
        indexed[key] = item
    expected_pairs = [(k, j) for k in range(6) for j in range(6 - k)]
    rows: list[dict[str, Any]] = []
    all_passed = not duplicate and set(indexed) == set(expected_pairs)
    for k, j in expected_pairs:
        item = indexed.get((k, j), {})
        node_id = item.get("node") if isinstance(item, Mapping) else None
        row = nodes[node_id] if isinstance(node_id, int) and 0 <= node_id < len(nodes) else {}
        checks: dict[str, Any] = {
            "coordinate_pair": item.get("Z") == k and item.get("angle") == j,
            "node_is_actual_graph_row": isinstance(node_id, int) and 0 <= node_id < len(nodes),
        }
        source_kind = "Phi" if j == 0 else "lambda"
        factor_id = None
        raw_lambda_id = None
        if j == 0:
            phi_id = phi_rows[k]
            if k < 2:
                checks["factorial_normalized_row"] = node_id == phi_id
            else:
                args = row.get("arguments", []) if isinstance(row, Mapping) else []
                factor_id = next((value for value in args
                                  if _exact_rational_row(nodes, value, 1, _factorial(k))), None)
                checks["factorial_normalized_row"] = (row.get("operation") == "product"
                                                       and factor_id is not None
                                                       and phi_id in args)
        else:
            args = row.get("arguments", []) if isinstance(row, Mapping) else []
            factor_id = next((value for value in args
                              if _exact_rational_row(nodes, value, 1, j)), None)
            if j == 1:
                # Multiplication by one is simplified by the source graph;
                # the row itself is the raw lambda coefficient.
                raw_lambda_id = node_id
                checks["one_normalized_lambda_row"] = factor_id is None
            else:
                raw_candidates = [value for value in args if value != factor_id]
                raw_lambda_id = raw_candidates[0] if len(raw_candidates) == 1 else None
                checks["one_over_j_normalized_lambda_row"] = (row.get("operation") == "product"
                                                               and factor_id is not None
                                                               and raw_lambda_id is not None)
            checks.setdefault("one_normalized_lambda_row", True)
        if isinstance(raw_lambda_id, int) and 0 <= raw_lambda_id < len(nodes):
            descendants = _descendants(nodes, raw_lambda_id)
            checks["lambda_K_descendant"] = any(
                isinstance(nodes[value], Mapping)
                and nodes[value].get("operation") == "positive_quotient"
                and nodes[value].get("numerator") == 1
                and "2*pi" in str(nodes[value].get("source_positive_certificate", ""))
                for value in descendants
            )
        elif j > 0:
            checks["lambda_K_descendant"] = False
        checks["source_generator_kind"] = source_kind in ("Phi", "lambda")
        passed = all(checks.values())
        all_passed = all_passed and passed
        rows.append({"Z": k, "angle": j, "node": node_id,
                     "source_generator": source_kind, "factor": factor_id,
                     "raw_lambda": raw_lambda_id, "checks": checks,
                     "passed": passed})
    mixed_lambda = _mixed_lambda_source_receipt(nodes, loop, roots)
    return {
        "passed": bool(all_passed and mixed_lambda["passed"]),
        "row_count": len(raw_rows),
        "expected_pairs": expected_pairs,
        "rows": rows,
        "mixed_lambda_source": mixed_lambda,
        "source_recurrence": {
            "angle_zero": "partials[(k,0)]=1/k!*Phi[k]",
            "angle_positive": "partials[(k,j)]=1/j*lambda[k,j-1]",
            "index_domain": "k+j<=5",
        },
        "source_executed": False,
    }


def _factorial(value: int) -> int:
    result = 1
    for item in range(2, value + 1):
        result *= item
    return result


def _phase_formula_receipt(archives: Mapping[str, Any]) -> dict[str, Any]:
    """Bind inverse, lambda, A_phi and B_phi to the actual C5 graph rows."""

    c5 = archives["c5"]["prescribed_loop_C5_bindings"]
    nodes = archives["c5"]["exact_graph_nodes"]
    charts: dict[str, Any] = {}
    all_passed = True
    for chart in CHARTS:
        data = c5["charts"][chart]
        loop = data["loop"]
        roots = data["roots"]
        phase_partial_graph = _phase_partial_graph_receipt(nodes, loop, roots)
        inverse_id = int(loop["inverse"][0])
        inverse = nodes[inverse_id]
        lambda_id = int(loop["phase_lambda_at_inverse"][0])
        lambda_row = nodes[lambda_id]
        direction_id = int(loop["direction_at_inverse"][0])
        direction = nodes[direction_id]
        inverse_ok = (inverse.get("operation") == "original_monotone_phase_inverse"
                      and inverse.get("phase_function") == loop["Phi"][0]
                      and inverse.get("target_phase") == data["actual_phase"]
                      and inverse.get("phase_endpoints") == [0, 1]
                      and inverse.get("branch_inverse_is_original_not_support_cap") is True
                      and inverse.get("true_prescribed_source_C5_implicit_inverse") is True)
        direction_ok = (direction.get("operation") == "substitute_original_inverse_angle"
                        and direction.get("inverse_angle") == inverse_id
                        and direction.get("angle_variable") == inverse.get("angle_variable"))
        lambda_body_id = lambda_row.get("body") if isinstance(lambda_row, Mapping) else None
        lambda_body = nodes[lambda_body_id] if isinstance(lambda_body_id, int) else {}
        lambda_args = list(lambda_body.get("arguments", [])) if isinstance(lambda_body, Mapping) else []
        inv_k_ids = [value for value in lambda_args if isinstance(value, int)
                     and isinstance(nodes[value], Mapping)
                     and nodes[value].get("operation") == "positive_quotient"
                     and nodes[value].get("numerator") == 1
                     and "2*pi" in str(nodes[value].get("source_positive_certificate", ""))]
        one_plus_ids = [value for value in lambda_args if isinstance(value, int)
                        and isinstance(nodes[value], Mapping)
                        and nodes[value].get("operation") == "sum"]
        direction_body = direction.get("body") if isinstance(direction, Mapping) else None
        t_squared_ok = False
        if len(one_plus_ids) == 1 and isinstance(direction_body, int):
            plus_args = list(nodes[one_plus_ids[0]].get("arguments", []))
            t_squared_ok = any(_exact_rational_row(nodes, value, 1) for value in plus_args) and any(
                isinstance(value, int) and isinstance(nodes[value], Mapping)
                and nodes[value].get("operation") == "product"
                and nodes[value].get("arguments") == [direction_body, direction_body]
                for value in plus_args)
        lambda_ok = (lambda_row.get("operation") == "substitute_original_inverse_angle"
                     and lambda_row.get("inverse_angle") == inverse_id
                     and lambda_row.get("angle_variable") == inverse.get("angle_variable")
                     and lambda_body.get("operation") == "product"
                     and len(inv_k_ids) == 1 and t_squared_ok)
        inv_lambda_id = None
        for root_id in (loop["A_phi"][0], loop["B_over_Pstar_phi"][0]):
            active_id = nodes[root_id].get("active_body") if isinstance(nodes[root_id], Mapping) else None
            for candidate in _descendants(nodes, active_id):
                row = nodes[candidate]
                if row.get("operation") == "positive_quotient" and row.get("numerator") == 1 and row.get("denominator") == lambda_id:
                    inv_lambda_id = candidate
            if inv_lambda_id is not None:
                break
        a_root = int(roots["a"][0])
        e_root = int(roots["E"][0])
        t_root = direction_body if isinstance(direction_body, int) else -1
        a_formula = _formula_core(nodes, int(loop["A_phi"][0]), a_root,
                                   inv_lambda_id if inv_lambda_id is not None else -1,
                                   inv_k_ids[0] if inv_k_ids else -1, "A_phi")
        b_formula = _formula_core(nodes, int(loop["B_over_Pstar_phi"][0]), e_root,
                                  inv_lambda_id if inv_lambda_id is not None else -1,
                                  inv_k_ids[0] if inv_k_ids else -1, "B_phi")
        lazy_rows_ok = all(isinstance(nodes[int(node_id)], Mapping)
                           and nodes[int(node_id)].get("operation") == "original_lazy_flat_branch"
                           and nodes[int(node_id)].get("Delta") == data["roots"]["Delta"][0]
                           and nodes[int(node_id)].get("eta") == 21
                           and nodes[int(node_id)].get("flat_value") == 0
                           and nodes[int(node_id)].get("active_body_not_evaluated_on_flat_branch") is True
                           for key in ("A_phi", "B_over_Pstar_phi") for node_id in loop[key])
        chart_passed = bool(inverse_ok and direction_ok and lambda_ok and a_formula["passed"]
                            and b_formula["passed"] and lazy_rows_ok
                            and phase_partial_graph["passed"])
        all_passed = all_passed and chart_passed
        charts[chart] = {
            "inverse": {"node": inverse_id, "phase_function": inverse.get("phase_function"),
                        "target_phase": inverse.get("target_phase"), "phase_endpoints": inverse.get("phase_endpoints"),
                        "passed": inverse_ok},
            "phase_lambda": {"node": lambda_id, "body": lambda_body_id, "K_factor": inv_k_ids[0] if inv_k_ids else None,
                              "direction_body": direction_body, "lambda_equals_K_times_one_plus_t_squared": lambda_ok,
                              "passed": bool(direction_ok and lambda_ok)},
            "A_phi": a_formula,
            "B_phi": b_formula,
            "lazy_flat_rows": {"passed": lazy_rows_ok, "count": len(loop["A_phi"]) + len(loop["B_over_Pstar_phi"])},
            "phase_partial_graph": phase_partial_graph,
            "passed": chart_passed,
        }
    return {
        "passed": all_passed,
        "charts": charts,
        "formula_identities": {
            "A_phi": "a/2*(1-1/(2*pi*lambda))",
            "B_phi": "E/2*(-a*t(psi)/(2*pi*lambda)-b)",
            "partial_psi_inverse": "D_target psi = 1/lambda after actual original inverse substitution",
        },
        "source_executed": False,
    }


def _parent_graph_wrapper_receipt(parent: Mapping[str, Any], actual: Mapping[str, Any]) -> dict[str, Any]:
    """Verify the saved 185690-node lower-wrapper substitution rows.

    The parent graph patch is consumed as a source receipt only.  No graph
    contraction or node rebuilding occurs here.  The two wrapper rows must be
    the actual finite-boundary lower R_Q_y rows substituted at coordinate zero,
    while the upper trace remains the standalone Q4 endpoint trace.
    """

    pressure = parent.get("actual_long_reference_pressure_source")
    patch = pressure.get("graph_patch") if isinstance(pressure, Mapping) else None
    joined = pressure.get("joined_pressure_source") if isinstance(pressure, Mapping) else None
    endpoint_receipt = pressure.get("endpoint_identity_receipts") if isinstance(pressure, Mapping) else None
    if not isinstance(patch, Mapping) or not isinstance(joined, Mapping):
        return {"passed": False, "reason": "accepted parent graph patch or joined source is missing", "rows": {}}
    expected_canonical_digest = "8422242e6c2bce6750323cf2aafd2a5a4bce2bf406f033ddd6a18ad61ae33d28"
    base = patch.get("base_node_count")
    suffix = patch.get("appended_nodes")
    metadata_ok = (patch.get("base_archive") == Q4_CACHE
                   and base == 185308
                   and isinstance(suffix, list) and len(suffix) == 382
                   and patch.get("base_prefix_unchanged") is True
                   and patch.get("new_node_start") == 185308
                   and patch.get("appended_node_count") == 382
                   and patch.get("full_node_count") == 185690
                   and patch.get("full_canonical_nodes_sha256") == expected_canonical_digest
                   and patch.get("joined_contracted") is None
                   and pressure.get("source_execution") is False)
    endpoint_handles = patch.get("endpoint_handles", {})
    expected_handles = {
        "post_power_upper_standalone": list(SEAM_DATA["post_power_upper_to_long_lower"]["upper_endpoint_R_Q_y"][:2]),
        "long_reshape_lower_q0": list(SEAM_DATA["post_power_upper_to_long_lower"]["graph_lower_wrappers"]),
        "long_reshape_upper_standalone": list(SEAM_DATA["long_upper_to_reference_lower"]["upper_endpoint_R_Q_y"][:2]),
        "reference_lower_q0": list(SEAM_DATA["long_upper_to_reference_lower"]["graph_lower_wrappers"]),
    }
    handle_checks = {name: endpoint_handles.get(name) == value for name, value in expected_handles.items()}
    rows: dict[str, Any] = {}
    all_rows = True
    for name, data in SEAM_DATA.items():
        lower_trace = actual["rqy_endpoint_traces"][name]["lower_chart_R_Q_y"]
        wrappers = list(data["graph_lower_wrappers"])
        upper = list(data["upper_endpoint_R_Q_y"][:2])
        seam_rows: list[dict[str, Any]] = []
        for order, (wrapper_id, lower_id) in enumerate(zip(wrappers, lower_trace[:2])):
            row = suffix[wrapper_id - base] if isinstance(wrapper_id, int) and base <= wrapper_id < base + len(suffix) else None
            lower_coordinate = 79966 if data["lower_chart"] == "long_reshape" else 81520
            expected = {"operation": "function_substitution", "expression": lower_id,
                        "variable": lower_coordinate, "value": 0,
                        "Z_independent_substitution": True}
            check = isinstance(row, Mapping) and all(row.get(key) == value for key, value in expected.items())
            seam_rows.append({"order": order, "wrapper": wrapper_id, "lower_R_Q_y": lower_id,
                              "row": dict(row) if isinstance(row, Mapping) else None,
                              "expected": expected, "passed": check})
        joined_upper_key = "post_power_upper_RQy" if name == "post_power_upper_to_long_lower" else "long_reshape_upper_RQy"
        joined_lower_key = "long_reshape_lower_RQy" if name == "post_power_upper_to_long_lower" else "reference_lower_RQy"
        joined_ok = (joined.get(joined_upper_key) == upper
                     and joined.get(joined_lower_key) == wrappers)
        # The parent endpoint receipt intentionally remains open for the full
        # modified seam.  Its presence identifies the accepted graph patch;
        # the actual wrapper rows and Q4 endpoint rows are checked here.
        seam_passed = bool(all(item["passed"] for item in seam_rows) and joined_ok
                           and isinstance(endpoint_receipt, Mapping))
        all_rows = all_rows and seam_passed
        rows[name] = {"wrappers": seam_rows, "joined_upper": joined.get(joined_upper_key),
                      "joined_lower": joined.get(joined_lower_key), "expected_upper": upper,
                      "expected_lower": wrappers, "joined_handles_match": joined_ok,
                      "passed": seam_passed}
    passed = bool(metadata_ok and all(handle_checks.values()) and all_rows)
    return {
        "passed": passed,
        "base_node_count": base,
        "full_node_count": patch.get("full_node_count"),
        "appended_node_count": len(suffix) if isinstance(suffix, list) else None,
        "full_canonical_nodes_sha256": patch.get("full_canonical_nodes_sha256"),
        "canonical_digest_expected": expected_canonical_digest,
        "metadata_ok": metadata_ok,
        "endpoint_handle_checks": handle_checks,
        "rows": rows,
        "source_execution": False,
    }


def _parent_endpoint_radius_receipt(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the saved common-phase/physical-radius endpoint theorem.

    The R_Q_y source rows carry a radius factor.  The parent loop receipt is
    therefore required for the endpoint step that identifies the upper and
    lower physical radii; a Q_y equality alone cannot discharge this gate.
    """

    loop = parent.get("actual_long_reference_loop_seams")
    receipt = loop.get("endpoint_identity_receipts") if isinstance(loop, Mapping) else None
    raw = receipt.get("seams") if isinstance(receipt, Mapping) else None
    pressure = parent.get("actual_long_reference_pressure_source")
    pressure_endpoint = (pressure.get("endpoint_identity_receipts")
                         if isinstance(pressure, Mapping) else None)
    common_N = (pressure_endpoint.get("common_N")
                if isinstance(pressure_endpoint, Mapping) else None)
    expected_names = (
        "post_power_upper_to_long_reshape_lower",
        "long_reshape_upper_to_reference_lower",
    )
    rows: dict[str, Any] = {}
    all_rows = isinstance(raw, list) and len(raw) == 2
    for index, seam in enumerate(SEAM_DATA):
        row = raw[index] if isinstance(raw, list) and index < len(raw) else None
        checks = {
            "name": isinstance(row, Mapping) and row.get("name") == expected_names[index],
            "left_fraction": isinstance(row, Mapping) and row.get("left_fraction") == 1,
            "right_fraction": isinstance(row, Mapping) and row.get("right_fraction") == 0,
            "direct_endpoint_radius_equality": isinstance(row, Mapping) and row.get("direct_endpoint_radius_equality") is True,
            "global_offset_equality": isinstance(row, Mapping) and row.get("global_offset_equality") is True,
            "actual_phase_equality_under_common_N": isinstance(row, Mapping) and row.get("actual_phase_equality_under_common_N") is True,
            "integer_offset_difference_not_assumed": isinstance(row, Mapping) and row.get("integer_offset_difference_assumed") is False,
            "source_radius_identity": isinstance(row, Mapping) and isinstance(row.get("source_radius_identity"), str),
        }
        passed = all(checks.values())
        all_rows = all_rows and passed
        rows[seam] = {"row": dict(row) if isinstance(row, Mapping) else None,
                      "checks": checks, "passed": passed}
    passed = bool(isinstance(receipt, Mapping)
                  and receipt.get("passed") is True
                  and receipt.get("source_executed") is False
                  and common_N == 18
                  and all_rows)
    return {
        "passed": passed,
        "common_N": common_N,
        "rows": rows,
        "source_executed": False,
    }


def _parent_axial_scope_receipt(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Check the inherited certified axial Z window without enlarging it."""

    expected = list(AXIAL_SCOPE["Z_interval"])
    sections = (
        "actual_micro_q_mixed_bounds",
        "actual_micro_phase_inverse_bounds",
        "actual_micro_b_p2_mixed_bounds",
        "actual_micro_mean_Z3_bounds",
        "actual_micro_Q_increment_bounds",
    )
    rows: dict[str, Any] = {}
    all_rows = True
    for name in sections:
        section = parent.get(name)
        domain = section.get("domain") if isinstance(section, Mapping) else None
        value = domain.get("Z") if isinstance(domain, Mapping) else None
        checks = {
            "section_passed": isinstance(section, Mapping) and section.get("passed") is True,
            "source_executed_false": isinstance(section, Mapping)
            and section.get("source_executed", False) is False
            and section.get("accepted_ancestors_replayed", False) is False,
            "Z_interval": value == expected,
            "global_five_moment_closure_false": not (isinstance(section, Mapping)
                                                       and section.get("global_five_moment_closure") is True),
        }
        passed = all(checks.values())
        rows[name] = {"Z": value, "checks": checks, "passed": passed}
        all_rows = all_rows and passed
    background = parent.get("actual_long_reference_background_I1_numeric")
    directions = background.get("directions") if isinstance(background, Mapping) else None
    background_domain = directions.get("Z_domain") if isinstance(directions, Mapping) else None
    background_ok = (isinstance(background, Mapping) and background.get("passed") is True
                     and background.get("source_executed") is False
                     and background_domain == "[37/100,93/250]")
    rows["actual_long_reference_background_I1_numeric"] = {
        "Z_domain": background_domain, "passed": background_ok,
    }
    all_rows = all_rows and background_ok
    return {
        "passed": all_rows,
        "Z_interval": expected,
        "rows": rows,
        "global_five_moment_closure": False,
        "source_executed": False,
    }


def _modified_my_source_receipt(archives: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the physical m_y identity to C5 transport and Q4 source rows."""

    c6 = archives["c6"]["prescribed_R100_Rm_C6_bindings"]
    c5 = archives["c5"]["prescribed_loop_C5_bindings"]
    q4 = archives["q4"]["modified_Q4_pressure_C3_bindings"]
    c5_nodes = archives["c5"]["exact_graph_nodes"]
    q4_nodes = archives["q4"]["exact_graph_nodes"]
    history = parent.get("actual_long_reference_history_seams", {})
    conditional = history.get("conditional_history_ODE", {}) if isinstance(history, Mapping) else {}
    q0 = conditional.get("q0", {}) if isinstance(conditional, Mapping) else {}
    q0_ok = (conditional.get("passed") is True and q0.get("correction") is True
             and q0.get("exp_zero_is_one") is True and q0.get("integral_zero_is_zero") is True)
    rows: dict[str, Any] = {}
    all_rows = True
    for chart in CHARTS:
        c5chart = c5["charts"][chart]
        c6chart = c6["charts"][C6_CHART[chart]]
        q4chart = q4["charts"][chart]
        row_receipts: list[dict[str, Any]] = []
        for order in ORDERS:
            s = c5chart["finite_N"]["delta_V"][order]
            m_increment = c5chart["finite_N"]["five_signed_density_increments"]["m"][order]
            v_n = c5chart["finite_N"]["V_N"][order]
            native_v = c6chart["fields"]["V"][order]
            native_m = c6chart["fields"]["m"][order]
            correction_m = c5chart["history"]["m"][order]
            modified_m = c5chart["modified_fields"]["m"][order]
            q4_radial_m_y = q4chart["radial"]["m_y"][order]
            q4_physical_m_y = q4chart["physical"]["m_y"][order]
            q4_physical_u = q4chart["physical"]["U"][order]
            q4_physical_m = q4chart["physical"]["m"][order]
            v_row = c5_nodes[v_n] if isinstance(v_n, int) else {}
            mod_m_row = c5_nodes[modified_m] if isinstance(modified_m, int) else {}
            radial_row = q4_nodes[q4_radial_m_y] if isinstance(q4_radial_m_y, int) else {}
            physical_row = q4_nodes[q4_physical_m_y] if isinstance(q4_physical_m_y, int) else {}
            radial_neg = q4_nodes[radial_row.get("arguments", [None, None])[1]] if isinstance(radial_row, Mapping) and len(radial_row.get("arguments", [])) == 2 else {}
            physical_neg = q4_nodes[physical_row.get("arguments", [None, None])[1]] if isinstance(physical_row, Mapping) and len(physical_row.get("arguments", [])) == 2 else {}
            check = {
                "signed_density_m_equals_B_over_Pstar_over_N": m_increment == s,
                "V_N_equals_V_plus_s": v_row.get("operation") == "sum" and v_row.get("arguments") == [native_v, s],
                "modified_m_equals_native_m_plus_transport": mod_m_row.get("operation") == "sum" and mod_m_row.get("arguments") == [native_m, correction_m],
                "radial_m_y_equals_V_N_minus_modified_m": radial_row.get("operation") == "sum" and radial_row.get("arguments") == [v_n, radial_row.get("arguments", [None, None])[1]] and radial_neg.get("operation") == "negative" and radial_neg.get("argument") == modified_m,
                "physical_m_y_equals_U_minus_m": physical_row.get("operation") == "sum" and physical_row.get("arguments") == [q4_physical_u, physical_row.get("arguments", [None, None])[1]] and physical_neg.get("operation") == "negative" and physical_neg.get("argument") == q4_physical_m,
            }
            row_passed = all(check.values())
            all_rows = all_rows and row_passed
            row_receipts.append({"order": order, "s": s, "m_increment": m_increment,
                                 "V_N": v_n, "native_V": native_v, "native_m": native_m,
                                 "transport_m": correction_m, "modified_m": modified_m,
                                 "radial_m_y": q4_radial_m_y, "physical_m_y": q4_physical_m_y,
                                 "checks": check, "passed": row_passed})
        rows[chart] = {"rows": row_receipts, "passed": all(item["passed"] for item in row_receipts)}
    try:
        import sympy as sp
        V, s, m, dm, Pstar = sp.symbols("V s m delta_m Pstar")
        native_identity = sp.simplify((V + s - (m + dm)) - ((V - m) + (s - dm)))
        physical_identity = sp.simplify(Pstar * (V + s - (m + dm))
                                        - (Pstar * (V + s) - Pstar * (m + dm)))
        algebra = {"native_residual": str(native_identity), "physical_residual": str(physical_identity),
                   "passed": native_identity == 0 and physical_identity == 0}
    except Exception as exc:  # pragma: no cover - dependency diagnostic
        algebra = {"passed": False, "error": f"SymPy unavailable: {exc}"}
    return {
        "passed": bool(q0_ok and all_rows and algebra.get("passed") is True),
        "q0_lemma": {"passed": q0_ok, "receipt": q0},
        "transport_identity": "(m+delta_m)_y=(V+s)-(m+delta_m), with s=B_over_Pstar/N",
        "physical_identity": "Pstar*((V+s)-(m+delta_m))=U-m",
        "algebra": algebra,
        "charts": rows,
        "source_executed": False,
    }


def _parent_lemma_gate(parent: Mapping[str, Any]) -> dict[str, Any]:
    sections = [parent[name] for name in ("actual_long_reference_native_seams", "actual_long_reference_history_seams", "actual_long_reference_loop_seams")]
    required_false = ("source_executed", "accepted_builders_executed", "accepted_checkers_executed",
                      "accepted_producers_executed", "owners_executed", "norms_executed",
                      "recurrences_executed", "quadrature_executed")
    safe = all(section.get("passed") is True and section.get("source_executed") is False
               and all(section.get(key) is not True for key in required_false if key != "source_executed")
               for section in sections)
    contract = parent.get("contract", {})
    source_bound = isinstance(contract, Mapping) and contract.get("source_bound_seam_prerequisites") is True
    no_old_claim = isinstance(contract, Mapping) and contract.get("modified_RQy_seams_proved") is False
    return {
        "passed": bool(safe and source_bound and no_old_claim),
        "sections_passed": safe,
        "source_bound_contract": source_bound,
        "old_modified_claim_absent": no_old_claim,
        "consumed_as_lemmas_only": True,
        "source_executed": False,
    }


def _json_digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _walk_mappings(value: Any) -> Sequence[Mapping[str, Any]]:
    found: list[Mapping[str, Any]] = []
    seen: set[int] = set()
    def visit(item: Any) -> None:
        if isinstance(item, Mapping):
            identity = id(item)
            if identity in seen:
                return
            seen.add(identity)
            found.append(item)
            for child in item.values():
                visit(child)
        elif isinstance(item, list):
            identity = id(item)
            if identity in seen:
                return
            seen.add(identity)
            for child in item:
                visit(child)
    visit(value)
    return found


def _find_alias(value: Any, aliases: Sequence[str]) -> list[Any]:
    alias_set = set(aliases)
    return [mapping[key] for mapping in _walk_mappings(value) for key in mapping.keys() if key in alias_set]


def _proved_field(value: Any, key: str) -> bool:
    if not isinstance(value, Mapping):
        return False
    item = value.get(key)
    if isinstance(item, Mapping):
        return item.get("proved") is True or item.get("passed") is True
    return item is True


def _packet_provenance(root: Path, packet: Any, label: str, expected_module: str) -> dict[str, Any]:
    if not isinstance(packet, Mapping):
        return {"provided": False, "passed": False, "reason": f"{label} packet not provided"}
    parent_ok = packet.get("parent_report_sha256") in (None, PARENT_SHA256)
    family_ok = packet.get("source_family") in (None, _load_source_family(root))
    source_ok = packet.get("source_executed") is False
    path = packet.get("module_path") or packet.get("module") or expected_module
    if isinstance(path, Mapping):
        path = path.get("path") or path.get("name")
    path = str(path)
    if path.endswith(".py") and not path.startswith("src/"):
        path = expected_module
    own_pin = None
    pins = packet.get("source_pins")
    if isinstance(pins, Mapping):
        own_pin = pins.get(expected_module)
        if own_pin is None:
            own_pin = next((v for k, v in pins.items() if str(k).endswith(Path(expected_module).name)), None)
    module_hash = packet.get("module_sha256") or packet.get("code_sha256") or packet.get("self_sha256")
    module_path = None
    try:
        module_path = _locate(root, expected_module)
    except FileNotFoundError:
        module_path = None
    actual_hash = _sha(module_path) if module_path is not None else None
    # Both upstream workers bind their own source through a self pin rather
    # than duplicating ``module_sha256`` at the packet root.  Use that pin as
    # the code hash when the root field is absent.
    if module_hash is None and isinstance(own_pin, Mapping):
        module_hash = own_pin.get("actual_sha256")
    hash_ok = bool(actual_hash and module_hash == actual_hash)
    if isinstance(own_pin, Mapping):
        hash_ok = hash_ok and own_pin.get("actual_sha256") == actual_hash and own_pin.get("expected_sha256") == actual_hash and own_pin.get("passed") is True and own_pin.get("source_executed") is False
    else:
        hash_ok = False
    parent_pin_ok = True
    if isinstance(pins, Mapping):
        p = pins.get(PARENT) or next((v for k, v in pins.items() if str(k).endswith(Path(PARENT).name)), None)
        parent_pin_ok = isinstance(p, Mapping) and p.get("expected_sha256") == PARENT_SHA256 and p.get("actual_sha256") == PARENT_SHA256 and p.get("passed") is True and p.get("source_executed") is False
        if p is None:
            validation = packet.get("parent_validation")
            parent_pin_ok = isinstance(validation, Mapping) and validation.get("passed") is True and validation.get("actual_sha256") == PARENT_SHA256
    else:
        validation = packet.get("parent_validation")
        parent_pin_ok = isinstance(validation, Mapping) and validation.get("passed") is True and validation.get("actual_sha256") == PARENT_SHA256
    # The worker's own source-family field is not enough to identify the
    # graph that supplied its rows.  Require the six actual C6/C5/Q4 source
    # and cache pins as well.  The derivative worker has no parent pin, so its
    # parent_validation receipt is handled above; both workers still have the
    # same cache-family pins here.
    required_source_pins = (C6_SOURCE, C6_CACHE, C5_SOURCE, C5_CACHE, Q4_SOURCE, Q4_CACHE)
    cache_pin_rows: dict[str, Any] = {}
    cache_pins_ok = isinstance(pins, Mapping)
    for required in required_source_pins:
        row = pins.get(required) if isinstance(pins, Mapping) else None
        if row is None and isinstance(pins, Mapping):
            row = next((value for key, value in pins.items()
                        if str(key).endswith(Path(required).name)), None)
        expected_hash = PINNED_HASHES[required]
        row_ok = (isinstance(row, Mapping)
                  and row.get("expected_sha256") == expected_hash
                  and row.get("actual_sha256") == expected_hash
                  and row.get("passed") is True
                  and row.get("source_executed") is False)
        cache_pins_ok = cache_pins_ok and row_ok
        cache_pin_rows[required] = {
            "present": isinstance(row, Mapping),
            "expected_sha256": expected_hash,
            "actual_sha256": row.get("actual_sha256") if isinstance(row, Mapping) else None,
            "passed": row_ok,
        }
    return {
        "provided": True,
        "passed": bool(parent_ok and family_ok and source_ok and hash_ok and parent_pin_ok and cache_pins_ok),
        "parent_report_sha256_ok": parent_ok,
        "source_family_ok": family_ok,
        "source_executed_false": source_ok,
        "module_path": expected_module,
        "module_sha256": actual_hash,
        "packet_module_sha256": module_hash,
        "module_pin_ok": hash_ok,
        "parent_pin_ok": parent_pin_ok,
        "source_pins_present": isinstance(pins, Mapping),
        "cache_family_pins_ok": cache_pins_ok,
        "cache_family_pins": cache_pin_rows,
    }


def _load_source_family(root: Path) -> Mapping[str, Any]:
    return _json(_locate(root, PARENT))["source_family"]


def _chart_packet(packet: Mapping[str, Any], chart: str) -> Mapping[str, Any] | None:
    aliases = (chart, C6_CHART[chart])
    for mapping in _walk_mappings(packet):
        for key in aliases:
            value = mapping.get(key)
            if isinstance(value, Mapping) and any(name in value for name in (
                    "row_checks", "root_receipts", "seed_map", "seed_bindings", "field_rows",
                    "phase_partial_rows", "firstjet_rows", "fields", "Q4_radial_background")):
                return value
    return None


def _packet_seed_rows(packet: Mapping[str, Any], chart: str, expected: Mapping[str, Any],
                      field_names: Sequence[str] = FIELDS) -> tuple[bool, dict[str, Any]]:
    candidate = _chart_packet(packet, chart)
    if candidate is None:
        for value in _find_alias(packet, ("seed_map", "source_seed_map", "native_seed_map", "seed_identity", "field_rows")):
            if isinstance(value, Mapping) and any(name in value for name in field_names):
                candidate = {"field_rows": value}
                break
    if candidate is None:
        return False, {"passed": False, "reason": f"{chart}: no source seed table"}
    # The derivative worker nests the exact table under ``seed_bindings``;
    # the control worker exposes the same source rows under source_graph as
    # ``fields`` and ``Q4_radial_background``.  Normalize both to the common
    # field_rows representation while retaining every actual handle.
    if isinstance(candidate.get("seed_bindings"), Mapping):
        candidate = candidate["seed_bindings"]
    if isinstance(candidate.get("fields"), Mapping) and isinstance(candidate.get("Q4_radial_background"), Mapping):
        candidate = {
            "field_rows": {
                field: [{"background_field_handle": source, "radial_background_handle": physical,
                          "ordinary_order": order, "J_handle": expected[field][order]["physical_J_node"]}
                         for order, (source, physical) in enumerate(zip(candidate["fields"][field], candidate["Q4_radial_background"][field]))]
                for field in field_names
            }
        }
    seed = candidate.get("seed_map", candidate)
    if isinstance(seed, Mapping) and isinstance(seed.get("field_rows"), Mapping):
        seed = seed["field_rows"]
    if not isinstance(seed, Mapping):
        return False, {"passed": False, "reason": f"{chart}: malformed seed table"}
    rows_receipt: dict[str, Any] = {}
    passed = True
    for field in field_names:
        raw = seed.get(field)
        if not isinstance(raw, list):
            passed = False
            rows_receipt[field] = {"passed": False, "reason": "missing rows"}
            continue
        normalized: list[dict[str, Any]] = []
        for order, row in enumerate(raw):
            if not isinstance(row, Mapping):
                passed = False
                continue
            source = row.get("source_field_node", row.get("source_node", row.get("field_root", row.get("background_field_handle"))))
            physical = row.get("physical_radial_row_node", row.get("actual_root_node", row.get("radial_root", row.get("radial_background_handle"))))
            ok = source == expected[field][order]["source_field_node"] if order < len(expected[field]) else False
            ok = ok and physical == expected[field][order]["physical_radial_row_node"]
            if "ordinary_order" in row:
                ok = ok and row.get("ordinary_order") == order
            if "physical_J_node" in row:
                ok = ok and row.get("physical_J_node") == expected[field][order]["physical_J_node"]
            if "J_handle" in row:
                ok = ok and row.get("J_handle") == expected[field][order]["physical_J_node"]
            passed = passed and ok
            normalized.append({"order": order, "source": source, "physical": physical, "passed": ok})
        if len(raw) != len(expected[field]):
            passed = False
        rows_receipt[field] = normalized
    return passed, {"passed": passed, "rows": rows_receipt}


def _packet_row(packet: Mapping[str, Any], row: Mapping[str, Any], expected_source: int, expected_actual: int, order: int, field: str) -> bool:
    if not isinstance(row, Mapping):
        return False
    source = row.get("source_node", row.get("source_root", row.get("native_root", row.get("source_handle"))))
    actual = row.get("actual_root_node", row.get("actual_node", row.get("physical_root", row.get("root_handle"))))
    if source != expected_source or actual != expected_actual:
        return False
    row_field = row.get("field", row.get("field_name", field))
    row_order = row.get("order", row.get("ordinary_order", order))
    if row_field != field or row_order != order:
        return False
    return row.get("passed") is True


def _derivative_packet_gate(packet: Any, expected: Mapping[str, Any], provenance: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(packet, Mapping):
        return {"passed": False, "provenance": provenance, "charts": {}, "normalized_roots": 0, "expected_roots": 36,
                "seed_tables": {}, "reason": "derivative packet not provided"}
    charts_receipt: dict[str, Any] = {}
    seed_receipt: dict[str, Any] = {}
    rows_total = 0
    rows_passed = 0
    seeds_passed = True
    phases_passed = True
    for chart in CHARTS:
        chart_expected = expected["charts"][chart]
        seed_ok, seed_detail = _packet_seed_rows(packet, chart, expected["seed_rows"][chart], FIELDS[:-1])
        seed_receipt[chart] = seed_detail
        seeds_passed = seeds_passed and seed_ok
        candidate = _chart_packet(packet, chart)
        row_checks = None
        if isinstance(candidate, Mapping):
            row_checks = candidate.get("row_checks", candidate.get("root_receipts"))
        if not isinstance(row_checks, list):
            row_checks = []
            for value in _find_alias(packet, ("row_checks", "normalized_rows", "derivative_rows")):
                if isinstance(value, list) and len(value) >= 12:
                    row_checks = value
                    break
        expected_rows = []
        for field in ("A", "B_over_Pstar"):
            source = expected["charts"][chart]["loop_source_roots"][field]
            actual = expected["charts"][chart]["slow_derivative_roots"][field]
            expected_rows.extend((field, order, int(s), int(a)) for order, (s, a) in enumerate(zip(source, actual)))
        normalized: list[dict[str, Any]] = []
        for index, (field, order, source, actual) in enumerate(expected_rows):
            row = row_checks[index] if index < len(row_checks) else None
            ok = _packet_row(packet, row, source, actual, order, field)
            normalized.append({"field": field, "order": order, "source_node": source, "actual_root_node": actual, "passed": ok})
            rows_total += 1
            rows_passed += int(ok)
        # The worker must issue all six phase-partial rows at each ordinary
        # order; a scalar phase-partial boolean is insufficient.
        phase_rows = None
        if isinstance(candidate, Mapping):
            phase_rows = candidate.get("phase_partial_rows", candidate.get("actual_phase_partials"))
        if not isinstance(phase_rows, list):
            phase_rows = next((value for value in _find_alias(packet, ("phase_partial_rows", "actual_phase_partials", "phase_partial_coefficients"))
                               if isinstance(value, list) and len(value) >= 21), None)
        phase_rows_from_packet = isinstance(phase_rows, list) and len(phase_rows) == 21
        phase_ok = phase_rows_from_packet
        if not phase_ok:
            # The derivative worker intentionally emits the actual 36 A/B
            # roots.  When it does not duplicate the 21 mixed coefficients,
            # use the independently checked C5 source-jet graph receipt,
            # never a row-count fallback.
            phase_ok = expected["charts"][chart].get("phase_partial_graph_bound") is True
        if phase_rows_from_packet:
            expected_phase = expected["charts"][chart]["phase_partial_coefficients"]
            # Require each issued (Z, angle) node if the worker exposed rows;
            # this binds the actual source partial graph rather than its count.
            for expected_row in expected_phase:
                if not any(isinstance(row, Mapping) and row.get("Z") == expected_row.get("Z")
                           and row.get("angle") == expected_row.get("angle")
                           and row.get("node") == expected_row.get("node")
                           for row in phase_rows):
                    phase_ok = False
                    break
        phases_passed = phases_passed and phase_ok
        charts_receipt[chart] = {"row_checks": normalized,
                                 "phase_partial_rows_present": phase_rows_from_packet,
                                 "phase_partial_graph_bound": phase_ok,
                                 "passed": bool(seed_ok and all(row["passed"] for row in normalized) and phase_ok)}
    dag_gates = packet.get("gates")
    semantic = (packet.get("conditional_passed") is True
                and isinstance(dag_gates, Mapping)
                and all(value is True for value in dag_gates.values()))
    # The derivative worker deliberately leaves this external seed theorem
    # false.  Preserve that fact in the receipt without using it as a
    # substitute for the actual 36 held-phase roots above.
    seed_semantic = packet.get("seed_identity_discharged") is True
    root_wrappers = _wrapper_mapping(packet)
    wrappers_ok = root_wrappers["passed"]
    return {
        "passed": bool(provenance.get("passed") and semantic and seeds_passed
                        and phases_passed and rows_total == 36 and rows_passed == 36),
        "provenance": _copy(provenance),
        "charts": charts_receipt,
        "seed_tables": seed_receipt,
        "normalized_roots": rows_passed,
        "expected_roots": 36,
        "conditional_native_derivative_graph_passed": semantic,
        "native_seed_identity_proved": seeds_passed,
        "native_seed_identity_worker_flag": seed_semantic,
        "actual_phase_partials_proved": phases_passed,
        "wrapper_mapping": root_wrappers,
    }


def _wrapper_mapping(packet: Mapping[str, Any]) -> dict[str, Any]:
    expected = {name: list(data["graph_lower_wrappers"]) for name, data in SEAM_DATA.items()}
    expected_upper = {name: list(data["upper_endpoint_R_Q_y"][:2]) for name, data in SEAM_DATA.items()}
    candidates = _find_alias(packet, ("graph_lower_wrapper_mapping", "lower_wrapper_mapping", "actual_graph_lower_wrappers"))
    mapping = candidates[0] if candidates and isinstance(candidates[0], Mapping) else None
    passed = True
    rows: dict[str, Any] = {}
    for name, wrappers in expected.items():
        raw = mapping.get(name) if isinstance(mapping, Mapping) else None
        if isinstance(raw, Mapping):
            actual_wrappers = raw.get("graph_lower_wrappers", raw.get("lower_wrappers", raw.get("wrappers")))
            actual_upper = raw.get("upper_endpoint_R_Q_y", raw.get("upper_R_Q_y"))
        else:
            actual_wrappers, actual_upper = None, None
        ok = actual_wrappers == wrappers and actual_upper == expected_upper[name]
        passed = passed and ok
        rows[name] = {"expected_graph_lower_wrappers": wrappers, "actual_graph_lower_wrappers": actual_wrappers,
                      "expected_upper_endpoint_R_Q_y": expected_upper[name], "actual_upper_endpoint_R_Q_y": actual_upper,
                      "passed": ok}
    return {"passed": passed, "rows": rows, "required": True}


def _control_packet_gate(packet: Any, expected: Mapping[str, Any], provenance: Mapping[str, Any], actual: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(packet, Mapping):
        return {"passed": False, "provenance": provenance, "seed_tables": {}, "reason": "control packet not provided",
                "native_generic_firstjets_C1": False, "positive_phase_controls": False,
                "implicit_function_theorem": False, "actual_memory_value_proved": False,
                "actual_phase_partials": False, "native_seed_theorem_binding_proved": False}
    seed_tables: dict[str, Any] = {}
    seeds_passed = True
    for chart in CHARTS:
        ok, detail = _packet_seed_rows(packet, chart, expected["seed_rows"][chart])
        seed_tables[chart] = detail
        seeds_passed = seeds_passed and ok
    native_seed_receipt = _control_native_seed_binding(packet, expected)
    native_seed = packet.get("native_seed_theorem_binding")
    # The control theorem must bind the actual C6/Q4 rows.  A free-standing
    # generic-first-jet boolean cannot discharge this gate.
    firstjet_formula = _control_firstjet_formula_receipt(packet)
    firstjet_semantic = bool(native_seed_receipt["passed"] and firstjet_formula["passed"])
    positive_packet = packet.get("positive_active_loop_IFT")
    positive_receipt = _positive_ift_binding(packet)
    hypotheses = positive_packet.get("hypotheses") if isinstance(positive_packet, Mapping) else None
    positivity = _proved_field(hypotheses, "lambda_positive")
    ift = _proved_field(hypotheses, "smooth_IFT")
    # The actual C5 phase-partial rows and inverse/lazy metadata are source
    # bound by this module's cache receipt.  The control worker's positive/IFT
    # packet need not duplicate those 21 graph rows.
    phase = bool(actual.get("phase_formula_receipt", {}).get("passed")
                 if isinstance(actual.get("phase_formula_receipt"), Mapping)
                 else actual.get("all_actual_phase_partial_rows_present"))
    memory = packet.get("actual_memory_value_proved") is True
    memory_rows = _find_alias(packet, ("actual_memory_rows", "memory_value_rows", "endpoint_memory_rows"))
    memory_rows_ok = False
    if memory_rows:
        # A value proof must carry the actual five rows on both joins.  Pointer
        # aliases alone do not qualify; the table is compared to the archive.
        memory_rows_ok = any(isinstance(row, Mapping) and all(name in row for name in CORRECTION_FIELDS)
                             for row in memory_rows)
    generic_rows_ok = _generic_rows_match(packet, expected)
    p2_receipt = _control_p2_naturality(packet, generic_rows_ok, parent)
    wrappers = _wrapper_mapping(packet)
    return {
        # Memory values and parent graph wrappers are composition-owned gates;
        # the worker packet need only discharge its native control theorem and
        # positive/IFT hypotheses here.
        "passed": bool(provenance.get("passed") and seeds_passed and firstjet_semantic and generic_rows_ok
                        and p2_receipt["passed"] and positive_receipt["passed"] and phase),
        "provenance": _copy(provenance),
        "seed_tables": seed_tables,
        "native_seed_theorem_binding": native_seed_receipt,
        "native_seed_theorem_binding_proved": native_seed_receipt["passed"],
        "control_firstjet_formula": firstjet_formula,
        "native_generic_firstjets_C1": bool(firstjet_semantic and generic_rows_ok),
        "p2_source_naturality": p2_receipt,
        "positive_phase_controls": positivity,
        "implicit_function_theorem": ift,
        "positive_IFT_binding": positive_receipt,
        "actual_phase_partials": phase,
        "actual_memory_value_proved": bool(memory and memory_rows_ok),
        "memory_rows_present": memory_rows_ok,
        "wrapper_mapping": wrappers,
    }


def _generic_rows_match(packet: Mapping[str, Any], expected: Mapping[str, Any]) -> bool:
    candidates = _find_alias(packet, ("firstjet_rows", "generic_firstjet_rows", "generic_rows", "actual_generic_firstjets", "source_graph"))
    if not candidates:
        return False
    # Require each chart to expose at least one actual C6 and Q4 row map.  The
    # detailed seed map above binds all field roots; this check binds the
    # universal generic functional maps used by the source theorem.
    for chart in CHARTS:
        chart_value = None
        for candidate in candidates:
            if isinstance(candidate, Mapping) and chart in candidate:
                chart_value = candidate[chart]
                break
            if isinstance(candidate, Mapping) and isinstance(candidate.get("charts"), Mapping) and chart in candidate["charts"]:
                chart_value = candidate["charts"][chart]
                break
        if not isinstance(chart_value, Mapping):
            return False
        if "C6" not in chart_value and "Q4" not in chart_value and "generic" in chart_value and "Q4_generic" in chart_value:
            chart_value = {"C6": chart_value.get("generic"), "Q4": chart_value.get("Q4_generic")}
        expected_c6 = expected["generic_rows"][chart]["C6"]
        expected_q4 = expected["generic_rows"][chart]["Q4"]
        for group, exp in (("C6", expected_c6), ("Q4", expected_q4)):
            got = chart_value.get(group)
            if not isinstance(got, Mapping):
                return False
            for key, rows in exp.items():
                value = got.get(key)
                if value != rows:
                    return False
    return True


def _control_firstjet_formula_receipt(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Check the issued universal a/b/phase functional map itself."""

    theorem = packet.get("native_seed_theorem_binding")
    jets = theorem.get("control_first_jets") if isinstance(theorem, Mapping) else None
    ordinary = jets.get("ordinary_Z_functions") if isinstance(jets, Mapping) else None
    physical = jets.get("physical_y_functions") if isinstance(jets, Mapping) else None
    endpoint = jets.get("endpoint_values") if isinstance(jets, Mapping) else None
    expected_ordinary = {
        "a": "(E - 2*E_y)/E",
        "b": "2*V_y/E",
        "t0": "-b/a",
        "Delta": "a - 2",
        "kappa": "a",
    }
    expected_physical = {
        "a_y": "(-2*E*E_yy + 2*E_y**2)/E**2",
        "b_y": "(2*E*V_yy - 2*E_y*V_y)/E**2",
    }
    expected_endpoint = {
        "a": "4/5", "a_y": "0", "b": "0", "b_y": "0",
        "t0": "0", "Delta": "-6/5", "kappa": "4/5",
    }
    checks = {
        "theorem_passed": isinstance(jets, Mapping) and jets.get("passed") is True,
        "ordinary_Z_functional_map": ordinary == expected_ordinary,
        "physical_y_functional_map": physical == expected_physical,
        "endpoint_control_rows": endpoint == expected_endpoint,
        "native_p2_uses_C5_generic": isinstance(jets, Mapping)
        and jets.get("native_p2") == "C5 generic_roots takes chart['generic']['p2']; Q4 modified p2 is not used",
    }
    return {"passed": all(checks.values()), "checks": checks,
            "ordinary_Z_functions": ordinary, "physical_y_functions": physical,
            "endpoint_values": endpoint, "source_values_selected": False}


def _control_p2_naturality(packet: Mapping[str, Any], generic_rows_ok: bool,
                           parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the native inlet p2 and p2_y source map to common FunctionJets.

    The control worker's first-jet table supplies the common C6/Q4 inputs;
    this receipt separately requires its issued rational p2 source and the
    complete physical-y expression, including the retained P0(Z), P0_Z(Z),
    and ordinary-Z terms.  No endpoint value is selected here.
    """

    source = packet.get("native_p2_source_identity")
    if not isinstance(source, Mapping):
        return {"passed": False, "reason": "native_p2_source_identity missing"}
    required_uses = {
        "E_y", "m_y", "m_zy", "e_y", "e_zy", "p_y", "p_zy", "P0_y=0", "P0_zy=0",
    }
    uses = source.get("p2_y_uses")
    p2_expression = source.get("p2_expression")
    p2_y_expression = source.get("p2_y_expression")
    expression_tokens = ("P0", "P0_z", "Pstar", "R", "E", "Z**2", "delta")
    y_expression_tokens = ("E_y", "V_y", "m_y", "m_zy", "e_y", "e_zy", "p_y", "p_zy", "P0_z")
    reachability = source.get("actual_graph_reachability")
    reach_rows = {
        chart: isinstance(reachability, Mapping)
        and isinstance(reachability.get(chart), Mapping)
        and reachability[chart].get("passed") is True
        and reachability[chart].get("source_executed") is False
        for chart in CHARTS
    }
    amplitude = parent.get("actual_micro_mean_Z3_bounds", {}).get("amplitude_provenance", {})
    background = parent.get("actual_long_reference_background_I1_numeric", {})
    radii = background.get("radii", {}) if isinstance(background, Mapping) else {}
    macro_background = parent.get("actual_macro_background_I1_numeric", {})
    geometric_measure = (macro_background.get("geometric_measure", {})
                         if isinstance(macro_background, Mapping) else {})
    pressure_contract = parent.get("actual_long_reference_pressure_source", {}).get("source_contract", {})
    composition = packet.get("native_p2_first_jet_endpoint_composition")
    composition_graphs = composition.get("per_chart_source_graphs") if isinstance(composition, Mapping) else None
    required_common_rows = [
        "E", "E_y", "m", "m_y", "m_Z", "m_yZ", "e", "e_y", "e_Z", "e_yZ",
        "p", "p_y", "p_Z", "p_yZ", "P0_Z",
    ]
    composition_chart_checks: dict[str, Any] = {}
    composition_graphs_ok = isinstance(composition_graphs, Mapping)
    for chart in ("long_reshape", "reference"):
        graph = composition_graphs.get(chart) if isinstance(composition_graphs, Mapping) else None
        row_ok = (isinstance(graph, Mapping)
                  and graph.get("passed") is True
                  and graph.get("source_executed") is False
                  and graph.get("p2_graph_reaches_P0_first_two") is True
                  and graph.get("P0_first_two_reachable") is True
                  and graph.get("p2_first_two_product_nodes") is not None)
        composition_chart_checks[chart] = {"passed": row_ok,
                                           "p2_graph_reaches_P0_first_two":
                                           graph.get("p2_graph_reaches_P0_first_two") if isinstance(graph, Mapping) else None,
                                           "P0_first_two_reachable":
                                           graph.get("P0_first_two_reachable") if isinstance(graph, Mapping) else None}
        composition_graphs_ok = composition_graphs_ok and row_ok
    composition_checks = {
        "receipt_present": isinstance(composition, Mapping),
        "receipt_passed": isinstance(composition, Mapping) and composition.get("passed") is True,
        "source_executed_false": isinstance(composition, Mapping) and composition.get("source_executed") is False,
        "same_rational_map": isinstance(composition, Mapping)
        and composition.get("same_rational_map") == "p2=R*(axial_linear + Pstar*axial_quadratic)/E",
        "Z_scope": isinstance(composition, Mapping)
        and composition.get("Z_scope") == {"lower": "37/100", "upper": "93/250", "source": "inherited parent scope"},
        "common_native_fields_and_first_y_verified": isinstance(composition, Mapping)
        and composition.get("common_native_fields_and_first_y_verified") is True,
        "ordinary_Z_source_rows_verified": isinstance(composition, Mapping)
        and composition.get("ordinary_Z_source_rows_verified") is True,
        "P0_Z_retained": isinstance(composition, Mapping) and composition.get("P0_Z_retained") is True,
        "P0_Zy_physical_y_zero": isinstance(composition, Mapping)
        and composition.get("P0_Zy_physical_y_zero") is True,
        "required_common_native_rows": isinstance(composition, Mapping)
        and composition.get("required_common_native_rows") == required_common_rows,
        "conditional_composed_jet_conclusion": isinstance(composition, Mapping)
        and composition.get("conditional_composed_jet_conclusion") is True,
        "endpoint_seam_equality_open": isinstance(composition, Mapping)
        and composition.get("endpoint_seam_equality_proved") is False
        and composition.get("parent_endpoint_seam_equality_proved") is False
        and composition.get("parent_endpoint_congruence_claimed") is False,
        "per_chart_graph_binding": composition_graphs_ok,
    }
    checks = {
        "source_identity_passed": source.get("passed") is True,
        "source_path": source.get("source") == INLET_SOURCE,
        "p2_formula": source.get("formula") == "p2 = R*(axial_linear + Pstar*axial_quadratic)/E",
        "p2_expression_complete": isinstance(p2_expression, str)
        and all(token in p2_expression for token in expression_tokens),
        "p2_y_expression_complete": isinstance(p2_y_expression, str)
        and all(token in p2_y_expression for token in y_expression_tokens),
        "p2_y_uses_complete": isinstance(uses, list) and required_uses.issubset(set(uses)),
        "P0_retained": source.get("full_P0_retained") is True,
        "ordinary_Z_FunctionJet_source": isinstance(source.get("ordinary_Z_function_jet_source"), str)
        and "FunctionJet" in source["ordinary_Z_function_jet_source"],
        "native_source_not_Q4_modified": source.get("native_generic_not_Q4_modified") is True,
        "no_fixture_subtraction": source.get("fixture_subtraction_used") is False,
        "common_functional_inputs": generic_rows_ok,
        "actual_graph_reachability": all(reach_rows.values()),
        # The rational map is only applied on the inherited source domain,
        # where its actual E and physical R factors are certified nonzero.
        "E_nonzero_on_inherited_domain": isinstance(amplitude, Mapping)
        and amplitude.get("passed") is True
        and amplitude.get("actual_positive_E0_is_the_same_native_E") is True
        and amplitude.get("F0_positive_and_at_most_one") is True
        and amplitude.get("native_phi_positive_and_nonincreasing") is True,
        "R_nonzero_on_inherited_domain": (
            isinstance(radii, Mapping) and radii.get("R_positive_for_0_lt_m_le_1") is True
        ) or (
            isinstance(geometric_measure, Mapping)
            and geometric_measure.get("R_positive_for_0_lt_m_le_1") is True
        ),
        "physical_R_source_binding": isinstance(pressure_contract, Mapping)
        and pressure_contract.get("physical_R_not_a_mean_scale") is True
        and pressure_contract.get("physical_slow_J_applied_once") is True,
        "native_p2_first_jet_endpoint_composition": all(composition_checks.values()),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "source": source.get("source"),
        "p2_formula": source.get("formula"),
        "p2_y_uses": uses,
        "common_functional_inputs": generic_rows_ok,
        "graph_reachability": reach_rows,
        "p2_first_jet_endpoint_composition": {
            "passed": all(composition_checks.values()),
            "checks": composition_checks,
            "charts": composition_chart_checks,
            "source_executed": False,
        },
        "source_values_selected": False,
    }


def _control_native_seed_binding(packet: Mapping[str, Any], expected: Mapping[str, Any]) -> dict[str, Any]:
    """Check the control theorem's field binding table against cache rows."""

    theorem = packet.get("native_seed_theorem_binding")
    bindings = theorem.get("field_bindings") if isinstance(theorem, Mapping) else None
    rows: dict[str, Any] = {}
    all_rows = isinstance(bindings, Mapping)
    for field in FIELDS:
        field_receipt: dict[str, Any] = {"passed": False, "charts": {}}
        item = bindings.get(field) if isinstance(bindings, Mapping) else None
        item_ok = (isinstance(item, Mapping) and item.get("fieldname") == field
                   and item.get("ordinaryorder") == list(range(7))
                   and item.get("ordinaryorder_rows") == list(range(7))
                   and item.get("source_executed") is False)
        for chart in CHARTS:
            seed_rows = expected["seed_rows"][chart][field]
            source_expected = [row["source_field_node"] for row in seed_rows]
            radial_expected = [row["physical_radial_row_node"] for row in seed_rows]
            refs = item.get("bgfieldref", {}).get(chart) if isinstance(item, Mapping) else None
            radial = item.get("radial_bgref", {}) if isinstance(item, Mapping) else {}
            chart_rows = radial.get("chart_rows", {}).get(chart) if isinstance(radial, Mapping) else None
            source_rows = refs.get("field_rows") if isinstance(refs, Mapping) else None
            source_chart = refs.get("source_chart") if isinstance(refs, Mapping) else None
            J = seed_rows[0]["physical_J_node"] if seed_rows else None
            # The frozen control worker keeps the actual Q4 rows under
            # ``chart_rows[chart].actual_cache_rows``.  Accept that issued
            # shape while retaining the compact direct shape used by older
            # review packets; both are compared against the same cache rows.
            row_payload = chart_rows.get("actual_cache_rows") if isinstance(chart_rows, Mapping) else None
            if not isinstance(row_payload, Mapping):
                row_payload = chart_rows
            radial_rows = row_payload.get("Q4_radial_background") if isinstance(row_payload, Mapping) else None
            row_J = chart_rows.get("Jhandle") if isinstance(chart_rows, Mapping) else None
            chart_ok = (isinstance(refs, Mapping)
                        and source_chart == C6_CHART[chart]
                        and refs.get("J_handle") == J
                        and source_rows == source_expected
                        and isinstance(chart_rows, Mapping)
                        and row_J == J
                        and radial_rows == radial_expected
                        and row_J == radial.get("chart_rows", {}).get(chart, {}).get("Jhandle"))
            field_receipt["charts"][chart] = {
                "source_chart": source_chart,
                "source_rows": source_rows,
                "expected_source_rows": source_expected,
                "radial_rows": radial_rows,
                "expected_radial_rows": radial_expected,
                "J_handle": refs.get("J_handle") if isinstance(refs, Mapping) else None,
                "expected_J_handle": J,
                "passed": chart_ok,
            }
            item_ok = item_ok and chart_ok
        field_receipt["passed"] = item_ok
        rows[field] = field_receipt
        all_rows = all_rows and item_ok
    theorem_ok = (isinstance(theorem, Mapping) and theorem.get("passed") is True
                  and isinstance(theorem.get("control_first_jets"), Mapping)
                  and theorem["control_first_jets"].get("passed") is True)
    return {"passed": bool(theorem_ok and all_rows), "theorem_flag": theorem_ok,
            "field_bindings": rows}


def _positive_ift_binding(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Require the finite source-parameter handles behind the control proof."""

    positive = packet.get("positive_active_loop_IFT")
    if not isinstance(positive, Mapping):
        return {"passed": False, "handles": None, "reason": "positive_active_loop_IFT missing"}
    # A nested hypothesis is not allowed to promote a packet whose own
    # source receipt is stale or incomplete.  This is deliberately checked
    # before inspecting the parameter handles below.
    packet_passed = positive.get("passed") is True
    handles = positive.get("source_parameter_handles", positive.get("source_parameters"))
    if not isinstance(handles, Mapping):
        return {"passed": False, "handles": None, "reason": "source eta/dstar handles missing"}
    handles_ok = (handles.get("N") == 18 and handles.get("eta") == 21
                  and handles.get("dstar") == 22)
    hypotheses = positive.get("hypotheses")
    lambda_ok = _proved_field(hypotheses, "lambda_positive")
    ift_ok = _proved_field(hypotheses, "smooth_IFT")
    # The finite p2 control is an issued source receipt, not a free symbol
    # assumption.  The current worker exposes it under hypotheses.finite_p2;
    # retain the concrete finite-log and accepted-parameter records as part
    # of the same receipt so a future packet cannot pass by changing only a
    # scalar boolean.
    finite_p2 = hypotheses.get("finite_p2") if isinstance(hypotheses, Mapping) else None
    finite_p2_ok = (isinstance(finite_p2, Mapping)
                    and (finite_p2.get("proved") is True or finite_p2.get("passed") is True)
                    and finite_p2.get("value_not_selected") is True
                    and isinstance(finite_p2.get("basis"), str)
                    and "p2" in finite_p2["basis"].lower())
    finite_logs = positive.get("finite_log_records")
    eta_log = finite_logs.get("eta") if isinstance(finite_logs, Mapping) else None
    dstar_log = finite_logs.get("dstar") if isinstance(finite_logs, Mapping) else None
    finite_log_ok = all(
        isinstance(row, Mapping)
        and row.get("finite") is True
        and row.get("directed") is True
        and row.get("exp_strictly_positive") is True
        and row.get("canonical_source_handle") == handle
        for row, handle in ((eta_log, 21), (dstar_log, 22))
    )
    accepted_records = positive.get("accepted_positive_records")
    parameter_record = (accepted_records.get("q_mixed_parameter_envelopes")
                        if isinstance(accepted_records, Mapping) else None)
    parameter_record_ok = (isinstance(parameter_record, Mapping)
                           and parameter_record.get("passed") is True
                           and parameter_record.get("source_executed") is False)
    native_p2 = positive.get("native_p2_finite_receipt")
    native_p2_charts = native_p2.get("charts") if isinstance(native_p2, Mapping) else None
    native_p2_chart_receipts: dict[str, Any] = {}
    native_p2_ok = (isinstance(native_p2, Mapping)
                    and native_p2.get("passed") is True
                    and native_p2.get("source_executed") is False
                    and native_p2.get("value_not_selected") is True
                    and isinstance(native_p2_charts, Mapping))
    # The worker emits only the two actual seam charts for this finite source
    # receipt.  Check every p2 cap and both positive-E denominator logs while
    # retaining only compact row metadata in our packet (the exact MP limbs
    # remain in the worker's source-bound receipt).
    for chart in ("long_reshape", "reference"):
        chart_receipt = native_p2_charts.get(chart) if isinstance(native_p2_charts, Mapping) else None
        caps = chart_receipt.get("p2_log_caps") if isinstance(chart_receipt, Mapping) else None
        positive_e = chart_receipt.get("positive_E_source") if isinstance(chart_receipt, Mapping) else None
        cap_checks = []
        if isinstance(caps, list):
            for order, row in enumerate(caps):
                cap_checks.append({
                    "ordinary_Z_order": row.get("ordinary_Z_order") if isinstance(row, Mapping) else None,
                    "finite": row.get("finite") is True if isinstance(row, Mapping) else False,
                    "directed": row.get("directed") is True if isinstance(row, Mapping) else False,
                    "value_not_selected": row.get("value_not_selected") is True if isinstance(row, Mapping) else False,
                    "order_matches": row.get("ordinary_Z_order") == order if isinstance(row, Mapping) else False,
                })
        e_checks = []
        if isinstance(positive_e, Mapping):
            for component in ("phi", "Utheta"):
                component_receipt = positive_e.get(component)
                component_ok = isinstance(component_receipt, Mapping)
                component_rows = {}
                for side in ("lower", "upper"):
                    row = component_receipt.get(side) if isinstance(component_receipt, Mapping) else None
                    row_ok = (isinstance(row, Mapping) and row.get("finite") is True
                              and row.get("directed") is True
                              and row.get("value_not_selected") is True)
                    component_rows[side] = row_ok
                    component_ok = component_ok and row_ok
                e_checks.append({"component": component, "passed": component_ok,
                                 "rows": component_rows})
        chart_ok = (isinstance(chart_receipt, Mapping)
                    and chart_receipt.get("source_executed") is False
                    and chart_receipt.get("value_not_selected") is True
                    and chart_receipt.get("ordinary_Z_orders") == [0, 1, 2, 3]
                    and isinstance(caps, list) and len(caps) == 4
                    and all(item.get("finite") is True
                            and item.get("directed") is True
                            and item.get("value_not_selected") is True
                            and item.get("order_matches") is True
                            for item in cap_checks)
                    and len(e_checks) == 2 and all(item.get("passed") is True for item in e_checks))
        native_p2_chart_receipts[chart] = {
            "p2_graph_rows": chart_receipt.get("p2_graph_rows") if isinstance(chart_receipt, Mapping) else None,
            "p2_cap_rows": cap_checks,
            "positive_E_rows": e_checks,
            "passed": chart_ok,
        }
        native_p2_ok = native_p2_ok and chart_ok
    # The control worker records a local IFT receipt in the smooth hypothesis;
    # accepting either spelling keeps the API stable without accepting a bare
    # top-level boolean.
    local_ift = (isinstance(hypotheses, Mapping)
                 and isinstance(hypotheses.get("smooth_IFT"), Mapping)
                 and (hypotheses["smooth_IFT"].get("local_IFT") is True
                      or hypotheses["smooth_IFT"].get("proved") is True))
    return {"passed": bool(packet_passed and handles_ok and lambda_ok and ift_ok
                            and local_ift and finite_p2_ok and finite_log_ok
                            and parameter_record_ok and native_p2_ok),
            "packet_passed": packet_passed,
            "handles": dict(handles), "handles_ok": handles_ok,
            "lambda_positive": lambda_ok, "smooth_IFT": ift_ok,
            "local_IFT": local_ift,
            "finite_p2_source_receipt": {
                "passed": bool(finite_p2_ok and native_p2_ok),
                "hypothesis": _copy(finite_p2),
                "value_not_selected": finite_p2.get("value_not_selected") if isinstance(finite_p2, Mapping) else None,
                "native_receipt_present": isinstance(native_p2, Mapping),
                "native_receipt_passed": native_p2_ok,
                "charts": native_p2_chart_receipts,
            },
            "finite_parameter_log_receipt": {
                "passed": finite_log_ok,
                "eta": _copy(eta_log),
                "dstar": _copy(dstar_log),
            },
            "accepted_parameter_envelope_receipt": {
                "passed": parameter_record_ok,
                "record": _copy(parameter_record),
            }}


def _memory_gate_from_cache(actual: Mapping[str, Any], parent: Mapping[str, Any],
                            archives: Mapping[str, Any]) -> dict[str, Any]:
    """Discharge endpoint memory values from the accepted FTC/q0 source lemma.

    A predecessor pointer alone is kept as a separate receipt.  The value
    identity needs the actual C5 transport source theorem and its q0 facts:
    ``exp(0)=1`` and the zero-length integral is zero.  Those facts make each
    lower history equal its incoming pointer, including the rate-zero p row.
    """

    history = parent.get("actual_long_reference_history_seams", {})
    conditional = history.get("conditional_history_ODE", {}) if isinstance(history, Mapping) else {}
    c5_source = history.get("c5_history_source", {}) if isinstance(history, Mapping) else {}
    q0 = conditional.get("q0", {}) if isinstance(conditional, Mapping) else {}
    q0_ok = (conditional.get("passed") is True and q0.get("correction") is True
             and q0.get("exp_zero_is_one") is True and q0.get("integral_zero_is_zero") is True)
    c5_ok = (c5_source.get("passed") is True
             and c5_source.get("p_rate_zero_memory_retained") is True)
    c5charts_parent = c5_source.get("charts", {}) if isinstance(c5_source, Mapping) else {}
    c5charts = archives["c5"]["prescribed_loop_C5_bindings"]["charts"]
    rows: dict[str, Any] = {}
    all_rows = True
    for chart in ("long_reshape", "reference"):
        parent_chart = c5charts_parent.get(chart) if isinstance(c5charts_parent, Mapping) else None
        cache_chart = c5charts.get(chart)
        history_rows = parent_chart.get("histories", {}) if isinstance(parent_chart, Mapping) else {}
        chart_rows: dict[str, Any] = {}
        chart_ok = isinstance(parent_chart, Mapping) and parent_chart.get("passed") is True
        for field in CORRECTION_FIELDS:
            p_hist = history_rows.get(field) if isinstance(history_rows, Mapping) else None
            c_in = cache_chart.get("incoming", {}).get(field) if isinstance(cache_chart, Mapping) else None
            c_hist = cache_chart.get("history", {}).get(field) if isinstance(cache_chart, Mapping) else None
            q = p_hist.get("q0", {}) if isinstance(p_hist, Mapping) else {}
            metadata = p_hist.get("integral_metadata", {}) if isinstance(p_hist, Mapping) else {}
            row_ok = (isinstance(p_hist, Mapping) and p_hist.get("incoming") == c_in
                      and p_hist.get("history") == c_hist
                      and q.get("exp_zero_is_one") is True
                      and q.get("integral_zero_is_zero") is True
                      and q.get("incoming_memory_retained") is True
                      and metadata.get("physical_Jacobian_applied_once") is True
                      and metadata.get("predecessor_memory_never_reset") is True)
            chart_ok = chart_ok and row_ok
            chart_rows[field] = {"incoming": c_in, "history": c_hist,
                                 "q0": q, "integral_metadata": metadata, "passed": row_ok}
        # p has rate zero and must retain the incoming memory rather than be
        # replaced by a fresh zero.  The source parent and cache both expose
        # its exact unit memory handle.
        p_memory = cache_chart.get("memory", {}).get("p") if isinstance(cache_chart, Mapping) else None
        chart_rows["p_rate_zero"] = {"memory_handle": p_memory,
                                      "parent_proof": c5_ok and p_memory == 1,
                                      "passed": c5_ok and p_memory == 1}
        chart_ok = chart_ok and chart_rows["p_rate_zero"]["passed"]
        rows[chart] = {"fields": chart_rows, "passed": chart_ok}
        all_rows = all_rows and chart_ok
    pointer_rows = bool(actual.get("all_nonreset_history_pointers_bound"))
    passed = bool(pointer_rows and q0_ok and c5_ok and all_rows)
    return {
        "passed": passed,
        "pointer_rows_bound": pointer_rows,
        "actual_value_proof": passed,
        "q0_lemma": {"passed": q0_ok, "receipt": q0},
        "c5_history_source": {"passed": c5_ok, "p_rate_zero_memory_retained": c5_source.get("p_rate_zero_memory_retained")},
        "rows": rows,
        "source_executed": False,
    }


def _symbolic_qy_check() -> dict[str, Any]:
    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover - dependency diagnostic
        return {"passed": False, "error": f"SymPy unavailable: {exc}"}
    z, delta, eps = sp.symbols("Z delta eps")
    U = sp.Function("U")(z)
    m = sp.Function("m")(z)
    Uy = sp.Function("Uy")(z)
    my = sp.Function("my")(z)
    L = 1 - delta * z**2
    numerator = lambda u, mean: 2 * z * u - ((1 - delta) * z * mean + (1 - z**2) * sp.diff(mean, z))
    q_eps = numerator(U + eps * Uy, m + eps * my) / L
    qy = sp.diff(q_eps, eps).subs(eps, 0)
    issued = numerator(Uy, my) / L
    quotient_residual = sp.factor(qy - issued)
    commutation_residual = sp.simplify(sp.diff(qy, z) - sp.diff(sp.diff(q_eps, z), eps).subs(eps, 0))
    return {
        "passed": quotient_residual == 0 and commutation_residual == 0,
        "L": "1-delta*Z**2",
        "issued_Qy": "(2*Z*Uy-((1-delta)*Z*my+(1-Z**2)*D_Z(my)))/L",
        "quotient_residual": str(quotient_residual),
        "ordinary_Z_commutation_residual": str(commutation_residual),
        "ordinary_Z_orders": [0, 1],
        "source_tool": "SymPy exact symbols; no archived source execution",
    }


def _field_identities(full: bool) -> dict[str, Any]:
    status = "proved under all source gates" if full else "conditional; endpoint value/jet gate remains open"
    return {
        "status": status,
        "U": "Delta U = Pstar*(Delta V_N)",
        "Uy": "Delta Uy = Pstar*Delta(V_y + slow_B_y/N + B_phi)",
        "m": "Delta m = Pstar*Delta(m_native + m_correction_transport)",
        "my": "Delta my = Delta U - Delta m",
        "Qy": "Delta Qy = Delta[(2*Z*Uy-((1-delta)*Z*my+(1-Z**2)*D_Z(my))/(1-delta*Z**2))]",
        "ordinary_Z1": "D_Z commutes with the displayed rational source operator",
        "unresolved_if_conditional": [] if full else ["native seed/value identities and worker source controls"],
    }


def _seams(full: bool, actual: Mapping[str, Any], wrapper_ok: bool,
           radius_ok: bool, symbolic: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    operator_residual = (symbolic.get("quotient_residual")
                         if symbolic.get("passed") is True else None)
    commutation_residual = (symbolic.get("ordinary_Z_commutation_residual")
                            if symbolic.get("passed") is True else None)
    for name, data in SEAM_DATA.items():
        trace = actual["rqy_endpoint_traces"][name]
        seam_full = bool(full and wrapper_ok and radius_ok
                         and trace["cache_endpoint_trace_matches_contract"])
        result[name] = {
            "upper_chart": data["upper_chart"],
            "lower_chart": data["lower_chart"],
            "ordinary_Z_orders": [0, 1],
            "upper_endpoint_R_Q_y": trace["upper_endpoint_R_Q_y"],
            "lower_chart_R_Q_y": trace["lower_chart_R_Q_y"],
            "graph_lower_wrappers": trace["graph_lower_wrappers"],
            "left_RQy_upper": trace["upper_endpoint_R_Q_y"][:2],
            "right_RQy_lower_wrapper": trace["graph_lower_wrappers"],
            "graph_lower_wrapper_to_upper_endpoint": list(zip(trace["graph_lower_wrappers"], trace["upper_endpoint_R_Q_y"][:2])),
            "field_identities": _field_identities(seam_full),
            "finite_N_Qy_identity_proven": seam_full,
            "physical_Qy_C1_equal": seam_full,
            "physical_RQy_C1_equal": seam_full,
            "Qy_operator_symbolic_residual": operator_residual,
            "ordinary_Z1_commutation_residual": commutation_residual,
            "source_values_selected": False,
            "interval_or_cap_values_used": False,
            "N_materialized": False,
            "passed": seam_full,
        }
    return result


def build(root: str | Path, parent: Mapping[str, Any] | str | Path | None = None,
          control_packet: Mapping[str, Any] | None = None,
          derivative_packet: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build the source-bound conditional full modified seam receipt.

    ``control_packet`` and ``derivative_packet`` are receipts from the two
    upstream source workers.  They are data only; this function never imports
    or calls their modules.  Missing packets intentionally produce an honest
    conditional result rather than a fabricated full seam.
    """
    root_path = Path(root)
    _, disk_parent = _load_parent(root_path, parent)
    family = disk_parent["source_family"]
    pins = _source_pins(root_path)
    guards = _source_guards(root_path)
    archives = _load_archives(root_path, family)
    actual = _graph_port_receipt(archives)
    qy_denominator = _qy_denominator_receipt(archives, disk_parent)
    phase_formula = _phase_formula_receipt(archives)
    actual["phase_formula_receipt"] = phase_formula
    actual["all_actual_phase_formula_rows_bound"] = phase_formula["passed"]
    wrapper_receipt = _parent_graph_wrapper_receipt(disk_parent, actual)
    radius_receipt = _parent_endpoint_radius_receipt(disk_parent)
    axial_scope_receipt = _parent_axial_scope_receipt(disk_parent)
    my_receipt = _modified_my_source_receipt(archives, disk_parent)
    parent_gate = _parent_lemma_gate(disk_parent)
    symbolic = _symbolic_qy_check()
    symbolic["actual_Qy_graph_rows_bound"] = actual.get("all_R_Q_y_source_rows_bound") is True
    symbolic["actual_pressure_operator_bound"] = actual.get("passed") is True
    derivative_provenance = _packet_provenance(root_path, derivative_packet, "derivative", "src/openai_ns_reconstruction/original_physical_y_dag_derivative.py")
    control_provenance = _packet_provenance(root_path, control_packet, "control", "src/openai_ns_reconstruction/original_long_reference_control_seams.py")
    # Worker packets are checked against the exact cache-derived seed maps and
    # phase rows.  The packet itself can add stronger conditions, but cannot
    # replace these data-bound checks with a free-standing boolean.
    expected_for_workers = {
        "charts": {},
        "seed_rows": actual["seed_rows"],
        "generic_rows": actual["generic_rows"],
    }
    for chart in CHARTS:
        c = actual["charts"][chart]
        c5chart = archives["c5"]["prescribed_loop_C5_bindings"]["charts"][chart]
        expected_for_workers["charts"][chart] = {
            "loop_source_roots": {"A": list(c5chart["loop"]["A"]), "B_over_Pstar": list(c5chart["loop"]["B_over_Pstar"])},
            "slow_derivative_roots": {"A": list(archives["q4"]["modified_Q4_pressure_C3_bindings"]["charts"][chart]["radial"]["slow_A_y"]),
                                       "B_over_Pstar": list(archives["q4"]["modified_Q4_pressure_C3_bindings"]["charts"][chart]["radial"]["slow_B_y"])},
            "phase_partial_coefficients": list(c["phase_partial_coefficients"]),
            "phase_partial_graph_bound": phase_formula["charts"][chart]["phase_partial_graph"]["passed"],
        }
    derivative_gate = _derivative_packet_gate(derivative_packet, expected_for_workers, derivative_provenance)
    control_gate = _control_packet_gate(control_packet, expected_for_workers, control_provenance, actual, disk_parent)
    # The accepted history source theorem discharges endpoint memory values
    # from q0/FTC; no worker memory boolean is required.  The saved parent
    # graph patch independently discharges the two lower-wrapper substitutions.
    memory_gate = _memory_gate_from_cache(actual, disk_parent, archives)
    memory_gate["worker_actual_value_proof"] = control_gate.get("actual_memory_value_proved") is True
    wrapper_ok = wrapper_receipt.get("passed") is True
    control_seed_rows_ok = bool(control_gate.get("seed_tables")
                                and all(detail.get("passed") is True
                                        for detail in control_gate["seed_tables"].values()))
    gate_map = {
        "parent_source_bound_lemmas": parent_gate["passed"],
        "source_guards_and_actual_cache_ports": actual["passed"] and all(pin["passed"] for pin in pins.values()) and bool(guards),
        "Qy_denominator_positive_source_binding": qy_denominator["passed"],
        "symbolic_Qy_operator_and_ordinary_Z_commutation": symbolic["passed"],
        # Keep the packet-level receipt itself in the admission graph.  The
        # individual semantic gates below remain separately visible, but a
        # stale/partially malformed packet must never be promoted by a few
        # top-level ``True`` aliases.
        "derivative_source_packet_binding": derivative_gate.get("passed") is True,
        "control_source_packet_binding": control_gate.get("passed") is True,
        "native_seed_map_correct": bool(derivative_gate.get("native_seed_identity_proved")
                                        and control_seed_rows_ok
                                        and control_gate.get("native_seed_theorem_binding_proved") is True
                                        and derivative_gate.get("provenance", {}).get("passed")
                                        and control_gate.get("provenance", {}).get("passed")),
        "actual_36_held_phase_derivative_roots": derivative_gate.get("passed") is True,
        "native_generic_firstjets_C1": control_gate.get("native_generic_firstjets_C1") is True,
        "control_p2_source_naturality": control_gate.get("p2_source_naturality", {}).get("passed") is True,
        "positive_phase_controls_and_IFT": bool(control_gate.get("positive_phase_controls")
                                                  and control_gate.get("implicit_function_theorem")
                                                  and control_gate.get("positive_IFT_binding", {}).get("passed") is True),
        "actual_phase_partials": bool(phase_formula.get("passed")
                                       and derivative_gate.get("actual_phase_partials_proved")
                                       and control_gate.get("actual_phase_partials")),
        "actual_nonreset_memory_values": memory_gate["passed"],
        "actual_modified_my_source_identity": my_receipt["passed"],
        "actual_graph_lower_wrapper_mapping": wrapper_ok,
        "same_physical_radius_endpoint_identity": radius_receipt["passed"],
        "certified_axial_Z_scope": axial_scope_receipt["passed"],
    }
    full = bool(all(gate_map.values()))
    unresolved = [name for name, passed in gate_map.items() if not passed]
    seams = _seams(full, actual, wrapper_ok, radius_receipt["passed"], symbolic)
    module_sha = _sha(_locate(root_path, MODULE))
    result: dict[str, Any] = {
        "passed": bool(gate_map["source_guards_and_actual_cache_ports"] and gate_map["symbolic_Qy_operator_and_ordinary_Z_commutation"]),
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": _copy(family),
        "source_pins": pins,
        "input_pins": pins,
        "source_guards": guards,
        "sourceguards_tuple": True,
        "accepted_ancestors_replayed": False,
        "accepted_builders_executed": False,
        "accepted_checkers_executed": False,
        "accepted_producers_executed": False,
        "owners_executed": False,
        "norms_executed": False,
        "recurrences_executed": False,
        "quadrature_executed": False,
        "ordinary_Z_orders": [0, 1],
        "scope": "full modified U, Uy, m, my, Qy ordinary-C1 endpoint identities on post_power/long_reshape/reference within Z in [37/100,93/250]",
        "axial_scope": _copy(AXIAL_SCOPE),
        "parent_lemma_gate": parent_gate,
        "actual_source_composition": actual,
        "qy_denominator_receipt": qy_denominator,
        "derivative_packet": derivative_gate,
        "control_packet": control_gate,
        "memory_value_gate": memory_gate,
        "modified_my_source_receipt": my_receipt,
        "graph_lower_wrapper_receipt": wrapper_receipt,
        "endpoint_radius_receipt": radius_receipt,
        "axial_scope_receipt": axial_scope_receipt,
        "symbolic_checks": symbolic,
        "gates": gate_map,
        "unresolved_gates": unresolved,
        "conditional_theorem": {
            "hypotheses": [
                "the exact C6/C5/Q4 source pins and graph prefixes remain unchanged",
                "the derivative worker proves all 36 actual held-phase A/B roots and the exact C6-to-Q4 seed maps",
                "the control worker proves the same universal generic first-jet map, positive phase lambda, and IFT smoothness",
                "the accepted C5 FTC/q0 history lemma proves the non-reset endpoint memory values",
                "the accepted 185690-node parent graph patch proves the exact lower-wrapper substitutions",
            ],
            "conclusion": "under all hypotheses on the certified axial Z window [37/100,93/250], both displayed seams have identical finite-N Qy and physical Qy/RQy ordinary-C1 endpoint jets",
            "ordinary_Z_orders": [0, 1],
            "axial_scope": _copy(AXIAL_SCOPE),
            "status": "proved" if full else "conditional; unresolved gates are listed explicitly",
            "field_identities": _field_identities(full),
        },
        "seams": seams,
        "q4_physical_port_roots": actual["physical_port_roots"],
        "RQy_lower_upper_traces": actual["rqy_endpoint_traces"],
        "full_modified_seam_proved": full,
        "full_modified_seam_claimed": full,
        "limits": {
            "axial_scope": _copy(AXIAL_SCOPE),
            "modified_Uy_Qy_closed": full,
            "natural_B_y_closed": full,
            "q4_modified_Qy_seam_proved": full,
            "full_modified_Qy_claimed": full,
            "source_values_selected": False,
            "interval_or_cap_values_used": False,
            "N_materialized": False,
            "accepted_builder_checker_owner_producer_norm_recurrence_quadrature_executed": False,
            "graph_contraction_performed": False,
            "raw_graph_retained": True,
        },
        "module_sha256": module_sha,
    }
    result["receipt_sha256"] = _json_digest(result)
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["build"]
